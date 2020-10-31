#!python
#cython: boundscheck=False
#cython: wraparound=True
#cython: initializedcheck=False
#cython: cdivision=True


#import pylab as plt
import numpy as np
cimport numpy as np
import netCDF4 as nc
from scipy.interpolate import pchip_interpolate
from libc.math cimport pow, cbrt, exp, fmin, fmax
# from thermodynamic_functions cimport cpm_c
include 'parameters.pxi'
from profiles import profile_data

cimport ReadProfiles
# cimport TimeStepping
# from NetCDFIO cimport NetCDFIO_Stats
# from mlm_thermodynamic_functions import *
# from mlm_thermodynamic_functions cimport *
import cPickle
import pickle as pkl
from cfsites_forcing_reader import cfreader


# Note: the RRTM modules are compiled in the 'RRTMG' directory:
cdef extern:
    void c_rrtmg_lw_init(double *cpdair)
    void c_mcica_subcol_lw (int *iplon, int *ncol, int *nlay, int *icld, int *permuteseed, int *irng, 
             double *play, double *cldfrac, double *ciwp, double *clwp, double *rei, double *rel, double *tauc, 
             double *cldfmcl, double *ciwpmcl, double *clwpmcl, double *reicmcl, double *relqmcl, double *taucmcl)
    void c_rrtmg_lw (int *ncol    ,int *nlay    ,int *icld    ,int *idrv    ,
             double *play    ,double *plev    ,double *tlay    ,double *tlev    ,double *tsfc    ,
             double *h2ovmr  ,double *o3vmr   ,double *co2vmr  ,double *ch4vmr  ,double *n2ovmr  ,double *o2vmr,
             double *cfc11vmr,double *cfc12vmr,double *cfc22vmr,double *ccl4vmr ,double *emis    ,
             int *inflglw ,int *iceflglw,int *liqflglw,double *cldfr   ,
             double *taucld  ,double *cicewp  ,double *cliqwp  ,double *reice   ,double *reliq   ,
             double *tauaer  ,
             double *uflx    ,double *dflx    ,double *hr      ,double *uflxc   ,double *dflxc,  double *hrc,
             double *duflx_dt,double *duflxc_dt)
    void c_rrtmg_sw_init(double *cpdair)
    void c_mcica_subcol_sw(int *iplon, int *ncol, int *nlay, int *icld, int *permuteseed, int *irng, 
             double *play, double *cldfrac, double *ciwp, double *clwp, double *rei, double *rel, 
             double *tauc, double *ssac, double *asmc, double *fsfc,
             double *cldfmcl, double *ciwpmcl, double *clwpmcl, double *reicmcl, double *relqmcl,
             double *taucmcl, double *ssacmcl, double *asmcmcl, double *fsfcmcl)
    void c_rrtmg_sw (int *ncol    ,int *nlay    ,int *icld    ,int *iaer    ,
             double *play    ,double *plev    ,double *tlay    ,double *tlev    ,double *tsfc    ,
             double *h2ovmr  ,double *o3vmr   ,double *co2vmr  ,double *ch4vmr  ,double *n2ovmr  ,double *o2vmr,
             double *asdir   ,double *asdif   ,double *aldir   ,double *aldif   ,
             double *coszen  ,double *adjes   ,int *dyofyr  ,double *scon    ,
             int *inflgsw ,int *iceflgsw,int *liqflgsw,double *cldfr   ,
             double *taucld  ,double *ssacld  ,double *asmcld  ,double *fsfcld  ,
             double *cicewp  ,double *cliqwp  ,double *reice   ,double *reliq   ,
             double *tauaer  ,double *ssaaer  ,double *asmaer  ,double *ecaer   ,
             double *swuflx  ,double *swdflx  ,double *swhr    ,double *swuflxc ,double *swdflxc ,double *swhrc)


cdef class Radiation:
    def __init__(self, namelist):

        # Required for surface energy budget calculations, can also be used for stats io
        self.srf_lw_down = 0.0
        self.srf_sw_down = 0.0
        self.srf_lw_up = 0.0
        self.srf_sw_down = 0.0

        try:
            self.profile_name = namelist['radiation']['profile_name']
        except:
            self.profile_name = 'cgils_ctl_s6'

        try:
            self.read_file = namelist['radiation']['read_file']
        except:
            self.read_file = False
        if self.read_file:
            self.file = str(namelist['radiation']['file'])
            self.site = namelist['radiation']['site']

        try:
            self.n_buffer = namelist['radiation']['n_buffer']
        except:
            self.n_buffer = 4

        try:
            self.stretch_factor = namelist['radiation']['stretch_factor']
        except:
            self.stretch_factor = 1.0

        try:
            self.patch_pressure = namelist['radiation']['patch_pressure']
        except:
            self.patch_pressure = 1000.00*100.0

        # Namelist options related to gas concentrations
        try:
            self.co2_factor = namelist['radiation']['co2_factor']
        except:
            self.co2_factor = 1.0
        try:
            self.h2o_factor = namelist['radiation']['h2o_factor']
        except:
            self.h2o_factor = 1.0

        # Namelist options related to insolation
        try:
            self.dyofyr = namelist['radiation']['dyofyr']
        except:
            self.dyofyr = 0
        try:
            self.adjes = namelist['radiation']['adjes']
        except:
            print('Insolation adjustive factor not set so RadiationRRTM takes default value: adjes = 0.5 (12 hour of daylight).')
            self.adjes = 0.5

        try:
            self.scon = namelist['radiation']['solar_constant']
        except:
            print('Solar Constant not set so RadiationRRTM takes default value: scon = 1360.0 .')
            self.scon = 1360.0
        try:
            self.toa_sw = namelist['radiation']['toa_sw']
        except:
            print('TOA shortwave not set so RadiationRRTM takes default value: toa_sw = 420.0 .')
            self.toa_sw = 420.0
        if self.read_file:
            rdr = cfreader(self.file, self.site)
            self.toa_sw = rdr.get_timeseries_mean('rsdt')

        try:
            self.coszen = namelist['radiation']['coszen']
        except:
            if (self.toa_sw > 0.0):
                self.coszen = self.toa_sw / (self.scon * self.adjes)
            else:
                print('Mean Daytime cos(SZA) not set so RadiationRRTM takes default value: coszen = 2.0/pi .')
                self.coszen = 2.0/pi

        try:
            self.adif = namelist['radiation']['adif']
        except:
            print('Surface diffusive albedo not set so RadiationRRTM takes default value: adif = 0.06 .')
            self.adif = 0.06

        try:
            self.adir = namelist['radiation']['adir']
        except:
            if (self.coszen > 0.0):
                self.adir = (.026/(self.coszen**1.7+.065) + (.15*(self.coszen-0.10)*(self.coszen-0.50)*(self.coszen-1.00)))
            else:
                self.adir = 0.0
            print('Surface direct albedo not set so RadiationRRTM computes value: adir = %5.4f .'%(self.adir))

        try:
            self.uniform_reliq = namelist['radiation']['uniform_reliq']
        except:
            print('uniform_reliq not set so RadiationRRTM takes default value: uniform_reliq = False.')
            self.uniform_reliq = False

        try:
            self.radiation_frequency = namelist['radiation']['frequency']
        except:
            print('radiation_frequency not set so RadiationRRTM takes default value: radiation_frequency = 0.0 (compute at every step).')
            self.radiation_frequency = 60.0

        self.next_radiation_calculate = 0.0

        try:
            self.IsdacCC_dT = namelist['initial']['dSST'] + namelist['initial']['dTi'] - 5.0
            print('IsdacCC case: RRTM profiles are shifted according to %2.2f temperature change.'%(self.IsdacCC_dT))
        except:
            self.IsdacCC_dT = 0.0

        # self.bl = ReadProfiles.boundary_layer_profiles(namelist)

        #self.out_file = str(namelist['output']['file'])

        return

    cpdef initialize(self, ReadProfiles.ReadProfiles pf):

        self.heating_rate = np.zeros((pf.nz,), dtype=np.double)
        self.net_lw_flux = np.zeros((pf.nz,), dtype=np.double)

        # NS.add_profile('net_lw_flux')
        # NS.add_profile('radiative_heating_rate')

        return

    cpdef initialize_profiles(self, ReadProfiles.ReadProfiles pf):

        cdef:
            # Py_ssize_t qv_shift = DV.get_varshift(Gr, 'qv')
            # Py_ssize_t t_shift = DV.get_varshift(Gr, 'temperature')
            # double [:,:] qv_pencils =  self.z_pencil.forward_double(&Gr.dims, Pa, &DV.values[qv_shift])
            # double [:,:] t_pencils =  self.z_pencil.forward_double(&Gr.dims, Pa, &DV.values[t_shift])

            Py_ssize_t nz = pf.nz
            Py_ssize_t i,k


        # pf.get_profiles(mlm_vars)

        if not pf.model=='clima':
            self.t_surface = pf.t_surface
        # Construct the extension of the profiles, including a blending region between the given profile and LES domain (if desired)
        if self.read_file:
            rdr = cfreader(self.file, self.site)
            pressures = rdr.get_profile_mean('pfull')
            temperatures = rdr.get_profile_mean('ta')
            vapor_mixing_ratios = rdr.get_profile_mean('hus')
            if pf.model=='clima':
                self.t_surface = rdr.get_timeseries_mean('ts')
        else:
            pressures = profile_data[self.profile_name]['pressure'][:]
            temperatures = profile_data[self.profile_name]['temperature'][:]
            vapor_mixing_ratios = profile_data[self.profile_name]['vapor_mixing_ratio'][:]
            #specific_humidity = profile_data[self.profile_name]['specific_humidity'][:]


        dp = np.abs(pf.pressure[-1] - pf.pressure[-2])
        self.patch_pressure = np.minimum(self.patch_pressure, pf.pressure[-1] - dp)

        # n_profile = len(pressures[pressures<=self.patch_pressure]) # nprofile = # of points in the fixed profile to use
        n_profile = 0
        for pressure in pressures:
            if pressure <= self.patch_pressure:
                n_profile += 1
        # print(n_profile)
        self.n_ext =  n_profile + self.n_buffer # n_ext = total # of points to add to LES domain (buffer portion + fixed profile portion)

        # Create the space for the extensions (to be tacked on to top of LES pencils)
        # we declare these as class members in case we want to modify the buffer zone during run time
        # i.e. if there is some drift to top of LES profiles

        self.p_ext = np.zeros((self.n_ext,),dtype=np.double)
        self.t_ext = np.zeros((self.n_ext,),dtype=np.double)
        self.rv_ext = np.zeros((self.n_ext,),dtype=np.double)

        cdef Py_ssize_t count = 0
        for k in xrange(len(pressures)-n_profile, len(pressures)):
            self.p_ext[self.n_buffer+count] = pressures[k]
            #qt_new = specific_humidity[k]
            self.t_ext[self.n_buffer+count] = temperatures[k]
            #self.rv_ext[self.n_buffer+count] = qt_new / (1.0 - qt_new)
            self.rv_ext[self.n_buffer+count] = vapor_mixing_ratios[k]
            count += 1


        # Now  create the buffer zone
        if self.n_buffer > 0:
            #dp = np.abs(Ref.p0_half_global[nz + gw -1] - Ref.p0_half_global[nz + gw -2])
            dp = np.abs(pf.pressure[-1] - pf.pressure[-2])
            #self.p_ext[0] = Ref.p0_half_global[nz + gw -1] - dp
            self.p_ext[0] = pf.pressure[-1] - dp

            # print(self.p_ext[0])
            for i in range(1,self.n_buffer):
                self.p_ext[i] = self.p_ext[i-1] - (i+1.0)**self.stretch_factor * dp

            # for i in xrange(self.n_ext):
                # print i, self.p_ext[i]

            # Pressures of "data" points for interpolation, must be INCREASING pressure
            xi = np.array([self.p_ext[self.n_buffer+1],self.p_ext[self.n_buffer],pf.pressure[-1],pf.pressure[-2] ],dtype=np.double)
            print(xi)


            # interpolation for temperature
            ti = np.array([self.t_ext[self.n_buffer+1],self.t_ext[self.n_buffer], pf.temperature[-1], pf.temperature[-2] ], dtype = np.double)
            # interpolation for vapor mixing ratio
            # rv_m2 = qv_pencils[0, nz-2]/ (1.0 - qv_pencils[0, nz-2])
            # rv_m1 = qv_pencils[0,nz-1]/(1.0-qv_pencils[0,nz-1])

            rv_m2 = pf.qv[nz-2]/(1.0 - pf.qv[nz-2])
            rv_m1 = pf.qv[nz-1]/(1.0 - pf.qv[nz-1])

            ri = np.array([self.rv_ext[self.n_buffer+1],self.rv_ext[self.n_buffer], rv_m1, rv_m2 ], dtype = np.double)

            for i in xrange(self.n_buffer):
                self.rv_ext[i] = pchip_interpolate(xi, ri, self.p_ext[i] )
                self.t_ext[i] = pchip_interpolate(xi, ti, self.p_ext[i])

        # # Plotting to evaluate implementation of buffer zone
        # plt.figure(1)
        # plt.scatter(self.rv_ext,self.p_ext, label='Ext')
        # plt.plot(specific_humidity, pressures, label='Profile.py')
        # plt.plot(pf.qv, pf.pressure, label='LES')
        # plt.legend()
        # # plt.savefig('rrtm_buffer_rv.png')
        # plt.figure(2)
        # plt.scatter(self.t_ext,self.p_ext, label='Ext')
        # plt.plot(temperatures,pressures, label='Profile.py')
        # plt.plot(pf.temperature, pf.pressure, label='LES')
        # plt.legend()
        # # plt.savefig('rrtm_buffer_t.png')
        # plt.figure(10)
        # plt.plot(pf.ql, pf.pressure)
        # plt.show()

        self.p_full = np.zeros((self.n_ext+nz,), dtype=np.double)
        self.pi_full = np.zeros((self.n_ext+1+nz,),dtype=np.double)

        self.p_full[0:nz] = pf.pressure #Ref.p0_half_global[gw:nz+gw] # at cell center
        self.p_full[nz:]=self.p_ext[:]
        # print('p_ext ', np.array(self.p_ext))

        # self.pi_full[0:nz] = Ref.p0_global[gw:nz+gw] # at cell interface
        # self.pi_full[0:nz] = pf.pressure_i[1:]
        self.pi_full[0:nz] = pf.pressure_i[:-1]
        # print('pi_full before all filled:', np.array(self.pi_full))
        # print('pf.pressure_i:', np.array(pf.pressure_i))

        for i in range(nz,self.n_ext+nz):
            self.pi_full[i] = (self.p_full[i] + self.p_full[i-1]) * 0.5
        self.pi_full[self.n_ext +  nz] = 2.0 * self.p_full[self.n_ext + nz -1 ] - self.pi_full[self.n_ext + nz -1]

        # print('pi_full:', np.array(self.pi_full))
        # print('p_full:', np.array(self.p_full))


        # try to get ozone
        try:
            o3_trace = profile_data[self.profile_name]['o3_vmr'][:]   # O3 VMR (from SRF to TOP)
            o3_pressure = profile_data[self.profile_name]['pressure'][:]/100.0       # Pressure (from SRF to TOP) in hPa
            # can't do simple interpolation... Need to conserve column path !!!
            use_o3in = True
        except:
            try:
                o3_trace = profile_data[self.profile_name]['o3_mr'][:]*28.97/47.9982   # O3 MR converted to VMR
                o3_pressure = profile_data[self.profile_name]['pressure'][:]/100.0       # Pressure (from SRF to TOP) in hPa
                # can't do simple interpolation... Need to conserve column path !!!
                use_o3in = True

            except:
                print('O3 profile not set so default RRTM profile will be used.')
                use_o3in = False

        #Initialize rrtmg_lw and rrtmg_sw
        cdef double cpdair = np.float64(cpd)
        c_rrtmg_lw_init(&cpdair)
        c_rrtmg_sw_init(&cpdair)

        # Read in trace gas data
        lw_input_file = './RRTMG/lw/data/rrtmg_lw.nc'
        lw_gas = nc.Dataset(lw_input_file,  "r")

        lw_pressure = np.asarray(lw_gas.variables['Pressure'])
        lw_absorber = np.asarray(lw_gas.variables['AbsorberAmountMLS'])
        lw_absorber = np.where(lw_absorber>2.0, np.zeros_like(lw_absorber), lw_absorber)
        lw_ngas = lw_absorber.shape[1]
        lw_np = lw_absorber.shape[0]

        # 9 Gases: O3, CO2, CH4, N2O, O2, CFC11, CFC12, CFC22, CCL4
        # From rad_driver.f90, lines 546 to 552
        trace = np.zeros((9,lw_np),dtype=np.double,order='F')
        for i in xrange(lw_ngas):
            gas_name = ''.join(lw_gas.variables['AbsorberNames'][i,:])
            if 'O3' in gas_name:
                trace[0,:] = lw_absorber[:,i].reshape(1,lw_np)
            elif 'CO2' in gas_name:
                trace[1,:] = lw_absorber[:,i].reshape(1,lw_np)*self.co2_factor
            elif 'CH4' in gas_name:
                trace[2,:] = lw_absorber[:,i].reshape(1,lw_np)
            elif 'N2O' in gas_name:
                trace[3,:] = lw_absorber[:,i].reshape(1,lw_np)
            elif 'O2' in gas_name:
                trace[4,:] = lw_absorber[:,i].reshape(1,lw_np)
            elif 'CFC11' in gas_name:
                trace[5,:] = lw_absorber[:,i].reshape(1,lw_np)
            elif 'CFC12' in gas_name:
                trace[6,:] = lw_absorber[:,i].reshape(1,lw_np)
            elif 'CFC22' in gas_name:
                trace[7,:] = lw_absorber[:,i].reshape(1,lw_np)
            elif 'CCL4' in gas_name:
                trace[8,:] = lw_absorber[:,i].reshape(1,lw_np)

        # From rad_driver.f90, lines 585 to 620
        trpath = np.zeros((nz + self.n_ext + 1, 9),dtype=np.double,order='F')
        # plev = self.pi_full[:]/100.0
        for i in xrange(1, nz + self.n_ext + 1):
            trpath[i,:] = trpath[i-1,:]
            if (self.pi_full[i-1]/100.0 > lw_pressure[0]):
                trpath[i,:] = trpath[i,:] + (self.pi_full[i-1]/100.0 - np.max((self.pi_full[i]/100.0,lw_pressure[0])))/g*trace[:,0]
            for m in xrange(1,lw_np):
                #print i, m
                plow = np.min((self.pi_full[i-1]/100.0,np.max((self.pi_full[i]/100.0, lw_pressure[m-1]))))
                pupp = np.min((self.pi_full[i-1]/100.0,np.max((self.pi_full[i]/100.0, lw_pressure[m]))))
                if (plow > pupp):
                    pmid = 0.5*(plow+pupp)
                    wgtlow = (pmid-lw_pressure[m])/(lw_pressure[m-1]-lw_pressure[m])
                    wgtupp = (lw_pressure[m-1]-pmid)/(lw_pressure[m-1]-lw_pressure[m])
                    trpath[i,:] = trpath[i,:] + (plow-pupp)/g*(wgtlow*trace[:,m-1]  + wgtupp*trace[:,m])
            if (self.pi_full[i]/100.0 < lw_pressure[lw_np-1]):
                trpath[i,:] = trpath[i,:] + (np.min((self.pi_full[i-1]/100.0,lw_pressure[lw_np-1]))-self.pi_full[i]/100.0)/g*trace[:,lw_np-1]

        tmpTrace = np.zeros((nz + self.n_ext,9),dtype=np.double,order='F')
        for i in xrange(9):
            for k in xrange(nz + self.n_ext):
                tmpTrace[k,i] = g*100.0/(self.pi_full[k]-self.pi_full[k+1])*(trpath[k+1,i]-trpath[k,i])

        if use_o3in == False:
            self.o3vmr  = np.array(tmpTrace[:,0],dtype=np.double, order='F')
        else:
            # o3_trace, o3_pressure
            trpath_o3 = np.zeros(nz + self.n_ext+1, dtype=np.double, order='F')
            # plev = self.pi_full/100.0
            o3_np = o3_trace.shape[0]
            for i in xrange(1, nz + self.n_ext+1):
                trpath_o3[i] = trpath_o3[i-1]
                if (self.pi_full[i-1]/100.0 > o3_pressure[0]):
                    trpath_o3[i] = trpath_o3[i] + (self.pi_full[i-1]/100.0 - np.max((self.pi_full[i]/100.0,o3_pressure[0])))/g*o3_trace[0]
                for m in xrange(1,o3_np):
                    #print i, m
                    plow = np.min((self.pi_full[i-1]/100.0,np.max((self.pi_full[i]/100.0, o3_pressure[m-1]))))
                    pupp = np.min((self.pi_full[i-1]/100.0,np.max((self.pi_full[i]/100.0, o3_pressure[m]))))
                    if (plow > pupp):
                        pmid = 0.5*(plow+pupp)
                        wgtlow = (pmid-o3_pressure[m])/(o3_pressure[m-1]-o3_pressure[m])
                        wgtupp = (o3_pressure[m-1]-pmid)/(o3_pressure[m-1]-o3_pressure[m])
                        trpath_o3[i] = trpath_o3[i] + (plow-pupp)/g*(wgtlow*o3_trace[m-1]  + wgtupp*o3_trace[m])
                if (self.pi_full[i]/100.0 < o3_pressure[o3_np-1]):
                    trpath_o3[i] = trpath_o3[i] + (np.min((self.pi_full[i-1]/100.0,o3_pressure[o3_np-1]))-self.pi_full[i]/100.0)/g*o3_trace[o3_np-1]
            tmpTrace_o3 = np.zeros( nz + self.n_ext, dtype=np.double, order='F')
            for k in xrange(nz + self.n_ext):
                tmpTrace_o3[k] = g *100.0/(self.pi_full[k]-self.pi_full[k+1])*(trpath_o3[k+1]-trpath_o3[k])
            self.o3vmr = np.array(tmpTrace_o3[:],dtype=np.double, order='F')

        self.co2vmr = np.array(tmpTrace[:,1],dtype=np.double, order='F')
        self.ch4vmr =  np.array(tmpTrace[:,2],dtype=np.double, order='F')
        self.n2ovmr =  np.array(tmpTrace[:,3],dtype=np.double, order='F')
        self.o2vmr  =  np.array(tmpTrace[:,4],dtype=np.double, order='F')
        self.cfc11vmr =  np.array(tmpTrace[:,5],dtype=np.double, order='F')
        self.cfc12vmr =  np.array(tmpTrace[:,6],dtype=np.double, order='F')
        self.cfc22vmr = np.array( tmpTrace[:,7],dtype=np.double, order='F')
        self.ccl4vmr  =  np.array(tmpTrace[:,8],dtype=np.double, order='F')


        #Initialize NetCDF4 for RRTM outputs
        root_grp = nc.Dataset(pf.out_file, 'w', format='NETCDF4')
        root_grp.createDimension('z', len(self.pi_full))
        root_grp.createDimension('t', None)
        pi_full_ref = root_grp.createVariable('pi_full_ref', 'f8', ('z'))
        pi_full_ref[:] = np.array(self.pi_full)
        root_grp.createVariable('t', 'f8', ('t'))

        root_grp.createVariable('pi_full', 'f8', ('t', 'z'))
        root_grp.createVariable('uflx_lw', 'f8', ('t', 'z'))
        root_grp.createVariable('uflxc_lw', 'f8', ('t', 'z'))
        root_grp.createVariable('dflx_lw', 'f8', ('t', 'z'))
        root_grp.createVariable('dflxc_lw', 'f8', ('t', 'z'))
        root_grp.createVariable('uflx_sw', 'f8', ('t', 'z'))
        root_grp.createVariable('uflxc_sw', 'f8', ('t', 'z'))
        root_grp.createVariable('dflx_sw', 'f8', ('t', 'z'))
        root_grp.createVariable('dflxc_sw', 'f8', ('t', 'z'))

        root_grp.close()

        return

    cpdef update(self, ReadProfiles.ReadProfiles pf):
        # if TS.rk_step == 0:
        #     if self.radiation_frequency <= 0.0:
        #         self.update_RRTM(pf)
        #     elif TS.t >= self.next_radiation_calculate:
        #         self.update_RRTM(pf)
        #         self.next_radiation_calculate = (TS.t//self.radiation_frequency + 1.0) * self.radiation_frequency

        self.update_RRTM(pf)

        return

    cdef update_RRTM(self, ReadProfiles.ReadProfiles pf):

        #self.coszen = pf.toa_sw/self.scon
        #self.adir = pf.albedo

        cdef:
            Py_ssize_t nz = pf.nz
            Py_ssize_t nz_full = self.n_ext + nz
            Py_ssize_t k
            Py_ssize_t n_pencils = 1
            Py_ssize_t ngptlw = 140
            Py_ssize_t ngptsw = 112

        cdef:
            double [:] rl_full = np.zeros((nz_full,), dtype=np.double, order='F')
            double [:] ri_full = np.zeros((nz_full,), dtype=np.double, order='F')
            double [:] play_in = np.zeros((nz_full,), dtype=np.double, order='F')
            double [:] plev_in = np.zeros((nz_full + 1), dtype=np.double, order='F')
            double [:] tlay_in = np.zeros((nz_full,), dtype=np.double, order='F')
            double [:] tlev_in = np.zeros((nz_full + 1), dtype=np.double, order='F')
            double [:] tsfc_in = np.ones((n_pencils),dtype=np.double,order='F') * self.t_surface
            double [:] h2ovmr_in = np.zeros((nz_full,),dtype=np.double,order='F')
            double [:] o3vmr_in  = np.zeros((nz_full,),dtype=np.double,order='F')
            double [:] co2vmr_in = np.zeros((nz_full,),dtype=np.double,order='F')
            double [:] ch4vmr_in = np.zeros((nz_full,),dtype=np.double,order='F')
            double [:] n2ovmr_in = np.zeros((nz_full,),dtype=np.double,order='F')
            double [:] o2vmr_in  = np.zeros((nz_full,),dtype=np.double,order='F')
            double [:] cfc11vmr_in = np.zeros((nz_full,),dtype=np.double,order='F')
            double [:] cfc12vmr_in = np.zeros((nz_full,),dtype=np.double,order='F')
            double [:] cfc22vmr_in = np.zeros((nz_full,),dtype=np.double,order='F')
            double [:] ccl4vmr_in = np.zeros((nz_full,),dtype=np.double,order='F')
            double [:,:] emis_in = np.ones((n_pencils,16),dtype=np.double,order='F') * 0.95
            double [:] cldfr_raw  = np.zeros((nz_full,),dtype=np.double,order='F')
            double [:,:] cldfr_lw_in  = np.zeros((ngptlw,nz_full,),dtype=np.double,order='F')
            double [:,:] cldfr_sw_in  = np.zeros((ngptsw,nz_full,),dtype=np.double,order='F')
            double [:] cicewp_raw = np.zeros((nz_full,),dtype=np.double,order='F')
            double [:,:] cicewp_lw_in = np.zeros((ngptlw,nz_full,),dtype=np.double,order='F')
            double [:,:] cicewp_sw_in = np.zeros((ngptsw,nz_full,),dtype=np.double,order='F')
            double [:] cliqwp_raw = np.zeros((nz_full,),dtype=np.double,order='F')
            double [:,:] cliqwp_lw_in = np.zeros((ngptlw,nz_full,),dtype=np.double,order='F')
            double [:,:] cliqwp_sw_in = np.zeros((ngptsw,nz_full,),dtype=np.double,order='F')
            double [:] reice_raw  = np.zeros((nz_full,),dtype=np.double,order='F')
            double [:] reice_lw_in  = np.zeros((nz_full,),dtype=np.double,order='F')
            double [:] reice_sw_in  = np.zeros((nz_full,),dtype=np.double,order='F')
            double [:] reliq_raw  = np.zeros((nz_full,),dtype=np.double,order='F')
            double [:] reliq_lw_in  = np.zeros((nz_full,),dtype=np.double,order='F')
            double [:] reliq_sw_in  = np.zeros((nz_full,),dtype=np.double,order='F')
            double [:] coszen_in = np.ones((n_pencils),dtype=np.double,order='F') *self.coszen
            double [:] asdir_in = np.ones((n_pencils),dtype=np.double,order='F') * self.adir
            double [:] asdif_in = np.ones((n_pencils),dtype=np.double,order='F') * self.adif
            double [:] aldir_in = np.ones((n_pencils),dtype=np.double,order='F') * self.adir
            double [:] aldif_in = np.ones((n_pencils),dtype=np.double,order='F') * self.adif
            double [:,:] taucld_lw_raw  = np.zeros((ngptlw,nz_full,),dtype=np.double,order='F')
            double [:,:] taucld_lw_in  = np.zeros((ngptlw,nz_full,),dtype=np.double,order='F')
            double [:,:] tauaer_lw_in  = np.zeros((nz_full,16),dtype=np.double,order='F')
            double [:,:] taucld_sw_raw  = np.zeros((ngptsw,nz_full,),dtype=np.double,order='F')
            double [:,:] taucld_sw_in  = np.zeros((ngptsw,nz_full,),dtype=np.double,order='F')
            double [:,:] ssacld_sw_raw  = np.zeros((ngptsw,nz_full,),dtype=np.double,order='F')
            double [:,:] ssacld_sw_in  = np.zeros((ngptsw,nz_full,),dtype=np.double,order='F')
            double [:,:] asmcld_sw_raw  = np.zeros((ngptsw,nz_full,),dtype=np.double,order='F')
            double [:,:] asmcld_sw_in  = np.zeros((ngptsw,nz_full,),dtype=np.double,order='F')
            double [:,:] fsfcld_sw_raw  = np.zeros((ngptsw,nz_full,),dtype=np.double,order='F')
            double [:,:] fsfcld_sw_in  = np.zeros((ngptsw,nz_full,),dtype=np.double,order='F')
            double [:,:] tauaer_sw_in  = np.zeros((nz_full,14),dtype=np.double,order='F')
            double [:,:] ssaaer_sw_in  = np.zeros((nz_full,14),dtype=np.double,order='F')
            double [:,:] asmaer_sw_in  = np.zeros((nz_full,14),dtype=np.double,order='F')
            double [:,:] ecaer_sw_in  = np.zeros((nz_full,6),dtype=np.double,order='F')

            # Output
            double[:] uflx_lw_out = np.zeros((nz_full +1),dtype=np.double,order='F')
            double[:] dflx_lw_out = np.zeros((nz_full +1),dtype=np.double,order='F')
            double[:] hr_lw_out = np.zeros((nz_full,),dtype=np.double,order='F')
            double[:] uflxc_lw_out = np.zeros((nz_full + 1),dtype=np.double,order='F')
            double[:] dflxc_lw_out = np.zeros((nz_full +1),dtype=np.double,order='F')
            double[:] hrc_lw_out = np.zeros((nz_full,),dtype=np.double,order='F')
            double[:] duflx_dt_out = np.zeros((nz_full +1),dtype=np.double,order='F')
            double[:] duflxc_dt_out = np.zeros((nz_full +1),dtype=np.double,order='F')
            double[:] uflx_sw_out = np.zeros((nz_full +1),dtype=np.double,order='F')
            double[:] dflx_sw_out = np.zeros((nz_full +1),dtype=np.double,order='F')
            double[:] hr_sw_out = np.zeros((nz_full,),dtype=np.double,order='F')
            double[:] uflxc_sw_out = np.zeros((nz_full +1),dtype=np.double,order='F')
            double[:] dflxc_sw_out = np.zeros((nz_full +1),dtype=np.double,order='F')
            double[:] hrc_sw_out = np.zeros((nz_full,),dtype=np.double,order='F')

            double rv_to_reff = np.exp(np.log(1.2)**2.0)*10.0*1000.0

        self.p_full[0:nz] = pf.pressure
        self.pi_full[0:nz] = pf.pressure_i[:-1]
        # with nogil:
        for k in xrange(nz, nz_full):
            tlay_in[k] = self.t_ext[k-nz]
            h2ovmr_in[k] = self.rv_ext[k-nz] * Rv/Rd * self.h2o_factor
                # Assuming for now that there is no condensate above LES domain!
        for k in xrange(nz):
            tlay_in[k] = pf.temperature[k]
            h2ovmr_in[k] = pf.qv[k]/ (1.0 - pf.qv[k])* Rv/Rd * self.h2o_factor
            rl_full[k] = (pf.ql[k])/ (1.0 - pf.qv[k])
            ri_full[k] = (pf.qi[k])/ (1.0 - pf.qv[k])
            #ri_full[k] = 0.0
            cliqwp_raw[k] = ((pf.ql[k])/ (1.0 - pf.qv[k])
                               *1.0e3*(self.pi_full[k] - self.pi_full[k+1])/g)
            cicewp_raw[k] = ((pf.qi[k])/ (1.0 - pf.qv[k])
                               *1.0e3*(self.pi_full[k] - self.pi_full[k+1])/g)
            # if pf.ql[k] + pf.qi[k] > ql_threshold:
            #if pf.ql[k] > ql_threshold:
            #    cldfr_in[k] = 1.0
            cldfr_raw[k] = pf.cf[k]

        with nogil:
            for k in xrange(nz_full):
                play_in[k] = self.p_full[k]/100.0
                o3vmr_in[k] = self.o3vmr[k]
                co2vmr_in[k] = self.co2vmr[k]
                ch4vmr_in[k] = self.ch4vmr[k]
                n2ovmr_in[k] = self.n2ovmr[k]
                o2vmr_in [k] = self.o2vmr[k]
                cfc11vmr_in[k] = self.cfc11vmr[k]
                cfc12vmr_in[k] = self.cfc12vmr[k]
                cfc22vmr_in[k] = self.cfc22vmr[k]
                ccl4vmr_in[k] = self.ccl4vmr[k]


                if self.uniform_reliq:
                    reliq_raw[k] = 14.0*pf.cf[k]
                else:
                    reliq_raw[k] = ((3.0*self.p_full[k]/Rd/tlay_in[k]*rl_full[k]/
                                        fmax(pf.cf[k],1.0e-6))/(4.0*pi*1.0e3*100.0))**(1.0/3.0)
                    reliq_raw[k] = fmin(fmax(reliq_raw[ k]*rv_to_reff, 2.5), 60.0)

                # Boudala et al. (2002) Eqn 10a, this is dge (generalized effective size),
                # and is what iceflglw=3 calls for. Will only work with iceflglw=iceflgsw=3!
                reice_raw[k] = 53.005 * ((self.p_full[k]/Rd/tlay_in[k]*ri_full[k]*1.0e3)/
                                            fmax(pf.cf[k],1.0e-6)) ** 0.06 \
                                      * exp(0.013*(tlay_in[k] - 273.16))
                reice_raw[k] = fmin(fmax(reice_raw[k], 5.0), 140.0) # Threshold from rrtmg sw instruction


            with gil:
                tlev_in[0] = self.t_surface
            plev_in[0] = self.pi_full[0]/100.0
            for k in xrange(1,nz_full):
                tlev_in[k] = 0.5*(tlay_in[k-1]+tlay_in[k])
                plev_in[k] = self.pi_full[k]/100.0
            tlev_in[nz_full] = 2.0*tlay_in[nz_full-1] - tlev_in[nz_full-1]
            plev_in[nz_full] = self.pi_full[nz_full]/100.0

        #========================================

        # # Test whether RRTM works
        # # Construct the extension of the profiles, including a blending region between the given profile and LES domain (if desired)
        # pressures = profile_data[self.profile_name]['pressure'][::-1]
        # temperatures = profile_data[self.profile_name]['temperature'][::-1]
        # vapor_mixing_ratios = profile_data[self.profile_name]['vapor_mixing_ratio'][:]
        # specific_humidity = profile_data[self.profile_name]['specific_humidity'][::-1]
        #
        # #Interpolate onto input grid
        # tlay_in = pchip_interpolate(pressures, temperatures, play_in)
        # h2ovmr_in = pchip_interpolate(pressures, vapor_mixing_ratios, play_in) * Rv/Rd * self.h2o_factor
        #
        # for k in xrange(nz_full):
        #     tlev_in[0] = pf.t_surface
        #     plev_in[0] = self.pi_full[0]/100.0
        # for k in xrange(1,nz_full):
        #     tlev_in[k] = 0.5*(tlay_in[k-1]+tlay_in[k])
        #     plev_in[k] = self.pi_full[k]/100.0
        # tlev_in[nz_full] = 2.0*tlay_in[nz_full-1] - tlev_in[nz_full-1]
        # plev_in[nz_full] = self.pi_full[nz_full]/100.0

        #========================================

        # # Plot the variables to check
        # plt.figure(3)
        # plt.subplot(121)
        # plt.plot(tlay_in, play_in, label='tlay')
        # plt.plot(tlev_in, plev_in, label='tlev')
        # plt.legend()
        # plt.xlabel('tlay_in tlev_in')
        # plt.subplot(122)
        # plt.plot(h2ovmr_in, self.p_full)
        # plt.xlabel('h2ovmr_in')
        # plt.figure(4)
        # plt.subplot(121)
        # plt.plot(rl_full, self.p_full)
        # plt.xlabel('rl_full')
        # plt.subplot(122)
        # plt.plot(cliqwp_in, self.p_full)
        # plt.plot(cicewp_in, self.p_full)
        # plt.xlabel('cliqwp_in & cicewp_in')
        # plt.figure(5)
        # plt.subplot(121)
        # plt.plot(cldfr_in, self.p_full)
        # plt.xlabel('cldfr_in')
        # plt.subplot(122)
        # plt.plot(reliq_in, self.p_full)
        # plt.xlabel('reliq_in')
        # plt.figure(6)
        # plt.subplot(122)
        # plt.plot(reice_in, self.p_full)
        # plt.xlabel('reice_in')
        # plt.show()

        cdef:
            int iplon = 1
            int ncol = n_pencils
            int nlay = nz_full
            int icld = 1
            int irng = 1
            int seedlw = 1
            int seedsw = 10000
            int idrv = 0
            int iaer = 0
            int inflglw = 2
            int iceflglw = 3
            int liqflglw = 1
            int inflgsw = 2
            int iceflgsw = 3
            int liqflgsw = 1

        # print('Begin RRTM calculations!')
        c_mcica_subcol_lw (&iplon, &ncol, &nlay, &icld, &seedlw, &irng, 
           &play_in[0], &cldfr_raw[0], &cicewp_raw[0], &cliqwp_raw[0], &reice_raw[0], &reliq_raw[0], &taucld_lw_raw[0,0],
           &cldfr_lw_in[0,0], &cicewp_lw_in[0,0], &cliqwp_lw_in[0,0], &reice_lw_in[0], &reliq_lw_in[0], &taucld_lw_in[0,0])

        c_rrtmg_lw (
             &ncol    ,&nlay    ,&icld    ,&idrv,
             &play_in[0]    ,&plev_in[0]    ,&tlay_in[0]    ,&tlev_in[0]    ,&tsfc_in[0]    ,
             &h2ovmr_in[0]  ,&o3vmr_in[0]   ,&co2vmr_in[0]  ,&ch4vmr_in[0]  ,&n2ovmr_in[0]  ,&o2vmr_in[0],
             &cfc11vmr_in[0],&cfc12vmr_in[0],&cfc22vmr_in[0],&ccl4vmr_in[0] ,&emis_in[0,0]    ,
             &inflglw ,&iceflglw,&liqflglw,&cldfr_lw_in[0,0]   ,
             &taucld_lw_in[0,0]  ,&cicewp_lw_in[0,0]  ,&cliqwp_lw_in[0,0]  ,&reice_lw_in[0]   ,&reliq_lw_in[0]   ,
             &tauaer_lw_in[0,0]  ,
             &uflx_lw_out[0]    ,&dflx_lw_out[0]    ,&hr_lw_out[0]      ,&uflxc_lw_out[0]   ,&dflxc_lw_out[0],  &hrc_lw_out[0],
             &duflx_dt_out[0],&duflxc_dt_out[0] )
        # print('Done RRTM LW!')
        c_mcica_subcol_sw (&iplon, &ncol, &nlay, &icld, &seedsw, &irng, 
           &play_in[0], &cldfr_raw[0], &cicewp_raw[0], &cliqwp_raw[0], &reice_raw[0], &reliq_raw[0], 
           &taucld_sw_raw[0,0], &ssacld_sw_raw[0,0], &asmcld_sw_raw[0,0], &fsfcld_sw_raw[0,0],
           &cldfr_sw_in[0,0], &cicewp_sw_in[0,0], &cliqwp_sw_in[0,0], &reice_sw_in[0], &reliq_sw_in[0],
           &taucld_sw_in[0,0], &ssacld_sw_in[0,0], &asmcld_sw_in[0,0], &fsfcld_sw_in[0,0])
        c_rrtmg_sw (
            &ncol, &nlay, &icld, &iaer, &play_in[0], &plev_in[0], &tlay_in[0], &tlev_in[0],&tsfc_in[0],
            &h2ovmr_in[0], &o3vmr_in[0], &co2vmr_in[0], &ch4vmr_in[0], &n2ovmr_in[0],&o2vmr_in[0],
             &asdir_in[0]   ,&asdif_in[0]   ,&aldir_in[0]   ,&aldif_in[0]   ,
             &coszen_in[0]  ,&self.adjes   ,&self.dyofyr  ,&self.scon   ,
             &inflgsw ,&iceflgsw,&liqflgsw,&cldfr_sw_in[0,0]   ,
             &taucld_sw_in[0,0]  ,&ssacld_sw_in[0,0]  ,&asmcld_sw_in[0,0]  ,&fsfcld_sw_in[0,0]  ,
             &cicewp_sw_in[0,0]  ,&cliqwp_sw_in[0,0]  ,&reice_sw_in[0]   ,&reliq_sw_in[0]   ,
             &tauaer_sw_in[0,0]  ,&ssaaer_sw_in[0,0]  ,&asmaer_sw_in[0,0]  ,&ecaer_sw_in[0,0]   ,
             &uflx_sw_out[0]    ,&dflx_sw_out[0]    ,&hr_sw_out[0]      ,&uflxc_sw_out[0]   ,&dflxc_sw_out[0], &hrc_sw_out[0])

        # print('Done RRTM SW!')
        # for k in xrange(nz):
        #     self.heating_rate[k] = (hr_lw_out[k] + hr_sw_out[k]) * pf.rho[k] * cpd/86400.0
        #     self.net_lw_flux[k] = uflx_lw_out[k] - dflx_lw_out[k]

        # print(np.array(uflx_lw_out))
        # print(np.array(dflx_lw_out))
        # print(np.array(uflx_sw_out))

        # print(dflx_sw_out[-1])

        # print(np.array(hr_lw_out))
        # print(np.array(hr_sw_out))

        # plt.figure(6)
        # plt.subplot(121)
        # plt.plot(uflx_lw_out, self.pi_full, 'b')
        # plt.plot(uflxc_lw_out, self.pi_full, 'b--')
        # plt.plot(dflx_lw_out, self.pi_full, 'r')
        # plt.plot(dflxc_lw_out, self.pi_full, 'r--')
        # plt.xlabel('lw_out')
        # plt.subplot(122)
        # plt.plot(uflx_sw_out, self.pi_full, 'b')
        # plt.plot(dflx_sw_out, self.pi_full, 'r')
        # plt.plot(uflxc_sw_out, self.pi_full, 'b--')
        # plt.plot(dflxc_sw_out, self.pi_full, 'r--')
        # plt.xlabel('sw_out')
        # plt.figure(7)
        # plt.subplot(121)
        # # plt.plot(self.heating_rate, pf.pressure)
        # plt.plot(hr_lw_out, self.p_full, 'r', label='LW')
        # plt.plot(hrc_lw_out, self.p_full, 'r--',label='LW')
        # plt.plot(hr_sw_out, self.p_full, label='SW')
        # plt.xlabel('heating rates')
        # plt.legend()
        # # plt.subplot(122)
        # # plt.plot(self.net_lw_flux, pf.pressure)
        # # plt.xlabel('net lw flux')
        # plt.show()

        #==============================================
        #Save the RRTM outputs

        # data_out = {}
        # data_out['p_full'] = np.array(self.p_full)
        # data_out['pi_full'] = np.array(self.pi_full)
        #
        # #LW fluxes
        # data_out['uflx_lw'] = np.array(uflx_lw_out)
        # data_out['uflxc_lw'] = np.array(uflxc_lw_out)
        # data_out['dflx_lw'] = np.array(dflx_lw_out)
        # data_out['dflxc_lw'] = np.array(dflxc_lw_out)
        #
        # #SW fluxes
        # data_out['uflx_sw'] = np.array(uflx_sw_out)
        # data_out['uflxc_sw'] = np.array(uflxc_sw_out)
        # data_out['dflx_sw'] = np.array(dflx_sw_out)
        # data_out['dflxc_sw'] = np.array(dflxc_sw_out)
        #
        # #Heating rates (K/day)
        # data_out['hr_lw'] = np.array(hr_lw_out)
        # data_out['hrc_lw'] = np.array(hrc_lw_out)
        # data_out['hr_sw'] = np.array(hr_sw_out)
        # data_out['hrc_sw'] = np.array(hrc_sw_out)
        #
        #
        # output_path = pf.path+self.out_file
        # fh = open(output_path, 'wb')
        # pkl.dump(data_out, fh)
        # fh.close()

        # NetCDF4 format
        root_grp = nc.Dataset(pf.out_file, 'r+', format='NETCDF4')
        #First write time
        rrtm_t = root_grp.variables['t']
        rrtm_t[rrtm_t.shape[0]] = pf.count

        var = None
        var = root_grp.variables['pi_full']
        var[-1, :] = np.array(self.pi_full)

        var = None
        var = root_grp.variables['uflx_lw']
        var[-1, :] = np.array(uflx_lw_out)

        var = None
        var = root_grp.variables['uflxc_lw']
        var[-1, :] = np.array(uflxc_lw_out)

        var = None
        var = root_grp.variables['dflx_lw']
        var[-1, :] = np.array(dflx_lw_out)

        var = None
        var = root_grp.variables['dflxc_lw']
        var[-1, :] = np.array(dflxc_lw_out)

        var = None
        var = root_grp.variables['uflx_sw']
        var[-1, :] = np.array(uflx_sw_out)

        var = None
        var = root_grp.variables['uflxc_sw']
        var[-1, :] = np.array(uflxc_sw_out)

        var = None
        var = root_grp.variables['dflx_sw']
        var[-1, :] = np.array(dflx_sw_out)

        var = None
        var = root_grp.variables['dflxc_sw']
        var[-1, :] = np.array(dflxc_sw_out)

        root_grp.close()

        # print("Finished saving RRTM output!")



        return

    # cpdef stats_io(self, NetCDFIO_Stats NS):
    #
    #     NS.write_profile('net_lw_flux', self.net_lw_flux)
    #     NS.write_profile('radiative_heating_rate', self.heating_rate)
    #
    #     return

# def get_humidity(temperature_old, qt_old, pressure, temperature_new):
#     pv_star_1 = saturation_vapor_pressure(temperature_old)
#     pv_1 = (pressure * qt_old) / (eps_v * (1.0 - qt_old) + qt_old)
#     rh_ = pv_1 / pv_star_1
#     pv_star_2 = saturation_vapor_pressure(temperature_new)
#     pv_2 = rh_ * pv_star_2
#     qt_new = 1.0/(eps_vi * (pressure - pv_2)/pv_2 + 1.0)
#     return qt_new
