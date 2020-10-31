#!python
# cython: boundscheck=False
# cython: wraparound=True
# cython: initializedcheck=False
# cython: cdivision=True

import numpy as np
#import matplotlib.pyplot as plt

# from mlm_thermodynamic_functions import *
# cimport mlm_thermodynamic_functions as mfun
# from mlm_thermodynamic_functions cimport *
cimport Radiation
# import Radiation
include 'parameters.pxi'
from variables import variable_name
# cimport TimeStepping
# from NetCDFIO cimport NetCDFIO_Stats
from libc.math cimport fmin, fabs
from scipy.integrate import odeint
import netCDF4 as nc

# cdef extern from "mlm_thermodynamic_functions.h":

cdef class ReadProfiles:
    def __init__(self, namelist):

        self.path_plus_file = str(namelist['input']['path'])+str(namelist['input']['file'])
        print(self.path_plus_file)

        self.path_plus_file_ref = self.path_plus_file

        try:
            self.path_plus_file_control = str(namelist['input']['control_path'])+str(namelist['input']['control_file'])
        except:
            self.path_plus_file_control = self.path_plus_file
        
        if self.path_plus_file_control == self.path_plus_file:
            self.out_file = str(namelist['input']['path'])+'RRTM_'+str(namelist['input']['file'])
        else:
            self.out_file = str(namelist['input']['path'])+'RRTM_control_'+str(namelist['input']['file'])
        #self.path = str(namelist['input']['path'])
        #self.albedo = namelist['input']['albedo']
        #self.toa_sw = namelist['input']['toa_sw']

        self.t1 = namelist['input']['t1']
        self.t2 = namelist['input']['t2']

        # self.root_grp = None
        self.profile_grp = None
        self.ref_grp = None
        self.ts_grp = None
        self.profile_grp2 = None
        self.ref_grp2 = None
        self.ts_grp2 = None

        self.count = 0
        self.average = namelist['input']['time_average']
        self.model = str(namelist['input']['model'])

        self.fix_T = namelist['input']['fix_T']
        self.fix_qv = namelist['input']['fix_qv']
        self.fix_cloud = namelist['input']['fix_cloud']
        self.fix_albedo = namelist['input']['fix_albedo']
        
        if self.average:
            self.out_file = str(namelist['input']['path'])+'RRTM_avg_'+str(namelist['input']['file'])

        if self.fix_T:
            self.out_file = str(namelist['input']['path'])+'RRTM_control_fixT_'+str(namelist['input']['file'])
        if self.fix_qv:
            self.out_file = str(namelist['input']['path'])+'RRTM_control_fixqv_'+str(namelist['input']['file'])
        if self.fix_T and self.fix_qv:
            self.out_file = str(namelist['input']['path'])+'RRTM_control_fixTqv_'+str(namelist['input']['file'])


    cpdef initialize(self):

        #Original profiles
        self.root_grp = nc.Dataset(self.path_plus_file, 'r')
        self.root_grp_ref = nc.Dataset(self.path_plus_file_ref, 'r')
        if self.model=='clima':
            self.profile_grp = self.root_grp
            self.ref_grp = self.root_grp
            self.ts_grp = self.root_grp
        else:
            self.profile_grp = self.root_grp.groups['profiles']
            self.ref_grp = self.root_grp_ref.groups['reference']
            self.ts_grp = self.root_grp.groups['timeseries']
        # self.root_grp.close()

        #Fixed profiles for PRP
        self.control_grp = nc.Dataset(self.path_plus_file_control, 'r')
        if self.model=='clima':
            self.profile_grp2 = self.control_grp
            self.ts_grp2 = self.control_grp
        else:
            self.profile_grp2 = self.control_grp.groups['profiles']
            #self.ref_grp2 = self.control_grp.groups['reference']
            self.ts_grp2 = self.control_grp.groups['timeseries']
        
        # #Albedo
        # self.albedo_ts = self.albedo_grp.groups['timeseries']['surface_albedo'][:]
        # # self.albedo_grp.close()
        
        if not self.model=='clima':
            self.pressure = self.ref_grp['p0'][:]
            self.rho = self.ref_grp['rho0'][:]

        if self.model=='clima':
            self.nz = len(self.profile_grp['z'][:])
            self.ntime = len(self.ts_grp['time'][:])
        else:
            self.nz = len(self.rho)
            self.ntime = len(self.ts_grp['t'][:])
        
        if not self.model=='clima':
            self.pressure_i = np.zeros(self.nz+1, dtype=np.double)
        
            for i in xrange(self.nz-1):
                self.pressure_i[i+1] = (self.pressure[i] + self.pressure[i+1])*0.5
            self.pressure_i[0] = self.pressure[0]*2 - self.pressure_i[1]#(self.pressure[0] - self.pressure[1]) + self.pressure[0]
            self.pressure_i[-1] = 2.0 * self.pressure[-1] - self.pressure_i[-2]

        return


    cpdef update(self,  Radiation.Radiation Ra):


        temp_name = variable_name['temperature'][self.model]
        qv_name = variable_name['qv'][self.model]
        ql_name = variable_name['ql'][self.model]
        qi_name = variable_name['qi'][self.model]
        cf_name = variable_name['cf'][self.model]
        pres_name = variable_name['pressure'][self.model]

        if self.average:
            self.temperature = np.mean(self.profile_grp[temp_name][self.t1:self.t2, :], axis=0)
            self.qv = np.mean(self.profile_grp[qv_name][self.t1:self.t2, :], axis=0)
            self.ql = np.mean(self.profile_grp[ql_name][self.t1:self.t2, :], axis=0)
            # TODO: add snow
            self.qi = np.mean(self.profile_grp[qi_name][self.t1:self.t2, :], axis=0)
            self.cf = np.mean(self.profile_grp[cf_name][self.t1:self.t2, :], axis=0)
            self.ntime = 1
            
            if self.model=='clima':
                self.pressure = np.mean(self.profile_grp[pres_name][self.t1:self.t2, :], axis=0)
            else:
                self.t_surface = np.mean(self.ts_grp['surface_temperature'][self.t1:self.t2])

            #self.toa_sw = np.mean(self.ts_grp['toa_sw_flux'][self.t1:self.t2])
            # print(self.toa_sw)

            #self.albedo = np.mean(self.ts_grp['surface_albedo'][self.t1:self.t2])
            # print(self.albedo)

        else:
            if self.fix_T:
                self.temperature = self.profile_grp2[temp_name][self.count, :]
                if self.model=='clima':
                    self.pressure = self.profile_grp2[pres_name][self.count, :]
                else:
                    self.t_surface = self.ts_grp2['surface_temperature'][self.count]
            else:
                self.temperature = self.profile_grp[temp_name][self.count, :]
                if self.model=='clima':
                    self.pressure = self.profile_grp[pres_name][self.count, :]
                else:
                    self.t_surface = self.ts_grp['surface_temperature'][self.count]

            if self.fix_qv:
                self.qv = self.profile_grp2[qv_name][self.count, :]
            else:
                self.qv = self.profile_grp[qv_name][self.count, :]

            #if self.fix_albedo:
            #    self.albedo = self.ts_grp2['surface_albedo'][self.count]
            #else:
            #    self.albedo = self.ts_grp['surface_albedo'][self.count]

            self.ql = self.profile_grp[ql_name][self.count, :]
            # TODO: add snow
            self.qi = self.profile_grp[qi_name][self.count, :]
            self.cf = self.profile_grp[cf_name][self.count, :]

            #self.toa_sw = self.ts_grp['toa_sw_flux'][self.count]
            # print(self.toa_sw)


            # self.albedo = self.albedo_ts[self.count]
            # print(self.albedo)
        if self.model=='clima':
            self.pressure_i = np.zeros(self.nz+1, dtype=np.double)
            for i in xrange(self.nz-1):
                self.pressure_i[i+1] = (self.pressure[i] + self.pressure[i+1])*0.5
            self.pressure_i[0] = self.pressure[0]*2 - self.pressure_i[1]#(self.pressure[0] - self.pressure[1]) + self.pressure[0]
            self.pressure_i[-1] = 2.0 * self.pressure[-1] - self.pressure_i[-2]

        self.count += 1
        print(self.count)

        # self.root_grp.close()
        
        return

    # cpdef stats_io(self, NetCDFIO_Stats NS):
    #
    #     NS.write_ts('zi', self.values[0])
    #     NS.write_ts('thetal_ml', self.values[1])
    #     NS.write_ts('qt_ml', self.values[2])
    #
    #     NS.write_profile('thetal', self.thetal)
    #     NS.write_profile('qt', self.qt)
    #     NS.write_profile('ql', self.ql)
    #     NS.write_profile('qi', self.qi)
    #     NS.write_profile('temperature', self.temperature)
    #     NS.write_profile('rho', self.rho)
    #     NS.write_profile('pressure', self.pressure)
    #
    #     cdef:
    #         Py_ssize_t kmin = 0
    #         Py_ssize_t kmax = self.nz
    #         Py_ssize_t k
    #         double cb
    #         double lwp
    #         double iwp
    #
    #     # Compute cloud bottom height
    #     cb = 99999.9
    #     with nogil:
    #         for k in xrange(kmin, kmax):
    #             if self.ql[k] > 0.0:
    #                 cb = fmin(cb, self.z[k])
    #
    #     NS.write_ts('cloud_base', cb)
    #
    #     # Compute liquid water path
    #     with nogil:
    #         for k in xrange(kmin, kmax):
    #             lwp += self.rho[k] * self.ql[k] * self.dz
    #
    #     NS.write_ts('lwp', lwp)
    #
    #     # Compute ice water path
    #     with nogil:
    #         for k in xrange(kmin, kmax):
    #             iwp += self.rho[k] * self.qi[k] * self.dz
    #
    #     NS.write_ts('iwp', iwp)
    #
    #     return
