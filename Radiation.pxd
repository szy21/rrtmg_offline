cimport ReadProfiles
# cimport TimeStepping
# from NetCDFIO cimport NetCDFIO_Stats


cdef class Radiation:
    cdef:
        double srf_lw_down
        double srf_lw_up
        double srf_sw_down
        double srf_sw_up
        str profile_name
        bint read_file
        str file
        int site
        Py_ssize_t n_buffer
        double stretch_factor
        double patch_pressure
        double t_surface
        double co2_factor
        double h2o_factor
        int dyofyr
        double scon
        double adjes
        double toa_sw
        double coszen
        double adif
        double adir
        bint uniform_reliq
        double radiation_frequency
        double next_radiation_calculate
        double [:] heating_rate
        public double [:] net_lw_flux
        Py_ssize_t n_ext
        double [:] p_ext
        double [:] t_ext
        double [:] rv_ext
        double [:] p_full
        double [:] pi_full
        double [:] o3vmr
        double [:] co2vmr
        double [:] ch4vmr
        double [:] n2ovmr
        double [:] o2vmr
        double [:] cfc11vmr
        double [:] cfc12vmr
        double [:] cfc22vmr
        double [:] ccl4vmr
        double IsdacCC_dT
        str out_file


    cpdef initialize(self, ReadProfiles.ReadProfiles pf)
    cpdef initialize_profiles(self, ReadProfiles.ReadProfiles pf)
    cpdef update(self, ReadProfiles.ReadProfiles pf, )
    cdef update_RRTM(self, ReadProfiles.ReadProfiles pf)
    # cpdef stats_io(self, NetCDFIO_Stats NS)
