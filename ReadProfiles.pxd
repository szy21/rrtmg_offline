cimport Radiation
# cimport TimeStepping
# from NetCDFIO cimport NetCDFIO_Stats

cdef class ReadProfiles:
    cdef:

        Py_ssize_t nz
        double t_surface
        double [:] pressure
        double [:] pressure_i
        double [:] temperature
        double [:] qv
        double [:] ql
        double [:] qi
        double [:] cf
        double [:] qt
        double [:] rho
        double toa_sw
        double albedo

        Py_ssize_t t1
        Py_ssize_t t2

        str path_plus_file
        str path_plus_file_control
        str path_plus_file_ref
        object root_grp
        object root_grp_ref
        object control_grp
        object profile_grp
        object ref_grp
        object ts_grp
        object profile_grp2
        object ref_grp2
        object ts_grp2
        str path

        bint average
        str out_file

        bint fix_T
        bint fix_qv
        bint fix_cloud
        bint fix_albedo

    cdef public:
        Py_ssize_t count
        Py_ssize_t ntime


    cpdef initialize(self)
    cpdef update(self, Radiation.Radiation Ra)
    # cpdef stats_io(self, NetCDFIO_Stats NS)
