cimport ReadProfiles
# cimport TimeStepping
# cimport NetCDFIO
cimport Radiation
# import time
# import numpy as np


class Simulation:

    def __init__(self):
        return

    def initialize(self, namelist):
        self.pf = ReadProfiles.ReadProfiles(namelist)
        # self.TS = TimeStepping.TimeStepping()
        # self.StatsIO = NetCDFIO.NetCDFIO_Stats()
        self.Ra = Radiation.Radiation(namelist)

        # self.StatsIO.initialize(namelist, self.pf)
        # self.TS.initialize(namelist, self.pf)
        self.Ra.initialize(self.pf)
        self.pf.initialize()

        return

    def run(self):

        self.pf.update(self.Ra)
        self.Ra.initialize_profiles(self.pf)
        self.pf.count = 0
        while self.pf.count < self.pf.ntime:
            self.Ra.update(self.pf)
            self.pf.update(self.Ra)
        # self.pf.update(self.TS, self.Ra)
                # self.TS.update_second(self.pf)
                # self.io()
            # time2 = time.time()
            # print('T = ' + str(self.TS.t) + ' dt = ' + str(self.TS.dt) + ' walltime = ' + str(time2 - time1))
        # self.pf.root_grp.close()
        print('Finished RRTM!')
        return

    # def io(self):
    #     cdef:
    #         double stats_dt = 0.0
    #
    #     if self.TS.t > 0.0 and self.TS.rk_step == self.TS.n_rk_steps - 1:
    #         stats_dt = self.StatsIO.last_output_time + self.StatsIO.frequency - self.TS.t
    #
    #         if self.StatsIO.last_output_time + self.StatsIO.frequency == self.TS.t:
    #             self.StatsIO.last_output_time = self.TS.t
    #             self.StatsIO.open_files()
    #             self.StatsIO.write_simulation_time(self.TS.t)
    #
    #             self.mlm.stats_io(self.StatsIO)
    #             self.Ra.stats_io(self.StatsIO)
    #
    #             self.StatsIO.close_files()
    #
    #     return

