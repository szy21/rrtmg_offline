import subprocess
import glob
import xarray as xr
import numpy as np

def main():
    # schemes = [1, 2, 3, 4, 5, 6, 7, 8, 9]
    tau0s = ['0_40x','0_60x','0_80x','1_00x','1_20x','1_40x','1_60x']
    forcing_fname = '/export/data1/zhaoyi/GCMForcedLES/forcing/new_0_80x_default.nc'
    input_data = xr.open_dataset(forcing_fname)
    lons = input_data.variables['lons'].values[::4]

    for tau0 in tau0s:
        for lon in lons[0:1]:
            lon_str = str(np.round(lon,1))
            casename = tau0+'_'+lon_str+'_250m'
            f = glob.glob('/export/data1/zhaoyi/GCMForcedLES/'+tau0+'/250m/Output.'+casename+'.4x/stats/Stats.*')
            f_rrtm = glob.glob('/export/data1/zhaoyi/GCMForcedLES/'+tau0+'/250m/Output.'+casename+'.4x/stats/RRTM_Stats*')
            if not f:
                print casename+': does not exist!'
            elif not f_rrtm:
            #if True:
                run_str = 'python main.py input/'+casename+'_nofix.in > log/'+casename+'_nofix.log'
                print run_str
                subprocess.call([run_str], shell=True)
            else:
                print casename+': rrtmg already done!'

    return

if __name__ == "__main__":
    main()
