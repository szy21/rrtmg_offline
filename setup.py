from distutils.core import setup
from Cython.Build import cythonize
from distutils.extension import Extension
import numpy as np
import mpi4py as mpi4py
import sys
import platform
import subprocess as sp
import os.path
import string

# Now get include paths from relevant python modules
include_path = [mpi4py.get_include()]
include_path += [np.get_include()]

if sys.platform == 'darwin':
    #Compile flags for MacOSX
    library_dirs = []
    libraries = []
    extensions = []
    extra_compile_args = []
    extra_compile_args += ['-O3', '-march=native', '-Wno-unused', '-Wno-#warnings','-fPIC']
    extra_objects=['./RRTMG/rrtmg_build/rrtmg_combined.o']
    netcdf_include = '/opt/local/include'
    netcdf_lib = '/opt/local/lib'
    f_compiler = 'gfortran'
    # f_compiler = 'gfortran-mp-5'
elif 'euler' in platform.node():
    #Compile flags for euler @ ETHZ
    library_dirs = ['/cluster/apps/openmpi/1.6.5/x86_64/gcc_4.8.2/lib/']
    libraries = []
    libraries.append('mpi')
    libraries.append('gfortran')
    extensions = []
    extra_compile_args=[]
    extra_compile_args+=['-std=c99', '-O3', '-march=native', '-Wno-unused',
                         '-Wno-#warnings', '-Wno-maybe-uninitialized', '-Wno-cpp', '-Wno-array-bounds','-fPIC']
    extra_objects=['./RRTMG/rrtmg_build/rrtmg_combined.o']
    netcdf_include = '/cluster/apps/netcdf/4.3.1/x86_64/gcc_4.8.2/openmpi_1.6.5/include'
    netcdf_lib = '/cluster/apps/netcdf/4.3.1/x86_64/gcc_4.8.2/openmpi_1.6.5/lib'
    f_compiler = 'gfortran'
elif platform.machine()  == 'x86_64':
    #Compile flags for fram @ Caltech
    library_dirs = string.split(os.environ['LD_LIBRARY_PATH'],':')
    libraries = []
    libraries.append('mpi')
    libraries.append('gfortran')
    extensions = []
    extra_compile_args=[]
    extra_compile_args+=['-std=c99', '-O3', '-march=native', '-Wno-unused',
                         '-Wno-#warnings', '-Wno-maybe-uninitialized', '-Wno-cpp', '-Wno-array-bounds','-fPIC']
    extra_objects=['./RRTMG/rrtmg_build/rrtmg_combined.o']
    netcdf_include = '/share/apps/software/rhel6/software/netCDF/4.4.0-foss-2016a/include'
    netcdf_lib = '/share/apps/software/rhel6/software/netCDF/4.4.0-foss-2016a/lib'
    f_compiler = 'gfortran'

else:
    print('Unknown system platform: ' + sys.platform + 'or unknown system name: ' + platform.node())
    sys.exit()

# _ext = Extension('mlm_thermodynamic_functions', ['mlm_thermodynamic_functions.pyx'], include_dirs=include_path,
#                  extra_compile_args=extra_compile_args, libraries=libraries, library_dirs=library_dirs,
#                  runtime_library_dirs=library_dirs, extra_objects=extra_objects)
# extensions.append(_ext)

_ext = Extension('ReadProfiles', ['ReadProfiles.pyx'], include_dirs=include_path,
                 extra_compile_args=extra_compile_args, libraries=libraries, library_dirs=library_dirs,
                 runtime_library_dirs=library_dirs, extra_objects=extra_objects)
extensions.append(_ext)

_ext = Extension('Radiation', ['Radiation.pyx'], include_dirs=include_path,
                 extra_compile_args=extra_compile_args, libraries=libraries, library_dirs=library_dirs,
                 runtime_library_dirs=library_dirs, extra_objects=extra_objects)
extensions.append(_ext)
#
# _ext = Extension('TimeStepping', ['TimeStepping.pyx'], include_dirs=include_path,
#                  extra_compile_args=extra_compile_args, libraries=libraries, library_dirs=library_dirs,
#                  runtime_library_dirs=library_dirs, extra_objects=extra_objects)
# extensions.append(_ext)

_ext = Extension('Simulation', ['Simulation.pyx'], include_dirs=include_path,
                 extra_compile_args=extra_compile_args, libraries=libraries, library_dirs=library_dirs,
                 runtime_library_dirs=library_dirs, extra_objects=extra_objects)
extensions.append(_ext)

# _ext = Extension('NetCDFIO', ['NetCDFIO.pyx'], include_dirs=include_path,
#                  extra_compile_args=extra_compile_args, libraries=libraries, library_dirs=library_dirs,
#                  runtime_library_dirs=library_dirs)
# extensions.append(_ext)

#Build RRTMG

rrtmg_compiled = os.path.exists('./RRTMG/rrtmg_build/rrtmg_combined.o')
if not rrtmg_compiled:
    run_str = 'cd ./RRTMG; '
    run_str += ('FC='+ f_compiler + ' LIB_NETCDF=' + netcdf_lib + ' INC_NETCDF='+
               netcdf_include + ' csh ./compile_RRTMG_combined.csh')
    print run_str
    sp.call([run_str], shell=True)
else:
    print("RRTMG Seems to be already compiled.")



setup(
    ext_modules=cythonize(extensions, verbose=1, include_path=include_path)
)

