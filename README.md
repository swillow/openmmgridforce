GridForce plugin for OpenMM
===========================

`openmmgridforce` is a plugin for the `OpenMM` toolkit for molecular simulations using a potential energy on grid points. 
Currently, the `openmmgridforce` plugin is available within the "reference implementation" of the `OpenMM` toolkit.

## How to install from source

### Requirements (Note: openmm used SWIG 3 for all builds)
* conda create -n omm3
* conda activate omm3
* conda install -c conda-forge openmm cmake netcdf4
* conda install -c conda-forge swig=3.0.12

#### on macOS-64
* conda install -c conda-forge clang_osx-64 clangxx_osx-64

#### on linux-64
* conda install -c conda-forge gcc_linux-64 gxx_linux-64


### (Compile and install C++ codes) ###
* The current directory is at ./openmmgridforce. Modify 'CMakeLists.txt' as follows:

---(update where OpenMM is installed)---

SET(OPENMM_DIR "/Users/willow/opt/miniconda3/envs/omm3" CACHE PATH "Where OpenMM is installed")

* make a build directory

mkdir ../openmmgridforce_build

cd ../openmmgridforce_build
* type cmake

cmake ../openmmgridforce
* type make

make
* install (files will be installed into CMAKE_INSTALL_PREFIX)

make install

### (Compile and install python wrapper) ###
* The current directory is at ./openmmgridforce/python
* swig: gridforceplugin.i --> GridForcePluginWrapper.cpp

(You need to provide a correct header directory where OpenMM and openmmgridforce are installed.)

swig -python -c++ -o GridForcePluginWrapper.cpp -I${OPENMM_DIR}/include gridforceplugin.i
* python setup.py build
* python setup.py install
* Now you can testify the Test*py.

python TestReferenceGridForce.py
