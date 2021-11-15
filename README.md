GridForce plugin for OpenMM
===========================

`openmmgridforce` is a plugin for the `OpenMM` toolkit for molecular simulations using a potential energy on grid points. 
Currently, the `openmmgridforce` plugin is available within the "reference implementation" of the `OpenMM` toolkit.

## How to install from source

### Requirements
conda create -n omm3
conda activate omm3
conda install -c conda-forge openmm

conda install -c conda-forge cmake
conda install -c conda-forge swig
conda install -c conda-forge netcdf4
#### macOS-64
conda install -c conda-forge clang_osx-64
conda install -c conda-forge clangxx_osx-64
#### linux-64
conda install -c conda-forge gcc_linux-64
conda install -c conda-forge gxx_linux-64

### (Compile and install C++ codes) ###
0) The current directory is at ./openmmgridforce
Modify 'CMakeLists.txt' as follows:
---(update where OpenMM is installed)---
SET(OPENMM_DIR "/Users/willow/opt/miniconda3/envs/omm3" CACHE PATH "Where OpenMM is installed")
---
1) make a build directory
mkdir ../openmmgridforce_build
cd ../openmmgridforce_build
2) type cmake
cmake ../openmmgridforce
4) type make
make
5) install (files will be installed into CMAKE_INSTALL_PREFIX)
make install

### (Compile and install python wrapper) ###
0) The current directory is at ./openmmgridforce/python
1) swig: gridforceplugin.i --> GridForcePluginWrapper.cpp
(You need to provide a correct header directory where OpenMM and openmmgridforce are installed.)
swig -python -c++ -o GridForcePluginWrapper.cpp -I${OPENMM_DIR}/include gridforceplugin.i
2) python setup.py build
3) python setup.py install
# Now you can testify the Test*py.
4) python TestReferenceGridForce.py
