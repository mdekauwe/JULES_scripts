#!/bin/bash

# Wrapper to build JULES locally outside of rose,cycl
# Needs to go in the top of the src directory you checkout.
# Works with gcc8/gfortran 8, not the newer versions
#
# Martin De Kauwe, 14 Oct 2020

rm -rf build extract preprocess

export JULES_NETCDF=netcdf
export JULES_REMOTE=local
export JULES_BUILD=normal
export JULES_COMPILER=gfortran
export JULES_OMP=noomp
export JULES_MPI=nompi
export JULES_FFLAGS_EXTRA=
export JULES_LDFLAGS_EXTRA=
export JULES_NETCDF_INC_PATH=/opt/local/include/
export JULES_NETCDF_LIB_PATH=/opt/local/lib/

# Need this for some issue with "rpath"
export ncdf_ldflags_dynamic=

fcm make -v -f etc/fcm-make/make.cfg --new
#fcm make -v -f etc/fcm-make/make.cfg $1
##fcm make -vv -f etc/fcm-make/make.cfg $1
