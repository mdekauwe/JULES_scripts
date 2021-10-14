#!/bin/bash

# Wrapper to build JULES locally outside of rose,cycl
# Needs to go in the top of the src directory you checkout.
#
#
# Martin De Kauwe, 14 Oct 2020

export JULES_COMPILER=gfortran
export JULES_BUILD=normal
export JULES_OMP=noomp
export JULES_MPI=nompi
export JULES_NETCDF=netcdf
export JULES_NETCDF_INC_PATH=/opt/local/include/
export JULES_NETCDF_LIB_PATH=/opt/local/lib/
export JULES_REMOTE=local
export JULES_FFLAGS_EXTRA=
export JULES_LDFLAGS_EXTRA=

fcm make -f etc/fcm-make/make.cfg $1
#fcm make -vv -f etc/fcm-make/make.cfg $1
