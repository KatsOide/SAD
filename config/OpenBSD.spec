# New Makefile -- OS Dependence Part for OpenBSD
#

# Setup GCC for packages
__GCC_DEFAULT=4.9.0

GCC_FC_GFORTRAN?=egfortran
GCC_CC?=egcc
GCC_CXX?=eg++

-include $(SADCONF)/4.4BSD-Lite.spec

SYS_LIBCRYPT=

# End of File
