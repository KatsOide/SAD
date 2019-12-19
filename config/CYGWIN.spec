# New Makefile -- OS Dependence Part for Cygwin/i386
#
SYS_LIBCRYPT=-lcrypt

# Disable DynamicLink functionality,
# because extension module building on cygwin is not supported.
UNLINK_SAD_FUNCTIONS+=DynamicLink
#OBJDYNL=sim/dynl-dl.o

_USE_CYGWIN=$(shell echo "$(OS_MAJOR_VERSION).000$(OS_MINOR_VERSION).000" | sed -e 's/\.[^.]*\([^.][^.][^.]\)/\1/g' -e 's/^0*\(.\)/\1/g')

ifeq ($(shell test $(_USE_CYGWIN) -lt 1007000 ; echo $$?), 0)
CONFIG_UNSUPPORTED=cygwin before version 1.7 does not support posix_openpt, RFC 3493 network API, etc...
CONFIG_RECOMMEND=Use cygwin 1.7.1 or later
endif

ifeq ($(COMPILER),System)
override COMPILER=GCC
endif

include $(SADCONF)/GCC.spec

ifeq ($(COMPILER),GCC)
# Use gcc4 compiler packages
__GCC_DEFAULT?=4.4
SYS_CC=gcc-4
SYS_FC=gfortran-4
SYS_CXX=g++-4
endif

# Override -std=c99, because newlib is not compatible with ISO C99.
SYS_COPT+=	-std=gnu99

X11_PREFIX=/usr/X11R6

TCLTK_CFLAGS+=-DSAD_FORCE_X11 -DHAVE_TM_ZONE
SYS_COPT_tfTkInter_+=-DSAD_FORCE_X11

# End of File
