# New Makefile -- OS Dependence Part for FreeBSD
#

# Setup GCC for ports-CURRENT
__GCC_DEFAULT=4.8.4

__GCC_SUFFIX=$(shell        echo $(_USE_GCC)        | sed -e 's/\(...\)...$$/.\1/' -e 's/\.0*//' -e 's/^0$$//')
__GCC_SUFFIX_AS_CC=$(shell  echo $(_USE_GCC_AS_CC)  | sed -e 's/\(...\)...$$/.\1/' -e 's/\.0*//' -e 's/^0$$//')
__GCC_SUFFIX_AS_CXX=$(shell echo $(_USE_GCC_AS_CXX) | sed -e 's/\(...\)...$$/.\1/' -e 's/\.0*//' -e 's/^0$$//')

GCC_FC_GFORTRAN?=gfortran$(__GCC_SUFFIX)
GCC_CC?=gcc$(__GCC_SUFFIX_AS_CC)
GCC_CXX?=g++$(__GCC_SUFFIX_AS_CXX)

ifeq ($(_COMPILER),System)
ifeq ($(shell test $(OS_MAJOR_VERSION) -ge 10 ; echo $$?), 0)
USE_GCC_AS_CC?=NO
endif
USE_GCC_AS_CXX?=NO
endif

-include $(SADCONF)/4.4BSD-Lite.spec

USE_TCLTK_TRIM_VERSION=YES

X11_PREFIX=/usr/local

##SYS_OBJS+=sim/linux_freebsd.o

ifeq ($(_COMPILER),System)
ifeq ($(OS_MAJOR_VERSION),9)
SYS_CXX=clang++
endif
endif

# End of File
