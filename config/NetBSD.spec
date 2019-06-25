# New Makefile -- OS Dependence Part for NetBSD
#

# Setup GCC for pkgsrc
__GCC_DEFAULT=4.9.1

__GCC_SUFFIX=$(shell        echo $(_USE_GCC)        | sed -e 's/\(...\)...$$/.\1/' -e 's/\.0*//' -e 's/^0$$//')
__GCC_SUFFIX_AS_CC=$(shell  echo $(_USE_GCC_AS_CC)  | sed -e 's/\(...\)...$$/.\1/' -e 's/\.0*//' -e 's/^0$$//')
__GCC_SUFFIX_AS_CXX=$(shell echo $(_USE_GCC_AS_CXX) | sed -e 's/\(...\)...$$/.\1/' -e 's/\.0*//' -e 's/^0$$//')

GCC_FC_GFORTRAN?=/usr/pkg/gcc$(__GCC_SUFFIX)/bin/gfortran
GCC_CC?=/usr/pkg/gcc$(__GCC_SUFFIX_AS_CC)/bin/gcc
GCC_CXX?=/usr/pkg/gcc$(__GCC_SUFFIX_AS_CXX)/bin/g++

-include $(SADCONF)/4.4BSD-Lite.spec

HAVE_MKSTEMPS=NO
INSTALL=install -c -r

X11_PREFIX=/usr/X11R7

# End of File
