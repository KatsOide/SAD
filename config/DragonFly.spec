# New Makefile -- OS Dependence Part for DragonFly BSD
#

# Setup GCC for dports
__GCC_DEFAULT=4.9.3

__GCC_SUFFIX=$(shell        echo $(_USE_GCC)        | sed -e 's/\(...\)...$$/.\1/' -e 's/\.0*//' -e 's/^0$$//')
__GCC_SUFFIX_AS_CC=$(shell  echo $(_USE_GCC_AS_CC)  | sed -e 's/\(...\)...$$/.\1/' -e 's/\.0*//' -e 's/^0$$//')
__GCC_SUFFIX_AS_CXX=$(shell echo $(_USE_GCC_AS_CXX) | sed -e 's/\(...\)...$$/.\1/' -e 's/\.0*//' -e 's/^0$$//')

GCC_FC_GFORTRAN=gfortran$(__GCC_SUFFIX)
GCC_CC=gcc$(__GCC_SUFFIX_AS_CC)
GCC_CXX=g++$(__GCC_SUFFIX_AS_CXX)

-include $(SADCONF)/4.4BSD-Lite.spec

X11_PREFIX=/usr/local

# End of File
