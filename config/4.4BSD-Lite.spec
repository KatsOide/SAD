# New Makefile -- OS Dependence Part for 4.4BSD-Lite based PC-UNIX
#

HAVE_MKSTEMPS=YES
SYS_LIBCRYPT=-lcrypt
FETCH=ftp
INSTALL=install -c -S

OBJDYNL=sim/dynl-dl.o

ifeq ($(COMPILER),System)
override COMPILER=GCC
endif

# Compiler Dependence
#
SETUP_OSTYPE_DEP_COMPILER=NO

ifeq ($(COMPILER),GCC)
include $(SADCONF)/GCC.spec
ifeq ($(_COMPILER),System)
ifeq ($(_USE_GCC_AS_CC),0)
SYS_CC=cc
endif
ifeq ($(_USE_GCC_AS_CXX),0)
SYS_CXX=c++
endif
endif

SETUP_OSTYPE_DEP_COMPILER=YES
endif

ifeq ($(COMPILER),ICC)
include $(SADCONF)/ICC.spec
# many 4.4BSD-Lite variant use GNU-binutils's ld as system linker
SYS_LDOPT_RUNTIME_PATH_OPT=-Wl,-rpath,
SYS_LDOPT_EXPORT_SYMBOL=-Wl,--export-dynamic

SETUP_OSTYPE_DEP_COMPILER=YES
endif

ifeq ($(COMPILER),NAG)
include $(SADCONF)/GCC.spec
include $(SADCONF)/NAG.spec

ifeq ($(_COMPILER),System)
SYS_CC=cc
SYS_CXX=c++
endif

SETUP_OSTYPE_DEP_COMPILER=YES
endif

X11_PREFIX=/usr/X11R6

# End of File
