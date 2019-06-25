# New Makefile -- OS Dependence Part for Linux
#
SYS_LIBCRYPT=-lcrypt

OBJDYNL=sim/dynl-dl.o

ifeq ($(COMPILER),System)
override COMPILER=GCC
endif

# Compiler Dependence
#
SETUP_OSTYPE_DEP_COMPILER=NO

ifeq ($(COMPILER),GCC)
__GCC_DEFAULT?=4.4.3

include $(SADCONF)/GCC.spec

SYS_FC=f95
SYS_CC=cc
SYS_CXX=c++
ifeq ($(shell test $(_USE_GCC) -ge 4005000 ; echo $$?), 0)
__GCC_SUFFIX=$(shell echo -$(_USE_GCC) | sed -e 's/\(...\)...$$/.\1/' -e 's/\.0*/./')
SYS_FC=gfortran$(__GCC_SUFFIX)
SYS_CC=gcc$(__GCC_SUFFIX)
SYS_CXX=g++$(__GCC_SUFFIX)
endif

SETUP_OSTYPE_DEP_COMPILER=YES
endif

ifeq ($(COMPILER),ICC)
include $(SADCONF)/ICC.spec

SYS_LDOPT_RUNTIME_PATH_OPT=-Wl,-rpath,
SYS_LDOPT_EXPORT_SYMBOL=-Wl,--export-dynamic

SETUP_OSTYPE_DEP_COMPILER=YES
endif

ifeq ($(COMPILER),PGI)
include $(SADCONF)/PGI.spec

SYS_LDOPT_RUNTIME_PATH_OPT=-Wl,-rpath,
SYS_LDOPT_EXPORT_SYMBOL=-Wl,--export-dynamic

SETUP_OSTYPE_DEP_COMPILER=YES
endif

ifeq ($(COMPILER),SunStudio)
include $(SADCONF)/SunStudio.spec

SYS_LDOPT_RUNTIME_PATH_OPT=-Wl,-rpath,
SYS_LDOPT_EXPORT_SYMBOL=-Wl,--export-dynamic

SETUP_OSTYPE_DEP_COMPILER=YES
endif

SYS_COPT+=-D_POSIX_C_SOURCE=200112L -D_XOPEN_SOURCE=600

# for NI_MAXHOST/NI_MAXSERV
SYS_COPT_tfNetworkIO_=-D_BSD_SOURCE

# for wait4
SYS_COPT_tfProcess_=-D_BSD_SOURCE

# for mkdtemp
SYS_COPT_sim/fortran_io_+=-D_BSD_SOURCE

# for MAP_ANONYMOUS/MAP_32BIT
SYS_COPT_sim/unix_memory_+=-D_SVID_SOURCE
SYS_COPT_sim/unix_memory8_+=-D_SVID_SOURCE

# for RTLD_DEFAULT
SYS_COPT_sim/dynl-dl+=-D_GNU_SOURCE

# for getmode/setmode(3)
SYS_COPT_sim/setmode=-DS_ISTXT=0001000 -Du_int='unsigned int'

X11_PREFIX=/usr/X11R6
SYS_LIBS=-ldl

# End of File
