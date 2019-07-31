# New Makefile -- OS Dependence Part for SunOS
#
SYS_LIBCRYPT=-lcrypt

OBJDYNL=sim/dynl-dl.o

ifeq ($(OS_MAJOR_VERSION),4)
USE_TCLTK_TRIM_VERSION=YES
endif

# Compiler Dependence
#
SETUP_OSTYPE_DEP_COMPILER=NO

ifeq ($(COMPILER),GCC)
include $(SADCONF)/GCC.spec
SYS_COPT=-DBSD_COMP
SYS_LDOPT_RUNTIME_PATH_OPT=-R

SETUP_OSTYPE_DEP_COMPILER=YES
endif

ifeq ($(COMPILER),System)
include $(SADCONF)/SunStudio.spec
SYS_LDOPT_RUNTIME_PATH_OPT=-R

SETUP_OSTYPE_DEP_COMPILER=YES
endif

X11_PREFIX=/usr/X11R6
SYS_LIBS=-lsocket -lnsl -ldl

# End of File
