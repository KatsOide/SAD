# Generic Object Makefile OS Compiler Support Part -- SunStudio
#

ifeq ($(DEBUG),SunStudio)
debug:
	@echo "SYS_FOPT_ENABLE_BACKSLASH_ESCAPE=$(SYS_FOPT_ENABLE_BACKSLASH_ESCAPE)"
	@echo "SYS_FOPT_DISABLE_BACKSLASH_ESCAPE=$(SYS_FOPT_DISABLE_BACKSLASH_ESCAPE)"
endif

# Define missing Fortran2008 standard math functions
FORTRAN_MISSING_DMATH_FUNCS+=	hypot
# REAL(8) instance of asinh/acosh/atanh of Fortran2008
FORTRAN_MISSING_DMATH_FUNCS+=	asinh acosh atanh

HAVE_F_FORMAT_Q_EDIT=YES
HAVE_F_FGETC=YES
HAVE_F_GETFD=YES

SYS_FOPT_ENABLE_BACKSLASH_ESCAPE=-f77=%backslash
SYS_FOPT_DISABLE_BACKSLASH_ESCAPE=#-f77=no%backslash

SYS_FC=sunf95
SYS_CC=suncc
SYS_CXX=sunCC

SYS_FOPT=-xtypemap=integer:32,real:32,double:64 \
	-ftrap=%none -stackvar \
	-xalias=dummy,overindex
SYS_COPT=-xc99=all
SYS_CXXOPT=

SYS_FOPT+=$(SYS_FOPT_DISABLE_BACKSLASH_ESCAPE)

#SYS_FOPT+=-m64 -xmodel=small

# For linking ISO C99 library by Fortran compiler
_SYS_EXTSIM_XPG6=
ifeq ($(OS_TYPE),SunOS)
_SYS_EXTSIM_XPG6=/usr/lib/values-xpg6.o
endif
SYS_EXTSIM+=$(_SYS_EXTSIM_XPG6)

# Default debug build options
LIB_COPT=	-g -xO1
LIB_CXXOPT=	-g -xO1
FOPT=		-g -xO1
COPT=		-g -xO1
CXXOPT=		-g -xO1

# End of File
