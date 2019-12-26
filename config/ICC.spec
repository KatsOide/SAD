# Generic Object Makefile OS Compiler Support Part -- ICC
#

ifeq ($(DEBUG),ICC)
debug:
	@echo "_USE_IFC=$(_USE_IFC)"
	@echo "_USE_ICC=$(_USE_ICC)"
	@echo "SYS_FOPT_ENABLE_BACKSLASH_ESCAPE=$(SYS_FOPT_ENABLE_BACKSLASH_ESCAPE)"
	@echo "SYS_FOPT_DISABLE_BACKSLASH_ESCAPE=$(SYS_FOPT_DISABLE_BACKSLASH_ESCAPE)"
endif

# _USE_* <- major * 10^6 + minor * 10^3 + patch-level[minor,patch-level < 10^3]
ifdef USE_IFC
_USE_IFC=$(shell echo "$(USE_IFC)" | sed -e 's/ *//g' -e 's/[^0-9.].*//' -e 's/\.*$$/.0.0.0./' -e 's/^\.*//' -e 's/^\([^.][^.]*\)\.\([^.][^.]*\)\.\([^.][^.]*\)\..*/\1.000\2.000\3/' -e 's/\.[^.]*\([^.][^.][^.]\)/\1/g' -e 's/^0*\(.\)/\1/g')
else
_USE_IFC=0
endif

ifdef USE_ICC
_USE_ICC=$(shell echo "$(USE_ICC)" | sed -e 's/ *//g' -e 's/[^0-9.].*//' -e 's/\.*$$/.0.0.0./' -e 's/^\.*//' -e 's/^\([^.][^.]*\)\.\([^.][^.]*\)\.\([^.][^.]*\)\..*/\1.000\2.000\3/' -e 's/\.[^.]*\([^.][^.][^.]\)/\1/g' -e 's/^0*\(.\)/\1/g')
else
_USE_ICC=0
endif

# Default Intel compiler version: 8.001.000
ifeq ($(_USE_IFC),0)
ifneq ($(_USE_ICC),0)
_USE_IFC=$(_USE_ICC)
else
_USE_IFC=8001000
endif
endif

ifeq ($(_USE_ICC),0)
_USE_ICC=$(_USE_IFC)
endif

_ICC_ARCH_SUPPORTED=NO
ifeq ($(CPU_ARCH),i386)
_ICC_ARCH_SUPPORTED=YES
endif

ifeq ($(CPU_ARCH),AMD64)
_ICC_ARCH_SUPPORTED=YES
endif

ifeq ($(_ICC_ARCH_SUPPORTED),NO)
$(error Architecture $(CPU_ARCH) is not supported by current ICC.spec!)
endif

# Define missing Fortran2008 standard math functions
ifeq ($(shell test $(_USE_IFC) -lt 99000000 ; echo $$?), 0)
FORTRAN_MISSING_DMATH_FUNCS+=	hypot
# REAL(8) instance of asinh/acosh/atanh of Fortran2008
FORTRAN_MISSING_DMATH_FUNCS+=	asinh acosh atanh
endif

HAVE_F_FORMAT_Q_EDIT=YES
HAVE_F_FGETC=YES
HAVE_F_PXFFILENO=YES

SYS_FC=ifort
ifeq ($(shell test $(_USE_IFC) -lt 8000000 ; echo $$?), 0)
SYS_FC=ifc
endif
SYS_CC=icc
SYS_CXX=icpc

SYS_FOPT_DYNL_PIC=   -fPIC
SYS_COPT_DYNL_PIC=   -fPIC
SYS_CXXOPT_DYNL_PIC= -fPIC

SYS_FOPT=   -falias -ffnalias -i4 -r8 -dps
SYS_COPT=   -falias -ffnalias -std=c99
SYS_CXXOPT= -falias -ffnalias

# IFC 8.0 or later treat backslash as normal character in default
SYS_FOPT_ENABLE_BACKSLASH_ESCAPE=-assume bscc
SYS_FOPT_DISABLE_BACKSLASH_ESCAPE=
ifeq ($(shell test $(_USE_IFC) -lt 8000000 ; echo $$?), 0)
# IFC 7.x or prior treat backslash as escape character in default
SYS_FOPT_ENABLE_BACKSLASH_ESCAPE=
SYS_FOPT_DISABLE_BACKSLASH_ESCAPE=-nbs
endif

SYS_FOPT+=$(SYS_FOPT_DISABLE_BACKSLASH_ESCAPE)

# Use octet-bytes alignment for double on IA-32 architecture,
# because of alignment assumption of SAD core code.(Ex. rlist variable)
SYS_FOPT+=   -align -Zp8 -pc80
SYS_COPT+=   -align -Zp8 -pc80
SYS_CXXOPT+= -align -Zp8 -pc80

ifeq ($(shell test $(_USE_IFC) -lt 9000000 ; echo $$?), 0)
SYS_LDOPT=-Vaxlib -posixlib
endif

# Object specific options
SYS_FOPT_tfconvstr=-auto

# Default build options
LIB_COPT=	-O2
LIB_CXXOPT=	-O2
FOPT=		-O3 -ip -unroll
COPT=		-O3 -ip -unroll
CXXOPT=		-O3 -ip -unroll

# Pentium 4/SSE3 optimization
_ARCH_OPTIMIZE_FOPT=	-mtune=pentium4 -xP
ifeq ($(shell test $(_USE_IFC) -lt 9000000 ; echo $$?), 0)
_ARCH_OPTIMIZE_FOPT=	-tpp7 -xP
endif

_ARCH_OPTIMIZE_COPT=	-mtune=pentium4 -xP
ifeq ($(shell test $(_USE_ICC) -lt 9000000 ; echo $$?), 0)
_ARCH_OPTIMIZE_COPT=	-tpp7 -xP
endif

FOPT+=	 $(_ARCH_OPTIMIZE_FOPT)
COPT+=	 $(_ARCH_OPTIMIZE_COPT)
CXXOPT+= $(_ARCH_OPTIMIZE_COPT)

# End of File
