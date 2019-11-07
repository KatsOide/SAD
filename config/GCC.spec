# Generic Object Makefile OS Compiler Support Part -- GCC
#

# Default GCC version
__GCC_DEFAULT?=4.4.7

# Recommended GCC version
__GCC_RECOMMEND=4.9.1
__GCC_RECOMMEND_AS_CC=4.9.1
__GCC_RECOMMEND_AS_CXX=4.9.1

ifeq ($(DEBUG),GCC)
debug:
	@echo "_USE_GCC_DEFAULT=$(_USE_GCC_DEFAULT)"
	@echo "_USE_GCC=$(_USE_GCC)"
	@echo "_USE_G95=$(_USE_G95)"
	@echo "_USE_GCC_AS_CC=$(_USE_GCC_AS_CC)"
	@echo "_USE_GCC_AS_CXX=$(_USE_GCC_AS_CXX)"
	@echo "_USE_GNU_FORTRAN=$(_USE_GNU_FORTRAN)"
	@echo "GCC_FC_G95=$(GCC_FC_G95)"
	@echo "GCC_FC_GFORTRAN=$(GCC_FC_GFORTRAN)"
	@echo "GCC_CC=$(GCC_CC)"
	@echo "GCC_CXX=$(GCC_CXX)"
	@echo "SYS_FOPT_ENABLE_BACKSLASH_ESCAPE=$(SYS_FOPT_ENABLE_BACKSLASH_ESCAPE)"
	@echo "SYS_FOPT_DISABLE_BACKSLASH_ESCAPE=$(SYS_FOPT_DISABLE_BACKSLASH_ESCAPE)"
	@echo "_GCC_RLIBDIR=$(_GCC_RLIBDIR)"
	@echo "_GCC_RLIBDIR_AS_CC=$(_GCC_RLIBDIR_AS_CC)"
	@echo "_GCC_RLIBDIR_AS_CXX=$(_GCC_RLIBDIR_AS_CXX)"
	@echo "OPT_RLIBDIR=$(OPT_RLIBDIR)"
	@echo "_SYS_GCC_ABIOPT=$(_SYS_GCC_ABIOPT)"
endif

# _USE_* <- major * 10^6 + minor * 10^3 + patch-level[minor,patch-level < 10^3]

# _USE_GCC_DEFAULT <- __GCC_DEFAULT
_USE_GCC_DEFAULT=$(shell echo "$(__GCC_DEFAULT)" | sed -e 's/ *//g' -e 's/[^0-9.].*//' -e 's/\.*$$/.0.0.0./' -e 's/^\.*//' -e 's/^\([^.][^.]*\)\.\([^.][^.]*\)\.\([^.][^.]*\)\..*/\1.000\2.000\3/' -e 's/\.[^.]*\([^.][^.][^.]\)/\1/g' -e 's/^0*\(.\)/\1/g')

_USE_GCC=$(_USE_GCC_DEFAULT)
ifdef USE_GCC
_USE_GCC=$(shell echo "$(USE_GCC)" | sed -e 's/ *//g' -e 's/[^0-9.].*//' -e 's/\.*$$/.0.0.0./' -e 's/^\.*//' -e 's/^\([^.][^.]*\)\.\([^.][^.]*\)\.\([^.][^.]*\)\..*/\1.000\2.000\3/' -e 's/\.[^.]*\([^.][^.][^.]\)/\1/g' -e 's/^0*\(.\)/\1/g')
endif

_USE_G95=0
ifdef USE_G95
_USE_G95=$(shell echo "$(USE_G95)" | sed -e 's/ *//g' -e 's/[^0-9.].*//' -e 's/\.*$$/.0.0.0./' -e 's/^\.*//' -e 's/^\([^.][^.]*\)\.\([^.][^.]*\)\.\([^.][^.]*\)\..*/\1.000\2.000\3/' -e 's/\.[^.]*\([^.][^.][^.]\)/\1/g' -e 's/^0*\(.\)/\1/g')
endif

_USE_GCC_AS_CC=$(_USE_GCC)
ifdef USE_GCC_AS_CC
_USE_GCC_AS_CC=$(shell echo "$(USE_GCC_AS_CC)" | sed -e 's/ *//g' -e 's/[^0-9.].*//' -e 's/\.*$$/.0.0.0./' -e 's/^\.*//' -e 's/^\([^.][^.]*\)\.\([^.][^.]*\)\.\([^.][^.]*\)\..*/\1.000\2.000\3/' -e 's/\.[^.]*\([^.][^.][^.]\)/\1/g' -e 's/^0*\(.\)/\1/g')
ifeq ($(USE_GCC_AS_CC),YES)
_USE_GCC_AS_CC=$(_USE_GCC)
endif
ifeq ($(USE_GCC_AS_CC),NO)
_USE_GCC_AS_CC=0
endif
endif	# defined(USE_GCC_AS_CC)

_USE_GCC_AS_CXX=$(_USE_GCC)
ifdef USE_GCC_AS_CXX
_USE_GCC_AS_CXX=$(shell echo "$(USE_GCC_AS_CXX)" | sed -e 's/ *//g' -e 's/[^0-9.].*//' -e 's/\.*$$/.0.0.0./' -e 's/^\.*//' -e 's/^\([^.][^.]*\)\.\([^.][^.]*\)\.\([^.][^.]*\)\..*/\1.000\2.000\3/' -e 's/\.[^.]*\([^.][^.][^.]\)/\1/g' -e 's/^0*\(.\)/\1/g')
ifeq ($(USE_GCC_AS_CXX),YES)
_USE_GCC_AS_CXX=$(_USE_GCC)
endif
ifeq ($(USE_GCC_AS_CXX),NO)
_USE_GCC_AS_CXX=0
endif
endif	# defined(USE_GCC_AS_CXX)

# Check ABI option
_SYS_GCC_ABIOPT=
ifeq ($(ABI),32)
_SYS_GCC_ABIOPT=-m32
endif
ifeq ($(ABI),64)
_SYS_GCC_ABIOPT=-m64
endif
SYS_CC_ABIOPT=$(_SYS_GCC_ABIOPT)
SYS_FC_ABIOPT=$(_SYS_GCC_ABIOPT)
SYS_CXX_ABIOPT=$(_SYS_GCC_ABIOPT)

## Check ISO C11 for SYS_CC
#ifeq ($(shell test $(_USE_GCC_AS_CC) -gt 0 -a $(_USE_GCC_AS_CC) -lt 4009001 ; echo $$?), 0)
#CONFIG_UNSUPPORTED=GNU C(gcc) 4.9.0 or prior does not support ISO C11 standard
#CONFIG_RECOMMEND=Use GCC $(__GCC_RECOMMEND_AS_CC) or later.
#endif

# Check GNU Fortran compiler type
_USE_GNU_FORTRAN=GFORTRAN
ifeq ($(shell test $(_USE_GCC) -lt 4000000 ; echo $$?), 0)
_USE_GNU_FORTRAN=G77
endif
ifeq ($(shell test $(_USE_G95) -gt 0 ; echo $$?), 0)
_USE_GNU_FORTRAN=G95
endif

# Check compiler dependent library location
ifeq ($(_USE_GNU_FORTRAN),GFORTRAN)
ifeq ($(shell test $(_USE_GCC) -gt 0 ; echo $$?), 0)
_GCC_RLIBDIR=$(foreach _dir,$(firstword $(filter-out .,$(foreach _file,libgfortran.so libgfortran.a,$(shell dirname `$(GCC_FC_GFORTRAN) -print-file-name=$(_file)`)))),-L$(realpath $(_dir)))
endif
endif
ifeq ($(shell test $(_USE_GCC_AS_CC) -gt 0 ; echo $$?), 0)
_GCC_RLIBDIR_AS_CC=$(foreach _dir,$(firstword $(filter-out .,$(foreach _file,libgcc_s.so libiberty.a,$(shell dirname `$(GCC_CC) -print-file-name=$(_file)`)))),-L$(realpath $(_dir)))
endif
ifeq ($(shell test $(_USE_GCC_AS_CXX) -gt 0 ; echo $$?), 0)
_GCC_RLIBDIR_AS_CXX=$(foreach _dir,$(firstword $(filter-out .,$(foreach _file,libstdc++.so libstdc++.a,$(shell dirname `$(GCC_CXX) -print-file-name=$(_file)`)))),-L$(realpath $(_dir)))
endif
OPT_RLIBDIR+= $(_GCC_RLIBDIR)

# Define missing Fortran2008 standard math functions
ifeq ($(shell test $(_USE_GCC) -lt 4004000 ; echo $$?), 0)
FORTRAN_MISSING_DMATH_FUNCS+=	hypot
# REAL(8) instance of asinh/acosh/atanh of Fortran2008
FORTRAN_MISSING_DMATH_FUNCS+=	asinh acosh atanh
endif

# Default setting for GNU Fortran compiler
HAVE_F_FORMAT_DOLLAR_EDIT=YES
HAVE_F_FSEEK=GNU
HAVE_F_FGETC=YES
HAVE_F_FNUM=YES

ifeq ($(_USE_GNU_FORTRAN),G77)
CONFIG_UNSUPPORTED=GNU Fortran77(g77) does not support Fortran95/2003/2008 standard
CONFIG_RECOMMEND=Use gfortran in GCC $(__GCC_RECOMMEND) or later instead of g77.
endif

ifeq ($(_USE_GNU_FORTRAN),GFORTRAN)
# GFORTRAN 4.2.x or prior is not supported
ifeq ($(shell test $(_USE_GCC) -lt 4003000 ; echo $$?), 0)
CONFIG_UNSUPPORTED=GNU Fortran(gfortran) 4.2.x or prior does not support FSEEK intrinsic and has known erratas.
CONFIG_RECOMMEND=Use GCC $(__GCC_RECOMMEND) or later.
endif

# Construct -fcheck option for gfortran
ifdef USE_GFORTRAN_CHECK
_USE_GFORTRAN_CHECK=$(shell echo $(USE_GFORTRAN_CHECK) | sed -e 's/[ 	][ 	]*/ /g' -e 's/^ *//' -e 's/ *$$//' -e 's/ /,/g')
SYS_FOPT_GFORTRAN_CHECK=
ifneq ($(_USE_GFORTRAN_CHECK),)
SYS_FOPT_GFORTRAN_CHECK=-fcheck=$(_USE_GFORTRAN_CHECK)
endif
endif

# -fno-backslash is default backslash handling in GFORTRAN 4.3.0 or later
SYS_FOPT_ENABLE_BACKSLASH_ESCAPE=-fbackslash
SYS_FOPT_DISABLE_BACKSLASH_ESCAPE=#-fno-backslash
endif

ifeq ($(_USE_GNU_FORTRAN),G95)
SYS_FOPT_tfprinta+= -fsloppy-char

# -fbackslash is default backslash handling in G95
SYS_FOPT_ENABLE_BACKSLASH_ESCAPE=-fbackslash
SYS_FOPT_DISABLE_BACKSLASH_ESCAPE=-fno-backslash
endif

GCC_FC_GFORTRAN?=gfortran
GCC_FC_G95?=g95
GCC_CC?=gcc
GCC_CXX?=g++

SYS_FC=$(GCC_FC_$(_USE_GNU_FORTRAN))
SYS_CC=$(GCC_CC)
SYS_CXX=$(GCC_CXX)

SYS_FOPT_DYNL_PIC=   -fPIC
SYS_COPT_DYNL_PIC=   -fPIC
SYS_CXXOPT_DYNL_PIC= -fPIC

GCC_FOPT_GFORTRAN=	-Wall -std=gnu \
			-Wno-unused-dummy-argument $(SYS_FOPT_GFORTRAN_CHECK)
GCC_FOPT_G95=		-Wall -traditional
GCC_COPT_GCC=		-Wall -std=gnu11     -pedantic-errors
GCC_CXXOPT_GCC=		-Wall -std=c++98   -pedantic-errors

SYS_FOPT=$(GCC_FOPT_$(_USE_GNU_FORTRAN)) $(SYS_FOPT_DISABLE_BACKSLASH_ESCAPE)
SYS_COPT=$(GCC_COPT_GCC)
SYS_CXXOPT=$(GCC_CXXOPT_GCC)

SYS_FOPT+= -fno-second-underscore -fdollar-ok -fargument-alias
# o Disable appending a second underscore to externals
# o Allow dollar signs(`$') in symbol names
# o Specify that arguments may alias each other and globals

# Alignment options
# SAD assumes that all data structures are aligned by octet-bytes
# in order to access via rlist(*) REAL*8 array index.

# Apply 8bytes alignment for function entry point for setfnp_()
SYS_FOPT+= -falign-functions=8

# Apply 16bytes alignment for double on IA-32 architecture.
# Note:
# Natural alignment of Pentium is 4bytes.
# MacOS X for Intel requires 16bytes alignment.
ifeq ($(CPU_ARCH),i386)
#SYS_FOPT+=   -mpreferred-stack-boundary=4
#SYS_COPT+=   -mpreferred-stack-boundary=4
#SYS_CXXOPT+= -mpreferred-stack-boundary=4
endif

# Use GCC plugin
ifneq ($(USE_GCC_PLUGIN),)
ifeq ($(shell test $(_USE_GCC) -ge 4005000 ; echo $$?), 0)
_USE_GCC_PLUGIN=$(shell $(SADSRCDIR)/bin/realpath $(USE_GCC_PLUGIN))
SYS_FOPT+=   -fplugin=$(_USE_GCC_PLUGIN)
SYS_COPT+=   -fplugin=$(_USE_GCC_PLUGIN)
SYS_CXXOPT+= -fplugin=$(_USE_GCC_PLUGIN)
endif
endif

SYS_LDOPT_RUNTIME_PATH_OPT=-Wl,-rpath,
SYS_LDOPT_EXPORT_SYMBOL=-Wl,--export-dynamic

# Object specific options
SYS_COPT_gdtoa/gdtoa=-fno-strict-aliasing
SYS_COPT_gdtoa/misc=-fno-strict-aliasing
SYS_COPT_tfProcess_=-fno-strict-aliasing
SYS_FOPT_tffsmatch=-fno-tree-vectorize

# Disable range-check for using unsigned integer literal[not portable]
SYS_FOPT_tfstk=-fno-range-check

# Disable range-check for using unsigned integer literal sum[not portable]
SYS_FOPT+=-fno-range-check

# Default debug build options
LIB_COPT=	-g -O1
LIB_CXXOPT=	-g -O1
FOPT=		-g -O1
COPT=		-g -O1
CXXOPT=		-g -O1

# End of File
