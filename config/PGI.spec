# Generic Object Makefile Compiler Support Part -- PGI
#

# _USE_* <- major * 10^6 + minor * 10^3 + patch-level[minor,patch-level < 10^3]
ifdef USE_PGI
_USE_PGI=$(shell echo "$(USE_PGI)" | sed -e 's/ *//g' -e 's/[^0-9.].*//' -e 's/\.*$$/.0.0.0./' -e 's/^\.*//' -e 's/^\([^.][^.]*\)\.\([^.][^.]*\)\.\([^.][^.]*\)\..*/\1.000\2.000\3/' -e 's/\.[^.]*\([^.][^.][^.]\)/\1/g' -e 's/^0*\(.\)/\1/g')
else
_USE_PGI=0
endif

# Default PGI version: 6.000.000
ifeq ($(_USE_PGI),0)
_USE_PGI=6000000
endif

HAVE_F_FORMAT_Q_EDIT=YES
HAVE_F_FGETC=YES

SYS_FC=pgf95
SYS_CC=pgcc
SYS_CXX=pgCC

SYS_FOPT_DYNL_PIC=-fPIC
SYS_COPT_DYNL_PIC=-fPIC
SYS_CXXOPT_DYNL_PIC=-fPIC

SYS_FOPT=-fPIC -Mrecursive -Mnosecond_underscore
# o Enable local auto variable
# o Disable appending a second underscore to externals

SYS_COPT=-fPIC

SYS_FOPT_ENABLE_BACKSLASH_ESCAPE=#-Mnobackslash
SYS_FOPT_DISABLE_BACKSLASH_ESCAPE=-Mbackslash

SYS_FOPT+=$(SYS_FOPT_DISABLE_BACKSLASH_ESCAPE)

# Default debug build options
LIB_COPT=	-g
LIB_CXXOPT=	-g
FOPT=		-g
COPT=		-g
CXXOPT=		-g

# End of File
