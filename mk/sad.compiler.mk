# New Makefile -- Compiler Directive Part
#

# Check CC/CXX/FC Variable
ifneq ($(origin CC),default)
USER_CC:=$(CC)
_USER_CC:=$(notdir $(CC))
endif

ifneq ($(origin CXX),default)
USER_CXX:=$(CXX)
_USER_CXX:=$(notdir $(CXX))
endif

ifneq ($(origin FC),default)
USER_FC:=$(FC)
_USER_FC:=$(notdir $(FC))
endif

# Guess COMPILER type by user specified compiler name...
ifdef _USER_FC

ifeq ($(shell echo $(_USER_FC) | grep ^g95 >/dev/null; echo $$?),0)
ifeq ($(shell echo $(_USER_FC) | grep '^g95-[0-9][0-9]*$$' >/dev/null; echo $$?),0)
_FC_VER=$(shell echo $(_USER_FC) | sed -e 's/.*-\([0-9]\)/\1./')
endif
ifeq ($(shell echo $(_USER_FC) | grep '^g95.*-[0-9][0-9]*\.[0-9.]*$$' >/dev/null; echo $$?),0)
_FC_VER=$(shell echo $(_USER_FC) | sed -e 's/^.*-//')
endif

COMPILER?=GCC
ifneq ($(_FC_VER),)
USE_G95?=$(_FC_VER)
else
USE_G95?=0.91
endif
endif # ^g95

ifeq ($(shell echo $(_USER_FC) | grep ^g77 >/dev/null; echo $$?),0)
ifeq ($(shell echo $(_USER_FC) | grep '^g77-[0-9][0-9]*$$' >/dev/null; echo $$?),0)
_FC_VER=$(shell echo $(_USER_FC) | sed -e 's/.*-\([0-9]\)/\1./')
endif
ifeq ($(shell echo $(_USER_FC) | grep '^g77.*-[0-9][0-9]*\.[0-9.]*$$' >/dev/null; echo $$?),0)
_FC_VER=$(shell echo $(_USER_FC) | sed -e 's/^.*-//')
endif

COMPILER?=GCC
ifneq ($(_FC_VER),)
USE_GCC?=$(_FC_VER)
else
USE_GCC?=3.4.6
endif
endif # ^g77

ifeq ($(shell echo $(_USER_FC) | grep ^gfortran >/dev/null; echo $$?),0)
ifeq ($(shell echo $(_USER_FC) | grep '^gfortran[0-9][0-9]*$$' >/dev/null; echo $$?),0)
_FC_VER=$(shell echo $(_USER_FC) | sed -e 's/^gfortran\([0-9]\)/\1./')
endif
ifeq ($(shell echo $(_USER_FC) | grep '^gfortran.*-[0-9][0-9]*\.[0-9.]*$$' >/dev/null; echo $$?),0)
_FC_VER=$(shell echo $(_USER_FC) | sed -e 's/^.*-//')
endif

COMPILER?=GCC
ifneq ($(_FC_VER),)
USE_GCC?=$(_FC_VER)
else
USE_GCC?=4.3.0
endif
endif # ^gfortran

ifeq ($(shell echo $(_USER_FC) | grep ^ifc >/dev/null; echo $$?),0)
_FC_VER=$(shell $(USER_FC) -V 2>&1 | grep Version | sed -e 's/.*Version[ 	]*\([0-9][0-9.]*\).*/\1/')

COMPILER?=ICC
ifneq ($(_FC_VER),)
USE_IFC?=$(_FC_VER)
endif
endif # ^ifc

ifeq ($(shell echo $(_USER_FC) | grep ^ifort >/dev/null; echo $$?),0)
_FC_VER=$(shell $(USER_FC) -V 2>&1 | grep Version | sed -e 's/.*Version[ 	]*\([0-9][0-9.]*\).*/\1/')

COMPILER?=ICC
ifneq ($(_FC_VER),)
USE_IFC?=$(_FC_VER)
endif
endif # ^ifort

endif # _USER_FC

ifdef _USER_CC

ifeq ($(shell echo $(_USER_CC) | grep ^gcc >/dev/null; echo $$?),0)
ifeq ($(shell echo $(_USER_CC) | grep '^gcc[0-9][0-9]*$$' >/dev/null; echo $$?),0)
_CC_VER=$(shell echo $(_USER_CC) | sed -e 's/^gcc\([0-9]\)/\1./')
endif
ifeq ($(shell echo $(_USER_CC) | grep '^gcc.*-[0-9][0-9]*\.[0-9.]*$$' >/dev/null; echo $$?),0)
_CC_VER=$(shell echo $(_USER_CC) | sed -e 's/^.*-//')
endif

COMPILER?=GCC
ifneq ($(_CC_VER),)
USE_GCC?=$(_CC_VER)
endif
endif # ^gcc

ifeq ($(shell echo $(_USER_CC) | grep ^icc >/dev/null; echo $$?),0)
_CC_VER=$(shell $(USER_CC) -V 2>&1 | grep Version | sed -e 's/.*Version[ 	]*\([0-9][0-9.]*\).*/\1/')

COMPILER?=ICC
ifneq ($(_CC_VER),)
USE_ICC?=$(_CC_VER)
endif
endif # ^icc

endif # _USER_CC

ifdef _USER_CXX

ifeq ($(shell echo $(_USER_CXX) | grep ^g++ >/dev/null; echo $$?),0)
ifeq ($(shell echo $(_USER_CXX) | grep '^g++[0-9][0-9]*$$' >/dev/null; echo $$?),0)
_CXX_VER=$(shell echo $(_USER_CXX) | sed -e 's/^g++\([0-9]\)/\1./')
endif
ifeq ($(shell echo $(_USER_CXX) | grep '^g++.*-[0-9][0-9]*\.[0-9.]*$$' >/dev/null; echo $$?),0)
_CXX_VER=$(shell echo $(_USER_CXX) | sed -e 's/^.*-//')
endif

COMPILER?=GCC
ifneq ($(_CXX_VER),)
USE_GCC?=$(_CXX_VER)
endif
endif # ^g++

ifeq ($(shell echo $(_USER_CXX) | grep ^icpc >/dev/null; echo $$?),0)
_CXX_VER=$(shell $(USER_CXX) -V 2>&1 | grep Version | sed -e 's/.*Version[ 	]*\([0-9][0-9.]*\).*/\1/')

COMPILER?=ICC
ifneq ($(_CXX_VER),)
USE_ICC?=$(_CXX_VER)
endif
endif # ^icpc

endif # _USER_CXX

ifndef COMPILER
COMPILER=System
endif

# Convert COMPILER symbol
ifeq ($(COMPILER),GNU)
override COMPILER=GCC
endif

ifeq ($(COMPILER),Intel)
override COMPILER=ICC
endif

ifeq ($(COMPILER),Portland)
override COMPILER=PGI
endif

ifeq ($(COMPILER),SUN)
override COMPILER=SunStudio
endif

# Set USE_*CC directive
ifneq ($(COMPILER),System)
KNOWN_COMPILER=NO

ifeq ($(COMPILER),GCC)
KNOWN_COMPILER=YES
endif

ifeq ($(COMPILER),ICC)
KNOWN_COMPILER=YES
endif

ifeq ($(COMPILER),PGI)
KNOWN_COMPILER=YES
endif

ifeq ($(COMPILER),NAG)
KNOWN_COMPILER=YES
endif

ifeq ($(COMPILER),SunStudio)
KNOWN_COMPILER=YES
endif

ifeq ($(KNOWN_COMPILER),NO)
$(warning Unkown compiler type is specified!)
endif

endif # COMPILER!=System

# Preserve $COMPILER into _COMPILER
_COMPILER:=$(COMPILER)

# End of File
