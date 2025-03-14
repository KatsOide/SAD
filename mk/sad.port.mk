# New Makefile -- Porting Rule
#
ifeq ($(BEFOREPORTMK),YES)

ifeq ($(DEBUG),PORT)
debug:
	@echo "STD_MODULES=$(STD_MODULES)"
	@echo "MOSTLYCLEAN_DIRS=$(MOSTLYCLEAN_DIRS)"
	@echo "OPT_RLIBDIR=$(OPT_RLIBDIR)"
	@echo "ROPT_ADD=$(ROPT_ADD)"
	@echo "SYS_LDOPT_RUNTIME_PATH=$(SYS_LDOPT_RUNTIME_PATH)"
endif

## SAD Environment
ifndef SOURCE_ARCHIVE_DIR
SOURCE_ARCHIVE_DIR=$(SADDIR)/distfiles
endif

ifdef SAD_ROOT
SAD_SHARE_ROOT=$(SAD_ROOT)/share
else
SAD_ROOT=$(SADDIR)
SAD_SHARE_ROOT=$(SAD_ROOT)/arch
endif

SAD_ARCH_ROOT=$(SAD_ROOT)/arch
SAD_ARCH_DIR=$(SAD_ARCH_ROOT)/$(MACH_ARCH)

ifndef MODULE_ROOT
MODULE_ROOT=$(SAD_SHARE_ROOT)/Extension
endif

SADEXE=sad1.exe
#SADEXE=sad11.exe
LIBSAD=libsad.a

## Check Compiler Setup
ifeq ($(SETUP_OSTYPE_DEP_COMPILER),NO)
$(error OS Type depend compiler setup is failed!)
endif
ifeq ($(SETUP_OSNAME_DEP_COMPILER),NO)
$(error OS Name depend compiler setup is failed!)
endif

## Check OBJDYNL
SYS_USE_DYNL_MODULE=NO
ifdef OBJDYNL
ifneq ($(OBJDYNL),sim/dynl-dummy.o)
SYS_USE_DYNL_MODULE=YES
endif
endif

## Compiler Environment
ifdef USER_CC
override CC=$(USER_CC)
else
ifdef SYS_CC
override CC=$(SYS_CC)
endif
endif

ifdef USER_CXX
override CXX=$(USER_CXX)
else
ifdef SYS_CXX
override CXX=$(SYS_CXX)
endif
endif

ifdef USER_FC
override FC=$(USER_FC)
else
ifdef SYS_FC
override FC=$(SYS_FC)
endif
endif

## System Environment
ifndef PATCH
PATCH=patch
endif

ifndef AR
AR=ar
endif

ifndef RANLIB
RANLIB=ranlib
endif

ifndef FETCH
FETCH=wget --passive-ftp
endif

ifndef INSTALL
INSTALL=install -c
endif

INSTALL_PROG=$(INSTALL) -m 755
INSTALL_DATA=$(INSTALL) -m 644

ifndef SHLIB_SUFFIX
SHLIB_SUFFIX=.so
endif

include $(SADMK)/sad.builtin.mk
include $(SADMK)/sad.fortran.mk

include $(SADMK)/sad.epics.mk
include $(SADMK)/sad.tcltk.mk
include $(SADMK)/sad.x11.mk
include $(SADMK)/sad.libtai.mk
include $(SADMK)/sad.libffi.mk
include $(SADMK)/sad.framework.mk

# Add Standard dynamic loadable modules
ifeq ($(SYS_USE_DYNL_MODULE),YES)
##STD_MODULES+=Parallel
##STD_MODULES+=SpaceCharge/Scheff
endif

endif # BEFOREPORTMK

ifeq ($(AFTERPORTMK),YES)

include $(SADMK)/sad.builtin.mk

include $(SADMK)/sad.tcltk.mk
include $(SADMK)/sad.libtai.mk
include $(SADMK)/sad.libffi.mk
include $(SADMK)/sad.framework.mk

ifneq ($(_SYS_IN_MODULE_MK),YES)
include $(SADMK)/sad.glue.mk
include $(SADMK)/sad.obj.mk
endif

## Compiler Flags for LTO
_FOPT_LTO=$(USE_LTO)
ifdef FOPT_LTO
_FOPT_LTO=$(FOPT_LTO)
endif
_COPT_LTO=$(USE_LTO)
ifdef COPT_LTO
_COPT_LTO=$(COPT_LTO)
endif
_CXXOPT_LTO=$(USE_LTO)
ifdef CXXOPT_LTO
_CXXOPT_LTO=$(CXXOPT_LTO)
endif

## Compiler Flags
FC_FLAGS=$(DOPT)  $(FOPT)   $(FOPT_ADD)   $(_FOPT_LTO)   $(SYS_FOPT) \
	-I$(SRCDIR) $(OPT_INCDIR) -I$(OBJDIR)
CC_FLAGS=$(DOPT)  $(COPT)   $(COPT_ADD)   $(_COPT_LTO)   $(SYS_COPT) \
	$(OPT_COPT) \
	$(SYS_PTHREAD_COPT)   \
	-I$(SRCDIR) $(OPT_INCDIR) -I$(OBJDIR) -I$(SRCDIR)/gdtoa
CXX_FLAGS=$(DOPT) $(CXXOPT) $(CXXOPT_ADD) $(_CXXOPT_LTO) $(SYS_CXXOPT) \
	$(SYS_PTHREAD_CXXOPT) \
	-I$(SRCDIR) $(OPT_INCDIR) -I$(OBJDIR)

ifndef SYS_LDOPT_RUNTIME_PATH_SEP
SYS_LDOPT_RUNTIME_PATH_SEP=:
endif

SYS_LDOPT_RUNTIME_PATH=
ifneq ($(SYS_LDOPT_RUNTIME_PATH_OPT),)
SYS_LDOPT_RUNTIME_PATH=$(SYS_LDOPT_RUNTIME_PATH_OPT)$(shell echo "$(foreach _DIR,$(OPT_RLIBDIR) $(ROPT_ADD),$(_DIR))" | sed -e 's@  *-L@$(SYS_LDOPT_RUNTIME_PATH_SEP)@g' -e 's/^-L//')
endif

LD_FLAGS=$(OPT_LIBDIR) $(SYS_LDOPT_RUNTIME_PATH) \
	$(OPT_LIBS) $(SYS_LIBS) $(LIBS) \
	$(SYS_PTHREAD_LIBS) $(SYS_LDOPT_EXPORT_SYMBOL) $(SYS_LDOPT)\
	$(LDOPT)

## Common Object Rule
.f.o: $(@:.o=.f)
	$(FC) $(SYS_FC_ABIOPT) -o $@ -c $(FC_FLAGS) $(SYS_COMMON_FOPT) $(SYS_FOPT_$(@:.o=)) $(SRCDIR)/$(@:.o=.f)

.c.o: $(@:.o=.c)
	$(CC) $(SYS_CC_ABIOPT) -o $@ -c $(CC_FLAGS) $(SYS_COMMON_COPT) $(SYS_COPT_$(@:.o=)) $(SRCDIR)/$(@:.o=.c)

.cc.o: $(@:.o=.cc)
	$(CXX) $(SYS_CXX_ABIOPT) -o $@ -c $(CXX_FLAGS) $(SYS_COMMON_COPT) $(SYS_CXXOPT_$(@:.o=)) $(SRCDIR)/$(@:.o=.cc)

## Dynamic Generated Objects
arith.h:	gdtoa/arithchk.c
	$(CC) $(SYS_CC_ABIOPT) -o arithchk $(CC_FLAGS) $< || $(CC) $(SYS_CC_ABIOPT) -o arithchk -DNO_LONG_LONG $(CC_FLAGS) $<
	./arithchk >$@
	rm -f arithchk arithchk.exe arithchk.o

gd_qnan.h:	gdtoa/qnan.c arith.h 
	$(CC) $(SYS_CC_ABIOPT) -o qnan $(CC_FLAGS) $<
	./qnan >$@
	rm -f qnan qnan.exe qnan.o

tfDefFuncs_.o:	genDefFuncs
	$(CC) $(SYS_CC_ABIOPT) -o $@ -c $(CC_FLAGS) $(SYS_COMMON_COPT) $(SYS_COPT_$(@:.o=)) $(@:.o=.c)

buildinfo-db.o:	genBuildInfo
	$(CC) $(SYS_CC_ABIOPT) -o $@ -c $(CC_FLAGS) $(SYS_COMMON_COPT) $(SYS_COPT_$(@:.o=)) $(@:.o=.c)

genDefFuncs:
	@rm -f tfDefFuncs_.c
	@date "+/* tfDefFunc_.c generated at %Y-%m-%d %H:%M:%S %z */" >>tfDefFuncs_.c
	@for func in $(FRAMEWORK_INIT_HOOKS); do \
	    echo "extern int init_framework_$${func}(void);" >>tfDefFuncs_.c; \
	done
	@for func in $(_SAD_FUNC_DEFS); do \
	    echo "extern int sadDefFunc_$${func}(void);"     >>tfDefFuncs_.c; \
	done
	@echo ""					>>tfDefFuncs_.c
	@echo "int tfdeffuncs_(void) {"			>>tfDefFuncs_.c
	@echo "/* Call init_frameworks */"		>>tfDefFuncs_.c
	@for func in $(FRAMEWORK_INIT_HOOKS); do \
	    echo "  init_framework_$${func}();"		>>tfDefFuncs_.c; \
	done
	@echo ""					>>tfDefFuncs_.c
	@echo "/* Call sadDefFuncs */"			>>tfDefFuncs_.c
	@for func in $(_SAD_FUNC_DEFS); do \
	    echo "  sadDefFunc_$${func}();"		>>tfDefFuncs_.c; \
	done
	@echo ""					>>tfDefFuncs_.c
	@echo "  return 0;"				>>tfDefFuncs_.c
	@echo "}"					>>tfDefFuncs_.c
	@echo ""					>>tfDefFuncs_.c
	@echo "/* End of File */"			>>tfDefFuncs_.c

genBuildInfo:
	@echo "Updating Build Information"; \
	release=`uname -r 2>/dev/null`; \
	user=`whoami`; hostname=`hostname`; objdir=`pwd`; \
	today=`date "+%Y-%m-%d %H:%M:%S %z"`; \
	$(BINDIR)/buildinfo \
	-pre \
	-str Repository:Owner    "$(REPOSITORY_OWNER)"     \
	-str Repository:Type     "$(REPOSITORY_TYPE)"      \
	-str Repository:Location "$(REPOSITORY_LOCATION)"  \
	-str Source:TreeName     "$(SOURCE_TREE_PET_NAME)" \
	-str Config:ABI    	 "${ABI}"  \
	-str Config:FC    	 "${FC}"  \
	-str Config:CC           "${CC}"  \
	-str Config:CXX          "${CXX}" \
	-str Config:SYS_FOPT     "$(subst ",\",$(subst \,\\\\,${SYS_FOPT}))"   \
	-str Config:SYS_COPT     "$(subst ",\",$(subst \,\\\\,${SYS_COPT}))"   \
	-str Config:SYS_CXXOPT   "$(subst ",\",$(subst \,\\\\,${SYS_CXXOPT}))" \
	-str Config:SYS_CC_ABIOPT     "$(subst ",\",$(subst \,\\\\,${SYS_FC_ABIOPT}))"   \
	-str Config:SYS_FC_ABIOPT     "$(subst ",\",$(subst \,\\\\,${SYS_CC_ABIOPT}))"   \
	-str Config:SYS_CXX_ABIOPT   "$(subst ",\",$(subst \,\\\\,${SYS_CXX_ABIOPT}))" \
	-str Config:FOPT         "$(subst ",\",$(subst \,\\\\,${FOPT}))"   \
	-str Config:COPT         "$(subst ",\",$(subst \,\\\\,${COPT}))"   \
	-str Config:CXXOPT       "$(subst ",\",$(subst \,\\\\,${CXXOPT}))" \
	-str Target:OS_TYPE	"$(OS_TYPE)" \
	-str Target:OS_NAME	"$(OS_NAME)" \
	-int Target:OS_MAJOR	 $(OS_MAJOR_VERSION) \
	-int Target:OS_MINOR	 $(OS_MINOR_VERSION) \
	-str Target:OS_RELEASE	"$${release}"  \
	-str Target:CPU_ARCH	"$(CPU_ARCH)"  \
	-str Target:MACH_ARCH	"$(MACH_ARCH)" \
	-str Target:SAD_ROOT		"$(SAD_ROOT)" \
	-str Target:SAD_ARCH_ROOT	"$(SAD_ARCH_ROOT)"  \
	-str Target:SAD_SHARE_ROOT	"$(SAD_SHARE_ROOT)" \
	-str Target:SAD_PKG_ROOT	"$(SAD_SHARE_ROOT)/Packages" \
	-str Target:SAD_MOD_ROOT	"$(MODULE_ROOT)" \
	-str Built:Location	"$${user}@$${hostname}:$${objdir}" \
	-str Built:Date		"$${today}" \
	-post >buildinfo-db.c

endif # AFTERPORTMK

ifndef BEFOREPORTMK
include $(SADMK)/sad.port.pre.mk
include $(SADMK)/sad.port.post.mk
endif
# End of File
