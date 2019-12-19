# Makefile Skeleton for SAD Dynamic Loading Extension/External Engine Module

SRCDIR?=.
SCRIPTDIR?=.
MODULE_TYPE?=Extension
_SYS_EXTENSION_SUFFIX=.so
_SYS_ENGINE_SUFFIX=.exe

ifeq ($(MODULE_TYPE),Extension)
MODULE_TARGET=	$(MODULE_NAME)$(_SYS_EXTENSION_SUFFIX)
endif
ifeq ($(MODULE_TYPE),Engine)
MODULE_TARGET=	$(MODULE_NAME)$(_SYS_ENGINE_SUFFIX)
endif

PATH:=$(PATH):$(SADSRCDIR)/bin

ifndef _MAKELEVEL
_MAKELEVEL=$(MAKELEVEL)
endif

_MAKE_SUB_MODULES=cat /dev/null;
ifdef SUB_MODULES
_MAKE_SUB_MODULES=env MAKELEVEL=0 SADSRCDIR=`$(SADSRCDIR)/bin/realpath $(SADSRCDIR)` sh -c 'for dir in /dev/null $(SUB_MODULES); do if [ -d $${dir} ]; then (cd $${dir}; $(MAKE) _MAKELEVEL=$(MAKELEVEL) $@); fi; done';
endif

ifndef OS_NAME
MAKE_ENV=\
	\`$(SADSRCDIR)/bin/HostArch -config\`

all module install install-module install-script:	top-gensrcfiles

mostlyclean:	top-distclean

top-gensrcfiles gensrcfiles depend mostlyclean:
	@(target=$@; \
	if [ "x$${target##top-}" = "x$${target}" ]; then \
	    $(_MAKE_SUB_MODULES) \
	fi; \
	env  MAKELEVEL=$$(( $(_MAKELEVEL) + 1 )) \
	     TOPDIR=`pwd` \
	    _SRCDIR=`$(SADSRCDIR)/bin/realpath $(SRCDIR)` \
	    _SCRIPTDIR=`$(SADSRCDIR)/bin/realpath $(SCRIPTDIR)` \
	$(SHELL) -c "$(MAKE) -f Makefile $${target##top-} $(MAKE_ENV)")

top-distclean clean distclean:
	@(target=$@; \
	if [ "x$${target##top-}" = "x$${target}" ]; then \
	    $(_MAKE_SUB_MODULES) \
	fi; \
	if [ -d obj/`$(SADSRCDIR)/bin/HostArch` ]; then \
	    env  MAKELEVEL=$$(( $(_MAKELEVEL) + 1 )) \
		 TOPDIR=`pwd` \
		_SRCDIR=`$(SADSRCDIR)/bin/realpath $(SRCDIR)` \
		_SCRIPTDIR=`$(SADSRCDIR)/bin/realpath $(SCRIPTDIR)` \
		$(SHELL) -c "cd obj/`$(SADSRCDIR)/bin/HostArch`; \
		$(MAKE) -f ../../Makefile $${target##top-} $(MAKE_ENV)"; \
	    if [ "x$${target##top-}" = x"distclean" ]; then \
		rm -fr obj/`$(SADSRCDIR)/bin/HostArch`; \
		rmdir obj 2>/dev/null; true ; \
	    fi; \
	fi)

all module debug install install-module install-script:
	@(target=$@; \
	if [ "x$${target##top-}" = "x$${target}" ]; then \
	    $(_MAKE_SUB_MODULES) \
	fi; \
	mkdir -p obj/`$(SADSRCDIR)/bin/HostArch`; \
	env  MAKELEVEL=$$(( $(_MAKELEVEL) + 1 )) \
	     TOPDIR=`pwd` \
	    _SRCDIR=`$(SADSRCDIR)/bin/realpath $(SRCDIR)` \
	    _SCRIPTDIR=`$(SADSRCDIR)/bin/realpath $(SCRIPTDIR)` \
	$(SHELL) -c "cd obj/`$(SADSRCDIR)/bin/HostArch`; \
	 $(MAKE) -f ../../Makefile $${target##top-} $(MAKE_ENV)")

ifdef MODULE_EXTRA_TARGET
$(MODULE_EXTRA_TARGET):
	@( \
	mkdir -p obj/`$(SADSRCDIR)/bin/HostArch`; \
	env  MAKELEVEL=$$(( $(_MAKELEVEL) + 1 )) \
	     TOPDIR=`pwd` \
	    _SRCDIR=`$(SADSRCDIR)/bin/realpath $(SRCDIR)` \
	    _SCRIPTDIR=`$(SADSRCDIR)/bin/realpath $(SCRIPTDIR)` \
	$(SHELL) -c "cd obj/`$(SADSRCDIR)/bin/HostArch`; \
	 $(MAKE) -f ../../Makefile $@ $(MAKE_ENV)")
endif

endif # !OS_NAME

ifdef OS_NAME
ifdef _SRCDIR
SRCDIR=$(_SRCDIR)
SCRIPTDIR=$(_SCRIPTDIR)
endif

VPATH=$(SRCDIR):$(SCRIPTDIR)

CONF_FILE?=sad.conf

SADDIR=$(SADSRCDIR)

BINDIR=$(SADDIR)/bin
OBJDIR=$(TOPDIR)/obj/$(MACH_ARCH)
SADCONF=$(SADDIR)/config
SADMK=$(SADDIR)/mk
FILES=$(SADDIR)/files
CONTRIB=$(TOPDIR)/contrib

SYS_SAVE_USES=X11 TCLTK EPICS LIBFFI LIBTAI
$(foreach _FUNC,$(SYS_SAVE_USES),$(eval USER_USE_$(_FUNC):=$(USE_$(_FUNC))))

-include $(SADSRCDIR)/$(CONF_FILE)
$(foreach _FUNC,$(SYS_SAVE_USES),$(eval override USE_$(_FUNC)=$(USER_USE_$(_FUNC))))

include $(SADMK)/sad.compiler.mk
# -- Begin initialize options --
SYS_OBJS= 
SYS_LIBS= 
# -- End   initialize options --

-include $(SADCONF)/$(OS_NAME).spec
-include $(SADSRCDIR)/$(CONF_FILE)
$(foreach _FUNC,$(SYS_SAVE_USES),$(eval override USE_$(_FUNC)=$(USER_USE_$(_FUNC))))

override _SYS_IN_MODULE_MK=YES
include $(SADMK)/sad.port.mk

-include $(SRCDIR)/.depend

FC_FLAGS+=$(SYS_FOPT_DYNL_PIC)
CC_FLAGS+=$(SYS_COPT_DYNL_PIC)

# for Fortran module files
FC_FLAGS+=-I$(SADSRCDIR)/obj/$(MACH_ARCH)

ifndef LD_SHARED
LD_SHARED=$(CC) $(SYS_CC_ABIOPT) -shared
endif

ifndef LD_ENGINE
LD_ENGINE=$(CC) $(SYS_CC_ABIOPT)
endif

all:	module

ifeq ($(USE_FRAMEWORK),YES)
GENINC_DIRS+=$(SRCDIR)/framework
endif

ifeq ($(MODULE_TYPE),Extension)
GENINC_DIRS+=$(SRCDIR)/inc $(SRCDIR)/sim
endif

gensrcfiles:	$(GENINC_DIRS) $(GENSRC_FILES)

MKDEP=$(SADSRCDIR)/bin/mkdep
MKDEP_OPT=--exclude=/usr/local/include --exclude=/usr/include

ifeq ($(MODULE_TYPE),Script)
depend:
else
depend:		gensrcfiles
	@rm -f $(SRCDIR)/.depend
	env PATH="$${PATH}:$(SADSRCDIR)/bin" BIN_SH=xpg4 \
	$(SHELL) -c "cd $(SRCDIR); \
	$(MKDEP) ${MKDEP_OPT}               $(FOPT_ADD)   *.f ; \
	$(MKDEP) ${MKDEP_OPT} $(OPT_INCDIR) $(COPT_ADD)   *.c ; \
	$(MKDEP) ${MKDEP_OPT} $(OPT_INCDIR) $(CXXOPT_ADD) *.cc" \
	    | sort >>$(SRCDIR)/.depend
endif

mostlyclean:
	@for f in /dev/null $(GENSRC_FILES); do \
	    if [ "x$${f}" != "x/dev/null" ]; then \
		if [ -x $${f} ]; then \
		    rm -fr $${f}; \
		else \
		    rm -f $${f}; \
		fi; \
	    fi; \
	done
	@rm -f $(SRCDIR)/inc $(SRCDIR)/sim $(MOSTLYCLEAN_FILES)
	@rm -fr $(MOSTLYCLEAN_DIRS) $(SRCDIR)/framework obj
	@rmdir contrib 2>/dev/null; true

module:		$(MODULE_TARGET)

clean:
	@rm -f $(MODULE_TARGET)
	@rm -f $(OBJS)

distclean:	clean

$(SRCDIR)/inc:	$(SADSRCDIR)/src/inc
	@rm -f $@
	@ln -s $(SADSRCDIR)/src/inc $@

$(SRCDIR)/sim:	$(SADSRCDIR)/src/sim
	@rm -f $@
	@ln -s $(SADSRCDIR)/src/sim $@

$(SRCDIR)/framework:
	@mkdir -p $(SRCDIR)/framework
	@(cd $(SRCDIR)/framework; \
	for f in /dev/null $(FRAMEWORK_HEADERS); do \
	    if [ "x$${f}" != "x/dev/null" ]; then \
		ln -s $(SADSRCDIR)/src/$${f} $${f}; \
	    fi; \
	done)

ifeq ($(MODULE_TYPE),Extension)
$(MODULE_NAME)$(_SYS_EXTENSION_SUFFIX):	$(MODULE_PRE_BUILD_TARGET) $(OBJS)
	$(LD_SHARED) -o $@ $(OBJS) $(LDOPT_ADD)
endif

ifeq ($(MODULE_TYPE),Engine)
$(MODULE_NAME)$(_SYS_ENGINE_SUFFIX):	$(MODULE_PRE_BUILD_TARGET) $(OBJS)
	$(LD_ENGINE) -o $@ $(OBJS) $(LD_FLAGS) $(LDOPT_ADD)
endif

ifdef MODULE_ROOT
MODULE_DIR=$(MODULE_ROOT)
ifdef MODULE_SUBDIR
MODULE_DIR=$(MODULE_ROOT)/$(MODULE_SUBDIR)
endif

install:	install-module install-script

ifdef MODULE_SCRIPT
install-script:	$(MODULE_SCRIPT)
	@mkdir -p $(MODULE_DIR)
	@printf '%s' "Installing script files:"
	@for i in /dev/null $(MODULE_SCRIPT); do \
	    if [ "x$${i}" != "x/dev/null" ]; then \
		printf ' %s' $${i}; \
		$(INSTALL_DATA) $(SCRIPTDIR)/$${i} $(MODULE_DIR)/$${i}; \
	    fi; \
	done
	@printf ' ...done\n'
else
install-script:
endif

ifeq ($(MODULE_TYPE),Script)
install-module:
endif

ifeq ($(MODULE_TYPE),Extension)
install-module: $(MODULE_NAME)$(_SYS_EXTENSION_SUFFIX)
	@mkdir -p $(MODULE_DIR)/$(MACH_ARCH)
	@printf '%s' "Installing extension module:"
	@for f in $(MODULE_NAME)$(_SYS_EXTENSION_SUFFIX); do \
	    printf ' %s' $$f; \
	    $(INSTALL_DATA) $$f $(MODULE_DIR)/$(MACH_ARCH)/$$i; \
	done
	@printf ' ...done\n'
endif

ifeq ($(MODULE_TYPE),Engine)
install-module: $(MODULE_NAME)$(_SYS_ENGINE_SUFFIX)
	@mkdir -p $(MODULE_DIR)/$(MACH_ARCH)
	@printf '%s' "Installing external engine:"
	@for f in $(MODULE_NAME)$(_SYS_ENGINE_SUFFIX); do \
	    printf ' %s' $$f; \
	    $(INSTALL_DATA) $$f $(MODULE_DIR)/$(MACH_ARCH)/$$i; \
	done
	@printf ' ...done\n'
endif

endif # MODULE_ROOT

ifndef MODULE_ROOT
install install-module install-script:
	@echo "Edit MODULE_ROOT in Makefile or Type $(MAKE) install with \`MODULE_ROOT=(Top directory of extension module library)'"
endif # !MODULE_ROOT

endif # OS_NAME

# End of File
