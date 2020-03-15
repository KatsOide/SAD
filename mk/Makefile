# New MakeFile -- Driver Part
#
SHELL=/bin/sh

include ./config.mk

CONF_FILE?=sad.conf

BINDIR=$(SADDIR)/bin
SRCDIR=$(SADDIR)/src
PKGDIR=$(SADDIR)/Packages
SADCONF=$(SADDIR)/config
SADMK=$(SADDIR)/mk
FILES=$(SADDIR)/files
CONTRIB=$(SADDIR)/contrib
OBJDIR=$(SADDIR)/obj/$(MACH_ARCH)

PRE_BUILD= setup-gdtoa
PRE_DEPEND=setup-gdtoa

all:

-include $(SADDIR)/$(CONF_FILE)

include $(SADMK)/sad.compiler.mk
# -- Begin initialize options --
SYS_OBJS=
SYS_LIBS=
# -- End   initialize options --

-include $(SADCONF)/$(OS_NAME).spec
-include $(SADDIR)/$(CONF_FILE)

include $(SADMK)/sad.repository.mk
include $(SADMK)/sad.port.mk
include $(SADMK)/sad.doc.mk
include $(SADMK)/sad.install.mk
include $(SADMK)/sad.depend.mk

VPATH=$(SRCDIR)

## Message Output Rule
_MSG_DEPRECATED=cat /dev/null
ifdef CONFIG_DEPRECATED
_MSG_EXPIRATION_DATE=
ifneq ($(CONFIG_EXPIRATION_DATE),)
_MSG_EXPIRATION_DATE=It it scheduled to be removed on or after $(CONFIG_EXPIRATION_DATE).
endif
_MSG_RECOMMEND=
ifneq ($(CONFIG_RECOMMEND),)
_MSG_RECOMMEND=Recommendation: $(CONFIG_RECOMMEND)
endif
_MSG_DEPRECATED=( \
	echo ""; \
	echo "Current configuration is deprecated"; \
	echo "You may wish to reconsider installing it:"; \
	echo "$(CONFIG_DEPRECATED)"; \
	[ ! -z "$(_MSG_EXPIRATION_DATE)" ] && echo "$(_MSG_EXPIRATION_DATE)"; \
	[ ! -z "$(_MSG_RECOMMEND)" ] && echo "$(_MSG_RECOMMEND)"; \
	echo ""; \
	sleep 5)
endif
ifdef DEPRECATED
_MSG_EXPIRATION_DATE=
ifneq ($(EXPIRATION_DATE),)
_MSG_EXPIRATION_DATE=It it scheduled to be removed on or after $(EXPIRATION_DATE).
endif
_MSG_RECOMMEND=
ifneq ($(RECOMMEND),)
_MSG_RECOMMEND=Recommendation: $(RECOMMEND)
endif
_MSG_DEPRECATED=( \
	echo ""; \
	echo "This architecture[$(MACH_ARCH)] is deprecated!"; \
	echo "You may wish to reconsider installing it:"; \
	echo "$(DEPRECATED)"; \
	[ ! -z "$(_MSG_EXPIRATION_DATE)" ] && echo "$(_MSG_EXPIRATION_DATE)"; \
	[ ! -z "$(_MSG_RECOMMEND)" ] && echo "$(_MSG_RECOMMEND)"; \
	echo ""; \
	sleep 5)
endif

_MSG_BROKEN=cat /dev/null
ifdef CONFIG_UNSUPPORTED
_MSG_RECOMMEND=
ifneq ($(CONFIG_RECOMMEND),)
_MSG_RECOMMEND=Recommendation: $(CONFIG_RECOMMEND)
endif
_MSG_BROKEN=( \
	echo ""; \
	echo "Current configuration is marked as UNSUPPORTED,"; \
	echo "because $(CONFIG_UNSUPPORTED)."; \
	[ ! -z "$(_MSG_RECOMMEND)" ] && echo "$(_MSG_RECOMMEND)"; \
	echo ""; \
	test 0 = 1)
endif
ifdef BROKEN
_MSG_BROKEN=( \
	echo ""; \
	echo "This architecture[$(MACH_ARCH)] is marked as BROKEN,"; \
	echo "because $(BROKEN)."; \
	echo ""; \
	test 0 = 1)
endif

check-port:
	@$(_MSG_DEPRECATED)
	@$(_MSG_BROKEN)

## Main Target Rule
ifeq ($(USE_LIBTAI),YES)
PRE_PATCH+= libtai-patch
PRE_BUILD+= $(LIBTAI_INSTALL_TARGET)
endif

ifeq ($(USE_LIBFFI),YES)
PRE_PATCH+= libffi-patch
PRE_BUILD+= $(LIBFFI_INSTALL_TARGET)
endif

ifeq ($(USE_TCLTK),YES)
PRE_PATCH+= tcltk-patch
PRE_BUILD+= $(TCLTK_INSTALL_TARGET)
endif

all:	check-port $(PRE_BUILD) exe doc

setup-gdtoa:	setup-gdtoa-src setup-gdtoa-gen

setup-gdtoa-src:	$(CONTRIB)/gdtoa
	[ -d $(SRCDIR)/gdtoa ] || mkdir -p $(SRCDIR)/gdtoa
	(cd $(SRCDIR)/gdtoa; $(BINDIR)/lndir $(CONTRIB)/gdtoa)

setup-gdtoa-gen:	arith.h gd_qnan.h

clean-gdtoa:
	@rm -f  gd_qnan.h arith.h

mostlyclean-gdtoa:
	@rm -fr $(SRCDIR)/gdtoa

patch:	$(PRE_PATCH)

clean:	clean-gdtoa clean-std-modules
	-cp $(PKGDIR)/init.local.n $(PKGDIR)/init.local.n.bak
	@rm -f  buildinfo-db.c tfDefFuncs_.c
	@rm -f  $(OBJ_LIBSAD)        $(OBJSIM)        $(OBJGLUE)        MAIN.o
	@rm -f  $(OBJ_LIBSAD:.o=.il) $(OBJSIM:.o=.il) $(OBJGLUE:.o=.il) MAIN.il
	@rm -f  *.mod
	@rm -f  $(LIBSAD) $(SADEXE)
#	@rm -f  $(BINDIR)/gs $(PKGDIR)/init.local.n
	@rm -f  $(SADDIR)/script/design_example.sad.result
#	@rm -rf $(SADDIR)/Documents

distclean:	clean distclean-std-modules
	@rm -rf $(OBJDIR)

mostlyclean:	distclean mostlyclean-gdtoa mostlyclean-std-modules
	@rm -fr $(MOSTLYCLEAN_DIRS)

clean-std-modules:
	@for mod in /dev/null $(STD_MODULES); do \
	    if [ -x $(SADDIR)/extensions/Standard/$${mod} ]; then \
		( cd $(SADDIR)/extensions/Standard/$${mod}; \
		 unset OS_NAME MAKEFLAGS; \
		 env SADSRCDIR=$(SADDIR) $(MAKE) clean); \
	    fi; \
	done

distclean-std-modules:
	@for mod in /dev/null $(STD_MODULES); do \
	    if [ -x $(SADDIR)/extensions/Standard/$${mod} ]; then \
		(cd $(SADDIR)/extensions/Standard/$${mod}; \
		 unset OS_NAME MAKEFLAGS; \
		 env SADSRCDIR=$(SADDIR) $(MAKE) distclean); \
	    fi; \
	done

mostlyclean-std-modules:
	@for mod in /dev/null $(STD_MODULES); do \
	    if [ -x $(SADDIR)/extensions/Standard/$${mod} ]; then \
		(cd $(SADDIR)/extensions/Standard/$${mod}; \
		 unset OS_NAME MAKEFLAGS; \
		 env SADSRCDIR=$(SADDIR) $(MAKE) mostlyclean); \
	    fi; \
	done

exe:	check-port $(SADEXE) std-modules $(PKGDIR)/init.local.n $(BINDIR)/gs

lib:	check-port $(LIBSAD)

std-modules:
	@for mod in /dev/null $(STD_MODULES); do \
	    if [ -x $(SADDIR)/extensions/Standard/$${mod} ]; then \
		(cd $(SADDIR)/extensions/Standard/$${mod}; \
		 unset OS_NAME MAKEFLAGS; \
		 env SADSRCDIR=$(SADDIR) $(MAKE) module); \
	    fi; \
	done

$(LIBSAD):	$(OBJ_LIBSAD)
	@rm -f $(LIBSAD)
	@obj_libsad=`echo $(OBJ_LIBSAD) | tr ' ' '\n'| sort -u`; \
	    $(AR) r $(LIBSAD) $${obj_libsad}
	sleep 2
	$(RANLIB) $(LIBSAD)

_OBJ_LIBSAD=$(LIBSAD)
__OBJ_LIBSAD=$(LIBSAD)
ifdef USE_LTO
_OBJ_LIBSAD=$(OBJ_LIBSAD)
__OBJ_LIBSAD=$(shell lorder $(OBJ_LIBSAD) | tsort 2>/dev/null)
endif

$(SADEXE): $(PRE_BUILD) $(_OBJ_LIBSAD) $(OBJSIM) $(OBJGLUE) MAIN.o
	$(FC) $(SYS_FC_ABIOPT) -o $@ $(FC_FLAGS) MAIN.o $(SYS_EXTSIM) $(OBJSIM) $(OBJGLUE) $(__OBJ_LIBSAD) $(LD_FLAGS)

$(PKGDIR)/init.local.n:	$(SRCDIR)/init.local.n.in
	-cp $@ $@.bak
	sed -e "s,xSADROOT,$(SADDIR),g" \
	    -e "s,xSAD_SHARE_ROOT,$(SADDIR),g" \
	    -e "s,xScriptDirectory,script,g" \
		$< >$@

$(BINDIR)/gs:	$(SRCDIR)/gs.in
	sed -e "s,xSADROOT,$(SADDIR),g" \
	    -e "s,xSAD_SHARE_ROOT,$(SADDIR),g" \
	    -e "s,xSAD_ARCH_ROOT,$(SAD_ARCH_ROOT),g" \
	    -e "s,xSAD_EXE_ROOT,$(SADDIR)/obj,g" \
	    -e "s,xSADEXE,$(SADEXE),g" \
	    -e "s,xSAD_ARCH_OR_OBJ,obj,g" \
	    -e "s,xSAD_PACKAGES,Packages,g" \
	    $< >$@ && chmod 755 $@

# End of File
