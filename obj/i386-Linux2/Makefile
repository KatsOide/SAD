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

PRE_DEPEND+= setup-gdtoa

all:

include $(SADMK)/sad.depend.mk

-include $(SADDIR)/$(CONF_FILE)

include $(SADMK)/sad.compiler.mk
include $(SADCONF)/default.spec
ifneq ($(OS_TYPE),)
-include $(SADCONF)/$(OS_TYPE).spec
endif
ifneq ($(OS_NAME),$(OS_TYPE))
-include $(SADCONF)/$(OS_NAME).spec
endif

-include $(SADDIR)/$(CONF_FILE)

USE_LIBTAI=YES
USE_LIBUTIL=YES
USE_LIBCRYPT=YES

include $(SADMK)/sad.repository.mk
include $(SADMK)/sad.port.mk
include $(SADMK)/sad.doc.mk
include $(SADMK)/sad.install.mk

VPATH=$(SRCDIR)

## Main Target Rule
PRE_BUILD+= setup-gdtoa

ifeq ($(USE_LIBTAI),YES)
PRE_PATCH+= libtai-patch
PRE_BUILD+= $(LIBTAI_INSTALL_TARGET)
endif

ifeq ($(USE_TCLTK),YES)
PRE_PATCH+= tcltk-patch
PRE_BUILD+= $(TCLTK_INSTALL_TARGET)
endif

ifeq ($(USE_PYTHON),YES)
PRE_PATCH+= python-patch
PRE_BUILD+= $(PYTHON_INSTALL_TARGET)
endif

ifneq ($(HAVE_RANLIB),YES)
RANLIB="echo ranlib"
endif

all:	$(PRE_BUILD) exe doc

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
	@rm -f  MAIN_VERSION.f tfBuildInfo_.c tfDefFuncs_.c
	@rm -f  $(OBJ_LIBSAD)        $(OBJSIM)        $(OBJGLUE)        MAIN.o
	@rm -f  $(OBJ_LIBSAD:.o=.il) $(OBJSIM:.o=.il) $(OBJGLUE:.o=.il) MAIN.il
	@rm -f  $(LIBSAD) $(SADEXE)
	@rm -f  $(BINDIR)/gs $(PKGDIR)/init.local.n
	@rm -f  $(SADDIR)/script/design_example.sad.result
	@rm -rf $(SADDIR)/Documents

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

exe:	$(SADEXE) std-modules $(PKGDIR)/init.local.n $(BINDIR)/gs

lib:	$(LIBSAD)

std-modules:
	@for mod in /dev/null $(STD_MODULES); do \
	    if [ -x $(SADDIR)/extensions/Standard/$${mod} ]; then \
		(cd $(SADDIR)/extensions/Standard/$${mod}; \
		 unset OS_NAME MAKEFLAGS; \
		 env SADSRCDIR=$(SADDIR) $(MAKE) module); \
	    fi; \
	done

ifdef BROKEN
$(SADEXE) $(LIBSAD):
	@echo ""
	@echo "---- SAD port is broken! ---"
	@echo "  $(BROKEN)"
	@echo "---- SAD port is broken! ---"
	@echo ""
else
$(LIBSAD):	$(OBJ_LIBSAD)
	rm -f $(LIBSAD)
	obj_libsad=`echo $(OBJ_LIBSAD) | tr ' ' '\n'| sort -u`; \
	    $(AR) rv $(LIBSAD) $${obj_libsad}
	$(RANLIB) $(RANLIB_FLAGS) $(LIBSAD)

ifeq ($(USE_LIBSAD),YES)
$(SADEXE):	$(PRE_BUILD)     $(LIBSAD) $(OBJSIM) $(OBJGLUE) MAIN.o
	$(FC) -o $@ $(FC_FLAGS) MAIN.o $(OBJSIM) $(OBJGLUE)     $(LIBSAD) $(LD_FLAGS)
else
$(SADEXE):	$(PRE_BUILD) $(OBJ_LIBSAD) $(OBJSIM) $(OBJGLUE) MAIN.o
	$(FC) -o $@ $(FC_FLAGS) MAIN.o $(OBJSIM) $(OBJGLUE) $(OBJ_LIBSAD) $(LD_FLAGS)
endif
endif

$(PKGDIR)/init.local.n:	$(SRCDIR)/init.local.n.in
	sed -e "s,xSADROOT,$(SADDIR),g" \
	    -e "s,xSAD_SHARE_ROOT,$(SADDIR),g" \
	    -e "s,xScriptDirectory,script,g" \
	    -e "s,xMODROOT,\"$(MODULE_ROOT)\",g" \
	    -e "s,xFONT_SYSTEM,\"$(TCLTK_FONT_SYSTEM)\",g" \
		$< >$@

$(BINDIR)/gs:	$(SRCDIR)/gs.in
	sed -e "s,xSADROOT,$(SADDIR),g" \
	    -e "s,xSAD_SHARE_ROOT,$(SADDIR),g" \
	    -e "s,xSAD_ARCH_ROOT,$(SAD_ARCH_ROOT),g" \
	    -e "s,xSAD_EXE_ROOT,$(SADDIR)/obj,g" \
	    -e "s,xSADEXE,$(SADEXE),g" \
	    $< >$@ && chmod 755 $@

# End of File
