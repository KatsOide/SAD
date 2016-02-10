# New Makefile -- Tcl/Tk Part
#
ifeq ($(BEFOREPORTMK),YES)

ifeq ($(DEBUG),TCLTK)
debug:
	@echo "USE_TCLTK=$(USE_TCLTK)"
	@echo "BUILD_TCLTK=$(BUILD_TCLTK)"
	@echo "TCLTK_VERSION=$(TCLTK_VERSION)"
	@echo "_TCLTK_VERSION=$(_TCLTK_VERSION)"
	@echo "TCLTK_FEATURE_VERSION=$(TCLTK_FEATURE_VERSION)"
	@echo "TCLTK_GUI_BACKEND=$(TCLTK_GUI_BACKEND)"
	@echo "TCLTK_FONT_SYSTEM=$(TCLTK_FONT_SYSTEM)"
	@echo "TCLTK_CONFIG_ARG=$(TCLTK_CONFIG_ARG)"
	@echo "TCLTK_TCL_CONFIG_ARG=$(TCLTK_TCL_CONFIG_ARG)"
	@echo "TCLTK_TK_CONFIG_ARG=$(TCLTK_TK_CONFIG_ARG)"
	@echo "TCL_INSTALL_TARGETS=$(TCL_INSTALL_TARGETS)"
	@echo "TK_INSTALL_TARGETS=$(TK_INSTALL_TARGETS)"
endif

# Check USE_TCLTK
ifndef USE_TCLTK
USE_TCLTK=NO
ifeq ($(SYS_USE_DYNL_MODULE),YES)
USE_TCLTK=MODULE
ifndef USE_X11
USE_X11=MODULE
endif
endif
endif

ifneq ($(SYS_USE_DYNL_MODULE),YES)
ifeq ($(USE_TCLTK),MODULE)
override USE_TCLTK=YES
ifeq ($(USE_X11),MODULE)
override USE_X11=YES
endif
endif
endif

ifndef BUILD_TCLTK
BUILD_TCLTK=YES
endif

# Tcl/Tk Environment
ifeq ($(USE_TCLTK),YES)

ifndef TCLTK_VERSION
TCLTK_VERSION=8.5.18
endif

ifndef TCLTK_FEATURE_VERSION
TCLTK_FEATURE_VERSION=8.6
ifneq ($(TCLTK_VERSION),cvs)
TCLTK_FEATURE_VERSION=$(shell echo $(TCLTK_VERSION) | sed -e 's/\([0-9][0-9]*\.[0-9][0-9]*\).*/\1/')
endif
endif

_TCLTK_VERSION=8006000
ifneq ($(TCLTK_VERSION),cvs)
_TCLTK_VERSION=$(shell echo "$(TCLTK_VERSION)" | sed -e 's/ *//g' -e 's/[^0-9.].*//' -e 's/\.*$$/.0.0.0./' -e 's/^\.*//' -e 's/^\([^.][^.]*\)\.\([^.][^.]*\)\.\([^.][^.]*\)\..*/\1.000\2.000\3/' -e 's/\.[^.]*\([^.][^.][^.]\)/\1/g' -e 's/^0*\(.\)/\1/g')
endif

TCL_PATCH_SET=	$(FILES)/tcl/set-$(TCLTK_VERSION)
TK_PATCH_SET=	$(FILES)/tk/set-$(TCLTK_VERSION)

TCL_PATCH_FILES=$(shell grep -v ^\# $(TCL_PATCH_SET) | sed -e s%^%$(FILES)/tcl/%)
TK_PATCH_FILES=	$(shell grep -v ^\#  $(TK_PATCH_SET) | sed -e s%^%$(FILES)/tk/%)

ifndef TCLTK_TYPE
TCLTK_TYPE=unix
endif

ifeq ($(TCLTK_GUI_BACKEND),)
TCLTK_GUI_BACKEND=X11
endif

TCLTK_GUI_BACKEND_SETUP=FALSE

ifeq ($(TCLTK_GUI_BACKEND),X11)
TCLTK_GUI_BACKEND_SETUP=TRUE

ifndef TCLTK_FONT_SYSTEM
TCLTK_FONT_SYSTEM=Xft
endif

ifeq ($(TCLTK_FEATURE_VERSION),8.4)
override TCLTK_FONT_SYSTEM=Xlib
endif

ifneq ($(TCLTK_FONT_SYSTEM),Xft)
ifneq ($(TCLTK_FONT_SYSTEM),Aqua)
override TCLTK_FONT_SYSTEM=Xlib
endif
endif

_TCLTK_CONFIG_ARG_XFT=--disable-xft
ifeq ($(TCLTK_FONT_SYSTEM),Xft)
_TCLTK_CONFIG_ARG_XFT=--enable-xft
endif
TCLTK_CONFIG_ARG+=$(_TCLTK_CONFIG_ARG_XFT)
TCLTK_CONFIG_ARG+=--x-includes=$(X11_INCDIR) --x-libraries=$(X11_LIBDIR)

ifneq ($(USE_X11),YES)
$(error USE_TCLTK flag is couppled with USE_X11 flag!)
endif
endif

ifdef TCLTK_CONFIG_ARG_AQUA
ifeq ($(TCLTK_GUI_BACKEND),AQUA)
TCLTK_GUI_BACKEND_SETUP=TRUE
TCLTK_CONFIG_ARG+=$(TCLTK_CONFIG_ARG_AQUA)
override TCLTK_FONT_SYSTEM=Aqua
endif
endif

ifneq ($(TCLTK_GUI_BACKEND_SETUP),TRUE)
$(error GUI backend[$(TCLTK_GUI_BACKEND)] is not supported)
endif

ifndef USE_TCLTK_SHARED_LIBRARY
USE_TCLTK_SHARED_LIBRARY=YES
endif

_TCLTK_CONFIG_ARG_SHARED=--enable-shared
ifneq ($(USE_TCLTK_SHARED_LIBRARY),YES)
_TCLTK_CONFIG_ARG_SHARED=--disable-shared
endif
TCLTK_CONFIG_ARG+=$(_TCLTK_CONFIG_ARG_SHARED)

TCL_SHVER=1
TK_SHVER=1

ifndef USE_TCLTK_TRIM_VERSION
USE_TCLTK_TRIM_VERSION=NO
endif

TCLTK_LIBRARY_VERSION=$(TCLTK_FEATURE_VERSION)
ifeq ($(USE_TCLTK_TRIM_VERSION),YES)
TCLTK_LIBRARY_VERSION=$(shell echo $(TCLTK_FEATURE_VERSION) | tr -d .)
endif

TCL_SHLIB=libtcl$(TCLTK_LIBRARY_VERSION)$(SHLIB_SUFFIX)
TK_SHLIB=libtk$(TCLTK_LIBRARY_VERSION)$(SHLIB_SUFFIX)

TCLTK_MIRROR=nchc
TCLTK_MASTERSITE=http://$(TCLTK_MIRROR).dl.sourceforge.net/sourceforge/tcl
TCLTK_REPOSITORY=:pserver:anonymous@tktoolkit.cvs.sourceforge.net:/cvsroot

TCL_DISTNAME=tcl$(TCLTK_VERSION)-src.tar.gz
TK_DISTNAME=tk$(TCLTK_VERSION)-src.tar.gz
TCL_ARCHIVE=$(SOURCE_ARCHIVE_DIR)/$(TCL_DISTNAME)
TK_ARCHIVE=$(SOURCE_ARCHIVE_DIR)/$(TK_DISTNAME)

TCL_INSTALL_TARGETS=install-binaries install-libraries
TK_INSTALL_TARGETS=install-binaries install-libraries
ifeq ($(shell test $(_TCLTK_VERSION) -lt 8006000 -a $(_TCLTK_VERSION) -ge 8005013 ; echo $$?), 0)
TK_INSTALL_TARGETS=install-binaries install-headers install-libraries
endif

_TCLTK_VER=
ifneq ($(TCLTK_VERSION),cvs)
_TCLTK_VER=$(TCLTK_VERSION)
endif

ifndef TCLTK_SRCBASE
TCLTK_SRCBASE=$(CONTRIB)
endif

TCL_SRCDIR=$(TCLTK_SRCBASE)/tcl$(_TCLTK_VER)
TCL_OBJDIR=$(OBJDIR)/tcl$(_TCLTK_VER)

TK_SRCDIR=$(TCLTK_SRCBASE)/tk$(_TCLTK_VER)
TK_OBJDIR=$(OBJDIR)/tk$(_TCLTK_VER)

ifndef TCLTK_PREFIX
TCLTK_PREFIX=/usr/local
ifeq ($(BUILD_TCLTK),YES)
TCLTK_PREFIX=$(SAD_ARCH_DIR)
endif
endif # TCLTK_PREFIX

ifeq ($(BUILD_TCLTK),YES)
override TCLTK_INCDIR=$(TCLTK_PREFIX)/include
override TCLTK_LIBDIR=$(TCLTK_PREFIX)/lib
endif

ifndef TCLTK_INCDIR
TCLTK_INCDIR=$(TCLTK_PREFIX)/include
endif

ifndef TCLTK_LIBDIR
TCLTK_LIBDIR=$(TCLTK_PREFIX)/lib
endif

TCLTK_IOPT=-I$(TCLTK_INCDIR)
TCLTK_ROPT=-L$(TCLTK_LIBDIR)
TCLTK_LIBS=-L$(TCLTK_LIBDIR) -ltk$(TCLTK_LIBRARY_VERSION) -ltcl$(TCLTK_LIBRARY_VERSION)
ifeq ($(TCLTK_GUI_BACKEND)/$(TCLTK_FONT_SYSTEM),X11/Xft)
TCLTK_LIBS+=$(shell xft-config --libs)
endif

TCL_INSTALLED_MARK=$(TCLTK_LIBDIR)/.install_tcl$(TCLTK_VERSION)_done
TK_INSTALLED_MARK=$(TCLTK_LIBDIR)/.install_tk$(TCLTK_VERSION)_done

ifeq ($(BUILD_TCLTK),YES)
TCLTK_PATCH_TARGET=$(TCL_SRCDIR)/.patch_done $(TK_SRCDIR)/.patch_done
TCLTK_BUILD_TARGET=$(TCL_OBJDIR)/.build_done $(TK_OBJDIR)/.build_done
TCLTK_INSTALL_TARGET=$(TCL_INSTALLED_MARK) $(TK_INSTALLED_MARK)
MOSTLYCLEAN_DIRS+=$(TCL_SRCDIR) $(TK_SRCDIR)
endif # BUILD_TCLTK

endif # USE_TCLTK==YES

# Setup TkInter extension
ifeq ($(USE_TCLTK)/$(SYS_USE_DYNL_MODULE),MODULE/YES)
STD_MODULES+=Tkinter
endif

ifneq ($(USE_TCLTK),YES)
SYS_UNLINK_SAD_FUNCTIONS+=TkInter
endif

endif # BEFOREPORTMK

ifeq ($(AFTERPORTMK),YES)

tcltk-patch: $(TCLTK_PATCH_TARGET)

tcltk-build: $(TCLTK_BUILD_TARGET)

tcltk: $(TCLTK_INSTALL_TARGET)

ifeq ($(USE_TCLTK)/$(BUILD_TCLTK),YES/YES)

TCLTK_CONFIG_ENV+=CC="$(CC) $(SYS_CC_ABIOPT)"
TCLTK_CONFIG_ENV+=CFLAGS="$(SYS_PTHREAD_COPT) $(TCLTK_CFLAGS)"
TCLTK_CONFIG_ENV+=LDFLAGS="$(SYS_PTHREAD_LIBS)"

# Fetch
$(SOURCE_ARCHIVE_DIR)/$(TCL_DISTNAME):
	mkdir -p $(SOURCE_ARCHIVE_DIR)
	(cd $(SOURCE_ARCHIVE_DIR) ; \
	$(FETCH) $(TCLTK_MASTERSITE)/$(TCL_DISTNAME))

$(SOURCE_ARCHIVE_DIR)/$(TK_DISTNAME):
	mkdir -p $(SOURCE_ARCHIVE_DIR)
	(cd $(SOURCE_ARCHIVE_DIR) ; \
	$(FETCH) $(TCLTK_MASTERSITE)/$(TK_DISTNAME))

# Tcl/Tk Extract
ifeq ($(TCLTK_VERSION),cvs)
$(TCL_SRCDIR)/.extract_done: $(TCL_PATCH_SET) $(TCL_PATCH_FILES)
	(mkdir -p $(TCLTK_SRCBASE); rm -fr $(TCL_SRCDIR); \
	cd $(TCLTK_SRCBASE); \
	(cvs -z3 -d $(TCLTK_REPOSITORY)/tcl co -P tcl) && touch $@)
else
$(TCL_SRCDIR)/.extract_done: $(TCL_PATCH_SET) $(TCL_PATCH_FILES) $(TCL_ARCHIVE)
	(mkdir -p $(TCLTK_SRCBASE); rm -fr $(TCL_SRCDIR); \
	cd $(TCLTK_SRCBASE); \
	(gzip -dc $(TCL_ARCHIVE) | tar xvf -)    && touch $@)
endif

ifeq ($(TCLTK_VERSION),cvs)
$(TK_SRCDIR)/.extract_done: $(TK_PATCH_SET) $(TK_PATCH_FILES)
	(mkdir -p $(TCLTK_SRCBASE); rm -fr $(TK_SRCDIR); \
	cd $(TCLTK_SRCBASE); \
	(cvs -z3 -d $(TCLTK_REPOSITORY)/tktoolkit co -P tk) && touch $@)
else
$(TK_SRCDIR)/.extract_done: $(TK_PATCH_SET) $(TK_PATCH_FILES) $(TK_ARCHIVE)
	(mkdir -p $(TCLTK_SRCBASE); rm -fr $(TK_SRCDIR); \
	cd $(TCLTK_SRCBASE); \
	(gzip -dc $(TK_ARCHIVE) | tar xvf -)          && touch $@)
endif

# Tcl Base Patch
$(TCL_SRCDIR)/.patch_done: $(TCL_SRCDIR)/.extract_done $(TCL_PATCH_SET) $(TCL_PATCH_FILES)
	(cd $(TCL_SRCDIR); \
	cat /dev/null $(TCL_PATCH_FILES) | $(PATCH) -p0 && touch $@)

# Tk Base Patch
$(TK_SRCDIR)/.patch_done: $(TK_SRCDIR)/.extract_done $(TK_PATCH_SET) $(TK_PATCH_FILES)
	(cd $(TK_SRCDIR); \
	cat /dev/null $(TK_PATCH_FILES) | $(PATCH) -p0 && touch $@)

# Tcl lndir
ifeq ($(TCL_SRCDIR),$(shell $(BINDIR)/realpath $(OBJDIR)/../../contrib/tcl$(TCLTK_VERSION)))
$(TCL_OBJDIR)/.lndir_done: $(TCL_SRCDIR)/.patch_done
	(rm -fr $(TCL_OBJDIR); mkdir -p $(TCL_OBJDIR); cd $(TCL_OBJDIR); \
	 $(BINDIR)/lndir  "../../../contrib/tcl$(TCLTK_VERSION)"; \
	touch $@)
else
$(TCL_OBJDIR)/.lndir_done: $(TCL_SRCDIR)/.patch_done
	(rm -fr $(TCL_OBJDIR); mkdir -p $(TCL_OBJDIR); cd $(TCL_OBJDIR); \
	$(BINDIR)/lndir  $(TCL_SRCDIR); \
	touch $@)
endif

# Tk lndir
ifeq ($(TK_SRCDIR),$(shell $(BINDIR)/realpath $(OBJDIR)/../../contrib/tk$(TCLTK_VERSION)))
$(TK_OBJDIR)/.lndir_done: $(TK_SRCDIR)/.patch_done
	(rm -fr $(TK_OBJDIR); mkdir -p $(TK_OBJDIR); cd $(TK_OBJDIR); \
	$(BINDIR)/lndir  "../../../contrib/tk$(TCLTK_VERSION)"; \
	touch $@)
else
$(TK_OBJDIR)/.lndir_done: $(TK_SRCDIR)/.patch_done
	(rm -fr $(TK_OBJDIR); mkdir -p $(TK_OBJDIR); cd $(TK_OBJDIR); \
	$(BINDIR)/lndir  $(TK_SRCDIR); \
	touch $@)
endif

# Tcl Build
$(TCL_OBJDIR)/.configure_done: $(TCL_OBJDIR)/.lndir_done
	(cd $(TCL_OBJDIR)/$(TCLTK_TYPE); \
	env $(TCLTK_CONFIG_ENV) \
	    $(TCLTK_TCL_CONFIG_ENV) ./configure --prefix=$(TCLTK_PREFIX) \
		$(TCLTK_CONFIG_ARG) $(TCLTK_TCLCONFIG_ARG) && \
	touch $@)

$(TCL_OBJDIR)/.build_done:	$(TCL_OBJDIR)/.configure_done
	(cd $(TCL_OBJDIR)/$(TCLTK_TYPE); $(MAKE) && touch $@)

$(TCL_INSTALLED_MARK):	$(TCL_OBJDIR)/.build_done
	@rm -f $(TCLTK_LIBDIR)/.install_tcl*_done
	(cd $(TCL_OBJDIR)/$(TCLTK_TYPE); mkdir -p $(TCLTK_PREFIX) ; \
	$(MAKE) $(TCL_INSTALL_TARGETS) ; \
	if [ ! -e $(TCLTK_LIBDIR)/$(TCL_SHLIB) \
	    -a -e $(TCLTK_LIBDIR)/$(TCL_SHLIB).$(TCL_SHVER) ]; then \
	    ln -s $(TCL_SHLIB).$(TCL_SHVER) $(TCLTK_LIBDIR)/$(TCL_SHLIB) ; \
	fi) && touch $@

# Tk Build
$(TK_OBJDIR)/.configure_done:	$(TK_OBJDIR)/.lndir_done \
				$(TCL_INSTALLED_MARK)
	(cd $(TK_OBJDIR)/$(TCLTK_TYPE); \
	env $(TCLTK_CONFIG_ENV) \
	    $(TCLTK_TK_CONFIG_ENV)  ./configure --prefix=$(TCLTK_PREFIX) \
	    	$(TCLTK_CONFIG_ARG) $(TCLTK_TK_CONFIG_ARG) && \
	touch $@)

$(TK_OBJDIR)/.build_done:	$(TK_OBJDIR)/.configure_done
	(cd $(TK_OBJDIR)/$(TCLTK_TYPE); $(MAKE) && touch $@)

$(TK_INSTALLED_MARK):	$(TK_OBJDIR)/.build_done
	@rm -f $(TCLTK_LIBDIR)/.install_tk*_done
	(cd $(TK_OBJDIR)/$(TCLTK_TYPE); mkdir -p $(TCLTK_PREFIX) ; \
	$(MAKE) $(TK_INSTALL_TARGETS) ; \
	if [ ! -e $(TCLTK_LIBDIR)/$(TK_SHLIB) \
	    -a -e $(TCLTK_LIBDIR)/$(TK_SHLIB).$(TK_SHVER) ]; then \
	    ln -s $(TK_SHLIB).$(TK_SHVER) $(TCLTK_LIBDIR)/$(TK_SHLIB) ; \
	fi) && touch $@

endif # USE_TCLTK&&BUILD_TCLTK

endif # AFTERPORTMK

# End of File
