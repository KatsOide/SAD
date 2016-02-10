# New Makefile -- libffi Part
#
ifeq ($(BEFOREPORTMK),YES)

ifeq ($(DEBUG),LIBFFI)
debug:
	@echo "USE_LIBFFI=$(USE_LIBFFI)"
	@echo "BUILD_LIBFFI(USER/SYS): $(BUILD_LIBFFI)/$(_SYS_BUILD_LIBFFI)"
	@echo "LIBFFI_INCDIR=$(LIBFFI_INCDIR)"
	@echo "LIBFFI_LIBDIR=$(LIBFFI_LIBDIR)"
	@echo "LIBFFI_LIBS=$(LIBFFI_LIBS)"
	@echo "_USE_LIBFFI=$(_USE_LIBFFI)"
	@echo "_LIBFFI_HEADER_CANDIDATES=$(_LIBFFI_HEADER_CANDIDATES)"
	@echo "_LIBFFI_LIBRARY_CANDIDATES=$(_LIBFFI_LIBRARY_CANDIDATES)"
	@echo "_LIBFFI_HEADER=$(_LIBFFI_HEADER)"
	@echo "_LIBFFI_LIBRARY_=$(_LIBFFI_LIBRARY)"
	@echo "_LIBFFI_HEADER_PREFIX=$(_LIBFFI_HEADER_PREFIX)"
	@echo "_LIBFFI_PREFIX=$(_LIBFFI_PREFIX)"
	@echo "_LIBFFI_LIBDIR=$(_LIBFFI_LIBDIR)"
endif

## libffi Environment

ifeq ($(USE_LIBFFI),YES)

# Current default version to build
ifndef LIBFFI_VERSION
LIBFFI_VERSION=3.0.9
endif

LIBFFI_SHVER=5

# Search libffi.h
_LIBFFI_HEADER_CANDIDATES=

_LIBFFI_HEADER_CANDIDATES+=$(shell ls -r $(SAD_ARCH_DIR)/lib*/libffi-*/include/ffi.h 2>/dev/null)

ifdef LIBFFI_INCDIR
_LIBFFI_HEADER_CANDIDATES+=$(shell ls $(LIBFFI_INCDIR)/ffi.h 2>/dev/null)
endif

ifdef LIBFFI_PREFIX
_LIBFFI_HEADER_CANDIDATES+=$(shell ls $(LIBFFI_PREFIX)/include/ffi.h 2>/dev/null)
_LIBFFI_HEADER_CANDIDATES+=$(shell ls $(LIBFFI_PREFIX)/lib*/libffi-*/include/ffi.h 2>/dev/null)
endif

_LIBFFI_HEADER_CANDIDATES+=$(shell ls /usr/local/lib*/libffi-*/include/ffi.h 2>/dev/null)

_LIBFFI_HEADER_CANDIDATES+=$(shell ls /usr/lib*/libffi-*/include/ffi.h 2>/dev/null)

_LIBFFI_HEADER_CANDIDATES+=$(shell ls -f /usr/local/include/ffi.h /usr/pkg/include/ffi.h /usr/include/ffi.h 2>/dev/null)

_LIBFFI_HEADER=$(shell for f in $(_LIBFFI_HEADER_CANDIDATES); do if grep 'libffi  *[0-9].*Copyright' $$f >/dev/null 2>/dev/null; then echo $$f; break; fi; done)

# _USE_* <- major * 10^6 + minor * 10^3 + patch-level[minor,patch-level < 10^3]
_USE_LIBFFI=0
ifneq ($(_LIBFFI_HEADER),)
_USE_LIBFFI=$(shell grep 'libffi  *[0-9].*Copyright' $(_LIBFFI_HEADER) | sed -e 's/^.*libffi  *\([^ ][^ ]*\).*/\1/' | sed -e 's/ *//g' -e 's/[^0-9.].*//' -e 's/\.*$$/.0.0.0./' -e 's/^\.*//' -e 's/^\([^.][^.]*\)\.\([^.][^.]*\)\.\([^.][^.]*\)\..*/\1.000\2.000\3/' -e 's/\.[^.]*\([^.][^.][^.]\)/\1/g' -e 's/^0*\(.\)/\1/g')
_LIBFFI_HEADER_PREFIX=$(subst /include/ffi.h,,$(_LIBFFI_HEADER))
_LIBFFI_PREFIX=$(shell echo $(_LIBFFI_HEADER_PREFIX) | sed -e 's@/lib[^/]*/libffi-[0-9][^/]*$$@@')
endif

# Search libffi library
_LIBFFI_LIBRARY_CANDIDATES=

ifneq ($(_LIBFFI_HEADER),)
ifneq ($(_LIBFFI_HEADER_PREFIX),$(_LIBFFI_PREFIX))
_LIBFFI_LIBRARY_CANDIDATES+=$(shell ls $(shell echo $(_LIBFFI_HEADER_PREFIX) | sed -e 's@/libffi-[0-9][^/]*$$@@')/libffi.* 2>/dev/null)
endif
_LIBFFI_LIBRARY_CANDIDATES+=$(shell ls $(_LIBFFI_PREFIX)/lib*/libffi.* 2>/dev/null)
endif

_LIBFFI_LIBRARY_CANDIDATES+=$(shell ls -r $(SAD_ARCH_DIR)/lib*/libffi.* 2>/dev/null)

ifdef LIBFFI_LIBDIR
_LIBFFI_LIBRARY_CANDIDATES+=$(shell ls $(LIBFFI_LIBDIR)/libffi.* 2>/dev/null)
endif

ifdef LIBFFI_PREFIX
_LIBFFI_LIBRARY_CANDIDATES+=$(shell ls -r $(LIBFFI_PREFIX)/lib*/libffi.* 2>/dev/null)
endif

_LIBFFI_LIBRARY=$(shell for f in $(_LIBFFI_LIBRARY_CANDIDATES); do echo $$f; break; done)

ifneq ($(_LIBFFI_LIBRARY),)
_LIBFFI_LIBDIR=$(shell echo $(_LIBFFI_LIBRARY) | sed -e 's@/libffi\.[^/]*@@')
endif

# Build libffi by framework?
_SYS_BUILD_LIBFFI=YES
ifeq ($(shell test $(_USE_LIBFFI) -ge 3000000 ; echo $$?),0)
ifneq ($(_LIBFFI_PREFIX),)
ifneq ($(realpath $(_LIBFFI_PREFIX)),$(realpath $(SAD_ARCH_DIR)))
_SYS_BUILD_LIBFFI=NO
endif # _LIBFFI_PREFIX!=SAD_ARCH_DIR
endif # !empty(_LIBFFI_PREFIX)
endif

ifndef BUILD_LIBFFI
BUILD_LIBFFI=$(_SYS_BUILD_LIBFFI)
endif

# Setup LIBFFI* prefix for build
__LIBFFI_PREFIX=$(_LIBFFI_PREFIX)
__LIBFFI_INCDIR=$(_LIBFFI_HEADER_PREFIX)/include
__LIBFFI_LIBDIR=$(_LIBFFI_LIBDIR)
ifeq ($(BUILD_LIBFFI),YES)
__LIBFFI_PREFIX=$(SAD_ARCH_DIR)
ifdef LIBFFI_PREFIX
__LIBFFI_PREFIX=$(LIBFFI_PREFIX)
endif
__LIBFFI_INCDIR=$(__LIBFFI_PREFIX)/lib/libffi-$(LIBFFI_VERSION)/include
ifdef LIBFFI_INCDIR
__LIBFFI_INCDIR=$(LIBFFI_INCDIR)
endif
__LIBFFI_LIBDIR=$(__LIBFFI_PREFIX)/lib
ifdef LIBFFI_LIBDIR
__LIBFFI_LIBDIR=$(LIBFFI_LIBDIR)
endif
endif	# BUILD_LIBFFI==YES

override LIBFFI_INCDIR:=$(__LIBFFI_INCDIR)
override LIBFFI_LIBDIR:=$(__LIBFFI_LIBDIR)

LIBFFI_IOPT=-I$(LIBFFI_INCDIR)
LIBFFI_LIBS=-L$(LIBFFI_LIBDIR) -lffi
ifeq ($(realpath $(__LIBFFI_PREFIX)),$(realpath $(SAD_ARCH_DIR)))
LIBFFI_ROPT=-L$(LIBFFI_LIBDIR)
endif

LIBFFI_MASTERSITE=ftp://sourceware.org/pub/libffi

LIBFFI_DISTNAME=libffi-$(LIBFFI_VERSION).tar.gz
LIBFFI_ARCHIVE=$(SOURCE_ARCHIVE_DIR)/$(LIBFFI_DISTNAME)

LIBFFI_PATCH_SET=$(FILES)/libffi/set-$(LIBFFI_VERSION)
LIBFFI_PATCH_FILES=$(shell grep -v ^\# $(LIBFFI_PATCH_SET) | sed -e s%^%$(FILES)/libffi/%)

ifndef LIBFFI_SRCBASE
LIBFFI_SRCBASE=$(CONTRIB)
endif

LIBFFI_SRCDIR=$(LIBFFI_SRCBASE)/libffi-$(LIBFFI_VERSION)
LIBFFI_OBJDIR=$(OBJDIR)/libffi-$(LIBFFI_VERSION)

LIBFFI_CONFIG_ENV+=CC="$(CC) $(SYS_CC_ABIOPT)"

LIBFFI_INSTALLED_MARK=$(LIBFFI_LIBDIR)/.install_libffi-$(LIBFFI_VERSION)_done

ifeq ($(BUILD_LIBFFI),YES)
LIBFFI_PATCH_TARGET=$(LIBFFI_SRCDIR)/.patch_done
LIBFFI_BUILD_TARGET=$(LIBFFI_OBJDIR)/.build_done
LIBFFI_INSTALL_TARGET=$(LIBFFI_INSTALLED_MARK)
MOSTLYCLEAN_DIRS+=$(LIBFFI_SRCDIR)
endif # BUILD_LIBFFI

endif # USE_LIBFFI

endif # BEFOREPORTMK

ifeq ($(AFTERPORTMK),YES)

libffi-patch: $(LIBFFI_PATCH_TARGET)

libffi-build: $(LIBFFI_BUILD_TARGET)

libffi: $(LIBFFI_INSTALL_TARGET)

ifeq ($(USE_LIBFFI)/$(BUILD_LIBFFI),YES/YES)

# Fetch
$(SOURCE_ARCHIVE_DIR)/$(LIBFFI_DISTNAME):
	mkdir -p $(SOURCE_ARCHIVE_DIR)
	(cd $(SOURCE_ARCHIVE_DIR) ; \
	$(FETCH) $(LIBFFI_MASTERSITE)/$(LIBFFI_DISTNAME))

# Extract
$(LIBFFI_SRCDIR)/.extract_done: $(LIBFFI_PATCH_SET) $(LIBFFI_PATCH_FILES) $(LIBFFI_ARCHIVE)
	(mkdir -p $(LIBFFI_SRCBASE); rm -fr $(LIBFFI_SRCDIR); \
	cd $(LIBFFI_SRCBASE); \
	(gzip -dc $(LIBFFI_ARCHIVE) | tar xvf -) && touch $@)

# Patch
$(LIBFFI_SRCDIR)/.patch_done: $(LIBFFI_SRCDIR)/.extract_done $(LIBFFI_PATCH_SET) $(LIBFFI_PATCH_FILES)
	(cd $(LIBFFI_SRCDIR); \
	cat /dev/null $(LIBFFI_PATCH_FILES) | $(PATCH) -p0 && touch $@)

# lndir
ifeq ($(shell test -r $(SADDIR)/$(CONF_FILE) ; echo $$?),0)
_LIBFFI_LINDIR_DEPEND=$(SADDIR)/$(CONF_FILE)
endif

ifeq ($(LIBFFI_SRCDIR),$(shell $(BINDIR)/realpath $(OBJDIR)/../../contrib/libffi-$(LIBFFI_VERSION)))
$(LIBFFI_OBJDIR)/.lndir_done: $(LIBFFI_SRCDIR)/.patch_done $(_LIBFFI_LINDIR_DEPEND)
	(rm -fr $(LIBFFI_OBJDIR); mkdir -p $(LIBFFI_OBJDIR); \
	cd $(LIBFFI_OBJDIR); $(BINDIR)/lndir "../../../contrib/libffi-$(LIBFFI_VERSION)"; \
	touch $@)
else
$(LIBFFI_OBJDIR)/.lndir_done: $(LIBFFI_SRCDIR)/.patch_done $(_LIBFFI_LINDIR_DEPEND)
	(rm -fr $(LIBFFI_OBJDIR); mkdir -p $(LIBFFI_OBJDIR); \
	cd $(LIBFFI_OBJDIR); $(BINDIR)/lndir $(LIBFFI_SRCDIR); \
	touch $@)
endif

# Configure
$(LIBFFI_OBJDIR)/.configure_done: $(LIBFFI_OBJDIR)/.lndir_done
	(cd $(LIBFFI_OBJDIR); \
	env $(LIBFFI_CONFIG_ENV) ./configure --prefix=$(_LIBFFI_PREFIX) \
	    $(LIBFFI_CONFIG_ARG) && touch $@)

# Build
$(LIBFFI_OBJDIR)/.build_done: $(LIBFFI_OBJDIR)/.configure_done
	(cd $(LIBFFI_OBJDIR); $(MAKE) && touch $@)

# Install
$(LIBFFI_INSTALLED_MARK): $(LIBFFI_OBJDIR)/.build_done
	(cd $(LIBFFI_OBJDIR); \
	mkdir -p $(LIBFFI_INCDIR); \
	for f in ffi.h ffitarget.h; do \
	    $(INSTALL_DATA) include/$$f	$(LIBFFI_INCDIR)/; \
	done; \
	mkdir -p $(LIBFFI_LIBDIR); \
	for f in libffi.a; do \
	    $(INSTALL_DATA)   .libs/$$f	$(LIBFFI_LIBDIR)/; \
	    $(RANLIB)			$(LIBFFI_LIBDIR)/$$f; \
	done; \
	f=libffi$(SHLIB_SUFFIX); if [ -r .libs/$$f.$(LIBFFI_SHVER) ]; then \
	    $(INSTALL_DATA)   .libs/$$f.$(LIBFFI_SHVER) $(LIBFFI_LIBDIR)/; \
	    ln -s -f		    $$f.$(LIBFFI_SHVER) $(LIBFFI_LIBDIR)/$$f; \
	fi; \
	touch $@)

endif # USE_LIBFFI&&BUILD_LIBFFI

endif # AFTERPORTMK

# End of File
