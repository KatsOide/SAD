# New Makefile -- libtai Part
#
ifeq ($(BEFOREPORTMK),YES)

ifeq ($(DEBUG),LIBTAI)
debug:
	@echo "USE_LIBTAI=$(USE_LIBTAI)"
	@echo "BUILD_LIBTAI(USER/SYS): $(BUILD_LIBTAI)/$(_SYS_BUILD_LIBTAI)"
	@echo "LIBTAI_INCDIR=$(LIBTAI_INCDIR)"
	@echo "LIBTAI_LIBDIR=$(LIBTAI_LIBDIR)"
	@echo "LIBTAI_LIBS=$(LIBTAI_LIBS)"
	@echo "_USE_LIBTAI=$(_USE_LIBTAI)"
	@echo "_LIBTAI_PREFIX=$(_LIBTAI_PREFIX)"
	@echo "_LIBTAI_HEADER_PREFIX=$(_LIBTAI_HEADER_PREFIX)"
	@echo "_LIBTAI_HEADER=$(_LIBTAI_HEADER)"
	@echo "_LIBTAI_HEADER_CANDIDATES=$(_LIBTAI_HEADER_CANDIDATES)"
endif

## libtai Environment

ifeq ($(USE_LIBTAI),YES)

ifndef LIBTAI_ADDITIONAL_LEAPSECONDS
LIBTAI_ADDITIONAL_LEAPSECONDS=	+2005-12-31 +2008-12-31 +2012-06-30
endif

LIBTAI_LAST_LEAPSECONDS=$(shell echo $(LIBTAI_ADDITIONAL_LEAPSECONDS) | tr ' ' '\n' | tail -n 1)

ifndef LIBTAI_VERSION
LIBTAI_VERSION=0.60
endif

# Search tai.h
_LIBTAI_HEADER_CANDIDATES=

ifdef LIBTAI_INCDIR
_LIBTAI_HEADER_CANDIDATES+=$(shell ls -r $(LIBTAI_INCDIR)/libtai/tai.h 2>/dev/null)
endif

ifdef LIBTAI_PREFIX
_LIBTAI_HEADER_CANDIDATES+=$(shell ls -r $(LIBTAI_PREFIX)/include/libtai/tai.h 2>/dev/null)
endif

_LIBTAI_HEADER_CANDIDATES+=/usr/local/include/libtai/tai.h

_LIBTAI_HEADER=$(shell for h in $(_LIBTAI_HEADER_CANDIDATES); do if grep 'tai_add' $${h} >/dev/null 2>/dev/null; then echo $${h}; break; fi; done)

ifneq ($(_LIBTAI_HEADER),)
_LIBTAI_HEADER_PREFIX=$(subst /include/libtai/tai.h,,$(_LIBTAI_HEADER))
_LIBTAI_PREFIX=$(_LIBTAI_HEADER_PREFIX)
endif

# Build libtai by framework?
_SYS_BUILD_LIBTAI=YES
ifneq ($(_LIBTAI_PREFIX),)
ifneq ($(realpath $(_LIBTAI_PREFIX)),$(realpath $(SAD_ARCH_DIR)))
_SYS_BUILD_LIBTAI=NO
endif # _LIBTAI_PREFIX!=SAD_ARCH_DIR
endif # !empty(_LIBTAI_PREFIX)

# Setup _LIBTAI*_PREFIX for build
_LIBTAI_BUILD_SETUP=$(_SYS_BUILD_LIBTAI)
ifeq ($(BUILD_LIBTAI),YES)
_LIBTAI_BUILD_SETUP=YES
endif

ifeq ($(_LIBTAI_BUILD_SETUP),YES)
_LIBTAI_PREFIX=$(SAD_ARCH_DIR)
ifdef LIBTAI_PREFIX
_LIBTAI_PREFIX=$(LIBTAI_PREFIX)
endif
_LIBTAI_HEADER_PREFIX=$(_LIBTAI_PREFIX)
endif	# _LIBTAI_BUILD_SETUP==YES

ifndef BUILD_LIBTAI
BUILD_LIBTAI=$(_SYS_BUILD_LIBTAI)
endif

ifeq ($(BUILD_LIBTAI),YES)
override LIBTAI_INCDIR=$(_LIBTAI_HEADER_PREFIX)/include
override LIBTAI_LIBDIR=$(_LIBTAI_PREFIX)/lib
endif

ifndef LIBTAI_INCDIR
LIBTAI_INCDIR=$(_LIBTAI_HEADER_PREFIX)/include
endif

ifndef LIBTAI_LIBDIR
LIBTAI_LIBDIR=$(_LIBTAI_PREFIX)/lib
endif

LIBTAI_IOPT=-I$(LIBTAI_INCDIR)
LIBTAI_LIBS=-L$(LIBTAI_LIBDIR) -ltai
ifeq ($(realpath $(_LIBTAI_PREFIX)),$(realpath $(SAD_ARCH_DIR)))
LIBTAI_ROPT=-L$(LIBTAI_LIBDIR)
endif

ifndef LIBTAI_LEAPSECONDS_PREFIX
LIBTAI_LEAPSECONDS_PREFIX=$(_LIBTAI_PREFIX)/etc/libtai
ifeq ($(realpath $(_LIBTAI_PREFIX)),$(realpath $(SAD_ARCH_DIR)))
LIBTAI_LEAPSECONDS_PREFIX=$(SAD_SHARE_ROOT)/libtai
endif
endif

ifeq ($(HAVE_LEAPSECONDS),YES)
_LIBTAI_LEAPSECONDS_PREFIX=$(LIBTAI_LEAPSECONDS_PREFIX)
else
_LIBTAI_LEAPSECONDS_PREFIX=
endif

LIBTAI_MASTERSITE=http://cr.yp.to/libtai

LIBTAI_DISTNAME=libtai-$(LIBTAI_VERSION).tar.gz
LIBTAI_ARCHIVE=$(SOURCE_ARCHIVE_DIR)/$(LIBTAI_DISTNAME)

LIBTAI_PATCH_SET=$(FILES)/libtai/set-$(LIBTAI_VERSION)
LIBTAI_PATCH_FILES=$(shell grep -v ^\# $(LIBTAI_PATCH_SET) | sed -e s%^%$(FILES)/libtai/%)

ifndef LIBTAI_SRCBASE
LIBTAI_SRCBASE=$(CONTRIB)
endif

LIBTAI_SRCDIR=$(LIBTAI_SRCBASE)/libtai-$(LIBTAI_VERSION)
LIBTAI_OBJDIR=$(OBJDIR)/libtai-$(LIBTAI_VERSION)

LIBTAI_INSTALLED_MARK=$(LIBTAI_LEAPSECONDS_PREFIX)/.install_libtai-$(LIBTAI_VERSION)$(LIBTAI_LAST_LEAPSECONDS)_done

ifeq ($(BUILD_LIBTAI),YES)
LIBTAI_PATCH_TARGET=$(LIBTAI_SRCDIR)/.post_patch_done
LIBTAI_BUILD_TARGET=$(LIBTAI_OBJDIR)/.post_build_done
LIBTAI_INSTALL_TARGET=$(LIBTAI_INSTALLED_MARK)
MOSTLYCLEAN_DIRS+=$(LIBTAI_SRCDIR)
endif # BUILD_LIBTAI

endif # USE_LIBTAI

endif # BEFOREPORTMK

ifeq ($(AFTERPORTMK),YES)

libtai-patch: $(LIBTAI_PATCH_TARGET)

libtai-build: $(LIBTAI_BUILD_TARGET)

libtai: $(LIBTAI_INSTALL_TARGET)

ifeq ($(USE_LIBTAI)/$(BUILD_LIBTAI),YES/YES)

# Fetch
$(SOURCE_ARCHIVE_DIR)/$(LIBTAI_DISTNAME):
	mkdir -p $(SOURCE_ARCHIVE_DIR)
	(cd $(SOURCE_ARCHIVE_DIR) ; \
	$(FETCH) $(LIBTAI_MASTERSITE)/$(LIBTAI_DISTNAME))

# Extract
$(LIBTAI_SRCDIR)/.extract_done: $(LIBTAI_PATCH_SET) $(LIBTAI_PATCH_FILES) $(LIBTAI_ARCHIVE)
	(mkdir -p $(LIBTAI_SRCBASE); rm -fr $(LIBTAI_SRCDIR); \
	cd $(LIBTAI_SRCBASE); \
	(gzip -dc $(LIBTAI_ARCHIVE) | tar xvf -) && touch $@)

# Patch
$(LIBTAI_SRCDIR)/.patch_done: $(LIBTAI_SRCDIR)/.extract_done $(LIBTAI_PATCH_SET) $(LIBTAI_PATCH_FILES)
	(cd $(LIBTAI_SRCDIR); \
	cat /dev/null $(LIBTAI_PATCH_FILES) | $(PATCH) -p0 && touch $@)

$(LIBTAI_SRCDIR)/.post_patch_done: $(LIBTAI_SRCDIR)/.patch_done
	(cd $(LIBTAI_SRCDIR); \
	for f in "# leapseconds after libtai release" \
			$(LIBTAI_ADDITIONAL_LEAPSECONDS) \
				"# End of File"; do \
	    echo "$$f" >>leapsecs.txt; \
	done; \
	for f in *.h; do \
	    cp -f $$f $$f.bak; \
	    sed -e 's|\(#include[ 	]*\)"\([^"][^"]*\)"|\1<libtai/\2>|' $$f.bak >$$f; \
	done; \
	touch $@)

# lndir
ifeq ($(shell test -r $(SADDIR)/$(CONF_FILE) ; echo $$?),0)
_LIBTAI_LINDIR_DEPEND=$(SADDIR)/$(CONF_FILE)
endif

ifeq ($(LIBTAI_SRCDIR),$(shell $(BINDIR)/realpath $(OBJDIR)/../../contrib/libtai-$(LIBTAI_VERSION)))
$(LIBTAI_OBJDIR)/.lndir_done: $(LIBTAI_SRCDIR)/.post_patch_done $(_LIBTAI_LINDIR_DEPEND)
	(rm -fr $(LIBTAI_OBJDIR); mkdir -p $(LIBTAI_OBJDIR); \
	cd $(LIBTAI_OBJDIR); $(BINDIR)/lndir "../../../contrib/libtai-$(LIBTAI_VERSION)"; \
	touch $@)
else
$(LIBTAI_OBJDIR)/.lndir_done: $(LIBTAI_SRCDIR)/.post_patch_done $(_LIBTAI_LINDIR_DEPEND)
	(rm -fr $(LIBTAI_OBJDIR); mkdir -p $(LIBTAI_OBJDIR); \
	cd $(LIBTAI_OBJDIR); $(BINDIR)/lndir $(LIBTAI_SRCDIR); \
	touch $@)
endif

# Configure
$(LIBTAI_OBJDIR)/.configure_done: $(LIBTAI_OBJDIR)/.lndir_done
	(cd $(LIBTAI_OBJDIR); \
	rm -f libtai conf-cc conf-ld; \
	ln -s . libtai; \
	echo "$(CC) $(SYS_CC_ABIOPT) -I. $(LIB_COPT)" >conf-cc; \
	echo "$(CC) $(SYS_CC_ABIOPT)" >conf-ld; \
	for f in leapsecs_read.c leapsecs.3; do \
	    cp -f $$f $$f.bak; \
	    sed -e 's|%%LEAPSECONDS_PREFIX%%|$(_LIBTAI_LEAPSECONDS_PREFIX)|' $$f.bak >$$f; \
	done; \
	touch $@)

# Build
$(LIBTAI_OBJDIR)/.build_done: $(LIBTAI_OBJDIR)/.configure_done
	(cd $(LIBTAI_OBJDIR); $(MAKE) it && touch $@)

$(LIBTAI_OBJDIR)/.post_build_done: $(LIBTAI_OBJDIR)/.build_done
	(cd $(LIBTAI_OBJDIR); rm -f leapsecs.dat; \
	./leapsecs <leapsecs.txt >leapsecs.dat && touch $@)

# Install
$(LIBTAI_INSTALLED_MARK): $(LIBTAI_OBJDIR)/.post_build_done
	(cd $(LIBTAI_OBJDIR); \
	mkdir -p $(LIBTAI_LEAPSECONDS_PREFIX); \
	for f in leapsecs.txt leapsecs.dat; do \
	    $(INSTALL_DATA) $$f $(LIBTAI_LEAPSECONDS_PREFIX)/$$f; \
	done; \
	mkdir -p $(_LIBTAI_PREFIX)/bin; \
	for f in leapsecs yearcal easter nowutc; do \
	    $(INSTALL_PROG) $$f $(_LIBTAI_PREFIX)/bin/$$f; \
	done; \
	mkdir -p $(LIBTAI_INCDIR)/libtai; \
	for f in *.h; do \
	    $(INSTALL_DATA) $$f $(LIBTAI_INCDIR)/libtai/$$f; \
	done; \
	mkdir -p $(LIBTAI_LIBDIR); \
	for f in libtai.a; do \
	    $(INSTALL_DATA) $$f $(LIBTAI_LIBDIR)/$$f; \
	    $(RANLIB) $(LIBTAI_LIBDIR)/$$f; \
	done; \
	touch $@)

endif # USE_LIBTAI&&BUILD_LIBTAI

endif # AFTERPORTMK

# End of File
