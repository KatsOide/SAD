# SADScript Dynamic Loading Extension Module for Mersenne Twister

MODULE_TYPE=Extension
MODULE_NAME=	MT19937
MODULE_SUBDIR=	Random
#MODULE_SCRIPT=
OBJS=driver.o
USE_FRAMEWORK=YES

COPT_ADD=	-DWITH_EXTENSION_MODULE

COMMON=$(TOPDIR)/../Common

GENSRC_FILES=	core.c

MASTER_SITE=http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/MT2002/CODES
MT_CORE_SRC=mt19937ar.c

MT_CORE_ARCHIVE=$(SOURCE_ARCHIVE_DIR)/$(MT_CORE_SRC)

#SADSRCDIR=$(error Set top directory of SAD source tree)

module:

# Module Local Rules


# Don't touch after here

ifdef SADSRCDIR
include $(SADSRCDIR)/mk/sad.module.mk

COPT_ADD+=	-I$(COMMON)

ifdef WITH_PARABOLA_TRANSFORM
COPT_ADD+=	-DUSE_PARABOLA_TRANSFORM
endif

ifndef WITHOUT_TRUE_53BIT_RES
COPT_ADD+=	-DUSE_TRUE_53BIT_RES
endif

# Fetch
$(SOURCE_ARCHIVE_DIR)/$(MT_CORE_SRC):
	mkdir -p $(SOURCE_ARCHIVE_DIR)
	(cd $(SOURCE_ARCHIVE_DIR);\
	$(FETCH) $(MASTER_SITE)/$(MT_CORE_SRC))

core.c:	$(MT_CORE_ARCHIVE)
	ln -s $(MT_CORE_ARCHIVE) $@

else
all depend module clean distclean mostlyclean install install-module install-script:
	@echo "Edit SADSRCDIR in Makefile or Type $(MAKE) with \`SADSRCDIR=(Top directory of SAD source tree)'"
	@echo "Type $(MAKE) module to build extension module"
	@echo "Type $(MAKE) install to install extension module"
endif

# End of File
