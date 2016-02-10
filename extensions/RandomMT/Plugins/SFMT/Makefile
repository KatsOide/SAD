# SADScript Dynamic Loading Extension Module for Mersenne Twister

MODULE_TYPE=Extension
MODULE_NAME=	SFMT$(SFMT_MEXP)
MODULE_SUBDIR=	Random
MODULE_SCRIPT=	SFMT.n
OBJS=driver_$(SFMT_MEXP).o
USE_FRAMEWORK=YES

COPT_ADD=	-DWITH_EXTENSION_MODULE
COPT_ADD+=	-DUSE_INLINE_GENRAND

COMMON=$(TOPDIR)/../Common

GENSRC_FILES=		$(SFMT_CORE) \
			driver_$(SFMT_MEXP).c

MOSTLYCLEAN_FILES=	driver_*.c

MASTER_SITE=http://www.math.sci.hiroshima-u.ac.jp/%7Em-mat/MT/SFMT
SFMT_CORE_VERSION=1.3.3

SFMT_CORE=SFMT-src-$(SFMT_CORE_VERSION)
SFMT_CORE_TGZ=$(SFMT_CORE).tar.gz
SFMT_ARCHIVE=$(SOURCE_ARCHIVE_DIR)/$(SFMT_CORE_TGZ)

ifndef SFMT_MEXP
SFMT_MEXP=19937
endif

#SADSRCDIR=$(error Set top directory of SAD source tree)

module:

# Module Local Rules


# Don't touch after here

ifdef SADSRCDIR
include $(SADSRCDIR)/mk/sad.module.mk

COPT_ADD+=	-I$(COMMON) -I$(TOPDIR)/$(SFMT_CORE)

COPT_ADD+=	-DMEXP=$(SFMT_MEXP) -DSTR_MEXP=\"$(SFMT_MEXP)\"

ifdef WITH_PARABOLA_TRANSFORM
COPT_ADD+=	-DUSE_PARABOLA_TRANSFORM
endif

ifndef WITHOUT_TRUE_53BIT_RES
COPT_ADD+=	-DUSE_TRUE_53BIT_RES
endif

ifdef USE_SSE2
COPT_ADD+=	-DHAVE_SSE2 $(USE_SSE2)
endif

ifdef USE_ALTIVEC
COPT_ADD+=	-DHAVE_ALTIVEC $(USE_ALTIVEC)
endif

# Fetch
$(SOURCE_ARCHIVE_DIR)/$(SFMT_CORE_TGZ):
	mkdir -p $(SOURCE_ARCHIVE_DIR)
	(cd $(SOURCE_ARCHIVE_DIR); \
	$(FETCH) $(MASTER_SITE)/$(SFMT_CORE_TGZ))

$(SFMT_CORE):	$(SFMT_ARCHIVE)
	gzip -dc $(SFMT_ARCHIVE) | tar xvf -

driver_$(SFMT_MEXP).c:	driver.c
	ln -s driver.c $@

else
all depend module clean distclean mostlyclean install install-module install-script:
	@echo "Edit SADSRCDIR in Makefile or Type $(MAKE) with \`SADSRCDIR=(Top directory of SAD source tree)'"
	@echo "Type $(MAKE) module to build extension module"
	@echo "Type $(MAKE) install to install extension module"
endif

# End of File
