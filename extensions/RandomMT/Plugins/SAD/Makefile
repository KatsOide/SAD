# SADScript Dynamic Loading Extension Module for Mersenne Twister

MODULE_TYPE=Extension
MODULE_NAME=	SAD
MODULE_SUBDIR=	Random
#MODULE_SCRIPT=
OBJS=driver.o
USE_FRAMEWORK=YES

COPT_ADD=	-DWITH_EXTENSION_MODULE

#SADSRCDIR=$(error Set top directory of SAD source tree)

module:

# Module Local Rules


# Don't touch after here

ifdef SADSRCDIR
include $(SADSRCDIR)/mk/sad.module.mk

ifdef WITH_PARABOLA_TRANSFORM
COPT_ADD+=	-DUSE_PARABOLA_TRANSFORM
endif

else
all depend module clean distclean mostlyclean install install-module install-script:
	@echo "Edit SADSRCDIR in Makefile or Type $(MAKE) with \`SADSRCDIR=(Top directory of SAD source tree)'"
	@echo "Type $(MAKE) module to build extension module"
	@echo "Type $(MAKE) install to install extension module"
endif

# End of File
