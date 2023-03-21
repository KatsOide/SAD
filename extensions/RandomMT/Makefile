# SADScript Dynamic Loading Extension Module for Mersenne Twister

SUB_MODULES=
SUB_MODULES+=Plugins/SAD
SUB_MODULES+=Plugins/MT19937
SUB_MODULES+=Plugins/SFMT

MODULE_TYPE=Script
MODULE_NAME=RandomMT
#MODULE_SUBDIR=
MODULE_SCRIPT=	RandomMT.n RandomSFMT.n

#SADSRCDIR=$(error Set top directory of SAD source tree)

module:

# Module Local Rules


# Don't touch after here

ifdef SADSRCDIR
include $(SADSRCDIR)/mk/sad.module.mk
else
all depend module clean distclean mostlyclean install install-module install-script:
	@echo "Edit SADSRCDIR in Makefile or Type $(MAKE) with \`SADSRCDIR=(Top directory of SAD source tree)'"
	@echo "Type $(MAKE) module to build extension module"
	@echo "Type $(MAKE) install to install extension module"
endif

# End of File
