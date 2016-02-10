# External Module Build Driver
SADSRCDIR?=../..

all:	module

depend module install install-module install-script clean distclean mostlyclean:
	@env SADSRCDIR=`$(SADSRCDIR)/bin/realpath $(SADSRCDIR)` sh -c ' \
	for dir in `ls`; do \
	    if [     -d $${dir} \
		-a   -f $${dir}/Makefile \
		-a ! -f $${dir}/.skip_build ]; then \
		(cd $${dir}; $(MAKE) $@); \
	    fi; \
	done'

# End of File
