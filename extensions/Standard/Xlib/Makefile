# SAD/Xlib module

MODULE_TYPE=Extension
MODULE_NAME=Xlib
MODULE_SUBDIR=SAD
MODULE_SCRIPT=Xlib.n
OBJS=tfXlib_.o

GENSRC_FILES=		tfXlib_.c

COPT_ADD=	-DWITH_EXTENSION_MODULE

#USE_TCLTK=YES
USE_X11=YES

#SADSRCDIR=$(error Set top directory of SAD source tree)

module:

# Module Local Rules

tfXlib_.c:	$(SADSRCDIR)/src/tfXlib_.c
	ln -s -f $< $@

# Don't touch after here

ifdef SADSRCDIR
include $(SADSRCDIR)/mk/sad.module.mk
COPT_ADD+=	$(X11_COPT) $(X11_IOPT)
ROPT_ADD+=	$(X11_ROPT)
LDOPT_ADD+=	$(SYS_LDOPT_RUNTIME_PATH) $(X11_LIBS)
else
all depend module clean distclean mostlyclean install install-module install-script:
	@echo "Edit SADSRCDIR in Makefile or Type $(MAKE) with \`SADSRCDIR=(Top directory of SAD source tree)'"
	@echo "Type $(MAKE) module to build extension module"
	@echo "Type $(MAKE) install to install extension module"
endif

# End of File
