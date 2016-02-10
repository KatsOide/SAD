# SAD/Tkiner module

MODULE_TYPE=Extension
MODULE_NAME=Tkinter
MODULE_SUBDIR=SAD
MODULE_SCRIPT=Tkinter.n
OBJS=tfTkInter_.o tfTclArg.o tfcanvasclip.o

GENSRC_FILES=		tfTkInter_.c tfTclArg.f tfcanvasclip.f

COPT_ADD=	-DWITH_EXTENSION_MODULE

USE_TCLTK=YES
USE_X11=YES

MODULE_PRE_BUILD_TARGET=$(TCLTK_INSTALL_TARGET)
MODULE_EXTRA_TARGET=tcltk-patch tcltk-build tcltk

#SADSRCDIR=$(error Set top directory of SAD source tree)

module:

# Module Local Rules

tfTkInter_.c:	$(SADSRCDIR)/src/tfTkInter_.c
	ln -s -f $< $@

tfTclArg.f:	$(SADSRCDIR)/src/tfTclArg.f
	ln -s -f $< $@

tfcanvasclip.f:	$(SADSRCDIR)/src/tfcanvasclip.f
	ln -s -f $< $@

# Don't touch after here

ifdef SADSRCDIR
include $(SADSRCDIR)/mk/sad.module.mk
COPT_ADD+=	$(TCLTK_COPT) $(X11_COPT) $(TCLTK_IOPT) $(X11_IOPT)
ROPT_ADD+=	$(TCLTK_ROPT)
LDOPT_ADD+=	$(SYS_LDOPT_RUNTIME_PATH) $(TCLTK_LIBS)
else
all depend module clean distclean mostlyclean install install-module install-script:
	@echo "Edit SADSRCDIR in Makefile or Type $(MAKE) with \`SADSRCDIR=(Top directory of SAD source tree)'"
	@echo "Type $(MAKE) module to build extension module"
	@echo "Type $(MAKE) install to install extension module"
endif

# End of File
