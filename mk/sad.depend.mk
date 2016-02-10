# New Makefile -- Fortran Code Dependence Part
#

SAD_DEPEND_C=$(SRCDIR)/.depend_c
SAD_DEPEND_F=$(SRCDIR)/.depend_f

.PHONY:	depend_f depend_c
depend: $(PRE_DEPEND) depend_f depend_c depend-std-modules

depend_f:
ifeq ($(OS_TYPE),CYGWIN)
	(cd $(SRCDIR) ; \
	env PATH="$${PATH}:$(BINDIR)" bash \
	    $(BINDIR)/mkdep *.f sim/*.f) > ${SAD_DEPEND_F}
else
	(cd $(SRCDIR) ; \
	env PATH="$${PATH}:$(BINDIR)" BIN_SH=xpg4 \
	    $(BINDIR)/mkdep *.f sim/*.f) > ${SAD_DEPEND_F}
endif

depend_c:
ifeq ($(OS_TYPE),CYGWIN)
	(cd $(SRCDIR) ; \
	env PATH="$${PATH}:$(BINDIR)" bash \
	    $(BINDIR)/mkdep -Igdtoa *.c sim/*.c gdtoa/*.c) > ${SAD_DEPEND_C}
else
	(cd $(SRCDIR) ; \
	env PATH="$${PATH}:$(BINDIR)" BIN_SH=xpg4 \
	    $(BINDIR)/mkdep -Igdtoa *.c sim/*.c gdtoa/*.c) > ${SAD_DEPEND_C}
endif

depend-std-modules:
	@for mod in /dev/null $(STD_MODULES); do \
	    if [ -x $(SADDIR)/extensions/Standard/$${mod} ]; then \
		(cd $(SADDIR)/extensions/Standard/$${mod}; \
		 unset OS_NAME MAKEFLAGS; \
		 env SADSRCDIR=$(SADDIR) $(MAKE) depend); \
	    fi \
	done

-include $(SRCDIR)/.depend_f
-include $(SRCDIR)/.depend_c

# End of File
