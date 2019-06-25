# Generic Object Makefile OS Compiler Support Part -- NAG
#

ifeq ($(DEBUG),NAG)
debug:
	@echo "SYS_FOPT_ENABLE_BACKSLASH_ESCAPE=$(SYS_FOPT_ENABLE_BACKSLASH_ESCAPE)"
	@echo "SYS_FOPT_DISABLE_BACKSLASH_ESCAPE=$(SYS_FOPT_DISABLE_BACKSLASH_ESCAPE)"
	@echo "OPT_RLIBDIR=$(OPT_RLIBDIR)"
endif

SYS_FC=f95

SYS_FOPT_DYNL_PIC=   -PIC

SYS_FOPT=	-kind=byte -maxcontin=40
#SYS_FOPT+=	-ieee=full
#SYS_FOPT+=	-fcfuns

# Object specific options

# Default debug build options
FOPT=		-g -O1 -gline -nan

# End of File
