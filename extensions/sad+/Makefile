# SAD+ Differential Algebra External Engine Module

MODULE_TYPE=Engine
MODULE_NAME=sadp
#MODULE_SUBDIR=
MODULE_SCRIPT=	TaylorMap.n
SRCDIR=src
OBJS=	sadp.o \
	ran.o splined.o errf.o teigen.o \
	da_init.o dalib.o p_dalib.o c_dalib.o c_p_dalib.o \
	map.o matrix.o lin_map.o \
	beam.o track.o lie_da.o hamilton.o \
	element.o bb.o impedance.o linelem.o ph_rot.o \
	drift.o bend.o quad.o sext.o thin.o wig.o solenoid.o cav.o

LD_ENGINE=$(CXX) $(SYS_CXX_ABIOPT)
CXXOPT_ADD=	-I$(TOPDIR)/include
OPT_LIBS+=	-lm

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
