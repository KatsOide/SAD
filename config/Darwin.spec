# New Makefile -- OS Dependence Part for Darwin
#
-include $(SADCONF)/4.4BSD-Lite.spec

SYS_LIBCRYPT=

RANLIB=ranlib -s
SHLIB_SUFFIX=.dylib

SYS_LDOPT_RUNTIME_PATH_OPT=
SYS_LDOPT_EXPORT_SYMBOL=

LD_SHARED=$(CC) $(SYS_CC_ABIOPT) -dynamiclib \
    -flat_namespace -undefined suppress -Wl,-single_module

TCLTK_CONFIG_ARG_AQUA=--enable-corefoundation --enable-aqua
ifneq ($(TCLTK_GUI_BACKEND),AQUA)
TCLTK_CONFIG_ARG+=--disable-corefoundation --disable-aqua
endif
ifeq ($(TCLTK_GUI_BACKEND),AQUA)
SYS_COPT_tfTkInter_+=-DAQUA_TK
endif

# End of File
