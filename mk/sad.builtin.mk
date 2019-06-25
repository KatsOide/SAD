# New Makefile -- Built-in Function Part
#

ifeq ($(BEFOREPORTMK),YES)
SYS_LINK_SAD_FUNCTIONS=
SYS_UNLINK_SAD_FUNCTIONS=

# ---- Begin SAD functinal module list ----

SYS_LINK_SAD_FUNCTIONS+=DynamicLink
OBJ_FUNC_DynamicLink=tfDynamicLink_.o


SYS_LINK_SAD_FUNCTIONS+=FeatureQ
OBJ_FUNC_FeatureQ=tfFeatureQ_.o


SYS_LINK_SAD_FUNCTIONS+=BuildInfo
OBJ_FUNC_BuildInfo=tfBuildInfo_.o


SYS_LINK_SAD_FUNCTIONS+=Process
OBJ_FUNC_Process=tfProcess_.o


SYS_LINK_SAD_FUNCTIONS+=Crypt
OBJ_FUNC_Crypt=tfCrypt_.o
LIBS_FUNC_Crypt=$(SYS_LIBCRYPT)


SYS_LINK_SAD_FUNCTIONS+=FileIO
OBJ_FUNC_FileIO=tfFileIO_.o
ifeq ($(HAVE_MKSTEMPS),YES)
SYS_COPT_tfFileIO_+=-DHAVE_MKSTEMPS
endif


SYS_LINK_SAD_FUNCTIONS+=NetworkIO
OBJ_FUNC_NetworkIO=tfNetworkIO_.o


SYS_LINK_SAD_FUNCTIONS+=NetSemaphore
OBJ_FUNC_NetSemaphore=tfNetSemaphore_.o netSemaphore.o


SYS_LINK_SAD_FUNCTIONS+=TkInter
OBJ_FUNC_TkInter=tfTkInter_.o tfTclArg.o tfcanvasclip.o
COPT_FUNC_TkInter=$(TCLTK_COPT)
IOPT_FUNC_TkInter=$(TCLTK_IOPT)
ROPT_FUNC_TkInter=$(TCLTK_ROPT)
LIBS_FUNC_TkInter=$(TCLTK_LIBS)


SYS_LINK_SAD_FUNCTIONS+=Xlib
OBJ_FUNC_Xlib=tfXlib_.o
COPT_FUNC_Xlib=$(X11_COPT)
IOPT_FUNC_Xlib=$(X11_IOPT)
ROPT_FUNC_Xlib=$(X11_ROPT)
LIBS_FUNC_Xlib=$(X11_LIBS)


SYS_LINK_SAD_FUNCTIONS+=DateTAI
OBJ_FUNC_DateTAI=tfDateTAI_.o
USE_FUNC_DateTAI=LIBTAI
COPT_FUNC_DateTAI=$(LIBTAI_COPT)
IOPT_FUNC_DateTAI=$(LIBTAI_IOPT)
ROPT_FUNC_DateTAI=$(LIBTAI_ROPT)
LIBS_FUNC_DateTAI=$(LIBTAI_LIBS)


SYS_LINK_SAD_FUNCTIONS+=Random
OBJ_FUNC_Random=tfRandom_.o


SYS_LINK_SAD_FUNCTIONS+=RenumberElement
OBJ_FUNC_RenumberElement=tfRenumberElement_.o tfrenumber_element.o


# ---- End   SAD functinal module list ----
endif # BEFOREPORTMK

# ------ Don't touch ------
# Collect initialization hook, objects and related library
$(foreach _FUNC,$(SYS_LINK_SAD_FUNCTIONS),$(eval _WITH_FUNC_$(_FUNC)=YES))
$(foreach _FUNC,$(SYS_UNLINK_SAD_FUNCTIONS),$(eval WITHOUT_FUNC_$(_FUNC)=YES))
$(foreach _FUNC,$(UNLINK_SAD_FUNCTIONS),$(eval WITHOUT_FUNC_$(_FUNC)=YES))

_LINK_SAD_FUNCTIONS=$(SYS_LINK_SAD_FUNCTIONS) $(shell echo '$(foreach _FUNC,$(LINK_SAD_FUNCTIONS),$(if $(_WITH_FUNC_$(_FUNC)),,$(_FUNC)))' | tr '\t' '\n' | tr ' ' '\n' | sort -u | tr '\n' ' ' | sed -e 's@  *@ @g' -e 's@ $$@@')

_LINK_SAD_FUNCS=$(foreach _FUNC,$(_LINK_SAD_FUNCTIONS),$(if $(WITHOUT_FUNC_$(_FUNC)),,$(_FUNC)))

ifeq ($(AFTERPORTMK),YES)
_SAD_FUNC_DEFS=$(foreach _FUNC,$(_LINK_SAD_FUNCS),$(_FUNC))
_SAD_FUNC_OBJS=$(foreach _FUNC,$(_LINK_SAD_FUNCS),$(OBJ_FUNC_$(_FUNC)))
_SAD_FUNC_COPT=$(foreach _FUNC,$(_LINK_SAD_FUNCS),$(COPT_FUNC_$(_FUNC)))
_SAD_FUNC_IOPT=$(foreach _FUNC,$(_LINK_SAD_FUNCS),$(IOPT_FUNC_$(_FUNC)))
_SAD_FUNC_ROPT=$(foreach _FUNC,$(_LINK_SAD_FUNCS),$(ROPT_FUNC_$(_FUNC)))
_SAD_FUNC_LIBS=$(foreach _FUNC,$(_LINK_SAD_FUNCS),$(LIBS_FUNC_$(_FUNC)))
endif # AFTERPORTMK
_SAD_FUNC_USES=$(foreach _FUNC,$(_LINK_SAD_FUNCS),$(USE_FUNC_$(_FUNC)))
$(foreach _MOD,$(_SAD_FUNC_USES),$(eval USE_$(_MOD)=YES))

ifeq ($(BEFOREPORTMK),YES)
ifeq ($(DEBUG),BUILTIN)
debug:
	@echo "SYS_LINK_SAD_FUNCTIONS=$(SYS_LINK_SAD_FUNCTIONS)"
	@echo "    LINK_SAD_FUNCTIONS=$(LINK_SAD_FUNCTIONS)"
	@echo "SYS_UNLINK_SAD_FUNCTIONS=$(SYS_UNLINK_SAD_FUNCTIONS)"
	@echo "    UNLINK_SAD_FUNCTIONS=$(UNLINK_SAD_FUNCTIONS)"
	@echo "_LINK_SAD_FUNCTIONS=$(_LINK_SAD_FUNCTIONS)"
	@echo "_SAD_FUNC_DEFS=$(_SAD_FUNC_DEFS)"
	@echo "_SAD_FUNC_OBJS=$(_SAD_FUNC_OBJS)"
	@echo "_SAD_FUNC_COPT=$(_SAD_FUNC_COPT)"
	@echo "_SAD_FUNC_IOPT=$(_SAD_FUNC_IOPT)"
	@echo "_SAD_FUNC_ROPT=$(_SAD_FUNC_ROPT)"
	@echo "_SAD_FUNC_LIBS=$(_SAD_FUNC_LIBS)"
	@echo "_SAD_FUNC_USES=$(_SAD_FUNC_USES)"
endif
endif # BEFOREPORTMK

# End of File
