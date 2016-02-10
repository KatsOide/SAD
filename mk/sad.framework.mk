# Framework component handling

ifeq ($(BEFOREPORTMK),YES)

ifeq ($(DEBUG),FRAMEWORK)
debug:
	@echo "FRAMEWORK_HEADERS=$(FRAMEWORK_HEADERS)"
	@echo "FRAMEWORK_OBJS=$(FRAMEWORK_OBJS)"
	@echo "FRAMEWORK_INIT_HOOKS=$(FRAMEWORK_INIT_HOOKS)"
endif

# Framework header/object/hook list
FRAMEWORK_HEADERS=
FRAMEWORK_OBJS=
FRAMEWORK_INIT_HOOKS=

# Internal SADScript function call table framework
FRAMEWORK_INIT_HOOKS+=	FuncTBL

# Internal Dynamic Loader MI helper framework
FRAMEWORK_INIT_HOOKS+=	DynlMI

# Feature framework
FRAMEWORK_HEADERS+=	feature.h
FRAMEWORK_OBJS+=	feature.o
FRAMEWORK_INIT_HOOKS+=	Feature

# BuildInfo framework
FRAMEWORK_HEADERS+=	buildinfo.h
FRAMEWORK_OBJS+=	buildinfo.o buildinfo-db.o buildinfo_.o
FRAMEWORK_INIT_HOOKS+=	BuildInfo

# Random driver framework
FRAMEWORK_HEADERS+=	random_driver.h
FRAMEWORK_OBJS+=	random_driver.o
FRAMEWORK_INIT_HOOKS+=	Random

# Random plugin(SAD)
FRAMEWORK_OBJS+=	random_plugin_sad.o
FRAMEWORK_INIT_HOOKS+=	RandomPlugin_SAD

endif # BEFOREPORTMK

ifeq ($(AFTERPORTMK),YES)

ifeq ($(USE_FRAMEWORK)/$(_SYS_IN_MODULE_MK),YES/YES)
OPT_INCDIR+=-I$(SRCDIR)/framework
endif # USE_FRAMEWORK && _SYS_IN_MODULE_MK

endif # AFTERPORTMK

# End of File
