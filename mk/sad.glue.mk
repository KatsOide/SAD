# New Makefile -- Library Glue Part
#

OBJGLUE_EPICS=dummyroutCA.o
ifeq ($(USE_EPICS),YES)
OBJGLUE_EPICS=tfefun2.o tfefun3ep.o DbStatic.o CaSearch2.o CaSearch_.o
#OBJGLUE_EPICS=tfefun2.o tfefun3ep.o CaSearch2.o CaSearch_.o
OPT_INCDIR+=$(EPICS_IOPT)
OPT_RLIBDIR+=$(EPICS_ROPT)
OPT_LIBS+=$(EPICS_LIBS)
endif

OBJGLUE= dummyrout.o
OBJGLUE+=$(OBJGLUE_EPICS)

# End of File
