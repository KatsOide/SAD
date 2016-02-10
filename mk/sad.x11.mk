# New Makefile -- X11 Part
#

# Check USE_X11
ifndef USE_X11
USE_X11=NO
ifeq ($(SYS_USE_DYNL_MODULE),YES)
USE_X11=MODULE
endif
endif

ifneq ($(SYS_USE_DYNL_MODULE),YES)
ifeq ($(USE_X11),MODULE)
override USE_X11=YES
endif
endif

# X11 Environment
ifeq ($(USE_X11),YES)

ifndef X11_PREFIX
X11_PREFIX=/usr/X11R6
endif

ifndef X11_INCDIR
X11_INCDIR=$(X11_PREFIX)/include
endif

ifndef X11_LIBDIR
X11_LIBDIR=$(X11_PREFIX)/lib
endif

X11_IOPT=-I$(X11_INCDIR)
#X11_ROPT=-L$(X11_LIBDIR)
X11_LIBS=-L$(X11_LIBDIR) -lX11

ifeq ($(COMPILER),GCC)
X11_COPT+=-Wno-variadic-macros
endif

endif # USE_X11==YES

# Setup Xlib extension
ifeq ($(USE_X11)/$(SYS_USE_DYNL_MODULE),MODULE/YES)
STD_MODULES+=Xlib
endif

ifneq ($(USE_X11),YES)
ifneq ($(USE_TCLTK)/$(TCLTK_GUI_BACKEND),YES/AQUA)
SYS_UNLINK_SAD_FUNCTIONS+=Xlib
endif
endif

# End of File
