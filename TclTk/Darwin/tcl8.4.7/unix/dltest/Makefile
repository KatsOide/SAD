# Generated automatically from Makefile.in by configure.
# This Makefile is used to create several test cases for Tcl's load
# command.  It also illustrates how to take advantage of configuration
# exported by Tcl to set up Makefiles for shared libraries.
# RCS: @(#) $Id: Makefile.in,v 1.11.2.1 2004/07/20 11:13:10 das Exp $

TCL_DBGX =		
CC = gcc -pipe
LIBS =			-L/Users/oide/SAD/oldsad/TclTk/Darwin/tcl8.4.7/unix -ltclstub8.4${TCL_DBGX}   -framework CoreFoundation 
AC_FLAGS =		 -DMAC_OSX_TCL=1 -DHAVE_CFBUNDLE=1 -DUSE_VFORK=1 -DTCL_DEFAULT_ENCODING=\"utf-8\" -DTCL_WIDE_INT_TYPE=long\ long -DWORDS_BIGENDIAN=1 -DHAVE_GETCWD=1 -DHAVE_OPENDIR=1 -DHAVE_STRSTR=1 -DHAVE_STRTOL=1 -DHAVE_STRTOLL=1 -DHAVE_STRTOULL=1 -DHAVE_TMPNAM=1 -DHAVE_WAITPID=1 -DNO_VALUES_H=1 -DHAVE_LIMITS_H=1 -DHAVE_UNISTD_H=1 -DHAVE_SYS_PARAM_H=1 -DUSE_TERMIOS=1 -DHAVE_SYS_TIME_H=1 -DTIME_WITH_SYS_TIME=1 -DHAVE_TM_ZONE=1 -DHAVE_GMTIME_R=1 -DHAVE_LOCALTIME_R=1 -DHAVE_TM_GMTOFF=1 -DHAVE_ST_BLKSIZE=1 -DSTDC_HEADERS=1 -DHAVE_SIGNED_CHAR=1 -DHAVE_PUTENV_THAT_COPIES=1 -DHAVE_LANGINFO=1 -DHAVE_SYS_IOCTL_H=1 -DHAVE_SYS_FILIO_H=1 
SHLIB_CFLAGS =		-fno-common
SHLIB_LD =		cc -dynamiclib ${LDFLAGS}
SHLIB_LD_LIBS =		${LIBS}
SHLIB_SUFFIX =		.dylib
SRC_DIR =		.
TCL_VERSION=		8.4

CFLAGS_DEBUG		= -g
CFLAGS_OPTIMIZE		= -Os

CFLAGS =  $(CFLAGS_OPTIMIZE)
CC_SWITCHES = $(CFLAGS) -I${SRC_DIR}/../../generic -DTCL_MEM_DEBUG \
	${SHLIB_CFLAGS} -DUSE_TCL_STUBS ${AC_FLAGS}

all: pkga${SHLIB_SUFFIX} pkgb${SHLIB_SUFFIX} pkgc${SHLIB_SUFFIX} pkgd${SHLIB_SUFFIX} pkge${SHLIB_SUFFIX}
	@touch ../dltest.marker

pkga${SHLIB_SUFFIX}: $(SRC_DIR)/pkga.c
	$(CC) -c $(CC_SWITCHES) $(SRC_DIR)/pkga.c
	${SHLIB_LD} -o pkga${SHLIB_SUFFIX} pkga.o ${SHLIB_LD_LIBS}

pkgb${SHLIB_SUFFIX}: $(SRC_DIR)/pkgb.c
	$(CC) -c $(CC_SWITCHES) $(SRC_DIR)/pkgb.c
	${SHLIB_LD} -o pkgb${SHLIB_SUFFIX} pkgb.o ${SHLIB_LD_LIBS}

pkgc${SHLIB_SUFFIX}: $(SRC_DIR)/pkgc.c
	$(CC) -c $(CC_SWITCHES) $(SRC_DIR)/pkgc.c
	${SHLIB_LD} -o pkgc${SHLIB_SUFFIX} pkgc.o ${SHLIB_LD_LIBS}

pkgd${SHLIB_SUFFIX}: $(SRC_DIR)/pkgd.c
	$(CC) -c $(CC_SWITCHES) $(SRC_DIR)/pkgd.c
	${SHLIB_LD} -o pkgd${SHLIB_SUFFIX} pkgd.o ${SHLIB_LD_LIBS}

pkge${SHLIB_SUFFIX}: $(SRC_DIR)/pkge.c
	$(CC) -c $(CC_SWITCHES) $(SRC_DIR)/pkge.c
	${SHLIB_LD} -o pkge${SHLIB_SUFFIX} pkge.o ${SHLIB_LD_LIBS}

clean:
	rm -f *.o *${SHLIB_SUFFIX} config.cache config.log config.status
	rm -f lib.exp ../dltest.marker

distclean: clean
	rm -f Makefile
