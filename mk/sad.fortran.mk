# New Makefile -- Fortran Support Part
#

# for missing Fortran standard math functions
_FORTRAN_MISSING_DMATH_FUNCS=$(shell echo '$(FORTRAN_MISSING_DMATH_FUNCS)' | tr '\t' '\n' | tr ' ' '\n' | sort -u | tr '\n' ' ' | sed -e 's@  *@ @g' -e 's@ $$@@')
OBJ_F_DMATH=$(shell echo '$(_FORTRAN_MISSING_DMATH_FUNCS)' | sed -e 's@[ 	]*\([_a-zA-Z][_a-zA-Z0-9]*\)[ 	]*@sim/fortran_math_\1_.o @g' -e 's/ $$//')

# Complex tangent function
OBJ_F_ZTAN=
ifneq ($(HAVE_F_ZTAN),YES)
OBJ_F_ZTAN=ztan.o
endif

# FSEEK
OBJ_F_FSEEK=
ifneq ($(HAVE_F_FSEEK),YES)
OBJ_F_FSEEK=sim/fseek_Dummy_.o
ifeq ($(HAVE_F_FSEEK),GNU)
OBJ_F_FSEEK=sim/fseek_subroutine.o sim/fseek_subroutine_.o
endif
endif

# File descriptor getter from LUN
OBJ_F_GETFD=
ifneq ($(HAVE_F_GETFD),YES)
OBJ_F_GETFD=sim/getfd_Dummy.o
ifeq ($(HAVE_F_PXFFILENO),YES)
OBJ_F_GETFD=sim/getfd_PXFFILENO.o
endif
ifeq ($(HAVE_F_FNUM),YES)
OBJ_F_GETFD=sim/getfd_FNUM.o
endif
endif

# Binary write routine
OBJ_F_WRITEB=writeb_f95.o
ifeq ($(HAVE_F_FORMAT_DOLLAR_EDIT),YES)
OBJ_F_WRITEB=writeb_dollar.o
endif
# Disabled on MAIN trunk(BFA not yet)
OBJ_F_WRITEB=

# itfgetbuf implementation/FGETC SIM
#
# src/itfgetbuf.f	by READ with Q edit descriptor
# src/itfgetbuf_PPC.f	by FGETC function
# src/itfgetbuf_.c	emulated by read(2)
#			Breaking Fortran I/O buffering
# src/itfgetbuf_f90.f	emulated by READ with ADVANCE option
#			Can't handle '\r' character correctly!
# Note:	FGETC is required for binary Read[] in src/tfwrite.f
#
OBJ_F_ITGETBUF=itfgetbuf_.o

ifeq ($(HAVE_F_FGETC),YES)
OBJ_F_ITGETBUF=itfgetbuf_PPC.o
ifeq ($(HAVE_F_FORMAT_Q_EDIT),YES)
OBJ_F_ITGETBUF=itfgetbuf.o
endif
endif

# Fortran I/O SIM for openbuf
OBJ_F_FortranIO=itopenbuf.o itfopenread.o
ifneq ($(HAVE_F_OPEN_DISP_OPT),YES)
OBJ_F_FortranIO=itopenbuf.o itfopenread_G77.o
endif
ifeq ($(USE_NEW_FORTRAN_OPEN_SIM),YES)
OBJ_F_FortranIO=sim/fortran_io_.o sim/fortran_io.o
endif

# Collet Fortran support objects
OBJ_FORTRAN=	\
	sim/unix_fortran_.o \
	$(OBJ_F_FortranIO) \
	$(OBJ_F_DMATH) $(OBJ_F_ZTAN) \
	$(OBJ_F_FSEEK) $(OBJ_F_GETFD) $(OBJ_F_WRITEB) $(OBJ_F_ITGETBUF)

# ------ Don't touch ------
PRE_DEPEND+=	genFortranMathInc
PRE_BUILD+=	genFortranMathInc
# Generate include for missing Fortran math functions
genFortranMathInc:
	@( \
	target=inc/MATH.inc; \
	temp=`mktemp -q $${target}.XXXXXXXXXXXXXXXX`; \
	if [ $$? -ne 0 ]; then \
	    echo "Can't create temporary file!"; \
	    exit 1; \
	else \
	    echo "c     for missing standard math functions"	  >> $${temp}; \
	    echo ""						  >> $${temp}; \
	    for func in $(_FORTRAN_MISSING_DMATH_FUNCS); do \
		echo "c     for src/sim/fortran_math_$${func}_.c" >> $${temp}; \
		echo "      real(kind=8), external :: $${func}"	  >> $${temp}; \
		echo ""						  >> $${temp}; \
	    done; \
	    echo "c     End of File"				  >> $${temp}; \
	    if cmp -s $${target} $${temp}; then \
		rm -f $${temp}; \
	    else \
		mv -f $${temp} $${target}; \
		rm -f `grep $${target} $(SRCDIR)/*.f 2>/dev/null | sed -e 's@^$(SRCDIR)/@@' -e 's@.f:.*@.o@'`; \
	    fi \
	fi)

ifeq ($(DEBUG),FORTRAN)
debug:
	@echo " FORTRAN_MISSING_DMATH_FUNCS=$(FORTRAN_MISSING_DMATH_FUNCS)"
	@echo "_FORTRAN_MISSING_DMATH_FUNCS=$(_FORTRAN_MISSING_DMATH_FUNCS)"
	@echo "HAVE_F_ZTAN=$(HAVE_F_ZTAN)"
	@echo "HAVE_F_FORMAT_Q_EDIT=$(HAVE_F_FORMAT_Q_EDIT)"
	@echo "HAVE_F_FORMAT_DOLLAR_EDIT=$(HAVE_F_FORMAT_DOLLAR_EDIT)"
	@echo "HAVE_F_FSEEK=$(HAVE_F_FSEEK)"
	@echo "HAVE_F_FGETC=$(HAVE_F_FGETC)"
	@echo "HAVE_F_GETFD=$(HAVE_F_GETFD)"
	@echo "HAVE_F_FNUM=$(HAVE_F_FNUM)"
	@echo "HAVE_F_PXFFILENO=$(HAVE_F_PXFFILENO)"
	@echo "OBJ_F_DMATH=$(OBJ_F_DMATH)"
	@echo "OBJ_FORTRAN=$(OBJ_FORTRAN)"
endif

# End of File
