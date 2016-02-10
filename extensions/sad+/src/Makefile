COMPFLAGS     =  

BLINCLINK     =

CXX	      = g++
LD	      =	$(CXX)

CFLAGS	      = -I../include $(BLINCLINK) $(COMPFLAGS)
PFLAGS	      = -I../include $(BLINCLINK) $(COMPFLAGS)
FFLAGS	      = -I../include $(BLINCLINK) $(COMPFLAGS)
CCFLAGS       = $(BLINCLINK) $(COMPFLAGS)
CXXFLAGS      = -I../include $(BLINCLINK) $(COMPFLAGS)
LDFLAGS	      = $(COMPFLAGS)

OBJS	      = beam.o \
		bb.o \
		ph_rot.o \
		bend.o \
		sadp.o \
		c_dalib.o \
		c_p_dalib.o \
		cav.o \
		da_init.o \
		dalib.o \
		drift.o \
		element.o \
		errf.o \
		hamilton.o \
		impedance.o \
		lie_da.o \
		linelem.o \
		lin_map.o \
		map.o \
		matrix.o \
		p_dalib.o \
		quad.o \
		ran.o \
		sext.o \
		solenoid.o \
		splined.o \
		teigen.o \
		thin.o \
		track.o \
		wig.o

PROGRAM       = sadp.exe

$(PROGRAM):           $(OBJS)
		@echo "Linking $(PROGRAM) ..."
		$(LD) -o $(PROGRAM) $(OBJS) -lm
		@echo "done"

.cc.o:	$(@:.o=.cc)
	$(CXX) -c $(CXXFLAGS) $(@:.o=.cc)

clean:
	rm -f $(OBJS) $(PROGRAM)

# copy from .depend
bb.o: bb.cc ../include/Complex.h ../include/bb.h ../include/c_da.h ../include/c_p_da.h ../include/dacpp.h ../include/element.h ../include/lin_map.h ../include/map_da.h ../include/map_double.h ../include/map_p_da.h ../include/matrix.h ../include/p_da.h ../include/phys_const.h ../include/track.h
beam.o: beam.cc ../include/Complex.h ../include/c_da.h ../include/dacpp.h ../include/element.h ../include/lin_map.h ../include/map_da.h ../include/map_double.h ../include/map_p_da.h ../include/matrix.h ../include/p_da.h ../include/phys_const.h ../include/track.h
bend.o: bend.cc ../include/Complex.h ../include/bend.h ../include/c_da.h ../include/dacpp.h ../include/element.h ../include/lin_map.h ../include/map_da.h ../include/map_double.h ../include/map_p_da.h ../include/matrix.h ../include/p_da.h ../include/phys_const.h ../include/track.h
c_dalib.o: c_dalib.cc ../include/Complex.h ../include/c_da.h ../include/dacpp.h ../include/phys_const.h
c_p_dalib.o: c_p_dalib.cc ../include/Complex.h ../include/c_da.h ../include/c_p_da.h ../include/dacpp.h ../include/p_da.h ../include/phys_const.h
cav.o: cav.cc ../include/Complex.h ../include/c_da.h ../include/cavity.h ../include/dacpp.h ../include/element.h ../include/lin_map.h ../include/map_da.h ../include/map_double.h ../include/map_p_da.h ../include/matrix.h ../include/p_da.h ../include/phys_const.h ../include/track.h
concate.o: concate.cc ../include/Complex.h ../include/c_da.h ../include/dacpp.h ../include/lin_map.h ../include/map_c_da.h ../include/map_da.h ../include/map_double.h ../include/map_p_da.h ../include/matrix.h ../include/p_da.h ../include/track.h
dalib.o: dalib.cc ../include/Complex.h ../include/c_da.h ../include/dacpp.h ../include/lin_map.h ../include/map_c_da.h ../include/map_da.h ../include/map_double.h ../include/map_p_da.h ../include/matrix.h ../include/p_da.h ../include/track.h concate.cc oldmul.cc perturb.cc
drift.o: drift.cc ../include/Complex.h ../include/c_da.h ../include/dacpp.h ../include/element.h ../include/lin_map.h ../include/map_da.h ../include/map_double.h ../include/map_p_da.h ../include/matrix.h ../include/p_da.h ../include/phys_const.h ../include/track.h
element.o: element.cc ../include/Complex.h ../include/bb.h ../include/bend.h ../include/c_da.h ../include/cavity.h ../include/dacpp.h ../include/element.h ../include/impedance.h ../include/lin_map.h ../include/linelem.h ../include/map_da.h ../include/map_double.h ../include/map_p_da.h ../include/matrix.h ../include/p_da.h ../include/ph_rot.h ../include/phys_const.h ../include/quad.h ../include/sext.h ../include/solenoid.h ../include/thin.h ../include/track.h ../include/wiggler.h
errf.o: errf.cc ../include/Complex.h
hamilton.o: hamilton.cc ../include/Complex.h ../include/c_da.h ../include/dacpp.h ../include/lie_da.h ../include/lin_map.h ../include/map_da.h ../include/map_double.h ../include/matrix.h ../include/track.h
impedance.o: impedance.cc ../include/Complex.h ../include/c_da.h ../include/dacpp.h ../include/element.h ../include/impedance.h ../include/lin_map.h ../include/map_da.h ../include/map_double.h ../include/map_p_da.h ../include/matrix.h ../include/p_da.h ../include/phys_const.h ../include/track.h
lie_da.o: lie_da.cc ../include/Complex.h ../include/c_da.h ../include/dacpp.h ../include/lie_da.h ../include/lin_map.h ../include/map_c_da.h ../include/map_da.h ../include/map_double.h ../include/matrix.h ../include/track.h
lin_map.o: lin_map.cc ../include/dacpp.h ../include/lin_map.h ../include/matrix.h ../include/phys_const.h
linelem.o: linelem.cc ../include/Complex.h ../include/c_da.h ../include/dacpp.h ../include/element.h ../include/lin_map.h ../include/linelem.h ../include/map_da.h ../include/map_double.h ../include/map_p_da.h ../include/matrix.h ../include/p_da.h ../include/phys_const.h ../include/track.h
map.o: map.cc ../include/Complex.h ../include/c_da.h ../include/dacpp.h ../include/lin_map.h ../include/map_c_da.h ../include/map_da.h ../include/map_double.h ../include/map_p_da.h ../include/matrix.h ../include/p_da.h ../include/track.h
matrix.o: matrix.cc ../include/matrix.h
p_dalib.o: p_dalib.cc ../include/dacpp.h ../include/p_da.h
perturb.o: perturb.cc ../include/Complex.h
ph_rot.o: ph_rot.cc ../include/Complex.h ../include/c_da.h ../include/dacpp.h ../include/element.h ../include/lin_map.h ../include/map_da.h ../include/map_double.h ../include/map_p_da.h ../include/matrix.h ../include/p_da.h ../include/ph_rot.h ../include/phys_const.h ../include/track.h
quad.o: quad.cc ../include/Complex.h ../include/c_da.h ../include/dacpp.h ../include/element.h ../include/lin_map.h ../include/map_da.h ../include/map_double.h ../include/map_p_da.h ../include/matrix.h ../include/p_da.h ../include/phys_const.h ../include/quad.h ../include/track.h
sadp.o: sadp.cc ../include/Complex.h ../include/c_da.h ../include/dacpp.h ../include/element.h ../include/lie_da.h ../include/lin_map.h ../include/map_c_da.h ../include/map_da.h ../include/map_double.h ../include/map_p_da.h ../include/matrix.h ../include/p_da.h ../include/phys_const.h ../include/sadplus.h ../include/track.h
sext.o: sext.cc ../include/Complex.h ../include/c_da.h ../include/c_p_da.h ../include/dacpp.h ../include/element.h ../include/lin_map.h ../include/map_da.h ../include/map_double.h ../include/map_p_da.h ../include/matrix.h ../include/p_da.h ../include/phys_const.h ../include/sext.h ../include/track.h
solenoid.o: solenoid.cc ../include/Complex.h ../include/c_da.h ../include/dacpp.h ../include/element.h ../include/lin_map.h ../include/map_da.h ../include/map_double.h ../include/map_p_da.h ../include/matrix.h ../include/p_da.h ../include/phys_const.h ../include/solenoid.h ../include/track.h
teigen.o: teigen.cc ../include/Complex.h
thin.o: thin.cc ../include/Complex.h ../include/c_da.h ../include/c_p_da.h ../include/dacpp.h ../include/element.h ../include/lin_map.h ../include/map_da.h ../include/map_double.h ../include/map_p_da.h ../include/matrix.h ../include/p_da.h ../include/phys_const.h ../include/thin.h ../include/track.h
track.o: track.cc ../include/Complex.h ../include/c_da.h ../include/dacpp.h ../include/element.h ../include/lin_map.h ../include/map_da.h ../include/map_double.h ../include/map_p_da.h ../include/matrix.h ../include/p_da.h ../include/phys_const.h ../include/track.h
wig.o: wig.cc ../include/Complex.h ../include/c_da.h ../include/dacpp.h ../include/element.h ../include/lin_map.h ../include/map_da.h ../include/map_double.h ../include/map_p_da.h ../include/matrix.h ../include/p_da.h ../include/phys_const.h ../include/track.h ../include/wiggler.h

# End of File
