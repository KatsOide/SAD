# Don't edit Makefile! Use conf-* for configuration.

SHELL=/bin/sh

default: it

auto-ccld.sh: \
conf-cc conf-ld warn-auto.sh
	( cat warn-auto.sh; \
	echo CC=\'`head -1 conf-cc`\'; \
	echo LD=\'`head -1 conf-ld`\' \
	) > auto-ccld.sh

hasclock_gettime.h: \
tryclock_gettime.c compile load
	( ( ./compile tryclock_gettime.c && ./load tryclock_gettime ) >/dev/null \
	2>&1 \
	&& echo \#define HASCLOCK_GETTIME 1 || exit 0 ) > hasclock_gettime.h
	rm -f tryclock_gettime.o

caldate.0: \
caldate.3
	nroff -man caldate.3 > caldate.0

caldate_fmjd.o: \
compile caldate_fmjd.c caldate.h
	./compile caldate_fmjd.c

caldate_fmt.o: \
compile caldate_fmt.c caldate.h
	./compile caldate_fmt.c

caldate_mjd.0: \
caldate_mjd.3
	nroff -man caldate_mjd.3 > caldate_mjd.0

caldate_mjd.o: \
compile caldate_mjd.c caldate.h
	./compile caldate_mjd.c

caldate_norm.o: \
compile caldate_norm.c caldate.h
	./compile caldate_norm.c

caldate_scan.o: \
compile caldate_scan.c caldate.h
	./compile caldate_scan.c

caldate_ster.o: \
compile caldate_ster.c caldate.h
	./compile caldate_ster.c

caltime.0: \
caltime.3
	nroff -man caltime.3 > caltime.0

caltime_fmt.o: \
compile caltime_fmt.c caldate.h caltime.h
	./compile caltime_fmt.c

caltime_scan.o: \
compile caltime_scan.c caltime.h caldate.h
	./compile caltime_scan.c

caltime_tai.0: \
caltime_tai.3
	nroff -man caltime_tai.3 > caltime_tai.0

caltime_tai.o: \
compile caltime_tai.c tai.h uint64.h leapsecs.h caldate.h caltime.h
	./compile caltime_tai.c

caltime_utc.o: \
compile caltime_utc.c tai.h uint64.h leapsecs.h caldate.h caltime.h
	./compile caltime_utc.c

check: \
load check.o libtai.a
	./load check libtai.a 

check.o: \
compile check.c tai.h uint64.h leapsecs.h caltime.h caldate.h
	./compile check.c

compile: \
make-compile warn-auto.sh systype
	( cat warn-auto.sh; ./make-compile "`cat systype`" ) > \
	compile
	chmod 755 compile

easter: \
load easter.o libtai.a
	./load easter libtai.a 

easter.o: \
compile easter.c caldate.h
	./compile easter.c

find-systype: \
find-systype.sh auto-ccld.sh
	cat auto-ccld.sh find-systype.sh > find-systype
	chmod 755 find-systype

it: \
man libtai.a leapsecs check easter yearcal nowutc

leapsecs: \
load leapsecs.o libtai.a
	./load leapsecs libtai.a 

leapsecs.0: \
leapsecs.3
	nroff -man leapsecs.3 > leapsecs.0

leapsecs.o: \
compile leapsecs.c tai.h uint64.h leapsecs.h caldate.h
	./compile leapsecs.c

leapsecs_add.o: \
compile leapsecs_add.c leapsecs.h tai.h uint64.h
	./compile leapsecs_add.c

leapsecs_init.o: \
compile leapsecs_init.c leapsecs.h
	./compile leapsecs_init.c

leapsecs_read.o: \
compile leapsecs_read.c tai.h uint64.h leapsecs.h
	./compile leapsecs_read.c

leapsecs_sub.o: \
compile leapsecs_sub.c leapsecs.h tai.h uint64.h
	./compile leapsecs_sub.c

libtai.a: \
makelib tai_add.o tai_now.o tai_pack.o tai_sub.o tai_unpack.o \
taia_add.o taia_approx.o taia_fmtfrac.o taia_frac.o taia_half.o \
taia_less.o taia_now.o taia_pack.o taia_sub.o taia_tai.o \
taia_unpack.o caldate_fmt.o caldate_scan.o caldate_fmjd.o \
caldate_mjd.o caldate_norm.o caldate_ster.o leapsecs_read.o \
leapsecs_init.o leapsecs_add.o leapsecs_sub.o caltime_fmt.o \
caltime_scan.o caltime_tai.o caltime_utc.o
	./makelib libtai.a tai_add.o tai_now.o tai_pack.o \
	tai_sub.o tai_unpack.o taia_add.o taia_approx.o \
	taia_fmtfrac.o taia_frac.o taia_half.o taia_less.o \
	taia_now.o taia_pack.o taia_sub.o taia_tai.o taia_unpack.o \
	caldate_fmt.o caldate_scan.o caldate_fmjd.o caldate_mjd.o \
	caldate_norm.o caldate_ster.o leapsecs_read.o \
	leapsecs_init.o leapsecs_add.o leapsecs_sub.o caltime_fmt.o \
	caltime_scan.o caltime_tai.o caltime_utc.o

load: \
make-load warn-auto.sh systype
	( cat warn-auto.sh; ./make-load "`cat systype`" ) > load
	chmod 755 load

make-compile: \
make-compile.sh auto-ccld.sh
	cat auto-ccld.sh make-compile.sh > make-compile
	chmod 755 make-compile

make-load: \
make-load.sh auto-ccld.sh
	cat auto-ccld.sh make-load.sh > make-load
	chmod 755 make-load

make-makelib: \
make-makelib.sh auto-ccld.sh
	cat auto-ccld.sh make-makelib.sh > make-makelib
	chmod 755 make-makelib

makelib: \
make-makelib warn-auto.sh systype
	( cat warn-auto.sh; ./make-makelib "`cat systype`" ) > \
	makelib
	chmod 755 makelib

man: \
tai.0 tai_now.0 tai_pack.0 taia.0 taia_now.0 taia_pack.0 caldate.0 \
caldate_mjd.0 caltime.0 caltime_tai.0 leapsecs.0

nowutc: \
load nowutc.o libtai.a
	./load nowutc libtai.a 

nowutc.o: \
compile nowutc.c leapsecs.h tai.h uint64.h taia.h caltime.h caldate.h
	./compile nowutc.c

shar: \
FILES BLURB README INSTALL TODO THANKS CHANGES FILES VERSION Makefile \
leapsecs.txt leapsecs.dat leapsecs.c check.in check.out check.c \
easter.c yearcal.c nowutc.c caldate.3 caldate_mjd.3 caldate.h \
caldate_fmt.c caldate_scan.c caldate_mjd.c caldate_fmjd.c \
caldate_norm.c caldate_ster.c tai.3 tai_now.3 tai_pack.3 tai.h \
tai_add.c tai_now.c tai_pack.c tai_sub.c tai_unpack.c leapsecs.3 \
leapsecs.h leapsecs_add.c leapsecs_sub.c leapsecs_init.c \
leapsecs_read.c caltime.3 caltime_tai.3 caltime.h caltime_fmt.c \
caltime_scan.c caltime_tai.c caltime_utc.c taia.3 taia_now.3 \
taia_pack.3 taia.h taia_add.c taia_approx.c taia_fmtfrac.c \
taia_frac.c taia_half.c taia_less.c taia_now.c taia_pack.c taia_sub.c \
taia_tai.c taia_unpack.c uint64.h conf-cc conf-ld find-systype.sh \
make-compile.sh make-load.sh make-makelib.sh trycpp.c warn-auto.sh
	shar -m `cat FILES` > shar
	chmod 400 shar

systype: \
find-systype trycpp.c
	./find-systype > systype

tai.0: \
tai.3
	nroff -man tai.3 > tai.0

tai_add.o: \
compile tai_add.c tai.h uint64.h
	./compile tai_add.c

tai_now.0: \
tai_now.3
	nroff -man tai_now.3 > tai_now.0

tai_now.o: \
compile tai_now.c tai.h uint64.h
	./compile tai_now.c

tai_pack.0: \
tai_pack.3
	nroff -man tai_pack.3 > tai_pack.0

tai_pack.o: \
compile tai_pack.c tai.h uint64.h
	./compile tai_pack.c

tai_sub.o: \
compile tai_sub.c tai.h uint64.h
	./compile tai_sub.c

tai_unpack.o: \
compile tai_unpack.c tai.h uint64.h
	./compile tai_unpack.c

taia.0: \
taia.3
	nroff -man taia.3 > taia.0

taia_add.o: \
compile taia_add.c taia.h tai.h uint64.h
	./compile taia_add.c

taia_approx.o: \
compile taia_approx.c taia.h tai.h uint64.h
	./compile taia_approx.c

taia_fmtfrac.o: \
compile taia_fmtfrac.c taia.h tai.h uint64.h
	./compile taia_fmtfrac.c

taia_frac.o: \
compile taia_frac.c taia.h tai.h uint64.h
	./compile taia_frac.c

taia_half.o: \
compile taia_half.c taia.h tai.h uint64.h
	./compile taia_half.c

taia_less.o: \
compile taia_less.c taia.h tai.h uint64.h
	./compile taia_less.c

taia_now.0: \
taia_now.3
	nroff -man taia_now.3 > taia_now.0

taia_now.o: \
compile taia_now.c taia.h tai.h uint64.h hasclock_gettime.h
	./compile taia_now.c

taia_pack.0: \
taia_pack.3
	nroff -man taia_pack.3 > taia_pack.0

taia_pack.o: \
compile taia_pack.c taia.h tai.h uint64.h
	./compile taia_pack.c

taia_sub.o: \
compile taia_sub.c taia.h tai.h uint64.h
	./compile taia_sub.c

taia_tai.o: \
compile taia_tai.c taia.h tai.h uint64.h
	./compile taia_tai.c

taia_unpack.o: \
compile taia_unpack.c taia.h tai.h uint64.h
	./compile taia_unpack.c

yearcal: \
load yearcal.o libtai.a
	./load yearcal libtai.a 

yearcal.o: \
compile yearcal.c caldate.h
	./compile yearcal.c
