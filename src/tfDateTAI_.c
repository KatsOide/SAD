#include <feature.h>

#include <sim/sad_api.h>
#include <sim/TFCODE.h>
#include <sim/TFCBK.h>
#include <sim/TFSTK.h>

#include <libtai/leapsecs.h>
#include <libtai/caltime.h>
#include <libtai/taia.h>

#include <limits.h>
#include <string.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <math.h>

#ifndef	WITH_DATE_ARG
#define	DISABLE_DATE_ARG 1
#else
#define	DISABLE_DATE_ARG 0
#endif

#ifdef	WITHOUT_DATESTRING_ARG
#define	DISABLE_DATESTRING_ARG 1
#else
#define	DISABLE_DATESTRING_ARG 0
#endif

#define TICK_DIGITS	9
#define TICK_MIN	1e-9
#define TICK_MAX	1.0

#define DATECMD_FROMDATE		0
#define DATECMD_TODATE7			1
#define DATECMD_TODATE			2
#define DATECMD_DATE7			3
#define DATECMD_DATE			4
#define DATECMD_DAY			5
#define DATECMD_DAYNAME			6
#define DATECMD_TODATESTRING		7
#define DATECMD_DATESTRING		8
#define DATECMD_FROMDATESTRING		9
#define DATECMD_FROMDATESTRING2DATE7	10
#define DATECMD_FROMDATESTRING2DATE	11

#define DATEFMT_SAD_DAYNAME	-5
#define DATEFMT_SAD_DAY		-4
#define DATEFMT_SAD_LIST7	-3
#define DATEFMT_SAD_LIST	-2
#define DATEFMT_SAD_EPOCH	-1
#define DATEFMT_SAD		0
#define DATEFMT_ISO		1
#define DATEFMT_ISO_SHORT	2

#define DATEOPT_TIMEZONE	0x0001UL
#define DATEOPT_TICK		0x0002UL
#define DATEOPT_FORMAT		0x0004UL
#define DATEOPT_COMPRESSED	0x0008UL
#define DATEOPT_DISABLE_TICK	0x8000UL
typedef struct {
  double	tick;		/* DATEOPT_TICK */
  long		offset;		/* DATEOPT_TIMEZONE */
  int		format;		/* DATEOPT_FORMAT */
  bool		compress;	/* DATEOPT_COMPRESSED */
  bool		apply_tick;	/* for FormatDate() */
} sad_date_option_t;

typedef struct {
  char *name;
  long offset;
} sad_timezone_t;

typedef struct {
  char *name;
  int format;
} sad_date_format_t;

typedef struct tai tai_t;

static void tai_load(tai_t *t, double t0) {
  uint64 delta;
  double ti = floor(t0);

  t->x = 0;
  if(ti < 0) {
    delta = -ti;
    t->x -= delta;
  } else {
    delta =  ti;
    t->x += delta;
  }
}

static void tai_addi(tai_t *t, long dt) {
  uint64 delta;
  if(dt < 0) {
    delta = -dt;
    t->x -= delta;
  } else {
    delta =  dt;
    t->x += delta;
  }
}

static double tai_sub_approx(tai_t *t, tai_t *epoch) {
  tai_t tmp;

  if(tai_less(epoch, t)) {
    tai_sub(&tmp, t, epoch);
    return  tai_approx(&tmp);
  } else {
    tai_sub(&tmp, epoch, t);
    return -tai_approx(&tmp);
  }
}

static tai_t next_leapsecs_read;
static tai_t tai_sad_epoch;
static tai_t tai_unix_epoch;
static integer8 itfTimeZoneOffset = 0;

static void initialize_epoch(void) {
  struct caltime epoch;

  /* Initialize SAD epoch:	1900-01-01 00:00:00 +0900 */
  epoch.date.year	= 1900;
  epoch.date.month	=    1;
  epoch.date.day	=    1;
  epoch.hour		=    0;
  epoch.minute		=    0;
  epoch.second		=    0;
  epoch.offset		= +540;
  caltime_tai(&epoch, &tai_sad_epoch);

  /* Initialize UNIX epoch:	1970-01-01 00:00:00 +0000 */
  epoch.date.year	= 1970;
  epoch.date.month	=    1;
  epoch.date.day	=    1;
  epoch.hour		=    0;
  epoch.minute		=    0;
  epoch.second		=    0;
  epoch.offset		=   +0;
  caltime_tai(&epoch, &tai_unix_epoch);

  leapsecs_read();				/* Initialize leapsecs */
  tai_now(&next_leapsecs_read);
  tai_addi(&next_leapsecs_read, 2592000);	/* Add 30day offset */
}

static void leapsecs_check(struct tai *t) {
  struct tai now;

  if(t != NULL) {
    now = *t;
  } else {
    tai_now(&now);
  }
  if(tai_less(&next_leapsecs_read, &now) && leapsecs_read() != -1) {
    tai_addi(&now, 2592000);			/* Add 30day offset */
    next_leapsecs_read = now;
  }
}

static void caltime_shift_offset(struct caltime *ct, long offset,
				 int *weekday, int *yearday) {
  long s;

  s = ct->hour * 60 + ct->minute - ct->offset + offset;

  ct->minute = s % 60; s /= 60; if(ct->minute < 0) { ct->minute += 60; s -= 1; }
  ct->hour   = s % 24; s /= 24; if(ct->hour   < 0) { ct->hour   += 24; s -= 1; }
  ct->offset = offset;

  ct->date.day += s;

  s = caldate_mjd(&ct->date);
  caldate_frommjd(&ct->date, s, weekday, yearday);
}

static bool DecodeTimeZone(sad_date_option_t *opt, const char *offset) {
#include "timezone_table.h"
  char *endptr;
  long num;
  const sad_timezone_t *zone;

  if(offset == NULL || offset[0] == '\0') return false;

  if(offset[0] == '+' || offset[0] == '-') {
    if(strlen(offset) != 5) return false;
    num = strtol(offset + 1, &endptr, 10);
    if(*endptr != '\0') return false;
    if((num / 100) > 23 || (num % 100) > 59) return false;
    opt->offset = (num / 100 * 60) + (num % 100);
    if(offset[0] == '-') opt->offset *= -1;
    return true;
  }

  for(zone = sad_timezone_list; zone->name != NULL; zone++)
    if(strcmp(zone->name, offset) == 0) {
      opt->offset = zone->offset;
      return true;
    }

  return false;
}

static bool DecodeFormatOption(sad_date_option_t *opt, const char *format) {
  sad_date_format_t *node;
  static sad_date_format_t sad_date_format_list[] = {
    {"SAD",	DATEFMT_SAD},
    {"ISO",	DATEFMT_ISO},
    {"ISOs",	DATEFMT_ISO_SHORT},
    {NULL,	DATEFMT_SAD}};

  for(node = sad_date_format_list; node->name != NULL; node++)
    if(strcmp(node->name, format) == 0) {
      opt->format = node->format;
      return true;
    }

  return true;
}

static int DecodeOption(integer4 *isp1, integer4 *irtc, sad_date_option_t *opt,
			unsigned long active) {
  static const char *date_options[] = {
    "Compressed", "Format", "Tick", "TimeZone", NULL};
  const char **options;
  integer8 ia, kt;
  integer4 isp0, ispopt;
  int optc;

  if(itfTimeZoneOffset == 0)
    itfTimeZoneOffset = ktfsymbolz("TimeZoneOffset") - 4;

  /* Load standard value to date_option_t */
  opt->tick     = 0.0;
  opt->compress = false;
  opt->offset   = rlist(itfTimeZoneOffset);	/* TimeZoneOffset */
  opt->format	= DATEFMT_SAD;

  /* Load standard output format mode */
  opt->apply_tick	= (active & DATEOPT_TICK)
    && !(active & DATEOPT_DISABLE_TICK);

  /* Setup date_options skipping */
  options = date_options;
  if(!(active & DATEOPT_COMPRESSED)) {
    options += 1;
    if(!(active & DATEOPT_FORMAT)) {
      options += 1;
      if(!(active & DATEOPT_TICK)) {
	options += 1;
      }
    }
  }

  *irtc = 0;
  isp0 = isp;
  ispopt = itfgetoptionstk(*isp1, options);

  /*  fprintf(stderr, "DecodeOption %d,%d\n",ispopt,ispopt-isp0);*/
  if(ispopt < 0) {
    isp = isp0;
    return isp - *isp1;
  }

  optc = 0;
  while(options[optc] != NULL) optc += 1;

  if(active & DATEOPT_TIMEZONE)		/* DATEOPT_TIMEZONE */
    if((ktrmask & ktastk(isp0 + optc)) != ktfnr){
      opt->offset = rtastk(isp0 + optc - 0);}
    else{
      kt = ktfmask & ktastk(isp0 + optc -0);
      if(kt == ktfref){
      }
      else if(kt == ktfstring){
        ia = (ktamask & ktastk(isp0 + optc - 0));
#if SAD_REQUIRE_STRING_TERMINATION
        jlist(ilist(1, ia) + 1, ia + 1) = '\0';
#endif
        if(DecodeTimeZone(opt, &jlist(1, ia + 1))) {}
        else{
        *irtc = itfmessage(9, "General::wrongtype",
                           "\"TimeZone -> timezone name or offset from UTC\"");
        isp = isp0;
        return -1;}
      }
      else{
        *irtc = itfmessage(9, "General::wrongtype",
                           "\"TimeZone -> timezone name or offset from UTC\"");
        isp = isp0;
        return -1;
      }
    }

  if(active & DATEOPT_TICK)		/* DATEOPT_TICK */
    if((ktrmask & ktastk(isp0+optc-1)) != ktfnr){
      opt->tick = rtastk(isp0 + optc - 1);
      if(!(opt->tick > TICK_MIN)) opt->tick = TICK_MIN;
      if(!(TICK_MAX > opt->tick)) opt->tick = TICK_MAX;}
    else{
      kt = ktfmask & ktastk(isp0 + optc - 1);
      if(kt == ktfref){
      }
      else{
        *irtc = itfmessage(9, "General::wrongtype",
                           "\"Tick -> Real number\"");
        isp = isp0;
        return -1;
      }
    }

  if(active & DATEOPT_FORMAT)		/* DATEOPT_FORMAT */
    {
      kt = ktfmask & ktastk(isp0 + optc - 2);
      if(kt == ktfref){
      }
      else if (kt == ktfstring){
        ia = (ktamask & ktastk(isp0 + optc - 2));
#if SAD_REQUIRE_STRING_TERMINATION
        jlist(ilist(1, ia) + 1, ia + 1) = '\0';
#endif
        if(DecodeFormatOption(opt, &jlist(1, ia + 1))) {}
        else{
          *irtc = itfmessage(9, "General::wrongtype",
                           "\"Format -> DateTime string format name\"");
          isp = isp0;
          return -1;}
      }
      else{
      *irtc = itfmessage(9, "General::wrongtype",
			 "\"Format -> DateTime string format name\"");
      isp = isp0;
      return -1;
      }
    }

  if(active & DATEOPT_COMPRESSED)	/* DATEOPT_COMPRESSED */
    if((ktrmask & ktastk(isp0 + optc - 3)) != ktfnr){
      opt->compress = (rtastk(isp0 + optc - 3) != 0.0);
    }
    else{
      kt = ktfmask & ktastk(isp0 + optc - 3);
      if(kt == ktfref){
      }        
      else{
        *irtc = itfmessage(9, "General::wrongtype",
                           "\"Compressed -> True or False\"");
        isp = isp0;
        return -1;
      }
    }

  isp = isp0;
  return ispopt - *isp1 - 1;
}

static bool DecodeDateString(const sad_date_option_t *opt0, char *date,
			     tai_t *t0, double *frac0) {
  sad_date_option_t opt;
  struct caltime lt;
  char *beginptr, *endptr;
  long lfrac;
  double frac;
  bool success;
  int i;

  success = false;
  beginptr = date;
  frac = 0.0;
  lt.offset = opt0->offset;

  while(isspace(*beginptr)) beginptr += 1;

  switch(opt0->format) {
  case DATEFMT_ISO_SHORT:
  case DATEFMT_ISO:	/* "Y-m-d H:M:S.tick [TZ]" */
    if(!isdigit(beginptr[0])) break;
    lt.date.year  = strtol(beginptr, &endptr, 10);
    if(endptr[0] == '\0' || endptr[0] != '-') break; beginptr = endptr + 1;

    if(!isdigit(beginptr[0])) break;
    lt.date.month = strtol(beginptr, &endptr, 10);
    if(endptr[0] == '\0' || endptr[0] != '-') break; beginptr = endptr + 1;

    if(!isdigit(beginptr[0])) break;
    lt.date.day   = strtol(beginptr, &endptr, 10);
    if(endptr[0] == '\0' || !isspace(endptr[0])) break; beginptr = endptr + 1;

    while(isspace(*beginptr)) beginptr += 1;

    if(!isdigit(beginptr[0])) break;
    lt.hour       = strtol(beginptr, &endptr, 10);
    if(endptr[0] == '\0' || endptr[0] != ':') break; beginptr = endptr + 1;

    if(!isdigit(beginptr[0])) break;
    lt.minute     = strtol(beginptr, &endptr, 10);
    if(endptr[0] == '\0' || endptr[0] != ':') break; beginptr = endptr + 1;

    if(!isdigit(beginptr[0])) break;
    lt.second     = strtol(beginptr, &endptr, 10);
    if(endptr[0] == '.' && isdigit(endptr[1])) {
      beginptr = endptr + 1;
      lfrac = strtol(beginptr, &endptr, 10);
      for(i = 0, frac = lfrac; i < endptr - beginptr; i++) frac /= 10.0;
    }

    beginptr = endptr;
    while(isspace(*endptr)) endptr += 1;
    if(endptr[0] == '\0') {
      success = true;
    } else if(endptr - beginptr > 0) {
      if(DecodeTimeZone(&opt, endptr)) {
	lt.offset = opt.offset;
	success = true;
      }
    }
    break;

  case DATEFMT_SAD:	/* "m/d/Y H:M:S.tick" */
  default:
    lt.date.month = strtol(beginptr, &endptr, 10);
    if(endptr[0] == '\0' || endptr[0] != '/') break; beginptr = endptr + 1;

    lt.date.day   = strtol(beginptr, &endptr, 10);
    if(endptr[0] == '\0' || endptr[0] != '/') break; beginptr = endptr + 1;

    lt.date.year  = strtol(beginptr, &endptr, 10);
    if(endptr[0] == '\0' || !isspace(endptr[0])) break; beginptr = endptr + 1;

    lt.hour       = strtol(beginptr, &endptr, 10);
    if(endptr[0] == '\0' || endptr[0] != ':') break; beginptr = endptr + 1;

    lt.minute     = strtol(beginptr, &endptr, 10);
    if(endptr[0] == '\0' || endptr[0] != ':') break; beginptr = endptr + 1;

    lt.second     = strtol(beginptr, &endptr, 10);
    if(endptr[0] == '.' && isdigit(endptr[1])) {
      beginptr = endptr + 1;
      lfrac = strtol(beginptr, &endptr, 10);
      for(i = 0, frac = lfrac; i < endptr - beginptr; i++) frac /= 10.0;
    }

    while(isspace(*endptr)) endptr += 1;
    if(endptr[0] == '\0') success = true;
    break;
  }

  
  if(success) {
    leapsecs_check(NULL);
    caltime_tai(&lt, t0);
    *frac0 = frac;
  }

  return success;
}

static int FormatDate(integer8 *kx, 
		      integer4 *irtc,
		      tai_t *t, const double *frac0,
		      const sad_date_option_t *opt) {
  static const char *dayname[] = {"Sunday", "Monday", "Tuesday", "Wednesday",
				  "Thursday", "Friday", "Saturday"};
  struct caltime lt;
  char buffer0[64], buffer[64];
  int i, weekday;
  long abs_offset;
  double frac, tick;
  real8 vx;
  integer8 kax;

  tick = opt->tick;
  frac = !opt->apply_tick ? *frac0
    : (tick > 0.0) ? floor(*frac0 / tick) * tick : 0.0;

  switch(opt->format) {
  case DATEFMT_SAD_EPOCH:
    /*    fprintf(stderr, "format-SAD offset: %d",opt->offset);*/
    vx   = tai_sub_approx(t, &tai_sad_epoch) + frac;
    *kx = kfromr(vx);
    *irtc = 0;
    return 0;
    break;

  default:
    /*    fprintf(stderr, "format-default offset: %d ",opt->offset); */
    /* TAI -> UTC -> Local Time */
    caltime_utc(&lt, t, NULL, NULL);
    caltime_shift_offset(&lt, opt->offset, &weekday, NULL);
    break;
  }

  switch(opt->format) {
  case DATEFMT_SAD_DAY:
    vx   = caldate_mjd(&lt.date);
    *kx = kfromr(vx);
    *irtc = 0;
    return 0;
    break;

  case DATEFMT_SAD_DAYNAME:
    *kx  = ktfstring + ktsalocb(-1, dayname[weekday]);
    *irtc = 0;
    return 0;
    break;


  case DATEFMT_SAD_LIST7:
    kax  = ktavaloc(-1, 7);
    ilist(2, kax - 3) = lconstlist;
    rlist(kax + 1) = lt.date.year;
    rlist(kax + 2) = lt.date.month;
    rlist(kax + 3) = lt.date.day;
    rlist(kax + 4) = lt.hour;
    rlist(kax + 5) = lt.minute;
    rlist(kax + 6) = lt.second + frac;
    rlist(kax + 7) = lt.offset;
    *kx = ktflist + kax;
    *irtc = 0;
    return 0;
    break;

  case DATEFMT_SAD_LIST:
    kax  = ktavaloc(-1, 6);
    ilist(2, kax - 3) = lconstlist;
    rlist(kax + 1) = lt.date.year;
    rlist(kax + 2) = lt.date.month;
    rlist(kax + 3) = lt.date.day;
    rlist(kax + 4) = lt.hour;
    rlist(kax + 5) = lt.minute;
    rlist(kax + 6) = lt.second + frac;
    *kx   = ktflist + kax;
    *irtc = 0;
    return 0;
    break;

  case DATEFMT_ISO:
  case DATEFMT_ISO_SHORT:
    sprintf(buffer0, "%04ld-%02d-%02d %02d:%02d:%02d.%09ld",
	    lt.date.year, lt.date.month, lt.date.day,
	    lt.hour, lt.minute, lt.second,
	    (long)(frac / TICK_MIN));
    i = -1;
    if(tick > 0.0)
      for(i = 0; i < TICK_DIGITS && tick < 1.0; i++)
	tick *= 10.0;
    buffer0[strlen(buffer0) - TICK_DIGITS + i] = '\0';
    if(opt->format == DATEFMT_ISO_SHORT) {
      sprintf(buffer, "%s", buffer0);
    } else {
      abs_offset = opt->offset < 0 ? -opt->offset : opt->offset;
      sprintf(buffer, "%s %c%02ld%02ld", buffer0,
	      opt->offset < 0 ? '-' : '+',
	      abs_offset / 60, abs_offset % 60);
    }
    break;

  case DATEFMT_SAD:
  default:
    sprintf(buffer, opt->compress
	    ? "%d/%d/%ld %d:%d:%d.%09ld"
	    : "%02d/%02d/%04ld %02d:%02d:%02d.%09ld",
	    lt.date.month, lt.date.day, lt.date.year,
	    lt.hour, lt.minute, lt.second,
	    (long)(frac / TICK_MIN));
    i = -1;
    if(tick > 0.0)
      for(i = 0; i < TICK_DIGITS && tick < 1.0; i++)
	tick *= 10.0;
    buffer[strlen(buffer) - TICK_DIGITS + i] = '\0';
    break;
  }

  *kx  = ktfstring + ktsalocb(-1, buffer);
  *irtc = 0;
  return 0;
}

static int DateCommand(integer4 *isp1,
		       integer8 *kx,
		       integer4 *irtc, int mode) {
  integer8 ia, kt;
  struct caltime lt;
  struct taia now;
  tai_t t, tmp;
  double frac;
  sad_date_option_t opt;
  bool without_argv, from_date_string;
  unsigned long dateopt;
  const char *msg_narg, *msg_wrongtype;

  without_argv = false;
  from_date_string = false;
  dateopt = DATEOPT_TICK | DATEOPT_TIMEZONE;
  switch(mode) {
  case DATECMD_DATE:
  case DATECMD_DATE7:
#if DISABLE_DATE_ARG
    without_argv = true;
#endif
  case DATECMD_TODATE:
  case DATECMD_TODATE7:
    break;

  case DATECMD_FROMDATE:
    dateopt |= DATEOPT_DISABLE_TICK;
    break;

  case DATECMD_DATESTRING:
#if DISABLE_DATESTRING_ARG
    without_argv = true;
#endif
  case DATECMD_TODATESTRING:
    dateopt |= DATEOPT_COMPRESSED | DATEOPT_FORMAT;
    break;

  case DATECMD_FROMDATESTRING:
  case DATECMD_FROMDATESTRING2DATE7:
  case DATECMD_FROMDATESTRING2DATE:
    from_date_string = true;
    dateopt &= !DATEOPT_TICK;
    dateopt |= DATEOPT_FORMAT;
    break;

  case DATECMD_DAY:
  case DATECMD_DAYNAME:
    dateopt &= !DATEOPT_TICK;
    break;

  default:
    break;
  }

  if(from_date_string) {
    msg_narg = "\"1 (+ option)\"";
    msg_wrongtype = "\"DateTime string\"";
  } else {
    msg_narg = without_argv ? "\"0 (+ option)\"" : "\"0 or 1 (+ option)\"";
    msg_wrongtype = "\"Null or Real number or {Y,m,d[,H[,M[,S[,offset]]]]}\"";
  }

  switch(DecodeOption(isp1, irtc, &opt, dateopt)) {
  case -1:
    return -1;
    break;

  case 0:
    if(from_date_string) {
      *irtc = itfmessage(9, "General::narg", msg_narg);
      return -1;
    }
    taia_now(&now);	/* Get current TAIA time */
    taia_tai(&now, &t); frac = taia_frac(&now);
    opt.apply_tick = true;
    if(mode != DATECMD_FROMDATE) leapsecs_check(&t);
    break;

  case 1:

    kt = (ktfmask & ktastk(*isp1 + 1));
    if( kt == ktfoper) {
      if((ktamask & ktastk(*isp1 + 1)) != mtfnull) {
	*irtc = without_argv ? itfmessage(9, "General::narg", msg_narg)
	  : itfmessage(9, "General::wrongtype", msg_wrongtype);
	return -1;
      }
      if(from_date_string) {
	*irtc = itfmessage(9, "General::wrongtype", msg_wrongtype);
	return -1;
      }
      taia_now(&now);	/* Get current TAIA time */
      taia_tai(&now, &t); frac = taia_frac(&now);
      opt.apply_tick = true;
      if(mode != DATECMD_FROMDATE) leapsecs_check(&t);
    }
    else if(kt == ktflist){
      if(without_argv) {
        *irtc = itfmessage(9, "General::narg", msg_narg);
        return -1;
      }
      if(from_date_string) {
        *irtc = itfmessage(9, "General::wrongtype", msg_wrongtype);
        return -1;
      }
      ia = ktfaddr(ktastk(*isp1 + 1));
      if((klist(ia) == ktfoper + mtflist)
         && ktfreallistq(ia)) {
#define GET_FRACTIONAL(store_i, store_f, source)        \
        store_i  = (source);                            \
	store_f  = (source) - store_i;                  \
	store_i += (int)floor(store_f);                 \
	store_f -=      floor(store_f)

	double fractional[5] = {0, 0, 0, 0, 0};
	long source_offset;
	lt.offset       = opt.offset;
	lt.second       = 0; frac = 0.0;
	lt.minute       = 0;
	lt.hour         = 0;
	switch(ilist(2, ia - 1)) {
	case 7:
	  lt.offset     = rlist(ia + 7);
	case 6:
	  GET_FRACTIONAL(lt.second,     frac,          rlist(ia + 6));
	case 5:
	  GET_FRACTIONAL(lt.minute,     fractional[4], rlist(ia + 5));
	case 4:
	  GET_FRACTIONAL(lt.hour,       fractional[3], rlist(ia + 4));
	case 3:
	  GET_FRACTIONAL(lt.date.day,   fractional[2], rlist(ia + 3));
	  GET_FRACTIONAL(lt.date.month, fractional[1], rlist(ia + 2));
	  GET_FRACTIONAL(lt.date.year,  fractional[0], rlist(ia + 1));
          break;
	default:
	  *irtc = itfmessage(9, "General::wrongtype", msg_wrongtype);
	  return -1;
	  break;
	}
	leapsecs_check(NULL);

	source_offset = lt.offset;
	caltime_tai(&lt, &t);

	/* Normalize fractional part */
	for(int i = 0; i < 5; i++)
	  if(fractional[i] > 0) {
	    caltime_utc(&lt, &t, NULL, NULL);
	    caltime_shift_offset(&lt, source_offset, NULL, NULL);
	    switch(i) {
	    case 0:
	      lt.date.year   += 1;
	      break;
	    case 1:
	      lt.date.month  += 1;
	      break;
	    case 2:
	      lt.date.day    += 1;
	      break;
	    case 3:
	      lt.hour        += 1;
	      break;
	    case 4:
	      lt.minute      += 1;
	      break;
	    }
	    caltime_tai(&lt, &tmp);
	    frac += tai_sub_approx(&tmp, &t) * fractional[i];
	    if(frac > 1) {
	      tai_addi(&t, (long)floor(frac)); frac -= floor(frac);
	    }
	  }
        break;
      }
      *irtc = itfmessage(9, "General::wrongtype", msg_wrongtype);
      return -1;
    }
    else if(kt == ktfstring){
      if(without_argv) {
	*irtc = itfmessage(9, "General::narg", msg_narg);
	return -1;
      }
      if(!from_date_string) {
	*irtc = itfmessage(9, "General::wrongtype", msg_wrongtype);
	return -1;
      }
      ia = ktfaddr(ktastk(*isp1 + 1));
#if SAD_REQUIRE_STRING_TERMINATION
      jlist(ilist(1, ia) + 1, ia + 1) = '\0';
#endif
      if(!DecodeDateString(&opt, &jlist(1, ia + 1), &t, &frac)) {
	switch(opt.format) {
	case DATEFMT_ISO_SHORT:
	case DATEFMT_ISO:
	  msg_wrongtype = "\"''CCYY-mm-dd HH:MM:SS[.TICKS] [TimeZone]''\"";
	  break;

	case DATEFMT_SAD:
	  msg_wrongtype = "\"''mm/dd/CCYY HH:MM:SS[.TICKS]''\"";
	  break;

	default:
	  break;
	}
	*irtc = itfmessage(9, "General::wrongtype",
			   msg_wrongtype);
	return -1;
      }
    }
    else if((ktrmask & ktastk(*isp1+1)) != ktfnr){
      if(without_argv) {
        *irtc = itfmessage(9, "General::narg", msg_narg);
        return -1;
      }
      if(from_date_string) {
        *irtc = itfmessage(9, "General::wrongtype", msg_wrongtype);
        return -1;
      }
      tai_load(&t, rtastk(*isp1 + 1));
      tai_add(&t, &t, &tai_sad_epoch);	/* Add SAD epoch */
      frac = rtastk(*isp1 + 1) - tai_sub_approx(&t, &tai_sad_epoch);
      tai_addi(&t, (long)floor(frac)); frac -= floor(frac);
      leapsecs_check(NULL);}
    else {
      *irtc = without_argv ? itfmessage(9, "General::narg", msg_narg)
        : itfmessage(9, "General::wrongtype", msg_wrongtype);
      return -1;
    }
    break;
  default:
    *irtc = itfmessage(9, "General::narg", msg_narg);
    return -1;
  }

  switch(mode) {
  case DATECMD_DATE7:
  case DATECMD_TODATE7:
  case DATECMD_FROMDATESTRING2DATE7:
    opt.format = DATEFMT_SAD_LIST7;
    break;

  case DATECMD_DATE:
  case DATECMD_TODATE:
  case DATECMD_FROMDATESTRING2DATE:
    opt.format = DATEFMT_SAD_LIST;
    break;

  case DATECMD_FROMDATE:
  case DATECMD_FROMDATESTRING:
    opt.format = DATEFMT_SAD_EPOCH;
    break;

  case DATECMD_DAY:
    opt.format = DATEFMT_SAD_DAY;
    break;

  case DATECMD_DAYNAME:
    opt.format = DATEFMT_SAD_DAYNAME;
    break;

  default:
    break;
  }

  return FormatDate(kx, irtc, &t, &frac, &opt);
}

static int Date7(integer4 *isp1,
		 integer8 *kx,
		 integer4 *irtc) {
  return DateCommand(isp1, kx, irtc, DATECMD_DATE7);
}

static int Date(integer4 *isp1,
		integer8 *kx,
		integer4 *irtc) {
  return DateCommand(isp1, kx, irtc, DATECMD_DATE);
}

static int ToDate7(integer4 *isp1,
		   integer8 *kx,
		   integer4 *irtc) {
  return DateCommand(isp1, kx, irtc, DATECMD_TODATE7);
}

static int ToDate(integer4 *isp1,
		  integer8 *kx,
		  integer4 *irtc) {
  return DateCommand(isp1, kx, irtc, DATECMD_TODATE);
}

static int FromDate(integer4 *isp1,
		    integer8 *kx,
		    integer4 *irtc) {
  return DateCommand(isp1, kx, irtc, DATECMD_FROMDATE);
}

static int DateString(integer4 *isp1,
		      integer8 *kx,
		      integer4 *irtc) {
  return DateCommand(isp1, kx, irtc, DATECMD_DATESTRING);
}

static int ToDateString(integer4 *isp1,
			integer8 *kx,
			integer4 *irtc) {
  return DateCommand(isp1, kx, irtc, DATECMD_TODATESTRING);
}

static int JulianDay(integer4 *isp1,
		     integer8 *kx,
		     integer4 *irtc) {
  return DateCommand(isp1, kx, irtc, DATECMD_DAY);
}

static int DayName(integer4 *isp1,
		   integer8 *kx,
		   integer4 *irtc) {
  return DateCommand(isp1, kx, irtc, DATECMD_DAYNAME);
}

static int FromDateString(integer4 *isp1,
			  integer8 *kx,
			  integer4 *irtc) {
  return DateCommand(isp1, kx, irtc, DATECMD_FROMDATESTRING);
}

static int FromDateString2Date7(integer4 *isp1,
				integer8 *kx,
				integer4 *irtc) {
  return DateCommand(isp1, kx, irtc, DATECMD_FROMDATESTRING2DATE7);
}

static int FromDateString2Date(integer4 *isp1,
			       integer8 *kx,
			       integer4 *irtc) {
  return DateCommand(isp1, kx, irtc, DATECMD_FROMDATESTRING2DATE);
}

static int SetTimeZone(integer4 *isp1,
		       integer8 *kx,
		       integer4 *irtc) {
  integer8 ka;
  sad_date_option_t opt;

  if(itfTimeZoneOffset == 0)
    itfTimeZoneOffset = ktfsymbolz("TimeZoneOffset") - 4;

  if(isp != *isp1 + 1) {
    if(isp == *isp1) {
      *kx   = klist(itfTimeZoneOffset);
      *irtc=0;
      return 0;}
    *irtc = itfmessage(9, "General::narg", "\"1\"");
    return -1;
  }

  if(ktastk(*isp1 + 1) == ktfoper){
    *kx   = klist(itfTimeZoneOffset);
    *irtc=0;
    return 0;}
  if((ktrmask & ktastk(*isp1 + 1)) != ktfnr) {
    rlist(itfTimeZoneOffset) = rtastk( *isp1 + 1);}
  else if((ktfmask & ktastk(*isp1 + 1)) == ktfstring){
      ka = (ktamask & ktastk(*isp1 + 1));
#if SAD_REQUIRE_STRING_TERMINATION
      jlist(ilist(1, ka) + 1, ka + 1) = '\0';
#endif
      if(DecodeTimeZone(&opt, &jlist(1, ka + 1))) {
        rlist(itfTimeZoneOffset) = opt.offset;
      }
    else {
      *irtc = itfmessage(9, "General::wrongtype",
                         "\"Timezone name or offset from UTC\"");
      return -1;
      }
  }

  *kx   = klist(itfTimeZoneOffset);
  *irtc = 0;
  return 0;
}

static int GetTimeOfDay(integer4 *isp1,
			integer8 *kx,
			integer4 *irtc) {
  struct taia now;
  tai_t t;
  double frac;
  real8 vx;

  if(isp != *isp1 + 1
     || ktastk(isp) != ktfoper + mtfnull) {
    *irtc = itfmessage(9, "General::narg", "\"0\"");
    return -1;
  }

  taia_now(&now);	/* Get current TAIA time */
  taia_tai(&now, &t); frac = taia_frac(&now);

  vx   = tai_sub_approx(&t, &tai_unix_epoch) + frac;
  *kx = (integer8) vx;
  *irtc = 0;
  return 0;
}

/* SADScript function registration of Random stuff */
#define REG	dlfunaloc
#define REG8	dlfunaloc8
#ifdef WITH_EXTENSION_MODULE
int dldeffun_(void) {
#else
int sadDefFunc_DateTAI(void) {
#endif /* WITH_EXTENSION_MODULE */
  if(!feature_provide("Date/TAI", FEATURE_VERSION(1, 0), FEATURE_VERSION(1, 0)))
    return -1;

  initialize_epoch();

#ifdef EXPORT_TAI_SYMBOL
  REG8("TaiGetTimeOfDay",	GetTimeOfDay,	0, NULL, NULL, 0);
  REG8("TaiFromDate",		FromDate,	1, NULL, NULL, 0);
  REG8("TaiToDate$",		ToDate7,	1, NULL, NULL, 0);
  REG8("TaiToDate",		ToDate,		1, NULL, NULL, 0);
  REG8("TaiDate$",		Date7,		1, NULL, NULL, 0);
  REG8("TaiDate",		Date,		1, NULL, NULL, 0);
  REG8("TaiDateString",		DateString,	1, NULL, NULL, 0);
  REG8("TaiToDateString",	ToDateString,	1, NULL, NULL, 0);
  REG8("TaiDay",			DayName,	1, NULL, NULL, 0);
  REG8("TaiFromDateString",	FromDateString,	1, NULL, NULL, 0);
#endif
#if !defined(WITH_EXTENSION_MODULE) || defined(OVERRIDE_DATE)
  REG8("GetTimeOfDay",		GetTimeOfDay,	0, NULL, NULL, 0);
  REG8("FromDate",		FromDate,	1, NULL, NULL, 0);
  REG8("ToDate$",		ToDate7,	1, NULL, NULL, 0);
  REG8("ToDate",			ToDate,		1, NULL, NULL, 0);
  REG8("Date$",			Date7,		1, NULL, NULL, 0);
  REG8("Date",			Date,		1, NULL, NULL, 0);
  REG8("DateString$",		DateString,	1, NULL, NULL, 0);
  REG8("ToDateString",		ToDateString,	1, NULL, NULL, 0);
  REG8("Day",			DayName,	1, NULL, NULL, 0);
  REG8("FromDateString",		FromDateString,	1, NULL, NULL, 0);
  REG8("FromDateString2Date$",	FromDateString2Date7, 1, NULL, NULL, 0);
  REG8("FromDateString2Date",	FromDateString2Date,  1, NULL, NULL, 0);
#endif
  REG8("JulianDay",		JulianDay,	1, NULL, NULL, 0);
  REG8("SetTimeZone",		SetTimeZone,	1, NULL, NULL, 0);

  return 0;
}

/* End of File */
