#include <sim/sad_tcltk.h>
#include <sim/sad_api.h>
#include <sim/TFCODE.h>

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <inttypes.h>

#include <cadef.h>

#define MAXCOMLEN 1024
#define CA_PEND_IO_TIME	0.0001

static int bDebug = 0;
static int bInCaPend = 0;

static real8 *cavalue_buffer = NULL;
static size_t cavalue_bufsz = 0;

static void cavalue_alloc(size_t n0) {
  const size_t page_size = 4096 / sizeof(real8);
  real8 *ptr;
  size_t n;

  if(n0 > cavalue_bufsz) {
    if(n0 > page_size) {
      n = page_size * ((n0 + 2048) / page_size + 1);
    } else {
      n = 1;
      while(n0 > n) n *= 2;
    }

    ptr = realloc(cavalue_buffer, sizeof(real8) * n);
    if(ptr == NULL) {
      fprintf(stderr, "cavalue: Cannot allocate memory(cavalue_buffer)\n");
      exit(1);
    }
    cavalue_buffer = ptr;
    cavalue_bufsz = n;
  }
}

void cavalue(struct event_handler_args arg)
{
  struct dbr_time_double *p = (struct dbr_time_double *)arg.dbr;
  int n = arg.count;
  real8 chid = pointer2double(arg.chid);
  integer4 st = p->status, sev = p->severity;
  integer4 type = arg.type - DBR_TIME_STRING;
  real8 x,t = p->stamp.secPastEpoch + p->stamp.nsec * 1e-9;

  integer4 length, mode = -1;
  integer8 *i8ptr = NULL;
  int i;

  dbr_string_t *sv;
  dbr_char_t *cv;
  dbr_short_t *shv;
  dbr_long_t *lv;
  dbr_float_t *fv;

  if(bDebug) {
    printf("cavalue[%s,%d]: %d\n", ca_name(arg.chid), n, ca_state(arg.chid));
    fflush(stdout);
  }

  switch(arg.type) {
  case DBR_TIME_STRING:
    if(bDebug) printf("string[%s]\n",
		      ((struct dbr_time_string *)arg.dbr)->value);
    sv = &(((struct dbr_time_string *)arg.dbr)->value);
    cavalue_alloc(n); 
    i8ptr = (integer8 *)cavalue_buffer;
    for(i = 0; i < n; i++) {
      length = 0;
      while(length < sizeof(dbr_string_t) && sv[i][length] != '\0')
	length += 1;
      i8ptr[i]   = ktfstring + ktsalocb_(&mode, sv[i], &length, length);
    }
    break;

  case DBR_TIME_CHAR:
    if(bDebug) printf("char[%c]\n",
		      ((struct dbr_time_char *)arg.dbr)->value);
    cv = &(((struct dbr_time_char *)(arg.dbr))->value);
    cavalue_alloc(n); i8ptr = (integer8 *)cavalue_buffer;
    for(i = 0; i < n; i++) {x=(real8) cv[i];i8ptr[i] = kfromr(x);}
    break;

  case DBR_TIME_SHORT:
    if(bDebug) printf("short[%" PRId16 "]\n",
		      ((struct dbr_time_short *)arg.dbr)->value);
    shv = &(((struct dbr_time_short *)(arg.dbr))->value);
    cavalue_alloc(n); i8ptr = (integer8 *)cavalue_buffer;
    for(i = 0; i < n; i++) {x=(real8) shv[i];i8ptr[i] = kfromr(x);}
    break;

  case DBR_TIME_LONG:
    if(bDebug) printf("long[%" PRId32 "]\n",
		      ((struct dbr_time_long *)arg.dbr)->value);
    lv = &(((struct dbr_time_long *)(arg.dbr))->value);
    cavalue_alloc(n); i8ptr = (integer8 *)cavalue_buffer;
    for(i = 0; i < n; i++) {x=(real8) lv[i];i8ptr[i] = kfromr(x);}
    break;

  case DBR_TIME_FLOAT:
    if(bDebug) printf("float[%f]\n",
		      ((struct dbr_time_float *)arg.dbr)->value);
    fv = &(((struct dbr_time_float *)(arg.dbr))->value);
    cavalue_alloc(n); i8ptr = (integer8 *)cavalue_buffer;
    for(i = 0; i < n; i++) {x=(real8)fv[i];i8ptr[i] = kfromr(x);}
    break;

  case DBR_TIME_DOUBLE:
    if(bDebug) printf("double[%22.16le]\n",
		      ((struct dbr_time_double *)arg.dbr)->value);
    i8ptr = &kfromr(((struct dbr_time_double *)(arg.dbr))->value);
    break;

  default:
    printf("cavalue: unknown type: %ld\n", arg.type);
    return;
  }

  tfepicsvaluecb_(&chid, &st, &sev, &t, &type, &n, i8ptr);
  /* tfcavaluecb_(&chid, &st, &sev, &t, &type, &n, i4prt, r8ptr); */
}

void cachange(struct connection_handler_args arg)
{
  real8 chan_id;
  integer4 status;

  chan_id = pointer2double(arg.chid);
  status = ca_state(arg.chid);

  if(bDebug) printf("cachange: ID=%p[%s] %d\n",
		    (void*)arg.chid, ca_name(arg.chid), status);
  
  tfepicsconstatcb_(&chan_id, &status);
}
/*    Oide Change  5/8/1999
      void cachange(struct connection_handler_args arg)
      {
      integer4 len, itx, iax, irtc;
      real8 vx;
      int status = ca_state(arg.chid);
      char str[MAXCOMLEN];
      size_t full_len;

      if(bDebug) printf("cachange: ID=%p[%s]: %d\n",
      arg.chid, ca_name(arg.chid), status);

      full_len = snprintf(str, MAXCOMLEN, "EPICS$ConStatCB[%td,%d]",
      (ptrdiff_t)arg.chid, status);
      len = strlen(str);
      if(full_len > len) {
      fprintf(stderr, "cachange[%s]: "
      "SAD evaluation buffer overflowed by %dbytes. "
      "Change event is discarded\n",
      ca_name(arg.chid), full_len - len);
      fflush(stderr);
      } else
      tfevalb_(str, &len, &itx, &iax, &vx, &irtc, MAXCOMLEN);
      }
*/

void caputcb(struct event_handler_args arg)
{
  integer8 kx;
  integer4 len, irtc;
  char str[MAXCOMLEN];
  size_t full_len;

  if(bDebug) printf("caputcb: ID=%p[%s]\n",
		    (void*)arg.chid, ca_name(arg.chid));

  full_len = snprintf(str, MAXCOMLEN, "EPICS$PutCB[%td]", (ptrdiff_t)arg.chid);
  len = strlen(str);
  if(full_len > len) {
    fprintf(stderr, "caputcb[%s]: "
	    "SAD evaluation buffer overflowed by %zubytes. "
	    "Callback event is discarded\n",
	    ca_name(arg.chid), full_len - len);
    fflush(stderr);
  } else
    tfevalb_(str, len, &kx, &irtc);
}

integer ecaaddarrayevent_(real8 *chd,
			  real8 *evd, integer4 *vtype, integer4 *evt)
{
  int status,vt;
  chid chan_id = (chid)double2pointer(*chd);
  evid event_id;
  chtype cht;

  if (*vtype>=0)
    vt = *vtype;
  else
    vt = ca_field_type(chan_id);
#ifndef NDEBUG
  /*  fprintf(stderr, "%s : field type id %d",,vt);*/
#endif
  switch(vt) {
  case DBF_STRING:
    cht = DBR_TIME_STRING;
    break;
  case DBF_CHAR:
    cht = DBR_TIME_CHAR;
    break;
  case DBF_SHORT:
    cht = DBR_TIME_SHORT;
    break;
  case DBF_LONG:
  case DBF_ENUM:
    cht = DBR_TIME_LONG;
    break;
  case DBF_FLOAT:
    cht = DBR_TIME_FLOAT;
    break;
  case DBF_DOUBLE:
    cht = DBR_TIME_DOUBLE;
    break;
  default:
    printf("ecaaddarrayevent: unknown ca field type: %d\n",ca_field_type(chan_id));
    /*cht = DBR_TIME_DOUBLE;*/
    return -1;
  }
  status = ca_add_masked_array_event(cht, 0, chan_id, cavalue, NULL,
				     0, 0, 0, &event_id, *evt);
  SEVCHK(status, "error in ca_add_masked_array_event");

  *evd = pointer2double(event_id);
  return status;
}

integer ecasearchandconnect_(const char *chname, real8 *chd, void *pilist)
{
  chid chan_id;
  int status;

  if (bDebug)
    printf("casearch: [%s]\n", chname);
  status = ca_search_and_connect(chname, &chan_id, cachange, pilist);
  SEVCHK(status,"error in ca_search");
  if (status!=ECA_NORMAL) {
    *chd = 0;
    return status;
  }
  *chd = pointer2double(chan_id);
  return status;
}

void ecaclearchannel_(real8 *chd)
{
  chid chan_id = (chid)double2pointer(*chd);
  int status;

  status = ca_clear_channel(chan_id);
  SEVCHK(status, "ca_clear_channel");
}

void ecaclearevent_(real8 *evd)
{
  evid event_id = (evid)double2pointer(*evd);
  int status;

  status = ca_clear_event(event_id);
  SEVCHK(status, "ca_clear_event");
}

void ecaput_(real8 *ch_id, logical4 *cb,
	     integer4 *mode, integer4 *nc, void *val, integer4 *irtc)
{
  chid chan_id = (chid)double2pointer(*ch_id);
  integer8 kx, *ia;
  integer4 i, len;
  void *ptr = NULL;
  dbr_string_t *sv = NULL;
  dbr_double_t *dv = NULL;
  int type, status, size = *nc > 0 ? *nc : 1;

  switch(*mode) {
  case 0: /* STRING */
    type = DBR_STRING;
    break;

  case 1: /* DOUBLE */
    type = DBR_DOUBLE;
    break;

  default:
    type = DBR_DOUBLE;
    break;
  }

  switch(type) {
  case DBR_DOUBLE:
    if(*nc > 0) { /* Vector API */
      dv = malloc(sizeof(dbr_double_t) * *nc);
      if(dv != NULL) {
	*irtc = 0;
	ia = (integer8 *)val;
	for(i = 1; i <= *nc; i++) {
          kx = klist(*ia + i);
	  if(ktfnonrealq(kx)) {
	    *irtc = -1;
	    break;
	  }
	  dv[i - 1] = rfromk(kx);
	}
      }
    } else { /* Scalar API */
      dv = malloc(sizeof(dbr_double_t));
      if(dv != NULL) {
	*irtc = 0;
	*dv = *((real8 *)val);
      }
    }
    ptr = dv;
    break;

  case DBR_STRING:
    if(*nc > 0) { /* Vector API */
      sv = malloc(sizeof(dbr_string_t) * *nc + 1);
      if(sv != NULL) {
	*irtc = 0;
	ia = (integer8 *)val;
	for(i = 1; i <= *nc; i++) {
          kx = klist(*ia + i);
	  if(ktfnonstringq(kx)) {
	    *irtc = -1;
	    break;
	  }
	  tfgetstr_((character)sv[i - 1], sizeof(dbr_string_t), &kx, &len);
	  ((char *)sv[i - 1])[len] = '\0';
	}
      }
    } else { /* Scalar API */
      sv = malloc(sizeof(dbr_string_t) + 1);
      if(sv != NULL) {
	*irtc = 0;
	ia = (integer8 *)val;
	tfgetstr_((character)sv[0], sizeof(dbr_string_t), ia, &len);
	((char *)sv[0])[len] = '\0';
      }
    }
    ptr = sv;
    break;

  default:
    ptr = NULL;
    *irtc = -1;
  }

  /* Check channel connection if without callback mode */
  if(*irtc == 0 && !(*cb) && ca_state(chan_id) != cs_conn) *irtc = 1;

  if(*irtc == 0) {
    if(bDebug) {
      printf("ecaput%s[%d]: ", *cb ? "cb" : "", size);
      switch(type) {
      case DBR_DOUBLE:
	printf("%lf\n", dv[0]);	break;
      case DBR_STRING:
	printf("%s\n",  sv[0]);	break;
      default:
	printf("Unkown type[%d]\n", type);
      }
    }
    if(*cb)
      status = ca_array_put_callback(type, size, chan_id, ptr, caputcb, NULL);
    else
      status = ca_array_put(type, size, chan_id, ptr);
    SEVCHK(status, "error in ca_array_put");
    status = ca_flush_io();
    SEVCHK(status, "ca_flush_io");
  }

  if(ptr != NULL) free(ptr);
}

integer ecapendio_(real4 *timeout)
{
  int status;

  if (bDebug)
	printf("pend_io\n");
  if (bInCaPend) {
	fprintf(stderr, "CaSad: try to call ca_pend_io in ca_pend_* !!");
	return 0;
  }
  bInCaPend = 1;
  status = ca_pend_io(*timeout);
  bInCaPend = 0;
  SEVCHK(status,"ca_pend_io");

  return status;
}

integer ecapendevent_(real4 *timeout)
{
  int status;

  if (bDebug)
	printf("pend_event\n");
  if (bInCaPend) {
	fprintf(stderr, "CaSad: try to call ca_pend_event in ca_pend_* !!");
	return 0;
  }
  bInCaPend = 1;
  status = ca_pend_event(*timeout);
  bInCaPend = 0;
  SEVCHK(status,"ca_pend_event");

  return status;
}

integer ecaflushio_()
{
  int status;

  status = ca_flush_io();
  SEVCHK(status,"ca_flush_io");

  return status;
}

static void ecafdcallback(void *clientData, int mask)
{
  if (bDebug)
	printf("pend_event_callback 1\n");
  if (bInCaPend) {
	fprintf(stderr, "CaSad: try to call ca_pend_event in ca_pend_* !!");
	return;
  }
  bInCaPend = 1;
  ca_pend_event(CA_PEND_IO_TIME);
  bInCaPend = 0;
  if (bDebug)
	printf("pend_event_callback 2\n");
}

static void ecafdregister(userarg, fd, opened)
void *userarg;
int	fd,opened;
{
  if (bDebug)
	printf("ecafdregister %d %d %p\n",fd,opened,userarg);
  if(opened)
	sadTk_CreateFileHandler(fd, SAD_TK_READABLE | SAD_TK_EXCEPTION,
				ecafdcallback, userarg);
  else
	sadTk_DeleteFileHandler(fd);
}

integer ecainit_()
{
  int status;

  if (bDebug)
	printf("ecainit\n");
  status = ca_task_initialize();
  SEVCHK(status,"ca_task_initialize");
  status = ca_add_fd_registration(ecafdregister, NULL);
  SEVCHK(status,"ca_add_fd_registration");

  return status;
}

void ecahostname_(integer8 *kx)
{
  chid chan_id = (chid)double2pointer(*kx);
  integer4 length, mode = -1;
  char *hostname;

  hostname = (char *)ca_host_name(chan_id);
  if(hostname != NULL) {
    *kx = ktfstring + ktsalocb_(&mode, hostname, &length, length);
  } else {
    *kx = ktfoper + mtfnull;
  }
}

integer ecafieldtype_(real8 *chd, integer4 *type)
{
  chid chan_id = (chid)double2pointer(*chd);

  *type = ca_field_type(chan_id);

  return 0;
}

integer ecaelementcount_(real8 *chd, integer4 *count)
{
  chid chan_id = (chid)double2pointer(*chd);

  *count = ca_element_count(chan_id);

  return 0;
}

integer4 ecadebugprint_(integer4 *bNewDebug)
{
  integer4 bOldDebug = bDebug;

  bDebug = *bNewDebug;
  return bOldDebug;
}
