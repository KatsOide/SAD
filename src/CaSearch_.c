#include <sim/sad_api.h>
#include <sim/TFCODE.h>

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <cadef.h>

#define PEND_IO_TIME_OUT 100
#define CA_PEND_EVENT_TIME	1e-6		/* formerly 0.0001    */
#define MAX_CA_TRY 1

/* semaphore */
static int _ca_semid=0;
#define LOCK  {if(_ca_semid >0) --_ca_semid;}
#define UNLOCK {++_ca_semid;}
#define LOCKED (_ca_semid == 0)

/*  for Tk */
#ifdef TK
#include <tk.h>
#endif

#ifndef FALSE
#define FALSE 0
#endif

#ifndef TRUE
#define TRUE 1
#endif

static int inited=FALSE;

#ifdef TK
static void ca_service();
static void ca_fd_register(void *pfdctx, int fd, int condition);
#endif

static integer4 NTFREAL   = ntfreal;
static integer4 NTFSTRING = ntfstring;
static integer4 NTFLIST   = ntflist;

void cainit_(void)
{


#ifdef TK
  SEVCHK(ca_task_initialize(), 
	"init_ca: C.A. initialize failure.\n");

  SEVCHK(ca_add_fd_registration(ca_fd_register, NULL),
	"init_ca:   Fd registration failed.\n");
#endif

  /* create semaphore */
  UNLOCK;/* initialize semaphore( _ca_semid) */

  inited=TRUE;
}

void casearch_(const char *pname, integer4 *nc,
	       real8 *ch_id, integer4 *irtc)
{	
  chid  chan_id;
  int   status;
  char *chname;
  
  *irtc = -1;
  if(!inited) cainit_();

  chname = malloc(*nc+1);
  if(chname == NULL) {
    fprintf(stderr, "Cannot allocate memory in caSearch\n");
    exit(1);
  }
  memcpy(chname, pname, *nc);
  chname[*nc] = '\0';
  LOCK
    status = ca_search(chname, &chan_id);
  UNLOCK;
  SEVCHK(status, "CaSearch");
  if(status == ECA_NORMAL) { 
    *ch_id = pointer2double(chan_id);
    *irtc = 0; /* normal completion */
  } else {
    *irtc = - CA_EXTRACT_MSG_NO(status);
  }
  free(chname);
  return;
}

void caconstatus_(real8 *ch_id, integer4 *ival)
{
  chid chan_id = (chid)double2pointer(*ch_id);
  enum channel_state chs;
  
  chs = ca_state(chan_id);
  *ival = chs - cs_conn;
  /* *ival : -2->never_conn, -1-> prev_conn, 0 -> connected ,1 -> closed */ 
}

void caclose_(real8 *ch_id, integer4 *irtc)
{
  chid chan_id = (chid)double2pointer(*ch_id);
  int count, status;

  *irtc = -1;
  count = 0;
  status = ca_clear_channel(chan_id);
  do{
    LOCK{
      status = ca_pend_io(PEND_IO_TIME_OUT);
    }  UNLOCK;
  }while(status != 0 && count++ < MAX_CA_TRY);
  
  *irtc = status - 1;
  /* *irtc : status == 1 -> 0, otherwise -> (status - 1) */
}

void capendio_(real8 *t, integer4 *irtc)
{
  int status;

  *irtc = -1;
  LOCK
    status = ca_pend_io(*t);
  UNLOCK;
  SEVCHK(status,"pend io");

  switch(status) {
  case 1:
    *irtc = 0;
    break;
  case 0:
    *irtc = -1;
    break;
  default:
    *irtc = status;
    break;
  }
}

void cagettype_(real8 *ch_id, integer4 *type, integer4 *count)
{
  chid chan_id = (chid)double2pointer(*ch_id);

  assert(count);
  assert(type);

  *count = ca_element_count(chan_id);
  *type = ca_field_type(chan_id);
  return;
}

void cagetreal_(real8 *ch_id, integer8 *kax, integer4 *irtc)
{
  chid chan_id = (chid)double2pointer(*ch_id);
  int count, status=-1;
  integer4 i;
  real8 ts, st, sv;
  struct dbr_time_double tvalue;

  if(ca_state(chan_id) != cs_conn) {
    *irtc = -1;
    return;
  }
  count = 0;
  do{
    LOCK{
      status = ca_get(DBR_TIME_DOUBLE, chan_id, &tvalue);
      status |=  ca_pend_io(PEND_IO_TIME_OUT);
    } UNLOCK;
  }while(status !=0 && count++ < MAX_CA_TRY);
  
  SEVCHK(status, "ca_get in CaGetReal");
  if(status == 1) {
    st = (double) tvalue.status;
    sv = (double) tvalue.severity;
    ts = (double) tvalue.stamp.secPastEpoch
      + ((double) tvalue.stamp.nsec) / 1e9;
    rlist(*kax + 1) = tvalue.value;
    rlist(*kax + 2) = st;
    rlist(*kax + 3) = sv;
    rlist(*kax + 4) = ts;
    /*
    i=1; tfsetlist_(&NTFREAL, 0, &tvalue.value, iax, &i);
    i=2; tfsetlist_(&NTFREAL, 0, &st, iax, &i);
    i=3; tfsetlist_(&NTFREAL, 0, &sv, iax, &i);
    i=4; tfsetlist_(&NTFREAL, 0, &ts, iax, &i);
    */
    *irtc = 0;
  } else
    *irtc = (status == 1 ? 0 : status);
  return;
}

void cagetstring_(real8 *ch_id, integer8 *kax, integer4 *irtc)
{
  chid chan_id = (chid)double2pointer(*ch_id);
  int status=-1;
  integer4 mode = 0, i, nc, ivx;
  real8 ts, st, sv;
  struct dbr_time_string tvalue;
  dbr_string_t *pv;

  if(ca_state(chan_id) != cs_conn) {
    *irtc = -1;
    return;
  }
  LOCK {
    status = ca_get(DBR_TIME_STRING, chan_id, &tvalue);
    SEVCHK(status, "ca_get");
    status |= ca_pend_io(PEND_IO_TIME_OUT);
  } UNLOCK;
  SEVCHK(status, "pend after get");

  st = (double) tvalue.status;
  sv = (double) tvalue.severity;
  ts = (double) tvalue.stamp.secPastEpoch
    + ((double) tvalue.stamp.nsec) / 1e9;
  pv = (dbr_string_t *) &(tvalue.value);
  nc = 0;
  while(nc < sizeof(dbr_string_t) && (*pv)[nc] != '\0') nc += 1;
  klist(*kax+1)=ktfstring + ktsalocb(0, *pv);
  rlist(*kax+2)=st;
  rlist(*kax+3)=sv;
  rlist(*kax+4)=ts;
  /*
  ivx = itsalocb_(&mode, *pv, &nc, nc);
  i=1; tfsetlist_(&NTFSTRING, &ivx,   0, iax, &i);
  i=2; tfsetlist_(&NTFREAL,      0, &st, iax, &i);
  i=3; tfsetlist_(&NTFREAL,      0, &sv, iax, &i);
  i=4; tfsetlist_(&NTFREAL,      0, &ts, iax, &i);*/
  *irtc = (status == 1 ? 0 : status);
  return;
}

void cagetrealarray_(real8 *ch_id, integer8 *kax, integer4 *irtc)
{
  chid chan_id = (chid)double2pointer(*ch_id);
  int count, status=-1;
  integer8 kvx;
  integer4 i;
  real8 ts, st, sv, vx;
  dbr_double_t *pv;
  struct dbr_time_double *pval;

  count = ca_element_count(chan_id);
  pval = malloc(sizeof(struct dbr_time_double)
		+ (count - 1) * sizeof(dbr_double_t));

  LOCK
    status = ca_array_get(DBR_TIME_DOUBLE, count, chan_id, pval);
  UNLOCK;
  SEVCHK(status, "ca_get");
  if(status != 1) {*irtc = -1; free(pval); return;}
  LOCK
    status = ca_pend_io(PEND_IO_TIME_OUT);
  UNLOCK;
  SEVCHK(status," pend after get");
  if(status != 1) {*irtc = -1; free(pval); return;}

  kvx = ktavaloc(0, count);
  st = (double) pval->status;
  sv = (double) pval->severity;
  ts = (double) pval->stamp.secPastEpoch
    + ((double) pval->stamp.nsec) / 1e9;
  for(i = 1, pv = (dbr_double_t *)&(pval->value); i <= count; i++) {
    klist(kvx+i) = pv[i - 1];
  }
  klist(*kax+1)=ktflist+kvx;
  rlist(*kax+2)=st;
  rlist(*kax+3)=sv;
  rlist(*kax+4)=ts;
  *irtc = (status == 1 ? 0 : status);
  free(pval);
  return;
}

void cagetstringarray_(real8 *ch_id,
		       integer8 *kax, integer4 *irtc)
{
  chid chan_id = (chid)double2pointer(*ch_id);
  int count, status=-1;
  integer8 kas, kvx;
  integer4 mode = 0, i, nc;
  real8 ts, st, sv;
  struct dbr_time_string *rval;
  dbr_string_t *pv;

  count = ca_element_count(chan_id);
  rval = malloc(sizeof(struct dbr_time_string)
		+ (count - 1) * sizeof(dbr_string_t));

  LOCK
    status = ca_array_get(DBR_TIME_STRING, count, chan_id, rval);
  UNLOCK;
  SEVCHK(status, "ca_get");
  if(status != 1) {*irtc=-1; free(rval); 
    printf("count: %d\n",count);
    return;}
  LOCK
    status = ca_pend_io(PEND_IO_TIME_OUT);
  UNLOCK;
  SEVCHK(status," pend after get");
  if(status != 1) {*irtc=-1; free(rval); return;};

  kvx = ktadaloc(0, count);
  st = (double) rval->status;
  sv = (double) rval->severity;
  ts = (double) rval->stamp.secPastEpoch
    + ((double) rval->stamp.nsec) / 1e9;
  for(i = 1, pv = (dbr_string_t *)&(rval->value); i <= count; i++) {
    nc = 0;
    while(nc < sizeof(dbr_string_t) && pv[i-1][nc] != '\0') nc += 1;
    klist(kvx+i) = ktfstring + ktsalocb(0, pv[i-1]);
  }
  klist(*kax+1)=ktflist+kvx;
  rlist(*kax+2)=st;
  rlist(*kax+3)=sv;
  rlist(*kax+4)=ts;
  *irtc = (status == 1 ? 0 : status);
  free(rval);
  return;
}

void caputreal_(real8 *ch_id, real8 *val, integer4 *irtc)
{
  chid chan_id = (chid)double2pointer(*ch_id);
  dbr_double_t fval = *val;

  LOCK{
    SEVCHK(ca_put(DBR_DOUBLE, chan_id, &fval),"ca_put");
    SEVCHK(ca_pend_io(PEND_IO_TIME_OUT), "pend after put");
    SEVCHK(ca_get(DBR_DOUBLE, chan_id, &fval), "ca_get");
    SEVCHK(ca_pend_io(PEND_IO_TIME_OUT)," pend after get");
  }UNLOCK;
  *val = fval;
  *irtc = 0;
  return;
}

void caputstring_(real8 *ch_id, const_character *str, integer4 *nc,
		  integer4 *irtc, ftnlen sl)
{
  chid chan_id = (chid)double2pointer(*ch_id);
  dbr_string_t pchr;

  memcpy(pchr, str, *nc > sizeof pchr ? sizeof pchr : *nc);
  if(*nc < sizeof pchr) pchr[*nc] = '\0';
  LOCK{
    SEVCHK(ca_put(DBR_STRING, chan_id, pchr), "ca_put");
    SEVCHK(ca_pend_io(0.1), "pend after put");
    SEVCHK(ca_get(DBR_STRING, chan_id, pchr), "ca_get");
    SEVCHK(ca_pend_io(0.1), "pend after get");
  }UNLOCK;
  memcpy(str, pchr, *nc > sizeof pchr ? sizeof pchr : *nc);
  *irtc = 0;
  return;
}

void caputrealarray_(real8 *ch_id,
		     integer4 *size, integer8 *ka, integer4 *irtc)
{
  chid chan_id = (chid)double2pointer(*ch_id);
  int status;
  integer4 i;
  real8 vx;
  dbr_double_t *val;

  val = malloc(sizeof(dbr_double_t) * *size);
  if(val == NULL) {
    exit(1);
  }

  for(i = 1; i <= *size ; i++) {
    val[i - 1] = rlist(*ka +i);
    if(ktfnonrealq((integer8) val[i - 1])) {
      *irtc = -1;
      free(val);
      return;
    }
  }
  LOCK{
    status = ca_array_put(DBR_DOUBLE, *size, chan_id, val);
    SEVCHK(status, "ca_array_put");
    SEVCHK(ca_pend_io(0.5), "pend after put");
    status = ca_array_get(DBR_DOUBLE, *size, chan_id, val);
    SEVCHK(status, "ca_array_get");
    SEVCHK(ca_pend_io(0.5)," pend after get");
  }UNLOCK;
  for(i = 1; i <= *size; i++) {
    rlist(*ka + i)=val[i - 1];
  }
  free(val);
  *irtc = 0;
  return;
}

void caputstringarray_(real8 *ch_id,
		       integer4 *size, integer8 *ka, integer4 *irtc)
{
  chid chan_id = (chid)double2pointer(*ch_id);
  int status;
  integer8 kx,kax;
  integer4 mode = 0, nc, i;
  real8 vx;
  dbr_string_t *val;

  val = malloc(sizeof(dbr_string_t) * *size + 1);
  if(val == NULL) {
    exit(1);
  }

  for(i = 1; i <= *size; i++) {
    kx = klist(*ka + i);
    if(ktfnonstringq(kx)) {
      *irtc = -1;
      free(val);
      return;
    } else {
      kax=ktfaddr(kx);
      tfgetstr_((char *)val[i-1], sizeof(dbr_string_t), &kax, &nc);
      ((char *)val[i-1])[nc] = '\0';
    }
  }
  LOCK{
    status = ca_array_put(DBR_STRING, *size, chan_id, val);
    SEVCHK(status, "ca_array_put");
    SEVCHK(ca_pend_io(0.5), "pend after put");
    status = ca_array_get(DBR_STRING, *size, chan_id, val);
    SEVCHK(status, "ca_array_get");
    SEVCHK(ca_pend_io(0.5)," pend after get");
  }UNLOCK;
  for(i = 1; i <= *size; i++) {
    nc = 0;
    while(nc < sizeof(dbr_string_t) && val[i-1][nc] != '\0') nc += 1;
    klist(*ka + i) = ktfstring + ktsalocb_(0, val[i-1], &nc, nc);
  }
  free(val);
  *irtc = 0;
  return;
}

void cagetrealnowait_(real8 *ch_id, integer8 *ka, integer4 *irtc)
{
  chid chan_id = (chid)double2pointer(*ch_id);
  int count, type, status = -1;
  integer4 i;
  real8 v;
  struct dbr_time_double *ptvalue;

  count = ca_element_count(chan_id);
  type = ca_field_type(chan_id);
  ptvalue = malloc(sizeof(struct dbr_time_double)
		   + (count - 1) * sizeof(dbr_double_t));

  LOCK
    status = ca_array_get(DBR_TIME_DOUBLE,count, chan_id, ptvalue);
  UNLOCK;
  SEVCHK(status, "ca_get in CaGetRealNoWait");

  klist(*ka+1)=0;
  if(status == 1) {
    rlist(*ka + 2)=pointer2double(ptvalue);;
    *irtc = 0;
  } else {
    rlist(*ka + 2) = 0.0;
    *irtc = status - 1;
  }
  rlist(*ka + 3)=*ch_id;
  klist(*ka + 4)=0;
  return;
}

void cagetrealfinish_(integer8 *ka, integer4 *irtc)
{
  int count;
  integer8 kv, kav;
  integer4 i;
  real8 v, ts, st, sv;
  struct dbr_time_double *ptvalue;
  dbr_double_t *pv;
  kv = klist(*ka + 2);
  if(ktfnonrealq(kv) || kv < 0) {
    *irtc = -1;
    return;
  }
  v = rlist(*ka + 2);
  ptvalue = (struct dbr_time_double *)double2pointer(v);
  assert(ptvalue);

  kv = klist(*ka + 3);
  if(ktfnonrealq(kv)) {
    fprintf(stderr, " invalid frame\n");
  }
  v=rlist(*ka + 3);
  {
    chid chan_id = (chid)double2pointer(v);
    count = ca_element_count(chan_id);
  }

  if( ktfnonrealq(kv) || count <= 0 ) {
    *irtc = -1;
    return;
  }

  if(count == 1) {
    v = (dbr_double_t) ptvalue->value;
    rlist(*ka + 1) = v;
  } else if(count > 1) {
    kav = ktavaloc(0 , count);
    if(kav <= 0) {
      exit(kav);
    }
    for(i = 1, pv = &(ptvalue->value); i <= count; i++) {
      v = pv[i - 1];
      rlist(kav + i)=v;
    }
    klist(*ka + 1)=ktflist + kav;
    ilist(2, *ka-3) = lnonreallist;
  }
  *irtc=0;
  st = (double) ptvalue->status;
  sv = (double) ptvalue->severity;
  ts = (double) ptvalue->stamp.secPastEpoch
    + ((double) ptvalue->stamp.nsec) / 1e9;
  rlist(*ka+2)=st;
  rlist(*ka+3)=sv;
  rlist(*ka+4)=ts;
  free(ptvalue);
  return;
}

void cagetstringnowait_(real8 *ch_id, integer8 *ka, integer4 *irtc)
{
  chid chan_id = (chid)double2pointer(*ch_id);
  int count, type, status = -1;
  integer4 i;
  real8 v;
  struct dbr_time_string *ptvalue;

  count = ca_element_count(chan_id);
  type = ca_field_type(chan_id);
  ptvalue = malloc(sizeof(struct dbr_time_string)
		   + (count - 1) * sizeof(dbr_string_t));
  LOCK
    status=ca_array_get(DBR_TIME_STRING,count, chan_id, ptvalue);
  UNLOCK;
  SEVCHK(status,"ca_get in CaGetStringNoWait");
  klist(*ka + 1)=0;
  if(status == 1) {
    rlist(*ka + 2) = pointer2double(ptvalue);
    *irtc = 0;
  } else {
    klist(*ka + 2) = 0;
    *irtc = status - 1;
  }
  rlist(*ka + 3)=*ch_id;
  klist(*ka + 4)=1;
  return;
}

void cagetstringfinish_(integer8 *ka, integer4 *irtc)
{
  int count;
  integer8 kv, kav;
  integer4 mode = 0, i, iv, ivx, nc;
  real8 v, ts, st, sv;
  struct dbr_time_string *ptvalue;
  dbr_string_t *pv;
  
  kv = klist(*ka + 2);
  if(ktfnonrealq(kv) || kv < 0) {
    *irtc = -1;
    return;
  }
  
  ptvalue = (struct dbr_time_string *)double2pointer(rlist(*ka +2));
  assert(ptvalue);

  kv = klist(*ka + 3);
  if(ktfnonrealq(kv)) {
    fprintf(stderr, " invalid frame\n");
  }

  {
    chid chan_id = (chid)double2pointer(rlist(*ka + 3));
    count = ca_element_count(chan_id);
  }

  if(ktfnonrealq(kv) || count <= 0) {
    *irtc = -1;
    return;
  }

  if(count == 1) {
    pv = &(ptvalue->value);
    nc = 0;
    while(nc < sizeof(dbr_string_t) && pv[nc] != '\0') nc += 1;
    klist(*ka + 1) = ktfstring + ktsalocb(0, *pv);
  } 
  else if(count > 1) {
    kav = ktadaloc(0 , count);
    if(kav <= 0) {
      exit(kav);
    }
    
    for(i = 1, pv = &(ptvalue->value); i <= count; i++) {
      nc = 0;
      while(nc < sizeof(dbr_string_t) && pv[i-1][nc] != '\0') nc += 1;
      i=0;
      klist(kav + i) = ktfstring + ktsalocb_(&i, pv[i-1], &nc, nc);
    }
    klist(*ka + 1)=ktflist + kav;
  }
  *irtc = 0;
  st = (double) ptvalue->status;
  sv = (double) ptvalue->severity;
  ts = (double) ptvalue->stamp.secPastEpoch
    + ((double) ptvalue->stamp.nsec) / 1e9;
  rlist(*ka + 2) = st;
  rlist(*ka + 3) = sv;
  rlist(*ka + 4) = ts;
  free(ptvalue);
  return;
}

void caputrealnowait_(real8 *ch_id, real8 *val, integer4 *irtc)
{
  chid chan_id = (chid)double2pointer(*ch_id);
  dbr_double_t fval = *val;

  LOCK
    SEVCHK(ca_put(DBR_DOUBLE, chan_id, &fval),"ca_put");
  UNLOCK;
  *irtc = 0;
  return;
}

void caputstringnowait_(real8 *ch_id, const_character *str, integer4 *nc,
			integer4 *irtc, ftnlen sl)
{
  chid chan_id = (chid)double2pointer(*ch_id);
  dbr_string_t pchr;

  memcpy(pchr, str, *nc > sizeof pchr ? sizeof pchr : *nc);
  if(*nc < sizeof pchr) pchr[*nc] = '\0';
  LOCK
    SEVCHK(ca_put(DBR_STRING, chan_id, pchr), "ca_put");
  UNLOCK;
  *irtc = 0;
  return;
}

void caputrealarraynowait_(real8 *ch_id,
			   integer4 *size, integer8 *ka, integer4 *irtc)
{
  chid chan_id = (chid)double2pointer(*ch_id);
  integer4 i;
  real8 vx;
  int status;
  dbr_double_t *val;

  val = malloc(sizeof(dbr_double_t) * *size);
  if(val == NULL) {
    exit(1);
  }

  for(i = 1; i <= *size; i++) {
    if(ktfnonrealq(klist(*ka + i))) {
      *irtc = -1;
      free(val);
      return;
    };
    val[i - 1] = rlist(*ka + i);
  }
  status = ca_array_put(DBR_DOUBLE, *size, chan_id, val);
  SEVCHK(status, "ca_array_put");
  free(val);
  *irtc = 0;
  return;
}

void caputstringarraynowait_(real8 *ch_id,
			     integer4 *size, integer8 *ka, integer4 *irtc)
{
  chid chan_id = (chid)double2pointer(*ch_id);
  integer8 kx, kax;
  integer4 i, nc;
  real8 vx;
  int status;
  dbr_string_t *val;

  val = malloc(sizeof(dbr_string_t) * *size + 1);
  if(val == NULL) {
    exit(1);
  }

  for(i = 1; i <= *size; i++) {
    kx = klist(*ka + i);
    if(ktfnonstringq(kx)) {
      *irtc = -1;
      free(val);
      return;
    } else {
      kax=ktfaddr(kx);
      tfgetstr_((char *)val[i-1], sizeof(dbr_string_t), &kax, &nc);
      ((char *)val[i-1])[nc] = '\0';
    }
  }
  status = ca_array_put(DBR_STRING, *size, chan_id, val);
  SEVCHK(status, "ca_array_put");
  free(val);
  *irtc = 0;
  return;
}

#ifdef TK
static void ca_fd_register(pfdctx, fd, condition)
void	*pfdctx;
int	fd;
int	condition;
{
	if(condition) {
		Tk_CreateFileHandler(fd, TK_READABLE, ca_service, pfdctx);
	} else {
		Tk_DeleteFileHandler(fd);
	}
}

/*
******************************************************
* NAME
*	ca_service()
* DESCRIPTION
*	1) call ca_pend_event to allow event handler execution
******************************************************
*/

static void ca_service()
{
#ifdef DEBUG
	struct timeval tp,tp1;
 	struct timezone tzp;

  	gettimeofday(&tp1, &tzp);
#endif
	LOCK
	  ca_pend_event(CA_PEND_EVENT_TIME);
	UNLOCK;
#ifdef DEBUG
 	gettimeofday(&tp1, &tzp);
	printf(" %d %d\n", tp.tv_sec, tp.tv_usec);
  	printf(" %d %d\n", tp1.tv_sec, tp1.tv_usec);
#endif
}

#endif
