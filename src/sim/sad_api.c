/* Extended Allocation API */
#include <sim/sad_api.h>
#include <sim/sad_memory.h>

#include <sim/TFCODE.h>
#include <sim/TFCBK.h>
#include <sim/TFSTK.h>

#include <sys/types.h>
#include <stdbool.h>
#include <stdarg.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <pwd.h>

/* asprintf() for SAD */
int sad_asprintf(char **ptr, const char *format, ...) {
  int count;
  va_list va0, va;

  va_start(va0, format);

  va_copy(va, va0);
  count = vsnprintf(NULL, 0, format, va);
  if(count <= 0) {
    *ptr = NULL;
    return count;
  }

  *ptr = malloc(count + 1);
  count = vsnprintf(*ptr, count + 1, format, va0);

  va_end(va0);

  return count;
}

/* Expand `~' prefix by shell rule and return new allocated string */
char* expand_tilde(const char *path) {
  int n;
  char *user, *home, *buf;
  struct passwd *pe;

  if(path[0] != '~')  return strdup(path);

  if(path[1] == '/' || path[1] == '\0') {
    home = getenv("HOME"); n = 0;
    if(home == NULL) {
      pe = getpwuid(getuid());
      if(pe == NULL) return NULL;
      home = pe->pw_dir;
    }
  } else {
    for(n = 1; path[1 + n] != '/' && path[1 + n] != '\0'; n++);

    user = malloc(n + 1);
    if(user == NULL) return NULL;
    memcpy(user, path + 1, n); user[n] = '\0';
    pe = getpwnam(user); free(user);
    if(pe == NULL) return NULL;
    home = pe->pw_dir;
  }

  if(home == NULL) return NULL;

  buf = malloc(strlen(home) + strlen(path) - n);
  strcpy(buf, home);
  strcat(buf, path + 1 + n);

  return buf;
}

/* Fortran2008 LEN_TRIM() for C */
integer len_trim(const_character str, ftnlen slen) {
  for(int i = slen - 1; i >= 0; i--)
    if(str[i] != ' ')
      return i + 1;

  return 0;
}

integer8 ktsalocbref_(integer4 *mode, const char **str, integer4 *length) {
  return ktsalocb_(mode, *str, length, *length);
}

integer8 ktsalocbcstrs_(integer4 *mode, const char **str, integer4 *index) {
  const char *string = str[*index];
  integer4 length;

  length = strlen(string);
  return ktsalocb_(mode, string, &length, length);
}

/* Wrapper function for SAD internal API */
void tfreadbuf(integer4 lfn, integer4 *ib, integer4 *nc) {
   __readbuf_MOD_tfreadbuf(&lfn, ib, nc);
  }

void trbinit(integer4 lfn, integer4 *ib) {
   __tfrbuf_MOD_trbinit(&lfn, ib);
  }

integer8 kfromr_(real8 *x){
  return kfromr(*x);
}

real8 rfromk_(integer8 *k){
  return rfromk(*k);
}

real8 rgetgl1(const char *var) {
  return rgetgl1_(var, strlen(var));
}

integer8 ktfcopy(integer8 k) {
  return ktfcopy_(&k);
    }

integer4 itfsyserr(integer4 level) {
  return itfsyserr_(&level);
}

integer4 itfmessage(integer4 level, const char *msg, const char *arg) {
  return itfmessage_(&level, msg, arg, strlen(msg), strlen(arg));
}

integer8 ktfsymbolz(const char *symbol) {
  integer4 length = strlen(symbol);

  return ktfsymbolz_(symbol, &length, length);
}

integer8 ktfsymbolf(const char *symbol, logical4 l) {
  integer4 length = strlen(symbol);

  return ktfsymbolf_(symbol, &length, &l, length);
}

void tfree(integer8 ia) {
  __tfmem_MOD_tfree(&ia);
}

integer8 ktsalocb(integer4 mode, const char *str) {
  integer4 length = strlen(str);

  return ktsalocb_(&mode, str, &length, length);
}

integer8 ktsalocbl(integer4 mode, const char *str, integer4 length) {
  return ktsalocb_(&mode, str, &length, length);
}

/*
integer4 italoc(integer4 nd) {
  return italoc_(&nd);
  }*/

integer8 ktavaloc(integer4 mode, integer4 nd) {
  return ktavaloc_(&mode, &nd);
}

integer8 ktadaloc(integer4 mode, integer4 nd) {
  return ktadaloc_(&mode, &nd);
}

integer8 ktfmaloc(integer8 k, integer4 *n, integer4 *m,
		  logical4 vec, logical4 trans, integer4 *irtc) {
  return ktfmaloc_(&k, n, m, &vec, &trans, irtc);
}

integer8 ktfmakelist(integer4 isp1) {
  return ktfmakelist_(&isp1);
}

void tfsetlist(integer8 kx, integer8 iax, integer4 i) {
  tfsetlist_(&kx, &iax, &i);
}

void tfmakerulestk(integer8 ias, integer8 kx) {
  /*  fprintf(stderr,"mkrs_c: %d %f\n",kx); */
  tfmakerulestk_(&ias,&kx);
}

void tfevals(const char *buf,
	     integer8 *kx, integer4 *irtc) {
  integer4 length = strlen(buf);
  tfevalb_(buf, length, kx, irtc);
}

void tfevalc(const char *buf) {
  integer4 length = strlen(buf);

  tfevalc_(buf, length);
}

void tfdeval(integer4 isp1, integer8 iad,
	     integer8 *kx,
	     integer4 next, 
	     logical4 def, integer4 *ev, integer4 *irtc) {
  tfdeval_(&isp1, &iad, kx, &next, &def, ev, irtc);
}

void tflocal(integer8 kx) {
  tflocal_(&kx);
}

void tflocal1(integer8 kx) {
  tflocal1_(&kx);
}

/*
int tfinitstk(tfstk_t *stack, integer4 size) {
  integer4 stk_size, stk_offset;

  stk_size = size > 1024 ? size : 1024;
  stk_offset = italoc(stk_size * 2);
  if(stk_offset < 0)
    return -1;

  stack->__isporg      = stk_offset - 1;
  stack->__ivstkoffset = stk_size;
  stack->__isp         = stack->__isporg;
  stack->__mstk        = stack->__isporg + stack->__ivstkoffset;
  stack->__ipurefp     = 0;
  stack->__napuref     = 0;

  return 0;
  }*/

static bool itfgetoption1(integer8 ia, integer8 kr,
			  integer8 *kx) {
  integer8 iar;
  int i;

  iar = ktamask & kr;
  if(klist(iar) == ktfoper + mtflist) {
    for(i = 1; i <= ilist(2, iar - 1); i++) {
      if(itfgetoption1(ia, klist(iar + i), kx))
	return true;
    }
  } else if((ktfmask & klist(iar + 1)) == ktfsymbol
	    && tfsamesymbolqk_(&ia, &klist(iar + 1))) {
    *kx = klist(iar + 2);
    return true;
  }

  return false;
}

integer4 itfgetoptionstk(integer4 isp1, const char **optv) {
  integer8 *iopt;
  integer4 isp0, ispopt;
  int optc, i, j;

  optc = 0;
  while(optv[optc] != NULL) optc += 1;
  if(optc < 1) return -1;

  iopt = malloc(sizeof(integer8) * optc);
  if(iopt == NULL) return -1;

  for(i = 0; i < optc; i++)
    iopt[i] = ktfsymbolz(optv[i]);

  isp0 = isp;
  for(ispopt = isp0; isp1 <= ispopt; ispopt--)
    if(!tfruleqk_(&ktastk(ispopt)))
	break;
  ispopt += 1;

  isp = isp0 + optc;
  for(i = 1; i <= optc; i++)
    ktastk(isp0 + i) = ktfref;

  if(ispopt <= isp0)
     for(i = 1; i <= optc; i++)
       for(j = ispopt; j <= isp0; j++)
	 if(itfgetoption1(iopt[i - 1], ktastk(j), & ktastk(isp0 +i)))
	   break;

  free(iopt);
  return ispopt;
}

/* Debug routines */
#ifdef	DEBUG_DUMPER
void dumpstk_(integer4 *isp1) {
  int i, ia;

  printf(" Index    jstk(1)  jstk(2)    istk(1)    istk(2)      vstk()\n");
  for(i = 0, ia = *isp1; ia <= isp; i ++, ia++)
    printf("stk[%3d]:  %6d   %6d  0x%08x %10d    %e \n",
	   i, jtastk(1, ia), jtastk(2, ia), itastk(1, ia), itastk(2, ia),
	   vstk(ivstkoffset + ia));
  printf("\n");
}

void dumptok_(const_character str, integer4 *istart, integer4 *istop,
	      integer4 *iti, integer4 *iai, real8 *vi,
	      ftnlen slen) {
  int i;

  printf(" tfetok -> (itt=%d ia=%d v=%e): cur=\"",
	 *iti, *iai, *vi);
  for(i = 0; i < *istop - *istart; i++)
    switch(str[i]) {
    case '\\':
      printf("\\\\");
      break;

    case '"':
      printf("\\\"");
      break;

    case '\t':
      printf("\\t");
      break;

    case '\r':
      printf("\\r");
      break;

    case '\n':
      printf("\\n");
      break;

    case '\0':
      printf("\\0");
      break;

    default:
      printf("%c", str[i]);
      break;
    }

  printf("\" next=\"");

  for(i = *istop - *istart; i < slen; i++)
    switch(str[i]) {
    case '\\':
      printf("\\\\");
      break;

    case '"':
      printf("\\\"");
      break;

    case '\t':
      printf("\\t");
      break;

    case '\r':
      printf("\\r");
      break;

    case '\n':
      printf("\\n");
      break;

    case '\0':
      printf("\\0");
      break;

    default:
      printf("%c", str[i]);
      break;
    }

  printf("\"\n");
}
#endif	/* DEBUG_DUMPER */

/* End of File */
