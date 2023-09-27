#include <sim/sad_f2c.h>
#include <sim/sad_api.h>
#include <sim/TFCBK.h>

#include <float.h>
#include <math.h>

#include <sys/types.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <regex.h>
#include <fcntl.h>
#include <signal.h>

/* for regerror(3) buffer */
#define REGEXP_ERRBUF_SIZE 80

/* Support NAN/INFINITY by using DBL_QNAN/INFINITY for non-C99 compiler */
#if !defined(NAN)	&& defined(DBL_QNAN)
#define NAN DBL_QNAN
#endif
#if !defined(INFINITY)	&& defined(DBL_INFINITY)
#define INFINITY DBL_INFINITY
#endif

/* Machine independent Infinity/NaN Real*8 constant initialization */
void tfinfinit_(void) {
  /* NAN & INFINITY float constant are supplied by math.h since C99 */
  dinfinity   = INFINITY; /* Convert float Infinity  constant to Real*8 */
  /* dnotanumber = 0.0/0.0;	  /* Convert float quite NaN constant to Real*8 */
  dnotanumber = NAN;	  /* Convert float quite NaN constant to Real*8 */
}

/* stub routines for Fortran function call 
integer4 isnan_(real8 *vx) {
  return isnan(*vx);
  }*/

real8 second_(void) {
  struct rusage u_ru;

  if(getrusage(RUSAGE_SELF, &u_ru))
    return 0.0;

  return u_ru.ru_utime.tv_sec + u_ru.ru_utime.tv_usec * 1e-6;
}

integer4 tpause_(integer4 *microsec) {
  return (*microsec > 0) ? usleep(*microsec) : 0;
}

integer4 tkill_(integer4 *pr) {
  return kill((pid_t) pr, SIGTERM);
}

integer4 unixclose_(integer4 *fd) {
  return close(*fd);
}

/* regexp(string, regular expression)
 * Result:         1: Match
 *                 0: No Match
 *         otherwise: something wrong
 * Note: string, regular expression arguments
 *        MUST be C style `\0' terminated string.
 */
integer4 regexp_(const char *string, const char *pat) {
  static char ebuf[REGEXP_ERRBUF_SIZE];
  static char *pattern = NULL;
  static int pSize = 0;
  static regex_t preg;

  regmatch_t pmatch;
  char *buf;
  int error, plen;

  /* check input pattern length */
  plen = strlen(pat);

  if(pattern && strlen(pattern) == plen) {
    if(strncmp(pattern, pat, plen) == 0) goto reg_exec;
    regfree(&preg);
  } else if(plen + 1 > pSize) {
    buf = realloc(pattern, plen + 1);
    if(!buf) {
      error = REG_ESPACE;
      goto reg_error;
    }
    if(pattern) regfree(&preg);
    pattern = buf; pSize = plen + 1;
  } else regfree(&preg);

  /* Copy pattern string and compile it */
  memcpy(pattern, pat, plen); pattern[plen] = '\0';

  error = regcomp(&preg, pattern, REG_EXTENDED | REG_NOSUB);
  if(error) {
    pattern[0] = '\0'; /* Clear pattern buffer if regcomp failed */
    goto reg_error;
  }

 reg_exec:
  error = regexec(&preg, string, 0, &pmatch, 0);

 reg_error:
  if(!error) return 1;
  if(error == REG_NOMATCH) return 0;
  if(error == REG_ESPACE) {
    fprintf(stderr, "regexp_: can't allocate buffer\n");
  } else {
    regerror(error, &preg, ebuf, sizeof(ebuf));
    fprintf(stderr, "regexp_: %s\n", ebuf);
  }

  return -2;
}

/* Expand `~' prefix and return number of un-written-back characters */
integer4 cfexptilde_(character path, ftnlen plen)
{
  int len;
  char *buf, *expanded;

  if(plen < 2 || path[0] != '~') {
    return 0;
  }

  len = len_trim(path, plen);

  buf = malloc(len + 1);
  if(buf == NULL) {
    return -1;
  }
  memcpy(buf, path, len); buf[len] = '\0';

  expanded = expand_tilde(buf); free(buf);
  if(expanded == NULL) {
    return -1;
  }

  len = strlen(expanded);
  memcpy(path, expanded, len < plen ? len : plen); free(expanded);
  while(len < plen) {
    path[len] = ' ';
    len += 1;
  }
  return (plen > len) ? (plen - len) : 0;
}

/* tempnam(3) wrapper for Fortran [unsafe] */
integer4 tftmpnam_(character s, ftnlen slen) {
  /*  char buf[] = "/tmp/fileXXXXXX"; */
  size_t blen;

  int fd = mkstemp(s);
  fcntl(fd, F_SETFD, FD_CLOEXEC);

  if(!fd) {
    return 0;
  }

  close(fd);

  blen = strlen(s);
  if(slen < blen) {
    return 0;
  }

  return blen;
}

/* End of File */
