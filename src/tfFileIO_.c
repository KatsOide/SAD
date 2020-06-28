#include <feature.h>

#include <sim/sad_api.h>
#include <sim/TFCODE.h>
#include <sim/TFCBK.h>
#include <sim/TFSTK.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <sys/param.h>
#include <errno.h>
#include <unistd.h>
#include <string.h>
#include <stdlib.h>

/* SADScript function definition of FileIO stuff */
static int ExpandTilde(integer4 *isp1,
		       integer8 *kx,
		       integer4 *irtc) {
  integer8 ia;
  char *expanded;

  if(isp != *isp1 + 1) {
    *irtc = itfmessage(9, "General::narg", "\"1\"");
    return -1;
  }

  if((ktfmask & ktastk(*isp1 + 1)) != ktfstring) {
    *irtc = itfmessage(9, "General::wrongtype", "\"String\"");
    return -1;
  }

  ia = (ktamask & ktastk(*isp1 + 1));
#if SAD_REQUIRE_STRING_TERMINATION
  jlist(ilist(1, ia) + 1, ia + 1) = '\0';
#endif

  expanded = expand_tilde(&jlist(1, ia + 1));
  if(expanded == NULL) {
    switch(errno) {
    case ENOMEM:
      *irtc = itfsyserr(9);
      break;

    default:
      *irtc = itfmessage(9, "System::error",
			 "\"No such home directory\"");
    }
    return 1;
  }

  *kx = ktfstring + ktsalocb(-1, expanded); free(expanded);
  *irtc = 0;
  return 0;
}

static int FromFileDate(integer4 *isp1,
			integer8 *kx,
			integer4 *irtc) {
  integer8 ia;
  real8 vx;
  char *expanded;
  struct stat buf;

  if(isp != *isp1 + 1) {
    *irtc = itfmessage(9, "General::narg", "\"1\"");
    return -1;
  }
  if((ktfmask & ktastk( *isp1 + 1)) != ktfstring) {
    *irtc = itfmessage(9, "General::wrongtype", "\"String\"");
    return -1;
  }

  ia = (ktamask & ktastk(*isp1 + 1));
#if SAD_REQUIRE_STRING_TERMINATION
  jlist(ilist(1, ia) + 1, ia + 1) = '\0';
#endif

  expanded = expand_tilde(&jlist(1, ia + 1));
  if(expanded == NULL) {
*kx = ktfoper + mtfnull;;
    *irtc = 0;
#ifdef SAD_THROW_EXCEPTION
    switch(errno) {
    case ENOMEM:
      *irtc = itfsyserr(9);
      break;

    default:
      *irtc = itfmessage(9, "System::error",
			 "\"No such home directory\"");
    }
#ifdef SAD_NOBREAK_EXCEPTION
    *irtc = 0;
#endif
#endif
    return 1;
  }

  if(stat(expanded, &buf) != 0) {
    free(expanded);
    *kx = ktfoper + mtfnull;
    *irtc = 0;
#ifdef SAD_THROW_EXCEPTION
    *irtc = itfsyserr(9);
#ifdef SAD_NOBREAK_EXCEPTION
    *irtc = 0;
#endif
#endif
    return 1;
  }
  free(expanded);

  vx  = buf.st_mtime + SAD_EPOCH_OFFSET;
  *kx = kfromr(vx);
  *irtc = 0;
  return 0;
}

static int ReadLink(integer4 *isp1,
		    integer8 *kx,
		    integer4 *irtc) {
  integer8 ia;
  integer4 nc;
  char *expanded, buf[MAXPATHLEN];

  if(isp != *isp1 + 1) {
    *irtc = itfmessage(9, "General::narg", "\"1\"");
    return -1;
  }
  if((ktfmask & ktastk( *isp1 + 1)) != ktfstring) {
    *irtc = itfmessage(9, "General::wrongtype", "\"String\"");
    return -1;
  }

  ia = (ktamask & ktastk(*isp1 + 1));
#if SAD_REQUIRE_STRING_TERMINATION
  jlist(ilist(1, ia) + 1, ia + 1) = '\0';
#endif

  expanded = expand_tilde(&jlist(1, ia + 1));
  if(expanded == NULL) {
    *kx = ktfoper + mtfnull;;
    *irtc = 0;
#ifdef SAD_THROW_EXCEPTION
    switch(errno) {
    case ENOMEM:
      *irtc = itfsyserr(9);
      break;

    default:
      *irtc = itfmessage(9, "System::error",
			 "\"No such home directory\"");
    }
#ifdef SAD_NOBREAK_EXCEPTION
    *irtc = 0;
#endif
#endif
    return 1;
  }

  nc = readlink(expanded, buf, sizeof(buf)); free(expanded);
  if(nc < 0) {
    *kx = ktfoper + mtfnull;;
    *irtc = 0;
#ifdef SAD_THROW_EXCEPTION
    *irtc = itfsyserr(9);
#ifdef SAD_NOBREAK_EXCEPTION
    *irtc = 0;
#endif
#endif
    return 1;
  }

  *kx = ktfstring + ktsalocbl(-1, buf, nc);
  *irtc = 0;
  return 0;
}

static int RealPath(integer4 *isp1,
		    integer8 *kx,
		    integer4 *irtc) {
  integer8 ia;
  char *expanded, *resolved, resolved_buf[MAXPATHLEN];

  if(isp != *isp1 + 1) {
    *irtc = itfmessage(9, "General::narg", "\"1\"");
    return -1;
  }
  if((ktfmask & ktastk(*isp1 + 1)) != ktfstring) {
    *irtc = itfmessage(9, "General::wrongtype", "\"String\"");
    return -1;
  }

  ia = (ktamask & ktastk(*isp1 + 1));
#if SAD_REQUIRE_STRING_TERMINATION
  jlist(ilist(1, ia) + 1, ia + 1) = '\0';
#endif

  expanded = expand_tilde(&jlist(1, ia + 1));
  if(expanded == NULL) {
    *kx = kxfailed;
    *irtc = 0;
#ifdef SAD_THROW_EXCEPTION
    switch(errno) {
    case ENOMEM:
      *irtc = itfsyserr(9);
      break;

    default:
      *irtc = itfmessage(9, "System::error",
			 "\"No such home directory\"");
    }
#ifdef SAD_NOBREAK_EXCEPTION
    *irtc = 0;
#endif
#endif
    return 1;
  }

  resolved = realpath(expanded, resolved_buf); free(expanded);
  if(resolved == NULL) {
    *kx = kxfailed;
    *irtc = 0;
#ifdef SAD_THROW_EXCEPTION
    *irtc = itfsyserr(9);
#ifdef SAD_NOBREAK_EXCEPTION
    *irtc = 0;
#endif
#endif
    return 1;
  }

  *kx = ktfstring + ktsalocb(-1, resolved);
  *irtc = 0;
  return 0;
}

static int MkSecureTemp(integer4 *isp1,
			integer8 *kx,
			integer4 *irtc) {
  integer8 ka1, ka2 = 0;
  int fd, nc, slen;
  char *template;

  if(!(*isp1 < isp && isp < *isp1 + 3)) {
    *irtc = itfmessage(9, "General::narg", "\"1 or 2\"");
    return -1;
  }
  if((ktfmask & ktastk(*isp1 + 1)) != ktfstring) {
    *irtc = itfmessage(9, "General::wrongtype", "\"String\"");
    return -1;
  }

  ka1 = ktamask & ktastk(*isp1 + 1);
#if SAD_REQUIRE_STRING_TERMINATION
  jlist(ilist(1, ka1) + 1, ka1 + 1) = '\0';
#endif
  nc = itastk(1, ka1);
  slen = 0;
  if(nc < 1) {
    *kx = kxfailed;
    *irtc = 0;
    return -1;
  }

  if(isp == *isp1 + 2) {
    if((ktfmask & ktastk(*isp1 + 2)) != ktfstring) {
      *irtc = itfmessage(9, "General::wrongtype", "\"String\"");
      return -1;
    }
    ka2 = (ktamask & ktastk(*isp1 + 2));
#if SAD_REQUIRE_STRING_TERMINATION
    jlist(ilist(1, ka2) + 1, ka2 + 1) = '\0';
#endif
    slen = ktastk(ka2);
  }

  template = malloc(nc + slen + 1);
  if(template == NULL) {
    *kx = kxfailed;
    *irtc = 0;
    return -1;
  }

  memcpy(template, &jlist(1, ka1 + 1), nc);
#ifdef HAVE_MKSTEMPS
  if(slen > 0) memcpy(template + nc, &jlist(1, ka2 + 1), slen);
#else
  slen = 0; /* Remove suffix */
#endif
  template[nc + slen] = '\0';

#ifdef HAVE_MKSTEMPS
  if(slen > 0) {
    fd = mkstemps(template, slen);
  } else
#endif
    fd = mkstemp(template);

  if(fd != -1) {
    close(fd);
    *kx = ktfstring + ktsalocb(-1, template);
  } else {
    *kx = kxfailed;
  }

  free(template);
  *irtc = 0;
  return 0;
}

static int Pipe(integer4 *isp1, integer8 *kx,
		integer4 *irtc) {
  integer4 isp0, luns[2];
  int fildes[2],moder=1,modew=2;

  if(isp != *isp1 + 1
     || ktastk(isp) != ktfoper + mtfnull) {
    *irtc = itfmessage(9, "General::narg","\"0\"");
    return -1;
  }

  if(pipe(fildes) == -1) {
    *irtc = itfsyserr(9);
    return 1;
  }

  luns[0] = itopenbuf_(&moder,irtc); if(*irtc != 0) return 1;
  luns[1] = itopenbuf_(&modew,irtc); if(*irtc != 0) return 1;

  if(dup2(fildes[0], getfd_(&luns[0])) == -1 ||
     dup2(fildes[1], getfd_(&luns[1])) == -1) {
    close(fildes[0]); close(fildes[1]);
    *irtc = itfsyserr(9);
    return 1;
  }

  close(fildes[0]); close(fildes[1]);

  isp0 = isp;
  isp += 1;
  rtastk(isp) = luns[0];
  isp += 1;
  rtastk(isp) = luns[1];

  *kx = ktflist + ktfmakelist(isp0); 
  isp = isp0;
  *irtc = 0;
  return 0;
}

static int SetLUN2FD(integer4 *isp1, integer8 *kx,
		     integer4 *irtc) {
  integer4 lun;
  int fd;

  if(isp != *isp1 + 2) {
    *irtc = itfmessage(9, "General::narg", "\"2\"");
    return -1;
  }
  if((ktrmask & ktastk(*isp1 + 1)) == ktfnr
     || rtastk(*isp1 + 1) < 0) {
    *irtc = itfmessage(9, "General::wrongtype",
		       "\"Logical unit number for #1\"");
    return -1;
  }
  if((ktrmask & ktastk(*isp1 + 2)) == ktfnr
     || rtastk(*isp1 + 2) < 0) {
    *irtc = itfmessage(9, "General::wrongtype",
		       "\"File descriptor for #2\"");
    return -1;
  }

  lun = rtastk(*isp1 + 1);
  fd  = rtastk(*isp1 + 2);

  /* fprintf(stderr,"SetLUN %d %d %d\n",lun,fd,getfd_(&lun)); */

  if(dup2(getfd_(&lun),fd) == -1) {
    *irtc = itfsyserr(9);
    return 1;
  }

  *kx = ktfoper + mtfnull;
  *irtc = 0;
  return 0;
}

/* SADScript function registration of FileIO stuff */
#define REG	dlfunaloc
#define REG8	dlfunaloc8
#ifdef WITH_EXTENSION_MODULE
int dldeffun_(void) {
#else
int sadDefFunc_FileIO(void) {
#endif /* WITH_EXTENSION_MODULE */
  if(!feature_provide("FileIO/API",
                      FEATURE_VERSION(1, 0), FEATURE_VERSION(1, 0)))
    return -1;

  REG8("ExpandTilde",	ExpandTilde,		1, NULL, NULL, 0);
  REG8("FromFileDate",	FromFileDate,		1, NULL, NULL, 0);
  REG8("ReadLink",	ReadLink,		1, NULL, NULL, 0);
  REG8("RealPath",	RealPath,		1, NULL, NULL, 0);
  REG8("MkSecureTemp$",	MkSecureTemp,		2, NULL, NULL, 0);
  REG8("Pipe",		Pipe,			0, NULL, NULL, 0);
  REG8("SetLUN2FD",	SetLUN2FD,		2, NULL, NULL, 0);

  return 0;
}

/* End of File */
