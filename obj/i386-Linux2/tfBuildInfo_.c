
#include <sim/sad_api.h>
#include <sim/TFCODE.h>
#include <sim/TFCBK.h>
#include <sim/TFSTK.h>

#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sysexits.h>

typedef enum {
  _NOP,
  _INT,
  _STR,
} info_type_t;

typedef struct {
  const char *tag;
  info_type_t type;
  int32_t number;
  const char *string;
} build_info_t;

static const char sad_pkg_root[] = "/Users/oide/SAD/oldsad//share/Packages";

static build_info_t build_info_table[] = {
  {"Repository:Owner",		_STR, 0, "sadist"},
  {"Repository:Type",		_STR, 0, "CVS"},
  {"Repository:Location",	_STR, 0, "cvs://acsad0.kek.jp/SAD/cvsroot/oldsad"},
  {"Source:TreeName",		_STR, 0, "MAIN trunk"},
#if 0
  {"Source:Revision",		_INT,	 %%SourceRevision%%, NULL},
  {"Source:Date",		_STR, 0, %%SourceDate%%},
#endif
  {"Config:FC",			_STR, 0, "gfortran"},
  {"Config:CC",			_STR, 0, "gcc-4"},
  {"Config:CXX",		_STR, 0, "c++"},
  {"Config:SYS_FOPT",		_STR, 0, "-Wall"},
  {"Config:SYS_COPT",		_STR, 0, "-Wall -std=gnu99"},
  {"Config:SYS_CXXOPT",		_STR, 0, "-Wall -std=c++98 -pedantic-errors -mpreferred-stack-boundary=4"},
  {"Config:FOPT",		_STR, 0, "-O3 -fno-second-underscore -fdollar-ok -fargument-alias -mpreferred-stack-boundary=4 -mfancy-math-387 -frecursive -fbackslash -std=legacy -fall-intrinsics"},
  {"Config:COPT",		_STR, 0, "-g -O1"},
  {"Config:CXXOPT",		_STR, 0, "-g -O1"},
  {"Target:OS_TYPE",		_STR, 0, "Linux"},
  {"Target:OS_NAME",		_STR, 0, "Linux"},
  {"Target:OS_MAJOR",		_INT,	 2, NULL},
  {"Target:OS_MINOR",		_INT,	 6, NULL},
  {"Target:OS_RELEASE",		_STR, 0, "2.6.32-27-generic"},
  {"Target:CPU_ARCH",		_STR, 0, "i386"},
  {"Target:MACH_ARCH",		_STR, 0, "i386-Linux2"},
  {"Target:SAD_ROOT",		_STR, 0, "/Users/oide/SAD/oldsad/"},
  {"Target:SAD_ARCH_ROOT",	_STR, 0, "/Users/oide/SAD/oldsad//arch"},
  {"Target:SAD_SHARE_ROOT",	_STR, 0, "/Users/oide/SAD/oldsad//share"},
  {"Target:SAD_PKG_ROOT",	_STR, 0, "/Users/oide/SAD/oldsad//share/Packages"},
  {"Target:SAD_MOD_ROOT",	_STR, 0, "/Users/oide/SAD/oldsad//share/Extension"},
  {"Built:Location",		_STR, 0, "oide@ubuntu:/Users/oide/SAD/oldsad/obj/i386-Linux2"},
  {"Built:Date",		_STR, 0, "2010-12-28 17:08:30 -0800"},
  {NULL,			_NOP, 0, NULL}};

static bool check_build_info(void) {
  build_info_t *p;
  bool success = true;

  for(p = build_info_table; p->tag != NULL; p++) {
    if(strlen(p->tag) < 1) {
      success = false;
      fprintf(stderr, "build_info_table contains null tag string!\n");
    }
    switch(p->type) {
    case _INT:
      break;

    case _STR:
      if(p->string == NULL) {
	success = false;
	fprintf(stderr,
		"build_info_table does not contain string data at tag[%s]!\n",
		p->tag);
      }
      break;

    default:
      success = false;
      fprintf(stderr,
	      "build_info_table contains wrong data type[%d] at tag[%s]!\n",
	      p->type, p->tag);
      break;
    }
  }

  return success;
}

/* SADScript function definition */
static int BuildInfo(integer4 *isp1,
		     integer4 *itx, integer4 *iax, real8 *vx,
		     integer4 *irtc) {
  integer4 isp0, ia;
  build_info_t *p;

  if(isp != *isp1 + 1) {
    *irtc = itfmessage(9, "General::narg", "\"0 or 1\"");
    return -1;
  }

  if(itastk(1, isp) == ntfoper && itastk(2, isp) == mtfnull) {
    isp0 = isp;
    for(p = build_info_table; p->tag != NULL; p++) {
      isp += 1;
      itastk(1, isp) = ntfstring;
      itastk(2, isp) = itsalocb(-1, p->tag);
    }
    *itx = ntflist;
    *iax = itflist(isp0);
    isp = isp0;
  } else if(itastk(1, *isp1 + 1) == ntfstring) {
    ia = itastk(2, *isp1 + 1);
#if SAD_REQUIRE_STRING_TERMINATION
    jlist(ilist(1, ia) + 1, ia + 1) = '\0';
#endif
    for(p = build_info_table; p->tag != NULL; p++)
      if(strcmp(p->tag, &jlist(1, ia + 1)) == 0) break;
    if(p->tag != NULL) {
      switch(p->type) {
      case _INT:
	*itx = ntfreal;
	*vx  = p->number;
	break;

      case _STR:
	*itx = ntfstring;
	*iax = itsalocb(-1, p->string);
	break;

      default:
	*itx = ntfoper;
	*iax = mtfnull;
	break;
      }
    } else {
      *itx = ntfoper;
      *iax = mtfnull;
    }
  } else {
    *irtc = itfmessage(9, "General::wrongtype",
                       "\"Null or Build Information tag-string\"");
    return -1;
  }

  *irtc = 0;
  return 0;
}

/* SADScript function registration */
#define REG	dlfunaloc
int sadDefFunc_BuildInfo(void) {
  int id;

  if(!check_build_info()) exit(EX_CONFIG);

  id = REG("BuildInfo",		BuildInfo,	1, NULL, NULL, 0);

  return 0;
}

/* Fortran query function for Target:SAD_PKG_ROOT */
void sadpkgdir_(character s, ftnlen slen) {
  size_t len = strlen(sad_pkg_root);
  int i;

  if(len > slen) return;

  strncpy(s, sad_pkg_root, slen);
  for(i = len; i < slen; i++) s[i] = ' ';
}

/* End of File */
