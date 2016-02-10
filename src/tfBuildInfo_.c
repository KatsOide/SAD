#include <buildinfo.h>
#include <feature.h>

#include <sim/sad_api.h>
#include <sim/TFCODE.h>
#include <sim/TFCBK.h>
#include <sim/TFSTK.h>

#include <stdlib.h>
#include <string.h>

/* SADScript function definition */
static void BuildInfoScanCB(const buildinfo_t* entry) {
  isp += 1;
  ktastk(isp) = ktfstring + ktsalocb(-1, buildinfo_extract_keyword(entry));
}

static int BuildInfo(integer4 *isp1,
		     integer8 *kx, 
		     integer4 *irtc) {
  integer8 ka;
  integer4 isp0;
  real8 vx;

  if(isp != *isp1 + 1) {
    *irtc = itfmessage(9, "General::narg", "\"0 or 1\"");
    return -1;
  }

  if(ktastk(isp) == ktfoper + mtfnull) {
    isp0 = isp;
    buildinfo_scan(BuildInfoScanCB);
    *kx = ktflist + ktfmakelist(isp0);
    isp = isp0;
  } else if((ktfmask & ktastk(*isp1 + 1)) == ktfstring) {
    ka = (ktamask & ktastk(*isp1 + 1));
#if SAD_REQUIRE_STRING_TERMINATION
    jlist(ilist(1, ka) + 1, ka + 1) = '\0';
#endif

    const buildinfo_t* entry = buildinfo_search(&jlist(1, ka + 1));
    if(buildinfo_integer_entryQ(entry)) {
      vx  = buildinfo_extract_integer(entry);
      *kx = (integer8) vx;
    } else if(buildinfo_string_entryQ(entry)) {
      *kx = ktfstring + ktsalocb(-1, buildinfo_extract_string(entry));
    } else {
      *kx = ktfoper + mtfnull;
    }
  } else {
    *irtc = itfmessage(9, "General::wrongtype",
                       "\"Null or Build Information keyword-string\"");
    return -1;
  }

  *irtc = 0;
  return 0;
}

/* SADScript function registration */
#define REG	dlfunaloc
#define REG8	dlfunaloc8
#ifdef WITH_EXTENSION_MODULE
int dldeffun_(void) {
#else
int sadDefFunc_BuildInfo(void) {
#endif /* WITH_EXTENSION_MODULE */
  if(!(buildinfo_init() &&
       feature_provide("BuildInfo/API",
		       FEATURE_VERSION(1, 0), FEATURE_VERSION(1, 0)))) {
    return -1;
  }

  REG8("BuildInfo",		BuildInfo,	1, NULL, NULL, 0);

  return 0;
}

/* End of File */
