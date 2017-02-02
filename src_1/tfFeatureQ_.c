#include <sim/sad_api.h>
#include <sim/TFCODE.h>
#include <sim/TFCBK.h>
#include <sim/TFSTK.h>

#include <feature.h>

static void FeatureScanCB(const char *entry) {
  isp += 1;
  ktastk(isp) = ktfstring + ktsalocb(-1, entry);
}

/* SADScript function definition */
static int FeatureQ(integer8 *isp1,
		    integer8 *kx,
		    integer4 *irtc) {
  integer8 ia;
  integer8 isp0;
  real8 vx;
  uint32_t version = 0UL;

  if(isp == *isp1 + 1 && ktastk(*isp1 + 1) == ktfoper + mtfnull){

    isp0 = isp;
    feature_scan(FeatureScanCB);

    *kx = ktflist + ktfmakelist(isp0);
    isp = isp0;
    *irtc = 0;
    return 0;
  }

  switch(isp - *isp1) {
  case 2:
    if(ktfnonrealq(ktastk(*isp1 + 2))
       || rtastk(*isp1 + 2) < 0) {
      *irtc = itfmessage(9, "General::wrongtype",
			 "\"Required version number for #2\"");
      return -1;
    }
    version = rtastk(*isp1 + 2);
  case 1:
    if(ktfnonstringq(ktastk( *isp1 + 1))) {
      *irtc = itfmessage(9, "General::wrongtype",
			 "\"Feature name string for #1\"");
      return -1;
    }
    ia = ktfaddr(ktastk(*isp1 + 1));
#if SAD_REQUIRE_STRING_TERMINATION
    jlist(ilist(1, ia) + 1, ia + 1) = '\0';
#endif
    break;
  default:
    *irtc = itfmessage(9, "General::narg", "\"1 or 2\"");
    return -1;
    break;
  }

  if(version == 0UL)
    vx = feature_version(&jlist(1, ia + 1));
  else
    vx = SAD_BOOLEAN(feature_require(&jlist(1, ia + 1), version));
  *kx = kfromr(vx);
  *irtc = 0;
  return 0;
}

/* SADScript function registration */
#define REG	dlfunaloc
#define REG8	dlfunaloc8
int sadDefFunc_FeatureQ(void) {
  REG8("FeatureQ",		FeatureQ,	1, NULL, NULL, 0);

  return 0;
}

/* End of File */
