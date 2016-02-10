#include <sim/sad_api.h>
#include <feature.h>

#include <stdlib.h>

/* SADScript function registration of RenumberElement stuff */
/* tfrenumber_element.f */
#define RenumberElement		tfrenumber_element_
extern DECLARE_SAD_FUNC8(tfrenumber_element_);

/* Module initializer */
#define REG	dlfunaloc
#define REG8	dlfunaloc8

#ifdef WITH_EXTENSION_MODULE
int dldeffun_(integer4 *isp1,
	      integer4 *itx, integer4 *iax, real8 *vx,
	      integer4 *irtc) {
#else
int sadDefFunc_RenumberElement(void) {
#endif /* WITH_EXTENSION_MODULE */
  if(!feature_provide("RenumberElement/API",
		      FEATURE_VERSION(1, 0), FEATURE_VERSION(1, 0)))
    return -1;

  REG8("RenumberElement",	RenumberElement,	2, NULL, NULL, 0);

  return 0;
}
/* End of File */
