#include <feature.h>

#include <sim/sad_api.h>
#include <sim/TFCODE.h>
#include <sim/TFCBK.h>
#include <sim/TFSTK.h>

#include <sim/dynl.h>

#include <stdlib.h>
#include <stdio.h>

/* SADScript function for dynamic loader API */
static int DynamicLink(integer8 *isp1,
		       integer8 *kx,
		       integer4 *irtc) {
  integer8 ia;
  real8 vx;
  int i, id, mode, mode0;

  if(isp < *isp1 + 1) {
    *irtc = itfmessage(9, "General::narg", "\"1 or more\"");
    return -1;
  }

  if(ktfnonstringq(ktastk(*isp1 + 1))) {
    *irtc = itfmessage(9, "General::wrongtype", "\"Character-string\"");
    return -1;
  }

  ia = ktfaddr(ktastk(*isp1 + 1));
#if SAD_REQUIRE_STRING_TERMINATION
  jlist(ilist(1, ia) + 1, ia + 1) = '\0';
#endif

  mode = DYNL_LAZY; mode0 = 0;
  for(i = 2; *isp1 + i <= isp; i++) {
    if(ktfrealq(ktastk(*isp1 + i))) {
      *irtc = itfmessage(9, "General::wrongtype",
                         "\"DYNL flag for ##2\"");
    }
    id = rtastk(*isp1 + i);
    mode0 |= id;
    if(id & DYNL_MODEMASK) mode = id & DYNL_MODEMASK;
  }

  if(mode0 & DYNL_GLOBAL) mode |= DYNL_GLOBAL;

  id = dl_link(&jlist(1, ia + 1), mode);
  if(id > 0) {
    vx = id;
    *kx = kfromr(vx);
  } else {
    *kx = kxfailed;
  }
  *irtc = 0;
  return 0;
}

static int DynamicUnlink(integer8 *isp1,
			 integer8 *kx,
			 integer4 *irtc) {
  real8 vx;
  int id, status;

  if(isp != *isp1 + 1) {
    *irtc = itfmessage(9, "General::narg", "\"1\"");
    return -1;
  }

  if(ktfrealq(ktastk(isp))) {
    *irtc = itfmessage(9, "General::wrongtype",
		       "\"Dynamic object handle\"");
    return -1;
  }

  id = rtastk(isp);
  status = dl_unlink(id);

  vx   = SAD_BOOLEAN(status == 0);
  *kx = kfromr(vx);
  *irtc = 0;
  return 0;
}

#define SAD_INVALID_IRTC	INT32_MIN
static int DynamicCall(integer8 *isp1,
		       integer8 *kx,
		       integer4 *irtc) {
  integer8 ia, kx0;
  integer8 isp0, isp2;
  integer8 isp10;
  char *msg;
  int id;
  sad_func_t_8 func;

  if(isp < *isp1 + 2) {
    *irtc = itfmessage(9, "General::narg", "\"2 or more\"");
    return -1;
  }

  if(ktfnonrealq(ktastk(*isp1 + 1))) {
    *irtc = itfmessage(9, "General::wrongtype",
		       "\"Dynamic object handle for #1\"");
    return -1;
  }

  if(ktfnonstringq(ktastk(*isp1 + 2))) {
    *irtc = itfmessage(9, "General::wrongtype",
		       "\"Symbol name-string for #2\"");
    return -1;
  }

  id = rtastk(*isp1 + 1);
  ia = ktfaddr(ktastk(*isp1 + 2));
#if SAD_REQUIRE_STRING_TERMINATION
  jlist(ilist(1, ia) + 1, ia + 1) = '\0';
#endif

  func = (sad_func_t_8) dynl_function(dl_ref_handle(id), &jlist(1, ia + 1));
  if(func == NULL) {
    if(asprintf(&msg, "\"%s\"", &jlist(1, ia + 1)) >= 0) {
      *irtc = itfmessage(9, "DynamicLoader::nosymbol", msg);
      free(msg);	   
    } else {
      *irtc = itfmessage(9, "Memory::alloc", " ");
    }
    return -1;
  }

  isp0 = isp;		/* Backup previous stack pointer */
  isp2 = *isp1 + 2;	/* Shift stack frame(skip 2 argument) */
  if(isp == isp2) {	/* Append ``Null'' for zero argement case */
    isp += 1;
    ktastk(isp) = ktfoper + mtfnull;
  }

  /* Backup copy of arguments */
  isp10 = isp2;
  kx0 = *kx;

  /* Set invalid irtc for SADScript function call API */
  *irtc = SAD_INVALID_IRTC;

  func(&isp2, kx, irtc);

  /* Check argument changes.. */
  if(isp10 == isp2 && kx0 == *kx
     && SAD_INVALID_IRTC == *irtc) { /* func might be (*)(void) */
    /* In SADScript funtion call API
     * irtc	 0: success
     *		>0: iftmessage
     *		-1: Error
     *		-2: [used]
     *		-4: [un-used]
     *		-3: [used]
     *		-5: [used]
     *		-6: [used]
     */

    *kx = ktfoper + mtfnull; /* Return `Null' for (*)(void) results */
  }

  isp = isp0;		/* Restore stack pointer */
  return 0;
}

#define REG	dlfunaloc
#define REG8	dlfunaloc8
int sadDefFunc_DynamicLink(void) {
  if(!(dynl_provide() && dynl_init()))
    return -1;

  if(!feature_provide("DynamicLink/API",
		      FEATURE_VERSION(1, 0), FEATURE_VERSION(1, 0)))
    return -1;

  REG8("DynamicLink",	DynamicLink,	1, NULL, NULL, 0);
  REG8("DynamicUnlink",	DynamicUnlink,	1, NULL, NULL, 0);
  REG8("DynamicCall",	DynamicCall,	2, NULL, NULL, 0);

  dynl_error_report_level += 1;

  return 0;
}

/* End of File */
