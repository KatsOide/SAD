#include <sim/sad_api.h>
#include <sim/TFCODE.h>
#include <sim/TFSTK.h>

#include "netSemaphore.h"

#include <string.h>

/*
 * SADScript function prototype of NetSemaphore
 *
 * NetSemTake[ServerName, SemaphoreID, TimeOut]
 * Return value:	 0(success)
 *			-1(failed)
 *
 * NetSemGive[ServerName, SemaphoreID]
 * Return value:	 0(success)
 *			-1(failed)
 *
 * NetSemInfo[ServerName, SemaphoreID, InformationType]
 * Return value:	String(semaphore information)
 *			Null(failed)
 */

/* SADScript function definition of NetSemaphore stuff */
static int tfNetSemTake(integer4 *isp1,
			integer8 *kx,
			integer4 *irtc) {
  integer8 ia;
  real8 vx;
  int status, semaphore, timeout;

  if(isp != *isp1 + 3) {
    *irtc = itfmessage(9, "General::narg", "\"3\"");
    return -1;
  }

  if(ktfnonstringq(ktastk( *isp1 + 1))) {
    *irtc = itfmessage(9, "General::wrongtype", "\"Character-string for #1\"");
    return -1;
  }

  if(ktfnonrealq(ktastk( *isp1 + 2))) {
    *irtc = itfmessage(9, "General::wrongtype", "\"Real Number for #2\"");
    return -1;
  }

  if(ktfnonrealq(ktastk( isp))) {
    *irtc = itfmessage(9, "General::wrongtype", "\"Real Number for #3\"");
    return -1;
  }

  ia = ktfaddr(ktastk(*isp1 + 1));
#if SAD_REQUIRE_STRING_TERMINATION
  jlist(ilist(1, ia) + 1, ia + 1) = '\0';
#endif
  semaphore = rtastk(*isp1 + 2);
  timeout   = rtastk( isp);

  status = netSemTake(&jlist(1, ia + 1), semaphore, timeout);

  vx  = status;
  *kx = kfromr(vx);
  *irtc = 0;
  return 0;
}

static int tfNetSemGive(integer4 *isp1,
			integer8 *kx,
			integer4 *irtc) {
  integer8 ia;
  real8 vx;
  int status, semaphore;

  if(isp != *isp1 + 2) {
    *irtc = itfmessage(9, "General::narg", "\"2\"");
    return -1;
  }

  if(ktfnonstringq(ktastk( *isp1 + 1))) {
    *irtc = itfmessage(9, "General::wrongtype", "\"Character-string for #1\"");
    return -1;
  }

  if(ktfnonrealq(ktastk( isp))) {
    *irtc = itfmessage(9, "General::wrongtype", "\"Real Number for #2\"");
    return -1;
  }

  ia = ktfaddr(ktastk(*isp1 + 1));
#if SAD_REQUIRE_STRING_TERMINATION
  jlist(ilist(1, ia) + 1, ia + 1) = '\0';
#endif
  semaphore = rtastk( isp);

  status = netSemGive(&jlist(1, ia + 1), semaphore);

  vx  = status;
  *kx = kfromr(vx);
  *irtc = 0;
  return 0;
}

static int tfNetSemInfo(integer4 *isp1,
			integer8 *kx,
			integer4 *irtc) {
  integer8 ia;
  real8 vx;
  int status, semaphore, select;
  char buffer[MSG_LEN];

  if(isp != *isp1 + 3) {
    *irtc = itfmessage(9, "General::narg", "\"3\"");
    return -1;
  }

  if(ktfnonstringq(ktastk( *isp1 + 1))) {
    *irtc = itfmessage(9, "General::wrongtype", "\"Character-string for #1\"");
    return -1;
  }

  if(ktfnonrealq(ktastk( *isp1 + 2))) {
    *irtc = itfmessage(9, "General::wrongtype", "\"Real Number for #2\"");
    return -1;
  }

  if(ktfrealq(ktastk(isp))){
    select = rtastk( isp);
  }
  else if(ktfstringq(ktastk(isp))){
    ia = ktfaddr(ktastk(isp));
#if SAD_REQUIRE_STRING_TERMINATION
    jlist(ilist(1, ia) + 1, ia + 1) = '\0';
#endif
    if(strcmp("Owner", &jlist(1, ia + 1)) == 0) {
      select = OWNER; goto l1;
    }

    if(strcmp("Waiter", &jlist(1, ia + 1)) == 0) {
      select = WAITER; goto l1;
    }

    *irtc = itfmessage(9, "General::wrongtype", "\"Character-string "
		       "\\\"Owner\\\" or \\\"Waiter\\\" for #3\"");
    return -1;
  }
  else{
    *irtc = itfmessage(9, "General::wrongtype",
		       "\"Real Number or Character-string for #3\"");
    return -1;
  }

 l1:  ia = ktfaddr(ktastk(*isp1 + 1));
#if SAD_REQUIRE_STRING_TERMINATION
  jlist(ilist(1, ia) + 1, ia + 1) = '\0';
#endif
  semaphore = rtastk(*isp1 + 2);

  status = netSemInfo(&jlist(1, ia + 1), semaphore, select, buffer);

  if(status == 0) {
    *kx = ktfstring +ktsalocb(-1, buffer);
  } else {
    *kx = ktfoper + mtfnull;
  }
  *irtc = 0;
  return 0;
}

/* SADScript function registration of NetSemaphore stuff */
#define REG	dlfunaloc
#define REG8	dlfunaloc8
int sadDefFunc_NetSemaphore(void) {
  REG8("NetSemTake",	tfNetSemTake,	3, NULL, NULL, 0);
  REG8("NetSemGive",	tfNetSemGive,	2, NULL, NULL, 0);
  REG8("NetSemInfo",	tfNetSemInfo,	3, NULL, NULL, 0);

  return 0;
}

/* End of File */
