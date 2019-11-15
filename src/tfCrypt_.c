#include <sim/sad_api.h>
#include <sim/TFCODE.h>
#include <sim/TFCBK.h>
#include <sim/TFSTK.h>

#include <sys/types.h>
#include <unistd.h>

#if defined __has_include
#  if __has_include (<crypt.h>)
#     include <crypt.h>
#  endif
#endif

#include <pwd.h>

/* SADScript function definition */
static int Crypt(integer4 *isp1,
		 integer8 *kx,
		 integer4 *irtc) {
  integer8 ia1, ia2;
  char *encrypt;

  if(isp != *isp1 + 2) {
    *irtc = itfmessage(9, "General::narg", "\"2\"");
    return -1;
  }

  if(ktfstringq(ktastk(*isp1 + 1))) {
    *irtc = itfmessage(9, "General::wrongtype",
		       "\"Key string for #1\"");
    return -1;
  }

  if(ktfstringq(ktastk(*isp1 + 2))) {
    *irtc = itfmessage(9, "General::wrongtype" ,
		       "\"Slat string for #2\"");
    return -1;
  }

  ia1 = ktfaddr(ktastk(*isp1 + 1));
  ia2 = ktfaddr(ktastk(*isp1 + 2));
#if SAD_REQUIRE_STRING_TERMINATION
  jlist(ilist(1, ia1) + 1, ia1 + 1) = '\0';
  jlist(ilist(1, ia2) + 1, ia2 + 1) = '\0';
#endif

  encrypt = crypt(&jlist(1, ia1 + 1), &jlist(1, ia2 + 1));
  if(encrypt != NULL) {
    *kx = ktfstring + ktsalocb(-1, encrypt);
  } else {
    *kx = ktfoper + mtfnull;
  }

  *irtc = 0;
  return 0;
}

static int GetUserPassword(integer4 *isp1,
			   integer8 *kx,
			   integer4 *irtc) {
  integer8 ia;
  struct passwd *pe;

  if(isp != *isp1 + 1) {
    *irtc = itfmessage(9, "General::narg", "\"1\"");
    return -1;
  }

  if(ktfstringq(ktastk(*isp1 + 1))) {
    *irtc = itfmessage(9, "General::wrongtype", "\"User name\"");
    return -1;
  }

  ia = ktfaddr(ktastk(*isp1 + 1));
#if SAD_REQUIRE_STRING_TERMINATION
  jlist(ilist(1, ia) + 1, ia + 1) = '\0';
#endif

  if(geteuid() == 0) {
    pe = getpwnam(&jlist(1, ia + 1));
    if(pe != NULL) {
      *kx = ktfstring + ktsalocb(-1, pe->pw_passwd);
    } else {
      *kx = ktfoper + mtfnull;
    }
  } else {
    *kx = ktfstring + ktsalocb(-1, "*");
  }

  *irtc = 0;
  return 0;
}

/* SADScript function registration of Process system call stuff */
#define REG	dlfunaloc
#define REG8	dlfunaloc8
int sadDefFunc_Crypt(void) {
  REG8("Crypt",			Crypt,		 2, NULL, NULL, 0);
  REG8("GetUserPassWord",	GetUserPassword, 1, NULL, NULL, 0);

  return 0;
}

/* End of File */
