#include <sim/sad_api.h>
#include <sim/TFCODE.h>
#include <sim/TFSTK.h>

#include "Python.h"

static char argv0[] = "oldsad";
static int initialized = 0;

/* SADScript function definition of PyInter stuff */
int tfPyEvalString(integer4 *isp1,
		   integer4 *kx,
		   integer4 *irtc) {
  integer8 ia;
  char *argv = argv0;

  if(!initialized) { /* Initialize the Python interpreter if required. */
    initialized = 1;
    Py_Initialize();
    PySys_SetArgv(1, &argv);
  }

  if(isp != *isp1 + 1) {
    *irtc = itfmessage(9, "General::narg", "\"1\"");
    return -1;
  }

  if((ktfmask & ktastk(isp)) != ktfstring) {
    *irtc = itfmessage(9, "General::wrongtype", "\"Character-string\"");
    return -1;
  }

  ia = (ktamask & ktastk(isp));
#if SAD_REQUIRE_STRING_TERMINATION
  jlist(ilist(1, ia) + 1, ia + 1) = '\0';
#endif

  /* Execute some Python statements (in module __main__) */
  PyRun_SimpleString(&jlist(1, ia + 1));
  *kx = ktfoper + mtfnull;
  *irtc = 0;
  return 0;
}


int tfPyShutDown(integer4 *isp1,
		 integer4 *kx,
		 integer4 *irtc) {
  if(isp != *isp1 + 1
     || ktastk(isp) != ktfoper + mtfnull) {
    *irtc = itfmessage(9, "General::narg","\"0\"");
    return -1;
  }

  if(initialized) { /* Exit, cleaning up the interpreter */
    Py_Exit(0);
    initialized = 0;
  }

  *kx = ktfoper + mtfnull;
  *irtc = 0;
  return 0;
}

/* SADScript function registration of PyInter stuff */
#define REG	dlfunaloc
#define REG8	dlfunaloc8
int sadDefFunc_PyInter(void) {
  int id;

  id = REG8("PyEvalString",	tfPyEvalString,	   1, NULL, NULL, 0);
  id = REG8("PyShutDown",	tfPyShutDown,	   0, NULL, NULL, 0);

  return 0;
}

/* End of File */
