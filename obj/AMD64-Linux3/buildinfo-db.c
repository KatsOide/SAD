#define _BUILDINFO_C_
#include <buildinfo.h>
#include <stdlib.h>

const buildinfo_t buildinfo_db[] = {
  {"Repository:Owner",	_STR, 0, "sadist"},
  {"Repository:Type",	_STR, 0, "CVS"},
  {"Repository:Location",	_STR, 0, "cvs://acsad0.kek.jp/SAD/cvsroot/oldsad"},
  {"Source:TreeName",	_STR, 0, "k64-1-6 branch"},
  {"Config:ABI",	_STR, 0, ""},
  {"Config:FC",	_STR, 0, "gfortran"},
  {"Config:CC",	_STR, 0, "gcc"},
  {"Config:CXX",	_STR, 0, "c++"},
  {"Config:SYS_FOPT",	_STR, 0, "-Wall -m64 -O3"},
  {"Config:SYS_COPT",	_STR, 0, "-Wall -std=c99     -pedantic-errors -D_POSIX_C_SOURCE=200112L -D_XOPEN_SOURCE=600"},
  {"Config:SYS_CXXOPT",	_STR, 0, "-Wall -std=c++98   -pedantic-errors"},
  {"Config:SYS_CC_ABIOPT",	_STR, 0, ""},
  {"Config:SYS_FC_ABIOPT",	_STR, 0, "-Wall -std=\"gnu11\" -m64 "},
  {"Config:SYS_CXX_ABIOPT",	_STR, 0, ""},
  {"Config:FOPT",	_STR, 0, "-g -O3 -fno-second-underscore -fdollar-ok -fargument-alias  -frecursive -fbackslash -std=legacy -fall-intrinsics -Wall -m64 -fno-range-check -fbacktrace -falign-functions=64  -falign-jumps -falign-loops  -falign-labels"},
  {"Config:COPT",	_STR, 0, "-g -O3 -Wall -std=\"gnu11\" -m64 "},
  {"Config:CXXOPT",	_STR, 0, "-g -O1"},
  {"Target:OS_TYPE",	_STR, 0, "Linux"},
  {"Target:OS_NAME",	_STR, 0, "Linux"},
  {"Target:OS_MAJOR",	_INT, 3, NULL},
  {"Target:OS_MINOR",	_INT, 13, NULL},
  {"Target:OS_RELEASE",	_STR, 0, "3.13.0-62-generic"},
  {"Target:CPU_ARCH",	_STR, 0, "AMD64"},
  {"Target:MACH_ARCH",	_STR, 0, "AMD64-Linux3"},
  {"Target:SAD_ROOT",	_STR, 0, "/Users/oide/SAD/oldsad/"},
  {"Target:SAD_ARCH_ROOT",	_STR, 0, "/Users/oide/SAD/oldsad//arch"},
  {"Target:SAD_SHARE_ROOT",	_STR, 0, "/Users/oide/SAD/oldsad//share"},
  {"Target:SAD_PKG_ROOT",	_STR, 0, "/Users/oide/SAD/oldsad//share/Packages"},
  {"Target:SAD_MOD_ROOT",	_STR, 0, "/Users/oide/SAD/oldsad//share/Extension"},
  {"Built:Location",	_STR, 0, "oide@ubuntu:/home/oide/SAD/oldsad/obj/AMD64-Linux3"},
  {"Built:Date",	_STR, 0, "2015-09-02 21:12:36 +0900"},
  {NULL,	_NOP, 0, NULL}};

/* End of File */
