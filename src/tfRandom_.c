#include <random_driver.h>
#include <feature.h>

#include <sim/sad_api.h>
#include <sim/TFCODE.h>
#include <sim/TFCBK.h>
#include <sim/TFSTK.h>

#include <string.h>
#include <stdlib.h>

static int RandomRC(int, integer4, integer4,
		    integer8 *, integer4 *);

#define DefineRandom(type, mode, dist) \
static int Random##type(integer4 *isp1, \
			integer8 *kx,    \
			integer4 *irtc) { \
  switch(random_check_mode(mode)) { \
  case 0: \
    if(RANDOM_GAUSS == mode) { \
      random_gauss_cut = rgetgl1("GCUT"); \
      if(random_gauss_cut < RANDOM_GCUT_MIN) \
	random_gauss_cut = RANDOM_GCUT_MIN; \
    } \
    return RandomRC(mode, *isp1, isp, kx, irtc); \
    break; \
  case -2: \
    *irtc = itfmessage(9, "System::error", \
		       "\"No random plugin attached\""); \
    return -1; \
    break; \
  default: \
    *irtc = itfmessage(9, "System::error", \
		       "\"Current PRNG plugin does not support " \
		       "to generate " dist " distribution.\""); \
  } \
  return -1; \
}

/* random_driver FFS API */
DefineRandom(Gauss,	RANDOM_GAUSS,		"Gaussian")
DefineRandom(Parabola,	RANDOM_PARABOLA,	"parabolic (-1, 1)")
DefineRandom(Uniform,	RANDOM_UNIFORM,		"uniform (0,1)")
DefineRandom(Uniform0,	RANDOM_UNIFORM0,	"uniform [0,1)")
DefineRandom(Uniform1,	RANDOM_UNIFORM1,	"uniform (0,1]")
DefineRandom(Uniform01,	RANDOM_UNIFORM01,	"uniform [0,1]")
/* CAUTION: Don't append `;' after DefineRandom() macro		*/
/* ISO C does not allow extra `;' outside of a function!	*/

static int SeedRandom(integer4 *isp1,
		      integer8 *kx,
		      integer4 *irtc) {
  integer4 isp0;
  integer8 ia, ki, kt, kti, kai;
  real8 vi;
  size_t i, j, n, n0;
  uint32_t *buffer;
  const char *name, *sequence;

  if(isp != *isp1 + 1) {
    *irtc = itfmessage(9, "General::narg", "\"0 or 1\"");
    return -1;
  }

  n0 = random_size(); n = 0;
  buffer = NULL; name = NULL; sequence = NULL;
  if(ktfrealq(ktastk(*isp1 + 1))){
    buffer = malloc(sizeof(uint32_t) * n0);
    if(buffer == NULL) {
      *irtc = itfsyserr(9);
      return -1;
    }
    buffer[0] = rtastk(*isp1 + 1); n = 1;
  }
  else {
    kt = (ktfmask & ktastk(*isp1 + 1));
    if(kt == ktfoper) {
      if(ktastk(*isp1 + 1) != ktfoper + mtfnull) {
      *irtc = itfmessage(9, "General::wrongtype",
			 "\"String or Real or List or Null\"");
      return -1;
      }
    }
    else if(kt == ktfstring){
      ia = (ktamask & ktastk(*isp1 + 1));
#if SAD_REQUIRE_STRING_TERMINATION
      jlist(ilist(1, ia) + 1, ia + 1) = '\0';
#endif
      if(random_select(&jlist(1, ia + 1), NULL) != 0) {
        *irtc = itfmessage(9, "System::error",
                           "\"No such random plugin\"");
        return -1;
      }
    }
    else if(kt == ktflist) {
      ia = (ktamask & ktastk(*isp1 + 1));
      n = ilist(2, ia - 1);
      if(n > n0) n0 = n;
      buffer = malloc(sizeof(uint32_t) * n0);
      if(buffer == NULL) {
        *irtc = itfsyserr(9);
        return -1;
      };
      for(i = 0, j = 0; i < n; i++) {
        ki = klist(ia + i + 1);
        if(ktfrealq(ki)){
          buffer[j] = vi; j += 1;
        }
        else{
          kti = (ktfmask & ki);
          if(kti == ktfstring){
            kai=(ktamask & ki);
            if(j == 0) {
              if(i == 0) {
#if SAD_REQUIRE_STRING_TERMINATION
                jlist(ilist(1, kai) + 1, kai + 1) = '\0';
#endif
                name = &jlist(1, kai + 1);
                break;
              }
              if(i == 1) {
#if SAD_REQUIRE_STRING_TERMINATION
                jlist(ilist(1, kai) + 1, kai + 1) = '\0';
#endif
                sequence = &jlist(1, kai + 1);
                break;
              }
            }
          }
          else{
            *irtc = itfmessage(9, "General::wrongtype",
                               "\"Real\"");
            free(buffer);
            return -1;
          }
        }
      }
      n = j;
    }
    else{
      *irtc = itfmessage(9, "General::wrongtype",
                         "\"String or Real or List or Null\"");
      return -1;
    }
  }
  if(name != NULL) {
    if(random_select(name, sequence) != 0) {	/* Select random plugin */
      *irtc = itfmessage(9, "System::error",
			 "\"No such random plugin\"");
      free(buffer);
      return -1;
    }
  } else {
    name = random_name();			/* Check current plugin name */
  }

  sequence = random_sequence();

  if(name == NULL) {
    *irtc = itfmessage(9, "System::error",
		       "\"No random plugin attached\"");
    free(buffer);
    return -1;
  };
  if(buffer != NULL) {		/* Restore internal state from buffer */
    if(n > 0) random_restore(buffer, n);
  } 
  else {
    buffer = malloc(sizeof(uint32_t) * n0);
    if(buffer == NULL) {
      *irtc = itfsyserr(9);
      return -1;
    }
  }

  n = random_dump(buffer, n0);	/* Dump internal state to buffer */

  isp0 = isp;

  isp += 1;
  ktastk(isp) = ktfstring + ktsalocb(-1, name);

  if(sequence != NULL) {
    isp += 1;
    ktastk(isp) = ktfstring + ktsalocb(-1, sequence);
  }

  for(i = 0; i < n; i++){
    isp += 1;
    rtastk(isp) = buffer[i];
  }
  free(buffer);

  *kx = ktflist + ktfmakelist(isp0);
  *irtc = 0;
  return 0;
}

void seedrandom_(integer4 *isp1,
                 integer8 *kx,
                        integer4 *irtc){
  SeedRandom(isp1,kx,irtc);
    }


static void push_plugin_name(random_plugin_t *plugin, void *param) {
  isp += 1;
  ktastk(isp) = ktfstring + ktsalocb(-1, plugin->name);
}

static void push_sequence_name(random_plugin_t *plugin, void *param) {
  if(strcmp(plugin->name, param) == 0 && plugin->sequence != NULL) {
    isp += 1;
  ktastk(isp) = ktfstring + ktsalocb(-1, plugin->sequence);
  }
}

static int ListRandom(integer4 *isp1,
		      integer8 *kx,
		      integer4 *irtc) {
  integer8 ia, kt;
  integer4 isp0, i, j;

  if(isp != *isp1 + 1) {
    *irtc = itfmessage(9, "General::narg", "\"0 or 1\"");
    return -1;
  }

  isp0 = isp;
  kt = (ktfmask & ktastk(*isp1 + 1));
  if(kt == ktfstring){
    ia = (ktamask & ktastk(*isp1 + 1));
#if SAD_REQUIRE_STRING_TERMINATION
    jlist(ilist(1, ia) + 1, ia + 1) = '\0';
#endif
    random_scan_plugin_list(push_sequence_name, &jlist(1, ia + 1));
  }
  else if(kt == ktfoper){
    if(ktastk(*isp1 + 1) == ktfoper + mtfnull) {
      random_scan_plugin_list(push_plugin_name, NULL);
    }
  }
  else {
    *irtc = itfmessage(9, "General::wrongtype",
		       "\"Null or Character-String\"");
    return -1;
  }

  /* Reverse plugin list */
  for(i = isp0 + 1, j = isp; i < j; i += 1, j -= 1) {
    ia = ktastk(i); ktastk(i) = ktastk(j); ktastk(j) = ia;
  }

  *kx = ktflist + ktfmakelist(isp0);
  isp = isp0;
  *irtc = 0;
  return 0;
}

/* random_driver FFS API general backend */
static int RandomRC(int mode, integer4 isp1, integer4 isp2,
		    integer8 *kx,
		    integer4 *irtc) {
  integer8 kax, kt;
  integer4 isp0;
  real8 vx;
  int i, n;

  if(isp2 > isp1) {
    if(ktfrealq(ktastk(isp1 + 1))){
      n = rtastk(isp1 + 1);
      if(n < 0) {
	*irtc = itfmessage(9, "General::wrongnum", "\"positive or 0\"");
	return -1;
      };
      if(isp2 > isp1 + 1) {
	isp0 = isp;
	for(i = 0; i < n; i++) {
	  isp = isp + 1;
	  RandomRC(mode, isp1 + 1, isp2, &ktastk(isp), irtc);
	  if(*irtc != 0) {
	    isp = isp0;
	    return -1;
	  }
	}
	*kx = ktflist + ktfmakelist(isp0);
	isp = isp0;
      } else {
	kax = ktavaloc(-1, n);
	ilist(2, kax - 3) = lconstlist;
	random_generate(mode, n, &rlist(kax + 1));
        *kx = ktflist + kax;
      }
      *irtc = 0;
      return 0;
    }
    else if(ktfoperq(ktastk(isp1 + 1))){
      if(isp2 == isp1 + 1 && ktastk(isp1 + 1) == ktfoper + mtfnull) {
	random_generate(mode, 1, &vx);
        *kx = kfromr(vx);
	*irtc = 0;
	return 0;
      }
    }
    else{
      /* Nothing to do */
    }

    *irtc = itfmessage(9, "General::wrongtype", "\"Real numbers\"");
    return -1;
  }

  random_generate(mode, 1, &vx);
  *kx = kfromr(vx);
  *irtc = 0;
  return 0;
}

extern int init_framework_Random(void);

/* SADScript function registration of Random stuff */
#define REG	dlfunaloc
#define REG8	dlfunaloc8
#ifdef WITH_EXTENSION_MODULE
int dldeffun_(void) {
#else
int sadDefFunc_Random(void) {
#endif /* WITH_EXTENSION_MODULE */
  if(!feature_require("Random/Framework", FEATURE_VERSION(1, 6)))
    return -1;

  if(!feature_provide("Random/API",
		      FEATURE_VERSION(1, 0), FEATURE_VERSION(1, 0)))
    return -1;

#ifdef WITH_EXTENSION_MODULE
  REG8("SeedRandomMT",		SeedRandom,		1, NULL, NULL, 0);
  REG8("ListRandomMT",		ListRandom,		1, NULL, NULL, 0);
  REG8("RandomMT",		RandomUniform0,		1, NULL, NULL, 0);
  REG8("GaussRandomMT",		RandomGauss,		1, NULL, NULL, 0);
  REG8("ParabolaRandomMT",	RandomParabola,		1, NULL, NULL, 0);
  REG8("UniformRandomMT",	RandomUniform,		1, NULL, NULL, 0);
  REG8("UniformRandomMT0",	RandomUniform0,		1, NULL, NULL, 0);
  REG8("UniformRandomMT1",	RandomUniform1,		1, NULL, NULL, 0);
  REG8("UniformRandomMT01",	RandomUniform01,	1, NULL, NULL, 0);
#else
  REG8("SeedRandom",		SeedRandom,		1, NULL, NULL, 0);
  REG8("ListRandom",		ListRandom,		1, NULL, NULL, 0);
  REG8("Random",			RandomUniform0,		1, NULL, NULL, 0);
  REG8("GaussRandom",		RandomGauss,		1, NULL, NULL, 0);
  REG8("ParabolaRandom",		RandomParabola,		1, NULL, NULL, 0);
  REG8("UniformRandom",		RandomUniform,		1, NULL, NULL, 0);
  REG8("UniformRandom0",		RandomUniform0,		1, NULL, NULL, 0);
  REG8("UniformRandom1",		RandomUniform1,		1, NULL, NULL, 0);
  REG8("UniformRandom01",	RandomUniform01,	1, NULL, NULL, 0);
#endif /* WITH_EXTENSION_MODULE */

  return 0;
}

/* End of File */
