#include <sim/sad_api.h>
#include <sim/dynl.h>

#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/* Base offset number for itfunaloc_ */
#define FUNC_ID_OFFSET	4001

/* SADScript function pointer table */
#ifndef MAX_TBL_FUNCS
#define MAX_TBL_FUNCS	1000
#endif
static sad_func_t dl_func_tbl[MAX_TBL_FUNCS];
static bool dl_func_tbl_initialized = false;

static void dl_func_tbl_init(void) {
  if(!dl_func_tbl_initialized) {
    for(int i = 0; i < MAX_TBL_FUNCS; i++)
      dl_func_tbl[i] = NULL;
    dl_func_tbl_initialized = true;
  }
}

static int dl_reg_func(sad_func_t func) {
  int i = 0;

  for(i = 0; i < MAX_TBL_FUNCS; i++)
    if(dl_func_tbl[i] == NULL) {
      dl_func_tbl[i] = func;
      return i;
    }
  return -1;
}

/* SADScript function eval interface invoked by tfefun_@src/tfefun.f */
int tfefunctbl_(integer8 *isp1, integer4 *id_ref,
		integer4 *itx, integer4 *iax, real8 *vx,
		integer4 *itrc) {
  sad_func_t func;
  int id = *id_ref - FUNC_ID_OFFSET;

  if(MAX_TBL_FUNCS > id && id >=0) {
    func = dl_func_tbl[id];
    if(func != NULL)
      return func(isp1, itx, iax, vx, itrc);
  }

  return -1;
}

/* SADScript function eval interface invoked by tfefun_@src/tfefun.f */
int tfefunctbl8_(integer8 *isp1, integer4 *id_ref,
                 integer8 *kx, integer4 *itrc) {
  sad_func_t_8 func;
  int id = *id_ref - FUNC_ID_OFFSET;

  if(MAX_TBL_FUNCS > id && id >=0) {
    func = (sad_func_t_8) dl_func_tbl[id];
    if(func != NULL)
      return func(isp1, kx, itrc);
  }

  return -1;
}

/* SADScript function registration interface */
int dlfunaloc(const char *name, sad_func_t func,
	      int narg, const int *map, const int *ieval, int immed) {
  integer4 id, nfun, narg_ = narg, immed_ = immed;
  integer4 *map_, *ieval_;
  int i;

  map_ = malloc(2 * sizeof(integer4) * (narg + 1));
  if(map_ == NULL) {
    fprintf(stderr, "dlfunaloc: map/ieval duplication is failed\n");
    exit(1);
  }

  if(map != NULL) {
    for(i = 0; i < narg + 1; i++) map_[i] = map[i];
  } else {
    for(i = 0; i < narg + 1; i++) map_[i] = 0;
  }

  ieval_ = map_ + (narg + 1);
  if(ieval != NULL) {
    for(i = 0; i < narg + 1; i++) ieval_[i] = ieval[i];
  } else {
    for(i = 0; i < narg + 1; i++) ieval_[i] = 0;
  }

  /* Assign SADScript function pointer table entry */
  id = dl_reg_func(func);
  if(id < 0) {
    fprintf(stderr, "dlfunaloc: dynamic function table is over flowed\n");
    exit(1);
  }

  /* Register to SAD interpreter */
  id += FUNC_ID_OFFSET;
  nfun = itfunaloc_(name, &id, &narg_, map_, ieval_, &immed_, strlen(name));

  free(map_);
  return nfun;
}

int dlfunaloc8(const char *name, sad_func_t_8 func,
	      int narg, const int *map, const int *ieval, int immed) {
  integer4 id, nfun, narg_ = narg, immed_ = immed;
  integer4 *map_, *ieval_;
  int i;

  map_ = malloc(2 * sizeof(integer4) * (narg + 1));
  if(map_ == NULL) {
    fprintf(stderr, "dlfunaloc: map/ieval duplication is failed\n");
    exit(1);
  }

  if(map != NULL) {
    for(i = 0; i < narg + 1; i++) map_[i] = map[i];
  } else {
    for(i = 0; i < narg + 1; i++) map_[i] = 0;
  }

  ieval_ = map_ + (narg + 1);
  if(ieval != NULL) {
    for(i = 0; i < narg + 1; i++) ieval_[i] = ieval[i];
  } else {
    for(i = 0; i < narg + 1; i++) ieval_[i] = 0;
  }

  /* Assign SADScript function pointer table entry */
  id = dl_reg_func((sad_func_t) func);
  if(id < 0) {
    fprintf(stderr, "dlfunaloc: dynamic function table is over flowed\n");
    exit(1);
  }

  /* Register to SAD interpreter */
  id += FUNC_ID_OFFSET;
  nfun = itfunaloc_(name, &id, &narg_, map_, ieval_, &immed_, strlen(name));

  free(map_);
  return nfun;
}

/* Fortran version is provided only for compatibility! */
integer4 dlfunaloc_(const_character name, sad_func_t func, integer4 *narg,
		    integer4 *map, integer4 *ieval, integer4 *immed,
		    ftnlen len_name) {
  integer4 id;

  /* Assign SADScript function pointer table entry */
  id = dl_reg_func(func);
  if(id < 0) {
    fprintf(stderr, "dlfunaloc: dynamic function table is over flowed\n");
    exit(1);
  }

  /* Register to SAD interpreter */
  id += FUNC_ID_OFFSET;
  return itfunaloc_(name, &id, narg, map, ieval, immed, len_name); 
}

/* Initialization hook */
int init_framework_FuncTBL(void) {
  dl_func_tbl_init();
  return 0;
}

/* End of File */
