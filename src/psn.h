#ifndef PSN_H
#define PSN_H
/* psn.h*/
void calc_potential_psn(double *gridcharge, double *gridphi);
void create_env_psn(int nx, int ny, double cx, double cy);
void clear_env_psn(void);
#endif
