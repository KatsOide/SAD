#include <sim/sad_signal.h>
#include <sim/sad_f2c.h>

#include <signal.h>
#include <stdio.h>

/* SIGFPE counting handler */
static volatile int count_SIGFPE = 0;	/* changed async. by signal handler */

void sad_SIGFPE_handler(int no, siginfo_t *info, void *uap) {
  count_SIGFPE += 1;
}

/* SIGHUP handler */
void sad_SIGHUP_handler(int no, siginfo_t *info, void *uap) {
  fprintf(stderr, "SIGHUP is received.(continue)\n");
}

/* SIGABRT handler */
void sad_SIGABRT_handler(int no, siginfo_t *info, void *uap) {
  fprintf(stderr, "SIGABRT is received.(continue)\n");
}

/* SIGFPE counting Fortran API */
void tclrfpe_(void) {
  count_SIGFPE = 0;
}

void tsetfpe_(integer4 *n) {
  count_SIGFPE = *n;
}

integer4 itgetfpe_(void) {
  return count_SIGFPE;
}

/* Install default signal handler for SAD */
void tsetupsig_(void) {
  struct sigaction sa;

  sa.sa_flags = SA_SIGINFO | SA_RESTART;
  sigemptyset(&sa.sa_mask);

  /* Initialize counter */
  count_SIGFPE = 0;

#if SAD_TRAP_SIGFPE	/* Install signal handler for SIGFPE */
  sa.sa_sigaction = (___sa_sigaction_t)&sad_SIGFPE_handler;
  sigaction(SIGFPE,  &sa, NULL);
#endif

#if SAD_TRAP_SIGHUP	/* Install signal handler for SIGHUP */
  sa.sa_sigaction = (___sa_sigaction_t)&sad_SIGHUP_handler;
  sigaction(SIGHUP,  &sa, NULL);
#endif

#if SAD_TRAP_SIGABRT	/* Install signal handler for SIGABRT */
  sa.sa_sigaction = (___sa_sigaction_t)&sad_SIGABRT_handler;
  sigaction(SIGABRT, &sa, NULL);
#endif
}

/* End of File */
