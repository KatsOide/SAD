#ifndef _SAD_SIGNAL_H_
#define _SAD_SIGNAL_H_

#include <sim/sad_f2c.h>
#include <signal.h>

#ifndef SAD_TRAP_SIGFPE
#define SAD_TRAP_SIGFPE		1
#endif
#ifndef SAD_TRAP_SIGHUP
#define SAD_TRAP_SIGHUP		0
#endif
#ifndef SAD_TRAP_SIGABRT
#define SAD_TRAP_SIGABRT	0
#endif

/* constat for SigAction */
#define SAD_SA_ONSTACK	  0x0001  /* take signal on signal stack */
#define SAD_SA_RESTART	  0x0002  /* restart system call on signal return */
#define SAD_SA_RESETHAND  0x0004  /* reset to SIG_DFL when taking signal */
#define SAD_SA_NOCLDSTOP  0x0008  /* do not generate SIGCHLD on child stop */
#define SAD_SA_NODEFER	  0x0010  /* don't mask the signal we're delivering */
#define SAD_SA_NOCLDWAIT  0x0020  /* don't keep zombies around */
#define SAD_SA_SIGINFO	  0x0040  /* signal handler with SA_SIGINFO args */

/* prototypr for strcut sigaction.sa_sigaction */
typedef void (*___sa_sigaction_t)(int, siginfo_t *, void *);

/* Maximum signal number */
#define SAD_MAXSIG	(sizeof(sigset_t) * 8)

/* For signal action on SADScript */
#define SAD_SIGACT_DFL	0

/* SAD internal signal handler */
extern void sad_SIGFPE_handler(int, siginfo_t *, void *);
extern void sad_SIGHUP_handler(int, siginfo_t *, void *);
extern void sad_SIGABRT_handler(int, siginfo_t *, void *);

/* SAD Fortran interfaces */
extern void tsetupsig_(void);
extern void tclrfpe_(void);
extern void tsetfpe_(integer4*);
extern integer4 itgetfpe_(void);

#endif /* _SAD_SIGNAL_H_ */
