/*
 * tfProcess_.c contains wrong code casting pointer between object pointer
 * and function pointer. It is used to dump and restore between a function
 * pointer and "0x*****" hexadecimal address string by using sprintf() and
 * sscanf(). Some code blocks ``Violate ISO C standard''.
 */

#include <sim/sad_signal.h>
#include <sim/sad_api.h>
#include <sim/TFCODE.h>
#include <sim/TFSTK.h>
#include <sim/TFCBK.h>

#include <sys/types.h>
#include <sys/wait.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <sys/param.h>
#include <string.h>
#include <limits.h>
#include <signal.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdarg.h>
#include <ctype.h>
#include <fcntl.h>
#include <errno.h>
#include <stdio.h>
#include <time.h>
#include <math.h>

#define	USE_NULL_STRING_RESCUE_HACK

#ifndef	SAD_BSHELL_PATH
# if defined(__FreeBSD__) || defined(__NetBSD__) || defined(__OpenBSD__) || defined(__APPLE_CC__)
#  include <paths.h>
#  define	SAD_BSHELL_PATH		_PATH_BSHELL
# else
#  if defined(__sun__)
#   include <userdefs.h>
#   define	SAD_BSHELL_PATH		DEFSHL
#  else
#   define	SAD_BSHELL_PATH		"/bin/sh"
#  endif
# endif
#endif	/* SAD_BSHELL_PATH */

#ifndef HAVE_BSD_SETENV
#if defined(__FreeBSD__) || defined(__NetBSD__) || defined(__OpenBSD__) || defined(__APPLE_CC__)
#define HAVE_BSD_SETENV 1
#else
#define HAVE_BSD_SETENV 0
#endif
#endif	/* HAVE_BSD_SETENV */

#define True	(0 == 0)
#define False	(0 != 0)

/* Signal action number database */
static volatile sig_atomic_t sad_signal_handler_depth = 0;
static volatile sig_atomic_t sad_signal_actions[SAD_MAXSIG];
static volatile sig_atomic_t sad_signal_flags[SAD_MAXSIG];

	/*				   111111111122222222 */
	/* Template String	 0123456789012345678901234567 */
#define SIG_ACT_COMMAND		"   SigAction$Command[%d]=."
#define SIG_ACT_APPENDTO	"AppendTo[SigAction$IDs,%d]"
#define	SIG_ACT_CMD		"SigAction$Cmd [%d]"
#define SIG_ACT_RESET_CMD	"SigAction$CmdR[%d]"
#define SIG_ACT_MARK_CMD	"SigAction$CmdM[%d]"
#define SIG_ACT_CMD_OFFSET	13
#define SIG_SCRIPT_ACTION	"$SignalAction"
#define SIG_MINIMUM_ACTION  SIG_SCRIPT_ACTION"[%d,Signal->%d]"
#define	SIG_SCRIPT_BUF_SIZE	96

static void init_signal_handler(void) {
  int i;

  for(i = 0; i < SAD_MAXSIG; i++) {
    sad_signal_actions[i] = SAD_SIGACT_DFL;
    sad_signal_flags[i] = 0;
  }
}

static integer4 save_stk_frame(void) {
  integer4 isp0 = isp;	/* Save current SAD stk frame pointer for rewinding */

  isp += 2;		/* Make dump area for current stk frame context */

  /* Dump current stk frame context */
  itastk(1, isp0 + 1) = isporg;
  itastk(2, isp0 + 1) = 0;
  itastk(1, isp0 + 2) = ipurefp;
  itastk(2, isp0 + 2) = napuref;

  isporg  = isp;	/* Setup origin of next SAD stack frame */
  ipurefp = 0;
  napuref = 0;

  return isp0;
}

static void restore_stk_frame(integer4 isp0) {
  /* Restore previous stk frame context */
  isporg  = itastk(1, isp0 + 1);
  ipurefp = itastk(1, isp0 + 2);
  napuref = itastk(2, isp0 + 2);
  isp = isp0;		/* Rewind SAD stk frame pointer */
}

/*
 * At first, buf_printf try to store formated string into buf.
 * If buf is too short to srore formated string,
 * buf_printf try to allocate string buffer and store formated string into it.
 * If buf_printf succeeded, buf_printf returns current buffer size
 * and current buffer point into *ret.
 * If buf_printf failed, buf_printf returns -1 and store NULL pointer into *ret.
 */
static int buf_printf(char **ret, char *buf, size_t len, const char *fmt, ...) {
  va_list ap;
  int slen,  blen;

  blen = len;

  va_start(ap, fmt);
  slen = vsnprintf(buf, blen, fmt, ap);
  va_end(ap);

  if(slen < blen) {
    *ret = buf;
    return blen;
  }

  blen = slen + 1;
  *ret = (char *)malloc(sizeof(char) * blen);
  if(*ret == NULL) return -1;

  va_start(ap, fmt);
  slen = vsnprintf(*ret, blen, fmt, ap);
  va_end(ap);

  if(slen >= blen) {
    free(*ret);
    *ret = NULL;
    return -1;
  }

  return blen;
}

/* SADScript signal handler */
static void sad_SIGACT_handler(int no, siginfo_t *info, void *uap) {
  integer8 kx;
  integer4 isp0, irtc;
  struct sigaction sa;
  sigset_t old, full;
  int sigactno = SAD_SIGACT_DFL, sigact_flags = 0, blen;
  char *buf, buf0[SIG_SCRIPT_BUF_SIZE];

  if(info->si_signo > 0 && info->si_signo <= SAD_MAXSIG) {
    sigactno = sad_signal_actions[info->si_signo - 1];
    sigact_flags = sad_signal_flags[info->si_signo - 1];
  }
  sigfillset(&full);

  switch(info->si_signo) {
  case SIGCHLD:
    blen = buf_printf(&buf, buf0, sizeof(buf0), SIG_SCRIPT_ACTION "[%d,"
		      "Signal->%d,Error->%d,Code->%d,"
		      "UID->%d,Status->%d]",
		      sigactno,
		      info->si_signo, info->si_errno, info->si_code,
		      info->si_uid, info->si_status);
    break;

  case SIGILL:
  case SIGFPE:
  case SIGBUS:
  case SIGSEGV:
    blen = buf_printf(&buf, buf0, sizeof(buf0), SIG_SCRIPT_ACTION "[%d,"
		      "Signal->%d,Error->%d,Code->%d,"
		      "Address->\"%p\"]",
		      sigactno,
		      info->si_signo, info->si_errno, info->si_code,
		      info->si_addr);
    break;

  default:
    if(info->si_code == 0)
      blen = buf_printf(&buf, buf0, sizeof(buf0), SIG_SCRIPT_ACTION "[%d,"
			"Signal->%d,Error->%d,Code->%d,"
			"Sender->{PID->%d,UID->%d}]",
			sigactno,
			info->si_signo, info->si_errno, info->si_code,
			info->si_pid, info->si_uid);
    else
      blen = buf_printf(&buf, buf0, sizeof(buf0), SIG_SCRIPT_ACTION "[%d,"
			"Signal->%d,Error->%d,Code->%d]",
			sigactno,
			info->si_signo, info->si_errno, info->si_code);
    break;
  }

  sigprocmask(SIG_BLOCK, &full, &old);	/* Enter Signal Critical Section */
  isp0 = save_stk_frame();	/* Save SAD stk frame */
  sad_signal_handler_depth += 1;
  if(sigact_flags & SA_RESETHAND)
    switch(info->si_signo) {
#if SAD_TRAP_SIGFPE	/* Re-install signal handler for SIGFPE */
    case SIGFPE:
      sa.sa_sigaction = (___sa_sigaction_t)&sad_SIGFPE_handler;
      sa.sa_flags = SA_SIGINFO | SA_RESTART;
      sigemptyset(&sa.sa_mask);
      sigaction(SIGFPE,  &sa, NULL);
      break;

#endif
#if SAD_TRAP_SIGHUP	/* Re-install signal handler for SIGHUP */
    case SIGHUP:
      sa.sa_sigaction = (___sa_sigaction_t)&sad_SIGHUP_handler;
      sa.sa_flags = SA_SIGINFO | SA_RESTART;
      sigemptyset(&sa.sa_mask);
      sigaction(SIGHUP,  &sa, NULL);
      break;

#endif
#if SAD_TRAP_SIGABRT	/* Re-install signal handler for SIGABRT */
    case SIGABRT:
      sa.sa_sigaction = (___sa_sigaction_t)&sad_SIGABRT_handler;
      sa.sa_flags = SA_SIGINFO | SA_RESTART;
      sigemptyset(&sa.sa_mask);
      sigaction(SIGABRT,  &sa, NULL);
      break;

#endif
    default:
      break;
    }
  sigprocmask(SIG_SETMASK, &old, NULL);	/* Leave Signal Critical Section */

  if(buf != NULL) {
    tfevals(buf, &kx, &irtc);
    if(sigact_flags & SA_RESETHAND)	/* Make SIGACT marking script */
      snprintf(buf, blen, SIG_ACT_MARK_CMD, sigactno);
  }

  sigprocmask(SIG_BLOCK, &full, &old);	/* Enter Signal Critical Section */
  if(buf != NULL && sigact_flags & SA_RESETHAND)
    tfevals(buf, &kx, &irtc); /* Mark unused SIGACT Number */
  sad_signal_handler_depth -= 1;
  restore_stk_frame(isp0);	/* Rewind SAD stk frame */
  sigprocmask(SIG_SETMASK, &old, NULL);	/* Leave Signal Critical Section */

  if(buf != buf0) free(buf);
}

/* SADScript function definition */
static int Fork(integer4 *isp1, integer8 *kx,
		integer4 *irtc) {
  real8 vx;

  if(isp != *isp1 + 1
     || ktastk(isp) != ktfoper + mtfnull) {
    *irtc = itfmessage(999, "General::narg", "\"0\"");
    return -1;
  }

  vx = fork_worker_();
  *kx = kfromr(vx);
  *irtc = 0;
  return 0;
}

static int Wait(integer4 *isp1, integer8 *kx,
		integer4 *irtc) {
  integer8 kax;
  pid_t pid;
  int status;

  if(isp != *isp1 + 1
     || ktastk(isp) != ktfoper + mtfnull) {
    *irtc = itfmessage(999, "General::narg", "\"0\"");
    return -1;
  }

  pid = wait(&status);

  if(pid > 0) {
    kax = ktavaloc(-1, 2);
    ilist(2, kax - 3) = lconstlist;
    rlist(kax + 1) = pid;
    rlist(kax + 2) = status;
    *kx = ktflist + kax;
    *irtc = 0;
  } else
    *irtc = itfsyserr(9);

  return 0;
}

static int Wait4(integer4 *isp1, integer8 *kx,
		 integer4 *irtc) {
  static integer8 itfContinued = 0, itfExited = 0, itfSignaled = 0,
    itfStopped = 0, itfSignal = 0, itfCore = 0, itfUserTime, itfSysTime;
  integer4 isp0;
  pid_t pid, wpid;
  int status, options, options0, i;
  struct rusage resource_usage;

  if(itfContinued == 0) {
    itfContinued = ktfsymbolf("Continued",  True);
    itfExited    = ktfsymbolf("Exited",     True);
    itfSignaled  = ktfsymbolf("Signaled",   True);
    itfStopped   = ktfsymbolf("Stopped",    True);
    itfSignal    = ktfsymbolf("Signal",     True);
    itfCore      = ktfsymbolf("Core",       True);
    itfUserTime  = ktfsymbolf("UserTime",   True);
    itfSysTime   = ktfsymbolf("SystemTime", True);
  }

  if(isp == *isp1 + 1
     && ktastk(isp) == ktfoper + mtfnull) {
    wpid = -1; /* wait any process: wait() */
  } else if(ktfrealq(ktastk(*isp1 + 1))) {
    wpid = rtastk(*isp1 + 1);
  } else {
    *irtc = itfmessage(9, "General::wrongtype",
		       "\"wpid for #1\"");
    return -1;
  }

  options0 = 0;
  for(i = 2; *isp1 + i <= isp; i++) {
    if(ktfnonrealq(ktastk(*isp1 + i))) {
      *irtc = itfmessage(9, "General::wrongtype",
			 "\"Wait flag for ##2\"");
      return -1;
    }
    options0 |= (int)rtastk(*isp1 + i);
  }

  options = 0;
  if(options0 & 0x0001) options |= WNOHANG;
  if(options0 & 0x0002) options |= WUNTRACED;
#ifdef WCONTINUED
  if(options0 & 0x0004) options |= WCONTINUED;
#endif

  pid = wait4(wpid, &status, options, &resource_usage);

  isp0 = isp;
  switch(pid) {
  case 0: /* No children to wait(WNOHANG) */
    isp += 1;
    rtastk(isp) = 0;
    isp += 1;
    rtastk(isp) = 0;
    break;

  case -1:
    switch(errno) {
    case ECHILD: /* No children to wait */
      isp += 1;
      rtastk(isp) = 0;
      isp += 1;
      rtastk(isp) = 0;
      break;

    case EINTR: /* wait4() is interrupted */
      isp += 1;
      rtastk(isp) = -1;
      isp += 1;
      rtastk(isp) = 0;
      break;

    default:
      *irtc = itfsyserr(9);
      return 1;
      break;
    }
    break;

  default:
    isp += 1;
    rtastk(isp) = pid;
    isp += 1;
    rtastk(isp) = WIFEXITED(status) ? WEXITSTATUS(status) : 0;

    /* Process status bits... */
#if defined(WCONTINUED) && defined(WIFCONTINUED)
    tfmakerulestk(ktfsymbol+itfContinued, SAD_BOOLEAN(WIFCONTINUED(status)));
#endif
    tfmakerulestk(ktfsymbol+itfExited,    SAD_BOOLEAN(WIFEXITED(status)));
    tfmakerulestk(ktfsymbol+itfSignaled,  SAD_BOOLEAN(WIFSIGNALED(status)));
    tfmakerulestk(ktfsymbol+itfStopped,   SAD_BOOLEAN(WIFSTOPPED(status)));

    if(WIFSIGNALED(status)) { /* if signaled */
#if defined(WTERMSIG)
      tfmakerulestk(ktfsymbol+itfSignal,  WTERMSIG(status));
#endif
#if defined(WCOREDUMP)
      tfmakerulestk(ktfsymbol+itfCore,    SAD_BOOLEAN(WCOREDUMP(status)));
#endif
    }

    if(WIFSTOPPED(status)) { /* if stopped */
#if defined(WSTOPSIG)
      tfmakerulestk(ktfsymbol+itfSignal,  WSTOPSIG(status));
#endif
    }

    /* Process resource usage summary...(see getrusage(2)) */
    tfmakerulestk(ktfsymbol+itfUserTime,  resource_usage.ru_utime.tv_sec
		  + resource_usage.ru_utime.tv_usec * 1e-6);
    tfmakerulestk(ktfsymbol+itfSysTime,   resource_usage.ru_stime.tv_sec
		  + resource_usage.ru_stime.tv_usec * 1e-6);

  }

  *kx = ktflist + ktfmakelist(isp0);
  isp = isp0;
  *irtc = 0;
  return 0;
}

static int Kill(integer4 *isp1, integer8 *kx,
		integer4 *irtc) {
  pid_t pid;
  int sig;

  if(isp != *isp1 + 1 && isp != *isp1 + 2) {
    *irtc = itfmessage(9, "General::narg", "\"1 or 2\"");
    return -1;
  }

  if(ktfnonrealq(ktastk(*isp1 + 1))) {
    *irtc = itfmessage(9, "General::wrongtype", "\"Process ID for #1\"");
    return -1;
  }
  pid = rtastk(*isp1 + 1);

  if(isp == *isp1 + 2) {
    if(ktfnonrealq(ktastk(*isp1 + 2))) {
      *irtc = itfmessage(9, "General::wrongtype", "\"Signal Number for #2\"");
      return -1;
    }
    sig = rtastk(*isp1 + 2);
  } else
    sig = SIGTERM;

  if(kill(pid, sig) != 0)
    switch(errno) {
    case EINVAL:
      *irtc = itfmessage(9, "General::wrongtype", "\"Signal Number for #2\"");
      return -1;
      break;;

    default:
      *irtc = itfsyserr(9);
      return 1;
    }

  *kx = ktfoper + mtfnull;
  *irtc = 0;
  return 0;
}

static int Sleep(integer4 *isp1, integer8 *kx,
		 integer4 *irtc) {
  int status;
  struct timespec time_to_sleep, time_slept;
  real8 timeout;

  if(isp != *isp1 + 1) {
    *irtc = itfmessage(9, "General::narg", "\"1\"");
    return -1;
  }

  if(ktfnonrealq(ktastk(isp))) {
    *irtc = itfmessage(9, "General::wrongtype", "\"Real number (seconds)\"");
    return -1;
  }

  timeout = rtastk(isp);

  while(timeout > 0.5e-9) {
    if(timeout > INT_MAX) {
      time_to_sleep.tv_sec  = INT_MAX;
      time_to_sleep.tv_nsec = 0;
      timeout -= INT_MAX;
    } else {
      time_to_sleep.tv_sec  = floor(timeout); timeout -= time_to_sleep.tv_sec;
      time_to_sleep.tv_nsec = ceil(1e9 * timeout);
      timeout = 0.0;
    }

    status = nanosleep(&time_to_sleep, &time_slept);
    if(status != 0) {
#ifdef COMPAT_SAD_SLEEP_OLDAPI
      *kx = ktfoper + mtfnull;
#else
      *kx  = 0;
#endif
      *irtc = 0;
      return 0;
    }
  }

#ifdef COMPAT_SAD_SLEEP_OLDAPI
  *kx = ktfoper + mtfnull;
#else
  *kx = kfromr(r_true);
#endif
  *irtc = 0;
  return 0;
}

static int SigPending(integer4 *isp1, integer8 *kx,
		      integer4 *irtc) {
  integer4 isp0;
  sigset_t pending;
  int signo, status;

  if(isp != *isp1 + 1
     || ktastk(isp) != ktfoper + mtfnull) {
    *irtc = itfmessage(9, "General::narg", "\"0\"");
    return -1;
  }

  if(sigpending(&pending) != 0) {
    *irtc = itfsyserr(9);
    return 1;
  }

  isp0 = isp;
  for(signo = 1; signo <= SAD_MAXSIG && (status = sigismember(&pending, signo)) >= 0; signo++)
    if(status != 0) {
      isp += 1;
      rtastk(isp) = signo;
    }

  *kx = ktflist + ktfmakelist(isp0);
  isp = isp0;
  *irtc = 0;
  return 0;
}

static int SigSuspend(integer4 *isp1, integer8 *kx,
		      integer4 *irtc) {
  integer8 ia;
  integer4 isp2;
  sigset_t mask;
  int signo, i;

  sigemptyset(&mask);

  if(ktfrealq(ktastk(*isp1 + 1))){
    signo = rtastk(*isp1 + 2);
    sigaddset(&mask, signo);
  }else if(!(ktastk(*isp1 + 1) == ktfoper + mtfnull && *isp1 + 1 == isp)) {
      *irtc = itfmessage(9, "General::wrongtype",
			 "\"Null or Signal-number List or Signal-number for #1\"");
      return -1;
  } else if(ktflistq(ktastk(*isp1 + 1))){
    ia = ktfaddr(ktastk(*isp1 + 1));
    if(klist(ia) == ktfoper + mtflist
       && ktfreallistq(ia)){
      for(i = 1; i <= ilist(2, ia - 1); i++) {
	signo = rlist(ia + i);
	sigaddset(&mask, signo);
      }
    } else {
      *irtc = itfmessage(9, "General::wrongtype",
			 "\"Null or Signal-number List or Signal-number for #1\"");
      return -1;
    }
  } else {
    *irtc = itfmessage(9, "General::wrongtype",
		       "\"Null of Signal-number List or Signal-number for #1\"");
    return -1;
  }

  for(isp2 = *isp1 + 2; isp2 <= isp; isp2++)
    if(ktfrealq(ktastk(isp2))) {
      signo = rtastk(isp2);
      sigaddset(&mask, signo);
    } else {
      *irtc = itfmessage(9, "General::wrongtype",
			 "\"Signal-number for ##2\"");
      return -1;
    }

  sigsuspend(&mask);

  *kx = ktfoper + mtfnull;
  *irtc = 0;
  return 0;
}

static int SigProcMask(integer4 *isp1, integer8 *kx,
		       integer4 *irtc) {
  static integer8 itfUnBlock = 0, itfAll = 0;
  integer4 isp0, isp2, i, iop; 
  integer8 ia;
  sigset_t old, new, *mask;
  int action, signo, status;

  if(itfUnBlock == 0) {
    itfUnBlock = ktfsymbolf("UnBlock",    True);
    itfAll     = ktfsymbolf("AllSignals", True);
  }

  sigemptyset(&new);
  mask = &new;;
  action = SIG_SETMASK;

  if(ktfoperq(ktastk(*isp1 + 1))) {
    iop=ktfaddr(ktastk(*isp1 + 1));
    switch(iop) {
    case mtfnull:
      if(*isp1 + 1 == isp) {
	action = SIG_SETMASK;
	mask = NULL;
      } else {
	*irtc = itfmessage(9, "General::wrongtype",
			   "\"Block, UnBlock and Set for #1\"");
	return -1;
      }
      break;

    case mtfset:	/* Set -> (ntfoper, mtfset) */
      action = SIG_SETMASK;
      break;

    default:
      if(iop <= mtfnopc) {
	*irtc = itfmessage(9, "General::wrongtype",
			   "\"Block, UnBlock and Set for #1\"");
	return -1;
      }
      switch(ilist(2, klist(ifunbase + iop) + 1)) {
      case nfunblock:	/* Block -> (ntfoper, ifunbase + ia > mtfnopc :> nfunblock) */
	action = SIG_BLOCK;
	break;

      default:
	*irtc = itfmessage(9, "General::wrongtype",
			   "\"Block, UnBlock and Set for #1\"");
	return -1;
	break;
      }
      break;
    }
  } else if (ktfsymbolq(ktastk(*isp1 + 1))){
    if(tfsamesymbolqk_(&ktastk(*isp1 + 1), &itfUnBlock))
      action = SIG_UNBLOCK;
    else {
      *irtc = itfmessage(9, "General::wrongtype",
			 "\"Block, UnBlock and Set for #1\"");
      return -1;
    }
  } else {
    *irtc = itfmessage(9, "General::wrongtype",
		       "\"Block, UnBlock and Set for #1\"");
    return -1;
  }

  if(*isp1 + 1 < isp) {
    if(ktfrealq(ktastk(*isp1 + 2))){
      signo = rtastk(*isp1 + 2);
      sigaddset(&new, signo);
    }else {
      ia = ktfaddr(ktastk(*isp1 + 2));
      if(ktflistq(ktastk(*isp1 + 2))){
        if(klist(ia) == ktfoper + mtflist
           && ktfreallistq(ia)) {
          for(i = 1; i <= ilist(2, ia - 1); i++) {
            signo = rlist(ia + i);
            sigaddset(&new, signo);
          }
        } else {
          *irtc = itfmessage(9, "General::wrongtype",
                             "\"Signal-number List or Signal-number or AllSignals for #2\"");
          return -1;
        }
      }  else if (ktfsymbolq(ktastk(*isp1 + 2))){
        if(tfsamesymbolqk_(&ia, &itfAll)) {
          sigfillset(&new);
        } else {
          *irtc = itfmessage(9, "General::wrongtype",
                             "\"Signal-number List or Signal-number or AllSignals for #2\"");
          return -1;
        }
      }else
        *irtc = itfmessage(9, "General::wrongtype",
                           "\"Signal-number List or Signal-number or AllSignals for #2\"");
      return -1;
    }

    for(isp2 = *isp1 + 3; isp2 <= isp; isp2++)
      if(ktfrealq(ktastk(isp2))) {
	signo = rtastk(isp2);
	sigaddset(&new, signo);
      } else {
	*irtc = itfmessage(9, "General::wrongtype",
			   "\"Signal-number for ##3\"");
	return -1;
      }
  }

  if(sigprocmask(action, mask, &old) != 0) {
    *irtc = itfsyserr(9);
    return 1;
  }

  isp0 = isp;
  for(signo = 1; signo <= SAD_MAXSIG && (status = sigismember(&old, signo)) >= 0; signo++)
    if(status != 0) {
      isp += 1;
      rtastk(isp) = signo;
    }

  *kx = ktflist + ktfmakelist(isp0);
  isp = isp0;
  *irtc = 0;
  return 0;
}

static int SigAction(integer4 *isp1, integer8 *kx,
		     integer4 *irtc) {
  static const char *options[] = {"Action", "Flag", "Mask", NULL};
  static integer8 itfAction = 0, itfFlag = 0, itfMask = 0,
    itfDefault = 0, itfIgnore = 0, itfAll = 0;
  integer8 ia;
  integer4 ispopt, isp0, isp2;
  struct sigaction next_sa, old_sa, *new_sa;
  sigset_t full_mask, old_mask;
  char *buf, buf0[32];
  int new_sigactno, old_sigactno, install_sigactno;
  int signo, sad_sa_flags, status, i;

  if(itfAction == 0) {
    itfAction  = ktfsymbolf("Action", True);
    itfFlag    = ktfsymbolf("Flag",   True);
    itfMask    = ktfsymbolf("Mask",   True);

    itfDefault = ktfsymbolf("Default", True);
    itfIgnore  = ktfsymbolf("Ignore",  True);

    itfAll     = ktfsymbolf("AllSignals", True);
  }

  install_sigactno = False;
  new_sigactno = SAD_SIGACT_DFL;
  sad_sa_flags = 0;
  next_sa.sa_handler = SIG_DFL;
  next_sa.sa_flags = 0;
  sigemptyset(&next_sa.sa_mask);

  new_sa = &next_sa;

  isp0 = isp;
  ispopt = itfgetoptionstk(*isp1, options);
  if(ispopt > 0) {
    if(ktfrealq(ktastk(isp0 + 1))){
      install_sigactno = True;
      new_sigactno = rtastk(isp0 + 1);
      next_sa.sa_sigaction = (___sa_sigaction_t)sad_SIGACT_handler;
      sad_sa_flags |= SAD_SA_SIGINFO;
      }
    else if (ktastk(isp0 + 1) == ktfref){
        new_sa = NULL;	/* Read-out handler if Action is not defined */
      }
    else if(ktfsymbolq(ktastk(isp0 + 1))){
      if(tfsamesymbolqk_(&ktastk(isp0 + 1), &itfDefault))
        next_sa.sa_handler = SIG_DFL;
      else if(tfsamesymbolqk_(&ktastk(isp0 + 1), &itfIgnore))
        next_sa.sa_handler = SIG_IGN;
      else {
        isp = isp0;
        *irtc = itfmessage(9, "General::wrongtype",
                           "\"Action -> Default, Ignore, Real-number and Pointer address string\"");
        return -1;
      }
    }
    else if(ktfstringq(ktastk(isp0 + 1))){
      ia = ktfaddr(ktastk(isp0 + 1));
#if SAD_REQUIRE_STRING_TERMINATION
      jlist(ilist(1, ia) + 1, ia + 1) = '\0';
#endif
      if(ilist(1, ia) > 2 &&
         jlist(1, ia + 1) == '0' && jlist(2, ia + 1) == 'x') {
        for(i = 3; i <= ilist(1, ia); i++)
          if(!isxdigit((unsigned char)jlist(i, ia + 1))) {
            isp = isp0;
            *irtc = itfmessage(9, "General::wrongtype",
                               "\"Action -> Default, Ignore, Real-number and Pointer address string\"");
            return -1;
          }
        if(sscanf(&jlist(1, ia + 1), "%p", /* Violate ISO C standard */
                  (void **)&next_sa.sa_sigaction) == 1) goto j1;
      }
      isp = isp0;
      *irtc = itfmessage(9, "General::wrongtype",
                         "\"Action -> Default, Ignore, Real-number and Pointer address string\"");
      return -1;
    }
    else{
      isp = isp0;
      *irtc = itfmessage(9, "General::wrongtype",
			 "\"Action -> Default, Ignore, Real-number and Pointer address string\"");
      return -1;
    };
      
  j1:    if(ktfrealq(ktastk(isp0 + 2))){
      signo = rtastk(isp0 + 2);
      sad_sa_flags |= signo;
    }
    else if(ktastk(isp0 + 2) == ktfref){
    }
    else if(ktflistq(ktastk(isp0 + 2))){
      ia = ktfaddr(ktastk(isp0 + 2));
      if(klist(ia) == ktfoper + mtflist
         && ktfreallistq(ia)) {
        for(i = 1; i <= ilist(2, ia - 1); i++) {
          signo = rlist(ia + i);
          sad_sa_flags |= signo;
        }
      } 
      else {
        isp = isp0;
        *irtc = itfmessage(9, "General::wrongtype",
			   "\"Flag -> {Signal-action-flags}\"");
        return -1;
      }
    }
    else {
      isp = isp0;
      *irtc = itfmessage(9, "General::wrongtype",
			 "\"Flag -> {Signal-action-flags}\"");
      return -1;
    };

    if(ktfrealq(ktastk(isp0 + 3))){
      signo = rtastk(isp0 + 3);
      sigaddset(&next_sa.sa_mask, signo);
    }
    else if(ktastk(isp0 + 3) == ktfref){
    }
    else if(ktflistq(ktastk(isp0 + 3))){
      ia = ktfaddr(ktastk(isp0 + 3));
      if(klist(ia) == ktfoper + mtflist
         && ktfreallistq(ia)) {
        for(i = 1; i <= ilist(2, ia - 1); i++) {
          signo = rlist(ia + i);
          sigaddset(&next_sa.sa_mask, signo);
        }
      } else {
        isp = isp0;
        *irtc = itfmessage(9, "General::wrongtype",
			   "\"Mask -> {Signal-numbers} or AllSignals\"");
        return -1;
      }
    }
    else if(ktfsymbolq(ktastk(isp0 + 3))){
      if(tfsamesymbolqk_(&ktastk(isp0 + 3), &itfAll))
        sigfillset(&next_sa.sa_mask);
      else {
        isp = isp0;
        *irtc = itfmessage(9, "General::wrongtype",
			   "\"Mask -> {Signal-numbers} or AllSignals\"");
        return -1;
      }
    }
    else{
      isp = isp0;
      *irtc = itfmessage(9, "General::wrongtype",
                         "\"Mask -> {Signal-numbers} or AllSignals\"");
      return -1;
        
    };
  } else
    ispopt = isp + 1;
  isp = isp0;
  if(*isp1 + 2 != ispopt) {
    *irtc = itfmessage(9, "General::narg", "\"1 + option\"");
    return -1;
  }

  if(ktfnonrealq(ktastk(*isp1 + 1))) {
    *irtc = itfmessage(9, "General::wrongtype",
		       "\"Signal-number for #1\"");
    return -1;
  }
  signo = rtastk(*isp1 + 1);
  if(signo < 0 || signo > SAD_MAXSIG) {
    *irtc = itfmessage(9, "General::wrongtype",
		       "\"Signal-number for #1\"");
    return -1;
  }

  if(next_sa.sa_handler == SIG_DFL)
    switch(signo) {
#if SAD_TRAP_SIGHUP
    case SIGHUP:
      next_sa.sa_sigaction = (___sa_sigaction_t)sad_SIGHUP_handler;
      sad_sa_flags |= SAD_SA_SIGINFO | SAD_SA_RESTART;
      break;
#endif

#if SAD_TRAP_SIGFPE
    case SIGFPE:
      next_sa.sa_sigaction = (___sa_sigaction_t)sad_SIGFPE_handler;
      sad_sa_flags |= SAD_SA_SIGINFO | SAD_SA_RESTART;
      break;
#endif

#if SAD_TRAP_SIGABRT
    case SIGABRT:
      next_sa.sa_sigaction = (___sa_sigaction_t)sad_SIGABRT_handler;
      sad_sa_flags |= SAD_SA_SIGINFO | SAD_SA_RESTART;
      break;
#endif

    default:
      break;
    }

  /* Translate SAD sa_flags -> sigaction.sa_flags */
  next_sa.sa_flags = 0;
#ifdef	SA_ONSTACK
  if(sad_sa_flags & SAD_SA_ONSTACK)	next_sa.sa_flags |= SA_ONSTACK;
#endif
  if(sad_sa_flags & SAD_SA_RESTART)	next_sa.sa_flags |= SA_RESTART;
  if(sad_sa_flags & SAD_SA_RESETHAND)	next_sa.sa_flags |= SA_RESETHAND;
  if(sad_sa_flags & SAD_SA_NOCLDSTOP)	next_sa.sa_flags |= SA_NOCLDSTOP;
  if(sad_sa_flags & SAD_SA_NODEFER)	next_sa.sa_flags |= SA_NODEFER;
#ifdef	SA_NOCLDWAIT
  if(sad_sa_flags & SAD_SA_NOCLDWAIT)	next_sa.sa_flags |= SA_NOCLDWAIT;
#endif
  if(sad_sa_flags & SAD_SA_SIGINFO)	next_sa.sa_flags |= SA_SIGINFO;


  sigemptyset(&full_mask);
  sigaddset(&full_mask, signo);
  sigprocmask(SIG_BLOCK, &full_mask, &old_mask);
  /* Enter critical section of sigactno */

  if(sigaction(signo, new_sa, &old_sa) != 0) {
    sigprocmask(SIG_SETMASK, &old_mask, NULL);
    /* Leave critical section of sigactno */
    *irtc = itfsyserr(9);
    return 1;
  }

  old_sigactno = sad_signal_actions[signo - 1];
  if(install_sigactno) {
    sad_signal_actions[signo - 1] = new_sigactno;
    sad_signal_flags[signo - 1] = next_sa.sa_flags;
  }

  sigprocmask(SIG_SETMASK, &old_mask, NULL);
  /* Leave critical section of sigactno */

  isp0 = isp;
  if((old_sa.sa_flags & SA_SIGINFO) != 0) {	/* sa_sigaction */
    if(       old_sa.sa_sigaction == (___sa_sigaction_t)SIG_DFL) { /* SIG_DFL */
      tfmakerulestk(ktfsymbol+itfAction, ktfsymbol+itfDefault);
    } else if(old_sa.sa_sigaction == (___sa_sigaction_t)SIG_IGN) { /* SIG_IGN */
      tfmakerulestk(ktfsymbol+itfAction, ktfsymbol + itfIgnore);
#if SAD_TRAP_SIGHUP
    } else if(signo == SIGHUP &&			/* SAD's HUP handler */
	      old_sa.sa_sigaction == (___sa_sigaction_t)sad_SIGHUP_handler) {
      tfmakerulestk(ktfsymbol+itfAction, ktfsymbol + itfDefault);
#endif
#if SAD_TRAP_SIGFPE
    } else if(signo == SIGFPE &&			/* SAD's FPE handler */
	      old_sa.sa_sigaction == (___sa_sigaction_t)sad_SIGFPE_handler) {
      tfmakerulestk(ktfsymbol+itfAction, ktfsymbol + itfDefault);
#endif
#if SAD_TRAP_SIGABRT
    } else if(signo == SIGABRT &&			/* SAD's ABRT handler */
	      old_sa.sa_sigaction == (___sa_sigaction_t)sad_SIGABRT_handler) {
      tfmakerulestk(ktfsymbol+itfAction, ktfsymbol + itfDefault);
#endif
    } else if(old_sa.sa_sigaction == (___sa_sigaction_t)sad_SIGACT_handler) {
      integer8 ki, iai;
      integer4 irtci;

      buf_printf(&buf, buf0, sizeof(buf0), SIG_ACT_CMD, old_sigactno);
      if(buf != NULL) {
	tfevals(buf, &ki, &irtci);
        iai = ktfaddr(ki);
	if(irtci == 0 && ktflistq(ki)
	   && klist(iai) == ktfoper + mtfruledelayed) {
	  isp += 1;
          ktastk(isp) = ktfcopy_(&ki);
          if(new_sa != NULL && (!install_sigactno || new_sigactno != old_sigactno)) {
            sigprocmask(SIG_BLOCK, &full_mask, &old_mask);
            /* Enter critical section of sigactno */
            if(sad_signal_handler_depth == 0) {
              buf[SIG_ACT_CMD_OFFSET] = 'R';
              tfevals(buf, &ki, &irtci);
              tfevals("SigAction$Clean[]", &ki, &irtci);
            } else {
              buf[SIG_ACT_CMD_OFFSET] = 'M';
              tfevals(buf, &ki, &irtci);
            }
            sigprocmask(SIG_SETMASK, &old_mask, NULL);
            /* Leave critical section of sigactno */
          }
	} else
	  tfmakerulestk(ktfsymbol+itfAction, old_sigactno);
	if(buf != buf0) free(buf);
      } else
	tfmakerulestk(ktfsymbol+itfAction, old_sigactno);
    } else {						/* Otherwise */
      buf_printf(&buf, buf0, sizeof(buf0), /* Violate ISO C standard */
		 "%p", *(void **)&old_sa.sa_sigaction);
      if(buf == NULL) {
	*irtc = itfsyserr(9);
	return 1;
      }
      ia = ktsalocb(0, buf);
      tfmakerulestk(ktfsymbol+itfAction, ktfstring + ia);
      if(buf != buf0) free(buf);
    }
  } else {					/* ANSI C/BSD style */
    if(       old_sa.sa_handler == SIG_DFL) {		/* SIG_DFL */
      tfmakerulestk(ktfsymbol+itfAction, ktfsymbol + itfDefault);
    } else if(old_sa.sa_handler == SIG_IGN) {		/* SIG_IGN */
      tfmakerulestk(ktfsymbol+itfAction, ktfsymbol + itfIgnore);
    } else {						/* Otherwise */
      buf_printf(&buf, buf0, sizeof(buf0), /* Violate ISO C standard */
		 "%p", *(void **)&old_sa.sa_handler);
      if(buf == NULL) {
	*irtc = itfsyserr(9);
	return 1;
      }
      ia = ktsalocb(0, buf);
      tfmakerulestk(ktfsymbol+itfAction, ktfstring + ia);
      if(buf != buf0) free(buf);
    }
  }

  isp2 = isp;
#ifdef	SA_ONSTACK
  if(old_sa.sa_flags & SA_ONSTACK) {
    isp += 1;
    rtastk(isp) = SAD_SA_ONSTACK;
  }
#endif
  if(old_sa.sa_flags & SA_RESTART) {
    isp += 1;
    rtastk(isp) = SAD_SA_RESTART;
  }
  if(old_sa.sa_flags & SA_RESETHAND) {
    isp += 1;
    rtastk(isp) = SAD_SA_RESETHAND;
  }
  if(old_sa.sa_flags & SA_NOCLDSTOP) {
    isp += 1;
    rtastk(isp) = SAD_SA_NOCLDSTOP;
  }
  if(old_sa.sa_flags & SA_NODEFER) {
    isp += 1;
    rtastk(isp) = SAD_SA_NODEFER;
  }
#ifdef	SA_NOCLDWAIT
  if(old_sa.sa_flags & SA_NOCLDWAIT) {
    isp += 1;
    rtastk(isp) = SAD_SA_NOCLDWAIT;
  }
#endif
  if(old_sa.sa_flags & SA_SIGINFO) {
    isp += 1;
    rtastk(isp) = SAD_SA_SIGINFO;
  }
  ia = ktfmakelist(isp2);
  isp = isp2;
  tfmakerulestk(ktfsymbol+itfFlag, ktflist + ia);

  isp2 = isp;
  for(signo = 1; signo <= SAD_MAXSIG && (status = sigismember(&old_sa.sa_mask, signo)) >= 0; signo++)
    if(status != 0) {
      isp += 1;
      rtastk(isp) = signo;
    }
  ia = ktfmakelist(isp2);
  isp = isp2;
  tfmakerulestk(ktfsymbol+itfMask, ktflist + ia);

  *kx = ktflist + ktfmakelist(isp0);
  isp = isp0;
  *irtc = 0;
  return 0;
}

static int Environments(integer4 *isp1, integer8 *kx,
			integer4 *irtc) {
  extern char **environ;
  integer4 isp0;
  char **ep;

  if(isp != *isp1 + 1 || ktastk(isp) != ktfoper + mtfnull) {
    *irtc = itfmessage(9, "General::narg", "\"0\"");
    return -1;
  }

  isp0 = isp;
  for(ep = environ; *ep; ep++) {
    isp += 1;
    ktastk(isp) = ktfstring + ktsalocb(-1, *ep);
  }

  if(isp == isp0) {
    isp += 1;
    ktastk(isp) = ktfoper + mtfnull;
  }

  *kx = ktflist + ktfmakelist(isp0);
  isp = isp0;
  *irtc = 0;
  return 0;
}

static int GetEnv(integer4 *isp1, integer8 *kx, 
		  integer4 *irtc) {
  integer8 ka;
  char *buf;

  if(isp != *isp1 + 1) {
    *irtc = itfmessage(9, "General::narg", "\"1\"");
    return -1;
  }

  if((ktfmask & ktastk(*isp1 + 1)) != ktfstring) {
    *irtc = itfmessage(9, "General::wrongtype", "\"Character-string\"");
    return -1;
  }

  ka = (ktamask & ktastk(*isp1 + 1));
#if SAD_REQUIRE_STRING_TERMINATION
  jlist(ilist(1, ka) + 1, ka + 1) = '\0';
#endif

  if(ilist(1, ka) == 0) {
    *irtc = itfmessage(9, "General::wrongtype", "\"Character-string\"");
    return -1;
  }

  buf = getenv(&jlist(1, ka + 1));

  if(buf != NULL)
    *kx = ktfstring + ktsalocb(-1, buf);
  else
    *kx = ktfstring + ktsalocb(-1, "");
  *irtc = 0;
  return 0;
}

static int SetEnv(integer4 *isp1, integer8 *kx,
		  integer4 *irtc) {
  integer8 ia1, ia2;
#if !HAVE_BSD_SETENV
  char *buf;
#endif

  if(isp != *isp1 + 2) {
    *irtc = itfmessage(9, "General::narg", "\"2\"");
    return -1;
  }

  if((ktfmask & ktastk(*isp1 + 1)) != ktfstring) {
    *irtc = itfmessage(9, "General::wrongtype",
		       "\"Character-string for #1\"");
    return -1;
  }

  if((ktfmask & ktastk(*isp1 + 2)) != ktfstring) {
    *irtc = itfmessage(9, "General::wrongtype",
		       "\"Character-string for #2\"");
    return -1;
  }

  ia1 = (ktamask & ktastk(*isp1 + 1));
#if SAD_REQUIRE_STRING_TERMINATION
  jlist(ilist(1, ia1) + 1, ia1 + 1) = '\0';
#endif

  if(ilist(1, ia1) == 0) {
    *irtc = itfmessage(9, "General::wrongtype",
		       "\"Character-string for #1\"");
    return -1;
  }

  ia2 = (ktamask & ktastk(*isp1 + 2));
#if SAD_REQUIRE_STRING_TERMINATION
  jlist(ilist(1, ia2) + 1, ia2 + 1) = '\0';
#endif

#if HAVE_BSD_SETENV
#ifdef USE_NULL_STRING_RESCUE_HACK
  if((ilist(1, ia2) != 0
      ? setenv(&jlist(1, ia1 + 1), &jlist(1, ia2 + 1), 1)
      : setenv(&jlist(1, ia1 + 1), "", 1)) < 0) {
    *irtc = itfsyserr(9);
    return 1;
  }
#else
  if(setenv(&jlist(1, ia1 + 1), &jlist(1, ia2 + 1), 1) < 0) {
    *irtc = itfsyserr(9);
    return 1;
  }
#endif
#else
  buf = malloc(ilist(1, ia1) + ilist(1, ia2) + 2);
  if(buf == NULL) {
    *irtc = itfsyserr(9);
    return 1;
  }
  memcpy(buf, &jlist(1, ia1 + 1), ilist(1, ia1));
  buf[ilist(1, ia1)] = '=';
  memcpy(buf + ilist(1, ia1) + 1, &jlist(1, ia2 + 1), ilist(1, ia2));
  buf[ilist(1, ia1) + ilist(1, ia2) + 1] = '\0';
  if(putenv(buf) < 0) {
    free(buf);
    *irtc = itfsyserr(9);
    return 1;
  }
  /*
   * CAUTION:
   *  Don't free buf pointer if putenv is succeeded,
   *  because putenv links buf pointer into environment.(SUSv2)
   */
#endif

  *kx = ktfoper + mtfnull;
  *irtc = 0;
  return 0;
}

static int UnsetEnv(integer4 *isp1, integer8 *kx,
		    integer4 *irtc) {
  integer8 ia;

  if(isp != *isp1 + 1) {
    *irtc = itfmessage(9, "General::narg", "\"1\"");
    return -1;
  }

  if(ktfnonstringq(ktastk( *isp1 + 1))) {
    *irtc = itfmessage(9, "General::wrongtype", "\"Character-string\"");
    return -1;
  }

  ia = ktfaddr(ktastk(*isp1 + 1));
#if SAD_REQUIRE_STRING_TERMINATION
  jlist(ilist(1, ia) + 1, ia + 1) = '\0';
#endif

  if(ilist(1, ia) == 0) {
    *irtc = itfmessage(9, "General::wrongtype", "\"Character-string\"");
    return -1;
  }

  unsetenv(&jlist(1, ia + 1));

  *kx = ktfoper + mtfnull;
  *irtc = 0;
  return 0;
}

static int GetDirectory(integer4 *isp1, integer8 *kx, 
			integer4 *irtc) {
  char *buf;

  if(isp != *isp1 + 1
     || ktastk(isp) != ktfoper + mtfnull) {
    *irtc = itfmessage(9, "General::narg", "\"0\"");
    return -1;
  }

  buf = getcwd(NULL, MAXPATHLEN);
  if(buf == NULL) {
    *irtc = itfsyserr(9);
    return 1;
  }

  *kx = ktfstring + ktsalocb(-1, buf); free(buf);
  *irtc = 0;
  return 0;
}

static int SetDirectory(integer4 *isp1, integer8 *kx,
			integer4 *irtc) {
  integer8 ia;
  char *expanded, *buf;

  if(isp != *isp1 + 1) {
    *irtc = itfmessage(9, "General::narg", "\"1\"");
    return -1;
  }

  if(ktfnonstringq(ktastk(*isp1 + 1))){
    *irtc = itfmessage(9, "General::wrongtype", "\"Character-string\"");
    return -1;
  }

  ia = ktfaddr(ktastk(*isp1 + 1));
#if SAD_REQUIRE_STRING_TERMINATION
  jlist(ilist(1, ia) + 1, ia + 1) = '\0';
#endif

  expanded = expand_tilde(&jlist(1, ia + 1));
  if(expanded == NULL) {
    switch(errno) {
    case ENOMEM:
      *irtc = itfsyserr(9);
      return 1;
      break;

    default:
      *irtc = itfmessage(9, "System::error",
			 "\"No such home directory\"");
    }
    return -1;
  }

  if(chdir(expanded) < 0) {
    free(expanded);
    *irtc = itfsyserr(9);
    return 1;
  }
  free(expanded);

  /* Read-back currect directory */
  buf = getcwd(NULL, MAXPATHLEN);
  if(buf == NULL) {
    *irtc = itfsyserr(9);
    return 1;
  }

  *kx = ktfstring + ktsalocb(-1, buf); free(buf);
  *irtc = 0;
  return 0;
}

static int Execve(integer4 *isp1, integer8 *kx,
		  integer4 *irtc) {
  static const char *options[] = {"ForceCloseOnExec", "CloseStd", NULL};
  integer8 ia, ia1, ia2, ia3;
  integer4 isp0, ispopt;
  int i, j;
  char **argv, **envp;
  bit_field_t *bits = NULL;
  int max_open_files, force_close_on_exec, close_std;

  max_open_files = getdtablesize();

  force_close_on_exec = False;
  close_std = False;
  isp0 = isp;
  ispopt = itfgetoptionstk(*isp1, options);
  if(ispopt > 0) {
    if(ktfrealq(ktastk(isp0 + 1))){
      force_close_on_exec = (rtastk(isp0 + 1) != 0.0);
    }
    else if(ktastk(isp0 + 1) == ktfref){
    }
    else {
      isp = isp0;
      *irtc = itfmessage(9, "General::wrongtype",
			 "\"ForceCloseOnExec -> True of False\"");
      return -1;
    };

    if(ktfrealq(ktastk(isp0 + 2))){
      close_std = (rtastk(isp0 + 2) != 0.0);
    }
    else if(ktastk(isp0 + 1) == ktfref){
    }
    else {
      isp = isp0;
      *irtc = itfmessage(9, "General::wrongtype",
			 "\"ForceCloseOnExec -> True of False\"");
      return -1;
    }

  } else
    ispopt = isp + 1;

  isp = isp0;
  if(*isp1 + 4 != ispopt) {
    *irtc = itfmessage(9, "General::narg", "\"3 (+ option)\"");
    return -1;
  }

  if(ktfnonstringq(ktastk( *isp1 + 1))) {
    *irtc = itfmessage(9, "General::wrongtype",
		       "\"Character-string for #1\"");
    return -1;
  }

  if(ktfnonlistq(ktastk( *isp1 + 2))) {
    *irtc = itfmessage(9, "General::wrongtype",
		       "\"List of Character-strings for #2\"");
    return -1;
  }

  if(ktfnonlistq(ktastk( *isp1 + 3))) {
    *irtc = itfmessage(9, "General::wrongtype",
		       "\"List of Character-strings for #3\"");
    return -1;
  }

  ia1 = ktfaddr(ktastk(*isp1 + 1));
#if SAD_REQUIRE_STRING_TERMINATION
  jlist(ilist(1, ia1) + 1, ia1 + 1) = '\0';
#endif

  ia2 = ktfaddr(ktastk(*isp1 + 2));
  if(itastk(2, ia2 - 1) < 1) {
    *irtc = itfmessage(9, "General::wrongtype",
		       "\"List of Character-strings for #2\"");
    return -1;
  }
  for(i = 0, j = 1; i < itastk(2, ia2 - 1); i++, j++)
    if(ktfstringq(klist(ia2 + j))) {
      ia = ilist(2, ia2 + j);
#if SAD_REQUIRE_STRING_TERMINATION
      jlist(ilist(1, ia) + 1, ia + 1) = '\0';
#endif
    } else {
      *irtc = itfmessage(9, "General::wrongtype",
			 "\"List of Character-strings for #2\"");
      return -1;
    }

  ia3 = ktfaddr(ktastk(*isp1 + 3));
  for(i = 0, j = 1; i < itastk(2, ia3 - 1); i++, j++)
    if(ktfstringq(klist(ia3 + j))) {
      ia = ilist(2, ia3 + j);
#if SAD_REQUIRE_STRING_TERMINATION
      jlist(ilist(1, ia) + 1, ia + 1) = '\0';
#endif
    } else {
      *irtc = itfmessage(9, "General::wrongtype",
			 "\"List of Character-strings for #3\"");
      return -1;
    }

  /* Scan file descriptor and mark !close-on-exec descriptor */
  if(force_close_on_exec) {
    while(max_open_files > 0) {
      j = fcntl(max_open_files - 1, F_GETFD, 1);
      if(j != -1 && (j & FD_CLOEXEC) == 0) break;
      max_open_files -= 1;
    }

    if(max_open_files > 0) {
      if(!BF_ALLOC(max_open_files, bits)) {
	*irtc = itfsyserr(9);
	return 1;
      }
      BF_ZERO(max_open_files, bits);
    }
  }

  /* Construct arguments for execve(2) */
  envp = malloc(sizeof(char *) * (ilist(2, ia3 - 1) + 1)
		+ sizeof(char *) * (ilist(2, ia2 - 1) + 2));
  if(envp != NULL) {
    argv = envp + (ilist(2, ia3 - 1) + 1);

    for(i = 0, j = 1; j <= ilist(2, ia3 - 1); i++, j++)
      if(ktfstringq(klist(ia3 + j))) {
	ia = ilist(2, ia3 + j);
	envp[i] = &jlist(1, ia + 1);
      }
    envp[i] = NULL;

    argv[0] = "sh";
    for(i = 1, j = 1; j <= itastk(2, ia2 - 1); i++, j++)
      if(ktfstringq(klist(ia2 + j))) {
	ia = ilist(2, ia2 + j);
	argv[i] = &jlist(1, ia + 1);
      }
    argv[i] = NULL;

    if(force_close_on_exec) /* Set close-on-exec flag if needed */
      for(i = close_std ? 0 : 3; i < max_open_files; i++) {
	j = fcntl(i, F_GETFD, 1);
	if(j != -1 && (j & FD_CLOEXEC) == 0) {
	  BF_SET(i, bits);
	  fcntl(i, F_SETFD, FD_CLOEXEC);
	}
      }

    execve(&jlist(1, ia1 + 1), argv + 1, envp);
    if(errno == ENOEXEC) {
      argv[1] = &jlist(1, ia1 + 1);
      execve(SAD_BSHELL_PATH, argv, envp);
    }
    free(envp);

    if(force_close_on_exec) /* Recover close-on-exec flag if need */
      for(i = close_std ? 0 : 3; i < max_open_files; i++)
	if(BF_ISSET(i, bits))
	  fcntl(i, F_SETFD, 0);
    BF_FREE(bits);
  }
  *irtc = itfsyserr(9);
  return 0;
}

static int System(integer4 *isp1, integer8 *kx,
		  integer4 *irtc) {
  integer8 ka;
  int status;

  if(isp != *isp1 + 1) {
    *irtc = itfmessage(9, "General::narg", "\"1\"");
    return -1;
  }

  if((ktfmask & ktastk(*isp1 + 1)) != ktfstring) {
    *irtc = itfmessage(9, "General::wrongtype", "\"Character-string\"");
    return -1;
  }

  ka = (ktamask & ktastk(*isp1 + 1));
#if SAD_REQUIRE_STRING_TERMINATION
  jlist(ilist(1, ka) + 1, ka + 1) = '\0';
#endif

  status = system(&jlist(1, ka + 1));
  if(status < 0) {
    *kx = 0;
    *irtc = itfsyserr(9);
    return 1;
  }

  {/* *kx = kfromr(r_true); */
	  real8 vx=status;
	  *kx=kfromr(vx);
	  *irtc = 0;
	  return 0;
  }
}

/* Process Group ID family */
static int GetPGID(integer4 *isp1, integer8 *kx,
		   integer4 *irtc) {
  real8 vx;
  pid_t  pid = 0;
  pid_t pgid;

  if(isp != *isp1 + 1) {
    *irtc = itfmessage(9, "General::narg", "\"0 or 1\"");
    return -1;
  }

  if(ktfrealq(ktastk(*isp1 + 1))){
    pid = rtastk( *isp1 + 1);
  }
  else if((ktfmask & ktastk(*isp1 + 1)) == ktfoper){
    if(!(ktamask & ktastk(*isp1 + 1)) == mtfnull && *isp1 + 1 == isp) {
      *irtc = itfmessage(9, "General::wrongtype",
			 "\"Null or PID for #1\"");
      return -1;
    }
    pid = 0;	/* Current process */
  }
  else{
    *irtc = itfmessage(9, "General::wrongtype",
		       "\"Null or PID for #1\"");
    return -1;
  }

  if((pgid = getpgid(pid)) == -1) {
    *irtc = itfsyserr(9);
    return 1;
  }

  vx = pgid;
  *kx = kfromr(vx);
  *irtc = 0;

  return 0;
}

static int SetPGID(integer4 *isp1, integer8 *kx,
		   integer4 *irtc) {
  pid_t  pid;
  pid_t pgid;

  if(isp != *isp1 + 2) {
    *irtc = itfmessage(9, "General::narg", "\"2\"");
    return -1;
  }

  if((ktrmask & ktastk(*isp1 + 1)) == ktfnr) {
    *irtc = itfmessage(9, "General::wrongtype", "\"Process ID for #1\"");
    return -1;
  }
  pid = rtastk(*isp1 + 1);

  if((ktrmask & ktastk(*isp1 + 2)) == ktfnr) {
    *irtc = itfmessage(9, "General::wrongtype", "\"Process Group ID for #2\"");
    return -1;
  }
  pgid = rtastk(*isp1 + 2);

  if(setpgid(pid, pgid) != 0) {
    *irtc = itfsyserr(9);
    return 1;
  }

  *kx = ktfoper + mtfnull;
  *irtc = 0;
  return 0;
}

/* Process ID get family template */
#define DefineGetIDs(type, method) \
static int Get##type(integer4 *isp1, integer8 *kx, \
		     integer4 *irtc) {			\
  real8 vx; \
  if(isp != *isp1 + 1 || ktastk(isp) != ktfoper + mtfnull) { \
    *irtc = itfmessage(9, "General::narg", "\"0\""); \
    return -1; \
  } \
  vx = method(); \
  *kx = kfromr(vx); \
  *irtc = 0; \
  return 0; \
}

DefineGetIDs(PID,  getpid)
DefineGetIDs(PPID, getppid)
DefineGetIDs(UID,  getuid)
DefineGetIDs(EUID, geteuid)
DefineGetIDs(GID,  getgid)
DefineGetIDs(EGID, getegid)
/* CAUTION: Don't append `;' after DefineGetIDs() macro		*/
/* ISO C does not allow extra `;' outside of a function!	*/

/* SADScript function registration of Process system call stuff */
#define REG	dlfunaloc
#define REG8	dlfunaloc8
int sadDefFunc_Process(void) {
  init_signal_handler();

  /* Process control */
  REG8("Fork",		Fork,		 0, NULL, NULL, 0);
  REG8("Wait",		Wait,		 0, NULL, NULL, 0);
  REG8("Wait4",		Wait4,		 1, NULL, NULL, 0);
  REG8("Kill",		Kill,		 1, NULL, NULL, 0);
  REG8("Sleep",		Sleep,		 1, NULL, NULL, 0);
  REG8("Pause",		Sleep,		 1, NULL, NULL, 0);

  /* Signal control */
  REG8("SigPending",	SigPending,	0, NULL, NULL, 0);
  REG8("SigSuspend",	SigSuspend,	1, NULL, NULL, 0);
  REG8("SigProcMask",	SigProcMask,	2, NULL, NULL, 0);
  REG8("SigAction$",	SigAction,	1, NULL, NULL, 0);

  /* Environment control */
  REG8("Environments",	Environments,	0, NULL, NULL, 0);
  REG8("GetEnv",	GetEnv,		1, NULL, NULL, 0);
  REG8("SetEnv",		SetEnv,		2, NULL, NULL, 0);
  REG8("UnsetEnv",	UnsetEnv,	1, NULL, NULL, 0);
  REG8("GetDirectory",	GetDirectory,	0, NULL, NULL, 0);
  REG8("SetDirectory",	SetDirectory,	1, NULL, NULL, 0);

  /* Exec family */
  REG8("Execve",		Execve,		3, NULL, NULL, 0);
  REG8("System",	System,		1, NULL, NULL, 0);

  /* Process ID family */
  REG8("GetPID",		GetPID,		0, NULL, NULL, 0);
  REG8("GetPPID",	GetPPID,	0, NULL, NULL, 0);
  REG8("GetUID",		GetUID,		0, NULL, NULL, 0);
  REG8("GetEUID",	GetEUID,	0, NULL, NULL, 0);
  REG8("GetGID",		GetGID,		0, NULL, NULL, 0);
  REG8("GetEGID",	GetEGID,	0, NULL, NULL, 0);

  /* Process Group ID family */
  REG8("GetPGID",	GetPGID,	1, NULL, NULL, 0);
  REG8("SetPGID",	SetPGID,	2, NULL, NULL, 0);

  return 0;
}

/* End of File */
