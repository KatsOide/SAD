#include <sim/sad_f2c.h>

#include <sys/resource.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <errno.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <omp.h>

/* Fortran extensions for unix */

/* CALL IDate(TARRAY)
   Fills TARRAY with the numerical values at the current local time of
   day, month (in the range 1-12), and year in elements 1, 2, and 3,
   respectively.  The year has four significant digits. */
void idate_(integer4 tarray[3]) {
  time_t now;
  struct tm local;

  time(&now);
  localtime_r(&now, &local);

  tarray[0] = local.tm_mday;
  tarray[1] = local.tm_mon;
  tarray[2] = local.tm_year;
}

/* Time()
   Returns the current time encoded as an integer (in the manner of the
   UNIX function `time(3)').  This value is suitable for passing to
   `CTIME', `GMTIME', and `LTIME'. */
integer4 time_() {
  return time(NULL);
}

/* GETPID		Process ID function
 * Syntax:		RESULT = GETPID()
 * Retuen:		INTEGER of default kind
 */
integer getpid_(void) {
  return getpid();
}

/* FORK_WORKER		Create a new worker process with low-priority
 * Syntax:		PID = FORK()
 * Return:		INTEGER of default kind
 */
integer fork_worker_(void) {
  const int id = fork();
  if(id == 0){
    setpriority(PRIO_PROCESS, 0, 19);
    /* omp_set_num_threads(1); */
    tfsavesharedmap_();}
  return id;
}

/* FORK			Create a new process
 * Syntax:		PID = FORK()
 * Return:		INTEGER of default kind
 */
integer fork_(void) {
  return fork();
}

/* WAIT			Wait any child process
 * Syntax:		PID = WAIT(STATUS)
 * Argument:
 *   STATUS 		Shall be of default INTEGER type
 * Return:		INTEGER of default kind
 */
integer wait_(integer *status_) {
  pid_t pid;
  int status;

  pid = wait(&status);

  *status_ = status;
  return pid;
}

/* WAITPID		Wait givend PID process
 * Syntax:		PID = WAITPID(WPID, STATUS)
 * Argument:
 *   WPID		Shall be of default INTEGER type
 *   STATUS 		Shall be of default INTEGER type
 * Return:		INTEGER of default kind
 */
integer waitpid_(integer *wpid_, integer *status_) {
  pid_t pid;
  int status;

  pid = waitpid(*wpid_, &status, 0);

  *status_ = status;
  return pid;
}

/* WAITPID_NOHANG	Wait givend PID process without blocking
 * Syntax:		PID = WAITPID(WPID, STATUS)
 * Argument:
 *   WPID		Shall be of default INTEGER type
 *   STATUS 		Shall be of default INTEGER type
 * Return:		INTEGER of default kind
 */
integer waitpid_nohang_(integer *wpid_, integer *status_) {
  pid_t pid;
  int status;

  pid = waitpid(*wpid_, &status, WNOHANG);

  *status_ = status;
  return pid;
}

/* EXIT_WITHOUT_HOOKS	Terminate process without atexit hooks
 * Syntax:		CALL EXIT_WITHOUT_HOOKS(STATUS)
 * Argument:
 *   STATUS 		Shall be of default INTEGER type
 * Return:		None
 */
void exit_without_hooks_(integer *status_) {
  _exit(*status_);
}

/* SYSTEM		Execute a shell command
 * Syntax:		STATUS = SYSTEM(COMMAND)
 * Argument:
 *   COMMAND		Shall be of default CHARACTER type
 * Return:		Exit status of child process
 */
integer system_(const_character cmd, ftnlen clen) {
  size_t length = len_trim(cmd, clen);
  int status;

  if(length < 1) {
    return system(NULL);
  }

  if(cmd[length - 1] == '\0') {	/* null terminated command string case */
    status = system(cmd);
    return status;
  } else {
    char *buffer = malloc(length + 1);

    if(buffer == NULL) return 127;

    memcpy(buffer, cmd, length); buffer[length] = '\0';
    status = system(buffer); free(buffer);
    return status;
  }
}

/* PERROR		Print system error message
 * Syntax:		CALL PERROR(MESSAGE)
 * Argument:
 *   MESSAGE		Shall be of default CHARACTER type
 */
void perror_(const_character msg, ftnlen mlen) {
  size_t length = len_trim(msg, mlen);

  if(length < 1) {
    perror("");
    return;
  }

  if(msg[length - 1] == '\0') {	/* null terminated message string case */
    perror(msg);
  } else {
    char buffer[length + 1];

    memcpy(buffer, msg, length); buffer[length] = '\0';
    perror(buffer);
  }
}

/* GERROR		Get last system error message
 * Syntax:		CALL GERROR(RESULT)
 * Argument:
 *   RESULT		Shall be of default CHARACTER type
 */
void gerror_(character buff, ftnlen blen) {
  size_t length;
  char *str;

  str = strerror(errno); length = strlen(str);
  if(length > blen) length = blen;
  memcpy(buff, str, length);
  if(blen > length) memset(buff + length, ' ', blen - length);
}

/* End of File */
