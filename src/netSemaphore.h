
/* port number of the server */

#define NETSEM_SERVER_PORT 17006



/* maximum number of semaphores */

#define MAX_NUM_SEM 256



/* message length */

#define HDR_LEN 6
#define MSG_LEN 256



/* operation codes */

#define TAKE 1
#define GIVE 2
#define INFO 3



/* information selector */

#define OWNER  0
#define WAITER 1



/* message structure */

#define OPCODE 0
#define SEMAID 1
#define TMO_HI 2
#define TMO_LO 3
#define PID_HI 4
#define PID_LO 5
#define INFO_0 HDR_LEN

/* netSemaphore API prototype */

extern int netSemTake( char * server, int, int );
extern int netSemGive( char * server, int );
extern int netSemInfo( char * server, int, int, char * );

/* End of File */
