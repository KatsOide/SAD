/*
 *	Author:	J.Odagiri
 *		KEKB acclerator control group
 *		0298-64-5315 
 *	Date:	7-99
 */

#include <sys/types.h>
#include <unistd.h>
#include <string.h>
#include <ctype.h>
/*for SAD/fortran interface*/
#include <stdio.h>

#ifdef CYGWIN
int netSemTake(char *server, int semid, int tmo){
  return 0;
}

int netSemGive(char *server, int semid){
  return 0;
}

int netSemInfo(char *server, int semid,int sw,char *buf){
  return 0;
}

#else

#include <sys/socket.h>
#include <netinet/in.h>
#include <netdb.h>


#include "netSemaphore.h"

#define BYTE_MASK 0xff;
#define NETSEM_SERVER_ADDR "172.19.51.15"


int netsemtake(char *server, int *semid, int *tmo){
  return netSemTake(server, *semid, *tmo);
}

int netsemgive_(char *server, int *semid){
  return netSemGive(server, *semid);
}

int netseminfo_(char *server, int *semid,int *sw,char *buf){
  return netSemInfo(server, *semid,*sw, buf);
}

static int connect_to_server( char * );
static int take_semaphore( int, int );
static int give_semaphore( int );
static int get_sem_info( int, int, char * );
static int recv_msg( unsigned char * );
static int send_msg( unsigned char * );
static int check_msg( int, unsigned char *, unsigned char * );

static int cli_pid;
static int sock_fd;
static int connected;



/*
 *
 *	netSemTake()
 *
 */
int netSemTake( char *hostname, int semId, int timeout )
{
  if ( !connected ) {
    if ( (cli_pid = getpid()) < 0 ) return( -1 );
    if ( connect_to_server(hostname)  < 0 ) return( -1 );
    connected = 1;
  }

  if ( take_semaphore( semId, timeout ) < 0 ) {
    return( -1 );
  }

  return( 0 );
}



/*
 *
 *	netSemGive()
 *
 */
int netSemGive( char *hostname, int semId )
{
  if ( !connected ) {
    if ( (cli_pid = getpid()) < 0 ) return( -1 );
    if ( connect_to_server(hostname)  < 0 ) return( -1 );
    connected = 1;
  }

  if ( give_semaphore( semId ) < 0 ) {
    return( -1 );
  }

  return( 0 );
}



/*
 *
 *	netSemInfo()
 *
 */
int netSemInfo( char *hostname, int semId, int select, char *buffer )
{
  if ( !connected ) {
    if ( (cli_pid = getpid()) < 0 ) return( -1 );
    if ( connect_to_server(hostname)  < 0 ) return( -1 );
    connected = 1;
  }

  if ( get_sem_info( semId, select, buffer ) < 0 ) {
    return( -1 );
  }

  return( 0 );
}

/*
 *
 *	connect_to_server()
 *
 */
static int connect_to_server( char *NShost)
{
  struct sockaddr_in serverAddr;	/* server's address */
  struct hostent *phostent;

  /*
   * Open the socket. Use ARPA Internet address format and stream
   * sockets. Format described in <sys/socket.h>.
   */
  if ( ( sock_fd = socket( AF_INET,
		       SOCK_STREAM,
		       0 ) ) < 0 ) {
    perror( "netSemTake: socket creation error" );
    return( -1 );
  }

  /*
   * set the sever sock_addr structure
   */
  if (NShost == NULL || *NShost == '\0')
    phostent=gethostbyname(NETSEM_SERVER_ADDR);
  else{
    register unsigned char *pc;
    for (pc=(unsigned char*)NShost; *pc; pc++){
      *pc=tolower(*pc);
    }
    phostent=gethostbyname(NShost);
    if (phostent == NULL){
      fprintf(stderr,"cannot fined hostname : %s\n",NShost);
      return (-1);
    }
  }
  memset( (char *)&serverAddr, 0, sizeof( serverAddr ) );
  serverAddr.sin_family      = AF_INET;
  serverAddr.sin_port        = htons( NETSEM_SERVER_PORT );
  memcpy((char *) &serverAddr.sin_addr, phostent->h_addr_list[0], phostent->h_length);

  /*
   * connect to the network semaphore server
   */
  if ( connect( sock_fd,
		(struct sockaddr *)&serverAddr,
		sizeof(serverAddr)) < 0 ) {
    perror( "netSemTake: connect error" );
    close( sock_fd );
    return( -1 );
  }

  return( 0 );
}



/*
 *
 *	take_semaphore()
 *
 */
static int take_semaphore( int semId, int timeout )
{
  static unsigned char send_buf[MSG_LEN];
  static unsigned char recv_buf[MSG_LEN];

  /*
   * set the message structure
   */
  send_buf[OPCODE] = (unsigned char) TAKE;
  send_buf[SEMAID] = (unsigned char) semId;
  send_buf[TMO_HI] = (unsigned char) ( timeout >> 8 ) & BYTE_MASK; 
  send_buf[TMO_LO] = (unsigned char) timeout & BYTE_MASK;
  send_buf[PID_HI] = (unsigned char) ( cli_pid >> 8 ) & BYTE_MASK; 
  send_buf[PID_LO] = (unsigned char) cli_pid & BYTE_MASK;

  /*
   * send the message to the sever
   */
  if ( send_msg( send_buf ) < 0 ) {
    perror( "netSemTake: send_mgs failed" );
    return( -1 );
  }
  /*
   * recieve the response from the sever
   */
  if ( recv_msg( recv_buf ) < 0 ) {
    perror( "netSemTake: recv_msg failed" );
    return( -1 );
  }
  /*
   * check the response
   */
  if ( check_msg( TAKE, send_buf, recv_buf ) < 0 ) {
    return( -1 );
  }
  return( 0 );
}



/*
 *
 *	give_semaphore()
 *
 */
static int give_semaphore( int semId )
{
  static unsigned char send_buf[MSG_LEN];
  static unsigned char recv_buf[MSG_LEN];

  /*
   * set the message structure
   */
  send_buf[OPCODE] = (unsigned char) GIVE;
  send_buf[SEMAID] = (unsigned char) semId;
  send_buf[TMO_HI] = (unsigned char) 0; 
  send_buf[TMO_LO] = (unsigned char) 0;
  send_buf[PID_HI] = (unsigned char) ( cli_pid >> 8 ) & BYTE_MASK; 
  send_buf[PID_LO] = (unsigned char) cli_pid & BYTE_MASK;

  /*
   * send the message to the sever
   */
  if ( send_msg( send_buf ) < 0 ) {
    perror( "netSemGive: send_mgs failed" );
    return( -1 );
  }

  /*
   * recieve the response from the sever
   */
  if ( recv_msg( recv_buf ) < 0 ) {
    perror( "netSemGive: recv_msg failed" );
    return( -1 );
  }

  /*
   * check the response
   */
  if ( check_msg( GIVE, send_buf, recv_buf ) < 0 ) {
    return( -1 );
  }

  return( 0 );
}



/*
 *
 *	get_sem_info()
 *
 */
static int get_sem_info( int semId, int select, char *buffer )
{
  static unsigned char send_buf[MSG_LEN];
  static unsigned char recv_buf[MSG_LEN];

  /*
   * set the message structure
   */
  send_buf[OPCODE] = (unsigned char) INFO;
  send_buf[SEMAID] = (unsigned char) semId;
  send_buf[TMO_HI] = (unsigned char) select; 
  send_buf[TMO_LO] = (unsigned char) 0;
  send_buf[PID_HI] = (unsigned char) ( cli_pid >> 8 ) & BYTE_MASK; 
  send_buf[PID_LO] = (unsigned char) cli_pid & BYTE_MASK;

  /*
   * send the message to the sever
   */
  if ( send_msg( send_buf ) < 0 ) {
    perror( "netSemInfo: send_mgs failed" );
    return( -1 );
  }

  /*
   * recieve the response from the sever
   */
  if ( recv_msg( recv_buf ) < 0 ) {
    perror( "netSemInfo: recv_msg failed" );
    return( -1 );
  }

  /*
   * check the response
   */
  if ( check_msg( INFO, send_buf, recv_buf ) < 0 ) {
    return( -1 );
  }

  recv_buf[MSG_LEN - 1] = '\0';
  strncpy( buffer, (char *)&recv_buf[INFO_0], MSG_LEN - HDR_LEN );

  return( 0 );
}



/*
 *
 *	send_msg()
 *
 */
static int send_msg( unsigned char *buf )
{
  int nleft = MSG_LEN;
  int nwritten;

  while ( nleft ) {

    nwritten = send( sock_fd, buf, nleft, 0 );

    if ( nwritten <= 0 ) {
      perror( "send_msg: send error" );
      return( -1 );
    }

    nleft -= nwritten;
    buf   += nwritten;

  }

  return( 0 );
}



/*
 *
 *	recv_msg()
 *
 */
static int recv_msg( unsigned char *buf )
{
  int nleft = MSG_LEN;
  int nread;

  while ( nleft ) {

    nread = recv( sock_fd, buf, nleft, 0 );

    if ( nread <= 0 ) {
      perror( "recv_msg: recieve error" );
      return( -1 );
    }

    nleft -= nread;
    buf   += nread;

  }

  return( 0 );
}



/*
 *
 *	check_msg()
 *
 */
static int check_msg( int opCode,
		      unsigned char *send_buf,
		      unsigned char *recv_buf )
{
  int i;
  int err_flag = 0;

  if ( recv_buf[0] == 0 )
    {
      err_flag = 1;
      recv_buf[0] = opCode;
    }

  for ( i=0; i < HDR_LEN; i++ )
    {
      if ( send_buf[i] != recv_buf[i] )
	{
	  switch ( opCode )
	    {
	    case TAKE:
	      perror( "netSemTake: got inconsistent response" );
	      break;
	    case GIVE:
	      perror( "netSemGive: got inconsistent response" );
	      break;
	    case INFO:
	      perror( "netSemInfo: got inconsistent response" );
	      break;
	    default:
	      break;
	    }
	  return( -1 );
	}
    }

  if ( err_flag  ) return( -1 );

  return( 0 );
}


#endif
