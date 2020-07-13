#include <sim/sad_api.h>
#include <sim/TFCODE.h>
#include <sim/TFCBK.h>
#include <sim/TFSTK.h>

#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <netdb.h>
#include <sys/time.h>
#include <sys/select.h>

#include <sys/ioctl.h>
#ifndef FIONREAD
#  include <sys/filio.h>
#endif

#include <sys/param.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>

#define NTOP_BUFFER_LEN		64

/* SADScript function helper of Network stuff */
extern integer4 itopenbuf_(integer4*,integer4*);

static int tfSocketOpen(integer4 *isp1,
			integer8 *kx,
			integer4 *irtc, int socktype) {
  integer8 ia1, ia2, iax;
  integer4 ispi, in1, in2;
  real8 vx, timeout;
  int fd1, fd2, fdn, err;
  int port0 = 0, port = 0;
  int moder=1;
  char service[NI_MAXSERV], *src_addr, *src_addr0, *src_port;
  struct addrinfo hints, *res, *res0;
  struct sockaddr_storage addr;
  socklen_t na;

  if(isp < *isp1 + 1 || isp > *isp1 + 3) {
    *irtc = itfmessage(9, "General::narg", "\"1, 2, or 3\"");
    return -1;
  }

  ispi = *isp1;
  if(ktfrealq(ktastk(ispi + 1))){
    src_addr0 = NULL;
  }else if(ktfstringq(ktastk(ispi + 1))){
    ia1 = ktfaddr(ktastk(*isp1 + 1));
#if SAD_REQUIRE_STRING_TERMINATION
    jlist(ilist(1, ia1) + 1, ia1 + 1) = '\0';
#endif
    src_addr0 = &jlist(1, ia1 + 1);
    ispi = *isp1 + 1;
  }else{
    *irtc = itfmessage(9, "General::wrongtype",
		       "\"Network location or Port number for #1\"");
    return -1;
  }

  if(isp > ispi) {
    if(ktfrealq(ktastk(ispi + 1))){
      port = rtastk(ispi + 1);
      snprintf(service, sizeof(service), "%d", port);
      src_port = service;
    }
    else if(ktfstringq(ktastk(ispi + 1))){
      ia2 = ktfaddr(ktastk(*isp1 + 1));
#if SAD_REQUIRE_STRING_TERMINATION
      jlist(ilist(1, ia2) + 1, ia2 + 1) = '\0';
#endif
      src_port = &jlist(1, ia2 + 1);
    }else{
      if(src_addr0 == NULL)
	*irtc = itfmessage(9, "General::wrongtype",
			   "\"Port number for #1\"");
      else
	*irtc = itfmessage(9, "General::wrongtype",
			   "\"Port number for #2\"");
      return -1;
    }
  } else {
    if(src_addr0 == NULL)
      *irtc = itfmessage(9, "General::narg", "\"1 or 2\"");
    else
      *irtc = itfmessage(9, "General::narg", "\"2 or 3\"");
    return -1;
  }

  if(isp > ispi + 1) {
    if(ktfnonrealq(ktastk(ispi + 2))){
      if(src_addr0 == NULL)
	*irtc = itfmessage(9, "General::wrongtype",
			   "\"Timeout(seconds) for #2\"");
      else
	*irtc = itfmessage(9, "General::wrongtype",
			   "\"Timeout(seconds) for #3\"");
      return -1;
    }
    timeout = rtastk(ispi + 2);
  } else {
    timeout = 0.0;
  }

  src_addr = src_addr0;
  if(src_addr
     && (src_addr[0] == '\0'
	 || !strcmp(src_addr, "0.0.0.0")
	 || !strcmp(src_addr, "::"))) src_addr = NULL;

  in1 = itopenbuf_(&moder,irtc); if(*irtc != 0) return -1;
  fd1 = getfd_(&in1);
  if(src_addr != NULL) {
    in2 = itopenbuf_(&moder,irtc); if(*irtc != 0) return -1;
    fd2 = getfd_(&in2);
  } else {
    fd2 = 0;
  }

  /* Parse address */
  if(src_addr == NULL) { /* listen */
    fdn = -1;

    memset(&hints, 0, sizeof(hints));
    hints.ai_family = PF_UNSPEC;
    hints.ai_socktype = socktype;
    hints.ai_flags = AI_PASSIVE;

    if(src_addr0) {
      if(!strcmp(src_addr0, "0.0.0.0"))	hints.ai_family = PF_INET;
      if(!strcmp(src_addr0, "::"))	hints.ai_family = PF_INET6;
    };

    if(!(err = getaddrinfo(NULL, src_port, &hints, &res0))) {
      for(res = res0; res; res = res->ai_next) {
	fdn = socket(res->ai_family, res->ai_socktype, res->ai_protocol);
	if(fdn == -1) continue;

	if(!(err = getnameinfo(res->ai_addr, res->ai_addrlen,
			       NULL, 0, service, sizeof(service),
			       NI_NUMERICSERV)))
	  port0 = atoi(service);

	if(bind(fdn, res->ai_addr, res->ai_addrlen) != 0) {
	  close(fdn);
	  fdn = -1;
	  continue;
	}

	na = sizeof(addr);
	if(getsockname(fdn, (struct sockaddr*)&addr, &na) == 0)
	  if(!(err = getnameinfo((struct sockaddr*)&addr, na,
				 NULL, 0, service, sizeof(service),
				 NI_NUMERICSERV)))
	    port = atoi(service);

	switch(res->ai_socktype) {
	case SOCK_STREAM:
	  listen(fdn, 1);
	  break;

	default:
	  break;
	}

	break; /* listen is succeeded */
      }
      freeaddrinfo(res0);
    }
  } else { /* connect */
    fdn = -1;

    memset(&hints, 0, sizeof(hints));
    hints.ai_family = PF_UNSPEC;
    hints.ai_socktype = socktype;

    if(!(err = getaddrinfo(src_addr, src_port, &hints, &res0))) {
      for(res = res0; res; res = res->ai_next) {
	fdn = socket(res->ai_family, res->ai_socktype, res->ai_protocol);
	if(fdn == -1) continue;

	if(connect(fdn, res->ai_addr, res->ai_addrlen) != 0) {
	  close(fdn);
	  fdn = -1;
	  continue;
	}

	break; /* connection is succeeded */
      }
      freeaddrinfo(res0);
    }

    if(fdn != -1 && dup2(fdn, fd2) == -1) {
      close(fdn); fdn = -1;
    }
  }

  if(fdn != -1 && dup2(fdn, fd1) == -1) {
    close(fdn); fdn = -1;
  }

  if(fdn == -1) {
    *kx = kxfailed;
    *irtc = 0;
    return -1;
  }

  if(src_addr != NULL) {
    iax = ktavaloc(-1, 2);
    rlist(iax + 1) = in1;
    rlist(iax + 2) = in2;
    *kx = ktflist + iax;
  } else {
    if(port0 == 0) {
      iax = ktavaloc(-1, 2);
      rlist(iax + 1) = in1;
      rlist(iax + 2) = port;
      *kx = ktflist + iax;
    } else {
      vx = in1;
      *kx = kfromr(vx);
    }
  }
  close(fdn);
  *irtc = 0;
  return 0;
}

/* SADScript function definition of Network stuff */
static int tfHostName(integer4 *isp1,
		      integer8 *kx,
		      integer4 *irtc) {
  char hostname[MAXHOSTNAMELEN];

  if(isp != *isp1 + 1 || ktastk(isp) != ktfoper + mtfnull) {
    *irtc = itfmessage(9, "General::narg","\"0\"");
    return -1;
  }

  if(!gethostname(hostname, MAXHOSTNAMELEN)) {
    *kx = ktfstring + ktsalocb(-1, hostname);
  } else {
    *kx = ktfoper + mtfnull;
  }

  *irtc = 0;
  return 0;
}

static int tfGetHostAddrByName(integer4 *isp1,
			       integer8 *kx,
			       integer4 *irtc) {
  integer8 ia, kax;
  int addresses, count, err, i;
  integer8 *address, *newaddress;
  struct addrinfo hints, *res, *res0;
  char buffer[NI_MAXHOST];

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

  address = NULL; addresses = 0; count = 0;

  memset(&hints, 0, sizeof(hints));
  hints.ai_family = PF_UNSPEC;
  hints.ai_socktype = SOCK_STREAM;
  if(!(err = getaddrinfo(&jlist(1, ia + 1), NULL, &hints, &res0))) {
    for(res = res0; res; res = res->ai_next)
      if(!(err = getnameinfo(res->ai_addr, res->ai_addrlen,
			     buffer, sizeof(buffer), NULL, 0,
			     NI_NUMERICHOST))) {
	if(!(addresses < count)) {
	  count += 8;
	  newaddress = realloc(address, count * sizeof(integer4));
	  if(newaddress == NULL) break;
	  address = newaddress;
	}
	address[addresses] = ktsalocb(0, buffer); addresses += 1;
      }
    freeaddrinfo(res0);
  }

  if(addresses > 0) {
    kax = ktadaloc(-1, addresses);
    for(i = 0; i < addresses; i++)
      klist(kax + i + 1)=ktfstring + address[i];
    *kx = ktflist + kax;
  } else {
    *kx = kxnulll;
  }
  free(address);
  *irtc = 0;
  return 0;
}

static int tfGetHostNameByAddr(integer4 *isp1,
			       integer8 *kx,
			       integer4 *irtc) {
  integer4 ia;
  int err;
  struct addrinfo hints, *res, *res0;
  char buffer[NI_MAXHOST];

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

  *kx = ktfoper + mtfnull;
  *irtc = 0;

  memset(&hints, 0, sizeof(hints));
  hints.ai_family = PF_UNSPEC;
  hints.ai_socktype = SOCK_STREAM;
  hints.ai_flags = AI_NUMERICHOST;
  if(!(err = getaddrinfo(&jlist(1, ia + 1), NULL, &hints, &res0))) {
    for(res = res0; res; res = res->ai_next)
      if(!(err = getnameinfo(res->ai_addr, res->ai_addrlen,
			     buffer, sizeof(buffer), NULL, 0,
			     NI_NAMEREQD))) {
	*kx = ktfstring + ktsalocb(-1, buffer);
	break;
      }
    freeaddrinfo(res0);
  }

  if(ktfnonstringq(*kx)) return -1;

  return 0;
}

static int tfTCPOpen(integer4 *isp1,
		     integer8 *kx,
		     integer4 *irtc) {
  return tfSocketOpen(isp1, kx, irtc, SOCK_STREAM);
}

static int tfUDPOpen(integer4 *isp1,
		     integer8 *kx,
		     integer4 *irtc) {
  return tfSocketOpen(isp1, kx, irtc, SOCK_DGRAM);
}

static int tfTCPAccept(integer4 *isp1,
		       integer8 *kx,
		       integer4 *irtc) {
  integer8 kax;
  integer4 ins, in1, in2;
  int fds, fd1, fd2, fdn;
  int moder=1,modew=2;
  double timeout;
  struct sockaddr_in addr;
  socklen_t addr_len = sizeof(struct sockaddr_in);

  if(isp != *isp1 + 1 && isp != *isp1 + 2) {
    *irtc = itfmessage(9, "General::narg", "\"1 or 2\"");
    return -1;
  }

  if(isp == *isp1 + 2) {
    if(ktfnonrealq(ktastk( *isp1 + 1))) {
      *irtc = itfmessage(9, "General::wrongtype",
			 "\"Timeout(seconds) for #2\"");
      return -1;
    }
    timeout = rtastk(isp);
  } else {
    timeout = 0.0;
  }

  if(ktfnonrealq(ktastk( *isp1 + 1))) {
    *irtc = itfmessage(9, "General::wrongtype",
		       "\"Opened TCP socket logical unit number for #1\"");
    return -1;
  }
  ins = rtastk(*isp1 + 1);
  fds = getfd_(&ins);
  if(fds == -1) {
    *irtc = itfmessage(9, "General::wrongtype",
		       "\"Opened TCP socket logical unit number for #1\"");
    return -1;
  }

  in1 = itopenbuf_(&moder,irtc); if(*irtc != 0) return -1;
  fd1 = getfd_(&in1);
  in2 = itopenbuf_(&modew,irtc); if(*irtc != 0) return -1;
  fd2 = getfd_(&in2);

  if((fdn = accept(fds, (struct sockaddr*)&addr, &addr_len)) != -1) {
    if(dup2(fdn, fd1) == -1 || dup2(fdn, fd2) == -1)
      fdn = -1;
    close(fdn);
  }

  if(fdn == -1) {
    *kx = kxfailed;
    *irtc = 0;
    return -1;
  }

  kax = ktavaloc(-1, 2);
  rlist(kax + 1) = in1;
  rlist(kax + 2) = in2;
  *kx = ktflist + kax;
  *irtc = 0;
  return 0;
}

static int tfSelectUnit(integer4 *isp1,
			integer8 *kx,
			integer4 *irtc) {
  integer8 kax, ia, ki;
  integer4 inr;
  real8 vx;
  int i, fdr;
  int nfds, rfds, *rfd, status, sb_cc;
  double timeout;
  fd_set readfds;
  struct timeval tm;

  if(isp != *isp1 + 1 && isp != *isp1 + 2) {
    *irtc = itfmessage(9, "General::narg", "\"1 or 2\"");
    return -1;
  }

  if(isp == *isp1 + 2) {
    if(ktfnonrealq(ktastk( *isp1 + 2))) {
      *irtc = itfmessage(9, "General::wrongtype",
			 "\"Timeout(seconds) for #2\"");
      return -1;
    }
    timeout = rtastk(isp);
  } else {
    timeout = -1.0;
  }

  if(ktfrealq(ktastk(*isp1 + 1))){
      inr = rtastk(*isp1 + 1);
      fdr = getfd_(&inr);
      rfd = &fdr; rfds = 1;
  }
  else if(ktflistq(ktastk(*isp1 + 1))){
    ia = ktfaddr(ktastk(*isp1 + 1));
    rfds = ilist(2, ia - 1);
    rfd  = malloc(rfds * sizeof(int));
    if(rfd == NULL) {
      *irtc = itfmessage(9, "General::memoryfull", " ");
      return -1;
    }
    for(i = 0; i < rfds; i++) {
      ki = klist(ia + i + 1);
      if(ktfrealq(ki)) {
        inr = kfromr(ki);
        rfd[i] = getfd_(&inr);
      } else
        rfd[i] = -1;
    }
  }else{
    *irtc = itfmessage(9, "General::wrongtype",
                       "\"Opened logical unit number"
                       " or List of logcal unit numbers for #1\"");
    return -1;
  }

  /* Setup read file descripter set */
  nfds = 0;
  FD_ZERO(&readfds);
  for(i = 0; i < rfds; i++)
    if(rfd[i] >= 0) {
      FD_SET(rfd[i], &readfds);
      if(rfd[i] > nfds)
        nfds = rfd[i];
    }

  /* Setup timeout */
  tm.tv_sec  =  timeout;
  tm.tv_usec = (timeout - tm.tv_sec) * 1e6;

  nfds += 1;
  status = select(nfds, &readfds, NULL, NULL, (timeout < 0) ? NULL : &tm);

  /* Check file descripter status(number of character in buffer) */
  for(i = 0; i < rfds; i++)
    if(rfd[i] >= 0) {
      if(FD_ISSET(rfd[i], &readfds)) {
        status = ioctl(rfd[i], FIONREAD, &sb_cc);
        rfd[i] = (status != -1) ? sb_cc : 0;
      } else
        rfd[i] = 0;
    } else
      rfd[i] = 0;

  if(ktfrealq(ktastk(*isp1 + 1))){
  }
  else if(ktflistq(ktastk(*isp1 + 1))){
    kax = ktavaloc(-1, rfds);
    for(i = 0; i < rfds; i++)
      rlist(kax + i + 1)=rfd[i];
    free(rfd);
    *kx = ktflist + kax;
  }
  else{
    vx = rfd[0];
    *kx = kfromr(vx);
  };

  *irtc = 0;
  return 0;
}

/* SADScript function registration of Network stuff */
#define REG	dlfunaloc
#define REG8	dlfunaloc8
int sadDefFunc_NetworkIO(void) {
  REG8("HostName",		tfHostName,		0, NULL, NULL, 0);

  /* Resolver handling functions */
  REG8("GetHostAddrByName",	tfGetHostAddrByName,	1, NULL, NULL, 0);
  REG8("GetHostNameByAddr",	tfGetHostNameByAddr,	1, NULL, NULL, 0);

  REG8("GetHostByAddr",		tfGetHostNameByAddr,	1, NULL, NULL, 0);

  /* Socket handling funcions */
  REG8("TCPOpen",		tfTCPOpen,		3, NULL, NULL, 0);
  REG8("UDPOpen",		tfUDPOpen,		3, NULL, NULL, 0);
  REG8("TCPAccept",		tfTCPAccept,		2, NULL, NULL, 0);
  REG8("SelectUnit",		tfSelectUnit,		2, NULL, NULL, 0);

  return 0;
}

/* End of File */
