#include <sim/sad_f2c.h>

#include <sys/types.h>
#include <sys/uio.h>
#include <unistd.h>

/*
 * CAUTION: Don't used both Fortran READ statement and following
 *          itfgetbuf_/fgetc_ with same logical unit number!
 */

/* C implementation of itfgetbuf_() */
integer4 itfgetbuf_(integer4 *unit,
		    character buf, integer4 *limit, integer4 *irtc) {
  const char eol = '\n';
  char *pbuf;
  int fd, nc;

  *irtc = 0;	/* Reset status code */

  fd = getfd_(unit);
  for(pbuf = buf, nc = 0; nc < *limit; pbuf += 1, nc += 1) {
    switch(read(fd, pbuf, 1)) {	/* Try to read 1-character */
    case -1:	/* read(2) is failed */
      *irtc = 1;	/* Set error status(1) */
      return -999;	/* Return error code(-999) */
      break;

    case  0:	/* read(2) detects End of File */
      if(nc == 0) {	/* If EoF is detected without buffered characters, */
	*irtc = -1;	/* set EoF status(-1) */
	return -99;	/* Return EoF code(-99) */
      }
      *pbuf = eol;	/* Add tailing EoL into buffer */
      return nc;	/* Return buffered length without EoL character */
      break;

    default:	/* read(2) is succeeded */
      if(*pbuf == eol) {/* If End of Line is detected, */
	if(nc > 0 && pbuf[-1] == '\r' ) {	/* "\r\n" support */
	  return nc-1;	/* return buffered length without EoL character */
	} else {
	  return nc;	/* return buffered length without EoL character */
	}
      }
      break;
    }
  }

  /* Return buffer length(*limit), because buffer is filled! */
  return nc;	/* nc MIGHT be *limit */
}

/* C implementation of fgetc_() */
integer4 fgetc_(integer4 *unit, character buf) {
  switch(read(getfd_(unit), buf, 1)) {
  case -1:	/* read(2) is failed */
    return  1;
    break;

  case  0:	/* read(2) detects End of File */
    return -1;
    break;

  default:	/* read(2) is succeeded */
    return 0;
    break;
  }
}

/* End of File */
