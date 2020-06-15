#include <sim/sad_api.h>
#include <sim/TFCODE.h>
#include <sim/TFSTK.h>
#include <sim/TFRBUF.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <fcntl.h>
#include <sys/mman.h>

#define	TEMP_TYPE_REGULAR	0
#define	TEMP_TYPE_SPECIAL	1
#define	TEMP_TYPE_TTY		2

#undef	DEBUG_TEMPORARY_FILES

/* Prototype for sim/fortran_io.f */
extern void close_lun_(integer4*);
extern integer4 open_read_(const_character, integer4*, ftnlen);
extern integer4 open_write_(const_character, integer4*, ftnlen);
extern integer4 open_append_(const_character, integer4*, ftnlen);

/* SAD open buffer API */
integer4 itopenbuf_(integer4 *mode, integer4 *irtc) {
  integer4 nc, in;
  integer4 argv[2] = {2, 0};
  char *template, template_buffer[32];
  int fd, dev_null = -1, type = TEMP_TYPE_REGULAR;

  strncpy(template_buffer,
	  "/tmp/openbuf.XXXXXXXXXXXXXXXXXXX", sizeof(template_buffer));
  template_buffer[sizeof(template_buffer) - 1] = '\0';

  template = template_buffer;
  nc = strlen(template);

#ifdef	USE_ITOPENBUF_MKFIFO
  type = TEMP_TYPE_SPECIAL;
#else
  type = TEMP_TYPE_TTY;
#endif

  fprintf(stderr,"itopenbuf %d",type);
  switch(type) {
  case TEMP_TYPE_REGULAR:
    fprintf(stderr,"REGULAR");
    fd = mkstemp(template);
    if(fd < 0) {
      *irtc = itfmessage(999, "System::error",
			 "\"Can not create temporary file\"");
      return 0;
    }
    close(fd);
    break;

  case TEMP_TYPE_SPECIAL:
    fprintf(stderr,"SPECIAL");
    template[nc - 5] = '\0';
    strcat(template, "/fifo");
    template[nc - 5] = '\0';

    if(mkdtemp(template) == NULL) {
      *irtc = itfmessage(999, "System::error",
			 "\"Can not create temporary file\"");
      return 0;
    }

    template[nc - 5] = '/';
    fd = mkfifo(template, S_IRUSR | S_IWUSR);
    if(fd == -1) {
      template[nc - 5] = '\0';
      rmdir(template);
      *irtc = itfmessage(999, "System::error",
			 "\"Can not create temporary file\"");
      return 0;
    }
    break;

  case TEMP_TYPE_TTY:
    fprintf(stderr,"TTY");
    dev_null = open("/dev/null", O_RDWR);
    if(dev_null == -1) {
      *irtc = itfmessage(999, "System::error",
			 "\"Can not create temporary file\"");
      return 0;
    }

    fd = posix_openpt(O_RDWR|O_NOCTTY);
    if(fd == -1) {
      close(dev_null);
      *irtc = itfmessage(999, "System::error",
			 "\"Can not create temporary file\"");
      return 0;
    }

    if(grantpt(fd)  != 0 ||
       unlockpt(fd) != 0 ||
       (template = ptsname(fd)) == NULL) {
      close(dev_null);
      close(fd);
      *irtc = itfmessage(999, "System::error",
			 "\"Can not create temporary file\"");
      return 0;
    }

    nc = strlen(template);
    break;

  default:
    *irtc = itfmessage(999, "System::error",
		       "\"Can not create temporary file\"");
    return 0;
  }

  in = open_read_(template, &nc, nc);
  switch(type) {
#ifndef	DEBUG_TEMPORARY_FILES
  case TEMP_TYPE_REGULAR:
    unlink(template);
    break;

  case TEMP_TYPE_SPECIAL:
    unlink(template);
    template[nc - 5] = '\0';
    rmdir(template);
    break;
#endif

  case TEMP_TYPE_TTY:
    if(in > 0 && dup2(dev_null, getfd_(&in)) == -1) {
      close_lun_(&in);
      in = -1;
    }
    close(dev_null);
    close(fd);
    break;

  default:
    break;
  }
  if(in < 0) {
    *irtc = itfmessage(999, "General::toomany",
		       "\"files opened\"");
    return 0;
  }

#ifdef DEBUG_ITOPENBUF
  fprintf(stderr, "itopenbuf(unit=%d, fd=%d, temp=\"%s\")\n", in, getfd_(&in), template);
#endif

  argv[0]=*mode;
  trbinit(in, &argv[0]);

  *irtc = 0;
  return in;
}

/* SAD itfopen* API family */
integer4 itfopenread_(integer8 *k, logical4 *disp,
		      integer4 *irtc) {
  integer8 ka;
  integer4 nc, in;
  integer4 argv[2] = {2, 0};
  char *expanded, *quote = NULL;

  if((*k & ktfmask) != ktfstring) {
    *irtc=itfmessage(9, "General::wrongval",
		     "\"Filename\",\"argument\"");
    return -1;
  }

  ka = (*k & ktamask);
#if SAD_REQUIRE_STRING_TERMINATION
  jlist(ilist(1, ka) + 1, ka + 1) = '\0';
#endif
  expanded = expand_tilde(&jlist(1, ka + 1));

  nc = strlen(expanded);
  in = open_read_(expanded, &nc, nc);
  if(in < 0) {
    asprintf(&quote, "\"%s\"", expanded); free(expanded);
    if(quote != NULL) {
      *irtc = itfmessage(999, "General::fileopen", quote);
      free(quote);
    } else {
      *irtc = itfmessage(999, "General::fileopen", "\"(unknown)\"");
    }
    return -2;
  }

  if(*disp)
    unlink(expanded);

  trbinit(in, &argv[0]);

  free(expanded);
  *irtc = 0;
  return in;
}

integer4 itfopenwrite_(integer4 *it, integer4 *ia, integer4 *irtc) {
  integer4 nc, in;
  char *expanded, *quote = NULL;

  if(*it != ntfstring) {
    *irtc=itfmessage(9, "General::wrongval",
		     "\"Filename\",\"argument\"");
    return -1;
  }

#if SAD_REQUIRE_STRING_TERMINATION
  jlist(ilist(1, *ia) + 1, *ia + 1) = '\0';
#endif
  expanded = expand_tilde(&jlist(1, *ia + 1));

  nc = strlen(expanded);
  in = open_write_(expanded, &nc, nc);
  if(in < 0) {
    asprintf(&quote, "\"%s\"", expanded); free(expanded);
    if(quote != NULL) {
      *irtc = itfmessage(999, "General::fileopen", quote);
      free(quote);
    } else {
      *irtc = itfmessage(999, "General::fileopen", "\"(unknown)\"");
    }
    return -2;
  }

  free(expanded);
  *irtc = 0;
  return in;
}

integer4 itfopenappend_(integer4 *it, integer4 *ia, integer4 *irtc) {
  integer4 nc, in;
  char *expanded, *quote = NULL;

  if(*it != ntfstring) {
    *irtc=itfmessage(9, "General::wrongval",
		     "\"Filename\",\"argument\"");
    return -1;
  }

#if SAD_REQUIRE_STRING_TERMINATION
  jlist(ilist(1, *ia) + 1, *ia + 1) = '\0';
#endif
  expanded = expand_tilde(&jlist(1, *ia + 1));

  nc = strlen(expanded);
  in = open_append_(expanded, &nc, nc);
  if(in < 0) {
    asprintf(&quote, "\"%s\"", expanded); free(expanded);
    if(quote != NULL) {
      *irtc = itfmessage(999, "General::fileopen", quote);
      free(quote);
    } else {
      *irtc = itfmessage(999, "General::fileopen", "\"(unknown)\"");
    }
    return -2;
  }

  free(expanded);
  *irtc = 0;
  return in;
}

/* End of File */
