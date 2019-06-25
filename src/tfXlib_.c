#include <sim/sad_xlib.h>
#include <sim/sad_api.h>
#include <sim/TFCODE.h>
#include <sim/TFCBK.h>
#include <sim/TFSTK.h>

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <strings.h>
#include <unistd.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xatom.h>

/* Maximum number of font names to be returned */
#define MAX_FONTLIST_HITS	10000

/*
 * This part contains most of the X calls called by SAD. Many of
 * these calls are just stubs and either don't make sense on the
 * Macintosh or thier implamentation just doesn't do anything.
 */

#if defined(AQUA_TK)

#include <tcl.h>
 
#define XFree(data) {if ((data) != NULL) ckfree((char *) (data));}

int XFlush(Display *display) {
  return False;
}

int XRestackWindows(Display *display, Window windows[], int nwindows) {
  return BadWindow;
}

char **XListFonts(Display *display, const char *pattern, int maxnames,
		  int *actual_count_return) {
  actual_count_return = 0;
  return NULL;
}

int XFreeFontNames(char **list) {
  return 1;
}

char **XGetFontPath(Display *display, int *npaths_return) {
  npaths_return = 0;
  return NULL;
}

int XFreeFontPath(char **list) {
  return 1;
}

Window
XDefaultRootWindow (dpy)
     Display *dpy;
{
  return (RootWindow(dpy,DefaultScreen(dpy)));
}

void
XSetTextProperty (dpy, w, tp, property)
     Display *dpy;
     Window w;
     Atom property;
     XTextProperty *tp;
{
  XChangeProperty (dpy, w, property, tp->encoding, tp->format,
                   PropModeReplace, tp->value, tp->nitems);
}

Status
XGetTextProperty (display, window, tp, property)
     Display *display;
     Window window;
     XTextProperty *tp;
     Atom property;
{
  Atom actual_type;
  int actual_format = 0;
  unsigned long nitems = 0L, leftover = 0L;
  unsigned char *prop = NULL;

  if (XGetWindowProperty (display, window, property, 0L, 1000000L, False,
                          AnyPropertyType, &actual_type, &actual_format,
                          &nitems, &leftover, &prop) == Success &&
      actual_type != None) {
    /* okay, fill it in */
    tp->value = prop;
    tp->encoding = actual_type;
    tp->format = actual_format;
    tp->nitems = nitems;
    return True;
  }

  tp->value = NULL;
  tp->encoding = None;
  tp->format = 0;
  tp->nitems = 0;
  return False;
}

Status
XGetWindowAttributes(dpy, w, attr)
     register Display *dpy;
     Window w;
     XWindowAttributes *attr;
{  
  return 0;
}

int
XDeleteProperty(dpy, window, property)
     Display *dpy;
     Window window;
     Atom property;
{
  return 0;
}

#endif /* AQUA_TK */

/* end of stub part */

/* Static variable to keep X server connection */
static Display *curDisplay = NULL;

static bool hasDisplay(void) {
  const char *display;

  if(curDisplay != NULL)
    return true;

  display = getenv("DISPLAY");
  if(display != NULL)
    curDisplay = XOpenDisplay(display);

  return curDisplay != NULL;
}

/* X server connection update API for SAD/Tkinter */
static bool setNewDisplay(void *new0) {
  Display *new = (Display *)new0;

  if(new == NULL) return false;

  if(curDisplay != NULL && curDisplay != new)
    XCloseDisplay(curDisplay);

  curDisplay = new;

  return true;
}

/* SADScript function definition of Xlib stuff */
int tfXServerVendor(integer4 *isp1,
		    integer8 *kx,
		    integer4 *irtc) {
  char *sv;

  if(isp != *isp1 + 1
     || ktastk(isp) != ktfoper + mtfnull) {
    *irtc = itfmessage(9, "General::narg","\"0\"");
    return -1;
  }

  if(!hasDisplay()) {
    *kx = kxfailed;
    *irtc = 0;
    return -1;
  }

  sv = ServerVendor(curDisplay);
  if(sv != NULL) {
    *kx = ktfstring + ktsalocb(-1, sv);
  } else {
    *kx = ktfoper + mtfnull;
  }

  *irtc = 0;
  return 0;
}

int tfXFlush(integer4 *isp1,
	     integer8 *kx,
	     integer4 *irtc) {

  if(isp != *isp1 + 1
     || ktastk(isp) != ktfoper + mtfnull) {
    *irtc = itfmessage(9, "General::narg","\"0\"");
    return -1;
  }

  if(!hasDisplay()) {
#ifdef COMPAT_XLIB_MAIN_TRUNK
    *kx = ktfoper + mtfnull;
#else
    *kx = kxfailed;
#endif
    *irtc = 0;
    return -1;
  }

  XFlush(curDisplay);
  *kx = ktfoper + mtfnull;
  *irtc = 0;
  return 0;
}

int tfXRaiseWindow(integer4 *isp1,
		   integer8 *kx,
		   integer4 *irtc) {
  Window w;

  if(isp != *isp1 + 1) {
    *irtc = itfmessage(9, "General::narg", "\"1\"");
    return -1;
  }

  if((ktrmask & ktastk(isp)) == ktfnr) {
    *irtc = itfmessage(9, "General::wrongtype", "\"Real Number\"");
    return -1;
  }


  if(!hasDisplay()) {
#ifdef COMPAT_XLIB_MAIN_TRUNK
    *kx = kfromr(-1.0);
#else
    *kx = kxfailed;
#endif
    *irtc = 0;
    return -1;
  }

  w = rtastk(isp);
  *kx  = (integer8) XRaiseWindow(curDisplay, w); XFlush(curDisplay);
  *irtc = 0;
  return 0;
}

int tfXLowerWindow(integer4 *isp1,
		   integer8 *kx,
		   integer4 *irtc) {
  Window w;

  if(isp != *isp1 + 1) {
    *irtc = itfmessage(9, "General::narg", "\"1\"");
    return -1;
  }

  if((ktrmask & ktastk(isp)) == ktfnr) {
    *irtc = itfmessage(9, "General::wrongtype", "\"Real Number\"");
    return -1;
  }


  if(!hasDisplay()) {
#ifdef COMPAT_XLIB_MAIN_TRUNK
    *kx  = kfromr(-1.0);
#else
    *kx = kxfailed;
#endif
    *irtc = 0;
    return -1;
  }

  w = rtastk(isp);

  *kx  = (integer8) XLowerWindow(curDisplay, w); XFlush(curDisplay);
  *irtc = 0;
  return 0;
}

int tfXRestackWindows(integer4 *isp1,
		      integer8 *kx,
		      integer4 *irtc) {
  integer8 ia, ki;
  Window *windows;
  int i, nwindows;

  if(isp != *isp1 + 1) {
    *irtc = itfmessage(9, "General::narg", "\"1\"");
    return -1;
  }

  if((ktfmask & ktastk(isp)) != ktflist) {
    *irtc = itfmessage(9, "General::wrongtype", "\"List\"");
    return -1;
  }

  if(!hasDisplay()) {
#ifdef COMPAT_XLIB_MAIN_TRUNK
    *kx = 0;
#else
    *kx = kxfailed;
#endif
    *irtc = 0;
    return -1;
  }

    ia = (ktamask & ktastk(isp));
  nwindows = ilist(2, ia - 1);

  windows = malloc(nwindows * sizeof(Window));
  if(windows == NULL) {
#ifdef COMPAT_XLIB_MAIN_TRUNK
    *kx = 0;
#else
    *kx = kxfailed;
#endif
    *irtc = 0;
    return -1;
  }

  /* Duplicate Window array */
  for(i = 0; i < nwindows; i++) {
    ki = klist(ia + i + 1);
    if(ktfrealq(ki)) {
      windows[i] = rfromk(ki);
    } else {
      free(windows);
      *irtc = itfmessage(9, "General::wrongtype", "\"List of Reals\"");
      return -1;
    }
  }

  XRestackWindows(curDisplay, windows, nwindows); free(windows);
  XFlush(curDisplay);

    *kx = 0;
    *irtc = 0;
    return 0;
}

int tfXInternAtom(integer4 *isp1,
		  integer8 *kx,
		  integer4 *irtc) {
  integer8 ia;
  Bool only_if_exists;

  if(isp != *isp1 + 2) {
    *irtc = itfmessage(9, "General::narg", "\"2\"");
    return -1;
  }

  if((ktfmask & ktastk(*isp1 + 1)) != ktfstring) {
    *irtc = itfmessage(9, "General::wrongtype", "\"Character-string for #1\"");
    return -1;
  }

  if((ktrmask & ktastk(isp)) == ktfnr) {
    *irtc = itfmessage(9, "General::wrongtype", "\"Real Number for #2\"");
    return -1;
  }

  if(!hasDisplay()) {
#ifdef COMPAT_XLIB_MAIN_TRUNK
    *kx =0;
#else
    *kx = kxfailed;
#endif 
    *irtc = 0;
    return -1;
  }

  ia = (ktamask & ktastk(*isp1 + 1));
#if SAD_REQUIRE_STRING_TERMINATION
  jlist(ilist(1, ia) + 1, ia + 1) = '\0';
#endif

  only_if_exists = rtastk(isp);

  real8 vx = XInternAtom(curDisplay, &jlist(1, ia + 1), only_if_exists);
    *kx  = kfromr(vx);
  *irtc = 0;
  return 0;
}

int tfXSetTextProperty(integer4 *isp1,
		       integer8 *kx,
		       integer4 *irtc) {
  integer8 ia;
  Window w;
  Atom property;
  XTextProperty text_prop;

  if(isp != *isp1 + 3) {
    *irtc = itfmessage(9, "General::narg", "\"3\"");
    return -1;
  }

  if((ktrmask & ktastk(*isp1 + 1)) == ktfnr) {
    *irtc = itfmessage(9, "General::wrongtype", "\"Real Number for #1\"");
    return -1;
  }

  if((ktrmask & ktastk(*isp1 + 2)) == ktfnr) {
    *irtc = itfmessage(9, "General::wrongtype", "\"Real Number for #2\"");
    return -1;
  }

  if((ktfmask & ktastk(isp)) != ktfstring) {
    *irtc = itfmessage(9, "General::wrongtype", "\"Character-string for #3\"");
    return -1;
  }

  if(!hasDisplay()) {
#ifdef COMPAT_XLIB_MAIN_TRUNK
    *kx =0;
#else
    *kx = kxfailed;
#endif 
    *irtc = 0;
    return -1;
  }

  w        = rtastk(*isp1 + 1);
  property = rtastk(*isp1 + 2);
  ia = (ktamask & ktastk(isp));
#if SAD_REQUIRE_STRING_TERMINATION
  jlist(ilist(1, ia) + 1, ia + 1) = '\0';
#endif

  text_prop.value    = (unsigned char *)&jlist(1, ia + 1);
  text_prop.encoding = XA_STRING;
  text_prop.format   = 8;
  text_prop.nitems   = ilist(1, ia); /* length of ntfstring */
  XSetTextProperty(curDisplay, w, &text_prop, property);
  XFlush(curDisplay);

  *kx = 0;
  *irtc = 0;
  return 0;
}

int tfXGetTextProperty(integer4 *isp1,
		       integer8 *kx,
		       integer4 *irtc) {
  Status status;
  Window w;
  Atom property;
  XTextProperty text_prop;

  if(isp != *isp1 + 2) {
    *irtc = itfmessage(9, "General::narg", "\"2\"");
    return -1;
  }

  if((ktrmask & ktastk(*isp1 + 1)) == ktfnr) {
    *irtc = itfmessage(9, "General::wrongtype", "\"Real Number for #1\"");
    return -1;
  }

  if((ktrmask & ktastk(isp)) == ktfnr) {
    *irtc = itfmessage(9, "General::wrongtype", "\"Real Number for #2\"");
    return -1;
  }

  if(!hasDisplay()) {
#ifdef COMPAT_XLIB_MAIN_TRUNK
    *kx = ktfoper + mtfnull;
#else
    *kx = kxfailed;
#endif 
    *irtc = 0;
    return -1;
  }

  w        = rtastk(*isp1 + 1);
  property = rtastk(isp);
  status = XGetTextProperty(curDisplay, w, &text_prop, property);

  if(status == True && text_prop.value != NULL) {
    *kx = ktfstring + ktsalocb(-1, (char *)text_prop.value);
  } else {
    *kx = ktfoper + mtfnull;
  }
  *irtc = 0;
  return 0;
}
int tfXChangeProperty(integer4 *isp1,
		      integer8 *kx,
		      integer4 *irtc) {
  integer8 ia;
  integer4 nc, itype, imode;
  Window w;
  Atom property, type;
  int format, mode;

  if(isp != *isp1 + 5) {
    *irtc = itfmessage(9, "General::narg", "\"5\"");
    return -1;
  }

  if((ktrmask & ktastk( *isp1 + 1)) == ktfnr) {
    *irtc = itfmessage(9, "General::wrongtype", "\"Real Number for #1\"");
    return -1;
  }

  if((ktrmask & ktastk( *isp1 + 2)) == ktfnr) {
    *irtc = itfmessage(9, "General::wrongtype", "\"Real Number for #2\"");
    return -1;
  }

  if((ktrmask & ktastk( *isp1 + 3)) == ktfnr) {
    *irtc = itfmessage(9, "General::wrongtype", "\"Real Number for #3\"");
    return -1;
  }

  if((ktrmask & ktastk( *isp1 + 4)) == ktfnr) {
    *irtc = itfmessage(9, "General::wrongtype", "\"Real Number for #4\"");
    return -1;
  }

  if(!hasDisplay()) {
#ifdef COMPAT_XLIB_MAIN_TRUNK
    *kx = 0;
#else
    *kx = kxfailed;
#endif 
    *irtc = 0;
    return -1;
  }

  w        = rtastk(*isp1 + 1);
  property = rtastk(*isp1 + 2);
  itype    = rtastk(*isp1 + 3);
  imode    = rtastk(*isp1 + 4);

  switch(itype) {
  case 1:
    type = XA_STRING;
    format = 8;
    if((ktfmask & ktastk(isp)) != ktfstring) {
      *irtc = itfmessage(9, "General::wrongtype",
			 "\"Character-string for #5\"");
      return -1;
    }
    break;

  case 2: /* need fix? */
    type = XA_INTEGER;
    format = 32;
    if((ktfmask & ktastk(isp)) != ktfstring) {
      *irtc = itfmessage(9, "General::wrongtype",
			 "\"Character-string for #5\"");
      return -1;
    }
    break;

  default:
    *irtc = itfmessage(9, "General::wrongtype", "\"1 or 2 for #3\"");
    return -1;
  }

  switch(imode) {
  case 1:
    mode = PropModeReplace;
    break;

  case 2:
    mode = PropModePrepend;
    break;

  case 3:
    mode = PropModeAppend;
    break;

  default:
    *irtc = itfmessage(9, "General::wrongtype", "\"1, 2, or 3 for #4\"");
    return -1;
  }

  ia = (ktamask & ktastk(isp));
  nc = ilist(1, ia);
#if SAD_REQUIRE_STRING_TERMINATION
  jlist(ilist(1, ia) + 1, ia + 1) = '\0';
#endif

  XChangeProperty(curDisplay, w, property, type, format, mode,
			   (unsigned char *)&jlist(1, ia + 1), nc);
  *kx = 0;
  *irtc = 0;
  return 0;
}

int tfXDeleteProperty(integer4 *isp1,
		      integer8 *kx,
		      integer4 *irtc) {
  Window w;
  Atom property;

  if(isp != *isp1 + 2) {
    *irtc = itfmessage(9, "General::narg", "\"2\"");
    return -1;
  }

  if((ktrmask & ktastk( *isp1 + 1)) == ktfnr) {
    *irtc = itfmessage(9, "General::wrongtype", "\"Real Number for #1\"");
    return -1;
  }

  if((ktrmask & ktastk( isp)) == ktfnr) {
    *irtc = itfmessage(9, "General::wrongtype", "\"Real Number for #2\"");
    return -1;
  }

  if(!hasDisplay()) {
#ifdef COMPAT_XLIB_MAIN_TRUNK
    *kx = 0;
#else
    *kx = kxfailed;
#endif 
    *irtc = 0;
    return -1;
  }

  w        = rtastk(*isp1 + 1);
  property = rtastk( isp);

  XDeleteProperty(curDisplay, w, property);
  *kx = 0;
  *irtc = 0;
  return 0;
}

int tfXGetWindowAttributes(integer4 *isp1,
			   integer8 *kx,
			   integer4 *irtc) {
  Status status;
  Window w;
  XWindowAttributes xwa;

  if(isp != *isp1 + 1) {
    *irtc = itfmessage(9, "General::narg", "\"1\"");
    return -1;
  }

  if((ktrmask & ktastk( isp)) == ktfnr) {
    *irtc = itfmessage(9, "General::wrongtype", "\"Real Number\"");
    return -1;
  }

  if(!hasDisplay()) {
#ifdef COMPAT_XLIB_MAIN_TRUNK
    *kx = ktfoper + mtfnull;
#else
    *kx = kxfailed;
#endif 
    *irtc = 0;
    return -1;
  }

  w = rtastk(isp);
  status = XGetWindowAttributes(curDisplay, w, &xwa);
  printf("XGetWindowAttributes: %d, %d, %d, %d\n", status,
	 xwa.backing_store, xwa.width, xwa.height);
  *kx = ktfoper + mtfnull;
  *irtc = 0;
  return 0;
}

int tfXGetWindowProperty(integer4 *isp1,
			 integer8 *kx,
			 integer4 *irtc) {
  integer2 *i2ptr;
  integer4 *i4ptr;
  integer8 kax;
  char      *cptr;
  int i, length, n, status;
  Window w;
  Atom property;
  Bool delete;
  Atom actual_type;
  int actual_format;
  unsigned long nitems, bytes_after;
  unsigned char *prop;

  if(isp != *isp1 + 3) {
    *irtc = itfmessage(9, "General::narg", "\"3\"");
    return -1;
  }

  if((ktrmask & ktastk( *isp1 + 1)) == ktfnr) {
    *irtc = itfmessage(9, "General::wrongtype", "\"Real Number for #1\"");
    return -1;
  }

  if((ktrmask & ktastk( *isp1 + 2)) == ktfnr) {
    *irtc = itfmessage(9, "General::wrongtype", "\"Real Number for #2\"");
    return -1;
  }

  if((ktrmask & ktastk( isp)) == ktfnr) {
    *irtc = itfmessage(9, "General::wrongtype", "\"Real Number for #3\"");
    return -1;
  }

  if(!hasDisplay()) {
#ifdef COMPAT_XLIB_MAIN_TRUNK
    *kx = ktfoper + mtfnull;
#else
    *kx = kxfailed;
#endif 
    *irtc = 0;
    return -1;
  }

  w        = rtastk(*isp1 + 1);
  property = rtastk(*isp1 + 2);
  delete   = rtastk( isp);

  status = XGetWindowProperty(curDisplay, w, property, 0, 1L<<16,
			      delete, AnyPropertyType, &actual_type,
			      &actual_format, &nitems,
			      &bytes_after, &prop);

  if(status != Success) {
    *kx = ktfoper + mtfnull;
    *irtc = 0;
    XFree(prop);
    return status;
  }

  switch(actual_format) {
  case 32:
    kax = ktavaloc(-1, nitems);
    i4ptr = (integer4 *)prop;
    for(i = 0; i < nitems; i++)
      rlist(kax + 1 + i) = i4ptr[i];
    *kx = ktflist + kax;
    break;

  case 16:
    kax = ktavaloc(-1, nitems);
    i2ptr = (integer2 *)prop;
    for(i = 0; i < nitems; i++)
      rlist(kax + 1 + i) = i2ptr[i];
    *kx = ktflist + kax;
    break;

  case 8:
    cptr = (char *)prop;
    if(nitems > 0 && cptr[nitems - 1] == '\0') {
      n = 0; i = 0;
      while(i < nitems) {
	i += strlen(cptr + i) + 1;
	n += 1;
      }
      kax = ktadaloc(-1, n);
      for(i = 0; i < n; i++) {
	length = strlen(cptr);
        klist(kax + i + 1)= ktfstring + ktsalocbl(0, cptr, length);
	cptr += length + 1;
      }
    } else {
      kax = ktadaloc(-1, 1);
      klist(kax + 1) = ktfstring + ktsalocbl(0, cptr, nitems);
    }
    *kx = ktflist + kax;
    break;

  default:
    *kx = ktfoper + mtfnull;
  }

  *irtc = 0;
  XFree(prop);
  return status;
}

int tfXQueryTree(integer4 *isp1,
		 integer8 *kx,
		 integer4 *irtc) {
  integer8 ka, kax;
  Status status;
  Window w, root, parent;
  Window *children;
  unsigned int nchildren;
  int i;

  if(isp != *isp1 + 1) {
    *irtc = itfmessage(9, "General::narg", "\"1\"");
    return -1;
  }

  if((ktrmask & ktastk( isp)) == ktfnr) {
    *irtc = itfmessage(9, "General::wrongtype", "\"Real Number\"");
    return -1;
  }

  if(!hasDisplay()) {
#ifdef COMPAT_XLIB_MAIN_TRUNK
    *kx - ktfoper + mtfnull;
#else
    *kx = kxfailed;
#endif 
    *irtc = 0;
    return -1;
  }

  w = rtastk(isp);
  if(w == 0)
    w = XDefaultRootWindow(curDisplay);

  status = XQueryTree(curDisplay, w, &root, &parent, &children, &nchildren);
  if(status != True) {
    *kx = ktfoper + mtfnull;
    *irtc = 0;
    return -1;
  }

  kax = ktadaloc(-1, 3);
  rlist(kax + 1) = root;
  rlist(kax + 2) = parent;
  if(nchildren > 0 && children != NULL) {
    ka = ktavaloc(0, nchildren);
    for(i = 0; i < nchildren; i++)
      rlist(ka + 1 + i) = children[i];
    klist(kax + 3) = ktflist + ka;
  } else {
    klist(kax + 3) = ktfcopy(kxnulll);
  }
  *kx = ktflist + kax;
  *irtc = 0;
  return 0;
}

int tfXListFonts(integer4 *isp1,
		 integer8 *kx,
		 integer4 *irtc) {
  integer8 ia, kax;
  integer4 n;
  char **pcl;
  int i;

  if(isp != *isp1 + 1) {
    *irtc = itfmessage(9, "General::narg", "\"1\"");
    return -1;
  }

  if(ktfnonstringq(ktastk(isp))) {
    *irtc = itfmessage(9, "General::wrongtype", "\"Character-string\"");
    return -1;
  }

  if(!hasDisplay()) {
    *kx = kxfailed;
    *irtc = 0;
    return -1;
  }

  ia = ktfaddr(ktastk(*isp1 + 1));
#if SAD_REQUIRE_STRING_TERMINATION
  jlist(ilist(1, ia) + 1, ia + 1) = '\0';
#endif

  pcl = XListFonts(curDisplay, &jlist(1, ia + 1), MAX_FONTLIST_HITS, &n);
  if(pcl != NULL) {
    kax = ktadaloc(-1, n);
    *irtc = 0;
    for(i = 0; i < n; i++)
      klist(kax + i + 1) = ktfstring + ktsalocb(0, pcl[i]);
    XFreeFontNames(pcl);
    *kx = ktflist + kax;
    return 0;
  }

  /* return zero length list */
  *kx = kxnulll;
  *irtc = 0;
  return 0;
}

int tfXGetFontPath(integer4 *isp1,
		   integer8 *kx,
		   integer4 *irtc) {
  integer8 kax;
  integer4 n;
  char **pcl;
  int i;

  if(isp != *isp1 + 1
     || ktastk(isp) != ktfoper + mtfnull) {
    *irtc = itfmessage(9, "General::narg","\"0\"");
    return -1;
  }

  if(!hasDisplay()) {
    *kx = kxfailed;
    *irtc = 0;
    return -1;
  }

  pcl = XGetFontPath(curDisplay, &n);
  if(pcl != NULL) {
    kax = ktadaloc(-1, n);
    *irtc = 0;
    for(i = 0; i < n; i++)
      klist(kax + i + 1) = ktfstring + ktsalocb(0, pcl[i]);
    *kx = ktflist + kax;
    XFreeFontPath(pcl);
    return 0;
  }

  /* return zero length list */
  *kx = kxnulll;
  *irtc = 0;
  return 0;
}

/* SADScript function registration of Xlib stuff */
#define REG	dlfunaloc
#define REG8	dlfunaloc8
#ifdef WITH_EXTENSION_MODULE
int dldeffun_(void) {
#else
int sadDefFunc_Xlib(void) {
#endif /* WITH_EXTENSION_MODULE */
  /* Linking SetNewDisplay handler */
  sadxlib_hooks.SetNewDisplay = setNewDisplay;

  REG8("XServerVendor",		tfXServerVendor,	0, NULL, NULL, 0);
  REG8("XFlush",			tfXFlush,		0, NULL, NULL, 0);
  REG8("XRaiseWindow",		tfXRaiseWindow,		1, NULL, NULL, 0);
  REG8("XLowerWindow",		tfXLowerWindow,		1, NULL, NULL, 0);
  REG8("XRestackWindows",	tfXRestackWindows,	1, NULL, NULL, 0);
  REG8("XInternAtom",		tfXInternAtom,		2, NULL, NULL, 0);
  REG8("XSetTextProperty",	tfXSetTextProperty,	3, NULL, NULL, 0);
  REG8("XGetTextProperty",	tfXGetTextProperty,	2, NULL, NULL, 0);
  REG8("XChangeProperty",	tfXChangeProperty,	5, NULL, NULL, 0);
  REG8("XDeleteProperty",	tfXDeleteProperty,	2, NULL, NULL, 0);
  REG8("XGetWindowAttributes",	tfXGetWindowAttributes,	1, NULL, NULL, 0);
  REG8("XGetWindowProperty",	tfXGetWindowProperty,	3, NULL, NULL, 0);
  REG8("XQueryTree",		tfXQueryTree,		1, NULL, NULL, 0);
  REG8("XListFonts",		tfXListFonts,		1, NULL, NULL, 0);
  REG8("XGetFontPath",		tfXGetFontPath,		0, NULL, NULL, 0);

  return 0;
}

/* End of File */
