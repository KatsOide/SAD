#include <sim/sad_tcltk.h>
#include <sim/sad_xlib.h>
#include <sim/sad_api.h>
#include <sim/TFCODE.h>
#include <sim/TFCBK.h>
#include <sim/TFSTK.h>

#include <sys/types.h>
#include <sysexits.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <string.h>
#include <unistd.h>

/* Define USE_INTERP_RESULT macro in order to use Tcl_Interp.result member
 * This macro is required 8.6a3 or later
 */
#define	USE_INTERP_RESULT

#include <tcl.h>
#include <tk.h>

/* Stub function/variable of Tk extensions for SAD */
#if (TK_MAJOR_VERSION >= 8) && (TK_MINOR_VERSION >= 1)
#  define SAD_TK_EXTENSION
#else
#  undef  SAD_TK_EXTENSION
#endif
#  define SAD_TK_EXTENSION

#ifdef SAD_TK_EXTENSION
extern int bTkCanvTextDebug; 
extern int ctkCanvGetType(const char*);
extern int ctkCanvGetTypes(void (*)(const char*));
extern int ctkCanvCreateItem(Tk_Window, int,
			     int, double*, int, int (*)(const char**));
extern int ctkCanvLastItem(Tk_Window);
extern int ctkCanvFindEnclosed(Tk_Window, const double[4], void (*)(int));
#endif

/* Pseudo Color Table */
static const unsigned char PesudeColorR[256]={
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 2, 3, 4, 4, 5, 6, 7, 8, 10, 11, 12, 13, 15, 16, 17, 19, 20, 22, 23, 25, 27, 29, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48, 51, 53, 55, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 164, 162, 160, 159, 157, 156, 154, 152, 151, 149, 147, 146, 145, 143, 141, 140, 138, 137, 135, 134, 132, 131, 129, 128, 127, 125, 124, 123, 121, 120, 119, 118, 117, 116, 115, 114, 114, 115, 116, 118, 119, 121, 123, 124, 126, 128, 129, 131, 132, 135, 136, 138, 140, 142, 144, 146, 148, 150, 152, 154, 156, 158, 160, 163, 165, 167, 170, 172, 175, 178, 180, 183, 186, 188, 192, 195, 199, 201, 205, 209, 214, 219, 224, 229, 237, 255
};

static const unsigned char PesudeColorG[256]={
0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 6, 6, 6, 6, 7, 7, 7, 7, 8, 8, 8, 9, 9, 9, 9, 10, 10, 10, 11, 11, 11, 11, 12, 12, 12, 13, 13, 14, 14, 14, 15, 15, 15, 16, 16, 17, 17, 17, 18, 19, 19, 20, 22, 23, 25, 26, 28, 29, 31, 33, 34, 36, 38, 39, 41, 42, 44, 46, 48, 50, 52, 53, 55, 57, 59, 61, 63, 65, 67, 69, 71, 72, 74, 76, 78, 80, 83, 84, 86, 89, 91, 93, 95, 97, 99, 101, 104, 106, 108, 110, 112, 114, 117, 119, 121, 123, 126, 128, 130, 133, 135, 137, 140, 142, 144, 146, 148, 151, 153, 155, 158, 160, 162, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239, 240, 240, 240, 240, 240, 240, 240, 240, 241, 242, 243, 245, 247, 255
};

static const unsigned char PesudeColorB[256]={
0, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 56, 55, 55, 54, 54, 53, 53, 52, 51, 51, 50, 49, 48, 47, 46, 45, 45, 44, 43, 41, 40, 39, 38, 37, 36, 34, 33, 32, 30, 29, 28, 26, 25, 23, 22, 20, 19, 19, 20, 20, 21, 21, 22, 22, 22, 23, 24, 24, 24, 25, 25, 26, 27, 27, 28, 28, 29, 29, 30, 30, 31, 32, 32, 33, 34, 34, 35, 35, 36, 36, 37, 38, 38, 39, 40, 41, 41, 42, 43, 44, 44, 45, 46, 47, 47, 48, 49, 50, 50, 51, 52, 53, 53, 54, 55, 56, 57, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 68, 69, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 85, 86, 87, 88, 89, 91, 92, 93, 94, 96, 97, 98, 99, 101, 102, 103, 105, 106, 107, 109, 110, 112, 114, 118, 122, 126, 130, 134, 138, 142, 145, 149, 153, 156, 160, 164, 168, 171, 175, 179, 182, 186, 189, 193, 196, 200, 204, 206, 210, 213, 216, 219, 223, 226, 229, 232, 235, 237, 239, 241, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 255
};

typedef struct {
  Tk_Window rw;
  char **_filehandler;
  int fd_setsize;
  integer8 _itfinterp;
  bool stoploop;
} TkInterState_t;

static TkInterState_t TkInterState = {NULL, NULL, 0, 0, false};
#define itfinterp	TkInterState._itfinterp
#define filehandler	TkInterState._filehandler

#define	FDSET_ALLOC_UNIT	16
static bool alloc_filehandler(int fd) {
  char **new_filehandler;
  int new_setsize, i;

  if(TkInterState.fd_setsize > fd) return true;

  new_setsize = FDSET_ALLOC_UNIT * (fd / FDSET_ALLOC_UNIT + 1);

  new_filehandler = realloc(TkInterState._filehandler,
			    sizeof(char *) * new_setsize);
  if(new_filehandler == NULL) return false;

  for(i = TkInterState.fd_setsize; i < new_setsize; i++)
    new_filehandler[i] = NULL;

  TkInterState._filehandler = new_filehandler;
  TkInterState.fd_setsize = new_setsize;

  return true;
}

static bool check_filehandler(int fd) {
  return TkInterState.fd_setsize > fd;
}

/* Event handler for Tcl/Tk */
static int tclsadeval(ClientData clientData, Tcl_Interp *interp,
		      int argc, CONST char *argv[]) {
  char *buffer;
  int i, offset;
  size_t len;

  if(argc < 1) return TCL_OK;	/* argv[0] shows tcl command */

  switch(argc) {
  case 2:
    tfevalc(argv[1]);
    break;

  default:
    for(i = 1, len = 0; i < argc; i++)
      len += strlen(argv[i]) + 1;
    if(len > 0) {
      buffer = malloc(sizeof(char) * (len + 1));
      if(buffer != NULL) {
	offset = snprintf(buffer, len + 1, "%s[", argv[1]);
	for(i = 2; i < argc && offset < len; i++)
	  offset += snprintf(buffer + offset, (len + 1) - offset,
			     i < argc - 1 ? "%s," : "%s]", argv[i]);
	buffer[len] = '\0';
	if(offset <= len)
	  tfevalc(buffer);
	free(buffer);
      }
    }
    break;
  }

  return TCL_OK;
}

static int tclsadcommand(ClientData clientData, Tcl_Interp *interp,
			 int objc, Tcl_Obj *CONST objv[]) {
  static integer8 itfWidgetCommand = 0;
  integer8 kx;
  integer4 isp0, isp1, irtc;
  logical4 ev;
  int result;

  if(objc < 2) return TCL_OK;	/* objv[1] contains WidgetCommand[] number */

  Tcl_GetIntFromObj(interp, objv[1], &result);

  if(itfWidgetCommand == 0) {
    itfWidgetCommand = ktfsymbolz("WidgetCommand");
  }

  levele += 1;
  isp0 = isp;

  isp += 1;
  ktastk(isp) = ktfsymbol + itfWidgetCommand;
  isp1 = isp;

  isp += 1;
  rtastk(isp) = result;

  tfdeval(isp1, itfWidgetCommand, &kx, 1, false, &ev, &irtc);

  if(irtc > 0 && ierrorprint != 0) {
    tfreseterror_();
  }

  isp = isp0;
  itfdownlevel();

  return TCL_OK;
}

static void sadFileHandleProc(ClientData clientData, int mask) {
  const char *fname = (const char*)clientData;

  if(fname == NULL) {
    fprintf(stderr, "** BadFileHandle trapped by SAD/Tkinter **\n");
    return;
  }

  tfevalc(fname);
}

static int sadBadWindowProc(ClientData clientData, XErrorEvent *errEventPtr) {
  fprintf(stderr, "** BadWindow trapped by SAD/Tkinter **\n");
  fprintf(stderr, "   Major opcode of failed request: %d\n",
	  errEventPtr->request_code);
  return 0;
}

static int sadBadFontProc(ClientData clientData, XErrorEvent *errEventPtr) {
  fprintf(stderr, "** BadFont   trapped by SAD/Tkinter **\n");
  fprintf(stderr, "   Major opcode of failed request: %d\n",
	  errEventPtr->request_code);
  return 0;
}

static int sadBadGCProc(ClientData clientData, XErrorEvent *errEventPtr) {
  fprintf(stderr, "** BadGC   trapped by SAD/Tkinter **\n");
  fprintf(stderr, "   Major opcode of failed request: %d\n",
	  errEventPtr->request_code);
  return 0;
}

/* Tcl interpreter support stuff */
static const char *tftcl(Tcl_Interp *interp, const char *cmd) {
  if(Tcl_EvalEx(interp, cmd, -1, 0) == TCL_ERROR) return NULL;
  Tcl_GetStringResult(interp);

  return interp->result;
}

/* Status control handler */
static void tclUpdate(int mode) {
  Tcl_Interp **interp;

  if(itfinterp == 0) return;

  interp = (Tcl_Interp **)&rlist(itfinterp);

  mode &= TCL_UPDATE_TASK_MASK;
  if(mode == 0) {
    tftcl(*interp, "update");
  } else {
    if(mode & TCL_UPDATE_TASK_FILE)   tftcl(*interp, "update filetasks");
    if(mode & TCL_UPDATE_TASK_TIMER)  tftcl(*interp, "update timertasks");
    if(mode & TCL_UPDATE_TASK_WINDOW) tftcl(*interp, "update windowtasks");
    if(mode & TCL_UPDATE_TASK_IDLE)   tftcl(*interp, "update idletasks");
  }
}

/* CreateFileHandler wrapper */
static void _sadTk_CreateFileHandler(int fd, int mask0,
				     sadTk_FileProc *proc, void *data) {
  int mask = 0;

  /* Translating mask bits */
  if(mask0 & SAD_TK_READABLE)	mask |= TK_READABLE;
  if(mask0 & SAD_TK_WRITABLE)	mask |= TK_WRITABLE;
  if(mask0 & SAD_TK_EXCEPTION)	mask |= TK_EXCEPTION;

  Tk_CreateFileHandler(fd, mask, proc, data);
}

/*
 * Tk_Window object wrapper functions
 * for hiding TkCanvasPointer representation
 */
static int return_Tk_Window(Tk_Window tkWin,
			    integer8 *kx,
			    integer4 *irtc) {
#ifdef ENCODE_POINTER_INTO_STRING
  uintptr_t uptr;
  unsigned char str[sizeof(Tk_Window) + 2];
  int i;
#endif
  intptr_t iptr, iptrR;
  real8 vx;

  if(sizeof(Tk_Window) * CHAR_BIT <= 53) {
    /* Lossless conversion: pointer <=> intptr_t <=> IEEE754 double */
    iptr = (intptr_t)tkWin;
    vx   = iptr;
    *kx = kfromr(vx);
    *irtc = 0;
    return 0;
  }

#ifdef ENCODE_POINTER_INTO_STRING
  else {	/* Encode into String object */
    uptr = (uintptr_t)tkWin;
    str[0] = '^';			/* Check character 1 */
    for(i = 0; i < sizeof(Tk_Window); i++) {
      str[i + 1] = uptr & 0x00ffUL;
      uptr >>= CHAR_BIT;
    }
    str[sizeof(Tk_Window) + 1] = '$';	/* Check character 2 */
    *kx  = ktfstring + ktsalocbl(-1, (char*)str, sizeof(Tk_Window) + 2);
    *irtc = 0;
    return 0;
  }
#endif

  /* Try lossy translation: pointer <=> intptr_t =-> IEEE754 double */
  iptr = (intptr_t)tkWin;
  vx = iptr;			/* Lossy conversion */
  iptrR = vx;
  if(iptr != iptrR) {		/* Check conversion loss */
    fprintf(stderr, "Can't map a Tk_Window object into SAD object\n");
    exit(EX_CONFIG);
  }
  *kx = kfromr(vx);
  *irtc = 0;
  return 0;
}

static Tk_Window refer_Tk_Window(integer4 ispA, integer4 *irtc) {
#ifdef ENCODE_POINTER_INTO_STRING
  integer8 ia;
  uintptr_t uptr;
  unsigned char *str;
  int i;
#endif
  intptr_t iptr;

  if(sizeof(Tk_Window) * CHAR_BIT <= 53) {
    /* Lossless conversion: pointer <=> intptr_t <=> IEEE754 double */
    if(ktfrealq(ktastk(ispA))) {
      iptr = rtastk(ispA);
      return (Tk_Window)iptr;
    }
    *irtc = itfmessage(9, "General::wrongtype",
		       "\"Real Number for #1\"");
    return NULL;
  }

#ifdef ENCODE_POINTER_INTO_STRING
  else {	/* Encode into String object */
    ia = (ktamask & ktastk(ispA));
    if(ktfstringq(ktastk(ispA))
       && ilist(1, ia) == sizeof(Tk_Window) + 2
       && jlist(1, ia + 1) == '^'
       && jlist(ilist(1, ia), ia + 1) == '$') {
      str = (unsigned char *)&jlist(2, ia + 1);
      for(uptr = 0UL, i = sizeof(Tk_Window) - 1; i >= 0; i--) {
	uptr <<= CHAR_BIT;
	uptr |= str[i];
      }
      return (Tk_Window)uptr;
    }
    *irtc = itfmessage(9, "General::wrongtype",
		       "\"Encoded-string for #1\"");
    return NULL;
  }
#endif

  /* Decode lossy translated pointer */
  if(ktfrealq(ktastk(ispA))) {
    iptr = rtastk(ispA);
    return (Tk_Window)iptr;
  }
  *irtc = itfmessage(9, "General::wrongtype",
		     "\"Real Number for #1\"");
  return NULL;
}

/* SADScript function definition of TkInter stuff */
#define tfGetCanvasArgStk_	tfgetcanvasargstk_
extern void tfgetcanvasargstk_(integer8*, integer4*); /* tfpyarg.f */


static int tfTclCreateInterp(integer4 *isp1,
			     integer8 *kx,
			     integer4 *irtc) {
  /*  Tcl_Interp **buffer = (Tcl_Interp **)vx; */
  Tcl_Interp *interp;
  Display *sadtkdisplay;

  if(TK_MAJOR_VERSION == 8 && TK_MINOR_VERSION > 0)
    Tcl_FindExecutable(NULL);
  interp = Tcl_CreateInterp();

  if(Tcl_Init(interp) == TCL_ERROR) {
    fprintf(stderr, "Tcl_Init error: %s\n", interp->result);
    *irtc = itfmessage(9, "General::ioerr", "\"initializing Tcl\"");
    return -1;
  }
  
  if(Tk_Init(interp)  == TCL_ERROR) {
    fprintf(stderr, "Tk_Init error: %s\n",  interp->result);
    *irtc = itfmessage(9, "General::ioerr", "\"initializing Tk\"");
    return -1;
  }
  
  TkInterState.rw = Tk_MainWindow(interp);
  sadtkdisplay = Tk_Display(TkInterState.rw);
  sadXlib_SetNewDisplay(sadtkdisplay);

  Tcl_CreateCommand(interp,	"sadeval",	tclsadeval,	0, NULL);
  Tcl_CreateObjCommand(interp,	"sadcommand",	tclsadcommand,	0, NULL);

  Tk_CreateErrorHandler(sadtkdisplay,
			BadWindow, -1, -1, sadBadWindowProc, NULL);
  Tk_CreateErrorHandler(sadtkdisplay,
			BadFont,   -1, -1, sadBadFontProc,   NULL);
  Tk_CreateErrorHandler(sadtkdisplay,
			BadGC,     -1, -1, sadBadGCProc,     NULL);

  *kx=(integer8) interp;
    /*  *buffer = interp;  store bit image to Real*8 *vx */
  *irtc = 0;
  return 0;
}

static int tfTcl(integer4 *isp1,
		 integer8 *kx, 
		 integer4 *irtc) {
  integer8 ia;
  Tcl_Interp **interp;
  int length;
  char *buffer;

  if(itfinterp == 0)
    itfinterp = ktfsymbolz("Tcl$Interp") - 4;

  interp = (Tcl_Interp **)&rlist(itfinterp);

  if(isp != *isp1 + 1 && isp != *isp1 + 2) {
    *irtc = itfmessage(9, "General::narg", "\"1 or 2\"");
    return -1;
  }

  if(ktfnonstringq(ktastk(*isp1 + 1))){
    if(isp > *isp1 + 1) {
      *irtc = itfmessage(9, "General::wrongtype",
			 "\"Character-string for #1\"");
    } else {
      *irtc = itfmessage(9, "General::wrongtype",
			 "\"Character-string\"");
    }
    return -1;
  }

  ia = (ktamask & ktastk(*isp1 + 1));
#if SAD_REQUIRE_STRING_TERMINATION
  jlist(ilist(1, ia) + 1, ia + 1) = '\0';
#endif

  if(Tcl_EvalEx(*interp, &jlist(1, ia + 1), ilist(1, ia), 0) == TCL_ERROR) {
    tftcl(*interp, "puts stdout $errorInfo");
    *kx = ktfoper + mtfnull;
    *irtc = itfmessage(9, "Tkinter::tclerror", "\"1501\"");
    return 0;
  }

  if(isp == *isp1 + 2) { /* return result string if narg == 2 */
    buffer = Tcl_GetStringFromObj(Tcl_GetObjResult(*interp), &length);
    *kx = ktfstring + ktsalocbl(-1, buffer, length);
  } else {
    *kx = ktfoper + mtfnull;
  }
  *irtc = 0;
  return 0;
}

static int tfTclSetResult(integer4 *isp1,
			  integer8 *kx, 
			  integer4 *irtc) {
  integer8 ia;
  Tcl_Interp **interp;

  if(isp != *isp1 + 1) {
    *irtc = itfmessage(9, "General::narg", "\"1\"");
    return -1;
  }

  if((ktfmask & ktastk(isp)) != ktfstring) {
    *irtc = itfmessage(9, "General::wrongtype", "\"Character-string\"");
    return -1;
  }

  if(itfinterp == 0)
    itfinterp = ktfsymbolz("Tcl$Interp") - 4;

  ia = (ktamask & ktastk(isp));
#if SAD_REQUIRE_STRING_TERMINATION
  jlist(ilist(1, ia) + 1, ia + 1) = '\0';
#endif

  interp = (Tcl_Interp **)&rlist(itfinterp);
  Tcl_SetResult(*interp, &jlist(1,ia+1), TCL_VOLATILE);

  *kx = ktfoper + mtfnull;
  *irtc = 0;
  return 0;
}

static int tfTkEventLoop(integer4 *isp1,
			 integer8 *kx, 
			 integer4 *irtc) {
  if(isp != *isp1 + 1 || ktastk(isp) != ktfoper + mtfnull) {
    *irtc = itfmessage(9, "General::narg","\"0\"");
    return -1;
  }

  TkInterState.stoploop = false;
  while(!TkInterState.stoploop)
    Tk_DoOneEvent(0);

  /* Clear stop loop request for nested event loop */
  TkInterState.stoploop = false;

  *kx = ktfoper + mtfnull;
  *irtc = 0;
  return 0;
}

static int tfTkQuitEventLoop(integer4 *isp1,
			     integer8 *kx,
			     integer4 *irtc) {
  if(isp != *isp1 + 1 || ktastk(isp) != ktfoper + mtfnull) {
    *irtc = itfmessage(9, "General::narg","\"0\"");
    return -1;
  }

  TkInterState.stoploop = true;

  *kx = ktfoper + mtfnull;
  *irtc = 0;
  return 0;
}

static int tfTkCreateFileHandler(integer4 *isp1,
				 integer8 *kx,
				 integer4 *irtc) {
  integer8 ia;
  integer4 lfno;
  int fd;
  void *ptr, *ptr0;

  if(isp != *isp1 + 2) {
    *irtc = itfmessage(9, "General::narg", "\"2\"");
    return -1;
  }

  if((ktrmask & ktastk(*isp1 + 1)) == ktfnr) {
    *irtc = itfmessage(9, "General::wrongtype", "\"Real Number for #1\"");
    return -1;
  }

  if((ktfmask & ktastk(isp)) != ktfstring) {
    *irtc = itfmessage(9, "General::wrongtype", "\"Character-string for #2\"");
    return -1;
  }

  lfno = rtastk(*isp1 + 1);
  ia = (ktamask & ktastk(isp));
#if SAD_REQUIRE_STRING_TERMINATION
  jlist(ilist(1, ia) + 1, ia + 1) = '\0';
#endif

  fd = getfd_(&lfno);

  if(!alloc_filehandler(fd) || (ptr = strdup(&jlist(1, ia + 1))) == NULL) {
    fprintf(stderr, "TkCreateFileHandler: malloc failed!\n");
    *irtc = itfmessage(9, "General::memoryfull", " ");
    return -1;
  }

  ptr0 = filehandler[fd]; filehandler[fd] = NULL;
  Tk_CreateFileHandler(fd, TK_READABLE, sadFileHandleProc, ptr);
  filehandler[fd] = ptr; if(ptr0 != NULL) free(ptr0);

  *kx = ktfoper + mtfnull;
  *irtc = 0;
  return 0;
}

static int tfTkDeleteFileHandler(integer4 *isp1,
				 integer8 *kx,
				 integer4 *irtc) {
  integer4 lfno;
  int fd;
  void *ptr0;

  if(isp != *isp1 + 1) {
    *irtc = itfmessage(9, "General::narg", "\"1\"");
    return -1;
  }

  if((ktrmask & ktastk(isp)) == ktfnr) {
    *irtc = itfmessage(9, "General::wrongtype", "\"Real Number\"");
    return -1;
  }

  lfno = rtastk( isp);
  fd = getfd_(&lfno);

  if(check_filehandler(fd)) {
    Tk_DeleteFileHandler(fd);
    ptr0 = filehandler[fd]; filehandler[fd] = NULL;
    if(ptr0 != NULL) free(ptr0);
  }

  *kx = ktfoper + mtfnull;
  *irtc = 0;
  return 0;
}

static int tfTkDefineBitmap(integer4 *isp1,
			    integer8 *kx,
			    integer4 *irtc) {
  integer8 ia,is,ki;
  int width, height, length;
  int i, status;
  char *buffer;
  Tcl_Interp **interp;

  if(itfinterp == 0) return -1;

  if(isp != *isp1 + 4) {
    *irtc = itfmessage(9, "General::narg", "\"4\"");
    return -1;
  }

  if((ktfmask & ktastk(*isp1 + 1)) != ktfstring) {
    *irtc = itfmessage(9, "General::wrongtype", "\"Character-string for #1\"");
    return -1;
  }

  if((ktfmask & ktastk(*isp1 + 2)) != ktflist) {
    *irtc = itfmessage(9, "General::wrongtype", "\"List of Reals for #2\"");
    return -1;
  }

  if((ktrmask & ktastk(*isp1 + 3)) == ktfnr) {
    *irtc = itfmessage(9, "General::wrongtype", "\"Real Number for #3\"");
    return -1;
  }

  if((ktrmask & ktastk(isp)) == ktfnr) {
    *irtc = itfmessage(9, "General::wrongtype", "\"Real Number for #4\"");
    return -1;
  }

  interp = (Tcl_Interp **)&rlist(itfinterp);

  is = (ktamask & ktastk(*isp1 + 1));
#if SAD_REQUIRE_STRING_TERMINATION
  jlist(ilist(1, is) + 1, is + 1) = '\0';
#endif
  ia = (ktamask & ktastk(*isp1 + 2));
  width  = rtastk(*isp1 + 3);
  height = rtastk( isp);

  length = ilist(2, ia - 1);
  buffer = malloc(length);
  if(buffer == NULL) {
    fprintf(stderr, "TkDefineBitmap: malloc failed!\n");
    *irtc = itfmessage(9, "General::memoryfull", " ");
    return -1;
  }

  for(i = 0; i < length; i++) {
    ki=klist(ia + i + 1);
    if((ktrmask & ki) == ktfnr) {
      free(buffer);
      *irtc = itfmessage(9, "General::wrongtype", "\"List of Reals for #2\"");
      return -1;
    }
    buffer[i] = rfromk(ki);
  }

  /* Don't free ``buffer'' if TK_DefineBitmap is succeeded */
  status = Tk_DefineBitmap(*interp, Tk_GetUid(&jlist(1, is + 1)),
			   buffer, width, height);

  if(status == TCL_OK) {
    *kx = kfromr(r_true);
  } else {
    *kx = 0;
    free(buffer);
  }
  *irtc = 0;
  return 0;
}

static int tfTkPhotoPutBlock(integer4 *isp1,
			     integer8 *kx,
			     integer4 *irtc) {
  integer8 is, ia;
  real8 *source;
  int type, x, y, width, height;
  int xzoom = 1, yzoom = 1, xdiv = 1, ydiv = 1;
  int npixel, depth, array_length;
  int i, j, index;
  Tcl_Interp **interp;
  Tk_PhotoHandle handle;
  Tk_PhotoImageBlock block;
  unsigned char *pixelBuf;

  if(isp != *isp1 + 7 && isp != *isp1 + 9 && isp != *isp1 + 11) {
    *irtc = itfmessage(9, "General::narg", "\"7, 9 or 11\"");
    return -1;
  }

  if((ktfmask & ktastk(*isp1 + 1)) != ktfstring) {
    *irtc = itfmessage(9, "General::wrongtype", "\"Character-string for #1\"");
    return -1;
  }

  if((ktfmask & ktastk(*isp1 + 2)) != ktflist) {
    *irtc = itfmessage(9, "General::wrongtype", "\"List of Reals for #2\"");
    return -1;
  }

  if((ktrmask & ktastk(*isp1 + 3)) == ktfnr) {
    *irtc = itfmessage(9, "General::wrongtype", "\"Real Number for #3\"");
    return -1;
  }

  if((ktrmask & ktastk(*isp1 + 4)) == ktfnr) {
    *irtc = itfmessage(9, "General::wrongtype", "\"Real Number for #4\"");
    return -1;
  }

  if((ktrmask & ktastk(*isp1 + 5)) == ktfnr) {
    *irtc = itfmessage(9, "General::wrongtype", "\"Real Number for #5\"");
    return -1;
  }

  if((ktrmask & ktastk(*isp1 + 6)) == ktfnr) {
    *irtc = itfmessage(9, "General::wrongtype",
		       "\"Positive Real Number for #6\"");
    return -1;
  }

  if((ktrmask & ktastk(*isp1 + 7)) == ktfnr) {
    *irtc = itfmessage(9, "General::wrongtype",
		       "\"Positive Real Number for #7\"");
    return -1;
  }

  if(isp > *isp1 + 7) {
    if((ktrmask & ktastk(*isp1 + 8)) == ktfnr) {
      *irtc = itfmessage(9, "General::wrongtype", "\"Real Number for #8\"");
      return -1;
    }

    if((ktrmask & ktastk(*isp1 + 9)) == ktfnr) {
      *irtc = itfmessage(9, "General::wrongtype", "\"Real Number for #9\"");
      return -1;
    }
  }

  if(isp > *isp1 + 9) {
    if((ktrmask & ktastk(*isp1 + 10)) == ktfnr) {
      *irtc = itfmessage(9, "General::wrongtype", "\"Real Number for #10\"");
      return -1;
    }

    if((ktrmask & ktastk(isp)) == ktfnr) {
      *irtc = itfmessage(9, "General::wrongtype", "\"Real Number for #11\"");
      return -1;
    }
  }

  if(itfinterp == 0) return -1;

  interp = (Tcl_Interp **)&rlist(itfinterp);

  is = (ktamask & ktastk(*isp1 + 1));
#if SAD_REQUIRE_STRING_TERMINATION
  jlist(ilist(1, is) + 1, is + 1) = '\0';
#endif
  ia = (ktamask & ktastk(*isp1 + 2));
  type   = rtastk(*isp1 + 3);
  x      = rtastk(*isp1 + 4);
  y      = rtastk(*isp1 + 5);
  width  = rtastk(*isp1 + 6);
  height = rtastk(*isp1 + 7);

  if(width <= 0) {
    *irtc = itfmessage(9, "General::wrongtype",
		       "\"Positive Real Number for #6\"");
    return -1;
  }

  if(height <= 0) {
    *irtc = itfmessage(9, "General::wrongtype",
		       "\"Positive Real Number for #7\"");
    return -1;
  }

  npixel = width * height;
  switch(type) {
  case 1:
    depth = 3;
    array_length = npixel;
    break;

  case 2:
  case 3:
    depth = 3;
    array_length = depth * npixel;
    break;

  default:
    depth = 1;
    array_length = npixel;
  }

  if((lnonreallist & ilist(2, ia - 3)) != 0) {
    *irtc = itfmessage(9, "General::wrongtype", "\"List of Reals for #2\"");
    return -1;
  }

  if(ilist(2, ia - 1) < array_length) {
    *irtc = itfmessage(9, "General::wrongtype",
		       "\"Array is too short for #2\"");
    return -1;
  }

  source = &rlist(ia + 1);

  if(isp > *isp1 + 7) {
    xzoom = rtastk(*isp1 +  8); if(xzoom < 1) xzoom = 1;
    yzoom = rtastk(*isp1 +  9); if(yzoom < 1) yzoom = 1;
  }

  if(isp > *isp1 + 9) {
    xdiv  = rtastk(*isp1 + 10); if(xdiv  < 1) xdiv  = 1;
    ydiv  = rtastk(*isp1 + 11); if(ydiv  < 1) xdiv  = 1;
  }

  /* Compose pixel buffer */
  if((pixelBuf = malloc(depth * npixel)) == NULL) {
    fprintf(stderr, "Memory allocation failed in TkPhotoPutBlock\n");
    *irtc = itfmessage(9, "General::memoryfull", " ");
    return -1;
  }

  switch(type) {
  case 1: /* Pseudo color */
    block.offset[0] = 0;
    block.offset[1] = 1 * npixel;
    block.offset[2] = 2 * npixel;
    block.offset[3] = 0;
    for(i = 0; i < npixel; i++) {
      index = source[i];
      if(index <   0) index = 0;
      if(index > 255) index = 255;
      pixelBuf[i             ] = PesudeColorR[index];
      pixelBuf[i +     npixel] = PesudeColorG[index];
      pixelBuf[i + 2 * npixel] = PesudeColorB[index];
    }
    break;

  case 2: /* Color format 1 (R...G...B...) */
    block.offset[0] = 0;
    block.offset[1] = 1 * npixel;
    block.offset[2] = 2 * npixel;
    block.offset[3] = 0;
    for(i = 0; i < depth * npixel; i++) pixelBuf[i] = source[i];
    break;

  case 3: /* Color format 2 (RGBRGB......) */
    block.offset[0] = 0;
    block.offset[1] = 1 * npixel;
    block.offset[2] = 2 * npixel;
    block.offset[3] = 0;
    for(i = j = 0; i < npixel; i++, j+= 3) {
      pixelBuf[i             ] = source[j    ];
      pixelBuf[i +     npixel] = source[j + 1];
      pixelBuf[i + 2 * npixel] = source[j + 2];
    }
    break;

  default:
    block.offset[0] = 0;
    block.offset[1] = 0;
    block.offset[2] = 0;
    block.offset[3] = 0;
    for(i = 0; i < depth * npixel; i++) pixelBuf[i] = source[i];
  }

  handle = Tk_FindPhoto(*interp, &jlist(1, is + 1));
  block.pixelPtr  = pixelBuf;
  block.width     = width;
  block.height    = height;
  block.pitch     = width;
  block.pixelSize = 1;

  width  = width  * xzoom / xdiv; if(width  < 1) width  = 1;
  height = height * yzoom / ydiv; if(height < 1) height = 1;

#if TK_MINOR_VERSION >= 5
#ifdef USE_COMPOSITELESS_PHOTO_PUT_BLOCK
  Tk_PhotoPutZoomedBlock(handle, &block, x, y, width, height,
			 xzoom, yzoom, xdiv, ydiv);
#else /* !USE_COMPOSITELESS_PHOTO_PUT_BLOCK */
#ifdef USE_PANIC_ON_PHOTO_ALLOC_FAILURE
  Tk_PhotoPutZoomedBlock(handle, &block, x, y, width, height,
			 xzoom, yzoom, xdiv, ydiv,
			 TK_PHOTO_COMPOSITE_SET);
#else
  Tk_PhotoPutZoomedBlock(*interp, handle, &block, x, y, width, height,
			 xzoom, yzoom, xdiv, ydiv,
			 TK_PHOTO_COMPOSITE_SET);
#endif
#endif
#else /* TK_MINOR_VERSION < 5 */
#ifdef USE_COMPOSITELESS_PHOTO_PUT_BLOCK
  Tk_PhotoPutZoomedBlock(handle, &block, x, y, width, height,
			 xzoom, yzoom, xdiv, ydiv);
#else /* !USE_COMPOSITELESS_PHOTO_PUT_BLOCK */
  Tk_PhotoPutZoomedBlock(handle, &block, x, y, width, height,
			 xzoom, yzoom, xdiv, ydiv,
			 TK_PHOTO_COMPOSITE_SET);
#endif
#endif
  free(pixelBuf);

  *kx = ktfoper + mtfnull;
  *irtc = 0;
  return 0;
}

static void ctkCanvGetTypeCB(const char *type) {
  isp += 1;
  ktastk(isp) = ktfstring + ktsalocb(-1, type);
}

static int tfTkCanvGetTypes(integer4 *isp1,
			    integer8 *kx,
			    integer4 *irtc) {
  integer4 isp0;

  if(isp != *isp1 + 1 || ktastk(isp) != ktfoper + mtfnull) {
    *irtc = itfmessage(9, "General::narg","\"0\"");
    return -1;
  }

  isp0 = isp;

#ifdef SAD_TK_EXTENSION
  ctkCanvGetTypes(ctkCanvGetTypeCB);
#endif

  if(isp == isp0) {
    isp += 1;
    ktastk(isp) = ktfoper + mtfnull;
  }

  *kx = ktflist + ktfmakelist(isp0);
  isp = isp0;
  *irtc = 0;
  return 0;
}

static int ctkCanvCreateItemCB(const char **buf) {
  integer8 ia;

  if(ktfnonstringq(ktastk(isp))) {
    *buf = NULL;
    return -1;			/* Notify no string in stack buffer */
  }

  ia = (ktamask & ktastk(isp));
  isp -= 1;			/* Step stack pointer backword!! */
  *buf = &jlist(1, ia + 1);	/* Save string pointer into *buf */
  return ilist(1, ia);		/* Return string length */
}

static int tfTkCanvCreateItem(integer4 *isp1,
			      integer8 *kx,
			      integer4 *irtc) {
  integer8 ia, kt;
  real8 vx;
  integer4 isp0, isp2, isp3;
  Tk_Window tkWin;
  int type, parameters, options;
  double *parameter;
  int i, j, id;

  if(isp < *isp1 + 3) {
    *irtc = itfmessage(9, "General::narg", "\"3 or more\"");
    return -1;
  }

  if((tkWin = refer_Tk_Window(*isp1 + 1, irtc)) == NULL)
    return -1;

  if(ktfnonrealq(ktastk(*isp1 + 2))) {
    *irtc = itfmessage(9, "General::wrongtype", "\"Real Number for #2\"");
    return -1;
  }

  type = rtastk(*isp1 + 2);

  isp0 = isp;	/* Backup current stack pointer */

  /* Extract Real parameter list block from *isp1 + 3 */
  /* isp2 is setuped to point next parameter block */
  if(ktfrealq(ktastk(*isp1 + 3))){
    for(isp2 = *isp1 + 4; isp2 <= isp; isp2++)
      if(ktfnonrealq(ktastk(isp2))) break;
    parameters = isp2 - (*isp1 + 3);
    parameter  = &rtastk(*isp1 + 3);
  }
  else if(ktflistq(ktastk(*isp1 + 3))){
    isp2 = *isp1 + 4;	/* Next block begins 4th argument */
    ia  = (ktamask & ktastk(*isp1 + 3));
    if((lnonreallist & ilist(2, ia - 3)) == 0) {
      parameters = ilist(2, ia - 1);
    } else {
      for(i = 0, j = 1; i < ilist(2, ia - 1); i++, j++)
	if(ktfnonrealq(klist(ia + j))) break;
      parameters = i;
    }
    parameter = &rlist(ia + 1);
  }
  else{
    *irtc = itfmessage(9, "General::wrongtype",
		       "\"Reals or List of Reals for parameter block\"");
    return -1;
  }

  /* Extract String parameter list block from isp2 */
  options = 0;
  if(isp >= isp2) {
    kt = (ktfmask & ktastk(isp2));
    if(kt == ktfstring){
      for(isp3 = isp2 + 1; isp3 <= isp; isp3++)
	if((ktfmask & ktastk(isp3)) != ktfstring) break;
      options = isp3 - isp2;

      /* Push non-string guard object[Null] */
      isp += 1;
      ktastk(isp) = ktfoper +  mtfnull;
          
      /* Push option string objects into stack backword */
      /* Object range in stack: [isp2, isp2 + options - 1 == isp3 - 1] */
      for(i = options - 1; i >= 0; i--) {
        isp += 1;
        ktastk(isp) = ktfstring + (ktamask & ktastk(isp2 + i));
      }
    }
    else if (kt == ktflist){
      ia  = (ktamask & ktastk(isp2));
      if((lnonreallist & ilist(2, ia - 3)) != 0) {
        
        for(i = 0, j = 1; i < ilist(2, ia - 1); i++, j++)
          if((ktfmask & klist(ia + j)) != ktfstring) break;

        options = i; if(options >= 1) {

      /* Push non-string guard object[Null] */
          isp += 1;
          ktastk(isp) = ktfoper + mtfnull;
        
      /* Push option string objects into stack backword */
      /* Object range in list: [iad + 1, iad + options] */
          for(i = options; i > 0; i--) {
            isp += 1;
            ktastk(isp) = ktfstring + (ktamask & klist(ia + i));
          }
        }
      }
      else{
        *irtc = itfmessage(9, "General::wrongtype",
			 "\"Strings or List of Strings for option block\"");
        return -1;
      }
    }
  }

*kx = ktfoper + mtfnull;

#ifdef SAD_TK_EXTENSION
  id = ctkCanvCreateItem(tkWin, type, parameters, parameter,
			 options, ctkCanvCreateItemCB);
  if(id >= 0) {
    vx = id;
    *kx = kfromr(vx);
  }
#endif

  *irtc = 0;
  isp = isp0;	/* Rewind stack pointer */
  return 0;
}

static int tfTkCanvCreateItemDirect(integer4 *isp1,
				    integer8 *kx,
				    integer4 *irtc) {
  integer8 ia;
  integer4 isp0, ispi;
  int status;

  if(isp < *isp1 + 2) {
    *irtc = itfmessage(9, "General::narg", "\"2 or more\"");
    return -1;
  }

  isp0 = isp;

  /* Copy first and second arguments */
  isp += 1;
  ktastk(isp) = ktastk(*isp1 + 1);

  isp += 1;
  ktastk(isp) = ktastk(*isp1 + 2);

  /* Scan following arguments */
  for(ispi = *isp1 + 3; ispi <= isp0; ispi++) {
    if((ktrmask & ktastk(ispi)) != ktfnr){
      isp += 1;
      rtastk(isp) = rtastk(ispi);
    }
    else if((ktfmask & ktastk(ispi)) == ktflist){
      ia = (ktamask & ktastk(ispi));
      if((ktfmask & klist(ia)) != ktfoper) continue; /* list has Non-operator head */

      tfGetCanvasArgStk_(&ia, irtc);
      if(*irtc != 0) {
	isp = isp0;
	return -1;
      }
      if(isp == isp0 + 2) { /* zero length list expansion?? */
	isp = isp0;
	return 0;
      }
    }
    else{
      /* Skip argument */
      break;
    }
  }
    status = tfTkCanvCreateItem(&isp0, kx, irtc);
  isp = isp0;
  return status;
}

static int tfTkCanvPtr(integer4 *isp1,
		       integer8 *kx,
		       integer4 *irtc) {
  integer8 ia;
  Tcl_Interp **interp;
  Tk_Window tkWin;

  if(itfinterp == 0) return -1;

  if(isp != *isp1 + 1) {
    *irtc = itfmessage(9, "General::narg", "\"1\"");
    return -1;
  }

  if((ktfmask & ktastk(isp)) != ktfstring) {
    *irtc = itfmessage(9, "General::wrongtype", "\"Character-string\"");
    return -1;
  }

  ia = (ktamask & ktastk(isp));
#if SAD_REQUIRE_STRING_TERMINATION
  jlist(ilist(1, ia) + 1, ia + 1) = '\0';
#endif

  interp = (Tcl_Interp **)&rlist(itfinterp);
  tkWin = Tk_NameToWindow(*interp, &jlist(1,ia+1), TkInterState.rw);

  return return_Tk_Window(tkWin, kx, irtc);
}

static int tfTkCanvLastItem(integer4 *isp1,
			    integer8 *kx,
			    integer4 *irtc) {
  Tk_Window tkWin;
  real8 vx;

  if(isp != *isp1 + 1) {
    *irtc = itfmessage(9, "General::narg", "\"1\"");
    return -1;
  }

  if((tkWin = refer_Tk_Window(*isp1 + 1, irtc)) == NULL)
    return -1;

#ifdef SAD_TK_EXTENSION
  vx  = ctkCanvLastItem(tkWin);
  *kx = kfromr(vx);
#else
  *kx = 0;
#endif

  *irtc = 0;
  return 0;
}

static void ctkCanvFindEnclosedCB(int id) {
  isp += 1;
  rtastk(isp) = id;
}

static int tfTkCanvFindEnclosed(integer4 *isp1,
				integer8 *kx,
				integer4 *irtc) {
  integer8 ia;
  integer4 isp0;
  Tk_Window tkWin;

  if(isp != *isp1 + 2) {
    *irtc = itfmessage(9, "General::narg", "\"2\"");
    return -1;
  }

  if((tkWin = refer_Tk_Window(*isp1 + 1, irtc)) == NULL)
    return -1;

  ia = (ktamask & ktastk(*isp1 + 2));
  if(!((ktfmask & ktastk(*isp1 + 2)) == ktflist
       && klist(ia) == ktfoper + mtflist
       && (lnonreallist & ilist(2,ia-3)) != 0 && ilist(2, ia - 1) == 4)) {
    *irtc = itfmessage(9, "General::wrongtype",
		       "\"List of four Reals for #2\"");
    return -1;
  }

  isp0 = isp;

#ifdef SAD_TK_EXTENSION
  ctkCanvFindEnclosed(tkWin, &rlist(ia + 1),
		      ctkCanvFindEnclosedCB);
#endif

  if(isp == isp0) {
    isp += 1;
    ktastk(isp) = ktfoper + mtfnull;
  }

  *kx = ktflist + ktfmakelist(isp0);
  isp = isp0;
  *irtc = 0;
  return 0;
}

static int tfTkCanvTextDebug(integer4 *isp1,
			     integer8 *kx,
			     integer4 *irtc) {
  int mode;

  if(isp != *isp1 + 1) {
    *irtc = itfmessage(9, "General::narg", "\"1\"");
    return -1;
  }

  if((ktrmask & ktastk(isp)) == ktfnr) {
    *irtc = itfmessage(9, "General::wrongtype", "\"Real Number\"");
    return -1;
  }

  mode = rtastk( isp);

#if (TK_MAJOR_VERSION>=8) && (TK_MINOR_VERSION>=1)
  bTkCanvTextDebug = mode;
#endif

  *kx = ktfoper + mtfnull;
  *irtc = 0;
  return 0;
}

/* SADScript function registration of TkInter stuff */
/* tfpyarg.f */
extern DECLARE_SAD_FUNC8(tftclarg_);
extern DECLARE_SAD_FUNC8(tftclarg1_);
extern DECLARE_SAD_FUNC8(tfuniconv_);

/* tfcanvasclip.f */
#define CanvasClip		tfcanvasclip_
#define CanvasSymbol		tfcanvassymbol_
#define CanvasSymbolDirect	tfcanvassymboldirect_
#define Canvas3DClipTriangle	tfcanvas3dcliptriangle_
#define Canvas3DLightTriangle	tfcanvas3dlighttriangle_
#define Canvas3DProjection	tfcanvas3dprojection_
extern DECLARE_SAD_FUNC8(tfcanvasclip_);
extern DECLARE_SAD_FUNC8(tfcanvassymbol_);
extern DECLARE_SAD_FUNC8(tfcanvassymboldirect_);
extern DECLARE_SAD_FUNC8(tfcanvas3dcliptriangle_);
extern DECLARE_SAD_FUNC8(tfcanvas3dlighttriangle_);
extern DECLARE_SAD_FUNC8(tfcanvas3dprojection_);

/* Namespace conversion SIM for tfcanvassymboldirect_ */
int tftkcanvcreateitem_(integer4 *isp1,
			integer8 *kx,
			integer4 *irtc) {
  return tfTkCanvCreateItem(isp1, kx, irtc);
}

#define REG	dlfunaloc
#define REG8	dlfunaloc8
#ifdef WITH_EXTENSION_MODULE
int dldeffun_(void) {
#else
int sadDefFunc_TkInter(void) {
#endif /* WITH_EXTENSION_MODULE */
  /* Linking TclUpdate handler */
  sadtk_hooks.TclUpdate = tclUpdate;

  /* Linking Create/DeleteFileHandler */
  sadtk_hooks.CreateFileHandler = _sadTk_CreateFileHandler;
  sadtk_hooks.DeleteFileHandler = Tk_DeleteFileHandler;

  /* Tcl interpreter stuff */
  REG8("TclCreateInterp",	tfTclCreateInterp,	 0, NULL, NULL, 0);
  REG8("Tcl",			tfTcl,			 1, NULL, NULL, 0);
  REG8("TclArg",			tftclarg_,		 1, NULL, NULL, 0);
  REG8("TclArg1",		tftclarg1_,		 1, NULL, NULL, 0);
  REG8("TclUnicode",		tfuniconv_,		 1, NULL, NULL, 0);
  REG8("TclSetResult",		tfTclSetResult,		 1, NULL, NULL, 0);

  /* Tk stuff*/
  REG8("TkEventLoop",		tfTkEventLoop,		 0, NULL, NULL, 0);
  REG8("TkQuitEventLoop",	tfTkQuitEventLoop,	 0, NULL, NULL, 0);
  REG8("TkCreateFileHandler",	tfTkCreateFileHandler,	 2, NULL, NULL, 0);
  REG8("TkDeleteFileHandler",	tfTkDeleteFileHandler,	 1, NULL, NULL, 0);
  REG8("TkDefineBitmap",		tfTkDefineBitmap,	 4, NULL, NULL, 0);
  REG8("TkPhotoPutBlock$",	tfTkPhotoPutBlock,	11, NULL, NULL, 0);

  /* TkCanvas stuff */
  REG8("TkCanvTextDebug",	tfTkCanvTextDebug,	 1, NULL, NULL, 0);

  REG8("TkCanvasGetTypes",	tfTkCanvGetTypes,	 0, NULL, NULL, 0);
  REG8("TkCanvasPointer",	tfTkCanvPtr,		 1, NULL, NULL, 0);
  REG8("TkCanvasLastItem",	tfTkCanvLastItem,	 1, NULL, NULL, 0);
  REG8("TkCanvasFindEnclosed",	tfTkCanvFindEnclosed,	 2, NULL, NULL, 0);

  REG8("TkCanvasCreateItem",	tfTkCanvCreateItem,	 3, NULL, NULL, 0);
  REG8("TkCanvasCreateItemDirect",
      tfTkCanvCreateItemDirect,				 2, NULL, NULL, 0);

  REG8("CanvasClipLine",		CanvasClip,		 2, NULL, NULL, 0);
  REG8("CanvasSymbol",		CanvasSymbol,		 1, NULL, NULL, 0);
  REG8("CanvasSymbolDirect",	CanvasSymbolDirect,	 1, NULL, NULL, 0);
  REG8("Canvas3DClipTriangle",	Canvas3DClipTriangle,	 1, NULL, NULL, 0);
  REG8("Canvas3DLightTriangle",	Canvas3DLightTriangle,	 2, NULL, NULL, 0);
  REG8("Canvas3DProjection",	Canvas3DProjection,	 4, NULL, NULL, 0);

  return 0;
}

/* End of File */
