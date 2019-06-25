#include <sim/sad_api.h>
#include <sim/TFCODE.h>
#include <sim/TFCBK.h>
#include <sim/TFSTK.h>

#include "Python.h"

#include <stdbool.h>
#include <stdio.h>

#include <tcl.h>
#include <tk.h>

extern char *Py_GetProgramName ();

/**** Tkapp Object Declaration ****/

staticforward PyTypeObject Tkapp_Type;

typedef struct
  {
    PyObject_HEAD
    Tcl_Interp *interp;
#ifdef NEED_TKCREATEMAINWINDOW
    Tk_Window tkwin;
#endif
  }
TkappObject;

#define Tkapp_Check(v) ((v)->ob_type == &Tkapp_Type)
#ifdef NEED_TKCREATEMAINWINDOW
#define Tkapp_Tkwin(v)  (((TkappObject *) (v))->tkwin)
#endif
#define Tkapp_Interp(v) (((TkappObject *) (v))->interp)

/* function prototypes */
static PyObject *sad_eval(PyObject *self, PyObject *args);
static PyObject *sad_evalcommand(PyObject *self, PyObject *args);
static PyObject *sad_interp(PyObject *self, PyObject *args);

/*  Method Table  */
static PyMethodDef SadMethods[]={
   {"eval",sad_eval,1}, /* name in Python,function, flag(1 always) */
   {"evalcommand",sad_evalcommand,1}, /* name in Python,function, flag(1 always) */
   {"interp",sad_interp,1}, /* name in Python,function, flag(1 always) */
   {NULL, NULL, 0}
 };

/* initialization routine */
/* module name and a initialization routine are registered in the inittab record
in PyConfig.c */
void initsad(){
  (void) Py_InitModule("_sad", SadMethods);
}

static PyObject *
sad_eval(PyObject *self, PyObject *args){
  integer4 nc;
  static char *command;
  int len, status = 0;

  if(!PyArg_ParseTuple(args, "z#", &command, &len))
    return NULL;

  /* evaluate command in SAD/FFS interpreter */
  nc = len;	/* command MAY contains nul character */
  tfevalc_(command, &nc, nc);

  return Py_BuildValue("i", status);
} 

static PyObject *
sad_evalcommand(PyObject *self, PyObject *args){
  static integer8 itfWidgetCommand = 0;
  integer8 kx;
  integer4 isp0, isp1, irtc;
  real8 vx;
  logical ev;
  static char *command;
  int len, status;

  if(!PyArg_ParseTuple(args, "z#", &command, &len))
    return NULL;

  if(itfWidgetCommand == 0) {
    itfWidgetCommand = ktfsymbolz("WidgetCommand");
  }

  /* evaluate WidgetCommand[command-id] in SAD/FFS interpreter */
  /*
   * >> tfevalwidgetcommand_(command, &nc, &kx, &irtc); 
   * Final version of this tfevalwidgetcommand_() proto-type came
   * came from CVS src/tfpyarg.f revision 1.33.
   * Reimplemented by using src/tfTkInter_.c code fragment.
   */

  levele += 1;
  isp0 = isp;

  isp += 1;
  ktastk(isp) = ktfsymbol + itfWidgetCommand;
  isp1 = isp;

  isp += 1;
  ktastk(isp) = ktfstring + ktsalocbl(-1, command, len);

  tfdeval(isp1, itfWidgetCommand, &kx, 1,
	  ktfsymbol + itfWidgetCommand, false, &ev, &irtc);

  if(irtc > 0 && ierrorprint != 0) tfreseterror_();

  isp = isp0;
  status = itfdownlevel();

  status = irtc;

  return Py_BuildValue("i", status);
} 


static PyObject *
sad_interp(PyObject *self,PyObject *args)
{
  PyObject *obj;
  if(!PyArg_ParseTuple(args, "O", &obj))
    return NULL;

  return Py_BuildValue("i",Tkapp_Interp(obj));
}
