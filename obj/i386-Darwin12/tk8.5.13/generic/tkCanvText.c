/*
 * tkCanvText.c --
 *
 *	This file implements text items for canvas widgets.
 *
 * Copyright (c) 1991-1994 The Regents of the University of California.
 * Copyright (c) 1994-1997 Sun Microsystems, Inc.
 *
 * See the file "license.terms" for information on usage and redistribution of
 * this file, and for a DISCLAIMER OF ALL WARRANTIES.
 */

#include <stdio.h>
#include "tkInt.h"
#include "tkCanvas.h"
#include "default.h"

EXTERN Tk_TextLayout	Tk_KGComputeTextLayout _ANSI_ARGS_((Tk_Font font,
			    CONST char *string, int numChars, int wrapLength,
			    Tk_Justify justify, int flags, int *widthPtr,
			    int *heightPtr, Tk_Font, Tk_Font, Tk_Font,
			    GC, GC, GC, GC, int));
EXTERN void             Tk_KGDrawTextLayout _ANSI_ARGS_((Display *display,
                            Drawable drawable, GC gc, Tk_TextLayout layout,
                            int  x, int y, int firstChar, int lastChar));


/*
 * The structure below defines the record for each text item.
 */

typedef struct TextItem {
    Tk_Item header;		/* Generic stuff that's the same for all
				 * types. MUST BE FIRST IN STRUCTURE. */
    Tk_CanvasTextInfo *textInfoPtr;
				/* Pointer to a structure containing
				 * information about the selection and
				 * insertion cursor. The structure is owned by
				 * (and shared with) the generic canvas
				 * code. */
    /*
     * Fields that are set by widget commands other than "configure".
     */

    double x, y;		/* Positioning point for text. */
    int insertPos;		/* Character index of character just before
				 * which the insertion cursor is displayed. */

    /*
     * Configuration settings that are updated by Tk_ConfigureWidget.
     */

    Tk_Anchor anchor;		/* Where to anchor text relative to (x,y). */
    Tk_TSOffset tsoffset;
    XColor *color;		/* Color for text. */
    XColor *activeColor;	/* Color for text. */
    XColor *disabledColor;	/* Color for text. */
    Tk_Font tkfont;		/* Font for drawing text. */
    Tk_Font tkaltfont;	/* Information about alt text font, or NULL. */
    Tk_Font tkscriptfont;	/* Information about script font, or NULL. */
    Tk_Font tkaltscriptfont;	/* Information about alt script font, or NULL. */
    Tk_Justify justify;		/* Justification mode for text. */
    Pixmap stipple;		/* Stipple bitmap for text, or None. */
    Pixmap activeStipple;	/* Stipple bitmap for text, or None. */
    Pixmap disabledStipple;	/* Stipple bitmap for text, or None. */
    char *text;			/* Text for item (malloc-ed). */
    int width;			/* Width of lines for word-wrap, pixels. Zero
				 * means no word-wrap. */
    int underline;		/* Index of character to put underline beneath
				 * or -1 for no underlining. */

    /*
     * Fields whose values are derived from the current values of the
     * configuration settings above.
     */

    int numChars;		/* Length of text in characters. */
    int numBytes;		/* Length of text in bytes. */
    Tk_TextLayout textLayout;	/* Cached text layout information. */
    int leftEdge;		/* Pixel location of the left edge of the text
				 * item; where the left border of the text
				 * layout is drawn. */
    int rightEdge;		/* Pixel just to right of right edge of area
				 * of text item. Used for selecting up to end
				 * of line. */
    GC gc;			/* Graphics context for drawing text. */
    GC altGc;			/* Graphics context for drawing alttext. */
    GC scriptGc;		/* Graphics context for drawing scripttext. */
    GC altScriptGc;		/* Graphics context for drawing altscripttext. */
    GC selTextGC;		/* Graphics context for selected text. */
    GC cursorOffGC;		/* If not None, this gives a graphics context
				 * to use to draw the insertion cursor when
				 * it's off. Used if the selection and
				 * insertion cursor colors are the same. */
    GC copyGC;			/* Used for copying information from an
				 * off-screen pixmap to the screen. */
    int rotation, dx, dy;
    int pixoffx, pixoffy;
    int usebg;
    XColor *bgcolor;
    GC bgGC;
    int bsmplw,brmdes;
} TextItem;

/*
 * Information used for parsing configuration specs:
 */

static Tk_CustomOption stateOption = {
    (Tk_OptionParseProc *) TkStateParseProc,
    TkStatePrintProc, (ClientData) 2
};
static Tk_CustomOption tagsOption = {
    (Tk_OptionParseProc *) Tk_CanvasTagsParseProc,
    Tk_CanvasTagsPrintProc, (ClientData) NULL
};
static Tk_CustomOption offsetOption = {
    (Tk_OptionParseProc *) TkOffsetParseProc,
    TkOffsetPrintProc, (ClientData) (TK_OFFSET_RELATIVE)
};

static Tk_ConfigSpec configSpecs[] = {
    {TK_CONFIG_COLOR, "-activefill", NULL, NULL,
	NULL, Tk_Offset(TextItem, activeColor), TK_CONFIG_NULL_OK},
    {TK_CONFIG_BITMAP, "-activestipple", NULL, NULL,
	NULL, Tk_Offset(TextItem, activeStipple), TK_CONFIG_NULL_OK},
    {TK_CONFIG_ANCHOR, "-anchor", NULL, NULL,
	"center", Tk_Offset(TextItem, anchor), TK_CONFIG_DONT_SET_DEFAULT},
    {TK_CONFIG_COLOR, "-disabledfill", NULL, NULL,
	NULL, Tk_Offset(TextItem, disabledColor), TK_CONFIG_NULL_OK},
    {TK_CONFIG_BITMAP, "-disabledstipple", NULL, NULL,
	NULL, Tk_Offset(TextItem, disabledStipple), TK_CONFIG_NULL_OK},
    {TK_CONFIG_COLOR, "-fill", NULL, NULL,
	"black", Tk_Offset(TextItem, color), TK_CONFIG_NULL_OK},
    {TK_CONFIG_FONT, "-font", NULL, NULL,
	DEF_CANVTEXT_FONT, Tk_Offset(TextItem, tkfont), 0},
    {TK_CONFIG_FONT, "-altfont", NULL, NULL,
	"Symbol -14", Tk_Offset(TextItem, tkaltfont), 0},
    {TK_CONFIG_FONT, "-scriptfont", NULL, NULL,
	"Helvetica -10", Tk_Offset(TextItem, tkscriptfont), 0},
    {TK_CONFIG_FONT, "-altscriptfont", NULL, NULL,
	"Symbol -10", Tk_Offset(TextItem, tkaltscriptfont), 0},
    {TK_CONFIG_JUSTIFY, "-justify", NULL, NULL,
	"left", Tk_Offset(TextItem, justify), TK_CONFIG_DONT_SET_DEFAULT},
    {TK_CONFIG_INT, "-removelastdescent", NULL, NULL,
	"0", Tk_Offset(TextItem, brmdes), 0},
    {TK_CONFIG_INT, "-rotation", NULL, NULL,
	"0", Tk_Offset(TextItem, rotation), 0},
    {TK_CONFIG_CUSTOM, "-offset", NULL, NULL,
	"0,0", Tk_Offset(TextItem, tsoffset),
	TK_CONFIG_DONT_SET_DEFAULT, &offsetOption},
    {TK_CONFIG_CUSTOM, "-state", NULL, NULL,
	NULL, Tk_Offset(Tk_Item, state), TK_CONFIG_NULL_OK, &stateOption},
    {TK_CONFIG_BITMAP, "-stipple", NULL, NULL,
	NULL, Tk_Offset(TextItem, stipple), TK_CONFIG_NULL_OK},
    {TK_CONFIG_CUSTOM, "-tags", NULL, NULL,
	NULL, 0, TK_CONFIG_NULL_OK, &tagsOption},
    {TK_CONFIG_STRING, "-text", NULL, NULL,
	"", Tk_Offset(TextItem, text), 0},
    {TK_CONFIG_COLOR, "-textbg", NULL, NULL,
	DEF_CANVAS_BG_COLOR, Tk_Offset(TextItem, bgcolor), 0},
    {TK_CONFIG_INT, "-usebg", NULL, NULL,
	"0", Tk_Offset(TextItem, usebg), 0},
    {TK_CONFIG_INT, "-underline", NULL, NULL,
	"-1", Tk_Offset(TextItem, underline), 0},
    {TK_CONFIG_PIXELS, "-width", NULL, NULL,
	"0", Tk_Offset(TextItem, width), TK_CONFIG_DONT_SET_DEFAULT},
    {TK_CONFIG_END, NULL, NULL, NULL, NULL, 0, 0}
};

/*
 * Prototypes for functions defined in this file:
 */

static void		ComputeTextBbox(Tk_Canvas canvas, TextItem *textPtr);
static int		ConfigureText(Tcl_Interp *interp,
			    Tk_Canvas canvas, Tk_Item *itemPtr, int argc,
			    Tcl_Obj *CONST objv[], int flags);
static int		CreateText(Tcl_Interp *interp,
			    Tk_Canvas canvas, struct Tk_Item *itemPtr,
			    int argc, Tcl_Obj *CONST objv[]);
static void		DeleteText(Tk_Canvas canvas,
			    Tk_Item *itemPtr, Display *display);
static void		DisplayCanvText(Tk_Canvas canvas,
			    Tk_Item *itemPtr, Display *display, Drawable dst,
			    int x, int y, int width, int height);
static int		GetSelText(Tk_Canvas canvas,
			    Tk_Item *itemPtr, int offset, char *buffer,
			    int maxBytes);
static int		GetTextIndex(Tcl_Interp *interp,
			    Tk_Canvas canvas, Tk_Item *itemPtr,
			    Tcl_Obj *obj, int *indexPtr);
static void		ScaleText(Tk_Canvas canvas,
			    Tk_Item *itemPtr, double originX, double originY,
			    double scaleX, double scaleY);
static void		SetTextCursor(Tk_Canvas canvas,
			    Tk_Item *itemPtr, int index);
static int		TextCoords(Tcl_Interp *interp,
			    Tk_Canvas canvas, Tk_Item *itemPtr,
			    int argc, Tcl_Obj *CONST objv[]);
static void		TextDeleteChars(Tk_Canvas canvas,
			    Tk_Item *itemPtr, int first, int last);
static void		TextInsert(Tk_Canvas canvas,
			    Tk_Item *itemPtr, int beforeThis, char *string);
static int		TextToArea(Tk_Canvas canvas,
			    Tk_Item *itemPtr, double *rectPtr);
static double		TextToPoint(Tk_Canvas canvas,
			    Tk_Item *itemPtr, double *pointPtr);
static int		TextToPostscript(Tcl_Interp *interp,
			    Tk_Canvas canvas, Tk_Item *itemPtr, int prepass);
static void		TranslateText(Tk_Canvas canvas,
			    Tk_Item *itemPtr, double deltaX, double deltaY);

int bTkCanvTextDebug=0;

/*
 * The structures below defines the rectangle and oval item types by means of
 * functions that can be invoked by generic item code.
 */

Tk_ItemType tkTextType = {
    "text",			/* name */
    sizeof(TextItem),		/* itemSize */
    CreateText,			/* createProc */
    configSpecs,		/* configSpecs */
    ConfigureText,		/* configureProc */
    TextCoords,			/* coordProc */
    DeleteText,			/* deleteProc */
    DisplayCanvText,		/* displayProc */
    TK_CONFIG_OBJS,		/* flags */
    TextToPoint,		/* pointProc */
    TextToArea,			/* areaProc */
    TextToPostscript,		/* postscriptProc */
    ScaleText,			/* scaleProc */
    TranslateText,		/* translateProc */
    (Tk_ItemIndexProc *) GetTextIndex,/* indexProc */
    SetTextCursor,		/* icursorProc */
    GetSelText,			/* selectionProc */
    TextInsert,			/* insertProc */
    TextDeleteChars,		/* dTextProc */
    NULL,			/* nextPtr */
};

/*
 *--------------------------------------------------------------
 *
 * CreateText --
 *
 *	This function is invoked to create a new text item in a canvas.
 *
 * Results:
 *	A standard Tcl return value. If an error occurred in creating the item
 *	then an error message is left in the interp's result; in this case
 *	itemPtr is left uninitialized so it can be safely freed by the caller.
 *
 * Side effects:
 *	A new text item is created.
 *
 *--------------------------------------------------------------
 */
static int		KBCreateText _ANSI_ARGS_((Tcl_Interp *interp,
			    Tk_Canvas canvas, struct Tk_Item *itemPtr,
			    int objc, Tcl_Obj *CONST objv[]));

static int
CreateText(
    Tcl_Interp *interp,		/* Interpreter for error reporting. */
    Tk_Canvas canvas,		/* Canvas to hold new item. */
    Tk_Item *itemPtr,		/* Record to hold new item; header has been
				 * initialized by caller. */
    int objc,			/* Number of arguments in objv. */
    Tcl_Obj *CONST objv[])	/* Arguments describing rectangle. */
{
    TextItem *textPtr = (TextItem *) itemPtr;
    int i;

    if (objc == 0) {
	Tcl_Panic("canvas did not pass any coords\n");
    }

    /*
     * Carry out initialization that is needed in order to clean up after
     * errors during the the remainder of this function.
     */

    textPtr->textInfoPtr = Tk_CanvasGetTextInfo(canvas);

    textPtr->insertPos	= 0;

    textPtr->anchor	= TK_ANCHOR_CENTER;
    textPtr->tsoffset.flags = 0;
    textPtr->tsoffset.xoffset = 0;
    textPtr->tsoffset.yoffset = 0;
    textPtr->color	= NULL;
    textPtr->activeColor = NULL;
    textPtr->disabledColor = NULL;
    textPtr->tkfont	= NULL;
    textPtr->tkaltfont	= NULL;
    textPtr->tkscriptfont	= NULL;
    textPtr->tkaltscriptfont	= NULL;
    textPtr->justify	= TK_JUSTIFY_LEFT;
    textPtr->stipple	= None;
    textPtr->activeStipple = None;
    textPtr->disabledStipple = None;
    textPtr->text	= NULL;
    textPtr->width	= 0;
    textPtr->underline	= -1;

    textPtr->numChars	= 0;
    textPtr->numBytes	= 0;
    textPtr->textLayout = NULL;
    textPtr->leftEdge	= 0;
    textPtr->rightEdge	= 0;
    textPtr->gc		= None;
    textPtr->altGc		= None;
    textPtr->scriptGc		= None;
    textPtr->altScriptGc		= None;
    textPtr->selTextGC	= None;
    textPtr->cursorOffGC = None;
    textPtr->copyGC = None;
    textPtr->rotation = 0;
    textPtr->usebg = 0;
    textPtr->bgcolor	= NULL;
    textPtr->bgGC = 0;
    textPtr->brmdes = 0;

    if (objc<0)
	return KBCreateText(interp, canvas, itemPtr, objc, objv);

    /*
     * Process the arguments to fill in the item record. Only 1 (list) or 2 (x
     * y) coords are allowed.
     */

    if (objc == 1) {
	i = 1;
    } else {
	char *arg = Tcl_GetString(objv[1]);

	i = 2;
	if ((arg[0] == '-') && (arg[1] >= 'a') && (arg[1] <= 'z')) {
	    i = 1;
	}
    }
    if ((TextCoords(interp, canvas, itemPtr, i, objv) != TCL_OK)) {
	goto error;
    }
    if (ConfigureText(interp, canvas, itemPtr, objc-i, objv+i, 0) == TCL_OK) {
	return TCL_OK;
    }

  error:
    DeleteText(canvas, itemPtr, Tk_Display(Tk_CanvasTkwin(canvas)));
    return TCL_ERROR;
}

static int
KBCreateText(interp, canvas, itemPtr, objc, objv)
    Tcl_Interp *interp;		/* Interpreter for error reporting. */
    Tk_Canvas canvas;		/* Canvas to hold new item. */
    Tk_Item *itemPtr;		/* Record to hold new item; header has been
				 * initialized by caller. */
    int objc;			/* Number of arguments in objv. */
    Tcl_Obj *CONST objv[];	/* Arguments describing rectangle. */
{
    TextItem *textPtr = (TextItem *) itemPtr;
    int nn;
    double *bb;

    objc = -objc;
    nn = *((int *)objv[objc-2]);
    bb = (double *)objv[objc-1];
    objc -= 2;

    /*
     * Process the arguments to fill in the item record.
     */

    if (nn != 2) {
	return TCL_ERROR;
    }
    textPtr->x = bb[0];
    textPtr->y = bb[1];
    ComputeTextBbox(canvas, textPtr);
    if (ConfigureText(interp, canvas, itemPtr, objc, objv, 0) == TCL_OK) {
	return TCL_OK;
    }

    /*error:*/
    DeleteText(canvas, itemPtr, Tk_Display(Tk_CanvasTkwin(canvas)));
    return TCL_ERROR;
}

/*
 *--------------------------------------------------------------
 *
 * TextCoords --
 *
 *	This function is invoked to process the "coords" widget command on
 *	text items. See the user documentation for details on what it does.
 *
 * Results:
 *	Returns TCL_OK or TCL_ERROR, and sets the interp's result.
 *
 * Side effects:
 *	The coordinates for the given item may be changed.
 *
 *--------------------------------------------------------------
 */

static int
TextCoords(
    Tcl_Interp *interp,		/* Used for error reporting. */
    Tk_Canvas canvas,		/* Canvas containing item. */
    Tk_Item *itemPtr,		/* Item whose coordinates are to be read or
				 * modified. */
    int objc,			/* Number of coordinates supplied in objv. */
    Tcl_Obj *CONST objv[])	/* Array of coordinates: x1, y1, x2, y2, ... */
{
    TextItem *textPtr = (TextItem *) itemPtr;

    if (objc == 0) {
	Tcl_Obj *obj = Tcl_NewObj();

	Tcl_Obj *subobj = Tcl_NewDoubleObj(textPtr->x);
	Tcl_ListObjAppendElement(interp, obj, subobj);
	subobj = Tcl_NewDoubleObj(textPtr->y);
	Tcl_ListObjAppendElement(interp, obj, subobj);
	Tcl_SetObjResult(interp, obj);
    } else if (objc < 3) {
	if (objc==1) {
	    if (Tcl_ListObjGetElements(interp, objv[0], &objc,
		    (Tcl_Obj ***) &objv) != TCL_OK) {
		return TCL_ERROR;
	    } else if (objc != 2) {
		char buf[64 + TCL_INTEGER_SPACE];

		sprintf(buf, "wrong # coordinates: expected 2, got %d", objc);
		Tcl_SetResult(interp, buf, TCL_VOLATILE);
		return TCL_ERROR;
	    }
	}
	if ((Tk_CanvasGetCoordFromObj(interp, canvas, objv[0],
		    &textPtr->x) != TCL_OK)
		|| (Tk_CanvasGetCoordFromObj(interp, canvas, objv[1],
		    &textPtr->y) != TCL_OK)) {
	    return TCL_ERROR;
	}
	ComputeTextBbox(canvas, textPtr);
    } else {
	char buf[64 + TCL_INTEGER_SPACE];

	sprintf(buf, "wrong # coordinates: expected 0 or 2, got %d", objc);
	Tcl_SetResult(interp, buf, TCL_VOLATILE);
	return TCL_ERROR;
    }
    return TCL_OK;
}

/*
 *--------------------------------------------------------------
 *
 * ConfigureText --
 *
 *	This function is invoked to configure various aspects of a text item,
 *	such as its border and background colors.
 *
 * Results:
 *	A standard Tcl result code. If an error occurs, then an error message
 *	is left in the interp's result.
 *
 * Side effects:
 *	Configuration information, such as colors and stipple patterns, may be
 *	set for itemPtr.
 *
 *--------------------------------------------------------------
 */

static int
ConfigureText(
    Tcl_Interp *interp,		/* Interpreter for error reporting. */
    Tk_Canvas canvas,		/* Canvas containing itemPtr. */
    Tk_Item *itemPtr,		/* Rectangle item to reconfigure. */
    int objc,			/* Number of elements in objv. */
    Tcl_Obj *CONST objv[],	/* Arguments describing things to configure. */
    int flags)			/* Flags to pass to Tk_ConfigureWidget. */
{
    TextItem *textPtr = (TextItem *) itemPtr;
    XGCValues gcValues;
    GC newGC, newSelGC;
    unsigned long mask = GCFont;
    Tk_Window tkwin;
    Tk_CanvasTextInfo *textInfoPtr = textPtr->textInfoPtr;
    XColor *selBgColorPtr;
    XColor *color;
    Pixmap stipple;
    Tk_State state;

    if (bTkCanvTextDebug)
      printf("TkCanvTextDebug: ConfigureText\n");

    tkwin = Tk_CanvasTkwin(canvas);
    if (TCL_OK != Tk_ConfigureWidget(interp, tkwin, configSpecs, objc,
	    (CONST char **) objv, (char *) textPtr, flags|TK_CONFIG_OBJS)) {
	return TCL_ERROR;
    }

    /*
     * A few of the options require additional processing, such as graphics
     * contexts.
     */

    state = itemPtr->state;

    if (textPtr->activeColor != NULL || textPtr->activeStipple != None) {
	itemPtr->redraw_flags |= TK_ITEM_STATE_DEPENDANT;
    } else {
	itemPtr->redraw_flags &= ~TK_ITEM_STATE_DEPENDANT;
    }

    if(state == TK_STATE_NULL) {
	state = ((TkCanvas *)canvas)->canvas_state;
    }

    color = textPtr->color;
    stipple = textPtr->stipple;
    if (((TkCanvas *)canvas)->currentItemPtr == itemPtr) {
	if (textPtr->activeColor!=NULL) {
	    color = textPtr->activeColor;
	}
	if (textPtr->activeStipple!=None) {
	    stipple = textPtr->activeStipple;
	}
    } else if (state==TK_STATE_DISABLED) {
	if (textPtr->disabledColor!=NULL) {
	    color = textPtr->disabledColor;
	}
	if (textPtr->disabledStipple!=None) {
	    stipple = textPtr->disabledStipple;
	}
    }

    newGC = newSelGC = None;
    if (textPtr->tkfont != NULL) {
	gcValues.font = Tk_FontId(textPtr->tkfont);
	mask = GCFont;
	if (color != NULL) {
	    gcValues.foreground = color->pixel;
	    mask |= GCForeground;
	    if (stipple != None) {
		gcValues.stipple = stipple;
		gcValues.fill_style = FillStippled;
		mask |= GCStipple|GCFillStyle;
	    }
	    newGC = Tk_GetGC(tkwin, mask, &gcValues);
	}
	mask &= ~(GCTile|GCFillStyle|GCStipple);
	if (stipple != None) {
	    gcValues.stipple = stipple;
	    gcValues.fill_style = FillStippled;
	    mask |= GCStipple|GCFillStyle;
	}
	if (textInfoPtr->selFgColorPtr != NULL) {
	    gcValues.foreground = textInfoPtr->selFgColorPtr->pixel;
	}
	newSelGC = Tk_GetGC(tkwin, mask|GCForeground, &gcValues);
    }
    if (textPtr->gc != None) {
	Tk_FreeGC(Tk_Display(tkwin), textPtr->gc);
    }
    textPtr->gc = newGC;
    if (textPtr->selTextGC != None) {
	Tk_FreeGC(Tk_Display(tkwin), textPtr->selTextGC);
    }
    textPtr->selTextGC = newSelGC;

    if (textPtr->copyGC == None) {
	textPtr->copyGC = Tk_GetGC(tkwin, 0, &gcValues);
    }


    gcValues.foreground = textPtr->color->pixel;
    gcValues.font = Tk_FontId(textPtr->tkaltfont);
    newGC = Tk_GetGC(tkwin, mask, &gcValues);
    if (textPtr->altGc != None) {
	Tk_FreeGC(Tk_Display(tkwin), textPtr->altGc);
    }
    textPtr->altGc = newGC;
    gcValues.font = Tk_FontId(textPtr->tkscriptfont);
    newGC = Tk_GetGC(tkwin, mask, &gcValues);
    if (textPtr->scriptGc != None) {
	Tk_FreeGC(Tk_Display(tkwin), textPtr->scriptGc);
    }
    textPtr->scriptGc = newGC;
    gcValues.font = Tk_FontId(textPtr->tkaltscriptfont);
    newGC = Tk_GetGC(tkwin, mask, &gcValues);
    if (textPtr->altScriptGc != None) {
	Tk_FreeGC(Tk_Display(tkwin), textPtr->altScriptGc);
    }
    textPtr->altScriptGc = newGC;
	
    gcValues.foreground = textPtr->bgcolor->pixel;
/*    printf("bgcolor = %X\n",gcValues.foreground);*/
    mask = GCForeground;
    newGC = Tk_GetGC(tkwin, mask, &gcValues);
    if (textPtr->bgGC != None) {
	Tk_FreeGC(Tk_Display(tkwin), textPtr->bgGC);
    }
    textPtr->bgGC = newGC;

    selBgColorPtr = Tk_3DBorderColor(textInfoPtr->selBorder);
    if (Tk_3DBorderColor(textInfoPtr->insertBorder)->pixel
	    == selBgColorPtr->pixel) {
	if (selBgColorPtr->pixel == BlackPixelOfScreen(Tk_Screen(tkwin))) {
	    gcValues.foreground = WhitePixelOfScreen(Tk_Screen(tkwin));
	} else {
	    gcValues.foreground = BlackPixelOfScreen(Tk_Screen(tkwin));
	}
	newGC = Tk_GetGC(tkwin, GCForeground, &gcValues);
    } else {
	newGC = None;
    }
    if (textPtr->cursorOffGC != None) {
	Tk_FreeGC(Tk_Display(tkwin), textPtr->cursorOffGC);
    }
    textPtr->cursorOffGC = newGC;


    /*
     * If the text was changed, move the selection and insertion indices to
     * keep them inside the item.
     */

    textPtr->numBytes = strlen(textPtr->text);
    textPtr->numChars = Tcl_NumUtfChars(textPtr->text, textPtr->numBytes);
    if (textInfoPtr->selItemPtr == itemPtr) {

	if (textInfoPtr->selectFirst >= textPtr->numChars) {
	    textInfoPtr->selItemPtr = NULL;
	} else {
	    if (textInfoPtr->selectLast >= textPtr->numChars) {
		textInfoPtr->selectLast = textPtr->numChars - 1;
	    }
	    if ((textInfoPtr->anchorItemPtr == itemPtr)
		    && (textInfoPtr->selectAnchor >= textPtr->numChars)) {
		textInfoPtr->selectAnchor = textPtr->numChars - 1;
	    }
	}
    }
    if (textPtr->insertPos >= textPtr->numChars) {
	textPtr->insertPos = textPtr->numChars;
    }

    ComputeTextBbox(canvas, textPtr);
    return TCL_OK;
}

/*
 *--------------------------------------------------------------
 *
 * DeleteText --
 *
 *	This function is called to clean up the data structure associated with
 *	a text item.
 *
 * Results:
 *	None.
 *
 * Side effects:
 *	Resources associated with itemPtr are released.
 *
 *--------------------------------------------------------------
 */

static void
DeleteText(
    Tk_Canvas canvas,		/* Info about overall canvas widget. */
    Tk_Item *itemPtr,		/* Item that is being deleted. */
    Display *display)		/* Display containing window for canvas. */
{
    TextItem *textPtr = (TextItem *) itemPtr;

    if (textPtr->color != NULL) {
	Tk_FreeColor(textPtr->color);
    }
    if (textPtr->activeColor != NULL) {
	Tk_FreeColor(textPtr->activeColor);
    }
    if (textPtr->disabledColor != NULL) {
	Tk_FreeColor(textPtr->disabledColor);
    }
    Tk_FreeFont(textPtr->tkfont);
    Tk_FreeFont(textPtr->tkaltfont);
    Tk_FreeFont(textPtr->tkscriptfont);
    Tk_FreeFont(textPtr->tkaltscriptfont);
    if (textPtr->stipple != None) {
	Tk_FreeBitmap(display, textPtr->stipple);
    }
    if (textPtr->activeStipple != None) {
	Tk_FreeBitmap(display, textPtr->activeStipple);
    }
    if (textPtr->disabledStipple != None) {
	Tk_FreeBitmap(display, textPtr->disabledStipple);
    }
    if (textPtr->text != NULL) {
	ckfree(textPtr->text);
    }

    Tk_FreeTextLayout(textPtr->textLayout);
    if (textPtr->gc != None) {
	Tk_FreeGC(display, textPtr->gc);
    }
    if (textPtr->altGc != None) {
	Tk_FreeGC(display, textPtr->altGc);
    }
    if (textPtr->scriptGc != None) {
	Tk_FreeGC(display, textPtr->scriptGc);
    }
    if (textPtr->altScriptGc != None) {
	Tk_FreeGC(display, textPtr->altScriptGc);
    }
    if (textPtr->selTextGC != None) {
	Tk_FreeGC(display, textPtr->selTextGC);
    }
    if (textPtr->cursorOffGC != None) {
	Tk_FreeGC(display, textPtr->cursorOffGC);
    }
    if (textPtr->copyGC != None) {
	Tk_FreeGC(display, textPtr->copyGC);
    }
    if (textPtr->bgGC != None) {
	Tk_FreeGC(display, textPtr->bgGC);
    }
}

static void rotatepoint(textPtr, x1, y1, dx1, dy1, px2, py2)
    TextItem *textPtr;
    int x1,y1,dx1,dy1,*px2,*py2;
{
    switch(textPtr->rotation) {
    case 90:
	*px2 =  (y1-textPtr->y) + textPtr->y;
	*py2 = -(x1-textPtr->x) + textPtr->x;
	break;

    case 180:
	*px2 = -(x1-textPtr->x) + textPtr->x;
	*py2 = -(y1-textPtr->y) + textPtr->y;
	break;

    case 270:
	*px2 = -(y1-textPtr->y) + textPtr->y;
	*py2 =  (x1-textPtr->x) + textPtr->x;
	break;
    }
    *px2 += dx1;
    *py2 += dy1;
}

/*
 *--------------------------------------------------------------
 *
 * ComputeTextBbox --
 *
 *	This function is invoked to compute the bounding box of all the pixels
 *	that may be drawn as part of a text item. In addition, it recomputes
 *	all of the geometry information used to display a text item or check
 *	for mouse hits.
 *
 * Results:
 *	None.
 *
 * Side effects:
 *	The fields x1, y1, x2, and y2 are updated in the header for itemPtr,
 *	and the linePtr structure is regenerated for itemPtr.
 *
 *--------------------------------------------------------------
 */

static void
ComputeTextBbox(
    Tk_Canvas canvas,		/* Canvas that contains item. */
    TextItem *textPtr)		/* Item whose bbox is to be recomputed. */
{
    Tk_CanvasTextInfo *textInfoPtr;
    int leftX, topY, width, height, fudge;
    Tk_State state = textPtr->header.state;

    if(state == TK_STATE_NULL) {
	state = ((TkCanvas *)canvas)->canvas_state;
    }

    if (bTkCanvTextDebug) {
	printf("TkCanvTextDebug: ComputeTextBbox(%s)\n",textPtr->text);
	printf("x,y: %f,%f\n", textPtr->x, textPtr->y);
    }

    Tk_FreeTextLayout(textPtr->textLayout);
    textPtr->bsmplw = (textPtr->text==NULL)||(strchr(textPtr->text,'`') ? 0 : 1);
    if (textPtr->bsmplw&&textPtr->rotation==0)
	textPtr->textLayout = Tk_ComputeTextLayout(textPtr->tkfont,
	    textPtr->text, textPtr->numChars, textPtr->width,
	    textPtr->justify, 0, &width, &height);
    else
	textPtr->textLayout = Tk_KGComputeTextLayout(textPtr->tkfont,
	    textPtr->text, textPtr->numChars, textPtr->width,
	    textPtr->justify, 0, &width, &height, textPtr->tkaltfont,
	    textPtr->tkscriptfont, textPtr->tkaltscriptfont, textPtr->gc,
	    textPtr->altGc, textPtr->scriptGc, textPtr->altScriptGc,
	    textPtr->brmdes);

    if (state == TK_STATE_HIDDEN || textPtr->color == NULL) {
	width = height = 0;
    }

    /*
     * Use overall geometry information to compute the top-left corner of the
     * bounding box for the text item.
     */

    leftX = (int) floor(textPtr->x + 0.5);
    topY = (int) floor(textPtr->y + 0.5);
    switch (textPtr->anchor) {
    case TK_ANCHOR_NW:
    case TK_ANCHOR_N:
    case TK_ANCHOR_NE:
	break;

    case TK_ANCHOR_W:
    case TK_ANCHOR_CENTER:
    case TK_ANCHOR_E:
	topY -= height / 2;
	break;

    case TK_ANCHOR_SW:
    case TK_ANCHOR_S:
    case TK_ANCHOR_SE:
	topY -= height;
	break;
    }
    switch (textPtr->anchor) {
    case TK_ANCHOR_NW:
    case TK_ANCHOR_W:
    case TK_ANCHOR_SW:
	break;

    case TK_ANCHOR_N:
    case TK_ANCHOR_CENTER:
    case TK_ANCHOR_S:
	leftX -= width / 2;
	break;

    case TK_ANCHOR_NE:
    case TK_ANCHOR_E:
    case TK_ANCHOR_SE:
	leftX -= width;
	break;
    }

    textPtr->leftEdge = leftX;
    textPtr->rightEdge = leftX + width;

    /*
     * Last of all, update the bounding box for the item. The item's bounding
     * box includes the bounding box of all its lines, plus an extra fudge
     * factor for the cursor border (which could potentially be quite large).
     */

    /*
    if (!textPtr->bsmplw &&
	  (textPtr->fontPtr->ascent < textPtr->altFontPtr->max_bounds.ascent)){
	textPtr->header.y1 -=
	    textPtr->altFontPtr->max_bounds.ascent-textPtr->fontPtr->ascent;
    }
    */

    if (textPtr->rotation != 0) {
	int currianch = 0, newanchor = 0, anchx1 = 0, anchy1 = 0, anchx2, anchy2
;
	static int lanchc[]={TK_ANCHOR_N,TK_ANCHOR_E,TK_ANCHOR_S,TK_ANCHOR_W};
	static int lanche[]={TK_ANCHOR_NW,TK_ANCHOR_NE,TK_ANCHOR_SE,TK_ANCHOR_SW};
	int x1,x2,y1,y2,nx1,nx2,ny1,ny2,dx,dy;

	textPtr->header.x1 = leftX;
	textPtr->header.y1 = topY;
	textPtr->header.x2 = leftX + width;
	textPtr->header.y2 = topY + height;

	switch(textPtr->anchor) {
	case TK_ANCHOR_CENTER:
	    break;
	case TK_ANCHOR_N:
	case TK_ANCHOR_NW:
	    currianch = 0;
	    break;
	case TK_ANCHOR_E:
	case TK_ANCHOR_NE:
	    currianch = 1;
	    break;
	case TK_ANCHOR_S:
	case TK_ANCHOR_SE:
	    currianch = 2;
	    break;
	case TK_ANCHOR_W:
	case TK_ANCHOR_SW:
	    currianch = 3;
	    break;
	}	    

	switch(textPtr->anchor) {
	case TK_ANCHOR_N:
	case TK_ANCHOR_E:
	case TK_ANCHOR_S:
	case TK_ANCHOR_W:
	    newanchor = lanchc[(currianch+(int)(textPtr->rotation/90))%4];
	    /*printf("newanchor = %d\n", newanchor);*/
	    break;
	
	case TK_ANCHOR_NW:
	case TK_ANCHOR_NE:
	case TK_ANCHOR_SE:
	case TK_ANCHOR_SW:
	    newanchor = lanche[(currianch+(int)(textPtr->rotation/90))%4];
	    /*printf("newanchor = %d\n", newanchor);*/
	    break;

	case TK_ANCHOR_CENTER:
	    newanchor = TK_ANCHOR_CENTER;
	    break;
	}

	switch (newanchor) {
	case TK_ANCHOR_NW:
	case TK_ANCHOR_N:
	case TK_ANCHOR_NE:
	    anchy1 = textPtr->header.y1;
	    break;

	case TK_ANCHOR_W:
	case TK_ANCHOR_CENTER:
	case TK_ANCHOR_E:
	    anchy1 = (textPtr->header.y1+textPtr->header.y2)/2.;
	    break;

	case TK_ANCHOR_SW:
	case TK_ANCHOR_S:
	case TK_ANCHOR_SE:
	    anchy1 = textPtr->header.y2;
	    break;
	}
	switch (newanchor) {
	case TK_ANCHOR_NW:
	case TK_ANCHOR_W:
	case TK_ANCHOR_SW:
	    anchx1 = textPtr->header.x1;
	    break;

	case TK_ANCHOR_N:
	case TK_ANCHOR_CENTER:
	case TK_ANCHOR_S:
	    anchx1 = (textPtr->header.x1+textPtr->header.x2)/2.;
	    break;

	case TK_ANCHOR_NE:
	case TK_ANCHOR_E:
	case TK_ANCHOR_SE:
	    anchx1 = textPtr->header.x2;
	    break;
	}

	rotatepoint(textPtr, anchx1, anchy1, 0, 0, &anchx2, &anchy2);
	textPtr->dx = textPtr->x - anchx2;
	textPtr->dy = textPtr->y - anchy2;
	/*printf("anch:%d %d %d %d\n",anchx1,anchy1,anchx2,anchy2);*/

	x1 = textPtr->header.x1;
	y1 = textPtr->header.y1;
	x2 = textPtr->header.x2;
	y2 = textPtr->header.y2;
	dx = textPtr->dx;
	dy = textPtr->dy;
	/*printf("x1,y1,x2,y2,dx,dy:%d %d %d %d %d %d\n",x1,y1,x2,y2,dx,dy);*/
	switch(textPtr->rotation) {
	case 90:
	    rotatepoint(textPtr, x2, y1, dx, dy, &nx1, &ny1);
	    rotatepoint(textPtr, x1, y2, dx, dy, &nx2, &ny2);
	    break;
	case 180:
	    rotatepoint(textPtr, x2, y2, dx, dy, &nx1, &ny1);
	    rotatepoint(textPtr, x1, y1, dx, dy, &nx2, &ny2);
	    break;
	case 270:
	    rotatepoint(textPtr, x1, y2, dx, dy, &nx1, &ny1);
	    rotatepoint(textPtr, x2, y1, dx, dy, &nx2, &ny2);
	    break;
	}
	/*printf("nx1,ny1,nx2,ny2:%d %d %d %d\n",nx1,ny1,nx2,ny2);*/
	textPtr->header.x1 = nx1;
	textPtr->header.y1 = ny1;
	textPtr->header.x2 = nx2;
	textPtr->header.y2 = ny2;

	textPtr->leftEdge  = nx1;
	textPtr->rightEdge = nx2;
	textPtr->pixoffx = textPtr->header.x1;
	textPtr->pixoffy = textPtr->header.y1;
    } else {
	textInfoPtr = textPtr->textInfoPtr;
	fudge = (textInfoPtr->insertWidth + 1) / 2;
	if (textInfoPtr->selBorderWidth > fudge) {
	    fudge = textInfoPtr->selBorderWidth;
	}
	textPtr->header.x1 = leftX - fudge;
	textPtr->header.y1 = topY;
	textPtr->header.x2 = leftX + width + fudge;
	textPtr->header.y2 = topY + height;
    }
}

/*
 *--------------------------------------------------------------
 *
 * DisplayCanvText --
 *
 *	This function is invoked to draw a text item in a given drawable.
 *
 * Results:
 *	None.
 *
 * Side effects:
 *	ItemPtr is drawn in drawable using the transformation information in
 *	canvas.
 *
 *--------------------------------------------------------------
 */

static void
DisplayCanvText(
    Tk_Canvas canvas,		/* Canvas that contains item. */
    Tk_Item *itemPtr,		/* Item to be displayed. */
    Display *display,		/* Display on which to draw item. */
    Drawable drawable,		/* Pixmap or window in which to draw item. */
    int x, int y, int width, int height)
				/* Describes region of canvas that must be
				 * redisplayed (not used). */
{
    TextItem *textPtr;
    Tk_CanvasTextInfo *textInfoPtr;
    int selFirstChar, selLastChar;
    short drawableX, drawableY;
    Pixmap stipple;
    Tk_State state = itemPtr->state;

    int xoff, yoff, xos, yos, xpos, ypos, w = 0, h = 0, wt = 0, ht = 0;
    int xsi = 0, ysi = 0, wsi, hsi;
    Pixmap pixmap;
    Tk_Window tkwin = Tk_CanvasTkwin(canvas);

    if (bTkCanvTextDebug)
      printf("TkCanvTextDebug: DisplayCanvText\n");

    textPtr = (TextItem *) itemPtr;
    textInfoPtr = textPtr->textInfoPtr;

    if(state == TK_STATE_NULL) {
	state = ((TkCanvas *)canvas)->canvas_state;
    }
    stipple = textPtr->stipple;
    if (((TkCanvas *)canvas)->currentItemPtr == itemPtr) {
	if (textPtr->activeStipple!=None) {
	    stipple = textPtr->activeStipple;
	}
    } else if (state==TK_STATE_DISABLED) {
	if (textPtr->disabledStipple!=None) {
	    stipple = textPtr->disabledStipple;
	}
    }

    if (textPtr->gc == None) {
	return;
    }

    if (textPtr->rotation != 0) {
	int i,j,dxp,dyp,dwp,dhp;
	unsigned int wd,hd,bd,dp;
	short tmpx, tmpy;
	XImage *ims, *imd;
	char *imddat;
	Window root;
	
	w = textPtr->header.x2-textPtr->header.x1;
	h = textPtr->header.y2-textPtr->header.y1;
	if ((w==0)||(h==0))
	    return;
	/*printf("w,h:%d,%d\n",w,h);*/
	xoff = textPtr->header.x1;
	yoff = textPtr->header.y1;
/*	printf("x,yoff:%d %d\n",xoff,yoff);*/
	Tk_CanvasDrawableCoords(canvas, (double) xoff,
		(double) yoff, &tmpx, &tmpy);
	xos = tmpx; yos = tmpy;
	xoff = textPtr->pixoffx;
	yoff = textPtr->pixoffy;
	Tk_CanvasDrawableCoords(canvas, (double) xoff,
		(double) yoff, &tmpx, &tmpy);
	xpos = tmpx; ypos = tmpy;
	/*printf("x,ypos: %d,%d\n", xpos, ypos);*/
/*	printf("x,yoff:%d %d x,yos: %d %d w,h:%d %d\n",xoff,yoff,xos,yos,w,h);*/
/*	XCopyArea(display, drawable, pixmap, textPtr->copyGC, xos, yos, w, h, 0, 0);*/
	XGetGeometry(display, drawable, &root, &i, &j, &wd, &hd, &bd, &dp);
/*	printf("i,j,wd,hd,bd,dp:%d %d %d %d %d %d\n",i,j,wd,hd,bd,dp);*/
	xsi = xos; ysi = yos; wsi = w; hsi = h;
/*	printf("xsi etc.:%d %d %d %d\n",xsi,ysi,wsi,hsi);*/
	dxp = 0; dyp = 0; dwp = 0; dhp = 0;
	if (xsi<0) {
	    wsi += xsi;
	    dxp = -xsi;
	    xsi = 0;
	}
	if (ysi<0) {
	    hsi += ysi;
	    dyp = -ysi;
	    ysi = 0;
	}
	if (xsi+wsi>wd) {
	    dwp = xsi+wsi - wd;
	    wsi -= dwp;
	}
	if (ysi+hsi>hd) {
	    dhp = ysi+hsi - hd;
	    hsi -= dhp;
	}
/*	printf("xsi etc.:%d %d %d %d %d %d\n",xsi,ysi,wsi,hsi,dxp,dyp);*/
	w = wsi; h = hsi;
	switch(textPtr->rotation) {
	case 90:
	    xpos += dhp;
	    ypos += dxp;
	    break;
	case 180:
	    xpos += dwp;
	    ypos += dhp;
	    break;
	case 270:
	    xpos += dyp;
	    ypos += dwp;
	    break;
	}
	if (textPtr->rotation == 180) {
	    wt = w;
	    ht = h;
	} else {
	    wt = h;
	    ht = w;
	}
	pixmap = Tk_GetPixmap(display, Tk_WindowId(tkwin),
	    wt, ht, Tk_Depth(tkwin));
	/*printf("pixmap w,h:%d,%d\n",wt,ht);*/
	if (textPtr->usebg) {
	    XFillRectangle(Tk_Display(tkwin), pixmap, textPtr->bgGC,
		0, 0, (unsigned int) wt, (unsigned int) ht);
	    /*Tk_Fill3DRectangle(tkwin, pixmap, textPtr->bgBorder, 0, 0, wt, ht,
	      0, TK_RELIEF_FLAT);*/
	} else {
	    ims = XGetImage(display, drawable, xsi, ysi, wsi, hsi,
		      AllPlanes, ZPixmap);
	    imddat = malloc(w*h*sizeof(int));
	    imd = XCreateImage(display, Tk_Visual(tkwin), Tk_Depth(tkwin),
		      ZPixmap, 0, imddat, wt, ht, ims->bitmap_pad, 0);
	    if (textPtr->rotation==270)
		for (i=0; i<wt; i++)
		    for (j=0; j<ht; j++)
			XPutPixel(imd, i, ht-j-1, XGetPixel(ims, j, i));
	    else if (textPtr->rotation==180)
		for (i=0; i<h; i++)
		    for (j=0; j<w; j++)
			XPutPixel(imd, w-j-1, h-i-1, XGetPixel(ims, j, i));
	    else if (textPtr->rotation==90)
		for (i=0; i<wt; i++)
		    for (j=0; j<ht; j++)
			XPutPixel(imd, wt-i-1, j, XGetPixel(ims, j, i));
	    TkPutImage(NULL, 0, display, pixmap, textPtr->copyGC, imd, 0, 0, 0, 0, 
		wt, ht);

	    XDestroyImage(ims);
	    XDestroyImage(imd);
	}
    } else {
	pixmap = drawable;
	xpos = 0;
	ypos = 0;
    }
    /*printf("x,ypos: %d,%d\n", xpos, ypos);*/

    /*
     * If we're stippling, then modify the stipple offset in the GC. Be sure
     * to reset the offset when done, since the GC is supposed to be
     * read-only.
     */

    if (stipple != None) {
	Tk_CanvasSetOffset(canvas, textPtr->gc, &textPtr->tsoffset);
    }

    selFirstChar = -1;
    selLastChar = 0;		/* lint. */

    if (textInfoPtr->selItemPtr == itemPtr) {
	selFirstChar = textInfoPtr->selectFirst;
	selLastChar = textInfoPtr->selectLast;
	if (selLastChar > textPtr->numChars) {
	    selLastChar = textPtr->numChars - 1;
	}
	if ((selFirstChar >= 0) && (selFirstChar <= selLastChar)) {
	    int xFirst, yFirst, hFirst;
	    int xLast, yLast, wLast;

	    /*
	     * Draw a special background under the selection.
	     */

	    Tk_CharBbox(textPtr->textLayout, selFirstChar, &xFirst, &yFirst,
		    NULL, &hFirst);
	    Tk_CharBbox(textPtr->textLayout, selLastChar, &xLast, &yLast,
		    &wLast, NULL);

	    /*
	     * If the selection spans the end of this line, then display
	     * selection background all the way to the end of the line.
	     * However, for the last line we only want to display up to the
	     * last character, not the end of the line.
	     */

	    x = xFirst;
	    height = hFirst;
	    for (y = yFirst ; y <= yLast; y += height) {
		if (y == yLast) {
		    width = xLast + wLast - x;
		} else {
		    width = textPtr->rightEdge - textPtr->leftEdge - x;
		}
		Tk_CanvasDrawableCoords(canvas,
			(double) (textPtr->leftEdge + x
				- textInfoPtr->selBorderWidth),
			(double) (textPtr->header.y1 + y),
			&drawableX, &drawableY);
		Tk_Fill3DRectangle(Tk_CanvasTkwin(canvas), pixmap,
			textInfoPtr->selBorder, drawableX, drawableY,
			width + 2 * textInfoPtr->selBorderWidth,
			height, textInfoPtr->selBorderWidth, TK_RELIEF_RAISED);
		x = 0;
	    }
	}
    }

    /*
     * If the insertion point should be displayed, then draw a special
     * background for the cursor before drawing the text. Note: if we're the
     * cursor item but the cursor is turned off, then redraw background over
     * the area of the cursor. This guarantees that the selection won't make
     * the cursor invisible on mono displays, where both are drawn in the same
     * color.
     */

    if ((textInfoPtr->focusItemPtr == itemPtr) && (textInfoPtr->gotFocus)) {
	if (Tk_CharBbox(textPtr->textLayout, textPtr->insertPos,
		&x, &y, NULL, &height)) {
	    Tk_CanvasDrawableCoords(canvas,
		    (double) (textPtr->leftEdge + x
			    - (textInfoPtr->insertWidth / 2)),
		    (double) (textPtr->header.y1 + y),
		    &drawableX, &drawableY);
	    Tk_SetCaretPos(Tk_CanvasTkwin(canvas), drawableX, drawableY,
		    height);
	    if (textInfoPtr->cursorOn) {
		Tk_Fill3DRectangle(Tk_CanvasTkwin(canvas), pixmap,
			textInfoPtr->insertBorder,
			drawableX, drawableY,
			textInfoPtr->insertWidth, height,
			textInfoPtr->insertBorderWidth, TK_RELIEF_RAISED);
	    } else if (textPtr->cursorOffGC != None) {
		/*
		 * Redraw the background over the area of the cursor, even
		 * though the cursor is turned off. This guarantees that the
		 * selection won't make the cursor invisible on mono displays,
		 * where both may be drawn in the same color.
		 */

		XFillRectangle(display, pixmap, textPtr->cursorOffGC,
			drawableX, drawableY,
			(unsigned) textInfoPtr->insertWidth,
			(unsigned) height);
	    }
	}
    }

    /*
     * If there is no selected text or the selected text foreground is the
     * same as the regular text foreground, then draw one text string. If
     * there is selected text and the foregrounds differ, draw the regular
     * text up to the selection, draw the selection, then draw the rest of the
     * regular text. Drawing the regular text and then the selected text over
     * it would causes problems with anti-aliased text because the two
     * anti-aliasing colors would blend together.
     */

    Tk_CanvasDrawableCoords(canvas, (double) textPtr->leftEdge,
	    (double) textPtr->header.y1, &drawableX, &drawableY);

    /*printf("drawableX,Y: %d,%d\n", drawableX, drawableY);*/
    if (textPtr->bsmplw)
	Tk_DrawTextLayout(display, pixmap, textPtr->gc, textPtr->textLayout,
	    drawableX-xpos, drawableY-ypos, 0, -1);
    else
	Tk_KGDrawTextLayout(display, pixmap, textPtr->gc, textPtr->textLayout,
	    drawableX-xpos, drawableY-ypos, 0, -1);
    if ((selFirstChar >= 0) && (textPtr->selTextGC != textPtr->gc)) {
	Tk_DrawTextLayout(display, pixmap, textPtr->selTextGC,
	    textPtr->textLayout, drawableX, drawableY, selFirstChar,
	    selLastChar + 1);
	Tk_DrawTextLayout(display, pixmap, textPtr->gc, textPtr->textLayout,
	    drawableX, drawableY, selLastChar + 1, -1);
    }
    Tk_UnderlineTextLayout(display, drawable, textPtr->gc, textPtr->textLayout,
	    drawableX, drawableY, textPtr->underline);

    if (stipple != None) {
	XSetTSOrigin(display, textPtr->gc, 0, 0);
    }

    /*
    if (xovl0 != 0) {
	int y = drawableY-ypos - textPtr->altFontPtr->max_bounds.ascent;
	if (xovl0>xovl1)
	    xovl1 = xStart-xpos;
	XDrawLine(display, pixmap, textPtr->gc, xovl0, y, xovl1, y);
    }
    */

    if (textPtr->rotation != 0) {
	int i,j;
	XImage *ims, *imd;
	char *imddat;

	ims = XGetImage(display, pixmap, 0, 0, wt, ht, 
	    AllPlanes, ZPixmap);
	imddat = malloc(w*h*sizeof(int));
	imd = XCreateImage(display, Tk_Visual(tkwin), Tk_Depth(tkwin),
            ZPixmap, 0, imddat, w, h, ims->bitmap_pad, 0);
	if (textPtr->rotation==90)
	    for (i=0; i<w; i++)
		for (j=0; j<h; j++)
		    XPutPixel(imd, i, h-j-1, XGetPixel(ims, j, i));
	else if (textPtr->rotation==180)
	    for (i=0; i<h; i++)
		for (j=0; j<w; j++)
		    XPutPixel(imd, w-j-1, h-i-1, XGetPixel(ims, j, i));
	else if (textPtr->rotation==270)
	    for (i=0; i<w; i++)
		for (j=0; j<h; j++)
		    XPutPixel(imd, w-i-1, j, XGetPixel(ims, j, i));
	/*printf("x,yoff:%d %d w,h:%d %d\n",xoff,yoff,w,h);*/
	TkPutImage(NULL, 0, display, drawable, textPtr->copyGC, imd, 0, 0, xsi, ysi, 
		w, h);

	Tk_FreePixmap(display, pixmap);
	XDestroyImage(ims);
	XDestroyImage(imd);
    }
}

/*
 *--------------------------------------------------------------
 *
 * TextInsert --
 *
 *	Insert characters into a text item at a given position.
 *
 * Results:
 *	None.
 *
 * Side effects:
 *	The text in the given item is modified. The cursor and selection
 *	positions are also modified to reflect the insertion.
 *
 *--------------------------------------------------------------
 */

static void
TextInsert(
    Tk_Canvas canvas,		/* Canvas containing text item. */
    Tk_Item *itemPtr,		/* Text item to be modified. */
    int index,			/* Character index before which string is to
				 * be inserted. */
    char *string)		/* New characters to be inserted. */
{
    TextItem *textPtr = (TextItem *) itemPtr;
    int byteIndex, byteCount, charsAdded;
    char *newStr, *text;
    Tk_CanvasTextInfo *textInfoPtr = textPtr->textInfoPtr;

    string = Tcl_GetStringFromObj((Tcl_Obj *) string, &byteCount);

    text = textPtr->text;

    if (index < 0) {
	index = 0;
    }
    if (index > textPtr->numChars) {
	index = textPtr->numChars;
    }
    byteIndex = Tcl_UtfAtIndex(text, index) - text;
    byteCount = strlen(string);
    if (byteCount == 0) {
	return;
    }

    newStr = (char *) ckalloc((unsigned) textPtr->numBytes + byteCount + 1);
    memcpy(newStr, text, (size_t) byteIndex);
    strcpy(newStr + byteIndex, string);
    strcpy(newStr + byteIndex + byteCount, text + byteIndex);

    ckfree(text);
    textPtr->text = newStr;
    charsAdded = Tcl_NumUtfChars(string, byteCount);
    textPtr->numChars += charsAdded;
    textPtr->numBytes += byteCount;

    /*
     * Inserting characters invalidates indices such as those for the
     * selection and cursor. Update the indices appropriately.
     */

    if (textInfoPtr->selItemPtr == itemPtr) {
	if (textInfoPtr->selectFirst >= index) {
	    textInfoPtr->selectFirst += charsAdded;
	}
	if (textInfoPtr->selectLast >= index) {
	    textInfoPtr->selectLast += charsAdded;
	}
	if ((textInfoPtr->anchorItemPtr == itemPtr)
		&& (textInfoPtr->selectAnchor >= index)) {
	    textInfoPtr->selectAnchor += charsAdded;
	}
    }
    if (textPtr->insertPos >= index) {
	textPtr->insertPos += charsAdded;
    }
    ComputeTextBbox(canvas, textPtr);
}

/*
 *--------------------------------------------------------------
 *
 * TextDeleteChars --
 *
 *	Delete one or more characters from a text item.
 *
 * Results:
 *	None.
 *
 * Side effects:
 *	Characters between "first" and "last", inclusive, get deleted from
 *	itemPtr, and things like the selection position get updated.
 *
 *--------------------------------------------------------------
 */

static void
TextDeleteChars(
    Tk_Canvas canvas,		/* Canvas containing itemPtr. */
    Tk_Item *itemPtr,		/* Item in which to delete characters. */
    int first,			/* Character index of first character to
				 * delete. */
    int last)			/* Character index of last character to delete
				 * (inclusive). */
{
    TextItem *textPtr = (TextItem *) itemPtr;
    int byteIndex, byteCount, charsRemoved;
    char *newStr, *text;
    Tk_CanvasTextInfo *textInfoPtr = textPtr->textInfoPtr;

    text = textPtr->text;
    if (first < 0) {
	first = 0;
    }
    if (last >= textPtr->numChars) {
	last = textPtr->numChars - 1;
    }
    if (first > last) {
	return;
    }
    charsRemoved = last + 1 - first;

    byteIndex = Tcl_UtfAtIndex(text, first) - text;
    byteCount = Tcl_UtfAtIndex(text + byteIndex, charsRemoved)
	- (text + byteIndex);

    newStr = (char *) ckalloc((unsigned) (textPtr->numBytes + 1 - byteCount));
    memcpy(newStr, text, (size_t) byteIndex);
    strcpy(newStr + byteIndex, text + byteIndex + byteCount);

    ckfree(text);
    textPtr->text = newStr;
    textPtr->numChars -= charsRemoved;
    textPtr->numBytes -= byteCount;

    /*
     * Update indexes for the selection and cursor to reflect the renumbering
     * of the remaining characters.
     */

    if (textInfoPtr->selItemPtr == itemPtr) {
	if (textInfoPtr->selectFirst > first) {
	    textInfoPtr->selectFirst -= charsRemoved;
	    if (textInfoPtr->selectFirst < first) {
		textInfoPtr->selectFirst = first;
	    }
	}
	if (textInfoPtr->selectLast >= first) {
	    textInfoPtr->selectLast -= charsRemoved;
	    if (textInfoPtr->selectLast < first - 1) {
		textInfoPtr->selectLast = first - 1;
	    }
	}
	if (textInfoPtr->selectFirst > textInfoPtr->selectLast) {
	    textInfoPtr->selItemPtr = NULL;
	}
	if ((textInfoPtr->anchorItemPtr == itemPtr)
		&& (textInfoPtr->selectAnchor > first)) {
	    textInfoPtr->selectAnchor -= charsRemoved;
	    if (textInfoPtr->selectAnchor < first) {
		textInfoPtr->selectAnchor = first;
	    }
	}
    }
    if (textPtr->insertPos > first) {
	textPtr->insertPos -= charsRemoved;
	if (textPtr->insertPos < first) {
	    textPtr->insertPos = first;
	}
    }
    ComputeTextBbox(canvas, textPtr);
    return;
}

/*
 *--------------------------------------------------------------
 *
 * TextToPoint --
 *
 *	Computes the distance from a given point to a given text item, in
 *	canvas units.
 *
 * Results:
 *	The return value is 0 if the point whose x and y coordinates are
 *	pointPtr[0] and pointPtr[1] is inside the text item. If the point
 *	isn't inside the text item then the return value is the distance from
 *	the point to the text item.
 *
 * Side effects:
 *	None.
 *
 *--------------------------------------------------------------
 */

static double
TextToPoint(
    Tk_Canvas canvas,		/* Canvas containing itemPtr. */
    Tk_Item *itemPtr,		/* Item to check against point. */
    double *pointPtr)		/* Pointer to x and y coordinates. */
{
    TextItem *textPtr;
    Tk_State state = itemPtr->state;
    double value;

    if (state == TK_STATE_NULL) {
	state = ((TkCanvas *)canvas)->canvas_state;
    }
    textPtr = (TextItem *) itemPtr;
    value = (double) Tk_DistanceToTextLayout(textPtr->textLayout,
	    (int) pointPtr[0] - textPtr->leftEdge,
	    (int) pointPtr[1] - textPtr->header.y1);

    if ((state == TK_STATE_HIDDEN) || (textPtr->color == NULL) ||
	    (textPtr->text == NULL) || (*textPtr->text == 0)) {
	value = 1.0e36;
    }
    return value;
}

/*
 *--------------------------------------------------------------
 *
 * TextToArea --
 *
 *	This function is called to determine whether an item lies entirely
 *	inside, entirely outside, or overlapping a given rectangle.
 *
 * Results:
 *	-1 is returned if the item is entirely outside the area given by
 *	rectPtr, 0 if it overlaps, and 1 if it is entirely inside the given
 *	area.
 *
 * Side effects:
 *	None.
 *
 *--------------------------------------------------------------
 */

static int
TextToArea(
    Tk_Canvas canvas,		/* Canvas containing itemPtr. */
    Tk_Item *itemPtr,		/* Item to check against rectangle. */
    double *rectPtr)		/* Pointer to array of four coordinates
				 * (x1,y1,x2,y2) describing rectangular
				 * area. */
{
    TextItem *textPtr;
    Tk_State state = itemPtr->state;

    if (state == TK_STATE_NULL) {
	state = ((TkCanvas *)canvas)->canvas_state;
    }

    textPtr = (TextItem *) itemPtr;
    return Tk_IntersectTextLayout(textPtr->textLayout,
	    (int) (rectPtr[0] + 0.5) - textPtr->leftEdge,
	    (int) (rectPtr[1] + 0.5) - textPtr->header.y1,
	    (int) (rectPtr[2] - rectPtr[0] + 0.5),
	    (int) (rectPtr[3] - rectPtr[1] + 0.5));
}

/*
 *--------------------------------------------------------------
 *
 * ScaleText --
 *
 *	This function is invoked to rescale a text item.
 *
 * Results:
 *	None.
 *
 * Side effects:
 *	Scales the position of the text, but not the size of the font for the
 *	text.
 *
 *--------------------------------------------------------------
 */

	/* ARGSUSED */
static void
ScaleText(
    Tk_Canvas canvas,		/* Canvas containing rectangle. */
    Tk_Item *itemPtr,		/* Rectangle to be scaled. */
    double originX, double originY,
				/* Origin about which to scale rect. */
    double scaleX,		/* Amount to scale in X direction. */
    double scaleY)		/* Amount to scale in Y direction. */
{
    TextItem *textPtr = (TextItem *) itemPtr;

    textPtr->x = originX + scaleX*(textPtr->x - originX);
    textPtr->y = originY + scaleY*(textPtr->y - originY);
    ComputeTextBbox(canvas, textPtr);
    return;
}

/*
 *--------------------------------------------------------------
 *
 * TranslateText --
 *
 *	This function is called to move a text item by a given amount.
 *
 * Results:
 *	None.
 *
 * Side effects:
 *	The position of the text item is offset by (xDelta, yDelta), and the
 *	bounding box is updated in the generic part of the item structure.
 *
 *--------------------------------------------------------------
 */

static void
TranslateText(
    Tk_Canvas canvas,		/* Canvas containing item. */
    Tk_Item *itemPtr,		/* Item that is being moved. */
    double deltaX, double deltaY)
				/* Amount by which item is to be moved. */
{
    TextItem *textPtr = (TextItem *) itemPtr;

    textPtr->x += deltaX;
    textPtr->y += deltaY;
    ComputeTextBbox(canvas, textPtr);
}

/*
 *--------------------------------------------------------------
 *
 * GetTextIndex --
 *
 *	Parse an index into a text item and return either its value or an
 *	error.
 *
 * Results:
 *	A standard Tcl result. If all went well, then *indexPtr is filled in
 *	with the index (into itemPtr) corresponding to string. Otherwise an
 *	error message is left in the interp's result.
 *
 * Side effects:
 *	None.
 *
 *--------------------------------------------------------------
 */

static int
GetTextIndex(
    Tcl_Interp *interp,		/* Used for error reporting. */
    Tk_Canvas canvas,		/* Canvas containing item. */
    Tk_Item *itemPtr,		/* Item for which the index is being
				 * specified. */
    Tcl_Obj *obj,		/* Specification of a particular character in
				 * itemPtr's text. */
    int *indexPtr)		/* Where to store converted character
				 * index. */
{
    TextItem *textPtr = (TextItem *) itemPtr;
    int length;
    int c;
    TkCanvas *canvasPtr = (TkCanvas *) canvas;
    Tk_CanvasTextInfo *textInfoPtr = textPtr->textInfoPtr;
    char *string = Tcl_GetStringFromObj(obj, &length);

    c = string[0];

    if ((c == 'e') && (strncmp(string, "end", (unsigned) length) == 0)) {
	*indexPtr = textPtr->numChars;
    } else if ((c == 'i')
	    && (strncmp(string, "insert", (unsigned) length) == 0)) {
	*indexPtr = textPtr->insertPos;
    } else if ((c == 's') && (length >= 5)
	    && (strncmp(string, "sel.first", (unsigned) length) == 0)) {
	if (textInfoPtr->selItemPtr != itemPtr) {
	    Tcl_SetResult(interp, "selection isn't in item", TCL_STATIC);
	    return TCL_ERROR;
	}
	*indexPtr = textInfoPtr->selectFirst;
    } else if ((c == 's') && (length >= 5)
	    && (strncmp(string, "sel.last", (unsigned) length) == 0)) {
	if (textInfoPtr->selItemPtr != itemPtr) {
	    Tcl_SetResult(interp, "selection isn't in item", TCL_STATIC);
	    return TCL_ERROR;
	}
	*indexPtr = textInfoPtr->selectLast;
    } else if (c == '@') {
	int x, y;
	double tmp;
	char *end, *p;

	p = string+1;
	tmp = strtod(p, &end);
	if ((end == p) || (*end != ',')) {
	    goto badIndex;
	}
	x = (int) ((tmp < 0) ? tmp - 0.5 : tmp + 0.5);
	p = end+1;
	tmp = strtod(p, &end);
	if ((end == p) || (*end != 0)) {
	    goto badIndex;
	}
	y = (int) ((tmp < 0) ? tmp - 0.5 : tmp + 0.5);
	*indexPtr = Tk_PointToChar(textPtr->textLayout,
		x + canvasPtr->scrollX1 - textPtr->leftEdge,
		y + canvasPtr->scrollY1 - textPtr->header.y1);
    } else if (Tcl_GetIntFromObj(NULL, obj, indexPtr) == TCL_OK) {
	if (*indexPtr < 0){
	    *indexPtr = 0;
	} else if (*indexPtr > textPtr->numChars) {
	    *indexPtr = textPtr->numChars;
	}
    } else {
	/*
	 * Some of the paths here leave messages in the interp's result, so we
	 * have to clear it out before storing our own message.
	 */

    badIndex:
	Tcl_SetResult(interp, NULL, TCL_STATIC);
	Tcl_AppendResult(interp, "bad index \"", string, "\"", NULL);
	return TCL_ERROR;
    }
    return TCL_OK;
}

/*
 *--------------------------------------------------------------
 *
 * SetTextCursor --
 *
 *	Set the position of the insertion cursor in this item.
 *
 * Results:
 *	None.
 *
 * Side effects:
 *	The cursor position will change.
 *
 *--------------------------------------------------------------
 */

	/* ARGSUSED */
static void
SetTextCursor(
    Tk_Canvas canvas,		/* Record describing canvas widget. */
    Tk_Item *itemPtr,		/* Text item in which cursor position is to be
				 * set. */
    int index)			/* Character index of character just before
				 * which cursor is to be positioned. */
{
    TextItem *textPtr = (TextItem *) itemPtr;

    if (index < 0) {
	textPtr->insertPos = 0;
    } else if (index > textPtr->numChars) {
	textPtr->insertPos = textPtr->numChars;
    } else {
	textPtr->insertPos = index;
    }
}

/*
 *--------------------------------------------------------------
 *
 * GetSelText --
 *
 *	This function is invoked to return the selected portion of a text
 *	item. It is only called when this item has the selection.
 *
 * Results:
 *	The return value is the number of non-NULL bytes stored at buffer.
 *	Buffer is filled (or partially filled) with a NULL-terminated string
 *	containing part or all of the selection, as given by offset and
 *	maxBytes.
 *
 * Side effects:
 *	None.
 *
 *--------------------------------------------------------------
 */

static int
GetSelText(
    Tk_Canvas canvas,		/* Canvas containing selection. */
    Tk_Item *itemPtr,		/* Text item containing selection. */
    int offset,			/* Byte offset within selection of first
				 * character to be returned. */
    char *buffer,		/* Location in which to place selection. */
    int maxBytes)		/* Maximum number of bytes to place at buffer,
				 * not including terminating NULL
				 * character. */
{
    TextItem *textPtr = (TextItem *) itemPtr;
    int byteCount;
    char *text;
    CONST char *selStart, *selEnd;
    Tk_CanvasTextInfo *textInfoPtr = textPtr->textInfoPtr;

    if ((textInfoPtr->selectFirst < 0) ||
	    (textInfoPtr->selectFirst > textInfoPtr->selectLast)) {
	return 0;
    }
    text = textPtr->text;
    selStart = Tcl_UtfAtIndex(text, textInfoPtr->selectFirst);
    selEnd = Tcl_UtfAtIndex(selStart,
	    textInfoPtr->selectLast + 1 - textInfoPtr->selectFirst);
    byteCount = selEnd - selStart - offset;
    if (byteCount > maxBytes) {
	byteCount = maxBytes;
    }
    if (byteCount <= 0) {
	return 0;
    }
    memcpy(buffer, selStart + offset, (size_t) byteCount);
    buffer[byteCount] = '\0';
    return byteCount;
}

/*
 *--------------------------------------------------------------
 *
 * TextToPostscript --
 *
 *	This function is called to generate Postscript for text items.
 *
 * Results:
 *	The return value is a standard Tcl result. If an error occurs in
 *	generating Postscript then an error message is left in the interp's
 *	result, replacing whatever used to be there. If no error occurs, then
 *	Postscript for the item is appended to the result.
 *
 * Side effects:
 *	None.
 *
 *--------------------------------------------------------------
 */

static int
TextToPostscript(
    Tcl_Interp *interp,		/* Leave Postscript or error message here. */
    Tk_Canvas canvas,		/* Information about overall canvas. */
    Tk_Item *itemPtr,		/* Item for which Postscript is wanted. */
    int prepass)		/* 1 means this is a prepass to collect font
				 * information; 0 means final Postscript is
				 * being created. */
{
    TextItem *textPtr = (TextItem *) itemPtr;
    int x, y;
    Tk_FontMetrics fm;
    char *justify;
    char buffer[500];
    XColor *color;
    Pixmap stipple;
    Tk_State state = itemPtr->state;

    if (state == TK_STATE_NULL) {
	state = ((TkCanvas *)canvas)->canvas_state;
    }
    color = textPtr->color;
    stipple = textPtr->stipple;
    if (state == TK_STATE_HIDDEN || textPtr->color == NULL ||
	    textPtr->text == NULL || *textPtr->text == 0) {
	return TCL_OK;
    } else if (((TkCanvas *)canvas)->currentItemPtr == itemPtr) {
	if (textPtr->activeColor!=NULL) {
	    color = textPtr->activeColor;
	}
	if (textPtr->activeStipple!=None) {
	    stipple = textPtr->activeStipple;
	}
    } else if (state==TK_STATE_DISABLED) {
	if (textPtr->disabledColor!=NULL) {
	    color = textPtr->disabledColor;
	}
	if (textPtr->disabledStipple!=None) {
	    stipple = textPtr->disabledStipple;
	}
    }

    if (Tk_CanvasPsFont(interp, canvas, textPtr->tkfont) != TCL_OK) {
	return TCL_ERROR;
    }
    if (prepass != 0) {
	return TCL_OK;
    }
    if (Tk_CanvasPsColor(interp, canvas, color) != TCL_OK) {
	return TCL_ERROR;
    }
    if (stipple != None) {
	Tcl_AppendResult(interp, "/StippleText {\n    ", NULL);
	Tk_CanvasPsStipple(interp, canvas, stipple);
	Tcl_AppendResult(interp, "} bind def\n", NULL);
    }

    sprintf(buffer, "%.15g %.15g [\n", textPtr->x,
	    Tk_CanvasPsY(canvas, textPtr->y));
    Tcl_AppendResult(interp, buffer, NULL);

    Tk_TextLayoutToPostscript(interp, textPtr->textLayout);

    x = 0;  y = 0;  justify = NULL;	/* lint. */
    switch (textPtr->anchor) {
    case TK_ANCHOR_NW:	   x = 0; y = 0; break;
    case TK_ANCHOR_N:	   x = 1; y = 0; break;
    case TK_ANCHOR_NE:	   x = 2; y = 0; break;
    case TK_ANCHOR_E:	   x = 2; y = 1; break;
    case TK_ANCHOR_SE:	   x = 2; y = 2; break;
    case TK_ANCHOR_S:	   x = 1; y = 2; break;
    case TK_ANCHOR_SW:	   x = 0; y = 2; break;
    case TK_ANCHOR_W:	   x = 0; y = 1; break;
    case TK_ANCHOR_CENTER: x = 1; y = 1; break;
    }
    switch (textPtr->justify) {
    case TK_JUSTIFY_LEFT:   justify = "0";   break;
    case TK_JUSTIFY_CENTER: justify = "0.5"; break;
    case TK_JUSTIFY_RIGHT:  justify = "1";   break;
    }

    Tk_GetFontMetrics(textPtr->tkfont, &fm);
    sprintf(buffer, "] %d ", fm.linespace);
    Tcl_AppendResult(interp, buffer, NULL);
    Tcl_PrintDouble(NULL, x / -2.0, buffer);
    Tcl_AppendResult(interp, buffer, NULL);
    Tcl_PrintDouble(NULL, y / 2.0, buffer);
    Tcl_AppendResult(interp, " ", buffer, NULL);
    sprintf(buffer, " %s %s DrawText\n",
	    justify, ((stipple == None) ? "false" : "true"));
    Tcl_AppendResult(interp, buffer, NULL);

    return TCL_OK;
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 4
 * fill-column: 78
 * End:
 */
