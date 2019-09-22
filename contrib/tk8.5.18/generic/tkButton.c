/*
 * tkButton.c --
 *
 *	This module implements a collection of button-like widgets for the Tk
 *	toolkit. The widgets implemented include buttons, checkbuttons,
 *	radiobuttons, and labels.
 *
 * Copyright (c) 1990-1994 The Regents of the University of California.
 * Copyright (c) 1994-1998 Sun Microsystems, Inc.
 *
 * See the file "license.terms" for information on usage and redistribution of
 * this file, and for a DISCLAIMER OF ALL WARRANTIES.
 */

#include "tkInt.h"
#include "tkButton.h"
#include "default.h"

typedef struct ThreadSpecificData {
    int defaultsInitialized;
} ThreadSpecificData;
static Tcl_ThreadDataKey dataKey;

/*
 * Class names for buttons, indexed by one of the type values defined in
 * tkButton.h.
 */

static const char *const classNames[] = {"Label", "Button", "Checkbutton", "Radiobutton"};

/*
 * The following table defines the legal values for the -default option. It is
 * used together with the "enum defaultValue" declaration in tkButton.h.
 */

static const char *const defaultStrings[] = {
    "active", "disabled", "normal", NULL
};

/*
 * The following table defines the legal values for the -state option.
 * It is used together with the "enum state" declaration in tkButton.h.
 */

static const char *const stateStrings[] = {
    "active", "disabled", "normal", NULL
};

/*
 * The following table defines the legal values for the -compound option.
 * It is used with the "enum compound" declaration in tkButton.h
 */

static const char *const compoundStrings[] = {
    "bottom", "center", "left", "none", "right", "top", NULL
};

char tkDefButtonHighlightWidth[TCL_INTEGER_SPACE] = DEF_BUTTON_HIGHLIGHT_WIDTH;
char tkDefButtonPadx[TCL_INTEGER_SPACE] = DEF_BUTTON_PADX;
char tkDefButtonPady[TCL_INTEGER_SPACE] = DEF_BUTTON_PADY;
char tkDefButtonBorderWidth[TCL_INTEGER_SPACE] = DEF_BUTTON_BORDER_WIDTH;
char tkDefLabelHighlightWidth[TCL_INTEGER_SPACE] = DEF_LABEL_HIGHLIGHT_WIDTH;
char tkDefLabelPadx[TCL_INTEGER_SPACE] = DEF_LABCHKRAD_PADX;
char tkDefLabelPady[TCL_INTEGER_SPACE] = DEF_LABCHKRAD_PADY;

/*
 * Information used for parsing configuration options.  There is a
 * separate table for each of the four widget classes.
 */

static const Tk_OptionSpec labelOptionSpecs[] = {
    {TK_OPTION_BORDER, "-activebackground", "activeBackground", "Foreground",
	DEF_BUTTON_ACTIVE_BG_COLOR, -1, Tk_Offset(TkButton, activeBorder),
	0, (ClientData) DEF_BUTTON_ACTIVE_BG_MONO, 0},
    {TK_OPTION_COLOR, "-activeforeground", "activeForeground", "Background",
	DEF_BUTTON_ACTIVE_FG_COLOR, -1, Tk_Offset(TkButton, activeFg),
	TK_OPTION_NULL_OK, (ClientData) DEF_BUTTON_ACTIVE_FG_MONO, 0},
    {TK_OPTION_ANCHOR, "-anchor", "anchor", "Anchor",
	DEF_BUTTON_ANCHOR, -1, Tk_Offset(TkButton, anchor), 0, 0, 0},
    {TK_OPTION_BORDER, "-background", "background", "Background",
	DEF_BUTTON_BG_COLOR, -1, Tk_Offset(TkButton, normalBorder),
	0, (ClientData) DEF_BUTTON_BG_MONO, 0},
    {TK_OPTION_SYNONYM, "-bd", NULL, NULL,
	NULL, 0, -1, 0, (ClientData) "-borderwidth", 0},
    {TK_OPTION_SYNONYM, "-bg", NULL, NULL,
	NULL, 0, -1, 0, (ClientData) "-background", 0},
    {TK_OPTION_BITMAP, "-bitmap", "bitmap", "Bitmap",
	DEF_BUTTON_BITMAP, -1, Tk_Offset(TkButton, bitmap),
	TK_OPTION_NULL_OK, 0, 0},
    {TK_OPTION_PIXELS, "-borderwidth", "borderWidth", "BorderWidth",
	tkDefButtonBorderWidth, Tk_Offset(TkButton, borderWidthPtr),
	Tk_Offset(TkButton, borderWidth), 0, 0, 0},
    {TK_OPTION_STRING_TABLE, "-compound", "compound", "Compound",
	 DEF_BUTTON_COMPOUND, -1, Tk_Offset(TkButton, compound), 0,
	 (ClientData) compoundStrings, 0},
    {TK_OPTION_CURSOR, "-cursor", "cursor", "Cursor",
	DEF_BUTTON_CURSOR, -1, Tk_Offset(TkButton, cursor),
	TK_OPTION_NULL_OK, 0, 0},
    {TK_OPTION_COLOR, "-disabledforeground", "disabledForeground",
	"DisabledForeground", DEF_BUTTON_DISABLED_FG_COLOR,
	-1, Tk_Offset(TkButton, disabledFg), TK_OPTION_NULL_OK,
	(ClientData) DEF_BUTTON_DISABLED_FG_MONO, 0},
    {TK_OPTION_SYNONYM, "-fg", "foreground", NULL,
	NULL, 0, -1, 0, (ClientData) "-foreground", 0},
    {TK_OPTION_FONT, "-font", "font", "Font",
	DEF_BUTTON_FONT, -1, Tk_Offset(TkButton, tkfont), 0, 0, 0},
    {TK_OPTION_FONT, "-altfont", "altfont", "Font",
	DEF_BUTTON_ALTFONT, -1, Tk_Offset(TkButton, tkaltfont), 0, 0, 0},
    {TK_OPTION_FONT, "-scriptfont", "scriptfont", "Font",
	DEF_BUTTON_SCRFONT, -1, Tk_Offset(TkButton, tkscriptfont), 0, 0, 0},
    {TK_OPTION_FONT, "-altscriptfont", "altscriptfont", "Font",
	DEF_BUTTON_ALTSCRFONT, -1, Tk_Offset(TkButton, tkaltscriptfont), 0, 0, 0},
    {TK_OPTION_COLOR, "-foreground", "foreground", "Foreground",
	DEF_BUTTON_FG, -1, Tk_Offset(TkButton, normalFg), 0, 0, 0},
    {TK_OPTION_STRING, "-height", "height", "Height",
	DEF_BUTTON_HEIGHT, Tk_Offset(TkButton, heightPtr), -1, 0, 0, 0},
    {TK_OPTION_BORDER, "-highlightbackground", "highlightBackground",
	"HighlightBackground", DEF_BUTTON_HIGHLIGHT_BG_COLOR,
	-1, Tk_Offset(TkButton, highlightBorder), 0,
	(ClientData) DEF_BUTTON_HIGHLIGHT_BG_MONO, 0},
    {TK_OPTION_COLOR, "-highlightcolor", "highlightColor", "HighlightColor",
	DEF_BUTTON_HIGHLIGHT, -1, Tk_Offset(TkButton, highlightColorPtr),
	0, 0, 0},
    {TK_OPTION_PIXELS, "-highlightthickness", "highlightThickness",
	"HighlightThickness", tkDefLabelHighlightWidth,
	Tk_Offset(TkButton, highlightWidthPtr),
	Tk_Offset(TkButton, highlightWidth), 0, 0, 0},
    {TK_OPTION_STRING, "-image", "image", "Image",
	DEF_BUTTON_IMAGE, Tk_Offset(TkButton, imagePtr), -1,
	TK_OPTION_NULL_OK, 0, 0},
    {TK_OPTION_JUSTIFY, "-justify", "justify", "Justify",
	DEF_BUTTON_JUSTIFY, -1, Tk_Offset(TkButton, justify), 0, 0, 0},
    {TK_OPTION_PIXELS, "-padx", "padX", "Pad",
	tkDefLabelPadx, Tk_Offset(TkButton, padXPtr),
	Tk_Offset(TkButton, padX), 0, 0, 0},
    {TK_OPTION_PIXELS, "-pady", "padY", "Pad",
	tkDefLabelPady, Tk_Offset(TkButton, padYPtr),
	Tk_Offset(TkButton, padY), 0, 0, 0},
    {TK_OPTION_RELIEF, "-relief", "relief", "Relief",
	DEF_LABCHKRAD_RELIEF, -1, Tk_Offset(TkButton, relief), 0, 0, 0},
    {TK_OPTION_INT, "-rotation", "rotation", "Rotation",
	"0", -1, Tk_Offset(TkButton, rotation), 0, 0, 0},
    {TK_OPTION_STRING_TABLE, "-state", "state", "State",
	DEF_BUTTON_STATE, -1, Tk_Offset(TkButton, state),
	0, (ClientData) stateStrings, 0},
    {TK_OPTION_STRING, "-takefocus", "takeFocus", "TakeFocus",
	DEF_LABEL_TAKE_FOCUS, Tk_Offset(TkButton, takeFocusPtr), -1,
	TK_OPTION_NULL_OK, 0, 0},
    {TK_OPTION_STRING, "-text", "text", "Text",
	DEF_BUTTON_TEXT, Tk_Offset(TkButton, textPtr), -1, 0, 0, 0},
    {TK_OPTION_STRING, "-textvariable", "textVariable", "Variable",
	DEF_BUTTON_TEXT_VARIABLE, Tk_Offset(TkButton, textVarNamePtr), -1,
	TK_OPTION_NULL_OK, 0, 0},
    {TK_OPTION_INT, "-underline", "underline", "Underline",
	DEF_BUTTON_UNDERLINE, -1, Tk_Offset(TkButton, underline), 0, 0, 0},
    {TK_OPTION_STRING, "-width", "width", "Width",
	DEF_BUTTON_WIDTH, Tk_Offset(TkButton, widthPtr), -1, 0, 0, 0},
    {TK_OPTION_PIXELS, "-wraplength", "wrapLength", "WrapLength",
	DEF_BUTTON_WRAP_LENGTH, Tk_Offset(TkButton, wrapLengthPtr),
	Tk_Offset(TkButton, wrapLength), 0, 0, 0},
    {TK_OPTION_END, NULL, NULL, NULL, NULL, 0, 0, 0, 0}
};

static const Tk_OptionSpec buttonOptionSpecs[] = {
    {TK_OPTION_BORDER, "-activebackground", "activeBackground", "Foreground",
	DEF_BUTTON_ACTIVE_BG_COLOR, -1, Tk_Offset(TkButton, activeBorder),
	0, (ClientData) DEF_BUTTON_ACTIVE_BG_MONO, 0},
    {TK_OPTION_COLOR, "-activeforeground", "activeForeground", "Background",
	DEF_BUTTON_ACTIVE_FG_COLOR, -1, Tk_Offset(TkButton, activeFg),
	TK_OPTION_NULL_OK, (ClientData) DEF_BUTTON_ACTIVE_FG_MONO, 0},
    {TK_OPTION_ANCHOR, "-anchor", "anchor", "Anchor",
	DEF_BUTTON_ANCHOR, -1, Tk_Offset(TkButton, anchor), 0, 0, 0},
    {TK_OPTION_BORDER, "-background", "background", "Background",
	DEF_BUTTON_BG_COLOR, -1, Tk_Offset(TkButton, normalBorder),
	0, (ClientData) DEF_BUTTON_BG_MONO, 0},
    {TK_OPTION_SYNONYM, "-bd", NULL, NULL,
	NULL, 0, -1, 0, (ClientData) "-borderwidth", 0},
    {TK_OPTION_SYNONYM, "-bg", NULL, NULL,
	NULL, 0, -1, 0, (ClientData) "-background", 0},
    {TK_OPTION_BITMAP, "-bitmap", "bitmap", "Bitmap",
	DEF_BUTTON_BITMAP, -1, Tk_Offset(TkButton, bitmap),
	TK_OPTION_NULL_OK, 0, 0},
    {TK_OPTION_PIXELS, "-borderwidth", "borderWidth", "BorderWidth",
	tkDefButtonBorderWidth, Tk_Offset(TkButton, borderWidthPtr),
	Tk_Offset(TkButton, borderWidth), 0, 0, 0},
    {TK_OPTION_STRING, "-command", "command", "Command",
	DEF_BUTTON_COMMAND, Tk_Offset(TkButton, commandPtr), -1,
	TK_OPTION_NULL_OK, 0, 0},
    {TK_OPTION_STRING_TABLE, "-compound", "compound", "Compound",
	 DEF_BUTTON_COMPOUND, -1, Tk_Offset(TkButton, compound), 0,
	 (ClientData) compoundStrings, 0},
    {TK_OPTION_CURSOR, "-cursor", "cursor", "Cursor",
	DEF_BUTTON_CURSOR, -1, Tk_Offset(TkButton, cursor),
	TK_OPTION_NULL_OK, 0, 0},
    {TK_OPTION_STRING_TABLE, "-default", "default", "Default",
        DEF_BUTTON_DEFAULT, -1, Tk_Offset(TkButton, defaultState),
	0, (ClientData) defaultStrings, 0},
    {TK_OPTION_COLOR, "-disabledforeground", "disabledForeground",
	"DisabledForeground", DEF_BUTTON_DISABLED_FG_COLOR,
	-1, Tk_Offset(TkButton, disabledFg), TK_OPTION_NULL_OK,
	(ClientData) DEF_BUTTON_DISABLED_FG_MONO, 0},
    {TK_OPTION_SYNONYM, "-fg", "foreground", NULL,
	NULL, 0, -1, 0, (ClientData) "-foreground", 0},
    {TK_OPTION_FONT, "-font", "font", "Font",
	DEF_BUTTON_FONT, -1, Tk_Offset(TkButton, tkfont), 0, 0, 0},
    {TK_OPTION_FONT, "-altfont", "altfont", "Font",
	DEF_BUTTON_ALTFONT, -1, Tk_Offset(TkButton, tkaltfont), 0, 0, 0},
    {TK_OPTION_FONT, "-scriptfont", "scriptfont", "Font",
	DEF_BUTTON_SCRFONT, -1, Tk_Offset(TkButton, tkscriptfont), 0, 0, 0},
    {TK_OPTION_FONT, "-altscriptfont", "altscriptfont", "Font",
	DEF_BUTTON_ALTSCRFONT, -1, Tk_Offset(TkButton, tkaltscriptfont), 0, 0, 0},
    {TK_OPTION_COLOR, "-foreground", "foreground", "Foreground",
	DEF_BUTTON_FG, -1, Tk_Offset(TkButton, normalFg), 0, 0, 0},
    {TK_OPTION_STRING, "-height", "height", "Height",
	DEF_BUTTON_HEIGHT, Tk_Offset(TkButton, heightPtr), -1, 0, 0, 0},
    {TK_OPTION_BORDER, "-highlightbackground", "highlightBackground",
	"HighlightBackground", DEF_BUTTON_HIGHLIGHT_BG_COLOR,
	-1, Tk_Offset(TkButton, highlightBorder), 0,
	(ClientData) DEF_BUTTON_HIGHLIGHT_BG_MONO, 0},
    {TK_OPTION_COLOR, "-highlightcolor", "highlightColor", "HighlightColor",
	DEF_BUTTON_HIGHLIGHT, -1, Tk_Offset(TkButton, highlightColorPtr),
	0, 0, 0},
    {TK_OPTION_PIXELS, "-highlightthickness", "highlightThickness",
	"HighlightThickness", tkDefButtonHighlightWidth,
	Tk_Offset(TkButton, highlightWidthPtr),
	Tk_Offset(TkButton, highlightWidth), 0, 0, 0},
    {TK_OPTION_STRING, "-image", "image", "Image",
	DEF_BUTTON_IMAGE, Tk_Offset(TkButton, imagePtr), -1,
	TK_OPTION_NULL_OK, 0, 0},
    {TK_OPTION_JUSTIFY, "-justify", "justify", "Justify",
	DEF_BUTTON_JUSTIFY, -1, Tk_Offset(TkButton, justify), 0, 0, 0},
    {TK_OPTION_RELIEF, "-overrelief", "overRelief", "OverRelief",
	 DEF_BUTTON_OVER_RELIEF, -1, Tk_Offset(TkButton, overRelief),
	 TK_OPTION_NULL_OK, 0, 0},
    {TK_OPTION_PIXELS, "-padx", "padX", "Pad",
	tkDefButtonPadx, Tk_Offset(TkButton, padXPtr),
	Tk_Offset(TkButton, padX), 0, 0, 0},
    {TK_OPTION_PIXELS, "-pady", "padY", "Pad",
	tkDefButtonPady, Tk_Offset(TkButton, padYPtr),
	Tk_Offset(TkButton, padY), 0, 0, 0},
    {TK_OPTION_RELIEF, "-relief", "relief", "Relief",
	DEF_BUTTON_RELIEF, -1, Tk_Offset(TkButton, relief),
	 0, 0, 0},
    {TK_OPTION_INT, "-repeatdelay", "repeatDelay", "RepeatDelay",
	 DEF_BUTTON_REPEAT_DELAY, -1, Tk_Offset(TkButton, repeatDelay),
	 0, 0, 0},
    {TK_OPTION_INT, "-repeatinterval", "repeatInterval", "RepeatInterval",
	 DEF_BUTTON_REPEAT_INTERVAL, -1, Tk_Offset(TkButton, repeatInterval),
	 0, 0, 0},
    {TK_OPTION_INT, "-rotation", "rotation", "Rotation",
	"0", -1, Tk_Offset(TkButton, rotation), 0, 0, 0},
    {TK_OPTION_STRING_TABLE, "-state", "state", "State",
	DEF_BUTTON_STATE, -1, Tk_Offset(TkButton, state),
	0, (ClientData) stateStrings, 0},
    {TK_OPTION_STRING, "-takefocus", "takeFocus", "TakeFocus",
	DEF_BUTTON_TAKE_FOCUS, Tk_Offset(TkButton, takeFocusPtr), -1,
	TK_OPTION_NULL_OK, 0, 0},
    {TK_OPTION_STRING, "-text", "text", "Text",
	DEF_BUTTON_TEXT, Tk_Offset(TkButton, textPtr), -1, 0, 0, 0},
    {TK_OPTION_STRING, "-textvariable", "textVariable", "Variable",
	DEF_BUTTON_TEXT_VARIABLE, Tk_Offset(TkButton, textVarNamePtr), -1,
	TK_OPTION_NULL_OK, 0, 0},
    {TK_OPTION_INT, "-underline", "underline", "Underline",
	DEF_BUTTON_UNDERLINE, -1, Tk_Offset(TkButton, underline), 0, 0, 0},
    {TK_OPTION_STRING, "-width", "width", "Width",
	DEF_BUTTON_WIDTH, Tk_Offset(TkButton, widthPtr), -1, 0, 0, 0},
    {TK_OPTION_PIXELS, "-wraplength", "wrapLength", "WrapLength",
	DEF_BUTTON_WRAP_LENGTH, Tk_Offset(TkButton, wrapLengthPtr),
	Tk_Offset(TkButton, wrapLength), 0, 0, 0},
    {TK_OPTION_END, NULL, NULL, NULL, NULL, 0, -1, 0, 0, 0}
};

static const Tk_OptionSpec checkbuttonOptionSpecs[] = {
    {TK_OPTION_BORDER, "-activebackground", "activeBackground", "Foreground",
	DEF_BUTTON_ACTIVE_BG_COLOR, -1, Tk_Offset(TkButton, activeBorder),
	0, (ClientData) DEF_BUTTON_ACTIVE_BG_MONO, 0},
    {TK_OPTION_COLOR, "-activeforeground", "activeForeground", "Background",
	DEF_CHKRAD_ACTIVE_FG_COLOR, -1, Tk_Offset(TkButton, activeFg),
	TK_OPTION_NULL_OK, (ClientData) DEF_BUTTON_ACTIVE_FG_MONO, 0},
    {TK_OPTION_ANCHOR, "-anchor", "anchor", "Anchor",
	DEF_BUTTON_ANCHOR, -1, Tk_Offset(TkButton, anchor), 0, 0, 0},
    {TK_OPTION_BORDER, "-background", "background", "Background",
	DEF_BUTTON_BG_COLOR, -1, Tk_Offset(TkButton, normalBorder),
	0, (ClientData) DEF_BUTTON_BG_MONO, 0},
    {TK_OPTION_SYNONYM, "-bd", NULL, NULL,
	NULL, 0, -1, 0, (ClientData) "-borderwidth", 0},
    {TK_OPTION_SYNONYM, "-bg", NULL, NULL,
	NULL, 0, -1, 0, (ClientData) "-background", 0},
    {TK_OPTION_BITMAP, "-bitmap", "bitmap", "Bitmap",
	DEF_BUTTON_BITMAP, -1, Tk_Offset(TkButton, bitmap),
	TK_OPTION_NULL_OK, 0, 0},
    {TK_OPTION_PIXELS, "-borderwidth", "borderWidth", "BorderWidth",
	tkDefButtonBorderWidth, Tk_Offset(TkButton, borderWidthPtr),
	Tk_Offset(TkButton, borderWidth), 0, 0, 0},
    {TK_OPTION_STRING, "-command", "command", "Command",
	DEF_BUTTON_COMMAND, Tk_Offset(TkButton, commandPtr), -1,
	TK_OPTION_NULL_OK, 0, 0},
    {TK_OPTION_STRING_TABLE, "-compound", "compound", "Compound",
	 DEF_BUTTON_COMPOUND, -1, Tk_Offset(TkButton, compound), 0,
	 (ClientData) compoundStrings, 0},
    {TK_OPTION_CURSOR, "-cursor", "cursor", "Cursor",
	DEF_BUTTON_CURSOR, -1, Tk_Offset(TkButton, cursor),
	TK_OPTION_NULL_OK, 0, 0},
    {TK_OPTION_COLOR, "-disabledforeground", "disabledForeground",
	"DisabledForeground", DEF_BUTTON_DISABLED_FG_COLOR,
	-1, Tk_Offset(TkButton, disabledFg), TK_OPTION_NULL_OK,
	(ClientData) DEF_BUTTON_DISABLED_FG_MONO, 0},
    {TK_OPTION_SYNONYM, "-fg", "foreground", NULL,
	NULL, 0, -1, 0, (ClientData) "-foreground", 0},
    {TK_OPTION_FONT, "-font", "font", "Font",
	DEF_BUTTON_FONT, -1, Tk_Offset(TkButton, tkfont), 0, 0, 0},
    {TK_OPTION_FONT, "-altfont", "altfont", "Font",
	DEF_BUTTON_ALTFONT, -1, Tk_Offset(TkButton, tkaltfont), 0, 0, 0},
    {TK_OPTION_FONT, "-scriptfont", "scriptfont", "Font",
	DEF_BUTTON_SCRFONT, -1, Tk_Offset(TkButton, tkscriptfont), 0, 0, 0},
    {TK_OPTION_FONT, "-altscriptfont", "altscriptfont", "Font",
	DEF_BUTTON_ALTSCRFONT, -1, Tk_Offset(TkButton, tkaltscriptfont), 0, 0, 0},
    {TK_OPTION_COLOR, "-foreground", "foreground", "Foreground",
	DEF_CHKRAD_FG, -1, Tk_Offset(TkButton, normalFg), 0, 0, 0},
    {TK_OPTION_STRING, "-height", "height", "Height",
	DEF_BUTTON_HEIGHT, Tk_Offset(TkButton, heightPtr), -1, 0, 0, 0},
    {TK_OPTION_BORDER, "-highlightbackground", "highlightBackground",
	"HighlightBackground", DEF_BUTTON_HIGHLIGHT_BG_COLOR,
	-1, Tk_Offset(TkButton, highlightBorder), 0,
	(ClientData) DEF_BUTTON_HIGHLIGHT_BG_MONO, 0},
    {TK_OPTION_COLOR, "-highlightcolor", "highlightColor", "HighlightColor",
	DEF_BUTTON_HIGHLIGHT, -1, Tk_Offset(TkButton, highlightColorPtr),
	0, 0, 0},
    {TK_OPTION_PIXELS, "-highlightthickness", "highlightThickness",
	"HighlightThickness", tkDefButtonHighlightWidth,
	Tk_Offset(TkButton, highlightWidthPtr),
	Tk_Offset(TkButton, highlightWidth), 0, 0, 0},
    {TK_OPTION_STRING, "-image", "image", "Image",
	DEF_BUTTON_IMAGE, Tk_Offset(TkButton, imagePtr), -1,
	TK_OPTION_NULL_OK, 0, 0},
    {TK_OPTION_BOOLEAN, "-indicatoron", "indicatorOn", "IndicatorOn",
	DEF_BUTTON_INDICATOR, -1, Tk_Offset(TkButton, indicatorOn), 0, 0, 0},
    {TK_OPTION_JUSTIFY, "-justify", "justify", "Justify",
	DEF_BUTTON_JUSTIFY, -1, Tk_Offset(TkButton, justify), 0, 0, 0},
    {TK_OPTION_RELIEF, "-offrelief", "offRelief", "OffRelief",
	 DEF_BUTTON_RELIEF, -1, Tk_Offset(TkButton, offRelief), 0, 0, 0},
    {TK_OPTION_STRING, "-offvalue", "offValue", "Value",
	DEF_BUTTON_OFF_VALUE, Tk_Offset(TkButton, offValuePtr), -1, 0, 0, 0},
    {TK_OPTION_STRING, "-onvalue", "onValue", "Value",
	DEF_BUTTON_ON_VALUE, Tk_Offset(TkButton, onValuePtr), -1, 0, 0, 0},
    {TK_OPTION_RELIEF, "-overrelief", "overRelief", "OverRelief",
	 DEF_BUTTON_OVER_RELIEF, -1, Tk_Offset(TkButton, overRelief),
	 TK_OPTION_NULL_OK, 0, 0},
    {TK_OPTION_PIXELS, "-padx", "padX", "Pad",
	tkDefLabelPadx, Tk_Offset(TkButton, padXPtr),
	Tk_Offset(TkButton, padX), 0, 0, 0},
    {TK_OPTION_PIXELS, "-pady", "padY", "Pad",
	tkDefLabelPady, Tk_Offset(TkButton, padYPtr),
	Tk_Offset(TkButton, padY), 0, 0, 0},
    {TK_OPTION_RELIEF, "-relief", "relief", "Relief",
	DEF_LABCHKRAD_RELIEF, -1, Tk_Offset(TkButton, relief), 0, 0, 0},
    {TK_OPTION_INT, "-rotation", "rotation", "Rotation",
	"0", -1, Tk_Offset(TkButton, rotation), 0, 0, 0},
    {TK_OPTION_BORDER, "-selectcolor", "selectColor", "Background",
	DEF_BUTTON_SELECT_COLOR, -1, Tk_Offset(TkButton, selectBorder),
	TK_OPTION_NULL_OK, (ClientData) DEF_BUTTON_SELECT_MONO, 0},
    {TK_OPTION_STRING, "-selectimage", "selectImage", "SelectImage",
	DEF_BUTTON_SELECT_IMAGE, Tk_Offset(TkButton, selectImagePtr), -1,
	TK_OPTION_NULL_OK, 0, 0},
    {TK_OPTION_STRING_TABLE, "-state", "state", "State",
	DEF_BUTTON_STATE, -1, Tk_Offset(TkButton, state),
	0, (ClientData) stateStrings, 0},
    {TK_OPTION_STRING, "-takefocus", "takeFocus", "TakeFocus",
	DEF_BUTTON_TAKE_FOCUS, Tk_Offset(TkButton, takeFocusPtr), -1,
	TK_OPTION_NULL_OK, 0, 0},
    {TK_OPTION_STRING, "-text", "text", "Text",
	DEF_BUTTON_TEXT, Tk_Offset(TkButton, textPtr), -1, 0, 0, 0},
    {TK_OPTION_STRING, "-textvariable", "textVariable", "Variable",
	DEF_BUTTON_TEXT_VARIABLE, Tk_Offset(TkButton, textVarNamePtr), -1,
	TK_OPTION_NULL_OK, 0, 0},
    {TK_OPTION_STRING, "-tristateimage", "tristateImage", "TristateImage",
	DEF_BUTTON_IMAGE, Tk_Offset(TkButton, tristateImagePtr), -1,
	TK_OPTION_NULL_OK, 0, 0},
    {TK_OPTION_STRING, "-tristatevalue", "tristateValue", "TristateValue",
	DEF_BUTTON_TRISTATE_VALUE, Tk_Offset(TkButton, tristateValuePtr), -1, 0, 0, 0},
    {TK_OPTION_INT, "-underline", "underline", "Underline",
	DEF_BUTTON_UNDERLINE, -1, Tk_Offset(TkButton, underline), 0, 0, 0},
    {TK_OPTION_STRING, "-variable", "variable", "Variable",
	DEF_CHECKBUTTON_VARIABLE, Tk_Offset(TkButton, selVarNamePtr), -1,
	TK_OPTION_NULL_OK, 0, 0},
    {TK_OPTION_STRING, "-width", "width", "Width",
	DEF_BUTTON_WIDTH, Tk_Offset(TkButton, widthPtr), -1, 0, 0, 0},
    {TK_OPTION_PIXELS, "-wraplength", "wrapLength", "WrapLength",
	DEF_BUTTON_WRAP_LENGTH, Tk_Offset(TkButton, wrapLengthPtr),
	Tk_Offset(TkButton, wrapLength), 0, 0, 0},
    {TK_OPTION_END, NULL, NULL, NULL, NULL, 0, -1, 0, 0, 0}
};

static const Tk_OptionSpec radiobuttonOptionSpecs[] = {
    {TK_OPTION_BORDER, "-activebackground", "activeBackground", "Foreground",
	DEF_BUTTON_ACTIVE_BG_COLOR, -1, Tk_Offset(TkButton, activeBorder),
	0, (ClientData) DEF_BUTTON_ACTIVE_BG_MONO, 0},
    {TK_OPTION_COLOR, "-activeforeground", "activeForeground", "Background",
	DEF_CHKRAD_ACTIVE_FG_COLOR, -1, Tk_Offset(TkButton, activeFg),
	TK_OPTION_NULL_OK, (ClientData) DEF_BUTTON_ACTIVE_FG_MONO, 0},
    {TK_OPTION_ANCHOR, "-anchor", "anchor", "Anchor",
	DEF_BUTTON_ANCHOR, -1, Tk_Offset(TkButton, anchor), 0, 0, 0},
    {TK_OPTION_BORDER, "-background", "background", "Background",
	DEF_BUTTON_BG_COLOR, -1, Tk_Offset(TkButton, normalBorder),
	0, (ClientData) DEF_BUTTON_BG_MONO, 0},
    {TK_OPTION_SYNONYM, "-bd", NULL, NULL,
	NULL, 0, -1, 0, (ClientData) "-borderwidth", 0},
    {TK_OPTION_SYNONYM, "-bg", NULL, NULL,
	NULL, 0, -1, 0, (ClientData) "-background", 0},
    {TK_OPTION_BITMAP, "-bitmap", "bitmap", "Bitmap",
	DEF_BUTTON_BITMAP, -1, Tk_Offset(TkButton, bitmap),
	TK_OPTION_NULL_OK, 0, 0},
    {TK_OPTION_PIXELS, "-borderwidth", "borderWidth", "BorderWidth",
	tkDefButtonBorderWidth, Tk_Offset(TkButton, borderWidthPtr),
	Tk_Offset(TkButton, borderWidth), 0, 0, 0},
    {TK_OPTION_STRING, "-command", "command", "Command",
	DEF_BUTTON_COMMAND, Tk_Offset(TkButton, commandPtr), -1,
	TK_OPTION_NULL_OK, 0, 0},
    {TK_OPTION_STRING_TABLE, "-compound", "compound", "Compound",
	 DEF_BUTTON_COMPOUND, -1, Tk_Offset(TkButton, compound), 0,
	 (ClientData) compoundStrings, 0},
    {TK_OPTION_CURSOR, "-cursor", "cursor", "Cursor",
	DEF_BUTTON_CURSOR, -1, Tk_Offset(TkButton, cursor),
	TK_OPTION_NULL_OK, 0, 0},
    {TK_OPTION_COLOR, "-disabledforeground", "disabledForeground",
	"DisabledForeground", DEF_BUTTON_DISABLED_FG_COLOR,
	-1, Tk_Offset(TkButton, disabledFg), TK_OPTION_NULL_OK,
	(ClientData) DEF_BUTTON_DISABLED_FG_MONO, 0},
    {TK_OPTION_SYNONYM, "-fg", "foreground", NULL,
	NULL, 0, -1, 0, (ClientData) "-foreground", 0},
    {TK_OPTION_FONT, "-font", "font", "Font",
	DEF_BUTTON_FONT, -1, Tk_Offset(TkButton, tkfont), 0, 0, 0},
    {TK_OPTION_FONT, "-altfont", "altfont", "Font",
	DEF_BUTTON_ALTFONT, -1, Tk_Offset(TkButton, tkaltfont), 0, 0, 0},
    {TK_OPTION_FONT, "-scriptfont", "scriptfont", "Font",
	DEF_BUTTON_SCRFONT, -1, Tk_Offset(TkButton, tkscriptfont), 0, 0, 0},
    {TK_OPTION_FONT, "-altscriptfont", "altscriptfont", "Font",
	DEF_BUTTON_ALTSCRFONT, -1, Tk_Offset(TkButton, tkaltscriptfont), 0, 0, 0},
    {TK_OPTION_COLOR, "-foreground", "foreground", "Foreground",
	DEF_CHKRAD_FG, -1, Tk_Offset(TkButton, normalFg), 0, 0, 0},
    {TK_OPTION_STRING, "-height", "height", "Height",
	DEF_BUTTON_HEIGHT, Tk_Offset(TkButton, heightPtr), -1, 0, 0, 0},
    {TK_OPTION_BORDER, "-highlightbackground", "highlightBackground",
	"HighlightBackground", DEF_BUTTON_HIGHLIGHT_BG_COLOR,
	-1, Tk_Offset(TkButton, highlightBorder), 0,
	(ClientData) DEF_BUTTON_HIGHLIGHT_BG_MONO, 0},
    {TK_OPTION_COLOR, "-highlightcolor", "highlightColor", "HighlightColor",
	DEF_BUTTON_HIGHLIGHT, -1, Tk_Offset(TkButton, highlightColorPtr),
	0, 0, 0},
    {TK_OPTION_PIXELS, "-highlightthickness", "highlightThickness",
	"HighlightThickness", tkDefButtonHighlightWidth,
	Tk_Offset(TkButton, highlightWidthPtr),
	Tk_Offset(TkButton, highlightWidth), 0, 0, 0},
    {TK_OPTION_STRING, "-image", "image", "Image",
	DEF_BUTTON_IMAGE, Tk_Offset(TkButton, imagePtr), -1,
	TK_OPTION_NULL_OK, 0, 0},
    {TK_OPTION_BOOLEAN, "-indicatoron", "indicatorOn", "IndicatorOn",
	DEF_BUTTON_INDICATOR, -1, Tk_Offset(TkButton, indicatorOn),
	0, 0, 0},
    {TK_OPTION_JUSTIFY, "-justify", "justify", "Justify",
	DEF_BUTTON_JUSTIFY, -1, Tk_Offset(TkButton, justify), 0, 0, 0},
    {TK_OPTION_RELIEF, "-offrelief", "offRelief", "OffRelief",
	 DEF_BUTTON_RELIEF, -1, Tk_Offset(TkButton, offRelief), 0, 0, 0},
    {TK_OPTION_RELIEF, "-overrelief", "overRelief", "OverRelief",
	 DEF_BUTTON_OVER_RELIEF, -1, Tk_Offset(TkButton, overRelief),
	 TK_OPTION_NULL_OK, 0, 0},
    {TK_OPTION_PIXELS, "-padx", "padX", "Pad",
	tkDefLabelPadx, Tk_Offset(TkButton, padXPtr),
	Tk_Offset(TkButton, padX), 0, 0, 0},
    {TK_OPTION_PIXELS, "-pady", "padY", "Pad",
	tkDefLabelPady, Tk_Offset(TkButton, padYPtr),
	Tk_Offset(TkButton, padY), 0, 0, 0},
    {TK_OPTION_RELIEF, "-relief", "relief", "Relief",
	DEF_LABCHKRAD_RELIEF, -1, Tk_Offset(TkButton, relief), 0, 0, 0},
    {TK_OPTION_INT, "-rotation", "rotation", "Rotation",
	"0", -1, Tk_Offset(TkButton, rotation), 0, 0, 0},
    {TK_OPTION_BORDER, "-selectcolor", "selectColor", "Background",
	DEF_BUTTON_SELECT_COLOR, -1, Tk_Offset(TkButton, selectBorder),
	TK_OPTION_NULL_OK, (ClientData) DEF_BUTTON_SELECT_MONO, 0},
    {TK_OPTION_STRING, "-selectimage", "selectImage", "SelectImage",
	DEF_BUTTON_SELECT_IMAGE, Tk_Offset(TkButton, selectImagePtr), -1,
	TK_OPTION_NULL_OK, 0, 0},
    {TK_OPTION_STRING_TABLE, "-state", "state", "State",
	DEF_BUTTON_STATE, -1, Tk_Offset(TkButton, state),
	0, (ClientData) stateStrings, 0},
    {TK_OPTION_STRING, "-takefocus", "takeFocus", "TakeFocus",
	DEF_BUTTON_TAKE_FOCUS, Tk_Offset(TkButton, takeFocusPtr), -1,
	TK_OPTION_NULL_OK, 0, 0},
    {TK_OPTION_STRING, "-text", "text", "Text",
	DEF_BUTTON_TEXT, Tk_Offset(TkButton, textPtr), -1, 0, 0, 0},
    {TK_OPTION_STRING, "-textvariable", "textVariable", "Variable",
	DEF_BUTTON_TEXT_VARIABLE, Tk_Offset(TkButton, textVarNamePtr), -1,
	TK_OPTION_NULL_OK, 0, 0},
    {TK_OPTION_STRING, "-tristateimage", "tristateImage", "TristateImage",
	DEF_BUTTON_IMAGE, Tk_Offset(TkButton, tristateImagePtr), -1,
	TK_OPTION_NULL_OK, 0, 0},
    {TK_OPTION_STRING, "-tristatevalue", "tristateValue", "TristateValue",
	DEF_BUTTON_TRISTATE_VALUE, Tk_Offset(TkButton, tristateValuePtr), -1, 0, 0, 0},
    {TK_OPTION_INT, "-underline", "underline", "Underline",
	DEF_BUTTON_UNDERLINE, -1, Tk_Offset(TkButton, underline), 0, 0, 0},
    {TK_OPTION_STRING, "-value", "value", "Value",
	DEF_BUTTON_VALUE, Tk_Offset(TkButton, onValuePtr), -1, 0, 0, 0},
    {TK_OPTION_STRING, "-variable", "variable", "Variable",
	DEF_RADIOBUTTON_VARIABLE, Tk_Offset(TkButton, selVarNamePtr), -1,
	0, 0, 0},
    {TK_OPTION_STRING, "-width", "width", "Width",
	DEF_BUTTON_WIDTH, Tk_Offset(TkButton, widthPtr), -1, 0, 0, 0},
    {TK_OPTION_PIXELS, "-wraplength", "wrapLength", "WrapLength",
	DEF_BUTTON_WRAP_LENGTH, Tk_Offset(TkButton, wrapLengthPtr),
	Tk_Offset(TkButton, wrapLength), 0, 0, 0},
    {TK_OPTION_END, NULL, NULL, NULL, NULL, 0, -1, 0, 0, 0}
};

/*
 * The following table maps from one of the type values defined in tkButton.h,
 * such as TYPE_LABEL, to the option template for that class of widgets.
 */

static const Tk_OptionSpec *const optionSpecs[] = {
    labelOptionSpecs,
    buttonOptionSpecs,
    checkbuttonOptionSpecs,
    radiobuttonOptionSpecs
};

/*
 * The following tables define the widget commands supported by each of the
 * classes, and map the indexes into the string tables into a single
 * enumerated type used to dispatch the widget command.
 */

static const char *commandNames[][8] = {
    {"cget", "configure", NULL},
    {"cget", "configure", "flash", "invoke", NULL},
    {"cget", "configure", "deselect", "flash", "invoke", "select",
	    "toggle", NULL},
    {"cget", "configure", "deselect", "flash", "invoke", "select", NULL}
};
enum command {
    COMMAND_CGET, COMMAND_CONFIGURE, COMMAND_DESELECT, COMMAND_FLASH,
    COMMAND_INVOKE, COMMAND_SELECT, COMMAND_TOGGLE
};
static enum command map[][8] = {
    {COMMAND_CGET, COMMAND_CONFIGURE},
    {COMMAND_CGET, COMMAND_CONFIGURE, COMMAND_FLASH, COMMAND_INVOKE},
    {COMMAND_CGET, COMMAND_CONFIGURE, COMMAND_DESELECT, COMMAND_FLASH,
	    COMMAND_INVOKE, COMMAND_SELECT, COMMAND_TOGGLE},
    {COMMAND_CGET, COMMAND_CONFIGURE, COMMAND_DESELECT, COMMAND_FLASH,
	    COMMAND_INVOKE, COMMAND_SELECT}
};

/*
 * Forward declarations for functions defined later in this file:
 */

static void		ButtonCmdDeletedProc(ClientData clientData);
static int		ButtonCreate(ClientData clientData,
			    Tcl_Interp *interp, int objc,
			    Tcl_Obj *const objv[], int type);
static void		ButtonEventProc(ClientData clientData,
			    XEvent *eventPtr);
static void		ButtonImageProc(ClientData clientData,
			    int x, int y, int width, int height,
			    int imgWidth, int imgHeight);
static void		ButtonSelectImageProc(ClientData clientData,
			    int x, int y, int width, int height,
			    int imgWidth, int imgHeight);
static void		ButtonTristateImageProc(ClientData clientData,
			    int x, int y, int width, int height,
			    int imgWidth, int imgHeight);
static char *		ButtonTextVarProc(ClientData clientData,
			    Tcl_Interp *interp, const char *name1,
			    const char *name2, int flags);
static char *		ButtonVarProc(ClientData clientData,
			    Tcl_Interp *interp, const char *name1,
			    const char *name2, int flags);
static int		ButtonWidgetObjCmd(ClientData clientData,
			    Tcl_Interp *interp, int objc,
			    Tcl_Obj *const objv[]);
static int		ConfigureButton(Tcl_Interp *interp, TkButton *butPtr,
			    int objc, Tcl_Obj *const objv[]);
static void		DestroyButton(TkButton *butPtr);
static char *		strchrpnl(char *);

/*
 *--------------------------------------------------------------
 *
 * Tk_ButtonCmd, Tk_CheckbuttonCmd, Tk_LabelCmd, Tk_RadiobuttonCmd --
 *
 *	These functions are invoked to process the "button", "label",
 *	"radiobutton", and "checkbutton" Tcl commands. See the user
 *	documentation for details on what they do.
 *
 * Results:
 *	A standard Tcl result.
 *
 * Side effects:
 *	See the user documentation. These functions are just wrappers; they
 *	call ButtonCreate to do all of the real work.
 *
 *--------------------------------------------------------------
 */

int
Tk_ButtonObjCmd(
    ClientData clientData,	/* Either NULL or pointer to option table. */
    Tcl_Interp *interp,		/* Current interpreter. */
    int objc,			/* Number of arguments. */
    Tcl_Obj *const objv[])	/* Argument values. */
{
    return ButtonCreate(clientData, interp, objc, objv, TYPE_BUTTON);
}

int
Tk_CheckbuttonObjCmd(
    ClientData clientData,	/* Either NULL or pointer to option table. */
    Tcl_Interp *interp,		/* Current interpreter. */
    int objc,			/* Number of arguments. */
    Tcl_Obj *const objv[])	/* Argument values. */
{
    return ButtonCreate(clientData, interp, objc, objv, TYPE_CHECK_BUTTON);
}

int
Tk_LabelObjCmd(
    ClientData clientData,	/* Either NULL or pointer to option table. */
    Tcl_Interp *interp,		/* Current interpreter. */
    int objc,			/* Number of arguments. */
    Tcl_Obj *const objv[])	/* Argument values. */
{
    return ButtonCreate(clientData, interp, objc, objv, TYPE_LABEL);
}

int
Tk_RadiobuttonObjCmd(
    ClientData clientData,	/* Either NULL or pointer to option table. */
    Tcl_Interp *interp,		/* Current interpreter. */
    int objc,			/* Number of arguments. */
    Tcl_Obj *const objv[])	/* Argument values. */
{
    return ButtonCreate(clientData, interp, objc, objv, TYPE_RADIO_BUTTON);
}

/*
 *--------------------------------------------------------------
 *
 * ButtonCreate --
 *
 *	This function does all the real work of implementing the "button",
 *	"label", "radiobutton", and "checkbutton" Tcl commands. See the user
 *	documentation for details on what it does.
 *
 * Results:
 *	A standard Tcl result.
 *
 * Side effects:
 *	See the user documentation.
 *
 *--------------------------------------------------------------
 */

static int
ButtonCreate(
    ClientData clientData,	/* NULL. */
    Tcl_Interp *interp,		/* Current interpreter. */
    int objc,			/* Number of arguments. */
    Tcl_Obj *const objv[],	/* Argument values. */
    int type)			/* Type of button to create: TYPE_LABEL,
				 * TYPE_BUTTON, TYPE_CHECK_BUTTON, or
				 * TYPE_RADIO_BUTTON. */
{
    TkButton *butPtr;
    Tk_OptionTable optionTable;
    Tk_Window tkwin;
    ThreadSpecificData *tsdPtr = (ThreadSpecificData *)
	    Tcl_GetThreadData(&dataKey, sizeof(ThreadSpecificData));

    if (!tsdPtr->defaultsInitialized) {
	TkpButtonSetDefaults();
	tsdPtr->defaultsInitialized = 1;
    }

    if (objc < 2) {
	Tcl_WrongNumArgs(interp, 1, objv, "pathName ?options?");
	return TCL_ERROR;
    }

    /*
     * Create the new window.
     */

    tkwin = Tk_CreateWindowFromPath(interp, Tk_MainWindow(interp),
	    Tcl_GetString(objv[1]), NULL);
    if (tkwin == NULL) {
	return TCL_ERROR;
    }

    /*
     * Create the option table for this widget class. If it has already been
     * created, the cached pointer will be returned.
     */

    optionTable = Tk_CreateOptionTable(interp, optionSpecs[type]);

    Tk_SetClass(tkwin, classNames[type]);
    butPtr = TkpCreateButton(tkwin);

    Tk_SetClassProcs(tkwin, &tkpButtonProcs, (ClientData) butPtr);

    /*
     * Initialize the data structure for the button.
     */

    butPtr->tkwin = tkwin;
    butPtr->display = Tk_Display(tkwin);
    butPtr->interp = interp;
    butPtr->widgetCmd = Tcl_CreateObjCommand(interp, Tk_PathName(tkwin),
	    ButtonWidgetObjCmd, (ClientData) butPtr, ButtonCmdDeletedProc);
    butPtr->type = type;
    butPtr->optionTable = optionTable;
    butPtr->textPtr = NULL;
    butPtr->underline = -1;
    butPtr->textVarNamePtr = NULL;
    butPtr->bitmap = None;
    butPtr->imagePtr = NULL;
    butPtr->image = NULL;
    butPtr->selectImagePtr = NULL;
    butPtr->selectImage = NULL;
    butPtr->tristateImagePtr = NULL;
    butPtr->tristateImage = NULL;
    butPtr->state = STATE_NORMAL;
    butPtr->normalBorder = NULL;
    butPtr->activeBorder = NULL;
    butPtr->borderWidthPtr = NULL;
    butPtr->borderWidth = 0;
    butPtr->relief = TK_RELIEF_FLAT;
    butPtr->highlightWidthPtr = NULL;
    butPtr->highlightWidth = 0;
    butPtr->highlightBorder = NULL;
    butPtr->highlightColorPtr = NULL;
    butPtr->inset = 0;
    butPtr->tkfont = NULL;
    butPtr->tkaltfont = NULL;
    butPtr->tkscriptfont = NULL;
    butPtr->tkaltscriptfont = NULL;
    butPtr->normalFg = NULL;
    butPtr->activeFg = NULL;
    butPtr->disabledFg = NULL;
    butPtr->normalTextGC = None;
    butPtr->normalAltTextGC = None;
    butPtr->normalScriptTextGC = None;
    butPtr->normalAltScriptTextGC = None;
    butPtr->activeTextGC = None;
    butPtr->activeAltTextGC = None;
    butPtr->activeScriptTextGC = None;
    butPtr->activeAltScriptTextGC = None;
    butPtr->disabledGC = None;
    butPtr->disabledAltGC = None;
    butPtr->disabledScriptGC = None;
    butPtr->disabledAltScriptGC = None;
    butPtr->stippleGC = None;
    butPtr->gray = None;
    butPtr->copyGC = None;
    butPtr->widthPtr = NULL;
    butPtr->width = 0;
    butPtr->heightPtr = NULL;
    butPtr->height = 0;
    butPtr->wrapLengthPtr = NULL;
    butPtr->wrapLength = 0;
    butPtr->padXPtr = NULL;
    butPtr->padX = 0;
    butPtr->padYPtr = NULL;
    butPtr->padY = 0;
    butPtr->anchor = TK_ANCHOR_CENTER;
    butPtr->justify = TK_JUSTIFY_CENTER;
    butPtr->indicatorOn = 0;
    butPtr->selectBorder = NULL;
    butPtr->textWidth = 0;
    butPtr->textHeight = 0;
    butPtr->textLayout = NULL;
    butPtr->indicatorSpace = 0;
    butPtr->indicatorDiameter = 0;
    butPtr->defaultState = DEFAULT_DISABLED;
    butPtr->selVarNamePtr = NULL;
    butPtr->onValuePtr = NULL;
    butPtr->offValuePtr = NULL;
    butPtr->tristateValuePtr = NULL;
    butPtr->cursor = None;
    butPtr->takeFocusPtr = NULL;
    butPtr->commandPtr = NULL;
    butPtr->flags = 0;
    butPtr->rotation = 0;
    butPtr->bspecial = 0;

    Tk_CreateEventHandler(butPtr->tkwin,
	    ExposureMask|StructureNotifyMask|FocusChangeMask,
	    ButtonEventProc, (ClientData) butPtr);

    if (Tk_InitOptions(interp, (char *) butPtr, optionTable, tkwin)
	    != TCL_OK) {
	Tk_DestroyWindow(butPtr->tkwin);
	return TCL_ERROR;
    }
    if (ConfigureButton(interp, butPtr, objc - 2, objv + 2) != TCL_OK) {
	Tk_DestroyWindow(butPtr->tkwin);
	return TCL_ERROR;
    }

    Tcl_SetStringObj(Tcl_GetObjResult(interp), Tk_PathName(butPtr->tkwin),
	    -1);
    return TCL_OK;
}

/*
 *--------------------------------------------------------------
 *
 * ButtonWidgetCmd --
 *
 *	This function is invoked to process the Tcl command that corresponds
 *	to a widget managed by this module. See the user documentation for
 *	details on what it does.
 *
 * Results:
 *	A standard Tcl result.
 *
 * Side effects:
 *	See the user documentation.
 *
 *--------------------------------------------------------------
 */

static int
ButtonWidgetObjCmd(
    ClientData clientData,	/* Information about button widget. */
    Tcl_Interp *interp,		/* Current interpreter. */
    int objc,			/* Number of arguments. */
    Tcl_Obj *const objv[])	/* Argument values. */
{
    TkButton *butPtr = (TkButton *) clientData;
    int index;
    int result;
    Tcl_Obj *objPtr;

    if (objc < 2) {
        Tcl_WrongNumArgs(interp, 1, objv, "option ?arg arg ...?");
	return TCL_ERROR;
    }
    result = Tcl_GetIndexFromObj(interp, objv[1], commandNames[butPtr->type],
	    "option", 0, &index);
    if (result != TCL_OK) {
	return result;
    }
    Tcl_Preserve((ClientData) butPtr);

    switch (map[butPtr->type][index]) {
    case COMMAND_CGET:
	if (objc != 3) {
	    Tcl_WrongNumArgs(interp, 1, objv, "cget option");
	    goto error;
	}
	objPtr = Tk_GetOptionValue(interp, (char *) butPtr,
		butPtr->optionTable, objv[2], butPtr->tkwin);
	if (objPtr == NULL) {
	    goto error;
	} else {
	    Tcl_SetObjResult(interp, objPtr);
	}
	break;

    case COMMAND_CONFIGURE:
	if (objc <= 3) {
	    objPtr = Tk_GetOptionInfo(interp, (char *) butPtr,
		    butPtr->optionTable, (objc == 3) ? objv[2] : NULL,
		    butPtr->tkwin);
	    if (objPtr == NULL) {
		goto error;
	    } else {
		Tcl_SetObjResult(interp, objPtr);
	    }
	} else {
	    result = ConfigureButton(interp, butPtr, objc-2, objv+2);
	}
	break;

    case COMMAND_DESELECT:
	if (objc > 2) {
	    Tcl_WrongNumArgs(interp, 1, objv, "deselect");
	    goto error;
	}
	if (butPtr->type == TYPE_CHECK_BUTTON) {
	    if (Tcl_ObjSetVar2(interp, butPtr->selVarNamePtr, NULL,
		    butPtr->offValuePtr, TCL_GLOBAL_ONLY|TCL_LEAVE_ERR_MSG)
		    == NULL) {
		goto error;
	    }
	} else if (butPtr->flags & SELECTED) {
	    if (Tcl_ObjSetVar2(interp, butPtr->selVarNamePtr, NULL,
		     Tcl_NewObj(), TCL_GLOBAL_ONLY|TCL_LEAVE_ERR_MSG) == NULL){
		goto error;
	    }
	}
	break;

    case COMMAND_FLASH:
	if (objc > 2) {
	    Tcl_WrongNumArgs(interp, 1, objv, "flash");
	    goto error;
	}
	if (butPtr->state != STATE_DISABLED) {
	    int i;

	    for (i = 0; i < 4; i++) {
		if (butPtr->state == STATE_NORMAL) {
		    butPtr->state = STATE_ACTIVE;
		    Tk_SetBackgroundFromBorder(butPtr->tkwin,
			    butPtr->activeBorder);
		} else {
		    butPtr->state = STATE_NORMAL;
		    Tk_SetBackgroundFromBorder(butPtr->tkwin,
			    butPtr->normalBorder);
		}
		TkpDisplayButton((ClientData) butPtr);

		/*
		 * Special note: must cancel any existing idle handler for
		 * TkpDisplayButton; it's no longer needed, and
		 * TkpDisplayButton cleared the REDRAW_PENDING flag.
		 */

		Tcl_CancelIdleCall(TkpDisplayButton, (ClientData) butPtr);
		XFlush(butPtr->display);
		Tcl_Sleep(50);
	    }
	}
	break;

    case COMMAND_INVOKE:
	if (objc > 2) {
	    Tcl_WrongNumArgs(interp, 1, objv, "invoke");
	    goto error;
	}
	if (butPtr->state != STATE_DISABLED) {
	    result = TkInvokeButton(butPtr);
	}
	break;

    case COMMAND_SELECT:
	if (objc > 2) {
	    Tcl_WrongNumArgs(interp, 1, objv, "select");
	    goto error;
	}
	if (Tcl_ObjSetVar2(interp, butPtr->selVarNamePtr, NULL,
		butPtr->onValuePtr, TCL_GLOBAL_ONLY|TCL_LEAVE_ERR_MSG)
		== NULL) {
	    goto error;
	}
	break;

    case COMMAND_TOGGLE:
	if (objc > 2) {
	    Tcl_WrongNumArgs(interp, 1, objv, "toggle");
	    goto error;
	}
	if (Tcl_ObjSetVar2(interp, butPtr->selVarNamePtr, NULL,
		(butPtr->flags & SELECTED) ? butPtr->offValuePtr
		: butPtr->onValuePtr, TCL_GLOBAL_ONLY|TCL_LEAVE_ERR_MSG)
		== NULL) {
	    goto error;
	}
	break;
    }
    Tcl_Release((ClientData) butPtr);
    return result;

  error:
    Tcl_Release((ClientData) butPtr);
    return TCL_ERROR;
}

/*
 *----------------------------------------------------------------------
 *
 * DestroyButton --
 *
 *	This function is invoked by ButtonEventProc to free all the resources
 *	of a button and clean up its state.
 *
 * Results:
 *	None.
 *
 * Side effects:
 *	Everything associated with the widget is freed.
 *
 *----------------------------------------------------------------------
 */

static void
DestroyButton(
    TkButton *butPtr)		/* Info about button widget. */
{
    butPtr->flags |= BUTTON_DELETED;
    TkpDestroyButton(butPtr);

    if (butPtr->flags & REDRAW_PENDING) {
	Tcl_CancelIdleCall(TkpDisplayButton, (ClientData) butPtr);
    }

    /*
     * Free up all the stuff that requires special handling, then let
     * Tk_FreeOptions handle all the standard option-related stuff.
     */

    Tcl_DeleteCommandFromToken(butPtr->interp, butPtr->widgetCmd);
    if (butPtr->textVarNamePtr != NULL) {
	Tcl_UntraceVar(butPtr->interp, Tcl_GetString(butPtr->textVarNamePtr),
		TCL_GLOBAL_ONLY|TCL_TRACE_WRITES|TCL_TRACE_UNSETS,
		ButtonTextVarProc, (ClientData) butPtr);
    }
    if (butPtr->image != NULL) {
	Tk_FreeImage(butPtr->image);
    }
    if (butPtr->selectImage != NULL) {
	Tk_FreeImage(butPtr->selectImage);
    }
    if (butPtr->tristateImage != NULL) {
	Tk_FreeImage(butPtr->tristateImage);
    }
    if (butPtr->normalTextGC != None) {
	Tk_FreeGC(butPtr->display, butPtr->normalTextGC);
    }
    if (butPtr->normalAltTextGC != None) {
	Tk_FreeGC(butPtr->display, butPtr->normalAltTextGC);
    }
    if (butPtr->normalScriptTextGC != None) {
	Tk_FreeGC(butPtr->display, butPtr->normalScriptTextGC);
    }
    if (butPtr->normalAltScriptTextGC != None) {
	Tk_FreeGC(butPtr->display, butPtr->normalAltScriptTextGC);
    }
    if (butPtr->activeTextGC != None) {
	Tk_FreeGC(butPtr->display, butPtr->activeTextGC);
    }
    if (butPtr->activeAltTextGC != None) {
	Tk_FreeGC(butPtr->display, butPtr->activeAltTextGC);
    }
    if (butPtr->activeScriptTextGC != None) {
	Tk_FreeGC(butPtr->display, butPtr->activeScriptTextGC);
    }
    if (butPtr->activeAltScriptTextGC != None) {
	Tk_FreeGC(butPtr->display, butPtr->activeAltScriptTextGC);
    }
    if (butPtr->disabledGC != None) {
	Tk_FreeGC(butPtr->display, butPtr->disabledGC);
    }
    if (butPtr->disabledAltGC != None) {
	Tk_FreeGC(butPtr->display, butPtr->disabledAltGC);
    }
    if (butPtr->disabledScriptGC != None) {
	Tk_FreeGC(butPtr->display, butPtr->disabledScriptGC);
    }
    if (butPtr->disabledAltScriptGC != None) {
	Tk_FreeGC(butPtr->display, butPtr->disabledAltScriptGC);
    }
    if (butPtr->stippleGC != None) {
	Tk_FreeGC(butPtr->display, butPtr->stippleGC);
    }
    if (butPtr->gray != None) {
	Tk_FreeBitmap(butPtr->display, butPtr->gray);
    }
    if (butPtr->copyGC != None) {
	Tk_FreeGC(butPtr->display, butPtr->copyGC);
    }
    if (butPtr->textLayout != NULL) {
	Tk_FreeTextLayout(butPtr->textLayout);
    }
    if (butPtr->selVarNamePtr != NULL) {
	Tcl_UntraceVar(butPtr->interp, Tcl_GetString(butPtr->selVarNamePtr),
		TCL_GLOBAL_ONLY|TCL_TRACE_WRITES|TCL_TRACE_UNSETS,
		ButtonVarProc, (ClientData) butPtr);
    }
    Tk_FreeConfigOptions((char *) butPtr, butPtr->optionTable,
	    butPtr->tkwin);
    butPtr->tkwin = NULL;
    Tcl_EventuallyFree((ClientData) butPtr, TCL_DYNAMIC);
}

/*
 *----------------------------------------------------------------------
 *
 * ConfigureButton --
 *
 *	This function is called to process an objc/objv list to set
 *	configuration options for a button widget.
 *
 * Results:
 *	The return value is a standard Tcl result. If TCL_ERROR is returned,
 *	then an error message is left in interp's result.
 *
 * Side effects:
 *	Configuration information, such as text string, colors, font, etc. get
 *	set for butPtr; old resources get freed, if there were any. The button
 *	is redisplayed.
 *
 *----------------------------------------------------------------------
 */

static int
ConfigureButton(
    Tcl_Interp *interp,		/* Used for error reporting. */
    register TkButton *butPtr,	/* Information about widget;  may or may
				 * not already have values for some fields. */
    int objc,			/* Number of arguments. */
    Tcl_Obj *const objv[])	/* Argument values. */
{
    Tk_SavedOptions savedOptions;
    Tcl_Obj *errorResult = NULL;
    int error, haveImage;
    Tk_Image image;

    /*
     * Eliminate any existing trace on variables monitored by the button.
     */

    if (butPtr->textVarNamePtr != NULL) {
	Tcl_UntraceVar(interp, Tcl_GetString(butPtr->textVarNamePtr),
		TCL_GLOBAL_ONLY|TCL_TRACE_WRITES|TCL_TRACE_UNSETS,
		ButtonTextVarProc, (ClientData) butPtr);
    }
    if (butPtr->selVarNamePtr != NULL) {
	Tcl_UntraceVar(interp, Tcl_GetString(butPtr->selVarNamePtr),
		TCL_GLOBAL_ONLY|TCL_TRACE_WRITES|TCL_TRACE_UNSETS,
		ButtonVarProc, (ClientData) butPtr);
    }

    /*
     * The following loop is potentially executed twice. During the first pass
     * configuration options get set to their new values. If there is an error
     * in this pass, we execute a second pass to restore all the options to
     * their previous values.
     */

    for (error = 0; error <= 1; error++) {
	if (!error) {
	    /*
	     * First pass: set options to new values.
	     */

	    if (Tk_SetOptions(interp, (char *) butPtr,
		    butPtr->optionTable, objc, objv,
		    butPtr->tkwin, &savedOptions, NULL) != TCL_OK) {
		continue;
	    }
	} else {
	    /*
	     * Second pass: restore options to old values.
	     */

	    errorResult = Tcl_GetObjResult(interp);
	    Tcl_IncrRefCount(errorResult);
	    Tk_RestoreSavedOptions(&savedOptions);
	}

	if ((butPtr->flags & BUTTON_DELETED)) {
	    /*
	     * Somehow button was deleted - just abort now. [Bug #824479]
	     */
	    return TCL_ERROR;
	}

	/*
	 * A few options need special processing, such as setting the
	 * background from a 3-D border, or filling in complicated defaults
	 * that couldn't be specified to Tk_SetOptions.
	 */

	if ((butPtr->state == STATE_ACTIVE)
		&& !Tk_StrictMotif(butPtr->tkwin)) {
	    Tk_SetBackgroundFromBorder(butPtr->tkwin, butPtr->activeBorder);
	} else {
	    Tk_SetBackgroundFromBorder(butPtr->tkwin, butPtr->normalBorder);
	}
	if (butPtr->borderWidth < 0) {
	    butPtr->borderWidth = 0;
	}
	if (butPtr->highlightWidth < 0) {
	    butPtr->highlightWidth = 0;
	}
	if (butPtr->padX < 0) {
	    butPtr->padX = 0;
	}
	if (butPtr->padY < 0) {
	    butPtr->padY = 0;
	}

	if (butPtr->type >= TYPE_CHECK_BUTTON) {
	    Tcl_Obj *valuePtr, *namePtr;

	    if (butPtr->selVarNamePtr == NULL) {
		butPtr->selVarNamePtr = Tcl_NewStringObj(
			Tk_Name(butPtr->tkwin), -1);
		Tcl_IncrRefCount(butPtr->selVarNamePtr);
	    }
	    namePtr = butPtr->selVarNamePtr;

	    /*
	     * Select the button if the associated variable has the
	     * appropriate value, initialize the variable if it doesn't exist,
	     * then set a trace on the variable to monitor future changes to
	     * its value.
	     */

	    valuePtr = Tcl_ObjGetVar2(interp, namePtr, NULL, TCL_GLOBAL_ONLY);
	    butPtr->flags &= ~SELECTED;
            butPtr->flags &= ~TRISTATED;
	    if (valuePtr != NULL) {
		const char *value = Tcl_GetString(valuePtr);
		if (strcmp(value, Tcl_GetString(butPtr->onValuePtr)) == 0) {
		    butPtr->flags |= SELECTED;
		} else if (strcmp(value,
                        Tcl_GetString(butPtr->tristateValuePtr)) == 0) {
		    butPtr->flags |= TRISTATED;

		    /*
		     * For checkbuttons if the tristate value is the
		     * same as the offvalue then prefer off to tristate
		     */

		    if (butPtr->offValuePtr
			&& strcmp(value,
			    Tcl_GetString(butPtr->offValuePtr)) == 0) {
			butPtr->flags &= ~TRISTATED;
		    }
                }
	    } else {
		if (Tcl_ObjSetVar2(interp, namePtr, NULL,
			(butPtr->type == TYPE_CHECK_BUTTON)
			? butPtr->offValuePtr : Tcl_NewObj(),
			TCL_GLOBAL_ONLY|TCL_LEAVE_ERR_MSG)
			== NULL) {
		    continue;
		}

		/*
		 * If a radiobutton has the empty string as value it should be
		 * selected.
		 */

 		if ((butPtr->type == TYPE_RADIO_BUTTON) &&
			(*Tcl_GetString(butPtr->onValuePtr) == 0)) {
		    butPtr->flags |= SELECTED;
		}
	    }
	}

	/*
	 * Get the images for the widget, if there are any. Allocate the new
	 * images before freeing the old ones, so that the reference counts
	 * don't go to zero and cause image data to be discarded.
	 */

	if (butPtr->imagePtr != NULL) {
	    image = Tk_GetImage(butPtr->interp, butPtr->tkwin,
		    Tcl_GetString(butPtr->imagePtr), ButtonImageProc,
		    (ClientData) butPtr);
	    if (image == NULL) {
		continue;
	    }
	} else {
	    image = NULL;
	}
	if (butPtr->image != NULL) {
	    Tk_FreeImage(butPtr->image);
	}
	butPtr->image = image;
	if (butPtr->selectImagePtr != NULL) {
	    image = Tk_GetImage(butPtr->interp, butPtr->tkwin,
		    Tcl_GetString(butPtr->selectImagePtr),
		    ButtonSelectImageProc, (ClientData) butPtr);
	    if (image == NULL) {
		continue;
	    }
	} else {
	    image = NULL;
	}
	if (butPtr->selectImage != NULL) {
	    Tk_FreeImage(butPtr->selectImage);
	}
	butPtr->selectImage = image;
	if (butPtr->tristateImagePtr != NULL) {
	    image = Tk_GetImage(butPtr->interp, butPtr->tkwin,
		    Tcl_GetString(butPtr->tristateImagePtr),
		    ButtonTristateImageProc, (ClientData) butPtr);
	    if (image == NULL) {
		continue;
	    }
	} else {
	    image = NULL;
	}
	if (butPtr->tristateImage != NULL) {
	    Tk_FreeImage(butPtr->tristateImage);
	}
	butPtr->tristateImage = image;

	haveImage = 0;
	if (butPtr->imagePtr != NULL || butPtr->bitmap != None) {
	    haveImage = 1;
	}
	if ((!haveImage || butPtr->compound != COMPOUND_NONE)
		&& (butPtr->textVarNamePtr != NULL)) {
	    /*
	     * The button must display the value of a variable: set up a trace
	     * on the variable's value, create the variable if it doesn't
	     * exist, and fetch its current value.
	     */

	    Tcl_Obj *valuePtr, *namePtr;

	    namePtr = butPtr->textVarNamePtr;
	    valuePtr = Tcl_ObjGetVar2(interp, namePtr, NULL, TCL_GLOBAL_ONLY);
	    if (valuePtr == NULL) {
		if (Tcl_ObjSetVar2(interp, namePtr, NULL, butPtr->textPtr,
			TCL_GLOBAL_ONLY|TCL_LEAVE_ERR_MSG)
			== NULL) {
		    continue;
		}
	    } else {
		if (butPtr->textPtr != NULL) {
		    Tcl_DecrRefCount(butPtr->textPtr);
		}
		butPtr->textPtr = valuePtr;
		Tcl_IncrRefCount(butPtr->textPtr);
	    }
	}

	if ((butPtr->bitmap != None) || (butPtr->imagePtr != NULL)) {
	    /*
	     * The button must display the contents of an image or bitmap.
	     */

	    if (Tk_GetPixelsFromObj(interp, butPtr->tkwin, butPtr->widthPtr,
		    &butPtr->width) != TCL_OK) {
	    widthError:
		Tcl_AddErrorInfo(interp, "\n    (processing -width option)");
		continue;
	    }
	    if (Tk_GetPixelsFromObj(interp, butPtr->tkwin, butPtr->heightPtr,
		    &butPtr->height) != TCL_OK) {
	    heightError:
		Tcl_AddErrorInfo(interp, "\n    (processing -height option)");
		continue;
	    }
	} else {
	    /*
	     * The button displays an ordinary text string.
	     */

	    if (Tcl_GetIntFromObj(interp, butPtr->widthPtr, &butPtr->width)
		    != TCL_OK) {
		goto widthError;
	    }
	    if (Tcl_GetIntFromObj(interp, butPtr->heightPtr, &butPtr->height)
		    != TCL_OK) {
		goto heightError;
	    }
	}
	break;
    }
    if (!error) {
	Tk_FreeSavedOptions(&savedOptions);
    }

    /*
     * Reestablish the variable traces, if they're needed.
     */

    if (butPtr->textVarNamePtr != NULL) {
	Tcl_TraceVar(interp, Tcl_GetString(butPtr->textVarNamePtr),
		TCL_GLOBAL_ONLY|TCL_TRACE_WRITES|TCL_TRACE_UNSETS,
		ButtonTextVarProc, (ClientData) butPtr);
    }
    if (butPtr->selVarNamePtr != NULL) {
	Tcl_TraceVar(interp, Tcl_GetString(butPtr->selVarNamePtr),
		TCL_GLOBAL_ONLY|TCL_TRACE_WRITES|TCL_TRACE_UNSETS,
		ButtonVarProc, (ClientData) butPtr);
    }

    TkButtonWorldChanged((ClientData) butPtr);
    if (error) {
	Tcl_SetObjResult(interp, errorResult);
	Tcl_DecrRefCount(errorResult);
	return TCL_ERROR;
    } else {
	return TCL_OK;
    }
}

/*
 *---------------------------------------------------------------------------
 *
 * TkButtonWorldChanged --
 *
 *      This function is called when the world has changed in some way and the
 *      widget needs to recompute all its graphics contexts and determine its
 *      new geometry.
 *
 * Results:
 *      None.
 *
 * Side effects:
 *      Button will be relayed out and redisplayed.
 *
 *---------------------------------------------------------------------------
 */

void
TkButtonWorldChanged(
    ClientData instanceData)	/* Information about widget. */
{
    XGCValues gcValues;
    GC newGC;
    unsigned long mask;
    TkButton *butPtr;
    int balt,bscr;
    char *text;

    butPtr = (TkButton *) instanceData;

    /*
     * Recompute GCs.
     */

    balt = 0;
    bscr = 0;
    /*printf("text: [%s]\n",butPtr->text);*/
    if (butPtr->textPtr!=NULL) {
	text = Tcl_GetString(butPtr->textPtr);
	if (*text!=0) {
	    if (strstr(text,"`f")!=NULL)
		balt = 1;
	    if ((strstr(text,"`u")!=NULL)||(strstr(text,"`d")!=NULL))
		bscr = 1;
	}
    }
    butPtr->bspecial = balt || bscr;
    /*printf("bspecial set to %d\n", butPtr->bspecial);*/

    gcValues.font = Tk_FontId(butPtr->tkfont);
    gcValues.foreground = butPtr->normalFg->pixel;
    gcValues.background = Tk_3DBorderColor(butPtr->normalBorder)->pixel;

    /*
     * Note: GraphicsExpose events are disabled in normalTextGC because it's
     * used to copy stuff from an off-screen pixmap onto the screen (we know
     * that there's no problem with obscured areas).
     */

    gcValues.graphics_exposures = False;
    mask = GCForeground | GCBackground | GCFont | GCGraphicsExposures;
    newGC = Tk_GetGC(butPtr->tkwin, mask, &gcValues);
    if (butPtr->normalTextGC != None) {
	Tk_FreeGC(butPtr->display, butPtr->normalTextGC);
    }
    butPtr->normalTextGC = newGC;

    if (butPtr->normalAltTextGC != None) {
	Tk_FreeGC(butPtr->display, butPtr->normalAltTextGC);
    }
    if (balt) {
	gcValues.font = Tk_FontId(butPtr->tkaltfont);
	newGC = Tk_GetGC(butPtr->tkwin,
			 GCForeground|GCBackground|GCFont, &gcValues);
	butPtr->normalAltTextGC = newGC;
    } else {
	butPtr->normalAltTextGC = None;
    }

    if (butPtr->normalScriptTextGC != None) {
	Tk_FreeGC(butPtr->display, butPtr->normalScriptTextGC);
    }
    if (bscr) {
	gcValues.font = Tk_FontId(butPtr->tkscriptfont);
	newGC = Tk_GetGC(butPtr->tkwin,
			 GCForeground|GCBackground|GCFont, &gcValues);
	butPtr->normalScriptTextGC = newGC;
    } else {
	butPtr->normalScriptTextGC = None;
    }

    if (butPtr->normalAltScriptTextGC != None) {
	Tk_FreeGC(butPtr->display, butPtr->normalAltScriptTextGC);
    }
    if (balt && bscr) {
	gcValues.font = Tk_FontId(butPtr->tkaltscriptfont);
	newGC = Tk_GetGC(butPtr->tkwin,
			 GCForeground|GCBackground|GCFont, &gcValues);
	butPtr->normalAltScriptTextGC = newGC;
    } else {
	butPtr->normalAltScriptTextGC = None;
    }

    if (butPtr->activeFg != NULL) {
	gcValues.font = Tk_FontId(butPtr->tkfont);
	gcValues.foreground = butPtr->activeFg->pixel;
	gcValues.background = Tk_3DBorderColor(butPtr->activeBorder)->pixel;
	mask = GCForeground | GCBackground | GCFont;
	newGC = Tk_GetGC(butPtr->tkwin, mask, &gcValues);
	if (butPtr->activeTextGC != None) {
	    Tk_FreeGC(butPtr->display, butPtr->activeTextGC);
	}
	butPtr->activeTextGC = newGC;

	if (butPtr->activeAltTextGC != None) {
	    Tk_FreeGC(butPtr->display, butPtr->activeAltTextGC);
	}
	if (balt) {
	    gcValues.font = Tk_FontId(butPtr->tkaltfont);
	    newGC = Tk_GetGC(butPtr->tkwin,
			     GCForeground|GCBackground|GCFont, &gcValues);
	    butPtr->activeAltTextGC = newGC;
	} else {
	    butPtr->activeAltTextGC = None;
	}

	if (butPtr->activeScriptTextGC != None) {
	    Tk_FreeGC(butPtr->display, butPtr->activeScriptTextGC);
	}
	if (bscr) {
	    gcValues.font = Tk_FontId(butPtr->tkscriptfont);
	    newGC = Tk_GetGC(butPtr->tkwin,
			     GCForeground|GCBackground|GCFont, &gcValues);
	    butPtr->activeScriptTextGC = newGC;
	} else {
	    butPtr->activeScriptTextGC = None;
	}

	if (butPtr->activeAltScriptTextGC != None) {
	    Tk_FreeGC(butPtr->display, butPtr->activeAltScriptTextGC);
	}
	if (balt && bscr) {
	    gcValues.font = Tk_FontId(butPtr->tkaltscriptfont);
	    newGC = Tk_GetGC(butPtr->tkwin,
			     GCForeground|GCBackground|GCFont, &gcValues);
	    butPtr->activeAltScriptTextGC = newGC;
	} else {
	    butPtr->activeAltScriptTextGC = None;
	}
    }

    gcValues.background = Tk_3DBorderColor(butPtr->normalBorder)->pixel;

    /*
     * Create the GC that can be used for stippling
     */

    if (butPtr->stippleGC == None) {
	gcValues.foreground = gcValues.background;
	mask = GCForeground;
	if (butPtr->gray == None) {
	    butPtr->gray = Tk_GetBitmap(NULL, butPtr->tkwin, "gray50");
	}
	if (butPtr->gray != None) {
	    gcValues.fill_style = FillStippled;
	    gcValues.stipple = butPtr->gray;
	    mask |= GCFillStyle | GCStipple;
	}
	butPtr->stippleGC = Tk_GetGC(butPtr->tkwin, mask, &gcValues);
    }

    /*
     * Allocate the disabled graphics context, for drawing text in its
     * disabled state.
     */

    mask = GCForeground | GCBackground | GCFont;
    if (butPtr->disabledFg != NULL) {
	gcValues.foreground = butPtr->disabledFg->pixel;
    } else {
	gcValues.foreground = gcValues.background;
    }
    newGC = Tk_GetGC(butPtr->tkwin, mask, &gcValues);
    if (butPtr->disabledGC != None) {
	Tk_FreeGC(butPtr->display, butPtr->disabledGC);
    }
    butPtr->disabledGC = newGC;

    if (butPtr->disabledAltGC != None) {
	Tk_FreeGC(butPtr->display, butPtr->disabledAltGC);
    }
    if (balt) {
	gcValues.font = Tk_FontId(butPtr->tkaltfont);
	newGC = Tk_GetGC(butPtr->tkwin,
			 GCForeground|GCBackground|GCFont, &gcValues);
	butPtr->disabledAltGC = newGC;
    } else {
	butPtr->disabledAltGC = None;
    }

    if (butPtr->disabledScriptGC != None) {
	Tk_FreeGC(butPtr->display, butPtr->disabledScriptGC);
    }
    if (bscr) {
	gcValues.font = Tk_FontId(butPtr->tkscriptfont);
	newGC = Tk_GetGC(butPtr->tkwin,
			 GCForeground|GCBackground|GCFont, &gcValues);
	butPtr->disabledScriptGC = newGC;
    } else {
	butPtr->disabledScriptGC = None;
    }

    if (butPtr->disabledAltScriptGC != None) {
	Tk_FreeGC(butPtr->display, butPtr->disabledAltScriptGC);
    }
    if (balt && bscr) {
	gcValues.font = Tk_FontId(butPtr->tkaltscriptfont);
	newGC = Tk_GetGC(butPtr->tkwin,
			 GCForeground|GCBackground|GCFont, &gcValues);
	butPtr->disabledAltScriptGC = newGC;
    } else {
	butPtr->disabledAltScriptGC = None;
    }

    if (butPtr->copyGC == None) {
	butPtr->copyGC = Tk_GetGC(butPtr->tkwin, 0, &gcValues);
    }

    TkpComputeButtonGeometry(butPtr);

    /*
     * Lastly, arrange for the button to be redisplayed.
     */

    if (Tk_IsMapped(butPtr->tkwin) && !(butPtr->flags & REDRAW_PENDING)) {
	Tcl_DoWhenIdle(TkpDisplayButton, (ClientData) butPtr);
	butPtr->flags |= REDRAW_PENDING;
    }
}

/*
 *--------------------------------------------------------------
 *
 * ButtonEventProc --
 *
 *	This function is invoked by the Tk dispatcher for various events on
 *	buttons.
 *
 * Results:
 *	None.
 *
 * Side effects:
 *	When the window gets deleted, internal structures get cleaned up. When
 *	it gets exposed, it is redisplayed.
 *
 *--------------------------------------------------------------
 */

static void
ButtonEventProc(
    ClientData clientData,	/* Information about window. */
    XEvent *eventPtr)		/* Information about event. */
{
    TkButton *butPtr = (TkButton *) clientData;
    if ((eventPtr->type == Expose) && (eventPtr->xexpose.count == 0)) {
	goto redraw;
    } else if (eventPtr->type == ConfigureNotify) {
	/*
	 * Must redraw after size changes, since layout could have changed and
	 * borders will need to be redrawn.
	 */

	goto redraw;
    } else if (eventPtr->type == DestroyNotify) {
	DestroyButton(butPtr);
    } else if (eventPtr->type == FocusIn) {
	if (eventPtr->xfocus.detail != NotifyInferior) {
	    butPtr->flags |= GOT_FOCUS;
	    if (butPtr->highlightWidth > 0) {
		goto redraw;
	    }
	}
    } else if (eventPtr->type == FocusOut) {
	if (eventPtr->xfocus.detail != NotifyInferior) {
	    butPtr->flags &= ~GOT_FOCUS;
	    if (butPtr->highlightWidth > 0) {
		goto redraw;
	    }
	}
    }
    return;

  redraw:
    if ((butPtr->tkwin != NULL) && !(butPtr->flags & REDRAW_PENDING)) {
	Tcl_DoWhenIdle(TkpDisplayButton, (ClientData) butPtr);
	butPtr->flags |= REDRAW_PENDING;
    }
}

/*
 *----------------------------------------------------------------------
 *
 * ButtonCmdDeletedProc --
 *
 *	This function is invoked when a widget command is deleted. If the
 *	widget isn't already in the process of being destroyed, this command
 *	destroys it.
 *
 * Results:
 *	None.
 *
 * Side effects:
 *	The widget is destroyed.
 *
 *----------------------------------------------------------------------
 */

static void
ButtonCmdDeletedProc(
    ClientData clientData)	/* Pointer to widget record for widget. */
{
    TkButton *butPtr = (TkButton *) clientData;

    /*
     * This function could be invoked either because the window was destroyed
     * and the command was then deleted or because the command was deleted,
     * and then this function destroys the widget. The BUTTON_DELETED flag
     * distinguishes these cases.
     */

    if (!(butPtr->flags & BUTTON_DELETED)) {
	Tk_DestroyWindow(butPtr->tkwin);
    }
}

/*
 *----------------------------------------------------------------------
 *
 * TkInvokeButton --
 *
 *	This function is called to carry out the actions associated with a
 *	button, such as invoking a Tcl command or setting a variable. This
 *	function is invoked, for example, when the button is invoked via the
 *	mouse.
 *
 * Results:
 *	A standard Tcl return value. Information is also left in the interp's
 *	result.
 *
 * Side effects:
 *	Depends on the button and its associated command.
 *
 *----------------------------------------------------------------------
 */

int
TkInvokeButton(
    TkButton *butPtr)		/* Information about button. */
{
    Tcl_Obj *namePtr = butPtr->selVarNamePtr;

    if (butPtr->type == TYPE_CHECK_BUTTON) {
	if (butPtr->flags & SELECTED) {
	    if (Tcl_ObjSetVar2(butPtr->interp, namePtr, NULL,
		    butPtr->offValuePtr, TCL_GLOBAL_ONLY|TCL_LEAVE_ERR_MSG)
		    == NULL) {
		return TCL_ERROR;
	    }
	} else {
	    if (Tcl_ObjSetVar2(butPtr->interp, namePtr, NULL,
		    butPtr->onValuePtr, TCL_GLOBAL_ONLY|TCL_LEAVE_ERR_MSG)
		    == NULL) {
		return TCL_ERROR;
	    }
	}
    } else if (butPtr->type == TYPE_RADIO_BUTTON) {
	if (Tcl_ObjSetVar2(butPtr->interp, namePtr, NULL, butPtr->onValuePtr,
		TCL_GLOBAL_ONLY|TCL_LEAVE_ERR_MSG)
		== NULL) {
	    return TCL_ERROR;
	}
    }
    if ((butPtr->type != TYPE_LABEL) && (butPtr->commandPtr != NULL)) {
	return Tcl_EvalObjEx(butPtr->interp, butPtr->commandPtr,
		TCL_EVAL_GLOBAL);
    }
    return TCL_OK;
}

/*
 *--------------------------------------------------------------
 *
 * ButtonVarProc --
 *
 *	This function is invoked when someone changes the state variable
 *	associated with a radio button. Depending on the new value of the
 *	button's variable, the button may be selected or deselected.
 *
 * Results:
 *	NULL is always returned.
 *
 * Side effects:
 *	The button may become selected or deselected.
 *
 *--------------------------------------------------------------
 */

	/* ARGSUSED */
static char *
ButtonVarProc(
    ClientData clientData,	/* Information about button. */
    Tcl_Interp *interp,		/* Interpreter containing variable. */
    const char *name1,		/* Name of variable. */
    const char *name2,		/* Second part of variable name. */
    int flags)			/* Information about what happened. */
{
    register TkButton *butPtr = (TkButton *) clientData;
    char *name, *value;
    Tcl_Obj *valuePtr;

    name = Tcl_GetString(butPtr->selVarNamePtr);

    /*
     * If the variable is being unset, then just re-establish the trace unless
     * the whole interpreter is going away.
     */

    if (flags & TCL_TRACE_UNSETS) {
	butPtr->flags &= ~SELECTED;
        butPtr->flags &= ~TRISTATED;
	if ((flags & TCL_TRACE_DESTROYED) && !(flags & TCL_INTERP_DESTROYED)) {
	    Tcl_TraceVar(interp, name,
		    TCL_GLOBAL_ONLY|TCL_TRACE_WRITES|TCL_TRACE_UNSETS,
		    ButtonVarProc, clientData);
	}
	goto redisplay;
    }

    /*
     * Use the value of the variable to update the selected status of the
     * button.
     */

    valuePtr = Tcl_GetVar2Ex(interp, name, NULL, TCL_GLOBAL_ONLY);
    if (valuePtr == NULL) {
	value = Tcl_GetString(butPtr->tristateValuePtr);
    } else {
	value = Tcl_GetString(valuePtr);
    }
    if (strcmp(value, Tcl_GetString(butPtr->onValuePtr)) == 0) {
	if (butPtr->flags & SELECTED) {
	    return NULL;
	}
	butPtr->flags |= SELECTED;
        butPtr->flags &= ~TRISTATED;
    } else if (butPtr->offValuePtr 
	&& strcmp(value, Tcl_GetString(butPtr->offValuePtr)) == 0) {
	if (!(butPtr->flags & (SELECTED | TRISTATED))) {
	    return NULL;
	}
	butPtr->flags &= ~(SELECTED | TRISTATED);
    } else if (strcmp(value, Tcl_GetString(butPtr->tristateValuePtr)) == 0) {
        if (butPtr->flags & TRISTATED) {
            return NULL;
        }
        butPtr->flags |= TRISTATED;
        butPtr->flags &= ~SELECTED;
    } else if (butPtr->flags & (SELECTED | TRISTATED)) {
	butPtr->flags &= ~(SELECTED | TRISTATED);
    } else {
	return NULL;
    }

  redisplay:
    if ((butPtr->tkwin != NULL) && Tk_IsMapped(butPtr->tkwin)
	    && !(butPtr->flags & REDRAW_PENDING)) {
	Tcl_DoWhenIdle(TkpDisplayButton, (ClientData) butPtr);
	butPtr->flags |= REDRAW_PENDING;
    }
    return NULL;
}

/*
 *--------------------------------------------------------------
 *
 * ButtonTextVarProc --
 *
 *	This function is invoked when someone changes the variable whose
 *	contents are to be displayed in a button.
 *
 * Results:
 *	NULL is always returned.
 *
 * Side effects:
 *	The text displayed in the button will change to match the variable.
 *
 *--------------------------------------------------------------
 */

	/* ARGSUSED */
static char *
ButtonTextVarProc(
    ClientData clientData,	/* Information about button. */
    Tcl_Interp *interp,		/* Interpreter containing variable. */
    const char *name1,		/* Not used. */
    const char *name2,		/* Not used. */
    int flags)			/* Information about what happened. */
{
    TkButton *butPtr = (TkButton *) clientData;
    char *name;
    Tcl_Obj *valuePtr;

    if (butPtr->flags & BUTTON_DELETED) {
	return NULL;
    }

    name = Tcl_GetString(butPtr->textVarNamePtr);

    /*
     * If the variable is unset, then immediately recreate it unless the whole
     * interpreter is going away.
     */

    if (flags & TCL_TRACE_UNSETS) {
	if ((flags & TCL_TRACE_DESTROYED) && !(flags & TCL_INTERP_DESTROYED)) {
	    Tcl_SetVar2Ex(interp, name, NULL, butPtr->textPtr,
		    TCL_GLOBAL_ONLY);
	    Tcl_TraceVar(interp, name,
		    TCL_GLOBAL_ONLY|TCL_TRACE_WRITES|TCL_TRACE_UNSETS,
		    ButtonTextVarProc, clientData);
	}
	return NULL;
    }

    valuePtr = Tcl_GetVar2Ex(interp, name, NULL, TCL_GLOBAL_ONLY);
    if (valuePtr == NULL) {
	valuePtr = Tcl_NewObj();
    }
    Tcl_DecrRefCount(butPtr->textPtr);
    butPtr->textPtr = valuePtr;
    Tcl_IncrRefCount(butPtr->textPtr);
    TkpComputeButtonGeometry(butPtr);

    if ((butPtr->tkwin != NULL) && Tk_IsMapped(butPtr->tkwin)
	    && !(butPtr->flags & REDRAW_PENDING)) {
	Tcl_DoWhenIdle(TkpDisplayButton, (ClientData) butPtr);
	butPtr->flags |= REDRAW_PENDING;
    }
    return NULL;
}

/*
 *----------------------------------------------------------------------
 *
 * ButtonImageProc --
 *
 *	This function is invoked by the image code whenever the manager for an
 *	image does something that affects the size or contents of an image
 *	displayed in a button.
 *
 * Results:
 *	None.
 *
 * Side effects:
 *	Arranges for the button to get redisplayed.
 *
 *----------------------------------------------------------------------
 */

static void
ButtonImageProc(
    ClientData clientData,	/* Pointer to widget record. */
    int x, int y,		/* Upper left pixel (within image) that must
				 * be redisplayed. */
    int width, int height,	/* Dimensions of area to redisplay (might be
				 * <= 0). */
    int imgWidth, int imgHeight)/* New dimensions of image. */
{
    register TkButton *butPtr = (TkButton *) clientData;

    if (butPtr->tkwin != NULL) {
	TkpComputeButtonGeometry(butPtr);
	if (Tk_IsMapped(butPtr->tkwin) && !(butPtr->flags & REDRAW_PENDING)) {
	    Tcl_DoWhenIdle(TkpDisplayButton, (ClientData) butPtr);
	    butPtr->flags |= REDRAW_PENDING;
	}
    }
}

/*
 *----------------------------------------------------------------------
 *
 * ButtonSelectImageProc --
 *
 *	This function is invoked by the image code whenever the manager for an
 *	image does something that affects the size or contents of the image
 *	displayed in a button when it is selected.
 *
 * Results:
 *	None.
 *
 * Side effects:
 *	May arrange for the button to get redisplayed.
 *
 *----------------------------------------------------------------------
 */

static void
ButtonSelectImageProc(
    ClientData clientData,	/* Pointer to widget record. */
    int x, int y,		/* Upper left pixel (within image) that must
				 * be redisplayed. */
    int width, int height,	/* Dimensions of area to redisplay (might be
				 * <= 0). */
    int imgWidth, int imgHeight)/* New dimensions of image. */
{
    register TkButton *butPtr = (TkButton *) clientData;

    /*
     * Don't recompute geometry: it's controlled by the primary image.
     */

    if ((butPtr->flags & SELECTED) && (butPtr->tkwin != NULL)
	    && Tk_IsMapped(butPtr->tkwin)
	    && !(butPtr->flags & REDRAW_PENDING)) {
	Tcl_DoWhenIdle(TkpDisplayButton, (ClientData) butPtr);
	butPtr->flags |= REDRAW_PENDING;
    }
}


/*
 *----------------------------------------------------------------------
 *
 * TkKGDisplayText --
 *
 *
 * Results:
 *	None.
 *
 * Side effects:
 *	I don't know. :-)
 *
 *----------------------------------------------------------------------
 */

void
TkKGDisplayText(display, drawable, butPtr, tkfont, tkaltfont, tkscriptfont,
                tkaltscriptfont, string, numChars, x, y, length, justify, underline, gc,
                altGc, scriptGc, altScriptGc, bIsMeas)
     Display *display;
     Drawable drawable;
     TkButton *butPtr;
     Tk_Font tkfont, tkaltfont, tkscriptfont, tkaltscriptfont;
     char *string;
     int numChars;
     int x, y;
     int length;
     Tk_Justify justify;
     int underline;
     GC gc, altGc, scriptGc, altScriptGc;
     int bIsMeas;
{
#define MAXFONTS 50
  char *p, *pnl, *pChars[MAXFONTS], *pnewl;
  int charsThisLine, lengthThisLine, xThisLine;
  GC currGc[MAXFONTS];
  int xStart[MAXFONTS], yScript[MAXFONTS], nChars[MAXFONTS];
  int iel, i, currNumChars;
  int bInScript, bInAltFont, x0, y0, w, h, j;
  Tk_Font pTkFont[MAXFONTS];
  Tk_FontMetrics tkfm,tkscriptfm;
  Tk_Window tkwin = butPtr->tkwin;
  Pixmap pixmap = 0;
  XImage *ims, *imd;
  char *imddat;
  unsigned long maxWidth=0, numLines=0, bSupScr=0, bSubScr=0;
  int dysc;
  int st, supscrpad=0;
    
  /*printf("TkKGDisplayText %d,%d\n",x,y);*/

  /*printf("text: [%s]\n", string);*/
  Tk_GetFontMetrics(tkfont, &tkfm);
  Tk_GetFontMetrics(tkscriptfont, &tkscriptfm);
  dysc = (tkscriptfm.descent-tkscriptfm.ascent)/2;
  if (length <= 0)
	length = INT_MAX;
  x0 = x; y0 = y; w = Tk_Width(tkwin)-2*x0; h = Tk_Height(tkwin)-2*y0;
  if (!bIsMeas) {
	if ((butPtr->rotation==90)||(butPtr->rotation==270)) {
      pixmap = Tk_GetPixmap(butPtr->display, Tk_WindowId(tkwin),
                            Tk_Height(tkwin), Tk_Width(tkwin), Tk_Depth(tkwin));
	  Tk_Fill3DRectangle(tkwin, pixmap, butPtr->normalBorder, x0, y0, h,
                         w, 0, TK_RELIEF_FLAT);
	} else {
      pixmap = drawable;
    }
  }

  p = string;
  pnl = strchrpnl(p);
  currNumChars = pnl - p;
  y += tkfm.ascent;
  currGc[0] = gc;
  pTkFont[0] = tkfont;
  yScript[0] = 0;
  iel = 0;
  while (numChars > 0) {
	/*if (iel>0) {
	  yScript[0] = yScript[iel-1];
	  pTkFont[0] = pTkFont[iel-1];
	  currGc[0] = currGc[iel-1];
	  }*/
	yScript[0] = 0;
	pTkFont[0] = tkfont;
	currGc[0] = gc;
	iel = 0;
	xStart[0] = 0;
	pChars[0] = p;
	charsThisLine = 0;
	lengthThisLine = 1;
	bInScript = 0;
	bInAltFont = 0;
	bSubScr = 0;
	while ((lengthThisLine > 0)&&(numChars > 0)) {
      while (currNumChars <= 0) {
		if (p[0] == '`') {
          switch (p[1]) {
          case 'n':
			pTkFont[iel] = tkfont;
			currGc[iel] = gc;
			yScript[iel] = 0;
			bInScript = 0;
			bInAltFont = 0;
			break;
          case 'f':
			bInAltFont = 1;
			if (bInScript) {
              pTkFont[iel] = tkaltscriptfont;
              currGc[iel] = altScriptGc;
			} else {
              pTkFont[iel] = tkaltfont;
              currGc[iel] = altGc;
			}
			break;
          case 'd':
			bSubScr = 1;
          case 'u':
			if (bInAltFont) {
              pTkFont[iel] = tkaltscriptfont;
              currGc[iel] = altScriptGc;
			}else {
              pTkFont[iel] = tkscriptfont;
              currGc[iel] = scriptGc;
			}
			yScript[iel] = dysc;
			if (p[1]=='u') {
              yScript[iel] += (int)(tkfm.ascent*0.8);
              bSupScr = 1;
			}
			bInScript = 1;
			break;
          }
          if (p[1]=='`') {
			p++; numChars--;
			pnl = strchrpnl(p+1);
          } else {
			p += 2; numChars -= 2;
			if (numChars <= 0)
              break;
			pnl = strchrpnl(p);
          }
		}
		pChars[iel] = p;
		currNumChars = pnl - p;
	    }
      /*printf("a1 %d,%d,%d,%d\n",currNumChars,numChars,xStart[iel],length);*/
      pnewl = strchr(p,'\n');
      if (pnewl!=NULL) {
		if (pnewl==p) {
          p++; numChars--; currNumChars--;
          break;
		} else if (pnewl-p<currNumChars)
          currNumChars = pnewl-p;
      }
      /*printf("a2 %d,%d,%d,%d\n",currNumChars,numChars,xStart[iel],length);*/
      nChars[iel] = Tk_MeasureChars(pTkFont[iel], p, currNumChars,
                                    length, TK_WHOLE_WORDS/*|TK_AT_LEAST_ONE*/,
                                    &lengthThisLine);
      charsThisLine += nChars[iel];
      p += nChars[iel];
      pChars[iel+1] = p;
      lengthThisLine += xStart[iel];
      xStart[iel+1] = lengthThisLine;
      currNumChars -= nChars[iel];
      numChars -= nChars[iel];
      pTkFont[iel+1] = pTkFont[iel];
      currGc[iel+1] = currGc[iel];
      yScript[iel+1] = yScript[iel];
      /*printf("b %d, %d, %d, %d, %d\n",iel,currNumChars,numChars,lengthThisLine,nChars[iel]);*/
      iel++;
	}
	/*printf("x %d %d %d %d\n", x, length, lengthThisLine, numChars);*/
	if ((numLines==0)&&bSupScr) {
      st = dysc+(int)(tkfm.ascent*0.8)+tkscriptfm.ascent;
      if (st > tkfm.ascent) {
        /*printf("st:%d ascent:%d y:%d\n", st, font->ascent,y);*/
		supscrpad = st - tkfm.ascent;
		y += supscrpad;
        /*printf("y:%d\n", y);*/
      }
	}
	numLines++;
	if (bIsMeas) {
      if (lengthThisLine > maxWidth)
		maxWidth = lengthThisLine;
	} else {
      if (justify == TK_JUSTIFY_LEFT) {
		xThisLine = x;
      } else if (justify == TK_JUSTIFY_CENTER) {
		xThisLine = x + (length - lengthThisLine)/2;
      } else {
		xThisLine = x + (length - lengthThisLine);
      }
      for (i=0; i<iel; i++) {
        /*printf("[%s]%d,%d\n",pChars[i],nChars[i],xThisLine+xStart[i]);*/
		Tk_DrawChars(display, pixmap, currGc[i], pTkFont[i], pChars[i],
                     nChars[i], xThisLine+xStart[i], y-yScript[i]);
        /*if (nChars[i]==0) return;*/
      }
      /*
	    if ((underline >= 0) && (underline < charsThisLine)) {
		TkUnderlineChars(display, pixmap, gc, tkfont, p,
        xThisLine, y, xThisLine, 0, underline, underline);
	    }
	    underline -= charsThisLine;
      */
      y += tkfm.ascent + tkfm.descent;
	}

	/*if (isspace(UCHAR(*p))) {
      p++;
      numChars--;
      currNumChars--;
      underline--;
      }*/
  }

  if (bIsMeas) {
	butPtr->textWidth = maxWidth;
	butPtr->textHeight = (numLines==0?1:numLines) * (tkfm.ascent + tkfm.descent) + supscrpad;
    /*printf("maxw&h:%d %d\n", maxWidth, butPtr->textHeight);*/
	if (bSubScr) {
      st = -dysc+tkscriptfm.descent;
      if (st > tkfm.descent)
		butPtr->textHeight += st - tkfm.descent;
	}
	return;
  }

  if ((butPtr->rotation==90)||(butPtr->rotation==180)||(butPtr->rotation==270)) {
	if (butPtr->rotation==180)
      ims = XGetImage(butPtr->display, pixmap, x0, y0, w, h, 
                      AllPlanes, ZPixmap);
	else
      ims = XGetImage(butPtr->display, pixmap, x0, y0, h, w,
                      AllPlanes, ZPixmap);
    /*printf("ims: %d %d %d %d %d\n", h, w, ims->xoffset, ims->bitmap_pad, ims->bytes_per_line);*/
	imddat = malloc(w*h*sizeof(int));
	imd = XCreateImage(butPtr->display, Tk_Visual(tkwin), Tk_Depth(tkwin),
                       ZPixmap, 0, imddat, w, h, ims->bitmap_pad, 0);
	if (butPtr->rotation==90)
      for (i=0; i<w; i++)
		for (j=0; j<h; j++)
          XPutPixel(imd, i, h-j-1, XGetPixel(ims, j, i));
	else if (butPtr->rotation==180)
      for (i=0; i<h; i++)
		for (j=0; j<w; j++)
          XPutPixel(imd, w-j-1, h-i-1, XGetPixel(ims, j, i));
	else if (butPtr->rotation==270)
      for (i=0; i<w; i++)
		for (j=0; j<h; j++)
          XPutPixel(imd, w-i-1, j, XGetPixel(ims, j, i));
	TkPutImage(NULL, 0, butPtr->display, drawable, butPtr->copyGC, imd, 0, 0, x0,
               y0, w, h);
	if (butPtr->rotation!=180)
      Tk_FreePixmap(butPtr->display, pixmap);
	XDestroyImage(ims);
	XDestroyImage(imd);
  }
}

static char *strchrpnl(p)
     char *p;
{
  char *x=p;

  while ((x = strchr(x,'`'))) {
	if ((x[1]=='f')||(x[1]=='n')||(x[1]=='d')||(x[1]=='u')||(x[1]=='`'))
      return x;
	x += 2;
  }
  return p + strlen(p);
}

/*
 *----------------------------------------------------------------------
 *
 * ButtonTristateImageProc --
 *
 *	This function is invoked by the image code whenever the manager for an
 *	image does something that affects the size or contents of the image
 *	displayed in a button when it is selected.
 *
 * Results:
 *	None.
 *
 * Side effects:
 *	May arrange for the button to get redisplayed.
 *
 *----------------------------------------------------------------------
 */

static void
ButtonTristateImageProc(
    ClientData clientData,	/* Pointer to widget record. */
    int x, int y,		/* Upper left pixel (within image) that must
				 * be redisplayed. */
    int width, int height,	/* Dimensions of area to redisplay (might be
				 * <= 0). */
    int imgWidth, int imgHeight)/* New dimensions of image. */
{
    register TkButton *butPtr = (TkButton *) clientData;

    /*
     * Don't recompute geometry: it's controlled by the primary image.
     */

    if ((butPtr->flags & TRISTATED) && (butPtr->tkwin != NULL)
	    && Tk_IsMapped(butPtr->tkwin)
	    && !(butPtr->flags & REDRAW_PENDING)) {
	Tcl_DoWhenIdle(TkpDisplayButton, (ClientData) butPtr);
	butPtr->flags |= REDRAW_PENDING;
    }
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 4
 * fill-column: 78
 * End:
 */
