#!/bin/sh
symbols=" \
	TclCreateInterp Tcl TclArg TclArg1 TclSetResult \
	TkEventLoop TkQuitEventLoop \
	TkCreateFileHandler TkDeleteFileHandler \
	TkDefineBitmap TkPhotoPutBlock\\$ TkPhotoPutZoomedBlock\\$ \
	TkCanvasCreateItem TkCanvasGetTypes TkCanvasPointer \
	TkCanvasLastItem TkCanvasFindEnclosed TkCanvasCreateItemDirect \
	TkCanvTextDebug \
	CanvasClipLine CanvasSymbol CanvasSymbolDirect \
	Canvas3DClipTriangle Canvas3DLightTriangle Canvas3DProjection \
"

for symbol in ${symbols}; do
    grep "${symbol}" Packages/*.n
done

# End of File
