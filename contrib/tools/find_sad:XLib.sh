#!/bin/sh
symbols=" \
	XServerVendor XFlush \
	XRaiseWindow XLowerWindow XRestackWindows \
	XInternAtom \
	XSetTextProperty XGetTextProperty XChangeProperty XDeleteProperty \
	XGetWindowAttributes XGetWindowProperty \
	XQueryTree XListFonts XGetFontPath \
"

for symbol in ${symbols}; do
    grep ${symbol} Packages/*.n
done

# End of File
