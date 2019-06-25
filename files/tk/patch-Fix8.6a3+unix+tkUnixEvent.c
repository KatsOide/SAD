--- unix/tkUnixEvent.c.orig	2008-10-06 03:21:58.000000000 +0900
+++ unix/tkUnixEvent.c	2010-04-23 10:53:53.000000000 +0900
@@ -291,7 +291,7 @@
 	if (event.type == GenericEvent) {
 	    xGenericEvent *xgePtr = (xGenericEvent *) &event;
 
-	    Tcl_Panic("Wild GenericEvent; panic! (extension=%d,evtype=%d)"
+	    Tcl_Panic("Wild GenericEvent; panic! (extension=%d,evtype=%d)",
 		    xgePtr->extension, xgePtr->evtype);
 	}
 #endif
