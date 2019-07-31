--- leapsecs.h.orig	Wed Oct 14 01:51:59 1998
+++ leapsecs.h	Fri Nov 11 20:03:33 2005
@@ -1,10 +1,12 @@
 #ifndef LEAPSECS_H
 #define LEAPSECS_H
 
-extern int leapsecs_init();
-extern int leapsecs_read();
+#include "tai.h"
 
-extern void leapsecs_add();
-extern int leapsecs_sub();
+extern int leapsecs_init(void);
+extern int leapsecs_read(void);
+
+extern void leapsecs_add(struct tai *, int);
+extern int leapsecs_sub(struct tai *);
 
 #endif
