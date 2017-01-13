--- caltime.h.orig	Wed Oct 14 01:51:59 1998
+++ caltime.h	Fri Nov 11 20:01:16 2005
@@ -2,6 +2,7 @@
 #define CALTIME_H
 
 #include "caldate.h"
+#include "tai.h"
 
 struct caltime {
   struct caldate date;
@@ -11,10 +12,10 @@
   long offset;
 } ;
 
-extern void caltime_tai();
-extern void caltime_utc();
+extern void caltime_tai(struct caltime *, struct tai *);
+extern void caltime_utc(struct caltime *, struct tai *, int*, int *);
 
-extern unsigned int caltime_fmt();
-extern unsigned int caltime_scan();
+extern unsigned int caltime_fmt(char *, struct caltime *);
+extern unsigned int caltime_scan(char *, struct caltime *);
 
 #endif
