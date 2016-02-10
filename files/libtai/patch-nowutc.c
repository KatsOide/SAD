--- nowutc.c~	1998-10-14 01:51:59.000000000 +0900
+++ nowutc.c	2012-05-14 13:36:15.334275650 +0900
@@ -1,4 +1,5 @@
 #include <stdio.h>
+#include <stdlib.h>
 #include "leapsecs.h"
 #include "tai.h"
 #include "taia.h"
@@ -10,7 +11,7 @@
 
 char x[TAIA_FMTFRAC];
 
-main()
+int main()
 {
   if (leapsecs_init() == -1) {
     fprintf(stderr,"utcnow: fatal: unable to init leapsecs\n");
@@ -23,7 +24,7 @@
   taia_tai(&now,&sec);
   caltime_utc(&ct,&sec,(int *) 0,(int *) 0);
 
-  printf("%d-%02d-%02d %02d:%02d:%02d.%s\n"
+  printf("%ld-%02d-%02d %02d:%02d:%02d.%s\n"
     ,ct.date.year
     ,ct.date.month
     ,ct.date.day
