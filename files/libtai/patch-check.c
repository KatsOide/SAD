--- check.c.orig	1998-10-14 01:51:59.000000000 +0900
+++ check.c	2012-05-14 13:49:05.643275172 +0900
@@ -1,4 +1,5 @@
 #include <stdio.h>
+#include <stdlib.h>
 #include <time.h>
 #include "tai.h"
 #include "leapsecs.h"
@@ -11,7 +12,7 @@
 char out[101];
 char x[TAI_PACK];
 
-main()
+int main()
 {
   struct tai t;
   struct tai t2;
@@ -36,7 +37,7 @@
       tai_sub(&t2,&t2,&t);
       packerr = tai_approx(&t2);
       for (i = 0;i < TAI_PACK;++i)
-        printf("%2.2x",(unsigned long) (unsigned char) x[i]);
+        printf("%2.2lx",(unsigned long) (unsigned char) x[i]);
       if (packerr)
         printf(" packerr=%f",packerr);
       printf(" %03d  %s",yearday,dayname[weekday]);
