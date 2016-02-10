--- leapsecs.c~	1998-10-14 01:51:59.000000000 +0900
+++ leapsecs.c	2012-05-14 13:36:54.101275297 +0900
@@ -1,4 +1,5 @@
 #include <stdio.h>
+#include <stdlib.h>
 #include "tai.h"
 #include "leapsecs.h"
 #include "caldate.h"
@@ -9,7 +10,7 @@
 
 char line[100];
 
-main()
+int main()
 {
   struct caldate cd;
   struct tai t;
