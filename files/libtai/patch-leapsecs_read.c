--- leapsecs_read.c.orig	1998-10-14 01:51:59.000000000 +0900
+++ leapsecs_read.c	2012-05-14 13:39:03.094275680 +0900
@@ -1,3 +1,5 @@
+#include <stdlib.h>
+#include <unistd.h>
 #include <sys/types.h>
 #include <sys/stat.h>
 #include <fcntl.h>
@@ -47,4 +49,5 @@
 
   leapsecs = t;
   leapsecs_num = n;
+  return 0;
 }
