--- easter.c~	1998-10-14 01:51:59.000000000 +0900
+++ easter.c	2012-05-14 13:44:30.275284236 +0900
@@ -1,11 +1,12 @@
 #include <stdio.h>
+#include <stdlib.h>
 #include "caldate.h"
 
 char *dayname[7] = { "Sun", "Mon", "Tue", "Wed", "Thu", "Fri", "Sat" } ;
 
 char out[101];
 
-main(argc,argv)
+int main(argc,argv)
 int argc;
 char **argv;
 {
@@ -23,7 +24,7 @@
       caldate_frommjd(&cd,day,&weekday,&yearday);
       if (caldate_fmt((char *) 0,&cd) + 1 >= sizeof out) exit(1);
       out[caldate_fmt(out,&cd)] = 0;
-      printf("%s %s  yearday %d  mjd %d\n",dayname[weekday],out,yearday,day);
+      printf("%s %s  yearday %d  mjd %ld\n",dayname[weekday],out,yearday,day);
     }
   }
   exit(0);
