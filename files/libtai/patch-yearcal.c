--- yearcal.c.orig	1998-10-14 01:51:59.000000000 +0900
+++ yearcal.c	2012-05-14 13:43:19.933277008 +0900
@@ -1,4 +1,5 @@
 #include <stdio.h>
+#include <stdlib.h>
 #include "caldate.h"
 
 char *montab[] = {
@@ -16,7 +17,7 @@
 , "December"
 } ;
 
-void main(argc,argv)
+int main(argc,argv)
 int argc;
 char **argv;
 {
