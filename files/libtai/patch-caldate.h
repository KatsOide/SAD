--- caldate.h.orig	Wed Oct 14 01:51:59 1998
+++ caldate.h	Fri Nov 11 19:59:56 2005
@@ -7,13 +7,13 @@
   int day;
 } ;
 
-extern unsigned int caldate_fmt();
-extern unsigned int caldate_scan();
+extern unsigned int caldate_fmt(char *, struct caldate *);
+extern unsigned int caldate_scan(char *, struct caldate *);
 
-extern void caldate_frommjd();
-extern long caldate_mjd();
-extern void caldate_normalize();
+extern void caldate_frommjd(struct caldate *, long, int *, int *);
+extern long caldate_mjd(struct caldate *);
+extern void caldate_normalize(struct caldate *);
 
-extern void caldate_easter();
+extern void caldate_easter(struct caldate *);
 
 #endif
