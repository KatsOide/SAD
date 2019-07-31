--- tai_now.c.orig	Wed Oct 14 01:51:59 1998
+++ tai_now.c	Wed Jul  9 20:46:47 2003
@@ -4,5 +4,5 @@
 void tai_now(t)
 struct tai *t;
 {
-  t->x = 4611686018427387914ULL + (uint64) time((long *) 0);
+  tai_unix(t,time((long *) 0));
 }
