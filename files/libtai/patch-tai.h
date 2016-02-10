--- tai.h.orig	Wed Oct 14 01:51:59 1998
+++ tai.h	Fri Nov 11 19:57:17 2005
@@ -7,16 +7,18 @@
   uint64 x;
 } ;
 
-extern void tai_now();
+#define tai_unix(t,u) ((void) ((t)->x = 4611686018427387914ULL + (uint64) (u)))
+
+extern void tai_now(struct tai *);
 
 #define tai_approx(t) ((double) ((t)->x))
 
-extern void tai_add();
-extern void tai_sub();
+extern void tai_add(struct tai *, struct tai *, struct tai *);
+extern void tai_sub(struct tai *, struct tai *, struct tai *);
 #define tai_less(t,u) ((t)->x < (u)->x)
 
 #define TAI_PACK 8
-extern void tai_pack();
-extern void tai_unpack();
+extern void tai_pack(char *, struct tai *);
+extern void tai_unpack(char *, struct tai *);
 
 #endif
