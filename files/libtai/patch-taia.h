--- taia.h.orig	Wed Oct 14 01:51:59 1998
+++ taia.h	Fri Nov 11 19:55:38 2005
@@ -9,23 +9,23 @@
   unsigned long atto; /* 0...999999999 */
 } ;
 
-extern void taia_tai();
+extern void taia_tai(struct taia *, struct tai *);
 
-extern void taia_now();
+extern void taia_now(struct taia *);
 
-extern double taia_approx();
-extern double taia_frac();
+extern double taia_approx(struct taia *);
+extern double taia_frac(struct taia *);
 
-extern void taia_add();
-extern void taia_sub();
-extern void taia_half();
-extern int taia_less();
+extern void taia_add(struct taia *, struct taia *, struct taia *);
+extern void taia_sub(struct taia *, struct taia *, struct taia *);
+extern void taia_half(struct taia *, struct taia *);
+extern int taia_less(struct taia *, struct taia *);
 
 #define TAIA_PACK 16
-extern void taia_pack();
-extern void taia_unpack();
+extern void taia_pack(char *, struct taia *);
+extern void taia_unpack(char *, struct taia *);
 
 #define TAIA_FMTFRAC 19
-extern unsigned int taia_fmtfrac();
+extern unsigned int taia_fmtfrac(char *, struct taia *);
 
 #endif
