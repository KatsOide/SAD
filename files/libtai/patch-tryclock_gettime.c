--- /dev/null	2008-03-09 17:00:01.000000000 +0100
+++ tryclock_gettime.c	2008-03-09 16:59:22.000000000 +0100
@@ -0,0 +1,8 @@
+#include <sys/time.h>
+
+main()
+{
+  struct timespec ts;
+
+  (void) clock_gettime(CLOCK_REALTIME, &ts);
+}
