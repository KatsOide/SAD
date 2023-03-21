Index: src/tfprint.f
===================================================================
--- src/tfprint.f	(revision 562)
+++ src/tfprint.f	(working copy)
@@ -259,6 +259,8 @@
       enddo
       write(*,*)'Buffer is damaged at unreadbuf. ',
      $     ipoint,lrecl,' ',word(1:l)
+c      buffer(lrecl:lrecl)=char(10)
+      ipoint=lrecl+1
       call skipline
       return
       end
