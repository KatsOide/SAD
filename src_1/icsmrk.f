      integer*4 function icsmrk()
      implicit none
      include 'inc/TFCSI.inc'
      icsmrk=ipoint
      return
      end

      integer*4 function icsstat()
      implicit none
      include 'inc/TFCSI.inc'
      icsstat=ios
      return
      end

      integer*4 function icslrecl()
      implicit none
      include 'inc/TFCSI.inc'
      icslrecl=lrecl
      return
      end

      integer*4 function icslfni()
      implicit none
      include 'inc/TFCSI.inc'
      icslfni=lfni
      return
      end

      integer*4 function icslfno()
      implicit none
      include 'inc/TFCSI.inc'
      icslfno=lfno
      return
      end

      integer*4 function icslfn1()
      implicit none
      include 'inc/TFCSI.inc'
      icslfn1=lfn1
      return
      end

      integer*4 function icslinep()
      implicit none
      include 'inc/TFCSI.inc'
      icslinep=linep
      return
      end

      logical*4 function csrec()
      implicit none
      include 'inc/TFCSI.inc'
      csrec=rec
      return
      end
