      subroutine cssetp(ip)
      implicit none
      include 'inc/TFCSI.inc'
      integer*4 ip
      ipoint=ip
      return
      end

      subroutine cssets(ip)
      implicit none
      include 'inc/TFCSI.inc'
      integer*4 ip
      ios=ip
      return
      end

      subroutine cssetl(ip)
      implicit none
      include 'inc/TFCSI.inc'
      integer*4 ip
      lrecl=ip
      return
      end

      subroutine cssetlfni(ip)
      implicit none
      include 'inc/TFCSI.inc'
      integer*4 ip
      lfni=ip
      return
      end

      subroutine cssetlfno(ip)
      implicit none
      include 'inc/TFCSI.inc'
      integer*4 ip
      lfno=ip
      return
      end

      subroutine cssetlfn1(ip)
      implicit none
      include 'inc/TFCSI.inc'
      integer*4 ip
      lfn1=ip
      return
      end

      subroutine cssetrec(f)
      implicit none
      include 'inc/TFCSI.inc'
      logical*4 f
      rec=f
      return
      end

      subroutine cssetlinep(ip)
      implicit none
      include 'inc/TFCSI.inc'
      integer*4 ip
      linep=ip
      return
      end
