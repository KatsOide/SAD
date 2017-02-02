      subroutine tfgetm(ndp,xa,ya)
      use tfstk
      use ffs
      use tffitcode
      implicit none
      integer*4 i,ndp(*)
      real*8 xa,ya,alphad,getva,xa1,ya1,xxa1,xya1,yya1
      integer*4 lfnig
      data lfnig /0/
      logical exist
c
      alphad=getva(exist)
      if(.not. exist)then
        alphad=1.d0
      endif
      call nfgetm(lfnig,xa1,ya1,xxa1,xya1,yya1)
c      print *,'data gotten from pipe:', xa1,ya1
c
      do i=1,flv%nfc
        if(flv%ifitp(i) .eq. flv%measp .and.
     $       flv%ifitp1(i) .eq. flv%measp)then
          if(flv%kfit(i) .eq. 32 .and. ndp(i) .ne. 0)then
            flv%fitval(i)=(xa1-xa)*.5d0*alphad+
     $           flv%fitval(i)*(1.d0-alphad)
          endif
          if(flv%kfit(i) .eq. 34 .and.
     $         ndp(i) .ne. 0)then
            flv%fitval(i)=(ya1-ya)*.5d0*alphad+
     $           flv%fitval(i)*(1.d0-alphad)
          endif
c          print *,'tfgetm:',i,flv%kfit(i),ndp(i),flv%fitval(i)
        endif
c          print *,'tfgetm:',i,flv%ifitp(i),flv%ifitp1(i),
c     $       flv%kfit(i),ndp(i),flv%fitval(i),flv%measp
      enddo
      return
c
      end
