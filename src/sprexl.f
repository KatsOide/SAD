      subroutine sprexl(idxl)
      use maccbk
      implicit real*8 (a-h,o-z)
      integer idxl
c
      include 'inc/MACCODE.inc'
      include 'inc/MACFILE.inc'
c for debug
      call ptrace('sprexl '//pname(idxl)//'!',1)
c end debug
      if (ilist(2,idval(idxl)) .le. 0) then
        call errmsg('prexln',
     &       pname(idxl)//' is expanded now!',0,0)
        call expnln(idxl)
      endif
      ip=ilist(2,idval(idxl))
      llen=ilist(1,ip)
      write(outfl,'(1H ,A8,'' expanded(length='',I4,'')!=('')')
     &             pname(idxl),llen
      do 1000 i=1,llen-mod(llen,8),8
         write(outfl,'(8(2X,A8))')
     &         (pname(ilist(1,ip+i+j)),j=0,7)
 1000 continue
      write(outfl,'(8(2X,:,A8))')
     &      (pname(ilist(1,ip+i)),i=llen-mod(llen,8)+1,llen)
      write (outfl,'('' *** End of '',A8)') pname(idxl)
c.....for debug
      call ptrace('sprexl '//pname(idxl)//'!',-1)
c.....end debug
      return
      end
