      Subroutine sprlin(idxl)
      use maccbk
      use mackw
      use macfile
      implicit none
      integer*8 ip
      integer*4 idxl,llen,i,j
c
c for debug
c     call ptrace('prexl '//pname(idxl)//'!',1)
c end debug
      ip=idval(idxl)
      llen=ilist(1,ip)
c     write(outfl,'(1H ,A8,'' '',I4,''!=('')')
c    &             pname(idxl),llen
      write(outfl,'(1H ,''LINE '',A8,''=('')')
     &             pname(idxl)
      do 1000 i=1,llen-mod(llen,5),5
         write(outfl,'(5(1X,I5,''*'',A8))')
     &         (ilist(1,ip+i+j),pname(ilist(2,ip+i+j)),j=0,4)
 1000 continue
      write(outfl,'(5(1X,:,I5,''*'',A8))')
     &  (ilist(1,ip+i),pname(ilist(2,ip+i)),
     &                 i=llen-mod(llen,5)+1,llen)
c     write (outfl,'('' *** End of '',A8,)') pname(idxl)
      write (outfl, *)');'
c for debug
c     call ptrace('prexl '//pname(idxl)//'!',-1)
c end debug
      return
      end
