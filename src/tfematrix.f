      subroutine tfematrix(kl1,kl2,kx,iopc1,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) kx,kx1,tfearray
      type (sad_dlist) kl1,kl2
      type (sad_dlist), pointer :: kind1,kind2,kb1,klx,klx1
      integer*4 iopc1,irtc
      call descr_sad(kl1%dbody(1),kind1)
      call descr_sad(kl2%dbody(1),kind2)
      if(tfsamelistqo(kind1,kind2))then
        kx1=tfearray(kl1%dbody(2),kl2%dbody(2),iopc1,irtc)
        if(irtc .ne. 0)then
          return
        endif
        call descr_sad(kl1%dbody(2),kb1)
        if(ktflistq(kx1,klx1) .and. klx1%nl .eq. kb1%nl)then
          kx%k=ktflist+ktadaloc(-1,2,klx)
          klx%head=dtfcopy1(kl1%head)
          klx%dbody(1)=dtfcopy1(kl1%dbody(1))
          klx%dbody(2)=dtfcopy1(kx1)
        else
          kx=kx1
        endif
        return
      endif
      irtc=-1
      return
      end







