      subroutine tfsetlist(k,kl,i)
      use tfstk
      implicit none
      type (sad_descriptor) ,intent(in):: k
      type (sad_descriptor) ,intent(out):: kl
      type (sad_dlist), pointer :: kll
      integer*4 ,intent(in):: i
      call descr_sad(kl,kll)
      kll%dbody(i)=k
      if(ktfnonrealq(k))then
        kll%attr=ior(kll%attr,lnonreallist)
      endif
      return
      end

      subroutine tfdimensions(isp1,kx,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) ,intent(out)::kx
      integer*4 ,intent(in):: isp1
      integer*4 ,intent(out):: irtc
      type (sad_rlist), pointer :: klr
      integer*4 isp0,m,itfmessage
      if(isp .ne. isp1+1)then
        irtc=itfmessage(9,'General::narg','"1"')
        return
      endif
      irtc=0
      if(ktfnonlistq(ktastk(isp)))then
        kx=dxnulll
      else
        isp0=isp
        call tfdimensionsstk(dtastk(isp),dtastk(isp))
        m=isp-isp0
        kx=kxavaloc(-1,m,klr)
c        do i=1,m
          klr%rbody(1:m)=dble(itastk(1,isp0+1:isp0+m))
c        enddo
      endif
      return
      end

      recursive subroutine tfdimensionsstk(k,ka0)
      use tfstk
      implicit none
      type (sad_descriptor) ,intent(in):: k,ka0
      type (sad_dlist), pointer :: list
      integer*4 isp1,isp2,i,j,m,n1
      if(ktfnonlistq(k,list))then
        return
      elseif(.not. tfsameheadq(k,ka0))then
        return
      else
        isp=isp+1
        m=list%nl
        itastk(1,isp)=m
        if(ktfnonreallistqo(list))then
          isp1=isp
          call tfdimensionsstk(list%dbody(1),ka0)
          if(isp .gt. isp1)then
            n1=isp-isp1
            isp2=isp
            do i=2,m
              call tfdimensionsstk(list%dbody(i),ka0)
              if(isp-isp2 .ne. n1)then
                isp=isp1
                return
              else
                do j=1,n1
                  if(itastk(1,isp1+j) .ne. itastk(1,isp2+j))then
                    isp=isp1
                    return
                  endif
                  isp=isp2
                enddo
              endif
            enddo
          endif
        endif
      endif
      return
      end

      integer*8 function ktfmakelist(isp0)
      use tfstk, kf=>ktfmakelist
      integer*4 ,intent(in):: isp0
      ktfmakelist=kf(isp0)
      return
      end

      subroutine tfsetlistdummy
      use mackw
      return
      end
