      subroutine tfsetlist(k,kl,i)
      use tfstk
      implicit none
      type (sad_descriptor) k,kl
      type (sad_list), pointer :: kll
      integer*4 i
      call descr_sad(kl,kll)
      kll%dbody(i)=k
      if(ktfnonrealq(k))then
        kll%attr=ior(kll%attr,lnonreallist)
      endif
      return
      end

      integer*8 function ktfcrelistr(m,ks,kh)
      use tfstk
      implicit none
      type (sad_descriptor) kh
      type (sad_list), pointer ::kl
      integer*4 m
      integer*8 ks(m)
      ktfcrelistr=ktaalocr(-1,m,kl)
      call tfcrelista(m,ks,kh,kl)
      return
      end

      subroutine tfcrelista(m,ks,kh,list)
      use tfstk
      implicit none
      type (sad_descriptor) kh
      type (sad_list) list
      integer*4 i,m
      integer*8 ks(m)
      logical*4 d,di
      list%dbody(0)=dtfcopy(kh)
      d=.false.
      do i=1,m
        list%body(i)=ktfcopyd(ks(i),di)
        d=d .or. di
      enddo
      if(d)then
        list%attr=ior(list%attr,lnonreallist)
      endif
      return
      end
      
      subroutine tfdimensions(isp1,kx,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) kx
      type (sad_rlist), pointer :: klr
      integer*4 isp1,irtc,isp0,i,m,itfmessage
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
        do i=1,m
          klr%rbody(i)=dble(itastk(1,isp0+i))
        enddo
      endif
      return
      end

      recursive subroutine tfdimensionsstk(k,ka0)
      use tfstk
      implicit none
      type (sad_descriptor) k,ka0
      type (sad_list), pointer :: list
      integer*4 isp1,isp2,i,j,m,n1
      logical*4 tfsameheadqk
      if(ktfnonlistqd(k,list))then
        return
      elseif(.not. tfsameheadqk(k%k,ka0%k))then
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
      integer*4 isp0
      ktfmakelist=kf(isp0)
      return
      end

      subroutine tfsetlistdummy
      use mackw
      return
      end
