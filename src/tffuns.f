      module funs
      use tfstk

      contains
      function tfsequence(isp1,isp2) result(kx)
      implicit none
      type (sad_descriptor) kx
      integer*4 ,intent(in):: isp1,isp2
      if(isp1 .ge. isp2)then
        kx=dxnull
      elseif(isp1+1 .eq. isp2)then
        kx=dtastk(isp2)
      else
        kx=kxcrelistm(isp2-isp1,ktastk(isp1+1:isp2),k_descr(ktfoper+mtfnull))
      endif
c      kx=merge(dxnull,merge(dtastk(isp2),
c     $     kxcrelistm(isp2-isp1,ktastk(isp1+1:isp2),
c     $     k_descr(ktfoper+mtfnull)),isp1+1 .eq. isp2),
c     $     isp1 .ge. isp2)
      return
      end

      function tfsameq1(isp1,iopc,irtc) result(kx)
      implicit none
      type (sad_descriptor) kx,k,k1
      integer*4 ,intent(in):: isp1,iopc
      integer*4 ,intent(out):: irtc
      integer*4 itfmessage
      kx=dxnullo
      if(isp .ne. isp1+2)then
        irtc=itfmessage(9,'General::narg','"2"')
        return
      endif
      k=dtastk(isp)
      k1=dtastk(isp1+1)
      kx%k=0
      if(k%k .eq. k1%k)then
        kx%k=ktftrue
      elseif(tfsameq(k,k1))then
        kx%k=ktftrue
      endif
      if(iopc .eq. mtfunsame)then
        kx%k=ktftrue-kx%k
      endif
      irtc=0
      return
      end

      function tfoverride(isp1,irtc) result(kx)
      implicit none
      type (sad_descriptor) kx
      type (sad_dlist), pointer :: kli
      integer*8 ki,k1
      integer*4 ,intent(in):: isp1
      integer*4 ,intent(out):: irtc
      integer*4 isp0,isp2,n,itfmessage,isp3,isp4,i,j
      if(isp1 .eq. isp)then
        kx=dxnulll
        irtc=0
        return
      elseif(isp1+1 .eq. isp)then
        if(ktastk(isp) .eq. ktfoper+mtfnull)then
          kx=dxnulll
          irtc=0
          return
        endif
      endif
      kx=dxnullo
      isp=isp+1
      isp0=isp
      do i=isp1+1,isp0-1
        ki=ktastk(i)
        if(ktflistq(ki,kli))then
          k1=kli%head%k
          if(k1 .eq. ktfoper+mtflist)then
            call tfgetllstkall(kli)
          elseif(k1 .eq. ktfoper+mtfrule .or.
     $           k1 .eq. ktfoper+mtfruledelayed)then
            isp=isp+1
            ktastk(isp)=ki
          else
            go to 9000
          endif
        else
          isp=isp+1
          ktastk(isp)=ki
        endif
      enddo
      isp2=isp
      do i=isp0+1,isp2
        isp=isp+1
        if(ktflistq(ktastk(i),kli))then
          dtastk(isp)=kli%dbody(1)
        else
          dtastk(isp)=dtastk(i)
        endif
c        dtastk(isp)=merge(kli%dbody(1),dtastk(i),
c     $       ktflistq(ktastk(i),kli))
      enddo
      n=isp-isp2
      isp3=isp
      call tfsortl(ktastk(isp2-3),.false.,n,2,ktfref,.true.,irtc)
      if(irtc .ne. 0)then
        isp=isp0-1
        return
      endif
      isp4=isp
      do i=1,n
        j=int(ktastk(isp3+i))
        if(j .ne. 0 .and.
     $       ktastk(isp0+j) .ne. ktfoper+mtfnull)then
          isp=isp+1
          ktastk(isp)=ktastk(isp0+j)
        endif
      enddo
      kx=kxmakelist(isp4)
      isp=isp0-1
      irtc=0
      return
 9000 irtc=itfmessage(9,'General::wrongtype','"List, Rule, Symbol, String, Real"')
      isp=isp0-1
      return
      end

      recursive subroutine tfgetstkstk(ks,rep)
      implicit none
      type (sad_descriptor) ,intent(in):: ks
      type (sad_descriptor) ki
      type (sad_dlist), pointer :: kl,kli
      integer*8 i,ka
      logical*4 ,intent(out):: rep
      logical*4 rep1
      if(ktfrefq(ks,ka))then
        rep=.true.
        if(ka > 3)then
          ka=ka-ispbase
          do i=ka+1,itastk2(1,ka)
            ki=dtastk(i)
            if(ktfrefq(ki))then
              call tfgetstkstk(ki,rep1)
            elseif(ktflistq(ki,kli))then
              if(kli%head%k == ktfoper+mtfnull)then
                call tfgetllstkall(kli)
              else
                isp=isp+1
                dtastk(isp)=ki
              endif
            else
              isp=isp+1
              dtastk(isp)=ki
            endif
          enddo
          return
        endif
        isp=isp+1
        dtastk(isp)=dxnull
        return
      elseif(ktflistq(ks,kl))then
        if(kl%head%k == ktfoper+mtfnull)then
          call tfgetllstkall(kl)
          rep=.true.
          return
        endif
      endif        
      rep=.false.
      isp=isp+1
      dtastk(isp)=ks
      return
      end

      subroutine tfgetllstk(list,i1,i2)
      implicit none
      type (sad_dlist) ,intent(in):: list
      type (sad_dlist), pointer :: kl
      integer*4 i,m
      integer*4 ,intent(in):: i1,i2
      if(i2 .ge. 0)then
        m=min(i2,list%nl)
      else
        m=list%nl+i2+1
      endif
      if(i1 .gt. m)then
        return
      endif
      do i=max(0,i1),m
        isp=isp+1
        dtastk(isp)=list%dbody(i)
        if(ktfsequenceq(ktastk(isp),kl))then
          isp=isp-1
          call tfgetllstkall(kl)
        endif
      enddo
      return
      end

      end module
