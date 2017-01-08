      subroutine tspect(isp1,n,cx)
      use tfstk
      use tmacro
      implicit real*8 (a-h,o-z)
      dimension cx(2)
      lpa=itastk(2,lspect+isp1)
      if(n .lt. ilist(1,lpa+2))then
        return
      endif
      lp=ilist(2,lpa)
      ipd=ilist(1,lp+1)+2
      ilist(1,lp+1)=ipd
      ioffa=ilist(1,lp+5)
      rlist(ioffa+ipd-2)=cx(1)
      rlist(ioffa+ipd-1)=cx(2)
      nd=ilist(1,lp)
      if(ipd .lt. nd*2)then
        return
      endif
      ilist(1,lp+1)=0
      call tcftr(rlist(ioffa),nd,.false.)
      ilist(1,lp+2)=0
      joff=itastk(1,lspect+isp1)+1
      do 50 j=0,nd-1
        rlist(joff+j)=rlist(joff+j)+
     1                (rlist(ioffa+j*2)**2+rlist(ioffa+j*2+1)**2)*
     1                rlist(lp+6)
c       write(*,*)isp1,j,joff+j,rlist(joff+j)
50    continue
      return
      end

      subroutine tsptrm
      use tfstk
      use tmacro
      implicit none
      integer*4 lp,lpa,isp1
      do 10 isp1=1,nspect
        lpa=itastk(2,lspect+isp1)
        lp=ilist(2,lpa)
        call tfree(int8(ilist(1,lp+5)))
        call tfree(int8(lp))
10    continue
      return
      end
