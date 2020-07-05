      subroutine pgsolvcond(isp1,kx,irtc)
      use tfstk
      use iso_c_binding
      implicit none
      logical*4 normalize,cond
      type (sad_descriptor) kx
      type (sad_rlist) , pointer :: kl
      real*8, pointer:: a(:,:),b(:),c(:,:),d(:)
      integer*8 ktfmaloc,k,kb,kap,kbp,kc,kcp,kd,kdp
      integer*4 isp1,irtc,narg,nb,mb,n,m,
     $     nc,mc,nd,md,nx,itfmessage
      real*8 eps
      narg=isp-isp1
      if(narg .lt. 2 .or. narg .gt. 7 .or. narg .eq. 6)then
        irtc=itfmessage(9,'General::narg','"2, 3, 4, 5, or 7"')
        return
      endif
      k=ktastk(isp1+1)
      kb=ktastk(isp1+2)
      eps=1.d-8
      cond=.false.
      normalize=.true.
      nx=0
      if(narg .ge. 3)then
        if(ktfnonrealq(ktastk(isp1+3),eps)) go to 99
        if(narg .ge. 4)then
          if(ktfnonrealq(ktastk(isp1+4),nx)) go to 99
          if(narg .ge. 5)then
            if(ktfnonrealq(ktastk(isp1+5))) go to 99
            normalize=rtastk(isp1+5) .ne. 0
            if(narg .ge. 6)then
              cond=.true.
              if(ktfnonlistq(ktastk(isp-1))
     $             .or. ktfnonlistq(ktastk(isp))) go to 99
            endif
          endif
        endif
      endif
c      write(*,*)'pgsolvcond-1 ',nx,normalize,cond,eps
      if(cond)then
        kc=ktastk(isp1+6)
        kd=ktastk(isp1+7)
        kcp=ktfmaloc(kc,nc,mc,.false.,.false.,irtc)
        if(irtc .ne. 0)then
          return
        endif
        kdp=ktfmaloc(kd,nd,md,.true.,.true.,irtc)
        if(irtc .ne. 0)then
          call tfree(kcp)
          return
        endif
        if(md .ne. 0)then
          irtc=itfmessage(9,'General::wrongtype','"vector for #7"')
          call tfree(kcp)
          call tfree(kdp)
          return
        endif
      else
        nc=0
        mc=0
        nd=0
        md=0
        kcp=0
        kdp=0
      endif
      kap=ktfmaloc(k,n,m,.false.,.false.,irtc)
      if(irtc .ne. 0) then
        go to 98
      endif
      kbp=ktfmaloc(kb,nb,mb,.true.,.true.,irtc)
      if(irtc .ne. 0)then
        go to 98
      endif
      if(mb .ne. 0)then
        irtc=itfmessage(9,'General::wrongtype','"vector for #2"')
        call tfree(kap)
        call tfree(kbp)
        go to 98
      endif
      if(n .ne. nb)then
        irtc=itfmessage(9,'General::equalleng','"#1 and #2"')
        call tfree(kap)
        call tfree(kbp)
        go to 98
      endif
      if(cond)then
        if(mc .ne. m .or. nd. ne. nc)then
          irtc=itfmessage(9,'General::equalleng',
     $         '"matrices"')
          call tfree(kcp)
          call tfree(kdp)
          call tfree(kap)
          call tfree(kbp)
          return
        endif
      endif
      kx=kxavaloc(-1,m,kl)
c      write(*,*)'pgsolvcond-3 ',m
      call c_f_pointer(c_loc(rlist(kap)),a,[n,m])
      call c_f_pointer(c_loc(rlist(kbp)),b,[n])
      call c_f_pointer(c_loc(rlist(kcp)),c,[nc,m])
      call c_f_pointer(c_loc(rlist(kdp)),d,[nc])
      call pgsolvcond1(a,b,kl%rbody(1:m),n,m,n,
     $     c,d,nc,cond,nx,normalize,eps,
     $     .false.)
c      write(*,*)'pgsolvcond-4 ',nc,nx
      call tfree(kap)
      call tfree(kbp)
      if(cond)then
        call tfree(kcp)
        call tfree(kdp)
      endif
      return
 98   print *,'irtc=',irtc
      if(cond)then
        call tfree(kcp)
        call tfree(kdp)
      endif
      return
 99   irtc=itfmessage(9,'General::wrongtype','"Real number"')
      return
      end
c
      subroutine pgsolvcond1(a,b,x,n,m,na,c,d,nc,
     $     cond,nx,norm,eps,svd)
      use maccbk
      implicit none
      logical cond,micado,norm,svd
      integer*4,intent(in):: nx,n,m,na,nc
      real*8 a(na,m),b(n),x(m),c(nc,m),d(nc),eps
      real*8 aa(na,m),xx(m),bb(n),cc(nc,m), tt(nc),uu(nc)
      integer*4 mm(m)
      integer*4 i,j, m1

      micado=nx .ne. 0
c      write(*,*)'pgsolvcond1 ',micado,cond,nc,nx
      if(cond)then
        if(micado) then
          call tmov(a,aa,na*m)
          call tmov(b,bb,n)
          call tmov(c,cc,nc*m)
          call pmicad(a,b,x,nx,n,m,na,norm,eps,svd,.false.)
          m1=0
          do i=1,m
            if(x(i).ne.0d0)then
              m1=m1+1
              mm(m1)=i
              call tmov(aa(1,i),a(1,m1),na)
              call tmov(cc(1,i),c(1,m1),nc)
            endif
          enddo
          call tmov(bb,b,n)
          call psvdcond(a,b,x,n,m1,na,c,d,nc,eps,xx,tt,uu,cc)
          do i=m1,1,-1
            j=mm(i)
            if(j.ne.i)then
              x(j)=x(i)
              x(i)=0d0
            endif
          enddo
        else
          call psvdcond(a,b,x,n,m,na,c,d,nc,eps,xx,tt,uu,cc)
        endif
      elseif(micado)then
        call pmicad(a,b,x,nx,n,m,na,norm,eps,svd,.true.)
      else
        call tsolva(a,b,x,n,m,na,eps)
      endif
      return
      end
c
      subroutine psvdcond(a,b,x,n,m,na,c,d,nc,eps,
     $     xc,t,u,c1)
      implicit none
      real*8 epsc
      parameter (epsc=1d-13)
      integer*4 n,m,na,nc
      real*8 a(na,m),b(n),x(m),c(nc,m),d(nc),xc(m),t(nc),u(nc),
     $     c1(nc,m),eps
      integer*4 i,j,k,l
      real*8 s
      do 19 j=1,m
        do 18 i=1,nc
          c1(i,j)=c(i,j)
 18     continue
 19   continue
      do 17 l=1,nc
        do 11 i=1,l-1
          s=0d0
          do 10 k=1,m
            s=s+c(l,k)*c(i,k)
 10       continue
          t(i)=s/u(i)
 11     continue
        do 13 k=1,m
          s=0d0
          do 12 i=1,l-1
            s=s+t(i)*c(i,k)
 12       continue
          c(l,k)=c(l,k)-s
 13     continue
        s=0d0
        do 15 k=1,m
          s=s+c(l,k)**2
 15     continue
        if(s.lt.m*epsc**2) then
          do 16 k=1,m
            c(l,k)=0d0
 16       continue
          u(l)=1d0
        else
          u(l)=s
        endif
 17   continue
      call tsolvg(c1,d,xc,nc,m,nc)
      do 23 l=1,n
        s=0d0
        do 22 k=1,m
          s=s+a(l,k)*xc(k)
 22     continue
        b(l)=b(l)-s
 23   continue
      do 28 l=1,n
        do 25 i=1,nc
          s=0d0
          do 24 k=1,m
            s=s+a(l,k)*c(i,k)
 24       continue
          t(i)=s/u(i)
 25     continue
        do 27 i=1,nc
          do 26 k=1,m
            a(l,k)=a(l,k)-t(i)*c(i,k)
 26       continue
 27     continue
c     b(l)=b(l)-t(i)*d(i)
 28   continue
      call tsolva(a,b,x,n,m,na,eps)
      do 29 i=1,m
        x(i)=x(i)+xc(i)
 29   continue
      return
      end

