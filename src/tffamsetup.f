      subroutine tffamsetup(ll,em)
      use tfstk
      use tffitcode
      use ffs, only:pi2
      use ffs_fit
      use ffs_pointer, only:elatt
      use calc,only:twmov
      use eeval
      use tfcsi,only:lfnm
      implicit none
      real*8 sqrt3
c      parameter (sqrt3=sqrt(3.d0))
      parameter (sqrt3=1.732050807568877d0)
      type (sad_dlist), pointer :: kl
      type (sad_rlist), pointer :: kli
      type (sad_descriptor) kx,ki,kxnfam
      type (sad_descriptor) , save :: kxmamp,kxnfamp
      integer*4, parameter :: ivoid=9999,maxnfp=12
      integer*4 ,intent(in):: ll
      real*8 ,intent(in):: em
      real*8 twissi(50),trans(4,4),dx(4,maxnfp*2),dxp(4,maxnfp*2),dpm(-nfr:nfr+1)
      integer*8 lp
      integer*4 m,nfa,i,j,irtc,l,nfp,k
      real*8 dpw,x0,px0,y0,py0,dpi,x,y,c,s
      data kxmamp%k /0/
      lp=elatt%comp(ll)
      kfam(-nfr:nfr)=0
      jfam(-nfr:nfr)=ivoid
      jfam(0)=0
      nfam=nfr
      if(kxmamp%k == 0)then
        kxmamp=kxsymbolz('`MatchingAmplitude',18)
        kxnfamp=kxsymbolz('`NFAMP',6)
      endif
      call tclrfpe
      levele=levele+1
      kx=tfeeval(kxmamp,.true.,irtc)
      if(irtc /= 0)then
        write(lfnm,*)'Error in MatchingAmplitude, code =',irtc
        go to 9000
      endif
      if(.not. tflistq(kx,kl))then
        go to 9100
      endif
      m=kl%nl
      if(m <= 0)then
        go to 9100
      endif
      kxnfam=tfeeval(kxnfamp,.true.,irtc)
      if(ktfnonrealq(kxnfam,nfp) .or. nfp <= 0)then
        call tfdebugprint(kxnfam,
     $       'Non-positive NFAMP - FAM is skipped:',1)
        go to 9100
      endif
      nfp=min(maxnfp,nfp)
      nfa=nfr+1
      dpw=dp(nfr)-dp(-nfr)
c      dp0=(dp(nfr)+dp(-nfr))*.5d0
      call twmov(ll,twissi,1,0,.true.)
      x0 =sqrt(twissi(mfitbx)*em)
      px0=x0/twissi(mfitbx)
      y0 =sqrt(twissi(mfitby)*em)
      py0=y0/twissi(mfitby)
      call tfnormaltophysical(twissi,trans)
      dx(:,1:nfp*2)=0.d0
      do i=1,nfp
        c=Cos(pi2/nfp*(i-1))
        s=Sin(pi2/nfp*(i-1))
        dx(1,i)=c*x0
        dx(2,i)=s*px0-twissi(mfitax)/twissi(mfitbx)*dx(1,i)
        dx(3,i+nfp)=c*y0
        dx(4,i+nfp)=s*py0-twissi(mfitay)/twissi(mfitby)*dx(3,i+nfp)
      enddo
      dxp(:,1:nfp*2)=matmul(trans,dx(:,1:nfp*2))
c      do i=1,nfp*2
c        dxp(1,i)=trans(1,1)*dx(1,i)+trans(1,2)*dx(2,i)
c     $       +trans(1,3)*dx(3,i)+trans(1,4)*dx(4,i)
c        dxp(2,i)=trans(2,1)*dx(1,i)+trans(2,2)*dx(2,i)
c     $       +trans(2,3)*dx(3,i)+trans(2,4)*dx(4,i)
c        dxp(3,i)=trans(3,1)*dx(1,i)+trans(3,2)*dx(2,i)
c     $       +trans(3,3)*dx(3,i)+trans(3,4)*dx(4,i)
c        dxp(4,i)=trans(4,1)*dx(1,i)+trans(4,2)*dx(2,i)
c     $       +trans(4,3)*dx(3,i)+trans(4,4)*dx(4,i)
c      enddo
      do i=-nfr,nfr-1
        dpm(i+1)=(dp(i)+dp(i+1))*.5d0
      enddo
      dpm(-nfr)=dp(-nfr)
      dpm(nfr+1)=dp(nfr)
      do i=1,m
        ki=kl%dbody(i)
        if(.not. tfreallistq(ki,kli) .or. kli%nl /= 3)then
          go to 9000
        endif
        dpi=kli%rbody(1)
        do k=-nfr,nfr
          if(dpi >= dpm(k) .and. dpi < dpm(k+1))then
            x=kli%rbody(2)
            y=kli%rbody(3)
            if(x /= 0.d0)then
              do j=1,nfp
                jfam(nfa)=k
                dfam(1:4,nfa)=dxp(1:4,j)*x
                dp(nfa)=dp(k)
                kfam(nfa)=j
                nfa=merge(-nfa,-nfa+1,nfa >= 0)
              enddo
            endif
            if(y /= 0.d0)then
              do j=nfp+1,nfp*2
                jfam(nfa)=k
                dfam(1:4,nfa)=dxp(1:4,j)*y
                dp(nfa)=dp(k)
                kfam(nfa)=nfp-j
                nfa=merge(-nfa,-nfa+1,nfa >= 0)
              enddo
            endif
            exit
          endif
        enddo
      enddo
      if(nfa < 0)then
        nfam=-nfa
        jfam(nfa)=ivoid
        kfam(nfa)=0
      else
        nfam=nfa-1
      endif
      go to 9100
 9000 write(lfnm,*)'MatchingAmplitude must be {{dp,xamp,yamp},...}.'
 9100 l=itfdownlevel()
      return
      end

      subroutine tfnormaltophysical(twissi,trans)
      use ffs
      use tffitcode
      implicit none
      real*8 ,intent(in):: twissi(*)
      real*8 ,intent(out):: trans(4,4)
      real*8 r1,r2,r3,r4,detr,cc
      r1=twissi(mfitr1)
      r2=twissi(mfitr2)
      r3=twissi(mfitr3)
      r4=twissi(mfitr4)
      detr=r1*r4-r2*r3
      cc=sqrt(1.d0-detr)
      if(twissi(mfitdetr) < 1.d0)then
        trans(1,1)=cc
        trans(1,2)=0.d0
        trans(1,3)= r4
        trans(1,4)=-r2
        trans(2,1)=0.d0
        trans(2,2)=cc
        trans(2,3)=-r3
        trans(2,4)= r1
        trans(3,1)=-r1
        trans(3,2)=-r2
        trans(3,3)=cc
        trans(3,4)=0.d0
        trans(4,1)=-r3
        trans(4,2)=-r4
        trans(4,3)=0.d0
        trans(4,4)=cc
      else
        trans(1,1)=-r1
        trans(1,2)=-r2
        trans(1,3)= cc
        trans(1,4)= 0.d0
        trans(2,1)=-r3
        trans(2,2)=-r4
        trans(2,3)=0.d0
        trans(2,4)= cc
        trans(3,1)= cc
        trans(3,2)=0.d0
        trans(3,3)= r4
        trans(3,4)=-r2
        trans(4,1)=0.d0
        trans(4,2)= cc
        trans(4,3)=-r3
        trans(4,4)= r1
      endif
      return
      end
