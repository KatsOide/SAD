      subroutine tffamsetup(lp,em)
      use tfstk
      use tffitcode
      use ffs_fit
      implicit none
      real*8 sqrt3
c      parameter (sqrt3=sqrt(3.d0))
      parameter (sqrt3=1.732050807568877d0)
      type (sad_list), pointer :: kl,kli
      integer*8 kx,ki
      integer*8 lp
      real*8 twissi(50),em,trans(4,4),dx(4,6),dxp(4,6)
      integer*4 m,nfa,i,j,irtc,itfdownlevel,l,k
      real*8 dpw,x0,px0,y0,py0,dpi,x,y,dp0
      type (sad_descriptor) kxmamp
      data kxmamp%k /0/
      integer*4 ivoid
      parameter (ivoid=9999)
      do i=-nfr,nfr
        kfam(i)=0
        jfam(i)=ivoid
      enddo
      jfam(0)=0
      nfam=nfr
      if(kxmamp%k .eq. 0)then
        kxmamp=kxsymbolz('`MatchingAmplitude',18)
      endif
      call tclrfpe
      levele=levele+1
      call tfeeval(kxmamp,kx,.true.,irtc)
      if(irtc .ne. 0)then
        write(*,*)'Error in MatchingAmplitude, code =',irtc
        go to 9100
      endif
c      call tfdebugprint(kx,'famsetup',1)
      if(.not. tflistqk(kx,kl))then
        go to 9100
      endif
      m=kl%nl
      if(m .le. 0)then
        go to 9100
      endif
      nfa=nfr+1
      dpw=dp(nfr)-dp(-nfr)
      dp0=(dp(nfr)+dp(-nfr))*.5d0
      call twmov(lp,twissi,1,0,.true.)
      x0=sqrt(twissi(mfitbx)*em)
      px0=sqrt(em/twissi(mfitbx))
      y0=sqrt(twissi(mfitby)*em)
      py0=sqrt(em/twissi(mfitby))
      call tfnormaltophysical(twissi,trans)
      dx=0.d0
      dx(1,1)=x0
      dx(1,2)=-.5d0*x0
      dx(2,2)=( px0*sqrt3-twissi(mfitax)/twissi(mfitbx)*(-x0))*.5d0
      dx(1,3)=-.5d0*x0
      dx(2,3)=(-px0*sqrt3-twissi(mfitax)/twissi(mfitbx)*(-x0))*.5d0
      dx(3,4)=y0
      dx(3,5)=-.5d0*y0
      dx(4,5)=( py0*sqrt3-twissi(mfitay)/twissi(mfitby)*(-y0))*.5d0
      dx(3,6)=-.5d0*y0
      dx(4,6)=(-py0*sqrt3-twissi(mfitay)/twissi(mfitby)*(-y0))*.5d0
      do i=1,6
        dxp(1,i)=trans(1,1)*dx(1,i)+trans(1,2)*dx(2,i)
     $       +trans(1,3)*dx(3,i)+trans(1,4)*dx(4,i)
        dxp(2,i)=trans(2,1)*dx(1,i)+trans(2,2)*dx(2,i)
     $       +trans(2,3)*dx(3,i)+trans(2,4)*dx(4,i)
        dxp(3,i)=trans(3,1)*dx(1,i)+trans(3,2)*dx(2,i)
     $       +trans(3,3)*dx(3,i)+trans(3,4)*dx(4,i)
        dxp(4,i)=trans(4,1)*dx(1,i)+trans(4,2)*dx(2,i)
     $       +trans(4,3)*dx(3,i)+trans(4,4)*dx(4,i)
      enddo
      do i=1,m
        ki=kl%body(i)
        if(.not. tfreallistq(ki,kli))then
          go to 9000
        endif
        if(kli%nl .ne. 3)then
          go to 9000
        endif
        dpi=kli%rbody(1)+dp0
        if(dpi .ge. dp(-nfr) .and. dpi .le. dp(nfr))then
          x=kli%rbody(2)
          y=kli%rbody(3)
          if(x .ne. 0.d0)then
            do j=1,3
              if(dpw .eq. 0.d0)then
                jfam(nfa)=0
              else
                jfam(nfa)=nint(2*nfr*(dpi-dp(-nfr))/dpw-nfr)
              endif
              do k=1,4
                dfam(k,nfa)=dxp(k,j)*x
              enddo
              dp(nfa)=dp(jfam(nfa))
              kfam(nfa)=j
              if(nfa .ge. 0)then
                nfa=-nfa
              else
                nfa=-nfa+1
              endif
            enddo
          endif
          if(y .ne. 0.d0)then
            do j=4,6
              if(dpw .eq. 0.d0)then
                jfam(nfa)=0
              else
                jfam(nfa)=nint(2*nfr*(dpi-dp(-nfr))/dpw-nfr)
              endif
              do k=1,4
                dfam(k,nfa)=dxp(k,j)*y
              enddo
              dp(nfa)=dp(jfam(nfa))
              kfam(nfa)=j
              if(nfa .ge. 0)then
                nfa=-nfa
              else
                nfa=-nfa+1
              endif
            enddo
          endif
        endif
      enddo
      if(nfa .lt. 0)then
        nfam=-nfa
        jfam(nfa)=ivoid
        kfam(nfa)=0
      else
        nfam=nfa-1
      endif
      go to 9100
 9000 write(*,*)'MatchingAmplitude should be {{dp,xamp,yamp},...}.'
 9100 l=itfdownlevel()
      return
      end

      subroutine tfnormaltophysical(twissi,trans)
      use ffs
      use tffitcode
      implicit none
      real*8 twissi(*),trans(4,4),r1,r2,r3,r4,detr,cc
      r1=twissi(mfitr1)
      r2=twissi(mfitr2)
      r3=twissi(mfitr3)
      r4=twissi(mfitr4)
      detr=r1*r4-r2*r3
      cc=sqrt(1.d0-detr)
      if(twissi(mfitdetr) .lt. 1.d0)then
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
