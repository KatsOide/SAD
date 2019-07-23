      subroutine tfdapert(isp1,kx,irtc)
      use tfstk
      use ffs_flag
      use tmacro
      use tparastat
      use tfcsi, only:icslfno
      implicit none
      type (sad_descriptor) :: kmode,kc,kx
      type (sad_dlist), pointer :: klc
      type (sad_rlist), pointer :: klar,klc1,klp,kml,klcx
      integer*4 isp1,irtc,narg,i,lfno,n1,itfmessage,np00,iv(3)
      real*8 range(3,3),damp,rgetgl1,phi(3),codsave(6),dampenough
      narg=isp-isp1
      if(narg .ne. 7)then
        irtc=itfmessage(9,'General::narg','"7"')
        return
      endif
      kmode=dtastk(isp1+7)
      if(.not. tfreallistq(kmode,kml))then
        irtc=itfmessage(9,'General::wrongtype','"{v1, v2} for #7"')
        return
      endif
      if(kml%nl .ne. 2)then
        irtc=itfmessage(9,'General::wrongtype','"{v1, v2} for #7"')
        return
      endif
      iv(1:2)=int(kml%rbody(1:2))
      iv(3)=6-iv(1)-iv(2)
      kc=dtastk(isp1+1)
      if(.not. tflistq(kc,klc))then
        irtc=itfmessage(9,'General::wrongtype','"List for #1"')
        return
      endif
      if(klc%nl .ne. 3)then
        irtc=itfmessage(9,'General::wrongleng','"#1","3"')
        return
      endif
      do i=2,3
        if(.not. tfreallistq(klc%dbody(iv(i)),klcx))then
          irtc=itfmessage(9,'General::wrongtype','Real List'," ")
          return
        endif
        if(klcx%nl .ne. 2)then
          irtc=itfmessage(9,'General::wrongtype','"{min, max}"," "')
          return
        endif
        range(1:2,i-1)=klcx%rbody(1:2)
      enddo
      if(.not. tfreallistq(klc%dbody(iv(1)),klc1))then
        irtc=itfmessage(9,'General::wrongtype','Real List'," ")
        return
      endif
      n1=min(200,klc1%nl)
      if(ktfnonrealq(dtastk(isp1+2),lfno))then
        return
      endif
      if(lfno .eq. -1)then
        lfno=icslfno()
      endif
      if(.not. tfreallistq(dtastk(isp1+4),klp))then
        irtc=itfmessage(9,'General::wrongtype',
     $       '"{PhaseX, PhaseY, PhaseZ} for #4"')
        return
      endif
      if(klp%nl .ne. 3)then
        irtc=itfmessage(9,'General::wrongtype',
     $       '"{PhaseX, PhaseY, PhaseZ} for #4"')
        return
      endif
      phi=klp%rbody(1:3)
      if(ktfnonrealq(dtastk(isp1+5),damp))then
        irtc=itfmessage(9,'General::wrongtype',
     $       '"Real for #5"')
        return
      endif
      if(ktfnonrealq(dtastk(isp1+6),dampenough))then
        irtc=itfmessage(9,'General::wrongtype',
     $       '"Real for #6"')
        return
      endif
      codsave=codin
      if(tfreallistq(dtastk(isp1+3),klar))then
        if(klar%nl .ne. 6 )then
          irtc=itfmessage(9,'General::wrongtype',
     $         '"{x0, px0, y0, py0, z0, dp0} for #3"')
          return
        endif
        codin=klar%rbody(1:6)
      endif
      np00=np0
      np0=max(2000,n1,n1*min(51,nint(rgetgl1('DAPWIDTH'))))
      call rsetgl1('NP',dble(np0))
      call tfdapert1(range,klc1%rbody(1),
     $     n1,kx%x(1),phi,damp,dampenough,iv(1),iv(2),lfno)
      np0=np00
      call rsetgl1('NP',dble(np0))
      irtc=0
      codin=codsave
      return
      end

      subroutine tfdapert1(range,r1,n1,trval,phi,damp,dampenough,
     $     ivar1,ivar2,lfno)
      use tfstk
      use ffs_pointer
      use ffs_flag
      use tmacro
      use tparastat
      implicit none
      integer*4 n1,itgetfpe,ivar1,ivar2,lfno
      real*8 x(np0),px(np0),y(np0),py(np0),z(np0),g(np0),dv(np0),pz(np0)
      real*8 trval,phi(3),range(3,3),r1(n1),damp,dampenough
      integer*4 kzx(2,np0)
      integer*4 kptbl(np0,6)
      integer*4 mturn(np0)
      logical*4 dapert0
c      write(*,*)'tfda1 ',phix,phiy,phiz,ivar1,ivar2
      call tfsetparam
      call tclrparaall
      dapert0=dapert
      dapert=.true.
      x(1:2)=range(1:2,1)
      y(1:2)=range(1:2,2)
      g(1:n1)=r1
      pz(1)=dble(n1)
      call trackd(latt,kptbl,x,px,y,py,z,g,dv,pz,
     1     kzx,mturn,trval,phi,damp,dampenough,ivar1,ivar2,lfno)
      dapert=dapert0
      if(itgetfpe() .ne. 0)then
        write(*,*)'DynamicApertureSurvey-FPE ',itgetfpe()
        call tclrfpe
      endif
      return
      end

      subroutine tftrad(trdtbl,trval,lfno)
      use tfstk
      use ffs
      implicit none
      real*8 trdtbl(3,6),trval,range(3,3),r1(3),phi(3)
      integer*4 lfno
      range(1:2,1)=trdtbl(1:2,1)
      range(1:2,2)=trdtbl(1:2,3)
      r1=trdtbl(1:3,6)
      phi(1:2)=0.d0
      phi(3)=-0.5d0*pi
      call tfdapert1(range,r1,3,trval,phi,0.d0,0.d0,3,1,lfno)
      return
      end
