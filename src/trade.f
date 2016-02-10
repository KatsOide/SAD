      subroutine trade(trans,beam,cod,bx,by,bz,br,
     $     bxx,bxy,byy,dldx,
     $     f1,f2,s,ala,al)
      use tfstk
      implicit none
      include 'inc/TMACRO1.inc'
      integer*4 i
      real*8 f1,al,dldx,bx,by,bz,
     $     bxx,bxy,byy,pr,pxi,pyi,dtdsc,ar,br,
     $     p1,brad,dtds,h1sq,alr,tlc,dp,alr1,tlc1,
     $     dp1,de,ala,s,sq,fact,fact1,f2,
     $     vx,vy,vz,qx,qy,qzi,pzi,htlc,htlc2
      real*8 trans(6,12),beam(42),cod(6)
      real*8 radi(6,6)
      pr=1.d0+cod(6)
      p1=p0*pr
      ar=.5d0*br
      pxi=(cod(2)+ar*cod(3))/pr
      pyi=(cod(4)-ar*cod(1))/pr
      sq=min(pxi**2+pyi**2,0.9999d0)
      pzi=sqrt(1.d0-sq)
      vx=pyi*bz-pzi*by
      vy=pzi*bx-pxi*bz
      vz=pxi*by-pyi*bx
      brad=(vx**2+vy**2+vz**2)
c      if(al .eq. 0.02d0)then
c      write(*,'(a,1p5g15.7)')
c     $       'trade ',pyi,cod(4),cod(1),ar,sqrt(brad)*al
c      endif
      call tradel(al,f1,f2,s,ala,alr,alr1)
      alr=alr+cod(1)*dldx
      alr1=alr1+cod(1)*dldx
      dtds=1.d0/pzi
      dtdsc=dtds*crad
      h1sq=1.d0+p1**2
      tlc=alr*dtdsc
      htlc=-h1sq*tlc
      htlc2=htlc*2.d0
      dp=htlc*brad
      u0=u0-dp
c      write(*,*)'trade ',u0,dp,brad,crad
      if(irad .gt. 6)then
        qzi=-sq/pzi
        qx=pyi*bz-qzi*by
        qy=qzi*bx-pxi*bz
        call tclr(radi,36)
        radi(6,1)=htlc2*
     $       (pzi*(-vx*bxy+vy*bxx)+vz*(pxi*bxy-pyi*bxx))
     1       -h1sq*brad*dtdsc*dldx
        radi(6,2)=(dp*dtds**2*pxi+htlc2*
     $       ((vx*by-vy*bx)*pxi/pzi-vy*bz+vz*by))/pr
        radi(6,3)=htlc2*
     $       (pzi*(vx*byy-vy*bxy)+vz*(pxi*byy-pyi*bxy))
        radi(6,4)=(dp*dtds**2*pyi+htlc2*
     $       ((-vx*by+vy*bx)*pyi/pzi+vx*bz-vz*bx))/pr
        radi(6,5)=0.d0
        radi(6,6)=-2.d0*p1*p0*tlc*brad
     1       -dp*dtds**2*sq/pr
     1       -htlc2*(qx*vx+qy*vy+vz**2)/pr
c        write(*,'(1p6g15.7)')(radi(6,i),i=1,6)
        call tradp(radi,pxi,pyi,dp,pr,ar)
        call tmuld(trans,radi)
        do i=1,6
          radi(i,i)=radi(i,i)+1.d0
        enddo
        tlc1=alr1*dtdsc*.5d0
        dp1=-h1sq*brad*tlc1
        de=-erad*dp1*h1sq*sqrt(h1sq*brad)/brhoz/pr
        beam(21)=beam(21)+de
        fact=1.d0+de/pr**2
        fact1=1.d0-de/pr**2
        beam(3)=beam(3)*fact
     $       +(ar*(2.d0*beam(5)+ar*beam(6))+pxi**2)*de
        beam(8)=beam(8)*fact
     $       +(ar*(beam(9)-beam(2)-ar*beam(4))+pxi*pyi)*de
        beam(10)=beam(10)*fact
     $       +(ar*(-2.d0*beam(7)+ar*beam(1))+pyi**2)*de
        beam(17)=beam(17)*fact1+(-ar*beam(18)/pr**2+pxi)*de
        beam(19)=beam(19)*fact1+( ar*beam(16)/pr**2+pyi)*de
        call tmulbs(beam,radi,.false.,.true.)
        beam(21)=beam(21)+de
        beam(3)=beam(3)*fact
     $       +(ar*(2.d0*beam(5)+ar*beam(6))+pxi**2)*de
        beam(8)=beam(8)*fact
     $       +(ar*(beam(9)-beam(2)-ar*beam(4))+pxi*pyi)*de
        beam(10)=beam(10)*fact
     $       +(ar*(-2.d0*beam(7)+ar*beam(1))+pyi**2)*de
        beam(17)=beam(17)*fact1+(-ar*beam(18)/pr**2+pxi)*de
        beam(19)=beam(19)*fact1+( ar*beam(16)/pr**2+pyi)*de
      endif
      if(radcod)then
        cod(6)=cod(6)+dp
        cod(2)=cod(2)+pxi*dp
        cod(4)=cod(4)+pyi*dp
        call tesetdv(cod(6))
      endif
      return
      end

      subroutine tradel(al,f1,f2,s,ala,alr,alr1)
      implicit none
      real*8 al,f1,f2,s,ala,alr,alr1,sc
      real*8 sa,sb,s1,s2,ra,rb,r1,r2,da1,d12,d2b,r0,dl
      alr=al
      alr1=al
      if(s .eq. 0.d0)then
        if(f1 .ne. 0.d0)then
          r0=min(.5d0+al/f1,1.d0)
          dl=max(al-f1*.5d0,0.d0)
          alr=r0**3/3.d0*f1+dl
          alr1=alr*r0*.75d0+dl
        endif
      elseif(s .eq. ala)then
        if(f2 .ne. 0.d0)then
          r0=min(.5d0+al/f2,1.d0)
          dl=max(al-f2*.5d0,0.d0)
          alr=r0**3/3.d0*f2+dl
          alr1=alr*r0*.75d0+dl
        endif
      else
        if(f1 .ne. 0.d0 .or. f2 .ne. 0.d0)then
          sa=s-al*.5d0
          sb=s+al*.5d0
          sc=f1*ala/(f1+f2)
          s1=max(sa,min(sc,sb,    f1*.5d0))
          s2=min(sb,max(sc,sa,ala-f2*.5d0))
          if(f1 .ne. 0.d0)then
            if(f2 .ne. 0.d0)then
              ra=min(.5d0+sa/f1,.5d0+(ala-sa)/f2,1.d0)
              rb=min(.5d0+sb/f1,.5d0+(ala-sb)/f2,1.d0)
              r1=min(.5d0+s1/f1,.5d0+(ala-s1)/f2,1.d0)
              r2=min(.5d0+s2/f1,.5d0+(ala-s2)/f2,1.d0)
            else
              ra=min(.5d0+sa/f1,1.d0)
              rb=min(.5d0+sb/f1,1.d0)
              r1=min(.5d0+s1/f1,1.d0)
              r2=min(.5d0+s2/f1,1.d0)
            endif
          else
            ra=min(.5d0+(ala-sa)/f2,1.d0)
            rb=min(.5d0+(ala-sb)/f2,1.d0)
            r1=min(.5d0+(ala-s1)/f2,1.d0)
            r2=min(.5d0+(ala-s2)/f2,1.d0)
          endif
          da1=max(s1-sa,0.d0)
          d12=max(s2-s1,0.d0)
          d2b=max(sb-s2,0.d0)
          alr=(da1*(ra*(ra+r1)+r1**2)+
     $         d12*(r1*(r1+r2)+r2**2)+
     $         d2b*(r2*(r2+rb)+rb**2))/3.d0
          alr1=(da1*(ra**2+r1**2)*(ra+r1)+
     $         d12*(r1**2+r2**2)*(r1+r2)+
     $         d2b*(r2**2+rb**2)*(r2+rb))/4.d0
        endif
      endif
      return
      end

      subroutine tradp(trans,pxi,pyi,dp,pr,ar)
      implicit none
      real*8 trans(6,6),pxi,pyi,dp,pr,ar
      trans(6,1)=trans(6,1)-ar*trans(6,4)
      trans(6,3)=trans(6,3)+ar*trans(6,2)
      trans(2,1)=pxi*trans(6,1)
      trans(2,2)=pxi*trans(6,2)
      trans(2,3)=pxi*trans(6,3)
      trans(2,4)=pxi*trans(6,4)
      trans(2,5)=pxi*trans(6,5)
      trans(2,6)=pxi*trans(6,6)
      trans(2,2)=trans(2,2)+dp/pr
      trans(2,6)=trans(2,6)-pxi*dp/pr
      trans(4,1)=pyi*trans(6,1)
      trans(4,2)=pyi*trans(6,2)
      trans(4,3)=pyi*trans(6,3)
      trans(4,4)=pyi*trans(6,4)
      trans(4,5)=pyi*trans(6,5)
      trans(4,6)=pyi*trans(6,6)
      trans(4,4)=trans(4,4)+dp/pr
      trans(4,6)=trans(4,6)-pyi*dp/pr
      return
      end

      subroutine trades(trans,beam,cod,bzs0,bzs1,f1,brhoz)
      implicit none
      real*8 trans(6,12),beam(42),cod(6),bxr,byr,bzs0,bzs1,
     $             bxx,f1,bzr,brhoz
      bzr=(bzs1-bzs0)*brhoz
      if(bzr .eq. 0.d0)then
        return
      endif
      bxx=-.5d0*bzr/f1
      bxr=bxx*cod(1)
      byr=bxx*cod(3)
      call trade(trans,beam,cod,bxr,byr,0.d0,bzs0,
     $     bxx,0.d0,bxx,0.d0,
     $     0.d0,0.d0,0.d0,0.d0,f1*.25d0)
      bxr=bxx*cod(1)
      byr=bxx*cod(3)
      call trade(trans,beam,cod,bxr,byr,0.d0,.5d0*(bzs0+bzs1),
     $     bxx,0.d0,bxx,0.d0,
     $     0.d0,0.d0,0.d0,0.d0,f1*.5d0)
      bxr=bxx*cod(1)
      byr=bxx*cod(3)
      call trade(trans,beam,cod,bxr,byr,0.d0,bzs1,
     $     bxx,0.d0,bxx,0.d0,
     $     0.d0,0.d0,0.d0,0.d0,f1*.25d0)
      return
      end
