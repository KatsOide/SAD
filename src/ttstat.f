      subroutine ttstat(np,x,px,y,py,z,g,dv,wp,
     1                  title,sa,ss,es,
     1                  lum,mat,lfno)
      use tfstk
      use ffs
      use tffitcode
      implicit none
      integer*4 maxbin,mbin
      parameter (maxbin=16,mbin=maxbin*2+1)
      integer*4 lenw,np,i,j,nbin,lfno
      real*8 x(np),px(np),y(np),py(np),z(np),g(np),dv(np),wp(np)
      real*8 s(8),sa(6),ss(6,6)
      real*8 param(21),param1(21),ss1(6,6),s1m,s2m,s3m,s4m,s5m,s6m,
     $     w,wn,s1,s2,s3,s4,s5,s6,s1w,s2w,s3w,s4w,s5w,s6w,
     $     sigz,ecorr,temix,temiy,temix1,temiy1,sigp,y4,y6,y8,
     $     es,dy,tlum,alum
      real*8 tgetgcut
      character*(*) title
      character*10 tag(6),autofg
      logical*4 lum,mat,waks
      data tag /'        x ','    px/p0 ',
     1          '        y ','    py/p0 ',
     1          '        z ','    dp/p0 '/
      s1m=0.d0
      s2m=0.d0
      s3m=0.d0
      s4m=0.d0
      s5m=0.d0
      s6m=0.d0
      call tclr(sa,6)
      call tclr(ss,36)
c      waks=twake .or. lwake
      waks=.false.
      wn=0.d0
      w=1.d0
      do 110 i=1,np
        if(waks)then
          w=wp(i)
        endif
        wn=wn+w
c        s6=g(i)*(2.d0+g(i))
        s6=g(i)
        s1=x(i)
        s2=px(i)*(1.d0+s6)
        s3=y(i)
        s4=py(i)*(1.d0+s6)
        s5=z(i)
        s1w=s1*w
        s2w=s2*w
        s3w=s3*w
        s4w=s4*w
        s5w=s5*w
        s6w=s6*w
        s1m=max(abs(s1),s1m)
        s2m=max(abs(s2),s2m)
        s3m=max(abs(s3),s3m)
        s4m=max(abs(s4),s4m)
        s5m=max(abs(s5),s5m)
        s6m=max(abs(s6),s6m)
        sa(1)=sa(1)+s1w
        ss(1,1)=ss(1,1)+s1*s1w
        ss(2,1)=ss(2,1)+s2*s1w
        ss(3,1)=ss(3,1)+s3*s1w
        ss(4,1)=ss(4,1)+s4*s1w
        ss(5,1)=ss(5,1)+s5*s1w
        ss(6,1)=ss(6,1)+s6*s1w
        sa(2)=sa(2)+s2w
        ss(2,2)=ss(2,2)+s2*s2w
        ss(3,2)=ss(3,2)+s3*s2w
        ss(4,2)=ss(4,2)+s4*s2w
        ss(5,2)=ss(5,2)+s5*s2w
        ss(6,2)=ss(6,2)+s6*s2w
        sa(3)=sa(3)+s3w
        ss(3,3)=ss(3,3)+s3*s3w
        ss(4,3)=ss(4,3)+s4*s3w
        ss(5,3)=ss(5,3)+s5*s3w
        ss(6,3)=ss(6,3)+s6*s3w
        sa(4)=sa(4)+s4w
        ss(4,4)=ss(4,4)+s4*s4w
        ss(5,4)=ss(5,4)+s5*s4w
        ss(6,4)=ss(6,4)+s6*s4w
        sa(5)=sa(5)+s5w
        ss(5,5)=ss(5,5)+s5*s5w
        ss(6,5)=ss(6,5)+s6*s5w
        sa(6)=sa(6)+s6w
        ss(6,6)=ss(6,6)+s6*s6w
110   continue
      do 140 i=1,6
        sa(i)=sa(i)/wn
140   continue
      do 150 i=1,6
        do 160 j=1,i
          ss(i,j)=ss(i,j)/wn-sa(i)*sa(j)
          ss(j,i)=ss(i,j)
160     continue
150   continue
      if(title .ne. ' ')then
        if(mat)then
          write(lfno,9001)
     1         ' Statistics at ',title(1:lenw(title)),':',
     1         '   particles =',np,
     $         ' RAD: ',rad,', RFSW: ',rfsw,', GAUSS: ',gauss,
     $         ', DP =',autofg(dpmax,'10.6'),
     $         ', DP0 =',autofg(dp0,'10.6'),
     1         ', GCUT =',autofg(tgetgcut(),'10.6')
9001      format(a,a,a,a,i5,/,3(a,l1),6a)
          write(lfno,9002)tag
9002      format(12x,6a)
        endif
        write(lfno,9003)'   C of M ',(sngl(sa(j)),j=1,6)
        if(mat)then
          do 210 i=1,6
            write(lfno,9003)tag(i),(sngl(ss(i,j)),j=1,i)
9003        format(a,': ',1p,6E10.3)
210       continue
        endif
        sigz=sqrt(ss(5,5))
        if(ss(5,5) .gt. 0.d0)then
          ecorr=ss(6,5)/sigz
        else
          ecorr=0.d0
        endif
        call tmov(ss,ss1,36)
        call tbdecoup(ss1,param,param1,temix,temiy,temix1,temiy1,sigp)
c Kikuchi added one line following. 9APR'93
ccc        call putsti(temix,temiy,0d0,0d0,0d0,0d0,0d0,0d0,.false.,.true.)
        write(lfno,9006)temix,param(2),param(1),param(7),param(8),
     1                  temiy,param(5),param(4),param(9),param(10),
     1                  temix1,param1(2),param1(1),param1(7),param1(8),
     1                  temiy1,param1(5),param1(4),param1(9),param1(10),
     1                  param1(11),param1(12),param1(13),param1(14),
     1                  param1(11)*param1(14)-param1(12)*param1(13),
     $       sqrt(ss(1,1)),sqrt(ss(3,3)),
     $       atan2(-2.d0*ss(3,1),ss(1,1)-ss(3,3))*.5d0,
     $       sqrt(ss(2,2)),sqrt(ss(4,4)),
     1                  sqrt(sigp),sigz,ecorr
9006    format(' x-y projected(coupled) parameters:'/
     1         ' emitx:',1pg11.4,' bx:',g11.4,
     1            ' ax:',  g11.4,' ex:',g11.4,' epx:',g11.4/
     1         ' emity:',  g11.4,' by:',g11.4,
     1            ' ay:',  g11.4,' ey:',g11.4,' epy:',g11.4/
     1         ' x-y decoupled parameters:'/
     1         ' emitu:',1pg11.4,' bu:',g11.4,
     1            ' au:',  g11.4,' eu:',g11.4,' epu:',g11.4/
     1         ' emitv:',  g11.4,' bv:',g11.4,
     1            ' av:',  g11.4,' ev:',g11.4,' epv:',g11.4/
     1         '    r1:',  g11.4,' r2:',g11.4,
     1            ' r3:',  g11.4,' r4:',g11.4,'detr:',g11.4/
     $         '   sigx:', g11.4,'  sigy:',g11.4,'   tilt:',g11.4/
     $         '  sigpx:', g11.4,' sigpy:',g11.4/
     1         ' sigp/p:', g11.4,'  sigz:',g11.4,' dp/p/z:',g11.4,
     $       '/sigz')
        if(twake .or. lwake)then
          ss(1,1)=temix**2
          ss(3,3)=temiy**2
        endif
        if(.not. lum)then
          return
        endif
        y4=0.d0
        y6=0.d0
        y8=0.d0
        do 310 i=1,np
          s(1)=x(i)
          s(2)=px(i)
          s(3)=y(i)
          s(4)=py(i)
          s(5)=z(i)
          s(6)=g(i)
          s(7)=dv(i)
          call tconv(s,s,1)
          dy=(s(3)-sa(3))**2/ss(3,3)
          y4=y4+dy**2
          y6=y6+dy**3
          y8=y8+dy**4
310     continue
        alum=tlum(np/2,x,y,sa(1),sa(3),sqrt(ss(1,1)),sqrt(ss(3,3)),nbin)
        write(lfno,9004)alum,nbin,nbin
9004    format(' ',/,'   Effective beam area:',1p,G10.4,' mc2',
     1         '     (',i3,'*',i3,' bin)')
        write(lfno,9005)y4/np,y6/np,y8/np
9005    format('   4-, 6-, 8th Vertical moments:',1p,3G12.4)
        es=alum
      endif
      return
      end
