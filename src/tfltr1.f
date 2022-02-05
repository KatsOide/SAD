      subroutine tfltr1(twissm,twisss,x,
     1                  mx,my,omegas,
     1                  coss,sins,cosphi,sinphi,ns, lfno)
      use ffs
      use ffs_pointer
      use tffitcode
      use ffs_fit, only:ffs_stat
      implicit none
      type (ffs_stat) optstat
      integer*4 ndp,maxdp,ntmax
      real*8 dpstep,em0
      parameter (ndp=30,maxdp=ndp*2,dpstep=.001d0,ntmax=6100,em0=1.d-8)
      real*8 twissm(-ndp:ndp,ntwissfun),twisss(-ndp:ndp,ntwissfun)
      real*8 x(4,maxdp,ns),dp(-ndp:ndp),work(-ndp:ndp),trans(4,5)
      real*8 omegas(ns),coss(ns),sins(ns),cosphi(ns),sinphi(ns),
     $     detr,x1,x2,x3,x4,as,bs,cos1,cs,dpi,ds
      integer*4 mx(maxdp,ns),my(maxdp,ns),i,j,mdp,mdpx,
     $     ns,lfno,i1,k,n,nl,nsc
      character rad62
      logical cell0,over
      cell0=cell
      cell=.true.
      call qcell(0,optstat,.false.)
      cell=.false.
      do 1 i=11,ntwissfun
        twiss(1,0,i)=0.d0
1     continue
      mdp=ndp
      do 10 i=-ndp,ndp
        dp(i)=i*dpstep
        call qtwiss(twiss,0,1,nlat,over)
c     1.d0+dp(i) ???
        do 20 j=1,ntwissfun
          twissm(i,j)=twiss(nlat,0,j)
20      continue
        twissm(i,2)=log(twissm(i,2))
        twissm(i,5)=log(twissm(i,5))
        twissm(i,1)=asinh(twissm(i,1))
        twissm(i,4)=asinh(twissm(i,4))
        detr=twiss(nlat,0,11)*twiss(nlat,0,14)
     1      -twiss(nlat,0,12)*twiss(nlat,0,13)
        if(over .or. detr .gt. 1.d0
     1     .or. abs(twissm(i,1)) .gt. 10.d0
     1     .or. abs(twissm(i,4)) .gt. 10.d0
     1     .or. abs(twissm(i,2)) .gt. 10.d0
     1     .or. abs(twissm(i,5)) .gt. 10.d0)then
          mdp=min(abs(i)-1,mdp)
        endif
10    continue
      do 110 i=1,ntwissfun
        call spline(2*mdp+1,dp(-mdp),
     1             twissm(-mdp,i),twisss(-mdp,i),work(-mdp),0,0.d0)
110   continue
      mdpx=mdp*2
c     call spldrw(2*mdp+1,dp(-mdp),twissm(-mdp, 1),twisss(-mdp, 1))
c     call spldrw(2*mdp+1,dp(-mdp),twissm(-mdp, 4),twisss(-mdp, 4))
c     call spldrw(2*mdp+1,dp(-mdp),twissm(-mdp, 3),twisss(-mdp, 3))
c     call spldrw(2*mdp+1,dp(-mdp),twissm(-mdp, 6),twisss(-mdp, 6))
c     call spldrw(2*mdp+1,dp(-mdp),twissm(-mdp,13),twisss(-mdp,13))
c     call spldrw(2*mdp+1,dp(-mdp),twissm(-mdp,14),twisss(-mdp,14))
      x1=sqrt(em0*twiss(1,0,2))
      x2=-twiss(1,0,1)/twiss(1,0,2)*x1
      x3=sqrt(em0*twiss(1,0,5))
      x4=-twiss(1,0,4)/twiss(1,0,5)*x3
      do 40 j=1,ns
        do 30 i=1,mdpx
          x(1,i,j)=x1
          x(2,i,j)=x2
          x(3,i,j)=x3
          x(4,i,j)=x4
          mx(i,j)=ntmax
          my(i,j)=ntmax
30      continue
        coss(j)=cos(omegas(j))
        sins(j)=sin(omegas(j))
        cosphi(j)=0.d0
        sinphi(j)=1.d0
40    continue
      do 210 n=1,ntmax
        do 410 k=1,ns
          cos1     =cosphi(k)*coss(k)-sinphi(k)*sins(k)
          sinphi(k)=cosphi(k)*sins(k)+sinphi(k)*coss(k)
          cosphi(k)=cos1
          do 220 i=1,mdpx
            if(mx(i,k) .lt. ntmax .and. my(i,k) .lt. ntmax)then
              go to 220
            endif
            dpi=dpstep*.5d0*(i-1)*cos1
            if(dpi .ge. 0.d0)then
              i1=int(dpi/dpstep)
            else
              i1=int(dpi/dpstep)-1
            endif
            as=(dp(i1+1)-dpi)/dpstep
            bs=(dpi-dp(i1))/dpstep
            cs=-bs*(as+1.d0)*as*dpstep**2
            ds=-as*(bs+1.d0)*bs*dpstep**2
            do 230 j=1,ntwissfun
              twiss(nlat,0,j)=as*twissm(i1,j)+bs*twissm(i1+1,j)
     1                       +cs*twisss(i1,j)+ds*twisss(i1+1,j)
230         continue
            twiss(nlat,0,2)=exp(twiss(nlat,0,2))
            twiss(nlat,0,5)=exp(twiss(nlat,0,5))
            twiss(nlat,0,1)=sinh(twiss(nlat,0,1))
            twiss(nlat,0,4)=sinh(twiss(nlat,0,4))
            call qgettr(trans,1,nlat,0,.false.,.false.)
c     1.d0+dpi ???
c           if(i .eq. 1)then
c             write(*,'(1p5g14.6)')((trans(k,l),l=1,5),k=1,4)
c           endif
            x1=trans(1,1)*x(1,i,k)+trans(1,2)*x(2,i,k)
     1        +trans(1,3)*x(3,i,k)+trans(1,4)*x(4,i,k)
     $           +twiss(nlat,0,mfitdx)
            x2=trans(2,1)*x(1,i,k)+trans(2,2)*x(2,i,k)
     1        +trans(2,3)*x(3,i,k)+trans(2,4)*x(4,i,k)
     $           +twiss(nlat,0,mfitdpx)
            x3=trans(3,1)*x(1,i,k)+trans(3,2)*x(2,i,k)
     1        +trans(3,3)*x(3,i,k)+trans(3,4)*x(4,i,k)
     $           +twiss(nlat,0,mfitdy)
            x4=trans(4,1)*x(1,i,k)+trans(4,2)*x(2,i,k)
     1        +trans(4,3)*x(3,i,k)+trans(4,4)*x(4,i,k)
     $           +twiss(nlat,0,mfitdpy)
            x(1,i,k)=min(1.d18,max(-1.d18,x1))
            x(2,i,k)=min(1.d18,max(-1.d18,x2))
            x(3,i,k)=min(1.d18,max(-1.d18,x3))
            x(4,i,k)=min(1.d18,max(-1.d18,x4))
            if(mx(i,k) .gt. n)then
              if(abs(x1) .gt. 1.d0)then
                mx(i,k)=n
              endif
            endif
            if(my(i,k) .gt. n)then
              if(abs(x3) .gt. 1.d0)then
                my(i,k)=n
              endif
            endif
220       continue
410     continue
210   continue
      nl=0
      nsc=ns/20+1
      nl=ns/nsc
      do 1010 i=1,ns
        if(mod(i-1,nl) .eq. 0)then
          write(lfno,'(3a)')
     1'   nus   ',
     1'0.-.-.-.-.:.-.-.-.-.1.-.-.-.-.:.-.-.-.-.2.-.-.-.-.:.-.-.-.-.3',
     1' min(X,Y)'
        endif
        write(lfno,9011)omegas(i)/pi2,
     1                  (rad62(min(mx(j,i),my(j,i))/100),j=1,mdpx)
9011    format(f8.5,2x,:60a)
1010  continue
      write(lfno,'(2a)')
     1'         ',
     1'0.-.-.-.-.:.-.-.-.-.1.-.-.-.-.:.-.-.-.-.2.-.-.-.-.:.-.-.-.-.3'
      do 1011 i=1,ns
        if(mod(i-1,nl) .eq. 0)then
          write(lfno,'(3a)')
     1'   nus   ',
     1'0.-.-.-.-.:.-.-.-.-.1.-.-.-.-.:.-.-.-.-.2.-.-.-.-.:.-.-.-.-.3',
     1'    X    '
        endif
        write(lfno,9011)omegas(i)/pi2,
     1                  (rad62(    mx(j,i)         /100),j=1,mdpx)
1011  continue
      write(lfno,'(2a)')
     1'         ',
     1'0.-.-.-.-.:.-.-.-.-.1.-.-.-.-.:.-.-.-.-.2.-.-.-.-.:.-.-.-.-.3'
      do 1012 i=1,ns
        if(mod(i-1,nl) .eq. 0)then
          write(lfno,'(3a)')
     1'   nus   ',
     1'0.-.-.-.-.:.-.-.-.-.1.-.-.-.-.:.-.-.-.-.2.-.-.-.-.:.-.-.-.-.3',
     1  '    Y    '
        endif
        write(lfno,9011)omegas(i)/pi2,
     1                  (rad62(    my(j,i)         /100),j=1,mdpx)
1012  continue
      write(lfno,'(2a)')
     1'         ',
     1'0.-.-.-.-.:.-.-.-.-.1.-.-.-.-.:.-.-.-.-.2.-.-.-.-.:.-.-.-.-.3'
      cell=cell0
      return
      end
