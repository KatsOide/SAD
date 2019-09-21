C created on HP by N. Yamamoto, Apr. 26,'93
C   08/12/92 212101617  MEMBER NAME  NDELW    *.FORT     M  E2FORT
      subroutine ndelw(wl,wa,wp,wd,
     &                 latt,icomp,couple,geo,pos,master)
      use tfstk
      use ffs
      use tffitcode
      use sad_main
      use ffs_pointer, only:elatt,idelc,idtypec
      real*8 wl,wa,wp,wd
      dimension geo(3,4,nlat),pos(nlat)
      integer*8 latt(nlat),le,idv
      dimension icomp(nlat),couple(nlat),master(nlat)
c
      real*8 kx,ky
      real*8 delz
      integer int
c
      delz(int)=wa*sin(kx*geo(1,4,int)+ky*geo(2,4,int)-wp)
c
      print *,' ndelw', wl,wp,wa,wd
      kx=pi2/wl*cos(wd)
      ky=pi2/wl*sin(wd)
      do 10 i=1,nlat-1
        id=idtypec(i)
        if(id .eq. 1 .or. (id .gt. 8 .and. id .ne. 20))then
          go to 10
        endif
        ie=master(i)
        if(ie .le. 0)then
          go to 10
        endif
        i0=i
        s=(pos(i0)+pos(ie+1))*.5d0
        delx=delz(ie)*geo(3,1,ie)
        dely=delz(ie)*geo(3,2,ie)
        if(id .eq. 20)then
          istep=ie-i0
        else
          istep=1
        endif
        do 110 j=i0,ie,istep
          if(master(j) .ne. 0)then
            le=latt(j)+5
            dx1=delx
            dy1=dely
            if(id .eq. 2)then
              le=le+4
            elseif(id .eq. 20)then
              le=le-2
              idv=idval(idcomp(elatt,j))
              dx1=rlist(idv+3)-delx
              dy1=rlist(idv+4)-dely
            endif
            rlist(le)=dx1
            rlist(le+1)=dy1
          endif
110     continue
10    continue
      return
      end
