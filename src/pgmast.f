      subroutine pgmast(isp1,kx,irtc)
      use tfstk
      use ffs
      use tffitcode
      use ffs_pointer
      implicit none
      integer*8 kx,kax
      integer*4 isp1,irtc,itfmessage
      integer*4 narg

      narg=isp-isp1
      if(narg .gt. 1) go to 9010
      irtc=0
      kax=ktavaloc(-1,nlat)
      call pgmast1(rlist(kax+1:kax+nlat),icomp(1:nlat),
     $     rlist(ifpos:ifpos+nlat-1))
      kx=ktflist+kax
      return
 9010 irtc=itfmessage(9,'General::narg','"1"')
      return
      end
c     
      subroutine pgmast1(rmaster,icomp,pos)
      use tfstk
      use ffs
      use tffitcode
      use ffs_pointer, only:idtypec
      implicit none
      integer*4 icomp(nlat)
      real*8 pos(nlat),rmaster(nlat)
      integer*4 i,ie,id,j
      do 10 i=1,nlat
        rmaster(i)=0.d0
 10   continue
      ie=0
      do 1010 i=1,nlat-1
        if(i .le. ie)then
          go to 1010
        endif
        id=idtypec(i)
c        if(id .eq. 20)then
c          if(rmaster(i) .eq. 0.d0)then
c            if(rlist(latt(2,i)+8) .ne. 0.d0)then
c              do 1030 j=i+1,nlat-1
c                if(idtype(latt(1,j)) .eq. 20)then
c                  if(rlist(latt(2,j)+8) .ne. 0.d0)then
c                    rmaster(i)=j
c                    rmaster(j)=-i
c                    ie=i
c                    write(*,*)'b ',i,j,rmaster(i),rmaster(j)
c                    go to 1010
c                  endif
c                endif
c 1030         continue
c              call termes(6,'?Missing end of solenoid',
c     1             pname(latt(1,i)))
c              write(*,*)'a ',i,rmaster(i)
c              ie=i
c            endif
c          endif
c          go to 1010
c        endif
c     if(id .gt. 8 .and. id .ne. 31 .and. id .eq. 22
c     1     .or. id .eq. 1)then
        if(id .gt. 8 .and. id .ne. 31 .and. id .ne. 22
     $       .and. id .ne. icSOL
     1       .or. id .eq. 1)then
          go to 1010
        endif
        ie=i
        do 1020 j=i+1,nlat-1
          if(icomp(j) .eq. icomp(i))then
            ie=j
            rmaster(j)=-i
          elseif(pos(j) .eq. pos(j+1))then
          else
            go to 1021
          endif
 1020   continue
 1021   rmaster(i)=ie
 1010 continue
      return
      end
