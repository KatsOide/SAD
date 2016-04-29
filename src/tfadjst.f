      subroutine tfadjst(latt,pos)
      use tfstk
      use ffs
      use tffitcode
      implicit none
      integer*4 latt(2,nlat)
      real*8 pos(nlat)
      integer*4 i0,i,id,ip,ld,idj,jp,jd,j
      real*8 s0,sx,sy
      i0=0
      do 10 i=1,nlat-1
        if(i .gt. i0)then
          id=idtype(latt(1,i))
          if((id .ge. 2 .and. id .le. 8) .or. id .eq. 22
     1       .or. id .eq. 31)then
            ip=latt(2,i)
            ld=idval(latt(1,i))
            s0=pos(i)
            sx=rlist(ip+kytbl(kwdx,id))-rlist(ld+kytbl(kwdx,id))
            sy=rlist(ip+kytbl(kwdy,id))-rlist(ld+kytbl(kwdy,id))
            do 20 j=i,nlat-1
              idj=idtype(latt(1,j))
              if(pos(j) .ne. s0)then
                go to 10
              endif
              if((idj .ge. 2 .and. idj .le. 8) .or. idj .eq. 22
     1           .or. idj .eq. 31)then
                i0=j
                s0=pos(j+1)
                jp=latt(2,j)
                jd=idval(latt(1,j))
                rlist(jp+kytbl(kwdx,idj))=sx+rlist(jd+kytbl(kwdx,idj))
                rlist(jp+kytbl(kwdy,idj))=sy+rlist(jd+kytbl(kwdy,idj))
              endif
20          continue
          endif
        endif
10    continue
      return
      end
