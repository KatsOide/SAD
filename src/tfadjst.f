      subroutine tfadjst
      use tfstk
      use ffs
      use ffs_pointer
      use tffitcode
      implicit none
      integer*8 ld,jd,ip,jp
      integer*4 i0,id,idj,i,j
      real*8 s0,sx,sy
      i0=0
      do 10 i=1,nlat-1
        if(i .gt. i0)then
          id=idtypec(i)
          if((id .ge. 2 .and. id .le. 8) .or. id == 22
     1       .or. id == 31)then
            ip=latt(i)
            ld=idvalc(i)
            s0=pos(i)
            sx=rlist(ip+kytbl(kwdx,id))-rlist(ld+kytbl(kwdx,id))
            sy=rlist(ip+kytbl(kwdy,id))-rlist(ld+kytbl(kwdy,id))
            do 20 j=i,nlat-1
              idj=idtypec(j)
              if(pos(j) .ne. s0)then
                go to 10
              endif
              if((idj .ge. 2 .and. idj .le. 8) .or. idj == 22
     1           .or. idj == 31)then
                i0=j
                s0=pos(j+1)
                jp=latt(j)
                jd=idvalc(j)
                rlist(jp+kytbl(kwdx,idj))=sx+rlist(jd+kytbl(kwdx,idj))
                rlist(jp+kytbl(kwdy,idj))=sy+rlist(jd+kytbl(kwdy,idj))
              endif
20          continue
          endif
        endif
10    continue
      return
      end
