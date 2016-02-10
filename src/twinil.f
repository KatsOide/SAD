      subroutine twinil(wakel,wakev,dbunch,sweight,nwx,
     1                  nb,ns)
      implicit none
      integer*4 nwx,nb,ns,i,j,i1,j1,koff,js,j1off,kds
      real*8 wakel(nwx),wakev(ns*nb),w
      real*8 dbunch(nb),sweight(ns)
      call tclr(wakev,nb*ns)
      do 10 i=1,nb
        do 20 j=1,ns
          w=sweight(j)*(1.d0+dbunch(i))
          do 30 i1=i,nb
            if(i1 .eq. i)then
              koff=-j+1
              js=j
            else
              koff=(i1-i)*(2*ns-1)-j+1
              js=1
            endif
            j1off=(i1-1)*ns
            do 40 j1=j1off+js,j1off+ns
              kds=koff+j1-j1off
              wakev(j1)=wakev(j1)+w*wakel(kds)
40          continue
30        continue
c        write(*,*)'twinil ',i,j,wakev((i-1)*ns+j)
20      continue
10    continue
      return
      end
