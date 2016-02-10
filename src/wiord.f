      subroutine wiord(maxord,iord,aord,rfsw)
      implicit none
      integer*4 iord(3,*),maxord,i,j,m,m1,l,jm,lcoeff
      real*8 aord(3,*)
      logical*4 rfsw
      iord(1,1)=0
      iord(2,1)=0
      iord(3,1)=0
      iord(1,2)=0
      iord(2,2)=0
      iord(3,2)=0
      if(rfsw)then
        iord(1,3)=0
        iord(2,3)=0
        iord(3,3)=0
        l=4
      else
        l=3
      endif
      do 10 m=1,maxord,2
        do 20 j=m,1-m,-1
          iord(1,l)=j
          iord(2,l)=m-abs(j)
          iord(3,l)=0
          l=l+1
20      continue
        if(rfsw)then
          do 40 i=1,m
            jm=m-i
            do 50 j=-jm,jm
              iord(1,l)=j
              iord(2,l)=jm-abs(j)
              iord(3,l)=i
              l=l+1
50          continue
            do 60 j=jm-1,1-jm,-1
              iord(1,l)=j
              iord(2,l)=abs(j)-jm
              iord(3,l)=i
              l=l+1
60          continue
40        continue
        endif
        m1=m+1
        if(m1 .gt. maxord)then
          go to 10
        endif
        if(rfsw)then
          do 70 i=m1,1,-1
            jm=m1-i
            do 80 j=1-jm,jm-1
              iord(1,l)=j
              iord(2,l)=abs(j)-jm
              iord(3,l)=i
              l=l+1
80          continue
            do 90 j=jm,-jm,-1
              iord(1,l)=j
              iord(2,l)=jm-abs(j)
              iord(3,l)=i
              l=l+1
90          continue
70        continue
        endif
        do 100 j=1-m1,m1
          iord(1,l)=j
          iord(2,l)=m1-abs(j)
          iord(3,l)=0
          l=l+1
100     continue
10    continue
      lcoeff=l-1
      do 110 i=1,3
        do 120 j=1,lcoeff
          aord(i,j)=iord(i,j)
120     continue
110   continue
c     do 310 i=1,lcoeff
c      write(*,9001)i,(iord(j,i),j=1,3)
c9001    format(4i5)
c310   continue
      return
      end
