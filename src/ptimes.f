      subroutine ptimes(a,b,n)
      real*8 a,b(n)
c     if(n .gt. 16)then
        do 10 i=1,n
          b(i)=b(i)*a
   10   continue
c     else
c*VOPTION NOVEC
c       do 20 i=1,n
c         b(i)=b(i)*a
c  20   continue
c     endif
      return
      end
