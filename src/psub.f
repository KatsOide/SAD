      subroutine psub(a,b,n)
      real*8 a(n),b(n)
c     if(n .gt. 16)then
        do 10 i=1,n
          b(i)=a(i)-b(i)
   10   continue
c     else
c*VOPTION NOVEC
c       do 20 i=1,n
c         b(i)=a(i)-b(i)
c  20   continue
c     endif
      return
      end
