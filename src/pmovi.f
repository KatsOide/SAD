      subroutine pmovi(ia,ib,n)
      integer*4 ia(n),ib(n)
c     if(n .gt. 16)then
        do 10 i=1,n
          ib(i)=ia(i)
10      continue
c     else
c*VOPTION NOVEC
c       do 20 i=1,n
c         ib(i)=ia(i)
c20      continue
c     endif
      return
      end
