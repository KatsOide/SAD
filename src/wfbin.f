      subroutine wfbin(f1,f2,f3,n1,n2,n3,fmax,n,iord)
      implicit none
      integer*4 n1,n2,n3,n,iord(3,*)
      real*8 f1,f2,f3,f,df,fmax
      f=n1*f1+n2*f2+n3*f3
      df=abs(f-int(f))
      if(df .gt. .5d0)then
        df=1.d0-df
      endif
      if(df .lt. fmax)then
        iord(1,n)=n1
        iord(2,n)=n2
        iord(3,n)=n3
        n=n+1
      endif
      return
      end
