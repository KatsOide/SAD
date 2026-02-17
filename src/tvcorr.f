      subroutine tvcorr(cv,x0,al,n)
      use macmath
      use ftr
      implicit none
      integer*4 n,i
      real*8 x,w,tgauss,s,x1,al,x0
      complex*16 cv(n),cd
      if(x0 .eq. 0.d0)then
        cv(1)=(0.d0,0.d0)
        w=sqrt(al/n)
        do 1010 i=2,n
          cv(i)=cv(i-1)+dcmplx(w*tgauss(),w*tgauss())
1010    continue
        cd=cv(n)/n
        do 1020 i=2,n
          cv(i)=cv(i)-cd*(i-1)
1020    continue
        return
      endif
      x=x0/al
      s=0.d0
      x1=.25d0*(x*pi2)**2
      do 10 i=2,n/2
        w=exp(-(i-1)**2*x1)
        cv(i)=dcmplx(w*tgauss(),w*tgauss())
        cv(n-i+1)=dcmplx(w*tgauss(),w*tgauss())
        s=s+w**2
10    continue
      cv(1)=dcmplx(tgauss(),tgauss())
      s=1.d0/sqrt(4.d0*s+2.d0)
      call tcftr(cv,n,.false.)
      do 20 i=1,n
        cv(i)=cv(i)*s
20    continue
      return
      end
