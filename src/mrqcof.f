      subroutine mrqcof(x,y,wei,ndata,a,ma,lista,mfit,alpha,beta,nalp,
     1                  chisq,funcs)
      implicit real*8 (a-h,o-z)
      external funcs
      parameter (nmax=100)
      dimension x(ndata),y(ndata),wei(ndata),a(ma),lista(mfit),
     1          alpha(nalp,nalp),beta(ma),
     1          dyda(nmax)
c
      do 12 j=1,mfit
        do 11 k=1,j
          alpha(j,k)=0d0
11      continue
        beta(j)=0d0
12    continue
      chisq=0d0
      do 15 i=1,ndata
        call funcs(x(i),a,ymod,dyda,ma)
        dy=y(i)-ymod
        do 14 j=1,mfit
          wt=dyda(lista(j))*wei(i)
          do 13 k=1,j
            alpha(j,k)=alpha(j,k)+wt*dyda(lista(k))
13        continue
          beta(j)=beta(j)+dy*wt
14      continue
        chisq=chisq+dy*dy*wei(i)
15    continue
      do 17 j=2,mfit
        do 16 k=1,j-1
          alpha(k,j)=alpha(j,k)
16      continue
17    continue
      return
      end
