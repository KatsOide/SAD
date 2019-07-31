      subroutine mrqmin(x,y,wei,ndata,a,ma,lista,mfit,covar,alpha,nca,
     1                  chisq,funcs,alamda)
      implicit real*8 (a-h,o-z)
      parameter (nmax=100)
      dimension x(ndata),y(ndata),wei(ndata),
     1          a(ma),lista(ma),
     1          covar(nca,nca),alpha(nca,nca),
     1          atry(nmax),beta(nmax),da(nmax)
c     begin initialize for preventing compiler warning
      ochisq=0d0
c     end   initialize for preventing compiler warning
      if(alamda.lt.0d0) then
        kk=mfit+1
        do 12 j=1,ma
          ihit=0
          do 11 k=1,mfit
            if(lista(k).eq.j) ihit=ihit+1
11        continue
          if(ihit.eq.0) then
            lista(kk)=j
            kk=kk+1
          elseif(ihit.gt.1) then
            print *,'Improper permutation in lista'
            return
          endif
12      continue
        if(kk.ne.ma+1) then
          print *,'Improper permutation in lista'
          return
        endif
        alamda=1d-3
        call mrqcof(x,y,wei,ndata,a,ma,lista,mfit,alpha,beta,nca,chisq,
     1              funcs)
        ochisq=chisq
        do 13 j=1,ma
          atry(j)=a(j)
13      continue
      endif
      do 15 j=1,mfit
        do 14 k=1,mfit
          covar(j,k)=alpha(j,k)
14      continue
        covar(j,j)=alpha(j,j)*(1d0+alamda)
        da(j)=beta(j)
15    continue
      call pgaussj(covar,mfit,nca,da,1,1)
      if(alamda.eq.0d0) then
        call mrqcov(covar,nca,ma,lista,mfit)
        return
      endif
      do 16 j=1,mfit
        atry(lista(j))=a(lista(j))+da(j)
16    continue
      call mrqcof(x,y,wei,ndata,atry,ma,lista,mfit,covar,da,nca,chisq,
     1            funcs)
      if(chisq.lt.ochisq) then
        alamda=1d-1*alamda
        ochisq=chisq
        do 18 j=1,mfit
          do 17 k=1,mfit
            alpha(j,k)=covar(j,k)
17        continue
          beta(j)=da(j)
          a(lista(j))=atry(lista(j))
18      continue
      else
        alamda=1d1*alamda
        chisq=ochisq
      endif
      return
      end
