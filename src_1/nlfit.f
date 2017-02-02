      subroutine nlfit(x,y,wei,ndata,a,ma,lista,covar,alpha,chisq,funcs)
      implicit real*8 (a-h,o-z)
      parameter (nmax=100,eps=1d-3)
      dimension x(ndata),y(ndata),wei(ndata),
     1          a(ma),lista(ma),
     1          covar(ma,ma),alpha(ma,ma)
      mfit=ma
      do i=1,ma
        if(lista(i).eq.0) mfit=mfit-1
      enddo
      alamda=-1.
c     begin initialize for preventing compiler warning
      ochisq=0.d0
c     end   initialize for preventing compiler warning
      do it=1,30
        call mrqmin(x,y,wei,ndata,a,ma,lista,mfit,covar,alpha,ma,
     1              chisq,funcs,alamda)
        if(it.eq.1) then
          ochisq=chisq
        else
          if(chisq.lt.ochisq) then
            if((ochisq-chisq)/dble(ndata-ma).lt.eps) then
              call mrqmin(x,y,wei,ndata,a,ma,lista,mfit,covar,alpha,
     1                    ma,chisq,funcs,0d0)
              return
            endif
            ochisq=chisq
          endif
        endif
      enddo
      print *,'Iteration failed (NLFIT)'
      return
      end
