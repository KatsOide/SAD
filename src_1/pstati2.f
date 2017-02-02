      subroutine pstati2(data,x,e,n,lfno)
      implicit real*8 (a-h,o-z)
      parameter (ndi=4)
      character*11 autofg,vout(2*ndi,ndi+1)
      dimension data(ndi,ndi,*),x(ndi,ndi),e(ndi,ndi),ceil(ndi,ndi),
     1          floor(ndi,ndi)
      call pclr(x,ndi*ndi)
      call pclr(e,ndi*ndi)
      do 11 i=1,ndi
        do 10 j=1,ndi
          ceil(i,j)=-1d31
          floor(i,j)=1d31
          do 20 k=1,n
            ceil(i,j)=max(ceil(i,j),data(i,j,k))
            floor(i,j)=min(floor(i,j),data(i,j,k))
 20       continue
          do 22 k=1,n
            x(i,j)=x(i,j)+data(i,j,k)
 22       continue
          x(i,j)=x(i,j)/n
          do 26 k=1,n
            e(i,j)=e(i,j)+(data(i,j,k)-x(i,j))**2
 26       continue
          e(i,j)=sqrt(e(i,j)/max(n-1,1))
 10     continue
 11   continue
      do 30 j=1,ndi
        vout(1,j)=autofg(    x(3,j),'11.8')
        vout(2,j)=autofg(    e(3,j),'11.8')
        vout(3,j)=autofg(floor(3,j),'11.8')
        vout(4,j)=autofg( ceil(3,j),'11.8')
        vout(5,j)=autofg(    x(1,j),'11.8')
        vout(6,j)=autofg(    e(1,j),'11.8')
        vout(7,j)=autofg(floor(1,j),'11.8')
        vout(8,j)=autofg( ceil(1,j),'11.8')
   30 continue
      vout(1,ndi+1)=autofg(    x(2,3),'11.8')
      vout(2,ndi+1)=autofg(    e(2,3),'11.8')
      vout(3,ndi+1)=autofg(floor(2,3),'11.8')
      vout(4,ndi+1)=autofg( ceil(2,3),'11.8')
      vout(5,ndi+1)=autofg(    x(2,1),'11.8')
      vout(6,ndi+1)=autofg(    e(2,1),'11.8')
      vout(7,ndi+1)=autofg(floor(2,1),'11.8')
      vout(8,ndi+1)=autofg( ceil(2,1),'11.8')
      write(lfno,1000) n,((vout(i,j),i=1,2*ndi),j=1,ndi+1)
 1000 format(1X,'Ndata=',I3,T19,'mean',T33,'sdev',T48,'min',T64,'max'/
     z       1X,'dx      rms  ',A,t27,A,T43,A,T59,A/
     z       1X,'       peak  ',A,t27,A,T43,A,T59,A/
     z       1X,'dy      rms  ',A,t27,A,T43,A,T59,A/
     z       1X,'       peak  ',A,t27,A,T43,A,T59,A/
     z       1X,'dEX     rms  ',A,t27,A,T43,A,T59,A/
     z       1X,'       peak  ',A,t27,A,T43,A,T59,A/
     z       1X,'dEY     rms  ',A,t27,A,T43,A,T59,A/
     z       1X,'       peak  ',A,t27,A,T43,A,T59,A/
     z       1X,'Steer   rms  ',A,t27,A,T43,A,T59,A/
     z       1X,' (k0)   max  ',A,t27,A,T43,A,T59,A)
      return
      end
