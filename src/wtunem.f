      subroutine wtunem(zbuf,r,lsp,f1,fa,nfre)
      use ftr
      implicit none
      integer*4 lsp,i,j,k,l,nfre,i1,i2
      real*8 a,ri,f1(3,nfre),fa(3,nfre),r(lsp),xi,b,fi,r1,r2,r3,rl
      complex*16 zbuf(lsp,3)
      do 110 k=1,3
        call tcftr(zbuf(:,k),lsp,.false.)
        r1=dble(zbuf(lsp,k))**2+imag(zbuf(lsp,k))**2
        r2=dble(zbuf(1  ,k))**2+imag(zbuf(1  ,k))**2
        rl=(dble(zbuf(lsp-1,k))**2+imag(zbuf(lsp-1,k))**2
     1          +r2)*.5d0+r1
        do 102 i=2,lsp
          r3=dble(zbuf(i,k))**2+imag(zbuf(i,k))**2
          r(i-1)=(r1+r3)*.5d0+r2
          r1=r2
          r2=r3
102     continue
        r(lsp)=rl
        do 120 j=1,nfre
          f1(k,j)=0.d0
          fa(k,j)=0.d0
120     continue
        do 130 i=1,lsp
          i1=mod(i+lsp-2,lsp)+1
          i2=mod(i,lsp)+1
          a=r(i1)+r(i2)-2.d0*r(i)
          if(a .lt. 0.d0)then
            xi=-(r(i2)-r(i1))*.5d0/a
            if(xi .ge. -.5d0 .and. xi .lt. .5d0)then
              ri=r(i)+xi*(r(i2)-r(i1))*.25d0
              fi=i-1+xi
              do 160 j=1,nfre
                if(abs(f1(k,j)-fi) .lt. 1.d0)then
                  if(ri .gt. fa(k,j))then
                    fa(k,j)=ri
                    f1(k,j)=fi
                  endif
                  go to 130
                endif
160           continue
              do 140 j=1,nfre
                if(ri .gt. fa(k,j))then
                  do 150 l=nfre,j+1,-1
                    fa(k,l)=fa(k,l-1)
                    f1(k,l)=f1(k,l-1)
150               continue
                  fa(k,j)=ri
                  f1(k,j)=(i-1)+xi
                  go to 130
                endif
140           continue
            endif
          endif
130     continue
        a=1.d0/lsp**2
        b=1.d0/lsp
        do 210 i=1,nfre
          fa(k,i)=fa(k,i)*a
          f1(k,i)=f1(k,i)*b
          if(f1(k,i) .gt. .5d0)then
            f1(k,i)=f1(k,i)-1.d0
          endif
210     continue
110   continue
      return
      end
