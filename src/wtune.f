      real*8 function wtune(zbuf,zwork,lsp)
      use ftr
      implicit none
      integer*4 lsp,i,im,im1,im2
      real*8 a,am,a1,a2
      complex*16 zbuf(lsp),zwork(lsp)
      call tmov(zbuf,zwork,lsp*2)
      call tcftr(zwork,lsp,.false.)
      im=0
      am=0.d0
      do 10 i=1,lsp
        a=dble(zwork(i))**2+imag(zwork(i))**2
        if(a .gt. am)then
          im=i
          am=a
        endif
10    continue
      im1=mod(im+lsp-2,lsp)+1
      im2=mod(im,lsp)+1
      a1=dble(zwork(im1))**2+imag(zwork(im1))**2
      a2=dble(zwork(im2))**2+imag(zwork(im2))**2
      wtune=(im-1+(a2-a1)/2.d0/(am*2.d0-a1-a2))/lsp
      if(wtune .gt. .5d0)then
        wtune=wtune-1.d0
      endif
      return
      end
