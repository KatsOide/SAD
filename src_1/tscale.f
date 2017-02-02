      subroutine tscale(nlist,scale,lfno)
      use ffs
      use tffitcode
      character*(*) nlist(mfit)
      dimension scale(mfit)
      write(lfno,9001)(nlist(i),scale(i),i=1,mfit)
9001  format(:3(a,1pg15.7,1x))
      return
      end
