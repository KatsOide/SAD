      subroutine tfzap(apert,pos,
     1                 alphax,betax,dx,dpx,
     1                 alphay,betay,lfno)
      use ffs
      use tffitcode
      dimension alphax(nlat),betax(nlat),dx(nlat),dpx(nlat),
     1          alphay(nlat),betay(nlat),pos(nlat)
      write(lfno,9002)nlat,pos(nlat)
9002  format(i6,',',f12.3)
      do 10 i=1,nlat
        write(lfno,9001)pos(i),betax(i),alphax(i),betay(i),alphay(i),
     1                    dx(i),dpx(i),apert
9001    format(f11.3,6f9.3,f7.3)
10    continue
      return
      end
