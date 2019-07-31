      subroutine tfojit(twiss,lfno)
      use ffs
      use tffitcode
      implicit real*8(a-h,o-z)
      real*8 twiss(nlat,-ndim:ndim,ntwissfun)
      logical exist
      xj=0.d0
      yj=0.d0
      call tfgetr(xj,1.d0,'ORBITJ_ITTER x y',lfno,exist)
      if(.not. exist)then
        return
      endif
      call tfgetr(yj,1.d0,'ORBITJ_ITTER x y',lfno,exist)
      if(.not. exist)then
        return
      endif
      do i=1,nlat
        twiss(i,0,mfitdx)=twiss(i,0,mfitdx)+xj*tgauss()
        twiss(i,0,mfitdy)=twiss(i,0,mfitdy)+yj*tgauss()
      enddo
      return
      end
