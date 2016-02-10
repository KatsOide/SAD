      subroutine tftrad(trdtbl,trval,lfno)
      use tfstk
      use ffs
      use tffitcode
      real*8 trdtbl(3,6)
      ix=italoc(np0)
      ipx=italoc(np0)
      iy=italoc(np0)
      ipy=italoc(np0)
      iz=italoc(np0)
      ig=italoc(np0)
      idv=italoc(np0)
      ipz=italoc(np0)
      imt=italoc((np0+1)/2)
      ikptbl=italoc(np0*3)
      call tclr(rlist(ikptbl),np0*3)
      ikzx=italoc(np0)
      call tspini(1,ilist(1,ikptbl),.false.)
      rlist(ix)=trdtbl(1,1)
      rlist(ix+1)=trdtbl(2,1)
      rlist(ix+2)=trdtbl(3,1)
      rlist(iy)=trdtbl(1,3)
      rlist(iy+1)=trdtbl(2,3)
      rlist(iy+2)=trdtbl(3,3)
      rlist(ig)=trdtbl(1,6)
      rlist(ig+1)=trdtbl(2,6)
      rlist(ig+2)=trdtbl(3,6)
      rlist(ipz)=0.d0
      call trackd(ilist(1,ilattp+1),ilist(1,ikptbl),
     $     rlist(ix),rlist(ipx),rlist(iy),rlist(ipy),
     $     rlist(iz),rlist(ig),rlist(idv),rlist(ipz),
     $     ilist(1,imt),ilist(1,ikzx),trval,
     $     0.d0,0.d0,0.d0,lfno)
      call tfree(int8(ix))
      call tfree(int8(ipx))
      call tfree(int8(iy))
      call tfree(int8(ipy))
      call tfree(int8(iz))
      call tfree(int8(ig))
      call tfree(int8(idv))
      call tfree(int8(ipz))
      call tfree(int8(imt))
      call tfree(int8(ikptbl))
      call tfree(int8(ikzx))
      return
      end
