      subroutine  prkick(latt,mult,istr,nstr,push,print,lfno)
      use tfstk
      use ffs
      use tffitcode

      logical push,print
      dimension latt(*),istr(*),mult(*)
      ix=italoc(nstr)
      call prkick1(latt,mult,istr,nstr,
     z             rlist(ix),push,print,lfno)
      call tfree(int8(ix))
      return
      end

