      subroutine pmbdata(latt,mult,iobs,nobs,observ,ndata,yhi,print,
     $     lfno)
      use ffs
      use tfstk
      use tffitcode
c----- Write table of beam size. -------      
      logical print
      character*(*) observ
      dimension latt(*),mult(*)
      dimension iobs(*)
      common /corbtune/ipyv,ippv,ipobs,ipiv,ipbckup1,ipemibuf,nemitd,
     $     nememax,iter
c
      call pmbdata1(latt,mult,rlist(ipemibuf),iobs,nobs,observ,ndata,
     $     yhi,print,lfno)
      return
      end
