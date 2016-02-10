      subroutine mstack1(iop,ax,latt,twiss,istr,nstr,imon,emon,nmon)
      use tfstk
      use ffs
c----- mstack1 is also called by mstore with PUSH mode-----------------
      use tffitcode
      parameter (kstack=3)
      parameter (mstkic=16)
      dimension latt(2,nlat),twiss(nlat,-ndim:ndim,ntwissfun)
      common /mcstack/icoma,ipnt(kstack),iistck(kstack),mstkmx(kstack)
      dimension istr(*),imon(*),emon(*)
      icom=icoma
      if(iop.eq.1)then
        if(ipnt(icom).eq.0) then
          if(iistck(icom) .ne. 0) call tfree(int8(iistck(icom)))
          iistck(icom)=italoc((mstkic+1)/2)
          mstkmx(icom)=mstkic
        elseif(ipnt(icom).eq.mstkmx(icom))then
          call palocx(iistck(icom),mstkmx(icom),mstkmx(icom)+mstkic)
          mstkmx(icom)=mstkmx(icom)+mstkic
        endif
      elseif(iop.eq.2 .or. iop.eq.3 .or. iop.eq.4 .or. iop.eq.8) then
        if(ipnt(icom).eq.0)then
          print *,' !!! Stack empty.'
          return
        endif
      endif
      call mstack2(iop,ax,latt,twiss,istr,nstr,imon,emon,nmon,
     1             rlist(iistck(icom)))
      if(ipnt(icom).eq.0) then
        if(iistck(icom).ne.0) then
          call tfree(int8(iistck(icom)))
          iistck(icom)=0
        endif
      endif
      return
      end
