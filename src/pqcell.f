      subroutine pqcell(latt,twiss,gammab,idp,ddp,stab)
      use ffs
      use tffitcode
      implicit real*8(a-h,o-z)
      logical hstab,vstab,stab,over
      integer*8 latt(nlat)
      dimension twiss(nlat,-ndim:ndim,ntwissfun),
     $     gammab(nlat)
      if(cell) then
        hstab=.true.
        vstab=.true.
        call qcell(idp,hstab,vstab,tracex,
     1             tracey,.false.,over)
        if( .not.hstab .or. .not.vstab ) then
          write(*,'(A,1PD12.5,A,2(A,0PF9.4))')
     z  ' @@@@ Unstable orbit  dp=',ddp,' :',
     z  ' Trx=',tracex,' Try=',tracey
        endif
        stab=hstab.and.vstab
        cellstab=stab
      else
        twiss(1,0,3)=0d0
        twiss(1,0,6)=0d0
        call qtwiss(twiss,idp,1,nlat,over)
      endif
      if(over) then
        write(*,'(A,1PD12.5)')
     z       ' @@@@ Optics over     dp=',ddp,' :'
      endif
      return
      end
