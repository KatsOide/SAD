      subroutine pqcell(latt,twiss,gammab,idp,ddp,stab)
      use ffs
      use ffs_fit, only: ffs_stat
      use tffitcode
      implicit real*8(a-h,o-z)
      type (ffs_stat) optstat
      integer*8 latt(nlat)
      dimension twiss(nlat,-ndim:ndim,ntwissfun),
     $     gammab(nlat)
      logical*4 stab
      if(cell) then
        optstat%stabx=.true.
        optstat%staby=.true.
        call qcell(idp,optstat,.false.)
        if( .not.optstat%stabx .or. .not.optstat%staby ) then
          write(*,'(A,1PD12.5,A,2(A,0PF9.4))')
     z  ' @@@@ Unstable orbit  dp=',ddp,' :',
     z  ' Trx=',optstat%tracex,' Try=',optstat%tracey
        endif
        stab=optstat%stabx.and.optstat%staby
        cellstab=stab
      else
        twiss(1,0,3)=0d0
        twiss(1,0,6)=0d0
        call qtwiss(twiss,idp,1,nlat,optstat%over)
      endif
      if(optstat%over) then
        write(*,'(A,1PD12.5)')
     z       ' @@@@ Optics over     dp=',ddp,' :'
      endif
      return
      end
