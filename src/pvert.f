      logical function pvert(latt,n)
      use tfstk
      use ffs
      use tffitcode
      use ffs_pointer, only:idelc,idtypec
      integer*8 latt(nlat)
c
      id=idtypec(n)
      if(id.eq.icbend) then
        lpr=5
      elseif(id.eq.icquad) then
        lpr=4
      elseif(id.eq.icsext) then
        lpr=4
      else
        pvert=.true.
        return
      endif
      phi=mod(max(rlist(latt(n)+lpr),-rlist(latt(n)+lpr)),pi)
      if(phi.gt.pi*0.25.and.phi.lt.pi*0.75) then
        pvert=.true.
      else
        pvert=.false.
      endif
      return
      end
