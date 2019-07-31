      subroutine putsti(emitx,emity,emitz,sigme,sigmz,btilt,sigmx,
     1                  sigmy,fintra,fcod)
      use tfstk
      use ffs
      use tffitcode
      parameter (ndl=16,nde=8,nd=20)
      logical fcod,fintra,start
      common /pstatis/na,nb,nc,is,isb,isc,start
      if(.not.start) return
      if(.not.fcod) return
      if(fintra) then
        if(nc.eq.0) then
          isc=italoc(nde*nd)
        elseif(mod(nc,nd).eq.0) then
          call palocx(isc,nde*nc,nde*(nc+nd))
        endif
      else
        if(nb.eq.0) then
          isb=italoc(nde*nd)
        elseif(mod(nb,nd).eq.0) then
          call palocx(isb,nde*nb,nde*(nb+nd))
        endif
      endif
      if(fintra) then
        ia1=isc+nc*nde
        nc=nc+1
      else
        ia1=isb+nb*nde
        nb=nb+1
      endif
      rlist(ia1  )=emitx
      rlist(ia1+1)=emity
      rlist(ia1+2)=emitz
      rlist(ia1+3)=sigme
      rlist(ia1+4)=sigmz
      rlist(ia1+5)=btilt
      rlist(ia1+6)=sigmx
      rlist(ia1+7)=sigmy
      return
      end
