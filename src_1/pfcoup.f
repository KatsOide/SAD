      subroutine  pfcoup(latt,twiss,idp,la,lb,coup)
      use tfstk
      use ffs
      use tffitcode
      implicit real*8(a-h,o-z)
      parameter (tol=1d-16)
      logical coup
      integer*8 latt(nlat),le,ld,lp
      dimension twiss(nlat,-ndim:ndim,ntwissfun),ls(2),lf(2)
      coup=.false.
      if(la+1.gt.lb) then
        ls(1)=la+1
        lf(1)=nlat-1
        ls(2)=1
        lf(2)=lb
      else
        ls(1)=la+1
        lf(1)=lb
        ls(2)=nlat-1
        lf(2)=1
      endif
      do 1999 i=1,2
        do 110 l=ls(i),lf(i)
          l1=l-1
          ltyp=idtype(ilist(2,latt(l1)))
          le=latt(l1)
          ld=idval(ilist(2,latt(l1)))
c         if(ideal)then
            lp=ld
c         else
c           lp=le
c         endif
          if(ltyp.gt.20) goto 110
          goto (110,120,110,140,110,160,110,110,110,110,
     1          110,110,110,110,110,160,110,110,110,200),ltyp
  120     continue
          coup=rlist(lp+5).ne.0d0.and.rlist(lp+2).ne.0d0
          goto 1000
  140     coup=rlist(lp+4).ne.0d0.and.rlist(lp+2).ne.0d0
          goto 1000
  160     if(rlist(lp+2).ne.0d0) then
            dx=abs(twiss(l1,idp,mfitdx)-rlist(lp+5))
            dy=abs(twiss(l1,idp,mfitdy)-rlist(lp+6))
            coup=rlist(lp+4).ne.0d0 .or. dx.gt.tol .or. dy.gt.tol
          endif
          goto 1000
  200     coup=rlist(lp+2).ne.0
 1000     if(coup) return
  110   continue
 1999 continue
      return
      end
