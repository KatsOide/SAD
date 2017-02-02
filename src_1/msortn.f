C   09/10/90 204301711  MEMBER NAME  MSORTN   *.FORT     M  E2FORT
      subroutine msortn(ar,br,cr,ndim,m,n)
c        ________
c        Heapsort       ref: W.H.Press,'Numerical Recipes',p229,
c        ~~~~~~~~             Cambridge University Press 1989.
      parameter (mmax=10)
      integer ar(n),rra
      real*8 br(n),cr(ndim,m),rrb,rrc(mmax)
      if(m.gt.mmax) then
        print *,' Too large dimension in MSORTN'
        return
      endif
      if(n.lt.2) return
      l=n/2+1
      ir=n
   10 continue
          if(l.gt.1) then
            l=l-1
            rra=ar(l)
            rrb=br(l)
            do 11 k=1,m
              rrc(k)=cr(l,k)
   11       continue
          else
            rra=ar(ir)
            rrb=br(ir)
            ar(ir)=ar(1)
            br(ir)=br(1)
            do 12 k=1,m
              cr(ir,k)=cr(1,k)
   12       continue
            ir=ir-1
            if(ir.eq.1) then
              ar(1)=rra
              br(1)=rrb
              do 13 k=1,m
                cr(1,k)=rrc(k)
   13         continue
              goto 1
            endif
          endif
          i=l
          j=l+l
   20     if(j.le.ir) then
            if(j.lt.ir) then
              if(ar(j).lt.ar(j+1)) j=j+1
            endif
            if(rra.lt.ar(j)) then
              ar(i)=ar(j)
              br(i)=br(j)
              do 21 k=1,m
                cr(i,k)=cr(j,k)
   21         continue
              i=j
              j=j+j
            else
              j=ir+1
            endif
          goto 20
          endif
          ar(i)=rra
          br(i)=rrb
          do 22 k=1,m
            cr(i,k)=rrc(k)
   22     continue
      goto 10
    1 continue
cc    exclude equal elements
c     ne=n
c     i=2
c  30 if(ar(i).eq.ar(i-1)) then
c       do 32 j=i,ne
c         ar(j-1)=ar(j)
c         br(j-1)=br(j)
c         do 31 k=1,m
c           cr(j-1,k)=cr(j,k)
c  31     continue
c  32   continue
c       ne=ne-1
c       goto 30
c     else
c       i=i+1
c       if(i.gt.ne) then
c         n=ne
c         return
c       else
c         goto 30
c       endif
c     endif
      end
