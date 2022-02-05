      subroutine msort(ar,n)
c        ________
c        Heapsort       ref: W.H.Press,'Numerical Recipes',p229,
c        ~~~~~~~~             Cambridge University Press 1989.
      integer ar(n),rra
      if(n.lt.2) return
      l=n/2+1
      ir=n
   10 continue
          if(l.gt.1) then
            l=l-1
            rra=ar(l)
          else
            rra=ar(ir)
            ar(ir)=ar(1)
            ir=ir-1
            if(ir.eq.1) then
              ar(1)=rra
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
              i=j
              j=j+j
            else
              j=ir+1
            endif
          goto 20
          endif
          ar(i)=rra
      goto 10
    1 continue
c     exclude equal elements
      ne=n
      i=2
   30 if(ar(i).eq.ar(i-1)) then
        do 32 j=i,ne
          ar(j-1)=ar(j)
   32   continue
        ne=ne-1
        goto 30
      else
        i=i+1
        if(i.gt.ne) then
          n=ne
          return
        else
          goto 30
        endif
      endif
      end
c
      real*8 function pselect(kk,arr,n)
      implicit none
      integer k,kk,n
      real*8 arr(n)
c     Returns the k-th smallest value in the array arr(1:n). The input
c     array willbe rearranged to have this value in location arr(k),
c     with all smaller elements moved to arr(1:k-1)(in arbutrary order)
c     and all larger elements in arr[k+1..n](also in arbitrary order).
      integer i,ir,j,l,mid
      real*8 a,temp
c
      k=min(max(1,kk),n)
      l=1
      ir=n
 1    if(ir-l.le.1)then
        if(ir-l.eq.1)then
          if(arr(ir).lt.arr(l))then
            temp=arr(l)
            arr(l)=arr(ir)
            arr(ir)=temp
          endif

        endif
        pselect=arr(k)
        return
      else
c       .not. (ir-l.le.1) => ir > l + 1
        mid=l+(ir-l)/2
        temp=arr(mid)
        arr(mid)=arr(l+1)
        arr(l+1)=temp
        if(arr(l+1).gt.arr(ir))then
          temp=arr(l+1)
          arr(l+1)=arr(ir)
          arr(ir)=temp
        endif
        if(arr(l).gt.arr(ir))then
          temp=arr(l)
          arr(l)=arr(ir)
          arr(ir)=temp
        endif
        if(arr(l+1).gt.arr(l))then
          temp=arr(l+1)
          arr(l+1)=arr(l)
          arr(l)=temp
        endif
        i=l+1
        j=ir
        a=arr(l)
 3      continue
        i=i+1
        if(arr(i).lt.a)go to 3
 4      continue
        j=j-1
        if(arr(j).gt.a)go to 4
        if(j.lt.i)go to 5
        temp=arr(i)
        arr(i)=arr(j)
        arr(j)=temp
        go to 3
 5      arr(l)=arr(j)
        arr(j)=a
        if(j.ge.k)ir=j-1
        if(j.le.k)l=i
      endif
      go to 1
      end
