      subroutine txcalc(icalc,ncalc,ip1,ip2,ic,ins,err)
      use tfstk
      use ffs
      use tffitcode
      implicit none
      integer*4 ncalc,ip1,ip2,ic,i,i1,j,k,l
      integer*4 icalc(3,*)
      logical*4 ins,err
      if(ins)then
        err=ncalc .ge. ndim1
        if(err)then
          return
        endif
        do 10 i=1,ncalc
          if(icalc(1,i) .gt. ip1)then
            i1=i
            go to 11
          elseif(icalc(1,i) .eq. ip1)then
            do 20 j=i,ncalc
              if(icalc(2,i) .gt. ip2)then
                i1=j
                go to 11
              elseif(icalc(2,i) .eq. ip2)then
                do 30 k=j,ncalc
                  if(icalc(3,i) .gt. ic)then
                    i1=k
                    go to 11
                  elseif(icalc(3,i) .eq. ic)then
                    return
                  endif
30              continue
              endif
20          continue
          endif
10      continue
        i1=ncalc+1
        go to 101
11      do 40 i=ncalc,i1,-1
          icalc(1,i+1)=icalc(1,i)
          icalc(2,i+1)=icalc(2,i)
          icalc(3,i+1)=icalc(3,i)
40      continue
101     icalc(1,i1)=min(ip1,ip2)
        icalc(2,i1)=max(ip1,ip2)
        icalc(3,i1)=ic
        ncalc=ncalc+1
        return
      else
        err=.false.
        do 210 i=1,ncalc
          if(icalc(1,i) .eq. ip1)then
            do 220 j=i,ncalc
              if(icalc(2,i) .eq. ip2)then
                do 230 k=j,ncalc
                  if(icalc(3,i) .eq. ic)then
                    do 240 l=k,ncalc-1
                      icalc(1,l)=icalc(1,l+1)
                      icalc(2,l)=icalc(2,l+1)
                      icalc(3,l)=icalc(3,l+1)
240                 continue
                    ncalc=ncalc-1
                    return
                  endif
230             continue
              endif
220         continue
          endif
210     continue
        return
      endif
      end

      subroutine tfinitcalc
      use tfstk
      use ffs
      use tffitcode
      implicit none
      logical*4 err
      flv%ncalc=0
      call txcalc(flv%icalc,flv%ncalc,nlat,nlat,mfitbx,.true.,err)
      call txcalc(flv%icalc,flv%ncalc,nlat,nlat,mfitby,.true.,err)
      call txcalc(flv%icalc,flv%ncalc,nlat,nlat,mfitax,.true.,err)
      call txcalc(flv%icalc,flv%ncalc,nlat,nlat,mfitay,.true.,err)
      call txcalc(flv%icalc,flv%ncalc,nlat,nlat,mfitnx,.true.,err)
      call txcalc(flv%icalc,flv%ncalc,nlat,nlat,mfitny,.true.,err)
      call txcalc(flv%icalc,flv%ncalc,nlat,nlat,mfitleng,.true.,err)
      return
      end
