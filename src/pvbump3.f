      subroutine pvbump3(kv,pat)
      use tfstk
      use ffs
c ....delete (kv,pat) from 'vbump' list, if it exists.          
      use tffitcode
      parameter (meminc=100)
      logical tmatch
      character*(*) pat
      character*12 vname
      common /cordefv/ ipvbmp,nvbmp
c
      nvbmpa=nvbmp
      do 14 i=nvbmpa,1,-1
        lp=ilist(mod(i-1,2)+1,ipvbmp+(i-1)/2)
        kvv=ilist(1,lp+1)
        call mcchar(ilist(2,lp+1),vname,3)
        lsep=ifany(vname,'_',1)-1
        if((kv.eq.kvv .or. kv.eq.0) .and. tmatch(vname(1:lsep),pat)
     $       ) then
          call tfree(int8(lp))
          do 13 j=i+1,nvbmp 
            ilist(mod(j-2,2)+1,ipvbmp+(j-2)/2)=
     $           ilist(mod(j-1,2)+1,ipvbmp+(j-1)/2)
 13       continue 
          nvbmp=nvbmp-1
          if(nvbmp.le.0) then
            if(ipvbmp.ne.0) then
              call tfree(int8(ipvbmp))
              ipvbmp=0
            endif 
          elseif(mod(nvbmp,meminc).eq.0) then
            call palocx(ipvbmp,nvbmp,nvbmp)
          endif
        endif 
 14   continue
      return
      end
