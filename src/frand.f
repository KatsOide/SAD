      real*8 function frand(idx)
      use maccbk
      implicit none
      integer idx
      integer ig0
      real*8 rgetgl1,gaussn
c
      ig0=nint(rgetgl1('$SEED'))
c     go to (1000) ilist(2,idx)   chnged as follows 04/09/88
      go to (1000) idval(ilist(2,idx))
      call errmsg('rand',
     &            'unsupported distribution type',0,0)
c for debug
c      print *,'rand ',idval(ilist(2,idx) )
c end debug
      frand=0.0d0
      return
c  type =1 normal(gaussian) distribution
 1000 continue
        if (ilist(1,idx) .eq. 1) then
          frand=rlist(idx+1)*gaussn(ig0)
        else
          frand=rlist(idx+1)+rlist(idx+2)*gaussn(ig0)
        endif
        call rsetgl1('$SEED',dble(ig0))
        return
      end
