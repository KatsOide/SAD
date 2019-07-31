      subroutine pgflag(isp1,kx,irtc)
      use tfstk
      use ffs
      use tffitcode
      implicit none
      integer*8 kx
      logical*4 exist,v,fl,pflogi
      integer*4 isp1,irtc,itfmessage
      integer*4 narg,nc
      character*16 name,tfgetstr

      narg=isp-isp1
      if(narg .gt. 2) go to 9010
      if(ktfnonstringq(ktastk(isp1+1)))then
        go to 9010
      endif
      irtc=0
      name=tfgetstr(ktastk(isp1+1),nc)
      if(narg .eq. 1)then
        v=pflogi(name,flags,fname,sino,nflag,exist,.false.,.false.)
        if(exist)then
          if(v)then
            kx=ktftrue
          else
            kx=0
          endif
        else
          print *,' ???Non-existent flag'
          go to 9010
        endif
      else
        if(ktfrealq(ktastk(isp)))then
          if(rtastk(isp) .eq. 0.d0)then
            fl=.false.
          else
            fl=.true.
          endif
          v=pflogi(name,flags,fname,sino,nflag,exist,.true.,fl)
          if(exist)then
            if(v)then
              kx=ktftrue
            else
              kx=0
            endif
          else
            print *,' ???Non-existent flag'
            go to 9010
          endif
        else
          go to 9010
        endif
      endif
      return
 9010 irtc=itfmessage(9,'FFS::undefflag',' ')
      return
      end

      logical*4 function pflogi(word,flg,fname,sino,nflag,exist,setmode,
     $v)
      implicit none
      integer*4 nflag,i
      character*(*) word,fname(nflag),sino(nflag)
      character*32 fn
      logical*4 flg(nflag),exist,setmode,v
      pflogi=.false.
      exist=.true.
      if(word .eq. ' ')then
        exist=.false.
        pflogi=.false.
      endif
      if(setmode) then
        pflogi=v
        do i=1,nflag
          fn=fname(i)
          if(word .eq. fn)then
            flg(i)=v
            return
          elseif(word .eq. 'NO'//fn)then
            flg(i)=.not.v
            return
          elseif(word .eq. sino(i))then
            flg(i)=.not.v
            return
          endif
        enddo
        exist=.false.
      else
        do i=1,nflag
          fn=fname(i)
          if(word .eq. fn)then
            pflogi=flg(i)
            return
          elseif(word .eq. 'NO'//fn)then
            pflogi=.not. flg(i)
            return
          elseif(word .eq. sino(i))then
            pflogi=.not. flg(i)
            return
          endif
        enddo
        exist=.false.
      endif
      return
      end
