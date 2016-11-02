      subroutine tftyp1(kx,l,lp,kp,lt,lfno,lpw)
      use tfstk
      use ffs, only:emx,emy,dpmax
      use ffs_pointer
      use tffitcode
      use mackw
      implicit none
      integer*8 lp
      integer*4 ioff,kx,l,kp,lt,lfno,lv,lene,lenw,lpw,lpname
      real*8 v
      character*32 autos
      character*132 vout
      character*(MAXPNAME) kw,tfkwrd
      character*8 unit
      logical*4 start
      external trim
      start=.true.
      do 10 ioff=1,100
        kw=tfkwrd(lt,ioff)
        if(kw .eq. ' ')then
          if(start)then
            call twbuf(pname(kp)(1:max(8,lpname(kp)))//'=()',
     $           lfno,7,lpw,7,1)
          else
            call twbuf(')',lfno,10,lpw,1,1)
          endif
          return
        elseif(kw .eq. '-')then
          go to 10
        endif
        if(lt .eq. icMARK)then
          if(ioff .le. ntwissfun)then
            if(l .eq. 1)then
              v=rlist(lp+ioff)
            else
              v=twiss(l,0,ioff)
              if(ioff .eq. mfitax .or. ioff .eq. mfitay
     $             .or. ioff .eq. mfitepx .or. ioff .eq. mfitepy
     $             .or. ioff .eq. mfitr2 .or. ioff .eq. mfitr3
     $             .or. ioff .eq. mfitdpx .or. ioff .eq. mfitdpy
     $             .or. ioff .eq. mfitaz .or. ioff .eq. mfitzpx
     $             .or. ioff .eq. mfitzpy)then
                v=v*direlc(l)
              endif
            endif
          elseif(ioff .eq. kytbl(kwEMIX,icMARK))then
            v=emx
          elseif(ioff .eq. kytbl(kwEMIY,icMARK))then
            v=emy
          elseif(ioff .eq. kytbl(kwDP,icMARK))then
            v=dpmax
          else
            v=rlist(lp+ioff)
          endif
        elseif(ioff .eq. ival(kx))then
          v=rlist(lp+ioff)/errk(1,l)
        else
          v=rlist(lp+ioff)
        endif
        if(v .ne. 0.d0 .or. ioff .eq. ival(kx))then
          if(kw .eq. 'ROTATE')then
            v=v*90.d0/asin(1.d0)
            unit=' DEG'
          else
            unit=' '
            if(kw .eq. 'SIGMAZ')then
              kw='SIGZ'
            endif
          endif
          vout=kw(1:lenw(kw))//' ='
     $         //autos(v)//unit(1:lene(unit))
          call trim(vout)
          if(abs(v) .gt. 1.d10 .and. index(vout,'.') .le. 0)then
            lv=lene(vout)
            vout(lv+1:lv+1)='.'
          endif
          if(start)then
            call twbuf(pname(kp)(1:max(8,lpname(kp)))//
     $           '=('//vout,lfno,7,lpw-2,7,1)
            start=.false.
          else
            call twbuf(vout,lfno,10,lpw-2,5,1)
          endif
        endif
10    continue
      return
      end
