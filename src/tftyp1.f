      subroutine tftyp1(kx,l,lp,kp,lt,lfno,lpw)
      use kyparam
      use tfstk
      use ffs, only:emx,emy,dpmax
      use ffs_pointer
      use tffitcode
      use sad_main
      use mackw
      implicit none
      type (sad_comp) , pointer ::cmp
      type (sad_rlist), pointer :: klv
      integer*8 lp
      integer*4 ioff,kx,l,kp,lt,lfno,lv,lene,lenw,lpw,lpname,i
      real*8 v,coeff
      character*32 autos
      character*132 vout
      character*(MAXPNAME) kw,tfkwrd
      character*8 unit
      logical*4 start,real,list
      external trim
      call loc_comp(lp,cmp)
      start=.true.
      do ioff=1,100
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
          cycle
        endif
        coeff=1.d0
        list=.false.
        if(lt .eq. icMARK)then
          if(ioff .le. ntwissfun)then
            if(l .eq. 1)then
              v=cmp%value(ioff)
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
          elseif(ioff .eq. ky_EMIX_MARK)then
            v=emx
          elseif(ioff .eq. ky_EMIY_MARK)then
            v=emy
          elseif(ioff .eq. ky_DP_MARK)then
            v=dpmax
          else
            v=cmp%value(ioff)
          endif
          real=.true.
        else
          real=ktfrealq(cmp%dvalue(ioff),v)
          if(ioff .eq. ival(kx))then
            coeff=1.d0/errk(1,l)
          endif
        endif
        if(kw .eq. 'ROTATE')then
          coeff=coeff*90.d0/asin(1.d0)
          unit=' DEG'
        else
          unit=' '
          if(kw .eq. 'SIGMAZ')then
            kw='SIGZ'
          endif
        endif
        if(real)then
          v=v*coeff
          vout=kw(1:lenw(kw))//' ='
     $         //autos(v)//unit(1:lene(unit))
          call trim(vout)
          if(v .ne. 0.d0 .or. ioff .eq. ival(kx))then
            if(abs(v) .gt. 1.d10 .and. index(vout,'.') .le. 0)then
              lv=lene(vout)
              vout(lv+1:lv+1)='.'
            endif
            if(start)then
              call twbuf(pname(kp)(1:max(8,lpname(kp)))//
     $             '=('//vout,lfno,7,lpw-2,7,1)
              start=.false.
            else
              call twbuf(vout,lfno,10,lpw-2,5,1)
            endif
          endif
        else
          if(tfreallistq(cmp%dvalue(ioff),klv))then
            do i=1,klv%nl
              if(i .eq. 1)then
                if(start)then
                  vout=pname(kp)(1:max(8,lpname(kp)))//
     $                 '=('//kw(1:lenw(kw))//' ={'//
     $                 autos(klv%rbody(i)*coeff)
                else
                  vout=kw(1:lenw(kw))//' ={'//
     $                 autos(klv%rbody(i)*coeff)
                endif
              else
                vout=autos(klv%rbody(i)*coeff)
              endif
              lv=lene(vout)
              if(i .ne. klv%nl)then
                vout(lv+1:lv+1)=','
              else
                vout(lv+1:lv+1)='}'
              endif
              if(start)then
                call twbuf(vout,lfno,7,lpw-2,7,1)
                start=.false.
              elseif(i .eq. 1)then
                call twbuf(vout,lfno,10,lpw-2,5,1)
              else
                call twbuf(vout,lfno,10,lpw-2,1,1)
              endif
            enddo
            call twbuf(unit(1:lene(unit)),lfno,10,lpw-2,1,1)
          endif
        endif
      enddo
      return
      end
