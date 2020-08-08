      subroutine tftyp1(kx,l,lp,kp,lt,lfno,lpw)
      use kyparam
      use tfstk
      use ffs, only:emx,emy,dpmax
      use ffs_pointer
      use tffitcode
      use sad_main
      use mackw
      use strbuf
      implicit none
      type (sad_comp) , pointer ::cmp
      type (sad_dlist), pointer :: klv
      type (sad_strbuf), pointer :: strb
      integer*8 lp
      integer*4 ioff,kx,l,kp,lt,lfno,lv,lene,lenw,lpw,irtc,nc,j,j1
      real*8 v,coeff
      character*32 autos1
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
            call twbuf(pname(kp)(1:max(8,lpname(kp))),'=()',
     $           lfno,7,lpw,7,1)
          else
            call twbuf(')','',lfno,10,lpw,1,1)
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
          elseif(ioff .eq. ky_EMIZ_MARK)then
            v=emz
          elseif(ioff .eq. ky_SIGZ_MARK)then
            v=sigzs
          elseif(ioff .eq. ky_SIGE_MARK)then
            v=sizedp
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
     $         //autos1(v)//unit(1:lene(unit))
          call trim(vout)
          if(v .ne. 0.d0 .or. ioff .eq. ival(kx))then
c            if(abs(v) .gt. 1.d10 .and. index(vout,'.') .le. 0
c     $           .and. v .ne. dinfinity)then
c              lv=lene(vout)
c              vout(lv+1:lv+1)='.'
c            endif
            if(start)then
              call twbuf(pname(kp)(1:max(8,lpname(kp))),
     $             '=('//vout,lfno,7,lpw-2,7,1)
              start=.false.
            else
              call twbuf(vout,'',lfno,10,lpw-2,5,1)
            endif
          endif
        elseif(ioff .eq. kytbl(kwPROF,lt))then
          if(tflistq(cmp%dvalue(ioff),klv))then
            if(start)then
              vout=pname(kp)(1:max(8,lpname(kp)))//
     $             '=('//kw(1:lenw(kw))//' ='
            else
              vout='  '//kw(1:lenw(kw))//' ='
            endif
            if(start)then
              call twbuf(vout,'',lfno,7,lpw-2,7,1)
              start=.false.
            else
              call twbuf(vout,'',lfno,10,lpw-2,1,1)
            endif
            call getstringbuf(strb,0,.true.)
            call tfconvstrb(strb,cmp%dvalue(ioff),nc,
     $           .false.,.false.,-1,'*',irtc)
            j=1
            do while(j .le. nc)
              j1=index(strb%str(j:nc),',')
              if(j1 .eq. 0)then
                call twbuf(strb%str(j:nc),'',lfno,10,lpw-2,1,1)
                j=nc+1
              else
                call twbuf(strb%str(j:j+j1-1),'',lfno,10,lpw-2,1,1)
                j=j+j1
              endif
            enddo                
            call tfreestringbuf(strb)
          endif
        endif
      enddo
      return
      end
