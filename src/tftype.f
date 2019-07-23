      subroutine tftype(lfno,word)
      use kyparam
      use tfstk
      use ffs_pointer
      use ffs_flag
      use tffitcode
      use mackw
      implicit none
      integer*4 ntyp
      parameter (ntyp=18)
      integer*8 kavl
      integer*4 lfno,k1,ltyp(ntyp)
      integer*4 j,kx,lt,kp,notchar,ifany,lpw,itfgetrecl,
     $     nl,kkk,irtc
      character*(*) word
      character*(MAXPNAME) name1,tfkwrd
      character*80 patt
      logical*4 exist,exist1,tmatch,start,mulc
      data ltyp /icDRFT,
     $     icBEND,icQUAD,icSEXT,icOCTU,icDECA,icDODECA,icMULT,icSOL,
     $     icCAVI,icTCAV,
     $     icMAP,icBEAM,icINS,icCOORD,icAPRT,icMONI,icMARK/
      exist1=.false.
      lpw=min(131,itfgetrecl())
      name1=' '
1     call getwdl(word)
      if(word .eq. 'ALL' .or. word .eq. ' ' .and. .not. exist1)then
        patt='*'
      elseif(word .eq. ' ')then
        patt=' '
        exist=.false.
        go to 21
      elseif(notchar(word,'*',1) .eq. 0)then
        patt='*'
        call getwdl(word)
      else
        patt=word
      endif
      exist=.false.
      if(patt .eq. ' ')then
        go to 21
      endif
      call tfgetlineps(patt,len_trim(patt),nl,kavl,1,irtc)
      if(irtc .ne. 0 .or. nl .le. 0)then
        go to 21
      endif
      do j=1,ntyp
        lt=ltyp(j)
        mulc=kytbl(kwMAX,lt) .gt. ky_MAX_DRFT
        start=.true.
        do kkk=1,nl
          k1=int(rlist(kavl+kkk))
          kx=klp(k1)
c     Note: Skip no-head multiple elements
c     *     klp(iele1(kx)) == kx if singlet or head of multipole elements
c          if(klp(iele1(kx)) .ne. kx)cycle
          kp=idelc(kx)
          if(tmatch(pname(kp),patt))then
            exist=.true.
            if(idtype(kp) .eq. lt)then
              if(start)then
                write(lfno,*)';'
                call twbuf(' ',lfno,1,lpw,0,0)
                call twbuf(tfkwrd(lt,0),lfno,1,lpw,7,1)
                start=.false.
              endif
              call tftyp1(iele1(kx),kx,latt(kx),kp,lt,lfno,lpw)
              if(mulc)then
                call twbuf(' ',lfno,10,lpw,0,-1)
              endif
            endif
          endif
        enddo
        if(.not. start .and. .not. mulc)then
          call twbuf(' ',lfno,10,lpw,0,-1)
        endif
      enddo
21    exist1=exist1 .or. exist
c      if(.not. exist1)then
c        patt='*'
c        go to 2
c      endif
      if(exist .and. notchar(patt,'*',1) .ne. 0)then
        go to 1
      else
        exist=ifany(patt,'*%{|',1) .ne. 0
        write(lfno,*)';'
        return
      endif
      end

      subroutine tftyp1(kx,l,lp,kp,lt,lfno,lpw)
      use kyparam
      use tfstk
      use ffs, only:emx,emy,dpmax
      use ffs_pointer
      use tffitcode
      use ffs_flag
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
          if(.not. k64)then
            if(ioff .eq. mfitaz .or. ioff .eq. mfitbz .or.
     $           ioff .eq. mfitnz)then
              v=0.d0
            endif
          endif
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
c          if(kw .eq. 'SIGMAZ')then
c            kw='SIGZ'
c          endif
        endif
        if(real)then
          v=v*coeff
          vout=kw(1:lenw(kw))//' ='
     $         //autos(v)//unit(1:lene(unit))
          call trim(vout)
          if(v .ne. 0.d0 .or. ioff .eq. ival(kx))then
            if(abs(v) .gt. 1.d10 .and. index(vout,'.') .le. 0
     $           .and. v .ne. dinfinity)then
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
        elseif(ioff .eq. kytbl(kwPROF,lt) .and. k64)then
          if(tflistq(cmp%dvalue(ioff),klv))then
            if(start)then
              vout=pname(kp)(1:max(8,lpname(kp)))//
     $             '=('//kw(1:lenw(kw))//' ='
            else
              vout='  '//kw(1:lenw(kw))//' ='
            endif
            if(start)then
              call twbuf(vout,lfno,7,lpw-2,7,1)
              start=.false.
            else
              call twbuf(vout,lfno,10,lpw-2,1,1)
            endif
            call getstringbuf(strb,0,.true.)
            call tfconvstrb(strb,cmp%dvalue(ioff),nc,
     $           .false.,.false.,-1,'*',irtc)
            j=1
            do while(j .le. nc)
              j1=index(strb%str(j:nc),',')
              if(j1 .eq. 0)then
                call twbuf(strb%str(j:nc),lfno,10,lpw-2,1,1)
                j=nc+1
              else
                call twbuf(strb%str(j:j+j1-1),lfno,10,lpw-2,1,1)
                j=j+j1
              endif
            enddo                
            call tfreestringbuf(strb)
          endif
        endif
      enddo
      return
      end
