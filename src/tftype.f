      subroutine tftype(lfno,word)
      use kyparam
      use tfstk
      use ffs, only:nelvx
      use ffs_pointer
      use ffs_flag
      use tffitcode
      use mackw
      implicit none
      integer*4 ntyp
      parameter (ntyp=18)
      integer*8 kavl
      integer*4 ,intent(in):: lfno
      integer*4 j,kx,lt,kp,notchar,ifany,lpw,itfgetrecl,
     $     nl,kkk,irtc,k1,ltyp(ntyp)
      character*(*) ,intent(out):: word
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
      if(word == 'ALL' .or. word == ' ' .and. .not. exist1)then
        patt='*'
      elseif(word == ' ')then
        patt=' '
        exist=.false.
        go to 21
      elseif(notchar(word,'*',1) == 0)then
        patt='*'
        call getwdl(word)
      else
        patt=word
      endif
      exist=.false.
      if(patt == ' ')then
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
          kx=nelvx(k1)%klp
c     Note: Skip no-head multiple elements
c     *     klp(iele1(kx)) == kx if singlet or head of multipole elements
c          if(klp(iele1(kx)) .ne. kx)cycle
          kp=idelc(kx)
          if(tmatch(pname(kp),patt))then
            exist=.true.
            if(idtype(kp) == lt)then
              if(start)then
                write(lfno,*)';'
                call twbuf(' ','',lfno,1,lpw,0,0)
                call twbuf(tfkwrd(lt,0),'',lfno,1,lpw,7,1)
                start=.false.
              endif
              call tftyp1(iele1(kx),kx,latt(kx),kp,lt,lfno,lpw)
              if(mulc)then
                call twbuf(' ','',lfno,10,lpw,0,-1)
              endif
            endif
          endif
        enddo
        if(.not. start .and. .not. mulc)then
          call twbuf(' ','',lfno,10,lpw,0,-1)
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
      use ffs, only:emx,emy,emz,dpmax,sizedp,sigzs,nelvx
      use ffs_pointer
      use tffitcode
      use ffs_flag
      use sad_main
      use mackw
      use strbuf
      use macmath,only:m_pi
      implicit none
      type (sad_comp) , pointer ::cmp
      type (sad_dlist), pointer :: klv
      type (sad_strbuf), pointer :: strb
      integer*8 lp
      integer*4 ioff,kx,l,kp,lt,lfno,lenw,lpw,irtc,nc,j,j1
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
        if(kw == ' ')then
          if(start)then
            call twbuf(pname(kp)(1:max(8,lpname(kp))),'=()',
     $           lfno,7,lpw,7,1)
          else
            call twbuf(')','',lfno,10,lpw,1,1)
          endif
          return
        elseif(kw == '-')then
          cycle
        endif
        coeff=1.d0
        list=.false.
        if(lt == icMARK)then
          if(ioff .le. ntwissfun)then
            if(l == 1)then
              v=cmp%value(ioff)
            else
              v=twiss(l,0,ioff)
              if(ioff == mfitax .or. ioff == mfitay
     $             .or. ioff == mfitepx .or. ioff == mfitepy
     $             .or. ioff == mfitr2 .or. ioff == mfitr3
     $             .or. ioff == mfitdpx .or. ioff == mfitdpy
     $             .or. ioff == mfitaz .or. ioff == mfitzpx
     $             .or. ioff == mfitzpy)then
                v=v*direlc(l)
              endif
            endif
          elseif(ioff == ky_EMIX_MARK)then
            v=emx
          elseif(ioff == ky_EMIY_MARK)then
            v=emy
          elseif(ioff == ky_EMIZ_MARK)then
            v=emz
          elseif(ioff == ky_SIGZ_MARK)then
            v=sigzs
          elseif(ioff == ky_SIGE_MARK)then
            v=sizedp
          elseif(ioff == ky_DP_MARK)then
            v=dpmax
          else
            v=cmp%value(ioff)
          endif
          real=.true.
          if(.not. k64)then
            if(ioff == mfitdetr .or. ioff == mfitbz .or.
     $           ioff == mfitnz .or. ioff == mfitdz
     $           .or. ioff == ky_EMIZ_MARK
     $           .or. ioff == ky_SIGZ_MARK
     $           .or. ioff == ky_SIGE_MARK
     $           )then
              v=0.d0
            endif
          endif
        else
          real=ktfrealq(cmp%dvalue(ioff),v)
          if(ioff == nelvx(kx)%ival)then
            coeff=1.d0/errk(1,l)
          endif
        endif
        if(kw(1:7 ) == 'ROTATE ')then
          coeff=coeff*180.d0/m_pi
          unit=' DEG'
        elseif(kw(1:7) == 'SIGMAZ ' .and. .not. k64)then
          kw(1:6)='SIGZ  '
        else
          unit=' '
        endif
        if(real)then
          v=v*coeff
          vout=kw(1:lenw(kw))//' ='//autos1(v)
          call trim(vout)
          if(v .ne. 0.d0 .or. ioff == nelvx(kx)%ival)then
c            if(abs(v) .gt. 1.d10 .and. index(vout,'.') .le. 0
c     $           .and. v .ne. dinfinity)then
c              lv=len_trim(vout)
c              vout(lv+1:lv+1)='.'
c            endif
            if(start)then
              call twbuf(pname(kp)(1:max(8,lpname(kp)))//
     $             '=('//vout,unit(:len_trim(unit)),lfno,7,lpw-2,7,1)
              start=.false.
            else
              call twbuf(vout,unit(:len_trim(unit)),lfno,10,lpw-2,5,1)
            endif
          endif
        elseif(ioff == kytbl(kwPROF,lt) .and. k64)then
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
              if(j1 == 0)then
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
