      subroutine tfdisp(word,idisp1,idisp2,dp00,lfno,exist)
      use tfstk
      use ffs
      use ffs_pointer
      use ffs_fit, only:scale
      use tffitcode
      implicit none
      integer*4 idisp1,idisp2,lfno,mdisp,icolm,ifany,id0,id3,idstep,
     $     lines,l,id,ielm,i
      real*8 dp00,dgam,bx0,by0,bx1,by1,bx2,by2,r,sigpp,tfchi,detr,
     $     etaxp,etapxp,sigxxp,sigxpxp,sigpxpxp,emixp,
     $     etayp,etapyp,sigyyp,sigypyp,sigpypyp,emiyp
      integer*4 lname,lb,lb1
      integer*4, parameter:: blen=ntwissfun*12+15
      real*8 pe(4)
      real*8 og(3,4)
      character*(*) word
      character*255 wordp,word1,name
      character*(blen) buff
      character*16 autofg,vout
      character*255 bname
      character*1 dir,hc
      logical*4 tfinsol
      logical exist,dpeak,seldis,abbrev,temat,mat
c     begin initialize for preventing compiler warning
      sigpp=0.d0
c     end   initialize for preventing compiler warning
      dgam=gammab(1)*dp00
      exist=.false.
      dpeak=.false.
      seldis=.false.
      mdisp=0
      icolm=0
      buff=' '
270   call getwdl2(word,wordp)
c      write(*,*)'tfdisp ',word,wordp
      if(word(1:3) .ne. '***' .and.
     1   ifany(word,'|<{*%',1) .gt. 0)then
        seldis=.true.
        word1=wordp
        go to 270
      elseif(abbrev(word,'ALL','_'))then
        idisp1=1
        idisp2=nlat
        go to 270
      elseif(abbrev(word,'A_CCELERATION','_'))then
        mdisp=5
        go to 270
      elseif(abbrev(word,'G_EOMETRY','_'))then
        mdisp=101
        go to 270
      elseif(abbrev(word,'OG_EOMETRY','_'))then
        mdisp=102
        go to 270
      elseif(abbrev(word,'O_RBIT','_'))then
        mdisp=3
        go to 270
      elseif(abbrev(word,'R_MATRIX','_'))then
        mdisp=2
        go to 270
      elseif(abbrev(word,'B_EAM','_'))then
        mdisp=4
        go to 270
      elseif(abbrev(word,'P_HYSICAL','_'))then
        mdisp=6
        go to 270
      elseif(abbrev(word,'D_UMPOPTICS','_'))then
        mdisp=10
        go to 270
      elseif(word .eq. 'Z')then
        mdisp=11
        go to 270
      elseif(abbrev(word,'E_XTREMUM','_'))then
        dpeak=.true.
        go to 270
      endif
      id0=ielm(wordp,exist)
      if(exist)then
        idisp1=id0
      else
        if(wordp .ne. ' ')then
          return
        endif
        go to 250
      endif
      call getwdl2(word,wordp)
      id0=ielm(wordp,exist)
      if(exist)then
        idisp2=id0
      else
        idisp2=idisp1
      endif
250   if(idisp1 .eq. idisp2)then
        id3=nlat
        idstep=max(nlat-idisp2,1)
      else
        id3=idisp2
        idstep=1
      endif
      bx0=twiss(max(1,idisp1-1),icolm,mfitbx)
      by0=twiss(max(1,idisp1-1),icolm,mfitby)
      bx1=twiss(idisp1,icolm,mfitbx)
      by1=twiss(idisp1,icolm,mfitby)
      lines=0
      do 200 l=idisp1,id3,idstep
        mat=temat(l,name,word1)
        lname=len_trim(name)
        if(seldis .and. .not. mat)then
          go to 200
        endif
        if(dpeak)then
          bx2=twiss(min(l+1,nlat),icolm,mfitbx)
          by2=twiss(min(l+1,nlat),icolm,mfitby)
          if(bx1 .ge. bx0 .and. bx1 .ge. bx2 .or.
     1       by1 .ge. by0 .and. by1 .ge. by2 .or.
     1       bx1 .le. bx0 .and. bx1 .le. bx2 .or.
     1       by1 .le. by0 .and. by1 .le. by2)then
            bx0=bx1
            by0=by1
            bx1=bx2
            by1=by2
          else
            bx0=bx1
            by0=by1
            bx1=bx2
            by1=by2
            go to 200
          endif
        endif
        if(l .eq. nlat)then
          id=41
        else
          id=idtypec(l)
        endif
        if(iele1(iele((l))) .gt. 0 .and.
     $       id .ne. icMARK .and. id .ne. 34)then
          if(ival(iele1(iele(l))) .gt. 0)then
            vout=autofg(rlist(latt(l)
     $           +ival(iele1(iele(l)))),'10.7')
          else
            vout=' 0'
          endif
        elseif(id .eq. icSOL)then
          vout=autofg(rlist(latt(l)+2),'10.7')
        else
          vout=' 0'
        endif
        if(mdisp .eq. 101 .or. mdisp .eq. 102)then
          if(mdisp .eq. 101)then
            hc=' '
          else
            hc='O'
          endif
          if(mod(lines,66) .eq. 0)then
            write(lfno,9022)hc,hc,hc,hc,hc,hc
9022        format(
     1           ' Element             ',a,
     $           'Gx        ',a,'Gy          ',a,'Gz    ',
     1           '        s        Length     Value',
     1           '       ',a,
     $           'Chi1       ',a,
     $           'Chi2       ',a,
     $           'Chi3      #'
     $           )
          endif
          buff(37:48)=autofg(pos(l)/scale(mfitleng),'12.6')
          if(kytbl(kwL,id) .eq. 0)then
            buff(49:58)=autofg(0.d0,'10.6')
          else
            buff(49:58)=autofg(
     $           rlist(latt(l)+kytbl(kwL,id)),'10.6')
          endif
          if(id .ne. 41 .and. id .ne. 42 .and. id .ne. 34)then
            buff(59:69)=' '//vout(1:10)
          else
            buff(59:69)=' 0'
          endif
          if(mdisp .eq. 101 .or. tfinsol(l))then
            do i=0,2
              buff(1+i*12:12+i*12)
     1           =autofg(geo(i+1,4,l)/scale(mfitgx+i),'12.6')
              buff(70+i*12:81+i*12)
     1             =autofg(tfchi(geo(1,1,l),i+1)/scale(i+mfitchi1),
     $             '12.6')
            enddo
          else
            call tforbitgeo(og,geo(1,1,l),
     $           twiss(l,icolm,mfitdx),
     $           twiss(l,icolm,mfitdpx),
     $           twiss(l,icolm,mfitdy),
     $           twiss(l,icolm,mfitdpy))
            do i=0,2
              buff(1+i*12:12+i*12)
     1           =autofg(og(i+1,4)/scale(mfitgx+i),'12.6')
              buff(70+i*12:81+i*12)
     1             =autofg(tfchi(og(1,1),i+1)/scale(i+mfitchi1),
     $             '12.6')
            enddo
          endif
          dir=' '
          if(l .ne. nlat)then
            if(direlc(l) .le. 0.d0)then
              dir='-'
            endif
          endif
          write(lfno,9024)dir,name(1:max(18,lname)),buff(1:105),
     $         mod(l-1,10000)+1
9024      format(a,a,a,1x,i5)
        elseif(mdisp .eq. 10)then
          if(l .eq. idisp1)then
            buff(1:15)=' '
            do i=1,ntwissfun
              buff((i-1)*12+16:i*12+15)=autofg(scale(i),'12.8')
            enddo
            write(lfno,'(a)')buff(1:ntwissfun*12+15)
          endif
          buff(1:12)=name
          write(buff(13:15),'(I3)')id
          do i=1,ntwissfun
            buff((i-1)*12+16:i*12+15)
     $           =autofg(twiss(l,icolm,i)/scale(i),'12.8')
          enddo
c$$$          buff((19-1)*12+16:19*12+15)=autofg(pos(l)/
c$$$     $         scale(mfitleng),'12.8')
c$$$          do i=0,2
c$$$            buff((19+i)*12+16:(20+i)*12+15)=
c$$$     $           autofg(geo(i+1,4,l)/scale(mfitgx+i),'12.8')
c$$$            buff((22+i)*12+16:(23+i)*12+15)=
c$$$     $           autofg(tfchi(geo(1,1,l),i+1)/scale(mfitchi1+i),'12.8')
c$$$          enddo
c$$$          buff((26-1)*12+16:26*12+15)=vout
          write(lfno,'(a)')buff(1:blen)
        else
          if(mod(lines,66) .eq. 0)then
            if(mdisp .eq. 2)then
              write(lfno,9021)
9021          format(
     1      '   AX      BX      NX      EX     EPX   ',
     1      ' Element    R1     R2     R3     R4   ',' ',
     1      '   AY      BY      NY      EY     EPY     DetR     #')
            elseif(mdisp .eq. 6)then
              write(lfno,9026)
 9026         format(
     1      '   AX      BX      NX      EX     EPX   ',
     1      ' Element    PEX    PEPX   PEY    PEPY ',' ',
     1      '   AY      BY      NY      EY     EPY     DetR     #')
            elseif(mdisp .eq . 3)then
              write(lfno,9041)
9041          format(
     1      '   AX      BX      NX      EX     EPX   ',
     1      ' Element    DX     DPX    DY     DPY  ',' ',
     1      '   AY      BY      NY      EY     EPY     DetR     #')
            elseif(mdisp .eq. 4)then
              write(lfno,9042)
9042          format(
     1      '  AXp     BXp     EMITXp  EXp    EPXp   ',
     1      ' Element   Sigx(mm)  Sigy(mm)  Rot(deg)',' ',
     1       ' AYp     BYp     EMITYp  EYp    EPYp     sigp     #')
            elseif(mdisp .eq. 5)then
              write(lfno,9043)
9043          format(
     1      '   AX      BX      NX      EX     EPX   ',
     1      ' Element    p(GeV)emitx(m)emity(m) DDP',' ',
     1      '   AY      BY      NY      EY     EPY      DZ      #')
            elseif(mdisp .eq. 11)then
              write(lfno,9032)
 9032         format(
     1      '   AZ      BZ      NZ                   ',
     1      ' Element   Length   Value      s(m)   ',' ',
     1      '   ZX      ZPX     ZY      ZPY            DetR     #')
            else
              write(lfno,9031)
9031          format(
     1      '   AX      BX      NX      EX     EPX   ',
     1      ' Element   Length   Value      s(m)   ',' ',
     1      '   AY      BY      NY      EY     EPY     DetR     #')
            endif
          endif
          if(mdisp .eq. 11)then
            buff( 1: 8)=autofg(twiss(l,icolm,mfitaz)/
     $           scale(mfitaz),'8.5')
            buff( 9:16)=autofg(twiss(l,icolm,mfitbz)/
     $           scale(mfitbz),'8.5')
            buff(17:24)=autofg(twiss(l,icolm,mfitnz)/
     $           scale(mfitnz),'8.5')
          elseif(mdisp .ne. 4)then
            buff( 1: 8)=autofg(twiss(l,icolm,mfitax)/
     $           scale(mfitax),'8.5')
            buff( 9:16)=autofg(twiss(l,icolm,mfitbx)/
     $           scale(mfitbx),'8.5')
            buff(17:24)=autofg(twiss(l,icolm,mfitnx)/
     $           scale(mfitnx),'8.5')
            buff(25:32)=autofg(twiss(l,icolm,mfitex)/
     $           scale(mfitex),'8.5')
            buff(33:40)=autofg(twiss(l,icolm,mfitepx)/
     $           scale(mfitepx),'8.5')
          else
            call ffs_init_sizep
            if(.not. updatesize .or. sizedp .ne. dpmax)then
              call tfsize
            endif
            sigpp=beamsize(21,l)
            etaxp=beamsize(16,l)/sigpp
            etapxp=beamsize(17,l)/sigpp
            sigxxp=beamsize(1,l)-etaxp**2*sigpp
            sigxpxp=beamsize(2,l)-etaxp*etapxp*sigpp
            sigpxpxp=beamsize(3,l)-etapxp**2*sigpp
            emixp=sqrt(sigxxp*sigpxpxp-sigxpxp**2)
            buff( 1: 8)=autofg(-sigxpxp/emixp/scale(mfitax),'8.5')
            buff( 9:16)=autofg(sigxxp/emixp/scale(mfitbx),'8.5')
            buff(17:24)=autofg(emixp,'8.5')
            buff(25:32)=autofg(etaxp/scale(mfitex),'8.5')
            buff(33:40)=autofg(etapxp/scale(mfitepx),'8.5')
          endif
          if(l .ne. nlat .and.
     $         direlc(l) .le. 0.d0)then
            bname(1:max(10,lname+2))=' -'//name(1:lname)
            lname=max(10,lname+2)
          else
            if(lname .le. 8)then
              bname(1:10)='  '//name(1:lname)
              lname=10
            else
              bname(1:lname+1)=' '//name(1:lname)
              lname=lname+1
            endif
          endif
          buff(79:79)=' '
          if(mdisp .eq. 2)then
            do i=0,3
              buff(51+i*7:57+i*7)
     $             =autofg(twiss(l,icolm,mfitr1+i)
     $             /scale(mfitr1+i),'7.4')
            enddo
          elseif(mdisp .eq. 6)then
            call tgetphysdisp(l,pe)
            do i=0,3
              buff(51+i*7:57+i*7)
     $             =autofg(pe(i+1)
     $             /scale(mfitpex+i),'7.4')
            enddo
          elseif(mdisp .eq. 3)then
            do 120 i=0,3
              buff(51+i*7:78+i*7)
     $             =autofg(twiss(l,icolm,mfitdx+i)/
     $             scale(mfitdx+i),'7.4')
120         continue
          elseif(mdisp .eq. 4)then
            buff(51:60)=autofg(sqrt(beamsize(1,l))*1.d3,'10.8')
            buff(61:70)=autofg(sqrt(beamsize(6,l))*1.d3,'10.8')
            buff(71:79)=autofg(
     1           atan2(-2.d0*beamsize(4,l),beamsize(1,l)-beamsize(6,l))
     $           *90.d0/pi,'9.4')
          elseif(mdisp .eq. 5)then
            buff(51:58)=autofg((gammab(l)+dgam)*amass/1.d9,'8.6')
            r=(gammab(1)+dgam)/(gammab(l)+dgam)
            buff(59:65)=autofg(emx*r,'7.5')
            buff(66:72)=autofg(emy*r,'7.5')
            buff(73:79)=autofg(twiss(l,icolm,mfitddp)
     $           /scale(mfitddp),'7.4')
          else
            if(kytbl(kwL,id) .eq. 0)then
              buff(51:58)=autofg(0.d0,'8.5')
            else
              buff(51:68)=autofg(
     $             rlist(latt(l)+kytbl(kwL,id)),'8.5')
            endif
            if(id .ne. 41 .and. id .ne. 42 .and. id .ne. 34)then
              buff(59:68)=vout(1:10)
            else
              buff(59:68)=' 0'
            endif
            buff(69:79)=autofg(pos(l)/scale(mfitleng),'11.6')
          endif
          if(mdisp .eq. 11)then
            buff(80:87)=autofg(twiss(l,icolm,mfitzx)/
     $           scale(mfitzx),'8.5')
            buff(88:95)=autofg(twiss(l,icolm,mfitzpx)/
     $           scale(mfitzpx),'8.5')
            buff(96:103)=autofg(twiss(l,icolm,mfitzy)/
     $           scale(mfitzy),'8.5')
            buff(104:111)=autofg(twiss(l,icolm,mfitzpy)/
     $           scale(mfitzpy),'8.5')
            buff(120:126)=autofg(twiss(l,icolm,mfitdetr)/
     $           scale(mfitdetr),'7.4')
          elseif(mdisp .ne. 4)then
            buff(80:87)=autofg(twiss(l,icolm,mfitay)/
     $           scale(mfitay),'8.5')
            buff(88:95)=autofg(twiss(l,icolm,mfitby)/
     $           scale(mfitby),'8.5')
            buff(96:103)=autofg(twiss(l,icolm,mfitny)/
     $           scale(mfitny),'8.5')
            buff(104:111)=autofg(twiss(l,icolm,mfitey)/
     $           scale(mfitey),'8.5')
            buff(112:119)=autofg(twiss(l,icolm,mfitepy)/
     $           scale(mfitepy),'8.5')
            if(mdisp .eq. 5)then
              buff(120:126)=autofg(twiss(l,icolm,mfitdz)/
     $             scale(mfitdz),'7.4')
            else
              detr=twiss(l,icolm,mfitdetr)
              buff(120:126)=autofg(detr,'7.4')
            endif
          else
            etayp=beamsize(18,l)/sigpp
            etapyp=beamsize(19,l)/sigpp
            sigyyp=beamsize(6,l)-etayp**2*sigpp
            sigypyp=beamsize(9,l)-etayp*etapyp*sigpp
            sigpypyp=beamsize(10,l)-etapyp**2*sigpp
            emiyp=sqrt(sigyyp*sigpypyp-sigypyp**2)
            buff(80:87)=autofg(-sigypyp/emiyp/scale(mfitay),'8.5')
            buff(88:95)=autofg(sigyyp/emiyp/scale(mfitby),'8.5')
            buff(96:103)=autofg(emiyp,'8.5')
            buff(104:111)=autofg(etayp/scale(mfitey),'8.5')
            buff(112:119)=autofg(etapyp/scale(mfitepy),'8.5')
            buff(120:126)=autofg(sqrt(sigpp),'7.4')
          endif
          write(buff(127:131),'(i5)')mod(l-1,100000)+1
          call tfsqueezespace(buff(51:131),lname-10,lb)
          call tfsqueezespace(buff(1:40),lb+lname-91,lb1)
          write(lfno,9020)buff(1:lb1),bname(1:lname),buff(51:50+lb)
9020      format(a,a,a)
        endif
        lines=lines+1
200   continue
      return
      end

      subroutine tfsqueezespace(str,n,ls)
      implicit none
      integer*4 n,i,l,nsq,i1,ls
      character*(*) str
      l=len(str)
      if(n .le. 0)then
        ls=l
        return
      endif
      nsq=0
      i1=1
      do i=2,l
        i1=i1+1
        if(str(i-1:i-1) .eq. ' ' .and.
     $       str(i:i) .eq. ' ' .and. nsq .lt. n)then
          nsq=nsq+1
          i1=i1-1
        endif
        str(i1:i1)=str(i:i)
      enddo
      do i=l-nsq+1,l
        str(i:i)=' '
      enddo
      ls=l-nsq
      return
      end
