      subroutine pwrlat(word,wordp,latt,mult,lfno)
      use tfstk
      use ffs
      use tffitcode
      character*(*) word,wordp
      character*(MAXPNAME) name
      character cline*130,patt*80
      logical psname,abbrev,tmatch,exist,patfit
      dimension latt(2,nlat),mult(*),ls(2),lf(2)
      data cline/' '/,ip/1/
c
      if(abbrev(word,'PS_NAME','_')) then
        psname=.true.
      else
        psname=.false.
      endif
      call getwdl2(word,wordp)
      if(word.eq.' ') then
        patt='*'
        go to 11
      endif
      do 10 i=1,nlat-1
        if( tmatch(pname(latt(1,i)),wordp) ) then
          patt=wordp
          call getwdl2(word,wordp)
          goto 11
        endif
   10 continue
      patt=' '
   11 continue
      ni=ielm(wordp,exist)
      if(exist) then
        call getwdl2(word,wordp)
        nf=ielm(wordp,exist)
        if(exist) then
          call getwdl2(word,wordp)
        else
          nf=nlat
        endif
      else
        ni=1
        nf=nlat
      endif
c     print *,ni,nf
      if(ni.gt.nf) then
        ls(1)=ni
        lf(1)=nlat
        ls(2)=1
        lf(2)=nf
      else
        ls(1)=ni
        lf(1)=nf
        ls(2)=nlat
        lf(2)=1
      endif
      noel=0
      do 21 l=1,2
        do 20 i=ls(l),lf(l)
          if(i.eq.nlat) then
            patfit='***'.eq.patt
     $           .or. patt .eq. '$$$' .or. patt.eq.'*'
          elseif(i .eq. 1)then
            patfit=tmatch(pname(latt(1,i)),patt) .or. patt .eq. '^^^'
          else
            patfit=tmatch(pname(latt(1,i)),patt)
          endif
          if( patfit ) then
            noel=noel+1
            if(i.eq.nlat) then
              name='$$$'
            else
              if(psname) then
                name=pname(latt(1,i))
              else
                call elname(i,name)
              endif
            endif
            cline(ip:ip+7)=name
            ip=ip+8
            if(ip.ge.123) then
              write(lfno,'('' '',a)') cline
              cline=' '
              ip=1
            endif
          endif
   20   continue
   21 continue
      if(cline.ne.' ') then
        write(lfno,'('' '',a)') cline
        cline=' '
        ip=1
      endif
      write(lfno,'(I6,2A)') noel,' elements match ',patt(1:lene(patt))
      return
      end
