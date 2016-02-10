      subroutine pkill(word,wordp,latt,mult,imon,nmon,istr,nstr,
     1                 lfno)
      use ffs
      use tffitcode
      character*(*) word,wordp
      logical mhogal,exist
      dimension latt(2,nlat),mult(*)
      dimension imon(nmona,4),istr(nstra,4)
      external pack
      include 'inc/common.inc'
c
      call getwdl2(word,wordp)
      if(word.eq.'M') then
        if(nmon.eq.0) return
        call getwdl2(word,wordp)
        ni=ielm(latt,wordp,1,mult,exist)
        if(exist) then
          call getwdl2(word,wordp)
          nf=ielm(latt,wordp,1,mult,exist)
          if(exist) then
            call getwdl2(word,wordp)
          else
            nf=ni
          endif
          do 10 i=1,nmon
            if(mhogal(ni,nf,imon(imon(i,2),1))) then
              imon(imon(i,2),3)=1
            endif
   10     continue
          nmo=nmon
          call pack(imon(1,2),imon(1,3),nmon,nmon)
          write(lfno,'(2(A,I4))')
     1              '  BPM available:',nmon,' in ',nmo
          nmonact=nmon
        else
          return
        endif
      elseif(word.eq.'C') then
        if(nmon.eq.0) return
        call getwdl2(word,wordp)
        ni=ielm(latt,wordp,1,mult,exist)
        if(exist) then
          call getwdl2(word,wordp)
          nf=ielm(latt,wordp,1,mult,exist)
          if(exist) then
            call getwdl2(word,wordp)
          else
            nf=ni
          endif
          do 20 i=1,nstr
            if(mhogal(ni,nf,istr(istr(i,2),1))) then
              istr(istr(i,2),3)=1
            endif
   20     continue
          nst=nstr
          call pack(istr(1,2),istr(1,3),nstr,nst)
          write(lfno,'(2(A,I4))')
     1             '  Correction dipole available :',nstr,' in ',nst
        else
          return
        endif
      endif
      return
      end
