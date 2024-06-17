      subroutine tfcalc(word,nlist,icalc,ncalc,mfpnta,mfpnta1,calc)
      use tfstk
      use ffs
      use tffitcode
      implicit none
      integer*4 ncalc,mfpnta,mfpnta1,i
      integer*4 icalc(3,*)
      character*8 nlist(mfit),word1
      character*(*) word
      logical calc,err
      calc=.true.
400   call getwdl(word)
      if(word(1:1) == ' ')then
        return
      endif
      if(index(word,'-') == len_trim(word))then
        word1=word(1:len_trim(word)-1)
        do 410 i=1,mfit
          if(nlist(i) == word1)then
            call txcalc(icalc,ncalc,mfpnta,mfpnta1,i,.false.,err)
            calc=.false.
            go to 400
          endif
410     continue
        return
      else
        do 510 i=1,mfit
          if(nlist(i) == word)then
            call txcalc(icalc,ncalc,mfpnta,mfpnta1,i,.true.,err)
            calc=.false.
            go to 400
          endif
510     continue
        return
      endif
      end
