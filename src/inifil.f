      subroutine inifil
      use maccbk
      use macfile
      use tfcsi
      implicit none
      integer*8 ktcaloc
c
         infl=STDIN
         outfl=STDOUT
         pltfl=STDPLT
         errfl=STDERR
         lstfl=STDLST
         inflpt=ktcaloc(5)
         klist(inflpt)=infl
         klist(inflpt+1)=0
         lfni=infl
         call csinit(0,1,'!')
      return
      end
