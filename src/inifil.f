      subroutine inifil
      use maccbk
      implicit none
      include 'inc/MACFILE.inc'
      integer*4 mctaloc
c
         infl=STDIN
         outfl=STDOUT
         pltfl=STDPLT
         errfl=STDERR
         lstfl=STDLST
         inflpt=mctaloc(4)
         ilist(1,inflpt)=infl
         ilist(2,inflpt)=0
      return
      end
