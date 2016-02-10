      subroutine title$(title,where)
      character*(*) title,where
      include 'inc/MACFILE.inc'
      integer*4 lene
c
      write(pltfl,*)'TITLE ',where,'''',title(:lene(title)),''''
      return
      end
