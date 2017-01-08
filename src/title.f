      subroutine title$(title,where)
      use macfile
      character*(*) title,where
      integer*4 lene
c
      write(pltfl,*)'TITLE ',where,'''',title(:lene(title)),''''
      return
      end
