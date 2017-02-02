      integer*8 ix
      byte b(8)
      equivalence (ix,b)
      ix=-1
      b(8)=iand(b(8),31)
      write(*,'(i20)')ix
      b(8)=z'80'
      write(*,*)b(8)
      end
