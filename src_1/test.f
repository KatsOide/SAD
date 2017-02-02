      type a
      integer*4, pointer :: ip
      end type

      type (a) x
      integer*4, target :: ia(10)
      x%ip=>ia(1)
      ia(1)=3
      write(*,*)x%ip,sizeof(x)

      stop
      end


