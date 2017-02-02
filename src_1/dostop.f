      subroutine dostop(dummy)
      implicit none
      real*8 dummy
      character*(*) msg
      parameter(msg='Good Bye!')
c
      call errmsg('dostop', msg,0,0)
      stop
      end
