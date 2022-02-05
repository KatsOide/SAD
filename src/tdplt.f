      subroutine tdplt$(lvar,bottom,left,top,right,limit)
      integer lvar
      character*(*) bottom,left,top,right
      logical limit
c      call ptrace('tdplt',1)
       call tdcmd$('SET SYMBOL .M; PLOT'
     &                   ,lvar,bottom,left,top,right,limit)
c      call ptrace('tdplt',-1)
      return
      end
