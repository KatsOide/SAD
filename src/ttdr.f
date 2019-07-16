      subroutine ttdr(lfno,err)
      use tfcsi,only:cssetp
      implicit none
      logical*4 err,exist
      character*80 word,display
      character*12 autofg
      real*8 getva
      integer*4 system,lfno,ir,next,lfn
      goto 9000
c$$$      lfn=getva(exist)
c$$$      if(exist)then
c$$$        write(word,'(a,I2.2)')'ftn',lfn
c$$$      else        
c$$$        call peekwd(word,next)
c$$$        if(word .eq. ' ')then
c$$$          call termes(lfno,'Missing file in TDR',' ')
c$$$          go to 9000
c$$$        endif
c$$$        call cssetp(next)
c$$$      endif
c$$$      call get_environment_variable('DISPLAY',display)
c$$$      if(display .ne. ' ')then
c$$$        ir=system('tdr -v X '//word(1:len_trim(word)))
c$$$      else
c$$$        ir=system('tdr '//word)
c$$$      endif
c$$$      if(ir .lt. 0)then
c$$$        call termes(lfno,'Error in tdr, code =',
c$$$     $       autofg(dble(ir),'S12.0'))
c$$$        err=.true.
c$$$      else
c$$$        err=.false.
c$$$      endif
c$$$      return
 9000 err=.true.
      call termes(lfno,
     $  'A message from SAD: TDR command is now obsolete.',
     $  ' ')
      return
      end

