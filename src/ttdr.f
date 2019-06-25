      subroutine ttdr(lfno,err)
      use tfcsi,only:cssetp
      implicit none
      logical*4 err,exist
      character*80 word,display
      character*12 autofg
      real*8 getva
      integer*4 system,lfno,ir,next,lfn
      lfn=getva(exist)
      if(exist)then
        write(word,'(a,I2.2)')'ftn',lfn
      else        
        call peekwd(word,next)
        if(word .eq. ' ')then
          call termes(lfno,'Missing file in TDR',' ')
          go to 9000
        endif
        call cssetp(next)
      endif
      call get_environment_variable('DISPLAY',display)
      if(display .ne. ' ')then
        ir=system('tdr -v X '//word(1:len_trim(word)))
      else
        ir=system('tdr '//word)
      endif
      if(ir .lt. 0)then
        call termes(lfno,'Error in tdr, code =',
     $       autofg(dble(ir),'S12.0'))
        err=.true.
      else
        err=.false.
      endif
      return
 9000 err=.true.
      call termes(lfno,'Usage: TDR {filename|file_number}',' ')
      return
      end

