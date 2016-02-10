      integer*4 function itopenbuf(irtc)
      use tfrbuf
      implicit none
      integer*4 irtc,in,itfmessage,nextfn
      in=nextfn(0)
      if(in .eq. 0)then
        itopenbuf=0
        irtc=itfmessage(999,'General::toomany','"files opened"')
        return
      endif
      open(in,file='/dev/null',status='UNKNOWN')
      call tfreadbuf(irbinit,in,int8(2),int8(0),0,' ')
      itopenbuf=in
      irtc=0
      return
      end
