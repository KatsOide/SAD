      logical*4 function skipch(char,token,slen,ttype,rval,ival)
c
      implicit none
      character*(*) char
      character*(*)  token
      integer*4 slen,ttype,ival
      real*8 rval
c
      if (token(:slen) == char) then
        call gettok(token,slen,ttype,rval,ival)
        skipch=.true.
      else
        skipch=.false.
      endif
      return
      end
