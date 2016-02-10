      logical function skipch(char,token,slen,ttype,rval,ival)
c
      implicit real*8 (a-h,o-z)
      character*(*) char
      character*(*)  token
      integer slen,ttype,ival
      real*8 rval
c
      if (token(:slen) .eq. char) then
        call gettok(token,slen,ttype,rval,ival)
        skipch=.true.
      else
        skipch=.false.
      endif
      return
      end
