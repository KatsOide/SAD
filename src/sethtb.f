      function sethtb(token,type,ival)
      use maccbk
      use mackw
      implicit none
      character*(*) token
      integer*4 sethtb,type,ival
      sethtb=sethtb8(token,type,int8(ival))
      return
      end

      function sethtbr(token,ival)
      use maccbk
      use mackw
      implicit none
      character*(*) token
      integer*8 ival
      integer*4 sethtbr
      sethtbr=sethtb8(token,icRSVD,ival)
      return
      end

