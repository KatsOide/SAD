      integer*4 function sethtbr(token,ival)
      use maccode
      use maccbk
      implicit none
      character*(*) token
      integer*8 ival
      sethtbr=sethtb8(token,icRSVD,ival)
      return
      end function
