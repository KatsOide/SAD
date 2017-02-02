      function sethtb(token,type,ival)
      use maccbk
      implicit none
      include 'inc/MACCODE.inc'
      character*(*) token
      integer*4 sethtb,type,ival,decloc
      sethtb=sethtb8(token,type,int8(ival))
      return
      end

      function sethtbr(token,ival)
      use maccbk
      implicit none
      include 'inc/MACCODE.inc'
      character*(*) token
      integer*8 ival
      integer*4 sethtbr
      sethtbr=sethtb8(token,icRSVD,ival)
      return
      end

