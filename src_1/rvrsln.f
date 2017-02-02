      subroutine rvrsln(idxl)
      use maccbk
      implicit none
      integer*4 idxl
c
      include 'inc/MACCODE.inc'
      include 'inc/MACKW.inc'
      include 'inc/MACTTYP.inc'
      include 'inc/MACFILE.inc'
      include 'inc/MACMISC.inc'
c
      integer*4 llen
      integer*8 oldpt,ioldpt
      integer*4 oldex,i
      real*8 work
c
      llen=ilist(1,idval(idxl))
      oldpt=idval(idxl)
      ioldpt=oldpt+llen+1
      do 100 i=1,llen/2
        work=rlist(oldpt+i)
        rlist(oldpt+i)=rlist(ioldpt-i)
        rlist(ioldpt-i)=work
        ilist(1,oldpt+i)=-ilist(1,oldpt+i)
        ilist(1,ioldpt-i)=-ilist(1,ioldpt-i)
 100  continue
      oldex=ilist(2,oldpt)
      if(oldex .eq. 0) return
      llen=ilist(1,oldex)
      ioldpt=oldex+llen+1
      do 200 i=1,llen/2
        work=rlist(oldex+i)
        rlist(oldex+i)=rlist(ioldpt-i)
        rlist(ioldpt-i)=work
 200  continue
      return
      end
