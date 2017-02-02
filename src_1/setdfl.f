      Subroutine setdfl(cmdidx,vidx,val)
      use maccbk
      implicit none
      integer*4 cmdidx,vidx
      real*8 val
c.....cmdidx brings pointer to command name
      include 'inc/MACMISC.inc'
      include 'inc/MACCODE.inc'
      include 'inc/MACTTYP.inc'
      include 'inc/MACVAR.inc'
      integer*8 argdfp
      integer*4 mtaloc,mctaloc
C
      integer*4 argp,i,offset
c
c.......for debug
c       print *,'setdfl>',pname(cmdidx),pname(vidx),val
c.......end debug
        argdfp=idval(cmdidx)
        if(ilist(2,argdfp) .le. 0) then
          argp=mctaloc(ilist(1,argdfp))
          ilist(2,argdfp)=argp
          ilist(1,argp)=ilist(1,argdfp)-1
          do 1100 i=1,ilist(1,argp)
            ilist(1,argp+i)=ilist(1,argdfp+i+1)
            ilist(2,argp+i)=0
 1100     continue
        else
          argp=ilist(2,argdfp)
        endif
        offset=0
        do 2100 i=1,ilist(1,argdfp)
          if(ilist(1,argp+i) .eq. vidx) then
             offset=i
             go to 2200
          endif
 2100   continue
c
 2200   continue
        if (idval(vidx) .eq. VarInt) then
          ilist(2,argp+offset)=int(val)
c.......for debug
c       print *,'setdfl<',pname(vidx),offset,ilist(2,argp+offset)
c.......end debug
        else if (idval(vidx) .eq. VarRl) then
          ilist(2,argp+offset)=mtaloc(4)
          rlist(ilist(2,argp+offset))=val
c.......for debug
c       print *,'setdfl<',pname(vidx),offset
c    &                   ,rlist(ilist(2,argp+offset))
c.......end debug
        else
          call errmsg('setdfl','Cannot initialize'//pname(vidx),0,4)
        endif
        return
        end
