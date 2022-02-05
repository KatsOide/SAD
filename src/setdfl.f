      Subroutine setdfl(cmdidx,vidx,val)
      use maccbk
      use mackw
      use macttyp
      use macvar
      use macmisc
      implicit none
      integer*4 cmdidx,vidx
      real*8 val
c.....cmdidx brings pointer to command name
      integer*8 argdfp,argp,ktcaloc
      integer*4 i,offset
c
c.......for debug
c       print *,'setdfl>',pname(cmdidx),pname(vidx),val
c.......end debug
      argdfp=idval(cmdidx)
      argp=klist(argdfp)
      if(argp .le. 0) then
        argp=ktcaloc((ilist(1,argdfp-1)-1)*2+1)
        klist(argdfp)=argp
        ilist(1,argp-1)=ilist(1,argdfp-1)-1
        do 1100 i=1,ilist(1,argp-1)
          ilist(1,argp+i*2-1)=ilist(1,argdfp+i+1)
 1100   continue
      endif
        offset=0
        do 2100 i=1,ilist(1,argp-1)
          if(ilist(1,argp+i*2-1) .eq. vidx) then
             offset=i
             go to 2200
          endif
 2100   continue
c
c 2200   write(*,*)'setdfl ',argp,offset,ilist(1,argp-1),
c     $       argdfp,klist(argdfp)
 2200   if(idval(vidx) .eq. VarInt) then
          ilist(1,argp+offset*2)=int(val)
c.......for debug
c       print *,'setdfl<',pname(vidx),offset,ilist(2,argp+offset)
c.......end debug
        else if (idval(vidx) .eq. VarRl) then
          rlist(argp+offset*2)=val
c.......for debug
c       print *,'setdfl<',pname(vidx),offset
c    &                   ,rlist(ilist(2,argp+offset))
c.......end debug
        else
          call errmsg('setdfl','Cannot initialize'//pname(vidx),0,4)
        endif
        return
        end
