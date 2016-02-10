      integer function rslvin(idx)
      use maccbk
      implicit real*8 (a-h,o-z)
      integer idx
      include 'inc/MACCODE.inc'
      integer llen,lptnew,lptold,direct,hsrchz
      integer STACSZ
      parameter (STACSZ=400)
      integer istack(STACSZ),pstack
c.....statment functions
      icar(i)=ilist(1,i)
      icdr(i)=ilist(2,i)
c
      lptold=idval(idx)
      llen=ilist(1,lptold)
      lptnew=mcfallo(llen+1)
      rslvin=lptnew
      ilist(1,lptnew)=llen
      ilist(2,lptnew)=0
c
      direct=1
      pstack=0
      i0=0
 1100 continue
        i0=i0+1
        if (i0 .gt. ilist(1,lptold)) go to 1400
c
        i= ((1-direct)*(ilist(1,lptold)+1)+ 2*direct*i0)/2
        if (idtype(ilist(2,lptold+i)) .lt. icMXEL) then
          if (idtype(icdr(lptold+i)) .eq. icNULL) then
            call errmsg('rslvin',
     &           pname(icdr(lptold+i))//' is not defined yet',0,0)
            call errmsg('rslvin',
     &           'warnig:unable to resolve inverted line',0,0)
          endif
          ilist(1,lptnew+i0)=ABS(icar(lptold+i))
          if(direct .ge. 0) then
            ilist(2,lptnew+i0)=ilist(2,lptold+i)
          else
            ilist(2,lptnew+i0)=lrflct(ilist(2,lptold+i))
          endif
        else if (idtype(icdr(lptold+i)) .eq. icLINE) then
          if(idx .eq. icdr(lptold+i)) then
            call errmsg('expnln',
     &           'line def. contains itself. unable to resolve',0,16)
          endif
c for debug
c         call ptrace('nested definition  '//
c    &                 pname(icdr(lptold+i))//'!'
c    &    ,1)
c end debug
          pstack= pstack+4
          if (pstack .gt. STACSZ) then
            call errmsg('rslvin',
     &           'stack over flow. too deep definition of line'
     &           ,0,16)
          endif
          istack(pstack)=lptold
          istack(pstack-1)=lptnew
          istack(pstack-2)=i0
          istack(pstack-3)=direct
c.........
          if (direct*ilist(1,lptold+i) .ge. 0) then
            ilist(1,lptnew+i0)=abs(icar(lptold+i))
            ilist(2,lptnew+i0)=hsrchz('+'//pname(icdr(lptold+i)))
            direct=1
          else
            ilist(1,lptnew+i0)=abs(icar(lptold+i))
            ilist(2,lptnew+i0)=hsrchz('-'//pname(icdr(lptold+i)))
            direct=-1
          endif
          idtype(icdr(lptnew+i0))=idtype(icdr(lptold+i))
          idval(icdr(lptnew+i0))=mcfallo(icar(idval(icdr(lptold+i)))+1)
          ilist(1,idval(icdr(lptnew+i0)))=
     &                            icar(idval(icdr(lptold+i)))
c.........
          lptold=idval(icdr(lptold+i))
          lptnew=idval(icdr(lptnew+i0))
          i0=0
          go to 1100
        else if (idtype(ilist(2,idx+i)) .gt. icMXEL) then
            call errmsg('rslvin',
     &           pname(ilist(2,idx+i))//' is not a element ',0,16)
        endif
      go to 1100
c
 1400 continue
      if ( pstack .gt. 0) then
        lptold=istack(pstack)
        lptnew=istack(pstack-1)
        i0=istack(pstack-2)
        direct=istack(pstack-3)
        pstack=pstack-4
c.......for debug
c         call ptrace('end of nest.  '//
c    &                 pname(icdr(lptnew+i0))//'!'
c    &    ,-1)
c.......end debug
        go to 1100
      endif
      return
      end
