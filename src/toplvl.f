      module trackbypass
      implicit none
      logical*4 :: bypasstrack=.false.
      integer*4 :: lattuse=0, lattredef=0
      end module

      subroutine toplvl
      use trackbypass, only: bypasstrack
      use tfrbuf
      use maccbk
      implicit none
      include 'inc/MACCODE.inc'
      include 'inc/MACKW.inc'
      include 'inc/MACTTYP.inc'
      include 'inc/MACMISC.inc'
      include 'inc/MACFILE.inc'
      integer*4 jslen,jttype,jival
      integer*4 slen,ttype,ival,hsrch,hsrchz,status
      integer*4 argp,index,IgetGL,idummy
      real*8  jrval,rval
      character*(MAXSTR) jtoken,token
c     
 1000 continue
      if (IgetGL('$CTIME$',idummy) .eq. FLAGON) call cputix
      call gettok(token,slen,ttype,rval,ival)
c     for debug
c       print *,token(:slen),slen,ttype
c     end debug
c     
 1100 continue
      if (ttype .eq. ttypEF) go to 9000
      if(ttype .eq. ttypID) then
c..........System defined name.
         index=hsrch(token(:slen))
         if (idtype(index) .eq. icDEF) then
            if (idval(index) .lt. icMXEL) then
              call tfinitstk
              call doelem(idval(index))
            else if((idval(index) .ge. icLINE)
     &              .and. (idval(index) .lt. icRSVD))  then
              call tfinitstk
              call doline(idval(index))
            else
              call errmsg('toplvl','system bug: check initbl',0,0)
            endif
            go to 1000
         else if (idtype(index) .eq. icACT) then
            if(bypasstrack)then
              return
            endif
            call tfinitstk
            call doACT(index)
            go to 1000
         else if (idtype(index) .eq. icRSVD) then
           call tfinitstk
           call funcall1(idval(index),argp)
           go to 1000
         endif
c..........User defined variable.
c     print *,token(:slen)
         index=hsrchz(token(:slen))
         if ((idtype(index) .lt. icMXEL) .and.
     &        (idtype(index) .ne. icNULL)) then
           call tfinitstk
           call doelm2(idtype(index),pname(index),slen,ttype)
         else if ((idtype(index) .eq. icLINE)) then
           call tfinitstk
           call dolin2(index,slen,ttype)
         else
           call dAssgn(token,slen,status)
           if (status .ne. 0) call errmsg('toplvl',
     &          'Unsupported function '//token(:slen)//'!',
     &          0,16)
c     for debug
c      print *,token(:slen),
c     &     index,pname(index),idval(index),idtype(index)
c     end debug
         endif
      else if(token(:slen) .eq. ';') then
         go to 1000
      else if(token(:slen) .eq. RCURL) then
         go to 1000
       elseif(token(:slen) .eq. char(13))then
         go to 1000
      else
         call errmsg('main'
     &        ,'syntax error: invalid input for toplvl '//
     &        '"'//token(:slen)//'"'
     &        ,0,0)
      endif
      go to 1000
c     
 9000 continue
cccccccccccc   K. Oide 11/22/1997
cccccccccccc   K. Oide 9/2/1999
      if(itbuf(infl) .ne. 0)then
        return
      endif
cccccccccccc   K. Oide end
      print *," SAD1 reads EOF."
      call errmsg("main","Stop execution.(READ EOF)" ,0,0)
c     
      stop
c.....Entry point for error handling.
      entry bigjmp(jtoken,jslen,jttype,jrval,jival)
c........for debug
c     print *,'Big Jump !!!'
c     
      token=jtoken
      slen=jslen
      ttype=jttype
      rval=jrval
      ival=jival
      go to 1100
      end
