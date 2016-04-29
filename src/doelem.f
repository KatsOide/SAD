      subroutine doelem(elmcd)
      use maccbk
      implicit none
c
      include 'inc/MACFILE.inc'
      include 'inc/MACMISC.inc'
      include 'inc/MACTTYP.inc'
      include 'inc/MACCODE.inc'
      include 'inc/MACKW.inc'
      integer*4 elmcd
      character*(MAXSTR) token,wtoken*(*)
      integer slen,ival,ttype,slen2,ttype2,idx
      integer*4 hsrchz,itcaloc
      real*8 rval
      logical skipch
c
      call gettok(token,slen,ttype,rval,ival)
      if (.not. skipch(LCURL,token,slen,ttype,rval,ival) ) then
        if (ttype .eq. ttypID) then
c         call errmsg('doelem',
c    &              'missing '''//LCURL//''' is assumed',0,0)
        else
          call errmsg('doelem',
     &                ' syntax error:'//LCURL//' is expected.',0,16)
        endif
      endif
      goto 2000
C
c strat element definition
      entry doelm2(elmcd,wtoken,slen2,ttype2)
      slen=slen2
      ttype=ttype2
      token=wtoken(:slen)
 2000 if (skipch(COMMA,token,slen,ttype,rval,ival) ) go to 2000
c for debug
c      print *,'read element name ',ttype,ttypID,token(:slen)
c end debug
      if((token(:slen) .eq. RCURL) .or.
     &   (token(:slen) .eq. SEMIC))    then
         return
      else if (ttype .eq. ttypID) then
        idx =hsrchz(token(:slen))
        if (idtype(idx) .eq. icNULL) then
          idtype(idx)=elmcd
          idval(idx)=itcaloc(kytbl(kwMAX,elmcd)+1)
          call rdkwdl(idx)
c for debug
c      print *,'doelem:',idtype(idx),elmcd,pname(kytbl(0,elmcd))
c    &         ,kytbl(kwMAX,elmcd)
c end debug
        else if(idtype(idx) .eq. elmcd) then
c         call errmsg('doelm',
c    &              'element: '//token(:slen)//' is redefined',0,0)
          call rdkwdl(idx)
        else if(idtype(idx) .lt. icMXEL) then
          call errmsg('doelem',
     &      'syntax error: '//token(:slen)//
     $         ' is already used as '//
     &      pname(kytbl(0,idtype(idx)))//'.'
     &        ,0,16)
        else if (idtype(idx) .gt. icLINE) then
          call errmsg('doelem',
     &        token(:slen)//'is used as a element name.'
     &        ,0,0)
          call errmsg('doelem',
     &        'syntax error:you can not use reserved words as a name.'
     &        ,0,16)
        else
          call errmsg('doelem',
     &      'syntax error: '//token(:slen)//' is already used as'//
     &      ' anohter element name.'
     &        ,0,16)
        endif
      else
        call errmsg('doelm',
     &       'syntax error:doelm{<line name>=(line definition)}',0,16)
      endif
c
      call gettok(token,slen,ttype,rval,ival)
      go to 2000
c
      end
