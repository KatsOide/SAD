      subroutine lread(idxl,token,slen,ttype,rval,ival)
      use maccbk
      use tfstk , only : forcesf
      implicit none
      include 'inc/CBKMAC.inc'
      integer*4 idxl,idx
      character*(MAXSTR) token,wkstr
      integer*4 slen,ival,ttype
      real*8 rval
c     
      integer*4 hsrchz,mtaloc
      integer*4 llen,mult,blksz,newblk,i
      logical skipch
c     
c     first character must be '='
      if (.not. skipch(EQCHR,token,slen,ttype,rval,ival)) then
         call errmsg('lread',
     &        'missing '''//EQCHR//''' is assumed.',0,0)
      endif
c     next is to be a '('
      if (.not. skipch(LPAR,token,slen,ttype,rval,ival)) then
         call errmsg('lread',
     &        'missing '''//LPAR//''' is assumed.',0,0)
      endif
c     intialize parameters
      blksz=pagesz/4
      idval(idxl)=mtaloc(blksz)
      if(idval(idxl) .eq. 0) then
         call errmsg('lread',' cannot allocate memory',0,0)
         call forcesf()
      end if
      llen=0
c     now read elemnt list
 2000 continue
      if (skipch(COMMA,token,slen,ttype,rval,ival)) go to 2000
      mult=1
 2100 continue
      if(skipch(PLUS,token,slen,ttype,rval,ival))then
         go to 2100
      else if(skipch(MINUS,token,slen,ttype,rval,ival)) then
         mult=-mult
         go to 2100
      else if (ttype .eq. ttypID) then
         idx=hsrchz(token(:slen))
         if ( (idtype(idx) .gt. icMXEL) .and.
     &        (idtype(idx) .ne. icLINE)      ) then
            wkstr(:slen)=token(:slen)
            call errmsg('lread',
     &           'syntax error:'//wkstr(:slen)
     &           //'is not line nor element.',0,16)
         else
            llen=llen+1
            if(llen+1 .gt. blksz ) then
               newblk=mtaloc(2*blksz)
               if (newblk .eq. 0) then
                  call errmsg('lread',
     &                 ' cannot extend working area.',32,0)
                  call forcesf()
               end if
               do i=1,llen-1
                  ilist(1,newblk+i)=ilist(1,idval(idxl)+i)
                  ilist(2,newblk+i)=ilist(2,idval(idxl)+i)
               end do
               call tfreem(idval(idxl),blksz)
c               call freeme(idval(idxl),blksz)
               blksz=2*blksz
               idval(idxl)=newblk
            endif
            ilist(1,idval(idxl)+llen)=mult
            ilist(2,idval(idxl)+llen)=idx
         endif
      else if(ttype .eq. ttypNM) then
         if (ival .eq. 0) then
            call errmsg('lread',
     &           'invalid muliplicity',0,16)
         else
            mult=ival
         endif
         call gettok(token,slen,ttype,rval,ival)
         if (.not. skipch(STAR,token,slen,ttype,rval,ival)) then
            if (ttype .eq. ttypID) then
               call errmsg('lread',
     &              'missing '''//STAR//''' is assumed',0,0)
            else
               call errmsg('lread',
     &              'syntax error:<line elemnt>:: <num>*!<ID>',0,16)
            endif
         endif
         go to 2100
      else if (ttype .eq. ttypDM) then
         if(.not. skipch(RPAR,token,slen,ttype,rval,ival))then
            call errmsg('lread',
     &           'ivalid delimiter',0,0)
         endif
         ilist(1,idval(idxl))=llen
         ilist(2,idval(idxl))=0
         if(llen+4 .gt. blksz) then
            call errmsg('lread','too long error list',0,32)
         endif
         call tfreem(idval(idxl)+llen+1,blksz-llen-1)
c         call freeme(idval(idxl)+llen+1,blksz-llen-1)
         return
      else
         call errmsg('lread',
     &        'syntax error: <num>*!<element name> is expected',0,16)
      endif
      call gettok(token,slen,ttype,rval,ival)
      go to 2000
c
      end
