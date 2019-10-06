      subroutine gettok(token,slen,type,rval,ival)
      use ISO_C_BINDING
      use macttyp
      use macfile
      use macmisc
      implicit none
      character*(*) token
      integer*4 slen,type,ival
      real*8 rval,atof
c
      integer*4 iost,lens,ik,ipak
      character*(*) , parameter ::TAB=C_HORIZONTAL_TAB,EOL=C_NEW_LINE,
     $     CR=c_carriage_return,
     $     SPACES=' '//TAB//EOL//CR,PERIOD='.',
     $     DELIMS=SPACES//'!"#%''()=!^-{`@}!*+;:/?><,~'
      integer*4, parameter::LLEN=MAXLLEN
c
      logical issign,isdgt,isalph
      integer fltnum
      character chr,quot,curc,nc
c     for ungettok
      logical*4 inbuf
      data inbuf /.false./
      character*(MAXLLEN) btoken
      integer*4 bslen,btype,bival
      real*8 brval
      save inbuf,btoken,bslen,btype,bival,brval

c macro functions
      issign(chr)=(chr .eq. '+') .or. (chr .eq. '-')
      isdgt(chr)=(LGE(chr,'0') .and. LLE(chr,'9'))
      isalph(chr)= ( (LGE(chr,'A') .and. LLE(chr,'Z'))  .or.
     $      (LGE(chr,'a') .and. LLE(chr,'z')) )
c     begin initialize for preventing compiler warning
      fltnum=0
      quot=' '
c     end   initialize for preventing compiler warning
      if (inbuf) then
        token=btoken
        slen=bslen
        type=btype
        rval=brval
        ival=bival
        inbuf=.false.
        return
      endif
 1000 continue
      call mygetchr(curc,iost,*999)
c skip spaces
      continue
      do while(index(SPACES//char(0),curc) .ne. 0)
        call mygetchr(curc,iost,*999)
      enddo
c skip comment
      if (curc .eq. '!' ) then
        call myfflush
        go to 1000
      endif
c
      slen=1
      token(1:1)=curc
      if (isalph(curc) ) then
        type=ttypID
      else if (isdgt(curc)) then
        type=ttypNM
        fltnum=0
      else if ( issign(curc) ) then
         type=ttypDM
      else if (curc .eq. '.') then
         call mycurchr(nc,iost)
         if ( isdgt(nc)) then
            type=ttypNM
            fltnum=1
         else
            type=ttypDM
         endif
c     for quoted string. april 27,'87 by N.Y.
      else if (curc .eq. '''') then
        type=ttypST
        quot=''''
        slen=0
      else if (token(1:1) .eq. '"') then
        type=ttypST
        quot='"'
        slen=0
      else
        type=ttypDM
      endif
c main routine of gettok
      call mycurchr(curc,iost)
      if (type .eq. ttypID)then
 100    continue
        if( index(char(0)//DELIMS, curc) .eq. 0) then
           call mygetchr(curc,iost,*900)
           slen=slen+1
           token(slen:slen)=curc
           call mycurchr(curc,iost)
           go to 100
        endif
        go to 900
      else if (type .eq. ttypNM) then
 200     continue
         if( isdgt(curc) ) then
           call mygetchr(curc,iost,*900)
           slen=slen+1
           if(slen .gt. LLEN)then
             go to 900
           endif
           token(slen:slen)=curc
           call mycurchr(curc,iost)
           go to 200
         endif
         if (fltnum .lt. 2)  then
           if ((fltnum .eq. 0) .and. (curc .eq. '.')) then
             fltnum=1
             call mygetchr(curc,iost,*900)
             slen=slen+1
             token(slen:slen)=curc
             call mycurchr(curc,iost)
             go to 200
           else if ( index('dDeE',curc) .ne. 0) then
             fltnum=2
             call mygetchr(curc,iost,*900)
             slen=slen+1
             token(slen:slen)=curc
             call mycurchr(curc,iost)
             if ( token(slen:slen) .ne. '.'
     &            .and. (issign(curc) )) then
                call mygetchr(curc,iost,*900)
                slen=slen+1
                token(slen:slen)=curc
                call mycurchr(curc,iost)
             endif
             go to 200
           endif
         endif
         if (fltnum .ne. 0) then
           lens=slen
           rval=atof(token,lens)
c          ival=int(rval)
         else
           ik=slen
           ival=ipak(token,ik)
           ik=slen
           rval=atof(token,ik)
         endif
      else if (type .eq. ttypST) then
 300  continue
        if(curc .eq. quot) then
           call mygetchr(curc,iost,*910)
           call mycurchr(curc,iost)
           if(curc .eq. quot) then
              call mygetchr(curc,iost,*910)
              slen=slen+1
              if(slen .gt. LLEN)then
                go to 900
              endif
              token(slen:slen)=quot
              call mycurchr(curc,iost)
              go to 300
           end if
           go to 900
        else 
           call mygetchr(curc,iost,*910)
           slen=slen+1
           if(slen .gt. LLEN)then
             go to 900
           endif
           token(slen:slen)=curc
           call mycurchr(curc,iost)
           go to 300
        end if
      end if
c   
 900  continue
      if (slen .gt. LLEN)
     &   call errmsg('gettok',
     &               'too long identifier: '//token(:slen),0,16)
c      write(*,*) "new_gettok:900: ",token(:slen)
      return
c
 910  continue
      call errmsg('gettok','unmatched quote',0,0)
c      write(*,*) "new_gettok:910: ",token(:slen)
      return
c
  999 continue
      call errmsg('gettok','EOF has been read.',0,0)
      token=' '
      type=ttypEF
      slen=0
c      write(*,*) "new_gettok:999: ",token(:slen)
      return
c     entry  ungettok(token,slen,type,rval,ival)
c      btoken=token
c      bslen=slen
c      btype=type
c      brval=rval
c      bival=ival
c      inbuf=.true.

c      return
c
      end
c
c****************** this routine uses an alternate return
      subroutine mygetchr(chr,iost,*)
      use ISO_C_BINDING
      use macttyp
      use macfile
      use macmisc
      implicit none
      character chr
      integer*4 iost
c
      integer LLEN
      parameter(LLEN=MAXLLEN)
c
      if (pbuf .le. 0) then
         call filbuf(iost)
         if (iost .eq. ttypEF) then
            chr=char(255)
            return 1
         endif
      end if
      if( pbuf .gt. LLEN) then
         chr=C_NEW_LINE
         pbuf=0
      else
         chr=buf(pbuf)
         pbuf=pbuf+1
      end if
      iost=0
      return
      end
c
      subroutine mycurchr(chr,iost)
      use ISO_C_BINDING
      use macttyp
      use macfile
      use macmisc
      implicit none
      character chr
      integer*4 iost
c
      integer LLEN
      parameter(LLEN=MAXLLEN)
c
      if(pbuf .le. 0) then
            chr=C_NULL_CHAR
            iost=0
            return
      end if
      if(pbuf .gt. LLEN) then
         chr=C_NEW_LINE
         iost=0
      else
         chr=buf(pbuf)
         iost=0
      end if
      return
c
      entry  mynxtchr(chr,iost)
      if(pbuf .le. 0) then
         chr=C_NULL_CHAR
         iost=0
      else if(pbuf .ge. LLEN) then
         chr=C_NEW_LINE
         iost=0
      else
         chr=buf(pbuf+1)
         iost=0
      end if
      return
c
      entry myfflush
      pbuf=0
      return
c
      end
