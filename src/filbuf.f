      subroutine filbuf(fst)
      use tfrbuf
      use maccbk
      use macttyp
      use macfile
      use macmisc
      implicit none
c
cccccc K. Oide 8/30/1999
c
      integer fst
      character*(MAXLLEN) str
      integer*4 ione,LL1, idummy,ii
      integer*4 lene,igetGL,irtc
      integer*8 flmgr
      logical at1st
      save at1st
      data at1st/.true./
cc
  10  continue
      if(at1st) then
        at1st=.false.
        LL1=MAXLLEN
        ione=1
        call JNLPRM(ione,buf,LL1,pbuf)
        if(pbuf .ne. 0) go to 11
      end if
c
c  2/15/2008 K. Oide 
      call readstr(infl,str,irtc)
      if(irtc .ne. 0)then
        if(irtc .lt. 0)then
          go to 999
        endif
        go to 998
      endif
c
c      read(infl,100,end=999,err=998,iostat=ios)str
c 100  format(A)
c
c     6/3/1998 K. Oide
c
      call tmovb(str,buf,MAXLLEN)
c      do i=1,MAXLLEN
c         buf(i)=str(i:i)
c      end do
  11  if (IgetGL('$ECHO$',idummy) .eq. FLAGON
     1     .and. (infl .ne. STDIN
     1     .or. igetGL('$LOG$',idummy) .eq. FLAGON) )
     & then
         write(errfl,'(a)')str(:lene(str))
      endif
      pbuf=1
      fst=pbuf
      return
c
 998  continue
      call errmsg('filbuf', ' file read error',0,0)
 999  continue
cccccccccccc   K. Oide 11/22/1997
cccccccccccc   K. Oide 8/30/1999
c        if (infl .eq. STDIN)then
        if (infl .eq. STDIN .or. itbuf(infl) .eq. 2)then
cccccccccccc   K. Oide end
          fst=ttypEF
          pbuf=MAXLLEN
          buf(pbuf)=' '
        else
          close(infl)
          ii=infl-infl
          infl=int(flmgr(ii))
c        write(*,*)'@ii',ii,infl
c         call errmsg('filbuf',
c    &         'input file is redirected to STDIN',0,0)
          GO TO 10
        endif
      return
      end
