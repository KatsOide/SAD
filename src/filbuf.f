      subroutine filbuf(fst)
      use tfrbuf
      use maccbk
      use macttyp
      use macfile
      use macmisc
      use ffsfile
      implicit none
c
cccccc K. Oide 8/30/1999
c
      integer fst
      character*(MAXLLEN) str
      integer*4 LL1,idummy
      integer*4 lene,igetGL,irtc
      logical at1st
      save at1st
      data at1st/.true./
cc
  10  continue
      if(at1st) then
        at1st=.false.
        LL1=MAXLLEN
        call JNLPRM(buf,LL1,pbuf)
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
  11  if (IgetGL('$ECHO$') .eq. FLAGON
     1     .and. (infl .ne. STDIN
     1     .or. IgetGL('$LOG$') .eq. FLAGON) )
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
      if (infl .eq. STDIN .or. itbuf(infl) .eq. modewrite)then
        fst=ttypEF
        pbuf=MAXLLEN
        buf(pbuf)=' '
      else
        call trbclose(infl)
        if(lfnp .gt. lfnbase)then
c          write(*,*)'filbu-EOF ',infl,lfnp,lfnstk(lfnp)
          lfnp=lfnp-1
          infl=lfnstk(lfnp)
          call trbassign(infl)
        elseif(lfnbase .eq. 0)then
c          write(*,*)'filbu-EOF1 ',infl
          infl=STDIN
c         call errmsg('filbuf',
c     &         'input file is redirected to STDIN',0,0)
          GO TO 10
        else
          fst=ttypEF
          pbuf=MAXLLEN
          buf(pbuf)=' '
        endif
      endif
      return
      end
