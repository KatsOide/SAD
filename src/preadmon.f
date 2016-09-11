      subroutine preadmon(word,latt,twiss,mult,imon,nmon,lfno)
      use ffs
      use tffitcode
      use tfcsi, only:cssetp
      logical lod,lex,stat,name
      character*(*) word
      character*80 outc,ename,en
      character*255 tfconvstr
C*DEC
      character*255 word2
C*DEC
      integer*8 latt(nlat)
      dimension twiss(nlat,-ndim:ndim,ntwissfun),mult(*),
     &          imon(nmona,4)
      external pack
      include 'inc/common.inc'
      data in/50/
c
c     call getwdlK(word)
c     call getwdl(word)
      itype=itfpeeko(ia,x,next)
      lfni1=x+.5d0
      if(itype.eq.1) then
        call cssetp(next)
        write(word,'(''ftn'',i2.2)') lfni1
      elseif(itype.eq.101) then
        call cssetp(next)
        word=tfconvstr(101,ia,x,nc,'*')
        if(word.eq.' ') then
          call permes('?Missing filename for MONREAD','.',' ',lfno)
          call getwdl(word)
          return
        endif
        call texpfn(word)
      endif
      inquire(unit=in,opened=lod,iostat=ios)
      if(lod) close(in)
C*DEC
      word2 = word
      inquire(iostat=ios,file=word2,exist=lex)
C*DEC
C*HP
C     inquire(iostat=ios,file=word,exist=lex)
C*HP
      if(.not.lex) then
        call permes('?File',word,' not found.',lfno)
        call getwdl(word)
        return
      endif
      open (in,file=word,iostat=ios,status='OLD')
      if(ios .ne. 0) then
        call permes('?Cannot open',word,'.',lfno)
        call getwdl(word)
        return
      endif
c
c     .... reset current specification of BPMs
      do 10 i=1,nmona
        imon(imon(i,2),3)=1
 10   continue
      stat=.true.
      name=.true.
      jm=0
 1    continue
      call preabuf(outc,in,stat)
      if(stat) goto 99
      if(name) then
        ename=outc
      else
        read(outc,*,err=98) va1
        call preabuf(outc,in,stat)
        if(stat) goto 99
        read(outc,*,err=98) va2
      endif
      name=.not.name
      if(.not.name) goto 1
      ls=jm
      do 20 i=ls+1,nmona
        call elnameK(imon(imon(i,2),1),en)
        if(ename.eq.en) then
          jm=i
          imon(imon(i,2),3)=0
          goto 22
        endif
 20   continue
      do 21 i=ls,1,-1
        call elnameK(imon(imon(i,2),1),en)
        if(ename.eq.en) then
          jm=i
          imon(imon(i,2),3)=0
          goto 22
        endif
 21   continue
      call permes(ename,' is not a BPM.',' ',lfno)
      goto 1
c
 22   continue
      twiss(imon(imon(jm,2),1),0,15)=va1
      twiss(imon(imon(jm,2),1),0,17)=va2
c     print *,ename(1:lene(ename)),va1,va2
      goto 1
 98   call permes('??',outc,' is not a number.',lfno)
 99   call pack(imon(1,2),imon(1,3),nmon,nmona)
c     errval(1)=0d0
c     errval(2)=0d0
c     errval(3)=0d0
c     errval(4)=0d0
      write(lfno,'(2(A,I4))')'  BPM available =',nmon,' in ',nmona
      nmonact=nmon
      call getwdl(word)
      close(in)
      return
      end
