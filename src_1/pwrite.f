      subroutine pwrite(word,latt,twiss,gammab,pos,imon,nmon)
      use tfstk
      use ffs
      use tffitcode
      implicit real*8(a-h,o-z)
      character*(*) word
C*DEC
      character*80  fnam, fnam2
C*DEC
C*HP
C     character*80  fnam
C*HP
      logical lod,lex,rmatq
      integer*8 latt(nlat)
      dimension twiss(nlat,-ndim:ndim,ntwissfun),gammab(nlat),pos(nlat)
      dimension imon(*)
      data io/50/
      if(word.eq.'RMATQ') then
        rmatq=.true.
      else
        rmatq=.false.
      endif
      call getwdl(word)
      if(word.eq.' ') return
      fnam=word
C*DEC
      fnam2=word
C*DEC
      call getwdl(word)
      inquire(unit=io,opened=lod,iostat=ios)
      if(lod) close(io)
      inquire(iostat=ios,file=fnam,exist=lex)
      if(lex) then
C*DEC
        open (io,file=fnam2,iostat=ios,err=200,status='OLD',
C*DEC
C*HP
C       open (io,file=fnam,iostat=ios,err=200,status='OLD',
C*HP
     1        access='SEQUENTIAL',form='UNFORMATTED'
     1        )
      else
C*DEC
        open (io,file=fnam2,iostat=ios,err=200,status='NEW',
C*DEC
C*HP
C       open (io,file=fnam,iostat=ios,err=200,status='NEW',
C*HP
     1        access='SEQUENTIAL',form='UNFORMATTED'
     1        )
      endif
      if(ios .ne. 0) then
C*DEC
        print *,'Cannot open ',fnam2
C*DEC
C*HP
C       print *,'Cannot open ',fnam
C*HP
        return
      endif
      if(.not.rmatq) goto 100
      lq=0
      do 10 l=1,nlat-1
        if(idtype(ilist(2,latt(l))).eq.icquad) then
          t=rlist(idval(ilist(2,latt(l)))+4)
c         .... reject skew quads ....
          if(abs(t-pi/4d0).lt.0.01 .or. abs(t+pi/4d0).lt.0.01) goto 10
          lq=lq+1
        endif
   10 continue
      ia=italoc(2*nmon*lq)
      ib=italoc(lq)
c
      call pwmatq(latt,twiss,gammab,pos,rlist(ia),rlist(ib),imon,nmon,
     1            lq,io)
c
      call tfree(int8(ia))
      call tfree(int8(ib))
  100 continue
  200 return
      end
