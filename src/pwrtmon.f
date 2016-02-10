      subroutine pwrtmon(word,latt,twiss,mult,imon,nmon,lfno)
      use ffs
      use tffitcode
      parameter (lio=50,lstdout=6)
      logical lod
      character*(*) word
      character*13 outc,autofg
      character*16 en
      dimension latt(2,nlat),twiss(nlat,-ndim:ndim,ntwissfun),mult(*),
     &          imon(nmona,4)
      dimension rms(2),amin(2),amax(2)
      include 'inc/common.inc'
c
c     call getwdlK(word)
      call getwdl(word)
      call texpfn(word)
      if(word.eq.' ') then
        call getwdl(word)
        io=lfno
      else
        io=lio
        inquire(unit=io,opened=lod,iostat=ios)
        if(lod) close(io)
        open (io,file=word,iostat=ios,status='UNKNOWN',err=99)
        if(ios .ne. 0) then
          call permes(' Cannot open',' ',word,lfno)
          return
        endif
      endif
      rms(1)=0d0
      rms(2)=0d0
      amin(1)=twiss(imon(imon(1,2),1),0,15)
      amin(2)=twiss(imon(imon(1,2),1),0,17)
      amax(1)=amin(1)
      amax(2)=amin(2)
      do 10 i=1,nmon
        call elnameK(latt,imon(imon(i,2),1),mult,en)
        call mbufw(en,.false.,io)
        outc=autofg(twiss(imon(imon(i,2),1),0,15),'S13.10')
        call mbufw(outc,.false.,io)
        outc=autofg(twiss(imon(imon(i,2),1),0,17),'S13.10')
        call mbufw(outc,.false.,io)
        rms(1)=rms(1)+twiss(imon(imon(i,2),1),0,15)**2
        rms(2)=rms(2)+twiss(imon(imon(i,2),1),0,17)**2
        amin(1)=min(amin(1),twiss(imon(imon(i,2),1),0,15))
        amin(2)=min(amin(2),twiss(imon(imon(i,2),1),0,17))
        amax(1)=max(amax(1),twiss(imon(imon(i,2),1),0,15))
        amax(2)=max(amax(2),twiss(imon(imon(i,2),1),0,17))
 10   continue
      rms(1)=sqrt(rms(1)/dble(max(1,nmon)))
      rms(2)=sqrt(rms(2)/dble(max(1,nmon)))
      call mbufw(' ',.true.,io)
c
      call mbufg('! rms',.false.,io) 
      outc=autofg(rms(1),'S13.10')
      call mbufw(outc,.false.,io)
      outc=autofg(rms(2),'S13.10')
      call mbufw(outc,.false.,io)
c
      outc=autofg(amin(1),'S13.10')
      call mbufw('min '//outc,.false.,io)
      outc=autofg(amin(2),'S13.10')
      call mbufw(outc,.false.,io)
c
      outc=autofg(amax(1),'S13.10')
      call mbufw('max '//outc,.false.,io)
      outc=autofg(amax(2),'S13.10')
      call mbufw(outc,.false.,io)
c
      call mbufw(' ',.true.,io)
      if(io.ne.lstdout) close(io)
 99   call getwdl(word)
      return
      end
