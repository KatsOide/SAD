      subroutine pwrtstr(word,latt,twiss,mult,istr,nstr,lfno)
      use tfstk
      use ffs
      use tffitcode
      implicit real*8(a-h,o-z)
      parameter (lio=50,stndout=6)
      logical lod,pvert
      character*(*) word
      character*13 outc,autofg
      character*16 en
      integer*8 latt(nlat)
      dimension twiss(nlat,-ndim:ndim,ntwissfun),mult(*),istr(nstra,4)
      dimension rms(2),amin(2),amax(2),ave(2)
      include 'inc/common.inc'
      include 'inc/coroper.inc'
c
      io=lio
c     call getwdlK(word)
      call getwdl(word)
      call texpfn(word)
      if(word.eq.' ') then
        call getwdl(word)
        io=lfno
      else
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
      amin(1)=1d20
      amin(2)=1d20
      amax(1)=-amin(1)
      amax(2)=-amin(2)
      ave(1)=0d0
      ave(2)=0d0
      nx=0
      ny=0
      if(istope.ne.0) then
        do 10 i=1,nstope
          j=ilist(mod(i-1,2)+1,istope+(i-1)/2)
          call elnameK(istr(j,1),en)
          call mbufw(en,.false.,io)
          outc=autofg(rlist(latt(istr(j,1))+11),'S13.10')
          call mbufw(outc,.false.,io)
          if(pvert(latt,istr(j,1))) then
            rms(2)=rms(2)+rlist(latt(istr(j,1))+11)**2
            amin(2)=min(amin(2),rlist(latt(istr(j,1))+11))
            amax(2)=max(amax(2),rlist(latt(istr(j,1))+11))
            ave(2)=ave(2)+rlist(latt(istr(j,1))+11)
            ny=ny+1
          else
            rms(1)=rms(1)+rlist(latt(istr(j,1))+11)**2
            amin(1)=min(amin(1),rlist(latt(istr(j,1))+11))
            amax(1)=max(amax(1),rlist(latt(istr(j,1))+11))
            ave(1)=ave(1)+rlist(latt(istr(j,1))+11)
            nx=nx+1
          endif
 10     continue
      else
        do 11 i=1,nstr
          j=istr(i,2)
          call elnameK(istr(j,1),en)
          call mbufw(en,.false.,io)
          outc=autofg(rlist(latt(istr(j,1))+11),'S13.10')
          call mbufw(outc,.false.,io)
          if(pvert(latt,istr(j,1))) then
            rms(2)=rms(2)+rlist(latt(istr(j,1))+11)**2
            amin(2)=min(amin(2),rlist(latt(istr(j,1))+11))
            amax(2)=max(amax(2),rlist(latt(istr(j,1))+11))
            ave(2)=ave(2)+rlist(latt(istr(j,1))+11)
            ny=ny+1
          else
            rms(1)=rms(1)+rlist(latt(istr(j,1))+11)**2
            amin(1)=min(amin(1),rlist(latt(istr(j,1))+11))
            amax(1)=max(amax(1),rlist(latt(istr(j,1))+11))
            ave(1)=ave(1)+rlist(latt(istr(j,1))+11)
            nx=nx+1
          endif
 11     continue
      endif
      call mbufw(' ',.true.,io)
      rms(1)=sqrt(rms(1)/dble(max(1,nx)))
      rms(2)=sqrt(rms(2)/dble(max(1,ny)))
      ave(1)=ave(1)/dble(max(1,nx))
      ave(2)=ave(2)/dble(max(1,ny))
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
      outc=autofg(ave(1),'S13.10')
      call mbufw('mean '//outc,.false.,io)
      outc=autofg(ave(2),'S13.10')
      call mbufw(outc,.false.,io)
c
      call mbufw(' ',.true.,io)
      if(io.ne.stndout) close(io)
 99   call getwdl(word)
      return
      end
