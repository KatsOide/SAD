      subroutine mfnst(latt,twiss,isb,ncb,istr,nstr,yplane,iret,lfno)
c  ----- find <ncb> correctors around the target <n> ---------
c        istr(i) is assumed in ascending order
c        ncb     is assumed to be even number
      use ffs
      use tffitcode
      logical yplane,extend,second,coup,pvert,mhogal
      dimension latt(*),twiss(*),isb(*),istr(nstra,4)
      data nbefor/1/
      include 'inc/common.inc'
      save nbefor
c
      extend=.false.
      iret=0
      nrmax=ncb/2
      nlmax=ncb/2
      second=.false.
    9 if(second)then
        is=1
        if=nbefor-1
      else
        is=nbefor
        if=nstr
      endif
      do 10 i=is,if
        j1=istr(i,2)
        j2=istr(mod(i+nstr,nstr)+1,2)
        if(mhogal(istr(j1,1),istr(j2,1),isb(1)))then
          il=i
          ir=mod(i+nstr,nstr)+1
c         print *,istr(istr(il,2),1),istr(istr(ir,2),1),isb(1),' mfnst'
          goto 1
        endif
   10 continue
      if(second)then
        iret=1
        goto 99
      endif
      second=.true.
      goto 9
c ..... search steerings .....
    1 continue
c ..... a gauche
      nl=0
      nx=0
      ny=0
      jl=il
      second=.false.
   20 jl=mod(jl-1+nstr,nstr)+1
c     print *,'jl=',jl,' extend=',extend
      if(nl.eq.nlmax) goto 21
      if(jl.eq.il .and. second) then
        iret=2
        goto 99
      endif
      j1=istr(jl,2)
      if( pvert(latt,istr(j1,1)) .eqv. yplane ) then
        if(extend) then
          ny=ny+1
          if(ny.le.nlmax-1) then
            nl=nl+1
            isb(3+nlmax-nl)=jl
          endif
        else
          nl=nl+1
          isb(3+nlmax-nl)=jl
        endif
c       print *,' pvert=',pvert(latt,istr(jl)),' ny=',ny,' nl=',nl
      elseif(extend) then
        nx=nx+1
        if(nx.lt.2) then
          nl=nl+1
          isb(3+nlmax-nl)=jl
        endif
c       print *,' pvert=',pvert(latt,istr(jl)),' nx=',nx,' nl=',nl
      endif
      jl=jl-1
      second=.true.
      goto 20
   21 continue
c ..... a droit
      nr=0
      nx=0
      ny=0
      jr=ir
      second=.false.
   30 jr=mod(jr-1+nstr,nstr)+1
      if(nr.eq.nrmax) goto 31
      if(jr.eq.ir .and. second) then
        iret=3
        goto 99
      endif
      j1=istr(jr,2)
      if( pvert(latt,istr(j1,1)) .eqv. yplane ) then
        if(extend) then
          ny=ny+1
          if(ny.le.nrmax-1) then
            nr=nr+1
            isb(2+nlmax+nr)=jr
          endif
        else
          nr=nr+1
          isb(2+nlmax+nr)=jr
        endif
      elseif(extend) then
        nx=nx+1
        if(nx.lt.2) then
          nr=nr+1
          isb(2+nrmax+nr)=jr
        endif
      endif
      jr=jr+1
      second=.true.
      goto 30
   31 continue
      if(.not.extend) then
        j1=istr(isb(3),2)
        j2=istr(isb(2+nlmax+nrmax),2)
        call pfcoup(latt,twiss,ndim,
     1              istr(j1,1),istr(j2,1)+1,coup)
        if(coup) then
          nlmax=nlmax+1
          nrmax=nrmax+1
          extend=.true.
          goto 1
        endif
      endif
      nbefor=isb(2+nlmax)
      isb(2)=nlmax+nrmax
      return
   99 call permes(' ','Correction dipole not found.',' ',lfno)
      return
      end
