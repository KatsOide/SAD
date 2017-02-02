      subroutine tdrwdt(kv,y,jp,ymin,ymax,ls,le,frs,fre,ale,np,
     $     only,monly,patt,
     1     latt,twiss,gammab,idp,pos,imon,emon,nmon)
      use tfstk
      use ffs
      use tffitcode
      implicit real*8 (a-h,o-z)
      logical*4 sol,over,fin,monel,only,monly,tmatch
      character*(*) patt
      real*8 y(2,np),twiss(nlat,-ndim:ndim,ntwissfun),pos(nlat),
     $     gammab(nlat),emon(nmona,4),twisss(ntwissfun),ymin,ymax,
     $     frs,fre,hbs
      integer*8 latt(nlat),ip
      integer*4 kv,ls,le,imon(nmona,4)
      real*8 beam(42),sv(5),vsave(100),ftwiss(ntwissfun)
      include 'inc/common.inc'
      theta=0d0
      altotal=pos(nlat)-pos(1)
      alen=pos(le)-pos(ls)
      if(alen .le. 0.d0)then
        alen=alen+altotal
      endif
      jp=0
      sol=.false.
      istrt=ls
      if(ls. gt. le)then
        istop=nlat-1
      else
        istop=le
      endif
      pa=0.d0
      fin=.false.
 111  LOOP_I: do i=istrt,istop-1
        if(only .and. .not. tmatch(pname(ilist(2,latt(i))),patt))then
        elseif(monly .and. .not. monel(imon,nmon,i,m1))then
        elseif(kv .gt. 32 .and. kv .le. 34) then
          if(monel(imon,nmon,i,m1)) then
            jp=jp+1
            y(1,jp)=pos(i)+pa
            ll=imon(m1,2)
            nq=imon(ll,4)
            y(2,jp)=twiss(i,idp,2*kv-51)-twiss(i,ndim,2*kv-51)
     1           -rlist(latt(nq)+kv-28)-emon(ll,kv-32)
          endif
        else
          id=idtype(ilist(2,latt(i)))
          if(id .ge. icMARK)then
            cycle LOOP_I
          endif
          if(i .eq. istrt .and. frs .ne. 0.d0)then
            f=1.d0-frs
            dpa=(pos(i+1)-pos(i))*f
          elseif(i .eq. istop-1 .and. fre .ne. 0.d0)then
            f=fre
            dpa=0.d0
          else
            f=1.d0
            dpa=0.d0
          endif
          if(kv .le. 27 .or. kv .eq. 35)then
            jp=jp+1
            y(1,jp)=pos(i)+pa+dpa
            y(2,jp)=drwfun(kv,i,twiss,gammab,idp,amass,nlat,ndim)
          elseif(kv .le. 32) then
            jp=jp+1
            y(1,jp)=pos(i)+pa+dpa
            if(trsize)then
              call tftrb(latt,i,sv)
              beam(1)=sv(3)
              beam(4)=sv(5)
              beam(6)=sv(4)
            else
              call tfbeam(i,theta,beam)
            endif
            if(kv.ne.32) then
              ia=(kv-27)*(kv-26)/2
              y(2,jp)=sqrt(beam(ia))
            else
              y(2,jp)=atan2(-2d0*beam(4),beam(1)-beam(6)) /2d0
            endif
          endif
          ip=latt(i)
          if(id .le. icdodeca .or. id .eq. icmult
     $         .or. id .eq. iccavi .or. id .eq. icTCAV)then
            al=rlist(ip+kytbl(kwL,id))
            ndiv=int(al*f/ale)
            if(ndiv .gt. 0)then
              nkey=kytbl(kwmax,id)-1
              call tmov(rlist(ip+1),vsave,nkey)
              gbs=gammab(i+1)
              hs=p2h(gbs)
c              hs=gbs*sqrt(1.d0+1.d0/gbs**2)
              ha=p2h(gammab(i))
c              ha=gammab(i)*sqrt(1.d0+1.d0/gammab(i)**2)
              dh=hs-ha
              do k=1,ntwissfun
                twisss(k)=twiss(i+1,idp,k)
              enddo
              if(kv .gt. 27 .and. kv .le. 32 .and. trsize)then
                call tracke(latt,i,sv,np,'STANDBY',' ',0)
                call tracke(latt,i,sv,np,'TRACK',' ',0)
                call tracke(latt,i,sv,np,'SAVE',' ',0)
              endif
              do n=1,ndiv
                rn=dble(n)/(ndiv+1)
                call qtwissfrac(ftwiss,i,rn*f,over)
                do ii=1,ntwissfun
                  twiss(i+1,idp,ii)=ftwiss(ii)
                enddo
                jp=jp+1
                y(1,jp)=pos(i)+pa+dpa+al*rn*f
                if(kv .le. 27 .or. kv .eq. 35) then
                  gammab(i+1)=h2p(rn*f*dh+ha)
c                  gammab(i+1)=sqrt((rn*f*dh+ha)**2-1.d0)
                  y(2,jp)=drwfun(kv,i+1,twiss,gammab,
     $                 idp,amass,nlat,ndim)
                  gammab(i+1)=gbs
                elseif(kv.le.32)then
                  if(trsize)then
                    call tpara(ilist(2,latt(i)))
                    call tracke(latt,i+1,sv,np,'TRACK',' ',0)
                    beam(1)=sv(3)
                    beam(4)=sv(5)
                    beam(6)=sv(4)
                    if(ilist(2,latt(i)) .ne. 0)then
                      call tfree(ilist(2,latt(i)))
                      ilist(2,latt(i))=0
                    endif
                  else
                    hbs=rn*f*dh+ha
                    gammab(i+1)=h2p(hbs)
c                    gammab(i+1)=hbs*sqrt(1.d0-1.d0/hbs**2)
                    call tfbeam(i+1,theta,beam)
                    gammab(i+1)=gbs
                  endif
                  if(kv.ne.32) then
                    ia=(kv-27)*(kv-26)/2
                    y(2,jp)=sqrt(beam(ia))
                  else
                    y(2,jp)=atan2(-2d0*beam(4),beam(1)-beam(6)) /2d0
                  endif
                endif
              enddo
              if(kv .gt. 27 .and. kv .le. 32 .and. trsize)then
                call tracke(latt,i+1,sv,np,'RESET',' ',0)
              endif
              call tmov(vsave,rlist(ip+1),nkey)
              gammab(i+1)=gbs
              do k=1,ntwissfun
                twiss(i+1,idp,k)=twisss(k)
              enddo
            endif
          endif
        endif
c        write(*,*)'tdrwdt ',i,id,jp
      enddo LOOP_I
      i=istop
      if(fre .ne. 0.d0)then
        if(i .eq. nlat-1)then
          dpa=0.d0
        else
          dpa=(pos(i)-pos(i-1))*(fre-1.d0)
        endif
      else
        dpa=0.d0
      endif
      if(only .and. .not. tmatch(pname(ilist(2,latt(i))),patt))then
      elseif(monly .and. .not. monel(imon,nmon,i,m1))then
      else
        if(kv .le. 27 .or. kv .eq. 35)then
          jp=jp+1
          y(1,jp)=pos(i)+pa+dpa
          y(2,jp)=drwfun(kv,i,twiss,gammab,idp,amass,nlat,ndim)
        elseif(kv .le. 32) then
          jp=jp+1
          y(1,jp)=pos(i)+pa+dpa
          if(trsize)then
            call tftrb(latt,i,sv)
            beam(1)=sv(3)
            beam(4)=sv(5)
            beam(6)=sv(4)
          else
            call tfbeam(i,theta,beam)
          endif
          if(kv.ne.32) then
            ia=(kv-27)*(kv-26)/2
            y(2,jp)=sqrt(beam(ia))
          else
            y(2,jp)=atan2(-2d0*beam(4),beam(1)-beam(6)) /2d0
          endif
        elseif(kv.le.34) then
          if(monel(imon,nmon,i,m1)) then
            jp=jp+1
            y(1,jp)=pos(i)+pa+dpa
            ll=imon(m1,2)
            nq=imon(ll,4)
            y(2,jp)=twiss(i,idp,2*kv-51)-twiss(i,ndim,2*kv-51)
     1           -rlist(latt(nq)+kv-28)-emon(ll,kv-32)
          endif
        endif
      endif
      if(ls .gt. le .and. .not. fin) then
        pa=pos(nlat)-pos(1)
        istrt=1
        istop=le
        fin=.true.
        goto 111
      endif
c      write(*,*)'tdrwdt ',jp,np
      ymin= 1d20
      ymax=-1d20
      do i=1,jp
        ymin=min(y(2,i),ymin)
        ymax=max(y(2,i),ymax)
      enddo
      if(kv .eq. mfitbx .or. kv .eq. mfitby
     $     .or. (kv .ge. 28 .and .kv .le. 31))then
        ymin=0d0
      endif
      return
      end

      real*8 function drwfun(kv,i,twiss,gammab,idp,
     $     amass,nlat,ndim)
      use tffitcode
      implicit none
      integer*4 kv,i,idp,nlat,ndim
      real*8 twiss(nlat,-ndim:ndim,ntwissfun),gammab(nlat),
     $     amass,u(4),v(4)
      if(kv .le. mfitdpy) then
        drwfun=twiss(i,idp,kv)
      elseif(kv .le. 22) then
c        u(1)=twiss(i,idp,mfitex)
c        u(2)=twiss(i,idp,mfitepx)
c        u(3)=twiss(i,idp,mfitey)
c        u(4)=twiss(i,idp,mfitepy)
c        call mc2to4(twiss,idp,i,u,v)
        call tgetphysdisp(i,v)
        drwfun=v(kv-18)
      elseif(kv .eq. 23) then
        drwfun=twiss(i,idp,mfitr1)*twiss(i,idp,mfitr4)
     z       -twiss(i,idp,mfitr2)*twiss(i,idp,mfitr3)
      elseif(kv .le. 27) then
        u(1)=twiss(i,idp,mfitdx)
        u(2)=twiss(i,idp,mfitdpx)
        u(3)=twiss(i,idp,mfitdy)
        u(4)=twiss(i,idp,mfitdpy)
        call mc4to2(twiss,idp,i,u,v)
        drwfun=v(kv-23)
      elseif(kv .eq. 35) then
        drwfun=amass*gammab(i)
      else
        drwfun=0.d0
      endif
      return
      end
