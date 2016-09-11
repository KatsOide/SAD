      subroutine pvbump(word,wordp,latt,twiss,mult,master,
     $     istr,nstr,nlist,kfit,ifitp,mfitp,fitval,nfc,lfno)
      use tfstk
      use ffs
      use tffitcode
      parameter (mfitc1=34,mfitc2=30,meminc=100,vinit0=1d-3)
      logical exist,delete,takemb,yplane
      character*(*) word,wordp
      character name*8,vname*12,autofg*11
      character*8  nlist(mfit1)
      integer*8 latt(nlat)
      dimension twiss(nlat,-ndim:ndim,ntwissfun),mult(*),master(*)
      dimension istr(*)
      dimension kfit(nfc),ifitp(nfc),mfitp(nfc),fitval(nfc)
      include 'inc/common.inc'
      common /cordefv/ ipvbmp,nvbmp
c     begin initialize for preventing compiler warning
      k=0
c     end   initialize for preventing compiler warning
      delete=.false.
      takemb=.false.
c
      ncor=getva(exist)
      if(exist) then
        if(ncor.lt.0) then
          delete=.true.
        elseif(ncor.eq.0) then
          takemb=.true.
        else
          ncor=2*((ncor+1)/2)
        endif
      else
        takemb=.true.
      endif
      call getwdl2(word,wordp)
      if(word.eq.' ') then
c ..... write status ....
        do 9 i=1,nvbmp
          lp=ilist(mod(i-1,2)+1,ipvbmp+(i-1)/2)
          nc=ilist(1,lp)
          lp1=lp+3
          lp2=lp+3+(nc+1)/2
          lp3=lp+3+2*((nc+1)/2)
          lv=ilist(2,lp)
          kv=ilist(1,lp+1)
          call mcchar(ilist(2,lp+1),vname,3)
          call elname(ilist(mod(lv-1,2)+1,lp2+(lv-1)/2),name)
          call mbufw(vname//': fit '//name,.false.,lfno)
          call mbufw(nlist(kv)(1:lene(nlist(k)))//
     $         autofg(rlist(lp3+lv-1),'S10.3'),.false.,lfno)
          do 8 j=1,nc
            if(j.ne.lv) then
              call elname(latt,ilist(mod(j-1,2)+1,lp2+(j-1)/2),mult,
     $             name)
              k=ilist(mod(j-1,2)+1,lp1+(j-1)/2)
              if(vname.eq.name) then
                call mbufw(nlist(k)(1:lene(nlist(k)))//
     $               autofg(rlist(lp3-1+j),'S3.1'),.false.,lfno)
              else
                call mbufw('fit '//name(1:lene(name))//' '
     $               //nlist(k)(1:lene(nlist(k)))//
     $               autofg(rlist(lp3-1+j),'S3.1'),.false.,lfno)
              endif
              vname=name
            endif
 8        continue
          call mbufw(' ',.true.,lfno)
 9      continue
        return
      endif
c.... read variable type ....      
      do 10 i=1,mfit1
        if(word.eq.nlist(i)) then
          kv=i
          goto 11
        endif
 10   continue
      if(delete) then
        kv=0
        goto 12
      endif
      return
 11   vinit=getva(exist)
      if(.not.exist) vinit=0d0
      call getwdl2(word,wordp)
      if(word.eq.' ') return
c.... void definition .....
 12   if(delete) then
        call pvbump3(kv,word)
        call getwdl2(word,wordp)
        return
      endif 
      exist=.false.
      if(takemb) then
c...... get conditions from matching-condition buffer ...
        nc=0
        ltarget=0
        do 20 i=1,nfc
          if(mfitp(i).ne.0) then
            if(kfit(i) .ge. mfitc1 .and. kfit(i) .le. mfitc1+3 .or.
     &         kfit(i) .ge. mfitc2 .and. kfit(i) .le. mfitc2+3  ) then
              nc=nc+1
              call elname(ifitp(i),name)
              if(wordp.eq.name) then
                if(kv.eq.kfit(i)) then
                  ltarget=nc
                  vname=name(1:lene(name))//'_'//nlist(kv)
                  exist=.true.
                endif
              endif
            endif
          endif
 20     continue
        if(nc.eq.0) return
        if(.not.exist) return
        nvbmp=nvbmp+1
        if(nvbmp.eq.1) then
          ipvbmp=italoc(meminc)
        elseif(mod(nvbmp,meminc).eq.1) then
          call palocx(ipvbmp,nvbmp-1,nvbmp-1+meminc)
        endif
        ip=italoc(3+2*((nc+1)/2)+nc)
        ilist(mod(nvbmp-1,2)+1,ipvbmp+(nvbmp-1)/2)=ip
        ilist(1,ip)=nc
        ilist(2,ip)=ltarget
        ilist(1,ip+1)=kv
c       print *,vname,nc,ltarget,kv,nvbmp,ipvbmp,ip
        call mcchar(vname,ilist(2,ip+1),3)
        call pvbump1(kfit,ifitp,mfitp,fitval,nfc,
     &       rlist(ip+3),rlist(ip+3+(nc+1)/2),rlist(ip+3+2*((nc+1)/2)))
      else
c...... read bump spec ....
c       print *,mfalloc(-1)
        call mbmpf(word,wordp,latt,mult,master,nb)
        isb=italoc((ncor+4)*nb)
        yplane=kv.eq.mfitc1+2 .or. kv.eq.mfitc1+3 .or.
     $         kv.eq.mfitc2+2 .or. kv.eq.mfitc2+3
        call mbmp(latt,twiss,mult,master,rlist(isb),nb,istr,nstr,ncor,
     $       yplane,iret,lfno)
c+++++++ debug ++++++++
c       do i=1,nb
c         j=(ncor+4)*(i-1)
c         k=ilist(mod(j+1,2)+1,isb+(j+1)/2)+2
c         write(*,'(10i4)')(ilist(mod(j+l,2)+1,isb+(j+l)/2),l=0,k-1)
c       enddo
c++++++++++++++++++++++
        do 30 i=1,nb
          j=(ncor+4)*(i-1)
          lt=ilist(mod(j,2)+1,isb+j/2)
          call elname(lt,name)
c ....... delete (kv,name) from 'vbump' list, if it exists.          
          call pvbump3(kv,name)
          vname=name(1:lene(name))//'_'//nlist(kv)
          nvbmp=nvbmp+1
          if(nvbmp.eq.1) then
            ipvbmp=italoc(meminc)
          elseif(mod(nvbmp,meminc).eq.1) then
            call palocx(ipvbmp,nvbmp-1,nvbmp-1+meminc)
          endif
          if(vinit.eq.0d0) then
            if(kv.eq.mfitc1+1 .or. kv.eq.mfitc2+1) then
              vinit=vinit0*(1d0-twiss(lt,0,1))/twiss(lt,0,2)
            elseif(kv.eq.mfitc1+3 .or. kv.eq.mfitc2+3) then
              vinit=vinit0*(1d0-twiss(lt,0,4))/twiss(lt,0,5)
            else
              vinit=vinit0
            endif
          endif
          call pvbump2(ilist(mod(j,2)+1,isb+j/2),ncor,kv,vinit,vname,
     $         rlist(ipvbmp),nvbmp,istr,yplane)
c         do ii=1,nb
c           jj=(ncor+4)*(ii-1)
c           kk=ilist(mod(jj+1,2)+1,isb+(jj+1)/2)+2
c           write(*,'(10i4)')(ilist(mod(jj+ll,2)+1,isb+(jj+ll)/2),ll=0,
c    $           kk-1)
c         enddo
 30     continue
        call tfree(int8(isb))
c       print *,mfalloc(-1)
        return
      endif 
c     print *,ilist(1,ip),ilist(2,ip),ilist(1,ip+1)
c     call mcchar(ilist(2,ip+1),name,2)
c     print *,name
c     write(*,*) (ilist(mod(i-1,2)+1,ip+3+(i-1)/2),i=1,nc)
c     write(*,*) (ilist(mod(i-1,2)+1,ip+3+(nc+1)/2+(i-1)/2),i=1,nc)
c     write(*,*) (rlist(ip+3+2*((nc+1)/2)+i-1),i=1,nc)
      call getwdl2(word,wordp)
      return
      end
