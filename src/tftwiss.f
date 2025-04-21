      module twissf
      use tfstk
      use ffs
      use tffitcode
      use ffs_fit, only:nlist
      use ffs_pointer, only:latt,idelc,idtypec,pnamec,compelc,direlc,twiss
      use sad_main
      use tflinepcom
      use gfun

      contains
      type (sad_descriptor) function tflinek(keyword,ia,cmp,ref) result(kx)
      implicit none
      integer*4 ,intent(out):: ia
      logical*4 ,intent(in):: ref
      character*(*),intent(in):: keyword
      type (sad_comp),pointer,intent(inout) ::cmp
      integer*4 nc
      integer*8 ip
      character*64 key1
      nc=len(keyword)
      if(keyword == '@GEO')then
        key1(1:3)='GEO'
        nc=3
      else
        key1(1:nc)=keyword(1:nc)
      endif
      if(ia < nlat)then
        kx=tfkeyv(int(ia),key1(1:nc),ip,cmp,ref,.false.)
        if(.not. ref)then
          cmp%update=cmp%nparam <= 0
          kx%k=ktfref+ip
        endif
        tparaed=.false.
      else
        kx%k=0
      endif
      return
      end function tflinek

      end module

      subroutine tftwiss(isp1,kx,ref,irtc)
      use twissf
      use geto
      implicit none
      type (sad_descriptor) ,intent(out):: kx
      type (sad_dlist), pointer :: klx
      type (sad_rlist), pointer :: ktl,kll
      integer*4 ,parameter ::nkey=mfito
      integer*8 kax,kaxi,itoff
      integer*4 ,intent(in):: isp1
      integer*4 ,intent(out):: irtc
      integer*4 narg,i,m,nc,isp0,nd,kt,itfmessage,lenw,icol
      real*8 ftwiss(ntwissfun),tfgettwiss
      logical*4 ,intent(in):: ref
      logical*4 over,dref
      character*(MAXPNAME+16) keyword,tfgetstrs
      narg=isp-isp1
      keyword=tfgetstrs(ktastk(isp1+1),nc)
      if(nc <= 0)then
        irtc=itfmessage(9,'General::wrongtype',
     $       '"Keyword-string for #1"')
        return
      endif
      irtc=0
      icol=0
      dref=.false.
      kx%k=0
      call capita(keyword(1:nc))
      if(keyword == 'LENGTH')then
        if(narg > 1)then
          irtc=itfmessage(9,'General::narg','"1"')
          return
        endif
        kx=dfromr(dble(nlat))
      elseif(keyword == '*' .or. keyword == 'ALL'
     $       .or. keyword == 'DALL' .or. keyword == 'RALL')then
        if(keyword == 'RALL')then
          icol=-1
        elseif(keyword == 'DALL')then
          dref=.true.
        endif
        if(narg == 1)then
          kax=ktadaloc(-1,nlat,klx)
          if(dref)then
            do i=1,nlat
              kaxi=ktatwissaloc(0,ktl)
              klx%dbody(i)%k=ktflist+kaxi
              call tfgetdref(twiss(i,0,1:ntwissfun),
     $             twiss(i,-1,1:ntwissfun),ktl%rbody(1:ntwissfun))
            enddo
          else
            do i=1,nlat
              kaxi=ktatwissaloc(0,ktl)
              klx%dbody(i)%k=ktflist+kaxi
              ktl%rbody(1:ntwissfun)=twiss(i,icol,1:ntwissfun)
            enddo
          endif
        elseif(narg == 2)then
          call tflinestk(dtastk(isp),narg,isp0,irtc)
          if(irtc /= 0)then
            return
          endif
          if(isp == isp0+1)then
            kax=ktatwissaloc(-1,ktl)
            if(vstk2(isp) == 0.d0)then
              if(dref)then
                call tfgetdref(twiss(itastk(2,isp),0,1:ntwissfun),
     $               twiss(itastk(2,isp),-1,1:ntwissfun),
     $               ktl%rbody(1:ntwissfun))
              else
                ktl%rbody(1:ntwissfun)
     $               =twiss(itastk(2,isp),icol,1:ntwissfun)
              endif
            elseif(icol /= 0 .or. dref)then
              go to 9000
            else
              call qtwissfrac(rlist(kax+1:kax+ntwissfun),itastk(2,isp),
     $             vstk2(isp),over)
            endif
          else
            m=isp-isp0
            kax=ktadaloc(-1,m,klx)
            do i=1,m
              kaxi=ktatwissaloc(0,ktl)
              klx%dbody(i)%k=ktflist+kaxi
              if(vstk2(isp0+i) == 0.d0)then
                if(dref)then
                  call tfgetdref(twiss(itastk(2,isp0+i),0,1:ntwissfun),
     $                 twiss(itastk(2,isp0+i),-1,1:ntwissfun),
     $                 ktl%rbody(1:ntwissfun))
                else
                  ktl%rbody(1:ntwissfun)=
     $                 twiss(itastk(2,isp0+i),icol,1:ntwissfun)
                endif
              elseif(icol /= 0 .or. dref)then
                go to 9000
              else
                call qtwissfrac(rlist(kaxi+1:kax+ntwissfun),
     $               itastk(2,isp0+i),vstk2(isp0+i),over)
              endif
            enddo
          endif
          isp=isp0
        else
          irtc=itfmessage(9,'General::narg','"1 or 2"')
          return
        endif
        kx%k=ktflist+kax
      elseif(keyword == 'FUNCTIONS')then
      else
        findkey: do
          do i=1,nkey
            if(keyword == nlist(i))then
              kt=i
              exit findkey
            endif
          enddo
          if(keyword(1:1) == 'D' .or. keyword (1:1) == 'R')then
            if(keyword(1:1) == 'R')then
              icol=-1
            else
              dref=.true.
            endif
            do i=1,nkey
              if(keyword(2:) == nlist(i))then
                kt=i
                exit findkey
              endif
            enddo
            icol=0
            dref=.false.
          endif
          if(keyword(1:3) == 'SIG' .or. keyword(1:4) == 'SIZE'
     $         .or. keyword == 'GAMMA'
     $         .or. keyword == 'GAMMABETA'
     $         .or. keyword == 'S')then
            call tfline(isp1,kx,ref,irtc)
          else
            irtc=itfmessage(9,'General::wrongval',
     $           '"#1 ('//keyword(1:lenw(keyword))//
     $           ') is","to be name of optical function"')
          endif
          return
        enddo findkey
        if(narg == 1)then
          kax=ktavaloc(-1,nlat,kll)
          if(kt <= ntwissfun)then
            if(dref)then
              select case (kt)
              case (mfitbx,mfitby,mfitbz)
                kll%rbody(1:nlat)=
     $               (twiss(1:nlat,0,kt)-twiss(1:nlat,-1,kt))
     $               /twiss(1:nlat,-1,kt)
              case default
                kll%rbody(1:nlat)=
     $               twiss(1:nlat,0,kt)-twiss(1:nlat,-1,kt)
              end select
            elseif(ref)then
              kll%rbody(1:nlat)=twiss(1:nlat,icol,kt)
            else
              itoff=((2*ndim+1)*(kt-1)+ndim*(icol+1))*nlat+iftwis
                klist(kax+1:kax+nlat)=
     $             (/(ktfref+itoff+i-1,i=1,nlat)/)
            endif
          elseif(kt .ge. mfitpex .and. kt <= mfitpepy .or.
     $           kt .ge. mfitpzx .and. kt <= mfito)then
            do i=1,nlat
              call tfgettwiss1(i,icol,kt,kll%dbody(i),dref,ref)
            enddo
          endif
          kx%k=ktflist+kax
        elseif(narg == 2)then
          call tflinestk(dtastk(isp),narg,isp0,irtc)
          if(irtc /= 0)then
            return
          endif
          if(isp == isp0+1)then
            if(vstk2(isp) == 0.d0)then
              call tfgettwiss1(itastk(2,isp),icol,kt,kx,dref,ref)
            elseif(icol /= 0 .or. dref)then
              go to 9000
            else
              call qtwissfrac(ftwiss,itastk(2,isp),
     $             vstk2(isp),over)
              kx=dfromr(tfgettwiss(kt,ftwiss))
c              write(*,*)'twiss ',rfromk(kx),kt,
c     $             itastk(2,isp),vstk2(isp)
            endif
          else
            m=isp-isp0
            kax=ktavaloc(-1,m,kll)
            kx%k=ktflist+kax
            if(kt <= ntwissfun)then
              if(ref)then
                do i=1,m
                  if(vstk2(isp0+i) == 0.d0)then
                    call tfgettwiss1(itastk(2,isp0+i),
     $                   icol,kt,kll%dbody(i),dref,ref)
c                    kll%rbody(i)=twiss(itastk(2,isp0+i),icol,kt)
                  elseif(icol /= 0 .or. dref)then
                    go to 9000
                  else
                    call qtwissfrac(ftwiss,itastk(2,isp0+i),
     $                   vstk2(isp0+i),over)
                    kll%rbody(i)=ftwiss(kt)
                  endif
                enddo
              else
                itoff=((2*ndim+1)*(kt-1)+ndim)*nlat+iftwis
                do i=1,m
                  if(vstk2(isp0+i) == 0.d0)then
                    klist(kax+i)=ktfref+itoff+itastk(2,isp0+i)-1
                  else
                    go to 9000
                  endif
                enddo
              endif
            elseif(kt .ge. mfitpex .and. kt. le. mfitpepy
     $             .or. kt .ge. mfitpzx .and. kt. le. mfito)then
              do i=1,m
                if(vstk2(isp0+i) == 0.d0)then
                  call tfgettwiss1(itastk(2,isp0+i),
     $                   icol,kt,kll%dbody(i),dref,ref)
                elseif(icol /= 0 .or. dref)then
                  go to 9000
                else
                  call qtwissfrac(ftwiss,itastk(2,isp0+i),
     $                 vstk2(isp0+i),over)
                  kll%rbody(i)=tphysdisp(kt,ftwiss)
                endif
              enddo
            endif
          endif
          isp=isp0
        elseif(narg == 3)then
          keyword=tfgetstrs(ktastk(isp-1),nc)
          if(nc <= 0)then
            irtc=itfmessage(9,'General::wrongtype',
     $           '"name of component for #2"')
            return
          endif
          call capita(keyword(1:nc))
          if(kt <= ntwissfun)then
            if(keyword == 'SET')then
              kx=dtastk(isp)
              if(ktflistq(kx,klx))then
                nd=min(klx%nl,nlat)
                twiss(1:nd,icol,kt)=klx%rbody(1:nd)
c                rlist(itoff:itoff+nd-1)=klx%rbody(1:nd)
                return
              endif
            endif
          endif
        endif
      endif
      return
 9000 irtc=itfmessage(9,'General::wrongval',
     $     '"Fractional component # not supported"')
      return
      end

      subroutine tfgetdref(ft,ft0,r)
      use ffs
      use tffitcode
      implicit none
      real*8 ,intent(in):: ft(ntwissfun),ft0(ntwissfun)
      real*8 ,intent(out):: r(ntwissfun)
      r=ft-ft0
      r(mfitbx)=r(mfitbx)/ft0(mfitbx)
      r(mfitby)=r(mfitby)/ft0(mfitby)
      r(mfitbz)=r(mfitbz)/ft0(mfitbz)
      return
      end

      subroutine tfgettwiss1(i,icol,kt,kx,dref,ref)
      use tfstk
      use ffs
      use ffs_pointer
      use tffitcode
      use gfun
      implicit none
      type (sad_descriptor) , intent(out):: kx
      integer*4 ,intent(in):: i,icol,kt
      integer*8 itoff
      logical*4 ,intent(in):: ref,dref
      real*8 pe(4),pe0(4)
      select case (kt)
      case (mfitbx,mfitby,mfitbz)
        if(dref)then
          kx=dfromr((twiss(i,0,kt)-twiss(i,-1,kt))/twiss(i,-1,kt))
        elseif(ref)then
          kx=dfromr(twiss(i,icol,kt))
        else
          itoff=((2*ndim+1)*(kt-1)+ndim*(icol+1))*nlat+iftwis
          kx%k=ktfref+itoff+i-1
        endif
      case (mfitpex,mfitpepx,mfitpey,mfitpepy)
        if(dref)then
          call tgetphysdispi(i,0,pe)
          call tgetphysdispi(i,-1,pe0)
          kx=dfromr(pe(kt-mfitpex+1)-pe0(kt-mfitpex+1))
        else
          call tgetphysdispi(i,icol,pe)
          kx=dfromr(pe(kt-mfitpex+1))
        endif
      case (mfitpzx,mfitpzpx,mfitpzy,mfitpzpy)
        if(dref)then
          call tgetphysdispzi(i,0,pe)
          call tgetphysdispzi(i,-1,pe0)
          kx=dfromr(pe(kt-mfitpzx+1)-pe0(kt-mfitpzx+1))
        else
          call tgetphysdispzi(i,icol,pe)
          kx=dfromr(pe(kt-mfitpzx+1))
        endif
      case (mfitgmx,mfitgmy,mfitgmz)
        if(dref)then
          kx=dfromr((tgetgm(kt,i,0)-tgetgm(kt,i,-1))/
     $         tgetgm(kt,i,-1))
        else
          kx=dfromr(tgetgm(kt,i,icol))
        endif
      case (mfitbmagx,mfitbmagy,mfitbmagz)
        kx=tgetbmag(i,kt)
      case default
        if(dref)then
          kx=dfromr(twiss(i,0,kt)-twiss(i,-1,kt))
        elseif(ref)then
          kx=dfromr(twiss(i,icol,kt))
        else
          itoff=((2*ndim+1)*(kt-1)+ndim*(icol+1))*nlat+iftwis
          kx%k=ktfref+itoff+i-1
        endif
      end select
      return
      end

      real*8 function tfgettwiss(kt,ftwiss)
      use tfstk
      use tffitcode
      use gfun
      implicit none
      integer*4 ,intent(in):: kt
      real*8 ,intent(in):: ftwiss(ntwissfun)
      tfgettwiss=0.d0
      if(kt <= ntwissfun)then
        tfgettwiss=ftwiss(kt)
      elseif(kt .ge. mfitpex .and. kt <= mfitpepy .or.
     $       kt .ge. mfitpzx .and. kt <= mfito)then
        tfgettwiss=tphysdisp(kt,ftwiss)
      endif
      if(ktfenanq(tfgettwiss))then
        tfgettwiss=rtfnan
      endif
      return
      end

      subroutine tfelement(isp1,kx,ref,irtc)
      use tfstk
      use ffs
      use tffitcode
      use repl, only:tfgetoption
      implicit none
      type (sad_descriptor) ,intent(out):: kx
      integer*4 ,intent(in):: isp1
      integer*4 ,intent(out):: irtc
      integer*4 i,narg,nc,isp0,m,itfmessage,ispa
      character*(MAXPNAME+16) keyword,tfgetstrs
      logical*4 ,intent(in):: ref
      logical*4 saved
      narg=isp-isp1
      irtc=0
      if(narg <= 0)then
        irtc=itfmessage(9,'General::narg','"1, 2, or 3"')
        return
      endif
      keyword=tfgetstrs(ktastk(isp1+1),nc)
      if(nc <= 0)then
        irtc=itfmessage(9,'General::wrongtype',
     $       '"Keyword for #1"')
        return
      endif
      call capita(keyword(1:nc))
      if(keyword == 'LENGTH')then
        if(narg > 1)then
          go to 9010
        endif
        kx=dfromr(dble(nele))
      elseif(keyword == 'EXPAND')then
        if(narg /= 1)then
          go to 9010
        endif
        call tffsadjust
        kx=dxnullo
      else
        if(narg > 2)then
          if(narg == 3)then
            call tfgetoption('Saved',dtastk(isp),kx,irtc)
            if(irtc /= 0)then
              return
            endif
            if(ktfrealq(kx))then
              saved=kx%k /= 0
              narg=2
            else
              go to 9010
            endif
          else
            go to 9010
          endif
        else
          saved=.false.
        endif
        call tfelementstk(dtastk(isp1+2),isp0,narg,irtc)
        if(irtc /= 0)then
          isp=isp0
          return
        endif
        if(isp == isp0+1)then
          call tfelement1(itastk(1,isp),itastk(2,isp),
     $         kx,keyword,saved,ref,irtc)
        else
          m=isp-isp0
          ispa=isp
          do i=1,m
            call tfelement1(itastk(1,isp0+i),itastk(2,isp0+i),
     $           kx,keyword,saved,ref,irtc)
            isp=isp+1
            dtastk(isp)=kx
          enddo
          kx=kxmakelist(ispa)
        endif
        isp=isp0
      endif
      return
 9010 irtc=6
      return
      end

      subroutine tfelement1(it,ia,kx,keyword,saved,ref,irtc)
      use tfstk
      use ffs
      use tffitcode
      use ffs_pointer, only:latt,idelc,idtypec,idvalc,sad_comp,compelc
      use tflinepcom
      use tparastat, only:ndivelm
      implicit none
      type (sad_descriptor) kx
      type (sad_comp), pointer :: cmp
      type (sad_rlist), pointer :: kl
      integer*8 iax
      integer*4 ,intent(out):: irtc
      integer*4 ,intent(in):: it,ia
      integer*4 id,lenw,iv,isps,l
      character*(*) ,intent(in):: keyword
      character*(MAXPNAME) key,tfkwrd
      logical*4 ,intent(in):: saved,ref
      irtc=0
      select case (keyword)
      case('NAME')
        id=idelc(ia)
        kx=kxsalocb(-1,pname(id),lpname(id))
      case('VALUE')
        iv=nelvx(it)%ival
        if(iv > 0)then
          if(saved)then
            iax=idvalc(ia)+iv
            if(ref)then
              kx=dlist(iax)
            else
              kx%k=ktfref+iax
            endif
c            kx%k=merge(klist(iax),ktfref+iax,ref)
          else
            iax=latt(ia)+iv
            if(ref)then
              kx=dfromr(rlist(iax)/rlist(iferrk+2*(ia-1)))
            else
              kx%k=ktfref+iax
              call compelc(ia,cmp)
              cmp%update=cmp%nparam <= 0
            endif
          endif
        else
          kx%k=0
        endif
      case('DEFAULT')
        iv=nelvx(it)%ival
        if(iv == 0)then
          key='                                '
        else
          key=tfkwrd(idtypec(ia),iv)
        endif
        kx=kxsalocb(-1,key,lenw(key))
      case('DEFAULT$SUM')
        iv=nelvx(it)%ival
        if(iv == 0)then
          key=' '
        else
          key=tfkwrd(idtypec(ia),iv)
          key=key(1:lenw(key))//"$SUM"
        endif
        kx=kxsalocb(-1,key,lenw(key))
      case('KEYWORDS','KEYWORDS_ALL')
        l=0
        id=idtypec(ia)
        isps=isp
        call tftypekeystk(id,keyword == 'KEYWORDS_ALL')
        kx=kxmakelist(isps)
        isp=isps
      case('TYPE')
        kx=dfromr(dble(idtypec(ia)))
      case('TYPENAME')
        key=pname(kytbl(0,idtypec(ia)))
        kx=kxsalocb(-1,key(2:),lenw(key)-1)
      case('POSITION')
        kx=dfromr(dble(it))
      case('COMPONENT')
        call elcompl(it,kl)
        kx=sad_descr(kl)
      case('NDIV')
        kx%x(1)=dble(ndivelm(ia))
      case default
        kx=tfkeyv(-it,keyword,iax,cmp,ref,saved)
        if(.not. ref)then
          kx%k=ktfref+iax
          if(.not. saved)then
            cmp%update=cmp%nparam <= 0
          endif
        endif
      end select
      return
      end

      subroutine tfelementstk(k,isp0,narg,irtc)
      use geto
      use ffs
      use tffitcode
      use ffs_pointer, only:idelc,pnamec
      implicit none
      type (sad_descriptor) ,intent(in):: k
      integer*4 ,intent(out):: isp0,irtc
      integer*4 ,intent(in):: narg
      integer*4 iv,nc,ifany1,i,itfmessage,j,ielmh
      character*32768 name
      logical*4 tmatch
      isp0=isp
      if(ktfrealq(k,iv) .and. narg == 2)then
        if(iv .lt. 0)then
          iv=nele+iv+1
        endif
        if(iv <= 0 .or. iv > nele)then
          irtc=itfmessage(9,'General::wrongnum',
     $         '"positive and less than length of beam line"')
          return
        endif
        isp=isp+1
        itastk(1,isp)=iv
        itastk(2,isp)=nelvx(iv)%klp
c        itastk(2,isp)=ilist(iv,ifklp)
        irtc=0
      else
        if(narg == 1)then
          name='*'
          nc=1
        else
          call tfgetstrns(k,name,nc)
          if(nc <= 0)then
            irtc=itfmessage(9,'General::wrongtype',
     $           '"Name of component"')
            return
          endif
          if(convcase)then
            call capita1(name(1:nc))
          endif
        endif
        irtc=0
        if(name(1:nc) /= '***' .and.
     $       ifany1(name(1:nc),nc,'*%{<|',1) > 0)then
          do i=1,nele
c            write(*,*)'elementstk',i,nele,pname(idelc(ilist(i,ifklp)))
            if(tmatch(pnamec(nelvx(i)%klp),
     $           name(1:nc)))then
              isp=isp+1
              itastk(1,isp)=i
              itastk(2,isp)=nelvx(i)%klp
            endif
          enddo
        else
          j=ielmh(name(1:nc),0)
          if(j /= 0)then
            i=ilist(j,ifele1)
            isp=isp+1
            itastk(1,isp)=i
            itastk(2,isp)=nelvx(i)%klp
            return
          endif
        endif
      endif
      return
      end

      subroutine tfline(isp1,kx,ref,irtc)
      use tfstk
      use ffs
      use tffitcode
      use ffs_pointer, only:latt,icomp
      implicit none
      type (sad_descriptor) ,intent(out):: kx
      integer*4 ,intent(in):: isp1
      integer*4 ,intent(out):: irtc
      integer*4 ie,iv,k,j,m,ispa
      integer*4 i,narg,nc,isp0,itfmessage
      character*(MAXPNAME+16) keyword,tfgetstrs
      logical*4 ,intent(in):: ref
      narg=isp-isp1
      keyword=tfgetstrs(ktastk(isp1+1),nc)
      if(nc <= 0)then
        irtc=itfmessage(9,'General::wrongtype',
     $       '"Keyword for #1"')
        return
      endif
      irtc=0
      call capita(keyword(1:nc))
      if(keyword == 'LENGTH')then
        if(narg > 1)then
          irtc=itfmessage(9,'General::narg','"1"')
        endif
        kx=dfromr(dble(nlat))
      elseif(keyword(1:nc) == 'EXPAND')then
        if(narg /= 1)then
          irtc=itfmessage(9,'General::narg','"1"')
        endif
        do i=1,nlat-1
          ie=ilist(icomp(i),ifele1)
          if(ie > 0)then
            iv=nelvx(ie)%ival
c            iv=ilist(ie,ifival)
            if(iv > 0)then
              k=nelvx(ie)%klp
c              k=ilist(ie,ifklp)
              rlist(latt(i)+iv)=
     $             rlist(ifcoup+i-1)*
     $             rlist(iferrk+i*2-2)/rlist(iferrk+k*2-2)*
     $             rlist(latt(k)+iv)
            endif
          endif
        enddo
        kx%k=ktfoper+mtfnull
      else
        call tflinestk(dtastk(isp),narg,isp0,irtc)
        if(irtc /= 0)then
          isp=isp0
          return
        endif
        m=isp-isp0
        if(m == 1)then
          call tfline1(isp,kx,keyword(1:nc),ref,irtc)
        else
          ispa=isp
          do j=1,m
            isp=isp+1
            call tfline1(isp0+j,dtastk(isp),keyword(1:nc),ref,irtc)
            if(irtc /= 0)then
              return
            endif
          enddo
          kx=kxmakelist(ispa)
        endif
        isp=isp0
      endif
      return
      end

      subroutine tfline1(isp1,kx,keyword,ref,irtc)
      use twissf
      use sad_main
      use tflinepcom
      use geolib
      use tparastat,only:ndivelm
      implicit none
      type (sad_descriptor) ,intent(out):: kx
      type (sad_comp), pointer ::cmp
      integer*8 kax,ktfgeol,kai,i,ip,j
      integer*4 ,intent(out):: irtc
      integer*4 ,intent(in):: isp1
      integer*4 lenw,lxp,ibz,ia
      real*8 v,beam(42),xp,fr,
     $     gv(3,4),ogv(3,4),cod(6),vtwiss(ntwissfun),tfbzs,tfbzt
      character*(*) ,intent(in):: keyword
      character*64 name
      integer*4 lv
      logical*4 ,intent(in):: ref
      logical*4 over
c      iaidx(m,n)=int(((m+n+abs(m-n))**2+2*(m+n)-6*abs(m-n))/8)
      irtc=0
      ip=itastk(1,isp1)
      ia=itastk(2,isp1)
      v=vstk2(isp1)
      select case (keyword)
      case ('S','LENG')
        kx=dfromr(rlist(ifpos+ia-1)*(1.d0-v)+
     $       rlist(ifpos+min(nlat-1,ia))*v)
      case('GAMMABETA')
        kx=dfromr(rlist(ifgamm+ia-1)*(1.d0-v)+
     $       rlist(ifgamm+min(nlat-1,ia))*v)
      case('GAMMA')
        kx=dfromr(sqrt(1.d0+(rlist(ifgamm+ia-1)*(1.d0-v)+
     $       rlist(ifgamm+min(nlat-1,ia))*v)**2))
      case('MULT')
        kx=dfromr(dble(ilist(ia,ifmult)))
      case('TYPE')
        if(ia == nlat)then
          kx=dxzero
        else
          kx=dfromr(dble(idtypec(ia)))
        endif
      case('TYPENAME')
        if(ia == nlat)then
          kx=dxnulls
        else
          name=pname(kytbl(0,idtypec(ia)))
          kx=kxsalocb(-1,name(2:),lenw(name)-1)
        endif
      case('NAME')
        call elname(ia,name)
        kx=kxsalocb(-1,name,lenw(name))
      case('ELEMENT')
        if(ia == nlat)then
          name='$$$'
        else
          name=pnamec(ia)
        endif
        kx=kxsalocb(-1,name,lenw(name))
      case('POSITION')
        kx=dfromr(dble(ia))
      case('GEO')
        xp=v+ia
        lxp=int(xp)
        fr=xp-lxp
        kx%k=ktflist+ktfgeol(tfgeofrac(lxp,fr,irtc))
      case('GX' ,'GY' ,'GZ' ,'GCHI1' ,'GCHI2' ,'GCHI3')
        xp=v+ia
        lxp=int(xp)
        fr=xp-lxp
        gv=tfgeofrac(lxp,fr,irtc)
        select case (keyword)
        case('GX')
          kx=dfromr(gv(1,4))
        case('GY')
          kx=dfromr(gv(2,4))
        case('GZ')
          kx=dfromr(gv(3,4))
        case('GCHI1')
          kx=dfromr(tfchi(gv,1))
        case('GCHI2')
          kx=dfromr(tfchi(gv,2))
        case('GCHI3')
          kx=dfromr(tfchi(gv,3))
        end select
      case('OGEO')
        xp=v+ia
        lxp=int(xp)
        fr=xp-lxp
        j=ifgeo+(lxp-1)*12
        if(fr == 0.d0)then
          cod=twiss(lxp,0,mfitdx:mfitddp)
          call tmov(rlist(j),gv,12)
        else
          levele=levele+1
          call qtwissfracgeo(vtwiss,gv,lxp,fr,.true.,over)
          cod=vtwiss(mfitdx:mfitddp)
          lv=itfdownlevel()
        endif
        kx%k=ktflist+ktfgeol(tforbitgeo(gv,cod))
      case('OGX' ,'OGY' ,'OGZ' ,'OCHI1' ,'OCHI2' ,'OCHI3')
        xp=v+ia
        lxp=int(xp)
        fr=xp-lxp
        j=ifgeo+(lxp-1)*12
        if(fr == 0.d0)then
          call tmov(rlist(j),gv,12)
          cod=twiss(lxp,0,mfitdx:mfitddp)
        else
          levele=levele+1
          call qtwissfracgeo(vtwiss,gv,lxp,fr,.true.,over)
          cod=vtwiss(mfitdx:mfitddp)
          lv=itfdownlevel()
        endif
        ogv=tforbitgeo(gv,cod)
c        write(*,'(a,i5,1p12g10.2)')'ogeo ',lxp,fr,gv(:,4),ogv(:,4),cod(1:4)
        select case (keyword)
        case('OGX')
          kx=dfromr(ogv(1,4))
        case('OGY')
          kx=dfromr(ogv(2,4))
        case('OGZ')
          kx=dfromr(ogv(3,4))
        case('OCHI1')
          kx=dfromr(tfchi(ogv,1))
        case('OCHI2')
          kx=dfromr(tfchi(ogv,2))
        case('OCHI3')
          kx=dfromr(tfchi(ogv,3))
        end select
      case('DIR')
        if(ia /= nlat)then
          if(ref)then
            kx=dfromr(direlc(ia))
          else
            kx%k=ktfref+latt(ia)+1
            call compelc(ia,cmp)
            cmp%update=cmp%nparam <= 0
          endif
        else
          kx%k=ktftrue
        endif
      case('BZS')
        if(ref)then
          kx=dfromr(tfbzs(ia,ibz))
        else
          kx=dxzero
        endif
      case('BZT')
        if(ref)then
          kx=dfromr(tfbzt(ia,ibz))
        else
          kx=dxzero
        endif
      case('UPDATE')
        if(ia .lt. nlat)then
          call compelc(ia,cmp)
          kx%k=merge(ktftrue,ktffalse,cmp%update)
        else
          kx%k=ktftrue
        endif
      case('DK')
        kax=iferrk-2+ia*2
        if(ref)then
          kx=dfromr(rlist(kax))
        else
          kx%k=ktfref+kax
          call compelc(ia,cmp)
          cmp%update=cmp%nparam <= 0
        endif
      case('NDIV')
        kx%x(1)=dble(ndivelm(ia))
      case default
        if(keyword(1:3) == 'SIG' .and. keyword(4:4) /= 'E' .and. keyword(1:5) /= 'SIGMA')then
          call tfbeamkey(keyword(4:),i,j,irtc)
          if(irtc == 0)then
            if(ifsize == 0)then
              call tfsize(.true.)
            endif
            call tfbeamfrac(ia,v,0.d0,beam)
            if(i == 0)then
              if(j == 0)then
                kax=ktadaloc(-1,6)
                do i=1,6
                  kai=ktavaloc(0,6)
                  klist(kax+i)=ktflist+kai
                  rlist(kai+1:kai+6)=beam(iaidx(i,1:6))
                enddo
                kx%k=ktflist+kax
              else
                kx=dfromr(sqrt(beam(iaidx(j,j))))
              endif
            else
              kx=dfromr(merge(sqrt(beam(iaidx(i,i))),beam(iaidx(i,j)),
     $             j == 0))
            endif
            return
          endif
        endif
        kx=tflinek(keyword,ia,cmp,ref)
      end select
      return
      end

      subroutine tfbeamkey(key,i,j,irtc)
      implicit none
      integer*8 ,intent(out):: i,j
      integer*4 ,intent(out):: irtc
      integer*4 lk,l1,k1,k,lenw,itfmessage
      character*(*) ,intent(in):: key
      character*2 key1
      character*2 ,parameter ::keyname(6)=['X ','PX','Y ','PY','Z ','DP']
      lk=lenw(key)
      if(lk == 0)then
        i=0
        j=0
        irtc=0
        return
      endif
      l1=max(1,lk-1)
      key1=key(max(1,lk-1):lk)
      do k=1,6
        if(key1 == keyname(k))then
          k1=k
          go to 1
        endif
      enddo
 2    l1=lk
      key1=key(lk:lk)
      do k=1,5,2
        if(key1 == keyname(k))then
          k1=k
          go to 1
        endif
      enddo
      go to 9000
 1    if(l1 == 1)then
        i=0
        j=k1
      else
        key1=key(1:l1-1)
        do k=1,6
          if(key1 == keyname(k))then
            i=k
            j=k1
            go to 100
          endif
        enddo
        if(lk == l1+1)then
          go to 2
        endif
        go to 9000
      endif
 100  irtc=0
      return
 9000 irtc=itfmessage(9,'General::wrongval',
     $     '"coordinate","X, PX, Y, PY, Z, DP"')
      return
      end

      subroutine tflinestk(k,narg,isp0,irtc)
      use geto
      use ffs
      use tffitcode
      use ffs_pointer, only:idelc,pnamec,ielma
      implicit none
      type (sad_descriptor) ,intent(in):: k
      integer*4 ,intent(in):: narg
      integer*4 ,intent(out):: irtc,isp0
      integer*4 nc,itfmessage,i
      real*8 r,v
      character*32768 name
      isp0=isp
      if(ktfrealq(k,v) .and. narg == 2)then
        i=floor(v)
        r=v-i
        i=ielma(i)
        isp=isp+1
        itastk(1,isp)=ilist(i,ifele1)
        itastk(2,isp)=i
        vstk2(isp)=r
        irtc=0
      else
        if(narg == 1)then
          name(1:1)='*'
          nc=1
        else
          call tfgetstrns(k,name,nc)
          if(nc <= 0)then
            irtc=itfmessage(9,'General::wrongtype',
     $           '"name of component"')
            return
          endif
          if(convcase)then
            call capita1(name(1:nc))
          endif
        endif
        call tflinenamestk(name(1:nc),narg,0,0.d0,isp0,irtc)
        do i=isp0+1,isp
          itastk(2,i)=ielma(itastk(2,i))
        enddo
      endif
      return
      end

      recursive subroutine tflinenamestk(name0,narg,ioff,fr,isp0,irtc)
      use tfstk
      use ffs
      use tffitcode
      use ffs_pointer, only:idelc,pnamec,iele1
      implicit none
      type (sad_descriptor) kx
      integer*8 kav,ka,j,jj
      integer*4 ,intent(in):: narg,ioff
      integer*4 ,intent(out):: irtc
      integer*4 isp0,nc,ifany1,ielmf,itehash,l,i,ipoff,io
      real*8 ,intent(in):: fr
      real*8 r
      character*(*) ,intent(in):: name0
      character*32768 name,name2
      character*(MAXPNAME+16) name1
      logical*4 exist,temat
      integer*4 nl
      nc=len(name0)
      name(1:max(1,nc))=name0
      if(name(1:1) == '@')then
        name(1:nc-1)=name(2:nc)
        nc=nc-1
      else
        ipoff=ifany1(name(1:nc),nc,'+-',1)
        if(ipoff == 0)then
          call tfgetlineps(name,nc,nl,kav,0,irtc)
          if(irtc /= 0)then
            return
          endif
          if(nl > 0)then
            itastk(2,isp+1:isp+nl)=int(rlist(kav+1:kav+nl))+ioff
            itastk(1,isp+1:isp+nl)=iele1(itastk(2,isp+1:isp+nl))
            vstk2(isp+1:isp+nl)=fr
            isp=isp+nl
            return
          endif
        else
          call tfevals(name(ipoff:nc),kx%k,irtc)
          if(irtc /= 0 .or. ktfnonrealq(kx,r))then
            if(irtc > 0 .and. ierrorprint /= 0)then
              call tfreseterror
            endif
            return
          endif
          io=floor(r)
          call tflinenamestk(name(1:ipoff-1),narg,
     $         io,r-dble(io),isp0,irtc)          
          return
        endif
      endif
      if(nc > 2 .and. name(nc-1:nc) == '.*' .and.
     $     ifany1(name(1:nc),nc-2,'*%{|',1) == 0)then
        name2(1:nc-2)=name(1:nc-2)
        ka=itehash(name2(1:nc-2),nc-2)*2
        j=klist(ielmhash+ka+2)
        if(j /= 0)then
          do jj=j,j+ilist(1,ielmhash+ka+1)-1
            l=ilist(1,jj)
            if(name2(1:nc-2) == pnamec(l))then
              isp=isp+1
              itastk(1,isp)=ilist(l+ioff,ifele1)
              itastk(2,isp)=l+ioff
              vstk2(isp)=fr
            endif
          enddo
        endif
        irtc=0
      elseif(name(1:nc) /= '***' .and. name(1:nc) /= '^^^' .and.
     $       ifany1(name(1:nc),nc,'*%{|',1) > 0)then
        do i=1,nlat
          if(temat(i,name1,name(1:nc)))then
            isp=isp+1
            itastk(1,isp)=ilist(i+ioff,ifele1)
            itastk(2,isp)=i+ioff
            vstk2(isp)=fr
          endif
        enddo
        irtc=0
      else
        i=ielmf(name(1:nc),r,exist,0)
        if(exist)then
          isp=isp+1
          itastk(1,isp)=ilist(i+ioff,ifele1)
          itastk(2,isp)=i+ioff
          vstk2(isp)=fr
          irtc=0
        else
          irtc=0
        endif
      endif
      return
      end

      subroutine tfclearlinep()
      use tflinepcom
      implicit none
      iflinep=0
      return
      end

      subroutine tfinitlinep(irtc)
      use tflinepcom
      use tfstk
      use maccbk
      use efun
      implicit none
      type (sad_descriptor) kx
      integer*4 ,intent(out):: irtc
      integer*4 isp0
      irtc=0
      if(ifinitlinep == 0)then
        ifinitlinep=ktfsymbol+
     $       ktfsymbolz('InitLINE$P',10)
        ksdumm=kxsalocb(0,'-',1)
      endif
      if(iflinep == 0)then
        iflinep=ktfsymbol+ktfsymbolz('LINE$P',6)
        ifelementp=ktfsymbol+ktfsymbolz('Element$P',9)
        ifelementkeyp=ktfsymbol+ktfsymbolz('Element$Key',11)
        iftypekeyp=ktfsymbol+ktfsymbolz('Type$Key',8)
        isp0=isp
        isp=isp+1
        ktastk(isp)=ifinitlinep
        kx=tfefunref(isp0+1,.false.,irtc)
        isp=isp0
      endif
      return
      end

      subroutine tfgetlinep(ks,nl,kax,mode,irtc)
      use tflinepcom
      use tfstk
      use efun
      implicit none
      type (sad_descriptor) ,intent(in):: ks
      type (sad_descriptor) kx
      type (sad_rlist), pointer :: klr
      integer*8 ,intent(out):: kax
      integer*4 ,intent(out):: irtc
      integer*4 ,intent(in):: mode
      integer*4 isp0,nl
      call tfinitlinep(irtc)
      if(irtc /= 0)then
        nl=0
        return
      endif
      isp0=isp
      isp=isp+1
      ktastk(isp)=merge(iflinep,ifelementp,mode == 0)
      isp=isp+1
      dtastk(isp)=ks
      kx=tfefunref(isp0+1,.false.,irtc)
      isp=isp0
      if(.not. tfreallistq(kx,klr))then
        nl=0
      else
        kax=ktfaddr(kx)
        nl=klr%nl
      endif
      return
      end

      subroutine tfgetlineps(name0,lname,nl,kax,mode,irtc)
      use tfstk
      use ffs_flag
      use tmacro
      implicit none
      type (sad_descriptor) ks
      integer*8 ,intent(out):: kax
      integer*4 ,intent(in):: lname,mode
      integer*4 ,intent(out):: irtc
      integer*4 nl
      character*(*) ,intent(in):: name0
      character*(lname) name
      name=name0(1:lname)
      if(convcase)then
        call capita1(name)
      endif
      ks=kxsalocb(-1,name,lname)
      call tfgetlinep(ks,nl,kax,mode,irtc)
      return
      end

      integer*4 function itfloc(k,irtc)
      use tfstk
      use ffs
      use tffitcode
      implicit none
      type (sad_descriptor) ,intent(in):: k
      integer*4 ,intent(out):: irtc
      integer*4 nc,ielm,itfmessage,itfmessagestr
      character*(MAXPNAME+16) tfgetstrs,name
      logical*4 exist
      irtc=0
      if(ktfrealq(k,itfloc))then
        if(itfloc <= 0 .or. itfloc > nlat)then
          irtc=itfmessage(9,'General::wrongval',
     $         '"Component number",'//
     $         '"positive and less than length of beam line"')
          return
        endif
      else
        name=tfgetstrs(k%k,nc)
        if(nc <= 0)then
          irtc=itfmessage(9,'General::wrongtype',
     $         '"name of component"')
          itfloc=0
          return
        endif
        itfloc=ielm(name(1:nc),exist)
        if(.not. exist)then
          irtc=itfmessagestr(9,'FFS::undefcomp',name(1:nc))
        endif
      endif
      return
      end

      subroutine tfclrtparaed
      use tmacro
      implicit none
      tparaed=.false.
      return
      end

      integer*8 function ktfgeol(geo)
      use tfstk
      use geolib
      implicit none
      type (sad_dlist), pointer :: kl
      type (sad_rlist), pointer :: klv,klv2
      integer*8 kax,kax1,kax2
      real*8 ,intent(in):: geo(3,4)
      kax=ktadaloc(-1,2,kl)
      kax1=ktavaloc(0,3,klv)
      klv%rbody(1:3)=geo(:,4)
      kl%body(1)=ktflist+kax1
      kax2=ktavaloc(0,3,klv2)
      klv2%rbody(1)=tfchi(geo,1)
      klv2%rbody(2)=tfchi(geo,2)
      klv2%rbody(3)=tfchi(geo,3)
      kl%body(2)=ktflist+kax2
      ktfgeol=kax
      return
      end
