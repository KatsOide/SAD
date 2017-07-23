      module tflinepcom
      use tfstk
      implicit none
      integer*8, save :: iflinep=0,ifinitlinep=0,ifelementp=0,
     $     ifelementkeyp=0,iftypekeyp=0
      type (sad_descriptor) ksdumm

      contains

      integer*4 function itftypekey(i,word1,lw1)
      use tfstk
      implicit none
      type (sad_descriptor) kx,ks
      integer*4 i,isp0,lw1,irtc
      character*(lw1) word1
      isp0=isp
      ktastk(isp0+1)=iftypekeyp
      rtastk(isp0+2)=dble(i)
      isp=isp0+3
      ks=kxsalocb(-1,word1,lw1)
      dtastk(isp)=ks
      levele=levele+1
      call tfefunref(isp0+1,kx,.false.,irtc)
      call tfconnect(kx,irtc)
      isp=isp0
      if(irtc .ne. 0)then
        call tfreseterror
        itftypekey=0
      elseif(.not. ktfrealq(kx))then
        itftypekey=0
      else
        itftypekey=ifromd(kx)
      endif        
      return
      end function

      integer*4 function itfelementkey(i,word1,lw1)
      use tfstk
      implicit none
      type (sad_descriptor) kx,ks
      integer*4 i,isp0,lw1,irtc
      character*(lw1) word1
      isp0=isp
      ktastk(isp0+1)=ifelementkeyp
      rtastk(isp0+2)=dble(i)
      isp=isp0+3
      ks=kxsalocb(-1,word1,lw1)
      dtastk(isp)=ks
      levele=levele+1
      call tfefunref(isp0+1,kx,.false.,irtc)
      call tfconnect(kx,irtc)
      isp=isp0
      if(irtc .ne. 0)then
        call tfreseterror
        itfelementkey=0
      elseif(.not. ktfrealq(kx))then
        itfelementkey=0
      else
        itfelementkey=ifromd(kx)
      endif        
      return
      end function

      subroutine tftypekey(isp1,kx,irtc)
      use tfstk
      use maccbk
      use maccode
      use mackw
      implicit none
      type (sad_descriptor) kx
      real*8 c
      integer*4 ic,isp0,isp1,irtc,itfmessage
      if(isp .ne. isp1+1)then
        irtc=itfmessage(9,'General::narg','"1"')
        return
      endif
      if(.not. ktfrealq(dtastk(isp),c))then
        irtc=itfmessage(9,'General::wrongtype',
     $       '"Typecode for #1"')
        return
      endif
      ic=int(c)
      if(ic .lt. 0 .or. ic .gt. icMXEL)then
        irtc=itfmessage(9,'General::wrongval',
     $       '"Typecode out of range"')
        return
      endif
      isp0=isp
      call tftypekeystk(ic,.true.)
      kx=kxmakelist(isp0)
      isp=isp0
      irtc=0
      return
      end subroutine

      subroutine tftypekeystk(ic,all)
      use tfstk
      use mackw
      implicit none
      type (sad_descriptor) kx
      integer*4 ic,lenw,i
      logical*4 all
      character*(MAXPNAME) key,tfkwrd,tfkwrd1
      do i=1,kytbl(kwMAX,ic)-1
        key=tfkwrd(ic,i)
        if(all .or. key .ne. '-')then
          isp=isp+1
          if(key .eq. '-')then
            dtastk(isp)=ksdumm
          else
            dtastk(isp)=kxsalocb(-1,key,lenw(key))
            if(kyindex1(i,ic) .ne. 0)then
              key=tfkwrd1(ic,i)
              isp=isp+1
              dtastk(isp)=kxsalocb(-1,key,lenw(key))
              kx=kxmakelist(isp-2)
              isp=isp-1
              dtastk(isp)=kx
            endif
          endif
        endif
      enddo
      return
      end subroutine

      type (sad_descriptor) function tfkeyv(i,keyword,ia,cmp,saved)
      use tfstk
      use ffs
      use tffitcode
      use ffs_pointer, only:idelc,elatt,idtypec,idvalc,sad_comp,
     $     compelc
      implicit none
      type (sad_rlist), pointer :: kld,klr
      integer*4 i,it,kl,l,j,lk,lenw
      integer*8 ia
      character*(*) keyword
      character*128 key
      logical*4 saved,sum
      real*8 s
      type (sad_comp), pointer :: cmp
c     begin initialize for preventing compiler warning
      kl=0
c     end   initialize for preventing compiler warning
      if(i .gt. 0)then
        call compelc(i,cmp)
      else
        kl=ilist(-i,ifklp)
        call compelc(kl,cmp)
      endif
      sum=.false.
      lk=lenw(keyword)
      if(lk .gt. 4 .and. keyword(lk-3:lk) .eq. '$SUM')then
        sum=.true.
        lk=lk-4
      endif
      key(1:lk)=keyword(1:lk)
      it=idtype(cmp%id)
      do l=1,kytbl(kwMAX,it)-1
        j=kyindex(l,it)
        if(j .gt. 0)then
          if(pname(kytbl(j,0))(2:) .eq. key(1:lk))then
            go to 1
          endif
          j=kyindex1(l,it)
          if(j .gt. 0)then
            if(pname(kytbl(j,0))(2:) .eq. key(1:lk))then
              go to 1
            endif
          endif
        endif
      enddo
      ia=0
      tfkeyv%k=0
      return
 1    if(i .gt. 0)then
        ia=elatt%comp(i)+l
        tfkeyv=dlist(ia)
      elseif(saved)then
        ia=idvalc(kl)+l
        tfkeyv=dlist(ia)
      else
        ia=elatt%comp(kl)+l
        if(l .eq. ilist(-i,ifival))then
          if(tfreallistq(dlist(ia),klr))then
            tfkeyv=kxavaloc(-1,klr%nl,kld)
            kld%rbody(1:klr%nl)=klr%rbody(1:klr%nl)
     $           /rlist(iferrk+(kl-1)*2)
c            write(*,*)'tfkeyv ',i,l,ia,l,kl,ilist(-i,ifival),
c     $           klr%nl,rlist(iferrk+(kl-1)*2)
c            call tfdebugprint(dlist(ia),'tfkeyv-s',1)
c            call tfdebugprint(tfkeyv,'tfkeyv-d',2)
c          tfkeyv=rlist(ia)/rlist(iferrk+(kl-1)*2)
          else
            tfkeyv=dlist(ia)
          endif
        else
          tfkeyv=dlist(ia)
        endif
      endif
      if(sum .and. tfreallistq(tfkeyv,klr))then
        s=klr%rbody(1)
        do i=2,klr%nl
          s=s+klr%rbody(i)
        enddo
        tfkeyv=dfromr(s)
      endif
      return
      end function

      end module

      subroutine tftwiss(isp1,kx,ref,irtc)
      use tfstk
      use ffs
      use tffitcode
      use ffs_fit, only:nlist
      implicit none
      type (sad_descriptor) kx
      type (sad_list), pointer :: klx
      integer*4 nkey
      parameter (nkey=mfitpzpy)
      integer*8 ktatwissaloc,kax,kaxi,itoff
      integer*4 isp1,irtc,narg,i,m,nc,j,
     $     isp0,nd,kt,itfmessage,lenw
      real*8 ftwiss(28),pe(4),tfgettwiss,tphysdisp,tphysdispz
      logical*4 over,ref
      character*(MAXPNAME+16) keyword,tfgetstrs
      narg=isp-isp1
      keyword=tfgetstrs(ktastk(isp1+1),nc)
      if(nc .le. 0)then
        irtc=itfmessage(9,'General::wrongtype',
     $       '"Keyword-string for #1"')
        return
      endif
      irtc=0
      kx%k=0
      call capita(keyword(1:nc))
      if(keyword .eq. 'LENGTH')then
        if(narg .gt. 1)then
          irtc=itfmessage(9,'General::narg','"1"')
          return
        endif
        kx=dfromr(dble(nlat))
      elseif(keyword .eq. '*' .or. keyword .eq. 'ALL')then
        if(narg .eq. 1)then
          kax=ktadaloc(-1,nlat,klx)
          do i=1,nlat
            kaxi=ktatwissaloc(0)
            klx%body(i)=ktflist+kaxi
            do j=1,ntwissfun
              itoff=((2*ndim+1)*(j-1)+ndim)*nlat+iftwis+i-1
              rlist(kaxi+i)=rlist(itoff)
            enddo
          enddo
        elseif(narg .eq. 2)then
          call tflinestk(dtastk(isp),narg,isp0,irtc)
          if(irtc .ne. 0)then
            return
          endif
          if(isp .eq. isp0+1)then
            kax=ktatwissaloc(-1)
            if(vstk2(isp) .eq. 0.d0)then
              do j=1,ntwissfun
                itoff=(2*ndim+1)*nlat*(j-1)+ndim*nlat+iftwis
     $               +itastk(2,isp)-1
                rlist(kax+j)=rlist(itoff)
              enddo
            else
              write(*,*)'tftwiss-1 '
              call qtwissfrac(rlist(kax+1),itastk(2,isp),
     $             vstk2(isp),over)
            endif
          else
            m=isp-isp0
            kax=ktadaloc(-1,m,klx)
            do i=1,m
              kaxi=ktatwissaloc(0)
              klx%body(i)=ktflist+kaxi
              if(vstk2(isp0+i) .eq. 0.d0)then
                do j=1,ntwissfun
                  itoff=(2*ndim+1)*nlat*(j-1)+ndim*nlat+iftwis
     $                 +itastk(2,isp0+i)-1
                  rlist(kaxi+j+1)=rlist(itoff)
                enddo
              else
                write(*,*)'tftwiss-2 '
                call qtwissfrac(rlist(kaxi+1),itastk(2,isp0+i),
     $               vstk2(isp0+i),over)
              endif
            enddo
          endif
          isp=isp0
        else
          irtc=itfmessage(9,'General::narg','"1 or 2"')
          return
        endif
        kx%k=ktflist+kax
      else
        do i=1,nkey
          if(keyword .eq. nlist(i))then
            kt=i
            go to 110
          endif
        enddo
        if(keyword(1:3) .eq. 'SIG' .or. keyword(1:4) .eq. 'SIZE'
     $       .or. keyword .eq. 'GAMMA'
     $       .or. keyword .eq. 'GAMMABETA'
     $       .or. keyword .eq. 'S')then
          call tfline(isp1,kx,ref,irtc)
        else
          irtc=itfmessage(9,'General::wrongval',
     $         '"#1 ('//keyword(1:lenw(keyword))//
     $         ') is","to be name of optical function"')
        endif
        return
 110    itoff=(2*ndim+1)*nlat*(kt-1)+ndim*nlat+iftwis
        if(narg .eq. 1)then
          kax=ktavaloc(-1,nlat)
          if(kt .le. ntwissfun)then
            if(ref)then
              rlist(kax+1:kax+nlat)=rlist(itoff:itoff+nlat-1)
            else
              do i=1,nlat
                klist(kax+i)=ktfref+itoff+i-1
              enddo
            endif
          elseif(kt .ge. mfitpex .and. kt .le. mfitpepy)then
            do i=1,nlat
              call tgetphysdisp(i,pe)
              rlist(kax+i)=pe(kt-mfitpex+1)
            enddo
          elseif(kt .ge. mfitpzx .and. kt .le. mfitpzpy)then
            do i=1,nlat
              call tgetphysdispz(i,pe)
              rlist(kax+i)=pe(kt-mfitpzx+1)
            enddo
          endif
          kx%k=ktflist+kax
        elseif(narg .eq. 2)then
          call tflinestk(dtastk(isp),narg,isp0,irtc)
          if(irtc .ne. 0)then
            return
          endif
          if(isp .eq. isp0+1)then
            if(vstk2(isp) .eq. 0.d0)then
              if(kt .le. ntwissfun)then
                if(ref)then
                  kx%k=klist(itoff+itastk(2,isp)-1)
                else
                  kx%k=ktfref+itoff+itastk(2,isp)-1
                endif
              elseif(kt .ge. mfitpex .and. kt .le. mfitpepy)then
                call tgetphysdisp(itastk(2,isp),pe)
                kx=dfromr(pe(kt-mfitpex+1))
              elseif(kt .ge. mfitpzx .and. kt .le. mfitpzpy)then
                call tgetphysdispz(itastk(2,isp),pe)
                kx=dfromr(pe(kt-mfitpzx+1))
              endif
c              call tfdebugprint(kx,'tftwiss',1)
            else
              call qtwissfrac(ftwiss,itastk(2,isp),
     $             vstk2(isp),over)
              kx=dfromr(tfgettwiss(kt,ftwiss))
c              write(*,*)'twiss ',rfromk(kx),kt,
c     $             itastk(2,isp),vstk2(isp)
            endif
          else
            m=isp-isp0
            kax=ktavaloc(-1,m)
            kx%k=ktflist+kax
            if(kt .le. ntwissfun)then
              if(ref)then
                do i=1,m
                  if(vstk2(isp0+i) .eq. 0.d0)then
                    rlist(kax+i)=rlist(itoff+itastk(2,isp0+i)-1)
                  else
                    write(*,*)'tftwiss-4 '
                    call qtwissfrac(ftwiss,itastk(2,isp0+i),
     $                   vstk2(isp0+i),over)
                    rlist(kax+i)=ftwiss(kt)
                  endif
                enddo
              else
                do i=1,m
                  if(vstk2(isp0+i) .eq. 0.d0)then
                    klist(kax+i)=ktfref+itoff+itastk(2,isp0+i)-1
                  else
                    write(*,*)'tftwiss-5 '
                    call qtwissfrac(ftwiss,itastk(2,isp0+i),
     $                   vstk2(isp0+i),over)
                    rlist(kax+i)=ftwiss(kt)
                  endif
                enddo
              endif
            elseif(kt .ge. mfitpex .and. kt. le. mfitpepy)then
              do i=1,m
                if(vstk2(isp0+i) .eq. 0.d0)then
                  call tgetphysdisp(itastk(2,isp0+i),pe)
                  rlist(kax+i)=pe(kt-mfitpex+1)
                else
                  write(*,*)'tftwiss-6 '
                  call qtwissfrac(ftwiss,itastk(2,isp0+i),
     $                 vstk2(isp0+i),over)
                  rlist(kax+i)=tphysdisp(kt,ftwiss)
                endif
              enddo
            elseif(kt .ge. mfitpzx .and. kt. le. mfitpzpy)then
              do i=1,m
                if(vstk2(isp0+i) .eq. 0.d0)then
                  call tgetphysdispz(itastk(2,isp0+i),pe)
                  rlist(kax+i)=pe(kt-mfitpzx+1)
                else
                  write(*,*)'tftwiss-7 '
                  call qtwissfrac(ftwiss,itastk(2,isp0+i),
     $                 vstk2(isp0+i),over)
                  rlist(kax+i)=tphysdispz(kt,ftwiss)
                endif
              enddo
            endif
          endif
          isp=isp0
        elseif(narg .eq. 3)then
          keyword=tfgetstrs(ktastk(isp-1),nc)
          if(nc .le. 0)then
            irtc=itfmessage(9,'General::wrongtype',
     $           '"name of component for #2"')
            return
          endif
          call capita(keyword(1:nc))
          if(kt .le. ntwissfun)then
            if(keyword .eq. 'SET')then
              kx=dtastk(isp)
              if(ktflistq(kx,klx))then
                nd=min(klx%nl,nlat)
                rlist(itoff:itoff+nd-1)=klx%rbody(1:nd)
                return
              endif
            endif
          endif
        endif
      endif
      return
      end

      integer*8 function ktatwissaloc(mode)
      use tfstk
      use ffs
      use tffitcode
      implicit none
      integer*8 kax
      integer*4 mode
      kax=ktavaloc(mode,28)
      rlist(kax+ntwissfun+1:kax+28)=0.d0
      rlist(kax+23)=1.d0
      ktatwissaloc=kax
      return
      end

      real*8 function tfgettwiss(kt,ftwiss)
      use tfstk
      use tffitcode
      implicit none
      integer*4 kt
      real*8 rfromk
      logical*4 isnan
      real*8 ftwiss(ntwissfun),tphysdisp,tphysdispz
      tfgettwiss=0.d0
      if(kt .le. ntwissfun)then
        tfgettwiss=ftwiss(kt)
      elseif(kt .ge. mfitpex .and. kt .le. mfitpepy)then
        tfgettwiss=tphysdisp(kt,ftwiss)
      elseif(kt .ge. mfitpzx .and. kt .le. mfitpzpy)then
        tfgettwiss=tphysdispz(kt,ftwiss)
      endif
      if(isnan(tfgettwiss))then
        tfgettwiss=rfromk(ktfnan)
      endif
      return
      end

      subroutine tfelement(isp1,kx,ref,irtc)
      use tfstk
      use ffs
      use tffitcode
      implicit none
      type (sad_descriptor) kx
      integer*4 isp1,irtc,i,narg,nc,isp0,m,itfmessage,ispa
      character*(MAXPNAME+16) keyword,tfgetstrs
      logical*4 saved,ref
      narg=isp-isp1
      irtc=0
      if(narg .le. 0)then
        irtc=itfmessage(9,'General::narg','"1, 2, or 3"')
        return
      endif
      keyword=tfgetstrs(ktastk(isp1+1),nc)
      if(nc .le. 0)then
        irtc=itfmessage(9,'General::wrongtype',
     $       '"List for #1"')
        return
      endif
      call capita(keyword(1:nc))
      if(keyword .eq. 'LENGTH')then
        if(narg .gt. 1)then
          go to 9010
        endif
        kx=dfromr(dble(nele))
      elseif(keyword .eq. 'EXPAND')then
        if(narg .ne. 1)then
          go to 9010
        endif
        call tffsadjust(flv%ntouch)
        kx%k=ktfoper+mtfnull
      else
        if(narg .gt. 2)then
          if(narg .eq. 3)then
            call tfgetoption('Saved',ktastk(isp),kx,irtc)
            if(irtc .ne. 0)then
              return
            endif
            if(ktfrealq(kx))then
              saved=kx%k .ne. 0
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
        if(irtc .ne. 0)then
          return
        endif
        if(isp .eq. isp0+1)then
          call tfelement1(itastk(1,isp),itastk(2,isp),
     $         kx,keyword,saved,ref,irtc)
        else
          m=isp-isp0
          ispa=isp
          do i=1,m
            isp=isp+1
            call tfelement1(itastk(1,isp0+i),itastk(2,isp0+i),
     $           dtastk(isp),keyword,saved,ref,irtc)
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
      use ffs_pointer, only:latt,idelc,idtypec,idvalc,sad_comp,
     $     compelc
      use tflinepcom
      implicit none
      type (sad_descriptor) kx
      type (sad_comp), pointer :: cmp
      integer*8 iax
      integer*4 irtc,id,lenw,it,ia,iv,isps,l,lpname
      character*(*) keyword
      character*(MAXPNAME) key,tfkwrd
      logical*4 saved,ref
      irtc=0
      if(keyword .eq. 'NAME')then
        id=idelc(ia)
        kx=kxsalocb(-1,pname(id),lpname(id))
      elseif(keyword .eq. 'VALUE')then
        iv=ilist(it,ifival)
        if(iv .gt. 0)then
          if(saved)then
            iax=idvalc(ia)+iv
            if(ref)then
              kx%k=klist(iax)
            else
              kx%k=ktfref+iax
            endif
          else
            iax=latt(ia)+iv
            if(ref)then
              kx=dfromr(rlist(iax)/rlist(iferrk+2*(ia-1)))
            else
              kx%k=ktfref+iax
              call compelc(ia,cmp)
              cmp%update=0
            endif
          endif
        else
          kx%k=0
        endif
      elseif(keyword .eq. 'DEFAULT')then
        iv=ilist(it,ifival)
        if(iv .eq. 0)then
          key=' '
        else
          key=tfkwrd(idtypec(ia),iv)
        endif
        Kx=kxsalocb(-1,key,lenw(key))
      elseif(keyword .eq. 'DEFAULT$SUM')then
        iv=ilist(it,ifival)
        if(iv .eq. 0)then
          key=' '
        else
          key=tfkwrd(idtypec(ia),iv)
          key=key(1:lenw(key))//"$SUM"
        endif
        Kx=kxsalocb(-1,key,lenw(key))
      elseif(keyword .eq. 'KEYWORDS' .or.
     $       keyword .eq. 'KEYWORDS_ALL')then
        l=0
        id=idtypec(ia)
        isps=isp
        call tftypekeystk(id,keyword .eq. 'KEYWORDS_ALL')
        kx=kxmakelist(isps)
        isp=isps
      elseif(keyword .eq. 'TYPE')then
        kx=dfromr(dble(idtypec(ia)))
      elseif(keyword .eq. 'TYPENAME')then
        key=pname(kytbl(0,idtypec(ia)))
        kx=kxsalocb(-1,key(2:),lenw(key)-1)
      elseif(keyword .eq. 'POSITION')then
        kx=dfromr(dble(it))
      else
        kx=tfkeyv(-it,keyword,iax,cmp,saved)
        if(.not. ref)then
          kx%k=ktfref+iax
          if(.not. saved)then
            cmp%update=0
          endif
        endif
      endif
      return
      end

      subroutine tfelementstk(k,isp0,narg,irtc)
      use tfstk
      use ffs
      use tffitcode
      use ffs_pointer, only:idelc,pnamec
      implicit none
      type (sad_descriptor) k
      integer*4 isp0,narg,irtc,iv,nc,ifany1,i,itfmessage,j,ielmh
      character*1024 name
      logical*4 tmatch
      isp0=isp
      if(ktfrealq(k,iv) .and. narg .eq. 2)then
        if(iv .lt. 0)then
          iv=nele+iv+1
        endif
        if(iv .le. 0 .or. iv .gt. nele)then
          irtc=itfmessage(9,'General::wrongnum',
     $         '"positive and less than length of beam line"')
          return
        endif
        isp=isp+1
        itastk(1,isp)=iv
        itastk(2,isp)=ilist(iv,ifklp)
        irtc=0
      else
        if(narg .eq. 1)then
          name='*'
          nc=1
        else
          call tfgetstrns(k,name,nc)
          if(nc .le. 0)then
            irtc=itfmessage(9,'General::wrongtype',
     $           '"Name of component"')
            return
          endif
          if(convcase)then
            call capita1(name(1:nc))
          endif
        endif
        irtc=0
        if(name(1:nc) .ne. '***' .and.
     $       ifany1(name(1:nc),nc,'*%{<|',1) .gt. 0)then
          do i=1,nele
c            write(*,*)'elementstk',i,nele,pname(idelc(ilist(i,ifklp)))
            if(tmatch(pnamec(ilist(i,ifklp)),
     $           name(1:nc)))then
              isp=isp+1
              itastk(1,isp)=i
              itastk(2,isp)=ilist(i,ifklp)
            endif
          enddo
        else
          j=ielmh(name(1:nc),0)
          if(j .ne. 0)then
            i=ilist(j,ifele1)
            isp=isp+1
            itastk(1,isp)=i
            itastk(2,isp)=ilist(i,ifklp)
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
      use ffs_pointer, only:latt
      implicit none
      type (sad_descriptor) kx
      integer*4 isp1,irtc,ie,iv,k,j,m,ispa
      integer*4 i,narg,nc,isp0,itfmessage
      character*(MAXPNAME+16) keyword,tfgetstrs
      logical*4 ref
      narg=isp-isp1
      keyword=tfgetstrs(ktastk(isp1+1),nc)
      if(nc .le. 0)then
        irtc=itfmessage(9,'General::wrongtype',
     $       '"Keyword for #1"')
        return
      endif
      irtc=0
      call capita(keyword(1:nc))
      if(keyword .eq. 'LENGTH')then
        if(narg .gt. 1)then
          irtc=itfmessage(9,'General::narg','"1"')
        endif
        kx=dfromr(dble(nlat))
      elseif(keyword .eq. 'EXPAND')then
        if(narg .ne. 1)then
          irtc=itfmessage(9,'General::narg','"1"')
        endif
        do i=1,nlat-1
          ie=ilist(ilist(i,ifele),ifele1)
          if(ie .gt. 0)then
            iv=ilist(ie,ifival)
            if(iv .gt. 0)then
              k=ilist(ie,ifklp)
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
        if(irtc .ne. 0)then
          isp=isp0
          return
        endif
        m=isp-isp0
        if(m .eq. 1)then
          call tfline1(isp,kx,keyword,ref,irtc)
        else
          ispa=isp
          do j=1,m
            isp=isp+1
            call tfline1(isp0+j,dtastk(isp),keyword,ref,irtc)
            if(irtc .ne. 0)then
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
      use tfstk
      use ffs
      use tffitcode
      use ffs_pointer, only:latt,idelc,idtypec,pnamec,compelc,direlc
      use sad_main
      use tflinepcom
      implicit none
      type (sad_descriptor) kx
      type (sad_comp), pointer ::cmp
      integer*8 kax,ktfgeol,kai,j,i,ip,itoff
      integer*4 irtc,lenw,lxp,nv,lt,kk,ia,isp1,ibz
      real*8 v,beam(42),xp,fr,geo1(12),pos0,
     $     gam0,gv(3,4),ogv(3,4),cod(4),vtwiss(27),tfchi,tfbzs
      type (sad_descriptor) dsave(kwMAX)
      character*(*) keyword
      character*64 name
      integer*4 iaa,lv,itfdownlevel
      integer*8 n,m
      logical*4 chg,over,ref
      iaa(m,n)=int(((m+n+abs(m-n))**2+2*(m+n)-6*abs(m-n))/8)
      irtc=0
      ip=itastk(1,isp1)
      ia=itastk(2,isp1)
      v=vstk2(isp1)
      if(keyword .eq. 'S' .or. keyword .eq. 'LENG')then
        kx=dfromr(rlist(ifpos+ia-1)*(1.d0-v)+
     $       rlist(ifpos+min(nlat-1,ia))*v)
      elseif(keyword .eq. 'GAMMABETA')then
        kx=dfromr(rlist(ifgamm+ia-1)*(1.d0-v)+
     $       rlist(ifgamm+min(nlat-1,ia))*v)
      elseif(keyword .eq. 'GAMMA')then
        kx=dfromr(sqrt(1.d0+(rlist(ifgamm+ia-1)*(1.d0-v)+
     $       rlist(ifgamm+min(nlat-1,ia))*v)**2))
      elseif(keyword(1:3) .eq. 'SIG' .and. keyword(1:5) .ne. 'SIGMA'
     $       .or. keyword(1:4) .eq. 'SIZE')then
        if(keyword(1:3) .eq. 'SIG')then
          call tfbeamfrac(ia,v,0.d0,beam)
          call tfbeamkey(keyword(4:),i,j,irtc)
        else
          call tfsizefrac(ia,v,beam)
          call tfbeamkey(keyword(5:),i,j,irtc)
        endif
        if(irtc .ne. 0)then
          return
        endif
        if(i .eq. 0)then
          if(j .eq. 0)then
            kax=ktadaloc(-1,6)
            do i=1,6
              kai=ktavaloc(0,6)
              klist(kax+i)=ktflist+kai
              do j=1,6
                rlist(kai+j)=beam(iaa(i,j))
              enddo
            enddo
            kx%k=ktflist+kax
          else
            kx=dfromr(sqrt(beam(iaa(j,j))))
          endif
        else
          if(j .eq. 0)then
            kx=dfromr(sqrt(beam(iaa(i,i))))
          else
            kx=dfromr(beam(iaa(i,j)))
          endif
        endif
      elseif(keyword .eq. 'MULT')then
        kx=dfromr(dble(ilist(ia,ifmult)))
      elseif(keyword .eq. 'TYPE')then
        if(ia .eq. nlat)then
          kx%k=0
        else
          kx=dfromr(dble(idtypec(ia)))
        endif
      elseif(keyword .eq. 'TYPENAME')then
        if(ia .eq. nlat)then
          kx=dxnulls
        else
          name=pname(kytbl(0,idtypec(ia)))
          kx=kxsalocb(-1,name(2:),lenw(name)-1)
        endif
      elseif(keyword .eq. 'NAME')then
        call elname(ia,name)
        kx=kxsalocb(-1,name,lenw(name))
      elseif(keyword .eq. 'ELEMENT')then
        if(ia .eq. nlat)then
          name='$$$'
        else
          name=pnamec(ia)
        endif
        kx=kxsalocb(-1,name,lenw(name))
      elseif(keyword .eq. 'POSITION')then
        kx=dfromr(dble(ia))
      elseif(keyword .eq. 'GEO')then
        xp=v+ia
        lxp=int(xp)
        fr=xp-lxp
        j=ifgeo+(lxp-1)*12
        if(fr .eq. 0.d0)then
          kax=ktfgeol(rlist(j))
        else
          lt=idtypec(lxp)
          call compelc(lxp,cmp)
          nv=kytbl(kwMAX,lt)-1
          dsave(1:nv)=cmp%dvalue(1:nv)
c          call tmov(rlist(latt(lxp)+1),dsave,nv)
          levele=levele+1
          call qfracseg(cmp,0.d0,fr,dsave,chg)
c          call qfraccomp(lxp,0.d0,fr,.false.,chg)
          if(chg)then
            geo1=rlist(j+1:j+12)
            pos0=rlist(ifpos+lxp)
            gam0=rlist(ifgamm+lxp)
            call tfgeo1(lxp,lxp+1,.true.,.false.)
            kax=ktfgeol(rlist(j+12))
            rlist(ifgamm+lxp)=gam0
            rlist(j+1:j+12)=geo1
            rlist(ifpos+lxp)=pos0
          else
            kax=ktfgeol(rlist(j+12))
          endif
          cmp%dvalue(1:nv)=dsave(1:nv)
          lv=itfdownlevel()
c          call tmov(vsave,rlist(latt(lxp)+1),nv)
        endif
        kx%k=ktflist+kax
      elseif(keyword .eq. 'GX' .or. keyword .eq. 'GY' .or.
     $       keyword .eq. 'GZ' .or. keyword .eq. 'GCHI1' .or.
     $       keyword .eq. 'GCHI2' .or. keyword .eq. 'GCHI3')then
        xp=v+ia
        lxp=int(xp)
        fr=xp-lxp
        j=ifgeo+(lxp-1)*12
        if(fr .eq. 0.d0)then
          call tmov(rlist(j),gv,12)
        else
          lt=idtypec(lxp)
          call compelc(lxp,cmp)
          nv=kytbl(kwmax,lt)-1
          dsave(1:nv)=cmp%dvalue(1:nv)
c          call tmov(rlist(latt(lxp)+1),dsave,nv)
          levele=levele+1
          call qfracseg(cmp,0.d0,fr,dsave,chg)
c          call qfraccomp(lxp,0.d0,fr,.false.,chg)
          if(chg)then
            call tmov(rlist(j+12),geo1,12)
            pos0=rlist(ifpos+lxp)
            gam0=rlist(ifgamm+lxp)
c            write(*,*)'tftwiss-gx ',lxp,v,ia
            call tfgeo1(lxp,lxp+1,.true.,.false.)
            call tmov(rlist(j+12),gv,12)
            rlist(ifgamm+lxp)=gam0
            call tmov(geo1,rlist(j+12),12)
            rlist(ifpos+lxp)=pos0
          else
            call tmov(rlist(j+12),gv,12)
          endif
          cmp%dvalue(1:nv)=dsave(1:nv)
          lv=itfdownlevel()
c          call tmov(vsave,rlist(latt(lxp)+1),nv)
        endif
        if(keyword .eq. 'GX')then
          kx=dfromr(gv(1,4))
        elseif(keyword .eq. 'GY')then
          kx=dfromr(gv(2,4))
        elseif(keyword .eq. 'GZ')then
          kx=dfromr(gv(3,4))
        elseif(keyword .eq. 'GCHI1')then
          kx=dfromr(tfchi(gv,1))
        elseif(keyword .eq. 'GCHI2')then
          kx=dfromr(tfchi(gv,2))
        elseif(keyword .eq. 'GCHI3')then
          kx=dfromr(tfchi(gv,3))
        endif
      elseif(keyword .eq. 'OGEO')then
        xp=v+ia
        lxp=int(xp)
        fr=xp-lxp
        j=ifgeo+(lxp-1)*12
        if(fr .eq. 0.d0)then
          do kk=mfitdx,mfitdpy
            itoff=(2*ndim+1)*nlat*(kk-1)+ndim*nlat+iftwis+lxp-1
            cod(kk-mfitdx+1)=rlist(itoff)
          enddo
          call tmov(rlist(j),gv,12)
        else
          lt=idtypec(lxp)
          call compelc(lxp,cmp)
          nv=kytbl(kwmax,lt)-1
          dsave(1:nv)=cmp%dvalue(1:nv)
c          call tmov(rlist(latt(lxp)+1),vsave,nv)
          levele=levele+1
          call qfracseg(cmp,0.d0,fr,dsave,chg)
c          call qfraccomp(lxp,0.d0,fr,.false.,chg)
          if(chg)then
            call tmov(rlist(j+12),geo1,12)
            pos0=rlist(ifpos+lxp)
            gam0=rlist(ifgamm+lxp)
            call tfgeo1(lxp,lxp+1,.true.,.false.)
            call tmov(rlist(j+12),gv,12)
            rlist(ifgamm+lxp)=gam0
            call tmov(geo1,rlist(j+12),12)
            rlist(ifpos+lxp)=pos0
          else
            call tmov(rlist(j+12),gv,12)
          endif
          cmp%dvalue(1:nv)=dsave(1:nv)
c          call tmov(vsave,rlist(latt(lxp)+1),nv)
          call qtwissfrac(vtwiss,lxp,fr,over)
          cod(1:4)=vtwiss(mfitdx:mfitdpy)
c          call tmov(vtwiss(mfitdx),cod,4)
          lv=itfdownlevel()
        endif
        call tforbitgeo(ogv,gv,cod(1),cod(2),cod(3),cod(4))
        kx%k=ktflist+ktfgeol(ogv)
      elseif(keyword .eq. 'OGX' .or. keyword .eq. 'OGY' .or.
     $       keyword .eq. 'OGZ' .or. keyword .eq. 'OCHI1' .or.
     $       keyword .eq. 'OCHI2' .or. keyword .eq. 'OCHI3')then
        xp=v+ia
        lxp=int(xp)
        fr=xp-lxp
        j=ifgeo+(lxp-1)*12
        if(fr .eq. 0.d0)then
          call tmov(rlist(j),gv,12)
          do kk=mfitdx,mfitdpy
            itoff=(2*ndim+1)*nlat*(kk-1)+ndim*nlat+iftwis
     $           +lxp-1
            cod(kk-mfitdx+1)=rlist(itoff)
          enddo
        else
          lt=idtypec(lxp)
          call compelc(lxp,cmp)
          nv=kytbl(kwmax,lt)-1
          dsave(1:nv)=cmp%dvalue(1:nv)
c          call tmov(rlist(latt(lxp)+1),vsave,nv)
          levele=levele+1
          call qfracseg(cmp,0.d0,fr,dsave,chg)
c          call qfraccomp(lxp,0.d0,fr,.false.,chg)
          if(chg)then
            call tmov(rlist(j+12),geo1,12)
            pos0=rlist(ifpos+lxp)
            gam0=rlist(ifgamm+lxp)
            call tfgeo1(lxp,lxp+1,.true.,.false.)
            call tmov(rlist(j+12),gv,12)
            rlist(ifgamm+lxp)=gam0
            call tmov(geo1,rlist(j+12),12)
            rlist(ifpos+lxp)=pos0
          else
            call tmov(rlist(j+12),gv,12)
          endif
          cmp%dvalue(1:nv)=dsave(1:nv)
c          call tmov(vsave,rlist(latt(lxp)+1),nv)
          call qtwissfrac(vtwiss,lxp,fr,over)
          cod(1:4)=vtwiss(mfitdx:mfitdpy)
c          call tmov(vtwiss(mfitdx),cod,4)
          lv=itfdownlevel()
        endif
        call tforbitgeo(ogv,gv,cod(1),cod(2),cod(3),cod(4))
        if(keyword .eq. 'OGX')then
          kx=dfromr(ogv(1,4))
        elseif(keyword .eq. 'OGY')then
          kx=dfromr(ogv(2,4))
        elseif(keyword .eq. 'OGZ')then
          kx=dfromr(ogv(3,4))
        elseif(keyword .eq. 'OCHI1')then
          kx=dfromr(tfchi(ogv,1))
        elseif(keyword .eq. 'OCHI2')then
          kx=dfromr(tfchi(ogv,2))
        elseif(keyword .eq. 'OCHI3')then
          kx=dfromr(tfchi(ogv,3))
        endif
      elseif(keyword .eq. 'DIR')then
        if(ia .ne. nlat)then
          if(ref)then
            kx=dfromr(direlc(ia))
          else
            kx%k=ktfref+latt(ia)+1
            call compelc(ia,cmp)
            cmp%update=0
          endif
        else
          kx%k=ktftrue
        endif
      elseif(keyword .eq. 'DK')then
        kax=iferrk-2+ia*2
        if(ref)then
          kx=dfromr(rlist(kax))
        else
          kx%k=ktfref+kax
          call compelc(ia,cmp)
          cmp%update=0
        endif
      elseif(keyword .eq. 'BZS')then
        if(ref)then
          kx=dfromr(tfbzs(ia,ibz))
        else
          kx%k=0
        endif
      else
        if(ia .lt. nlat)then
          kx=tfkeyv(int(ia),keyword,ip,cmp,.false.)
          if(.not. ref)then
            cmp%update=0
            kx%k=ktfref+ip
          endif
          tparaed=.false.
        else
          kx%k=0
        endif
      endif
      return
      end

      subroutine tfbeamkey(key,i,j,irtc)
      implicit none
      integer*8 i,j
      integer*4 irtc,lk,l1,k1,k,lenw,itfmessage
      character*(*) key
      character*2 key1,keyname(6)
      data keyname /'X ','PX','Y ','PY','Z ','DP'/
      lk=lenw(key)
      if(lk .eq. 0)then
        i=0
        j=0
        irtc=0
        return
      endif
      l1=max(1,lk-1)
      key1=key(max(1,lk-1):lk)
      do k=1,6
        if(key1 .eq. keyname(k))then
          k1=k
          go to 1
        endif
      enddo
 2    l1=lk
      key1=key(lk:lk)
      do k=1,5,2
        if(key1 .eq. keyname(k))then
          k1=k
          go to 1
        endif
      enddo
      go to 9000
 1    if(l1 .eq. 1)then
        i=0
        j=k1
      else
        key1=key(1:l1-1)
        do k=1,6
          if(key1 .eq. keyname(k))then
            i=k
            j=k1
            go to 100
          endif
        enddo
        if(lk .eq. l1+1)then
          go to 2
        endif
        go to 9000
      endif
 100  irtc=0
      return
 9000 irtc=itfmessage(9,'General::wrongval',
     $     '"coordinae","X, PX, Y, PY, Z, DP"')
      return
      end

      subroutine tflinestk(k,narg,isp0,irtc)
      use tfstk
      use ffs
      use tffitcode
      use ffs_pointer, only:idelc,pnamec
      implicit none
      type (sad_descriptor) k
      integer*8 kav,ka,j,jj
      integer*4 narg,isp0,irtc,nc,ifany1,ielmf,itfmessage,
     $     itehash,l,i
      real*8 r
      character*1024 name,name2
      character*(MAXPNAME+16) name1
      logical*4 exist,temat
      integer*4 nl
      isp0=isp
      if(ktfrealq(k) .and. narg .eq. 2)then
        i=int(rtastk(isp))
        if(i .lt. 0)then
          r=1.d0+rtastk(isp)-i
          if(r .ne. 1.d0)then
            i=i-1
          endif
          i=nlat+i+1
        else
          r=rtastk(isp)-i
        endif
        if(i .le. 0 .or. i .gt. nlat)then
          irtc=itfmessage(9,'General::wrongnum',
     $         '"positive and less than length of beam line"')
        else
          if(i .eq. nlat)then
            r=0.d0
          endif
          isp=isp+1
          itastk(1,isp)=ilist(i,ifele1)
          itastk(2,isp)=i
          vstk2(isp)=r
          irtc=0
        endif
      else
        if(narg .eq. 1)then
          name(1:1)='*'
          nc=1
        else
          call tfgetstrns(k,name,nc)
          if(nc .le. 0)then
            irtc=itfmessage(9,'General::wrongtype',
     $           '"name of component"')
            return
          endif
          if(convcase)then
            call capita1(name(1:nc))
          endif
        endif
        if(name(1:1) .eq. '@')then
          name(1:nc-1)=name(2:nc)
          nc=nc-1
        elseif(ifany1(name(1:nc),nc,'+-',1) .eq. 0)then
          call tfgetlineps(name,nc,nl,kav,0,irtc)
          if(irtc .ne. 0)then
            return
          endif
          if(nl .gt. 0)then
            do i=1,nl
              l=int(rlist(kav+i))
              isp=isp+1
              itastk(1,isp)=ilist(l,ifele1)
              itastk(2,isp)=l
              vstk2(isp)=0.d0
            enddo
            return
          endif
        endif
        if(nc .gt. 2 .and. name(nc-1:nc) .eq. '.*' .and.
     $       ifany1(name(1:nc),nc-2,'*%{|',1) .eq. 0)then
          name2(1:nc-2)=name(1:nc-2)
          ka=itehash(name2(1:nc-2),nc-2)*2
          j=klist(ielmhash+ka+2)
          if(j .ne. 0)then
            do jj=j,j+ilist(1,ielmhash+ka+1)-1
              l=ilist(1,jj)
              if(name2(1:nc-2) .eq. pnamec(l))then
                isp=isp+1
                itastk(1,isp)=ilist(l,ifele1)
                itastk(2,isp)=l
                vstk2(isp)=0.d0
              endif
            enddo
          endif
          irtc=0
        elseif(name(1:nc) .ne. '***' .and. name(1:nc) .ne. '^^^' .and.
     $       ifany1(name(1:nc),nc,'*%{|',1) .gt. 0)then
          do i=1,nlat
            if(temat(i,name1,name(1:nc)))then
              isp=isp+1
              itastk(1,isp)=ilist(i,ifele1)
              itastk(2,isp)=i
              vstk2(isp)=0.d0
            endif
          enddo
          irtc=0
        else
          i=ielmf(name(1:nc),r,exist)
          if(exist)then
            isp=isp+1
            itastk(1,isp)=ilist(i,ifele1)
            itastk(2,isp)=i
            vstk2(isp)=r
c            write(*,*)'linestk ',name(1:nc),r
            irtc=0
          else
            irtc=0
          endif
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
      implicit none
      type (sad_descriptor) kx
      integer*4 irtc,isp0
      irtc=0
      if(ifinitlinep .eq. 0)then
        ifinitlinep=ktfsymbol+
     $       ktfsymbolz('InitLINE$P',10)
        ksdumm=kxsalocb(0,'-',1)
      endif
      if(iflinep .eq. 0)then
        iflinep=ktfsymbol+ktfsymbolz('LINE$P',6)
        ifelementp=ktfsymbol+ktfsymbolz('Element$P',9)
        ifelementkeyp=ktfsymbol+ktfsymbolz('Element$Key',11)
        iftypekeyp=ktfsymbol+ktfsymbolz('Type$Key',8)
        isp0=isp
        isp=isp+1
        ktastk(isp)=ifinitlinep
        call tfefunref(isp0+1,kx,.false.,irtc)
        isp=isp0
      endif
      return
      end

      subroutine tfgetlinep(ks,nl,kax,mode,irtc)
      use tflinepcom
      use tfstk
      implicit none
      type (sad_descriptor) kx,ks
      type (sad_rlist), pointer :: klr
      integer*8 kax
      integer*4 irtc,isp0,nl,mode
      call tfinitlinep(irtc)
      if(irtc .ne. 0)then
        nl=0
        return
      endif
      isp0=isp
      isp=isp+1
      if(mode .eq. 0)then
        ktastk(isp)=iflinep
      else
        ktastk(isp)=ifelementp
      endif
      isp=isp+1
      dtastk(isp)=ks
      call tfefunref(isp0+1,kx,.false.,irtc)
      isp=isp0
      if(.not. tfreallistqd(kx,klr))then
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
      integer*8 kax
      integer*4 irtc,lname,nl,mode
      character*(*) name0
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
      type (sad_descriptor) k
      integer*4 irtc,nc,ielm,itfmessage
      character*(MAXPNAME+16) tfgetstrs,name
      logical*4 exist
      irtc=0
      if(ktfrealq(k,itfloc))then
        if(itfloc .le. 0 .or. itfloc .gt. nlat)then
          irtc=itfmessage(9,'General::wrongval',
     $         '"Component number",'//
     $         '"positive and less than length of beam line"')
          return
        endif
      else
        name=tfgetstrs(k,nc)
        if(nc .le. 0)then
          irtc=itfmessage(9,'General::wrongtype',
     $         '"name of component"')
          itfloc=0
          return
        endif
        itfloc=ielm(name(1:nc),exist)
        if(.not. exist)then
          irtc=itfmessage(9,'FFS::undefcomp',
     $         '"'//name(1:nc)//'"')
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
      implicit none
      integer*8 kax,kax1,kax2
      real*8 geo(3,4),tfchi
      kax=ktadaloc(-1,2)
      kax1=ktavaloc(0,3)
      rlist(kax1+1:kax1+3)=geo(1:3,4)
c      call tmov(geo(1,4),rlist(kax1+1),3)
      klist(kax+1)=ktflist+kax1
      kax2=ktavaloc(0,3)
      rlist(kax2+1)=tfchi(geo,1)
      rlist(kax2+2)=tfchi(geo,2)
      rlist(kax2+3)=tfchi(geo,3)
      klist(kax+2)=ktflist+kax2
      ktfgeol=kax
      return
      end
