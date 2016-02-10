      subroutine tfsetv(latt,ivarele,ivvar,ivcomp,valvar,
     $     nvar,nele,
     $     klp,ival,couple,errk,iele,iele1,nlat)
      use tfstk
      implicit none
      include 'inc/MACCODE.inc'
      include 'inc/MACKW.inc'
      integer*4 nvar,nele,nlat
      integer*4 ivarele(nvar),ivvar(nvar),klp(nele),
     $     iele(nlat),iele1(nlat),latt(2,nlat),ival(nele),ivcomp(nvar)
      real*8 valvar(nvar),couple(nlat),errk(2,nlat)
      integer*4 i,ie,iv,j,ii,ie1
      if(nvar .gt. 0)then
        call tfclrtparaed
        do i=1,nlat-1
          ii=iele(i)
          ie=iele1(ii)
          ie1=iele1(i)
          do j=1,nvar
            iv=ivvar(j)
            if(iv .eq. ival(ie) .and. ivarele(j) .eq. ie
     $           .and. (ivcomp(j) .eq. 0 .or. ivcomp(j) .eq. ii))then
              rlist(latt(2,i)+iv)=valvar(j)*errk(1,i)*couple(i)
            elseif(iv .ne. ival(ie) .and. ivarele(j) .eq. ie1
     $             .and. (ivcomp(j) .eq. 0 .or. ivcomp(j) .eq. i))then
              rlist(latt(2,i)+iv)=valvar(j)
            endif
            if(ivarele(j) .gt. ie)then
              exit
            endif
          enddo
        enddo
        call tfinitvar(ivarele,ivvar,ivcomp,valvar,nvar,
     $       latt,klp,ival,errk,nlat,nele)
      endif
      return
      end

      subroutine tfinitvar(ivarele,ivvar,ivcomp,valvar,nvar,
     $     latt,klp,ival,errk,nlat,nele)
      use tfstk
      implicit none
      integer*4 nvar,nlat,nele,ivarele(nvar),ivvar(nvar),
     $     latt(2,nlat),klp(nele),ival(nele),ivcomp(nvar),k,
     $     irtc,icslfno
      real*8 valvar(nvar),errk(2,nlat)
      integer*4 i,ie,iv
      call tffscoupledvar(irtc)
      if(irtc .ne. 0)then
        call termes(icslfno(),'?Error in CoupledVariables',' ')
      endif
      do i=1,nvar
        ie=ivarele(i)
        iv=ivvar(i)
        k=ivcomp(i)
        if(k .eq. 0)then
          k=klp(ie)
        endif
        if(iv .eq. ival(ie))then
          valvar(i)=rlist(latt(2,k)+iv)/errk(1,k)
        else
          valvar(i)=rlist(latt(2,k)+iv)
        endif
      enddo
      return
      end

      subroutine tfsavevar(ie,itouchele,itouchv,ntouch,latt,klp,
     $     nlat,nele)
      use tfstk
      implicit none
      integer*4 ntouch,nlat,nele,itouchele(ntouch),itouchv(ntouch),
     $     latt(2,nlat),klp(nele),i,ie
      do i=1,ntouch
        if(itouchele(i) .eq. ie)then
          rlist(idval(latt(1,klp(ie)))+itouchv(i))=
     $         rlist(latt(2,klp(ie))+itouchv(i))
        endif
      enddo
      return
      end

      subroutine tffsadjust(itouchele,itouchv,latt,errk,
     $     couple,iele,iele1,klp,ival,ntouch)
      use tfstk
      use ffs
      use tffitcode
      implicit none
      integer*4 itv(ntouch+1),ite(nele+1)
      integer*4 ntouch,itouchele(ntouch),itouchv(ntouch),
     $     latt(2,nlat),klp(nele),i,ie,ival(nele),iele(nlat),j,
     $     iele1(nlat),irtc,icslfno,ie1,
     $     ntv,k
      real*8 errk(2,nlat),couple(nlat)
      call tfclrtparaed
      call tclr(ite,(nele+1)/2)
      ntv=0
      do j=1,ntouch
        if(ival(itouchele(j)) .ne. itouchv(j))then
          ntv=ntv+1
          itv(ntv)=j
          ite(itouchele(j))=1
        endif
      enddo
      if(ntv .eq. 0)then
        do i=1,nlat-1
          ie=iele1(iele(i))
          if(ie .ne. 0 .and. ival(ie) .ne. 0)then
            rlist(latt(2,i)+ival(ie))
     $           =rlist(latt(2,iele(i))+ival(ie))
     $           /errk(1,iele(i))*errk(1,i)*couple(i)
          endif
        enddo
      else
        do i=1,nlat-1
          ie=iele1(iele(i))
          ie1=iele1(i)
          if(ie .ne. 0 .and. ival(ie) .ne. 0)then
            rlist(latt(2,i)+ival(ie))
     $           =rlist(latt(2,iele(i))+ival(ie))
     $           /errk(1,iele(i))*errk(1,i)*couple(i)
          endif
          if(ite(ie1) .ne. 0)then
            do k=1,ntv
              j=itv(k)
              if(itouchele(j) .eq. ie1)then
                rlist(latt(2,i)+itouchv(j))
     $               =rlist(latt(2,klp(ie1))+itouchv(j))
              endif
            enddo
          endif
        enddo
      endif
      call tffscoupledvar(irtc)
      if(irtc .ne. 0)then
        call termes(icslfno(),'?Error in CoupledVariables',' ')
      endif
      return
      end

      subroutine tffscoupledvar(irtc)
      use tfstk
      implicit none
      integer*4 ia,irtc,icslfno,itfdownlevel,irtc1
      type (sad_list), pointer :: kl
      type (sad_descriptor) ifcoupv,ifsetcoup,k,kx
      type (sad_symdef), pointer, save :: symdcoupv
      data ifcoupv%k,ifsetcoup%k /0,0/
      if(ifcoupv%k .eq. 0)then
        ifcoupv  =kxsymbolz('`EVList',7)
        call descr_sad(ifcoupv,symdcoupv)
        ifsetcoup=kxsymbolz('`SetCoupledElements',19)
      endif
      irtc=0
      k=symdcoupv%value
      if(tflistqd(k,kl))then
        if(kl%nl .eq. 0)then
          return
        endif
      endif
      levele=levele+1
      call tfsyeval(ifsetcoup,kx,irtc)
      ia=itfdownlevel()
      if(irtc .ne. 0)then
        if(irtc .gt. 0)then
          if(ierrorprint .ne. 0)then
            ierrorf=0
            irtc1=irtc
            call tfaddmessage(
     $           'SetCoupledElements',18,icslfno())
          endif
        endif
        return
      endif
      return
      end
