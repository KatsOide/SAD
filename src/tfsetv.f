      subroutine tfsetv(nvar)
      use tfstk
      use ffs, only:nlat,nvevx
      use ffs_pointer
      use mackw
      use ffs_seg
      implicit none
      type (sad_comp), pointer :: cmp
      integer*4 nvar
      integer*4 i,ie,iv,j,ii,ie1
      if(nvar .gt. 0)then
        do i=1,nlat-1
          call compelc(i,cmp)
          ii=icomp(i)
          ie=iele1(ii)
          ie1=iele1(i)
          do j=1,nvar
            iv=nvevx(j)%ivvar
            if(iv .eq. ival(ie) .and. nvevx(j)%ivarele .eq. ie
     $           .and. (nvevx(j)%ivcomp .eq. 0 .or.
     $           nvevx(j)%ivcomp .eq. ii))then
              cmp%value(iv)=nvevx(j)%valvar*errk(1,i)*couple(i)
              cmp%update=iand(cmp%update,2)
            elseif(iv .ne. 0 .and. iv .ne. ival(ie) .and.
     $             nvevx(j)%ivarele .eq. ie1
     $             .and. (nvevx(j)%ivcomp .eq. 0 .or.
     $             nvevx(j)%ivcomp .eq. i))then
              cmp%value(iv)=nvevx(j)%valvar
              cmp%update=iand(cmp%update,2)
            endif
            if(nvevx(j)%ivarele .gt. ie)then
              exit
            endif
          enddo
        enddo
        call tfinitvar
      endif
      return
      end

      subroutine tffsadjustvar
      use tfstk
      use ffs
      use tffitcode
      implicit none
      call tffsadjust
      call tfinitvar
      return
      end

      subroutine tfinitvar
      use tfstk
      use ffs, only: flv,nvevx
      use ffs_pointer
      use tfcsi, only:icslfno
      implicit none
      integer*4 k,irtc
      integer*4 i,ie,iv
      call tffscoupledvar(irtc)
      if(irtc .ne. 0)then
        call termes(icslfno(),'?Error in CoupledVariables',' ')
      endif
      do i=1,flv%nvar
        ie=nvevx(i)%ivarele
        iv=nvevx(i)%ivvar
        k=nvevx(i)%ivcomp
        if(k .eq. 0)then
          k=klp(ie)
        endif
        if(iv .eq. ival(ie))then
          nvevx(i)%valvar=tfvalvar(k,iv)/errk(1,k)
        else
          nvevx(i)%valvar=tfvalvar(k,iv)
        endif
      enddo
      return
      end

      subroutine tffsadjust
      use tfstk
      use ffs
      use ffs_pointer
      use tffitcode
      use tfcsi, only:icslfno
      use ffs_seg
      implicit none
      integer*4 itv(flv%ntouch),ite(nve*2)
      integer*4 i,ie,j,irtc,ie1,ntv,k
      ite=0
      ntv=0
      do j=1,flv%ntouch
        if(ival(nvevx(j)%itouchele) .ne. nvevx(j)%itouchv)then
          ntv=ntv+1
          itv(ntv)=j
          ite(nvevx(j)%itouchele)=1
        endif
      enddo
      if(ntv .eq. 0)then
        do i=1,nlat-1
          ie=iele1(icomp(i))
          if(icomp(i) .ne. i .and. ie .ne. 0 .and. ival(ie) .ne. 0)then
            call tfvcopy(icomp(i),i,ival(ie),
     $           errk(1,i)*couple(i)/errk(1,icomp(i)))
          endif
        enddo
      else
        do i=1,nlat-1
          ie=iele1(icomp(i))
          ie1=iele1(i)
          if(ie .ne. 0 .and. i .ne. icomp(i) .and. ival(ie) .ne. 0)then
            call tfvcopy(icomp(i),i,ival(ie),
     $           errk(1,i)*couple(i)/errk(1,icomp(i)))
          endif
          if(ite(ie1) .ne. 0)then
            do k=1,ntv
              j=itv(k)
              if(nvevx(j)%itouchele .eq. ie1 .and. klp(ie1) .ne. i)then
                call tfvcopy(klp(ie1),i,nvevx(j)%itouchv,1.d0)
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

      subroutine tffsadjust1(isp0,var,nv)
      use tfstk
      use ffs
      use ffs_pointer
      use tflinepcom
      use tfcsi
      use ffs_seg
      implicit none
      type (sad_rlist), pointer ::kl
      integer*4 i,j,isp0,isp1,l,irtc,k,k1,ie1
      logical*4 var,nv
      if(nv)then
        do i=1,nele
          k1=klp(i)
          do j=isp0+1,isp
            if(itastk(1,j) .eq. i)then
              call elcompl(i,kl)
              do l=1,kl%nl
                k=int(kl%rbody(l))
                if(k .ne. k1 .and. icomp(k) .eq. k1
     $               .and. itastk(2,j) .ne. ival(i))then
                  call tfvcopy(k1,k,itastk(2,j),1.d0)
                endif
              enddo
            endif
          enddo
        enddo
      endif
      if(var)then
        isp1=isp
        isp=isp+1
        do l=1,nlat-1
          ie1=iele1(icomp(l))
          itastk(1,isp)=ie1
          itastk(2,isp)=ival(ie1)
          do j=isp0+1,isp1
            if(ktastk(j) .eq. ktastk(isp))then
              call tfvcopy(icomp(l),l,ival(ie1),couple(l))
            endif
          enddo
        enddo
        isp=isp1
      endif
      call tffscoupledvar(irtc)
      if(irtc .ne. 0)then
        call termes(icslfno(),'?Error in CoupledVariables',' ')
      endif
      return
      end

      subroutine tffscoupledvar(irtc)
      use tfstk
      use tfcsi, only:icslfno
      implicit none
      integer*4 ia,irtc,itfdownlevel,irtc1
      type (sad_dlist), pointer :: kl
      type (sad_descriptor) ifcoupv,ifsetcoup,k,kx
      type (sad_symdef), pointer, save :: symdcoupv
      data ifcoupv%k,ifsetcoup%k /0,0/
      if(ifcoupv%k .eq. 0)then
        ifcoupv  =kxsymbolz('`EVList',7)
        call descr_sad(ifcoupv,symdcoupv)
        ifsetcoup=kxsymbolz('SetCoupledElements',18)
      endif
      irtc=0
      k=symdcoupv%value
      if(tflistq(k,kl))then
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

      subroutine tfsavevar(ie,ntou)
      use tfstk
      use ffs, only:nvevx
      use ffs_pointer
      use sad_main
      use ffs_seg
      implicit none
      type (sad_comp), pointer :: cmps,cmpd
      integer*4 ntou,i,ie
      do i=1,ntou
        if(nvevx(i)%itouchele .eq. ie)then
          call loc_comp(latt(klp(ie)),cmps)
          call loc_comp(idvalc(klp(ie)),cmpd)
          call tfvcopycmp(cmps,cmpd,nvevx(i)%itouchv,1.d0)
        endif
      enddo
      return
      end

