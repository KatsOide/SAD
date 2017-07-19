      subroutine tfsetv(nvar)
      use tfstk
      use ffs, only:nlat
      use ffs_pointer
      use mackw
      implicit none
      type (sad_comp), pointer :: cmp
      integer*4 nvar
      integer*4 i,ie,iv,j,ii,ie1
      if(nvar .gt. 0)then
        do i=1,nlat-1
          call compelc(i,cmp)
          ii=iele(i)
          ie=iele1(ii)
          ie1=iele1(i)
          do j=1,nvar
            iv=ivvar(j)
            if(iv .eq. ival(ie) .and. ivarele(j) .eq. ie
     $           .and. (ivcomp(j) .eq. 0 .or. ivcomp(j) .eq. ii))then
              call tfsetcmp(valvar(j)*errk(1,i)*couple(i),cmp,iv)
c              cmp%value(iv)=valvar(j)*errk(1,i)*couple(i)
              cmp%update=0
            elseif(iv .ne. 0 .and. iv .ne. ival(ie) .and.
     $             ivarele(j) .eq. ie1
     $             .and. (ivcomp(j) .eq. 0 .or. ivcomp(j) .eq. i))then
              call tfsetcmp(valvar(j),cmp,iv)
c              cmp%value(iv)=valvar(j)
              cmp%update=0
            endif
            if(ivarele(j) .gt. ie)then
              exit
            endif
          enddo
        enddo
        call tfinitvar(nvar)
      endif
      return
      end

      subroutine tfinitvar(nvar)
      use tfstk
      use ffs_pointer
      use tfcsi, only:icslfno
      implicit none
      integer*4 nvar,k,irtc
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
          valvar(i)=tfvalvar(k,iv)/errk(1,k)
        else
          valvar(i)=tfvalvar(k,iv)
        endif
      enddo
      return
      end

      subroutine tfsavevar(ie,ntouch)
      use tfstk
      use ffs_pointer
      use sad_main
      implicit none
      type (sad_comp), pointer :: cmps,cmpd
      integer*4 ntouch,i,ie
      do i=1,ntouch
        if(itouchele(i) .eq. ie)then
          call loc_comp(latt(klp(ie)),cmps)
          call loc_comp(idvalc(klp(ie)),cmpd)
          call tfvcopycmp(cmps,cmpd,itouchv(i),1.d0)
        endif
      enddo
      return
      end

      subroutine tffsadjust(ntouch)
      use tfstk
      use ffs
      use ffs_pointer
      use tffitcode
      use tfcsi, only:icslfno
      implicit none
      integer*4 itv(ntouch+1),ite(nele+1)
      integer*4 ntouch,i,ie,j,irtc,ie1,ntv,k
      ite=0
      ntv=0
c      write(*,*)'tffsadjust ',ntouch
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
          if(iele(i) .ne. i .and. ie .ne. 0 .and. ival(ie) .ne. 0)then
            call tfvcopy(iele(i),i,ival(ie),
     $           errk(1,i)*couple(i)/errk(1,iele(i)))
          endif
        enddo
      else
        do i=1,nlat-1
          ie=iele1(iele(i))
          ie1=iele1(i)
          if(ie .ne. 0 .and. i .ne. iele(i) .and. ival(ie) .ne. 0)then
            call tfvcopy(iele(i),i,ival(ie),
     $           errk(1,i)*couple(i)/errk(1,iele(i)))
          endif
          if(ite(ie1) .ne. 0)then
            do k=1,ntv
              j=itv(k)
              if(itouchele(j) .eq. ie1 .and. klp(ie1) .ne. i)then
                call tfvcopy(klp(ie1),i,itouchv(j),1.d0)
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
      use tfcsi, only:icslfno
      implicit none
      integer*4 ia,irtc,itfdownlevel,irtc1
      type (sad_list), pointer :: kl
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
