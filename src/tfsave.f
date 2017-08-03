      subroutine tfsave(word,cmd0,ntou)
      use kyparam
      use tfstk
      use ffs
      use ffs_pointer
      use tffitcode
      use tfcsi,only:cssetp
      use ffs_seg
      implicit none
      type (sad_comp), pointer ::cmp,cmpd
      integer*8 j
      integer*4 i,l,lt,next,itfdownlevel
      integer*4 ntou
      character*(*) word
      logical*4 cmd,cmd0,exist,tmatch,exist1,all
      call tfsetparam
      cmd=cmd0
      exist=.false.
      all=.false.
 1    exist1=.false.
      if(cmd)then
        call peekwd(word,next)
        if(word .eq. ' ')then
          if(exist)then
            return
          endif
        elseif(word .eq. 'ALL')then
          call cssetp(next)
          all=.true.
          go to 1
        endif
      endif
      levele=levele+1
      do i=1,nele
        l=klp(i)
        if(.not. cmd .or. tmatch(pnamec(l),word))then
          if(cmd .and. .not. exist1)then
            call cssetp(next)
            exist1=.true.
            exist=.true.
          endif
          j=idvalc(l)
          lt=idtypec(l)
          call compelc(l,cmp)
          call loc_comp(j,cmpd)
          if(all)then
            call tfvcopycmpall(cmp,cmpd,kytbl(kwMAX,lt)-1)
          endif
          call tfsavevar(i,ntou)
          if(lt .eq. icMARK)then
            if(l .eq. 1)then
c              do 30 k=1,ntwissfun
              rlist(j+1:j+ntwissfun)=cmp%value(1:ntwissfun)
c     $             rlist(latt(1)+1:latt(1)+ntwissfun)
c                rlist(j+k)=rlist(latt(1)+k)
c 30           continue
            else
c              do 31 k=1,ntwissfun
              rlist(j+1:j+ntwissfun)=twiss(l,0,1:ntwissfun)
c                rlist(j+k)=twiss(l,0,k)
c 31           continue
c              write(*,*)'tfsave ',l,
c     $             ilist(1,latt(l)+ilist(1,latt(l))),
c     $             rlist(j+1),rlist(j+4)
              if(direlc(l) .lt. 0.d0)then
                rlist(j+mfitax)=-rlist(j+mfitax)
                rlist(j+mfitay)=-rlist(j+mfitay)
                rlist(j+mfitaz)=-rlist(j+mfitaz)
                rlist(j+mfitepx)=-rlist(j+mfitepx)
                rlist(j+mfitepy)=-rlist(j+mfitepy)
                rlist(j+mfitzpx)=-rlist(j+mfitzpx)
                rlist(j+mfitzpy)=-rlist(j+mfitzpy)
                rlist(j+mfitr2)=-rlist(j+mfitr2)
                rlist(j+mfitr3)=-rlist(j+mfitr3)
                rlist(j+mfitdpx)=-rlist(j+mfitdpx)
                rlist(j+mfitdpy)=-rlist(j+mfitdpy)
              endif
            endif
            rlist(j+ky_EMIX_MARK)=emx
            rlist(j+ky_EMIY_MARK)=emy
            rlist(j+ky_DP_MARK)=dpmax
          elseif(ival(i) .gt. 0)then
            call tfvcopycmp(cmp,cmpd,ival(i),1.d0/errk(1,l))
c            rlist(j+ival(i))=cmp%value(ival(i))/errk(1,l)
          endif
        endif
      enddo
      l=itfdownlevel()
      if(cmd)then
        if(exist)then
          if(exist1)then
            go to 1
          else
            return
          endif
        else
          if(exist1)then
            exist=exist1
          else
            cmd=.false.
          endif
          go to 1
        endif
      endif
      return
      end

      subroutine tfrst(word,cmd0)
      use tfstk
      use ffs
      use ffs_pointer
      use tffitcode
      use tfcsi,only:cssetp
      use ffs_seg
      implicit none
      type (sad_comp), pointer :: cmp,cmps
      integer*8 j
      integer*4 l,i,lt,k,next,itfdownlevel
      character*(*) word
      logical*4 cmd,cmd0,exist,tmatch,exist1,all
      call tfsetparam
      cmd=cmd0
      all=.false.
      exist=.false.
 1    exist1=.false.
      if(cmd)then
        call peekwd(word,next)
        if(word .eq. ' ')then
          if(exist)then
            go to 9000
          endif
        elseif(word .eq. 'ALL')then
          call cssetp(next)
          all=.true.
          go to 1
        elseif(word .eq. 'TOTAL')then
          call cssetp(next)
          call tffsresetall
          go to 9000
        endif
      endif
      levele=levele+1
      do i=1,nele
        l=klp(i)
        if(.not. cmd .or. tmatch(pnamec(l),word))then
          if(cmd .and. .not. exist1)then
            call cssetp(next)
            exist1=.true.
            exist=.true.
          endif
          lt=idtypec(l)
          call compelc(l,cmp)
          cmp%update=0
          j=idvalc(l)
          call loc_comp(j,cmps)
          if(lt .eq. icMARK)then
            if(l .eq. 1)then
              cmp%value(1:ntwissfun)=cmps%value(1:ntwissfun)
            endif
          else
            if(all)then
              call tfvcopycmpall(cmps,cmp,kytbl(kwMAX,lt)-1)
            endif
            do k=1,flv%nvar
              if(ivarele(k) .eq. i)then
                valvar(k)=tfvalvar(l,ivvar(k))
c                valvar(k)=rlist(j+ivvar(k))
              endif
            enddo
            do k=1,flv%ntouch
              if(itouchele(k) .eq. i)then
                call tfvcopycmp(cmps,cmp,itouchv(k),1.d0)
c                cmp%value(itouchv(k))=rlist(j+itouchv(k))
              endif
            enddo
            if(ival(i) .gt. 0)then
              call tfvcopycmp(cmps,cmp,ival(i),errk(1,l))
c              cmp%value(ival(i))=rlist(j+ival(i))*errk(1,l)
            endif
          endif
        endif
      enddo
      l=itfdownlevel()
      if(cmd)then
        if(exist)then
          if(exist1)then
            go to 1
          else
            go to 9000
          endif
        else
          if(exist1)then
            exist=exist1
          else
            cmd=.false.
          endif
          go to 1
        endif
      endif
 9000 call tffsadjust
      call tfinitvar
      return
      end

      subroutine tffsresetall
      use tfstk
      use ffs
      use ffs_pointer
      use tffitcode
      use ffs_seg
      implicit none
      type (sad_comp), pointer :: cmp,cmps
      integer*8 j
      integer*4 l,lt
      do l=1,nlat-1
        j=idvalc(l)
        lt=idtypec(l)
        call compelc(l,cmp)
        call loc_comp(j,cmps)
        call tfvcopycmpall(cmps,cmp,kytbl(kwMAX,lt)-1)
        cmp%update=0
        if(ival(l) .gt. 0)then
          cmp%value(ival(l))=cmp%value(ival(l))*errk(1,l)
        endif
c        cmp%value(1:kytbl(kwMAX,lt)-1)=
c     $       rlist(j+1:j+kytbl(kwMAX,lt)-1)
      enddo
      return
      end
