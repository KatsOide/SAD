      subroutine tfsave(word,cmd0,ntou)
      use kyparam
      use tfstk
      use ffs
      use ffs_pointer
      use tffitcode
      use tfcsi,only:ipoint
      use ffs_seg
      implicit none
      type (sad_comp), pointer ::cmp,cmpd
      integer*8 j
      integer*4 i,l,lt,next,itfdownlevel
      integer*4 ,intent(in):: ntou
      character*(*) ,intent(in):: word
      logical*4 ,intent(in):: cmd0
      logical*4 cmd,exist,tmatch,exist1,all
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
          ipoint=next
          all=.true.
          go to 1
        endif
      endif
      levele=levele+1
      do i=1,nele
        l=nelvx(i)%klp
        if(.not. cmd .or. tmatch(pnamec(l),word))then
          if(cmd .and. .not. exist1)then
            ipoint=next
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
              rlist(j+1:j+ntwissfun)=cmp%value(1:ntwissfun)
            else
              rlist(j+1:j+ntwissfun)=twiss(l,0,1:ntwissfun)
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
            rlist(j+ky_EMIZ_MARK)=emz
            rlist(j+ky_SIGZ_MARK)=sigzs
            rlist(j+ky_SIGE_MARK)=sizedp
            rlist(j+ky_DP_MARK)=dpmax
          elseif(nelvx(i)%ival .gt. 0)then
            call tfvcopycmp(cmp,cmpd,nelvx(i)%ival,1.d0/errk(1,l))
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
      use tfcsi,only:ipoint
      use ffs_seg
      implicit none
      type (sad_comp), pointer :: cmp,cmps
      integer*8 j
      integer*4 l,i,lt,k,next,itfdownlevel
      character*(*) ,intent(in):: word
      logical*4 ,intent(in):: cmd0
      logical*4 cmd,exist,tmatch,exist1,all
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
          ipoint=next
          all=.true.
          go to 1
        elseif(word .eq. 'TOTAL')then
          ipoint=next
          call tffsresetall
          go to 9000
        endif
      endif
      levele=levele+1
      do i=1,nele
        l=nelvx(i)%klp
        if(.not. cmd .or. tmatch(pnamec(l),word))then
          if(cmd .and. .not. exist1)then
            ipoint=next
            exist1=.true.
            exist=.true.
          endif
          lt=idtypec(l)
          call compelc(l,cmp)
          cmp%update=cmp%nparam .le. 0
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
              if(nvevx(k)%ivarele .eq. i)then
                nvevx(k)%valvar=tfvalvar(l,nvevx(k)%ivvar)
c                valvar(k)=rlist(j+ivvar(k))
              endif
            enddo
            do k=1,flv%ntouch
              if(nvevx(k)%itouchele .eq. i)then
c                write(*,*)'tfrst-itouch ',k,i,itouchv(k)
                call tfvcopycmp(cmps,cmp,nvevx(k)%itouchv,1.d0)
c                cmp%value(itouchv(k))=rlist(j+itouchv(k))
              endif
            enddo
            if(nelvx(i)%ival .gt. 0)then
              call tfvcopycmp(cmps,cmp,nelvx(i)%ival,errk(1,l))
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
        cmp%update=cmp%nparam .le. 0
        if(nelvx(l)%ival .gt. 0)then
          cmp%value(nelvx(l)%ival)=cmp%value(nelvx(l)%ival)*errk(1,l)
        endif
      enddo
      return
      end
