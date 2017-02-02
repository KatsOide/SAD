      subroutine tfsave(word,cmd0,ntouch)
      use tfstk
      use ffs
      use ffs_pointer
      use tffitcode
      implicit none
      integer*8 j
      integer*4 i,l,lt,k,next
      integer*4 ntouch
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
      do i=1,nele
        l=klp(i)
        if(.not. cmd .or. tmatch(pname(ilist(2,latt(l))),word))then
          if(cmd .and. .not. exist1)then
            call cssetp(next)
            exist1=.true.
            exist=.true.
          endif
          j=idval(ilist(2,latt(l)))
          lt=idtype(ilist(2,latt(l)))
          if(all)then
            do k=1,kytbl(kwMAX,lt)-1
              rlist(j+k)=rlist(latt(l)+k)
            enddo
          endif
          call tfsavevar(i,ntouch)
          if(lt .eq. icMARK)then
            if(l .eq. 1)then
              do 30 k=1,ntwissfun
                rlist(j+k)=rlist(latt(1)+k)
 30           continue
            else
              do 31 k=1,ntwissfun
                rlist(j+k)=twiss(l,0,k)
 31           continue
c              write(*,*)'tfsave ',l,
c     $             ilist(1,latt(l)+ilist(1,latt(l))),
c     $             rlist(j+1),rlist(j+4)
              if(rlist(latt(l)+ilist(1,latt(l))) .lt. 0.d0)then
                rlist(j+mfitax)=-rlist(j+mfitax)
                rlist(j+mfitay)=-rlist(j+mfitay)
                rlist(j+mfitepx)=-rlist(j+mfitepx)
                rlist(j+mfitepy)=-rlist(j+mfitepy)
                rlist(j+mfitr2)=-rlist(j+mfitr2)
                rlist(j+mfitr3)=-rlist(j+mfitr3)
                rlist(j+mfitdpx)=-rlist(j+mfitdpx)
                rlist(j+mfitdpy)=-rlist(j+mfitdpy)
              endif
            endif
            rlist(j+kytbl(kwEMIX,icMARK))=emx
            rlist(j+kytbl(kwEMIY,icMARK))=emy
            rlist(j+kytbl(kwDP,icMARK))=dpmax
          elseif(ival(i) .gt. 0)then
            rlist(j+ival(i))=rlist(latt(l)+ival(i))
     $           /errk(1,l)
          endif
        endif
      enddo
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

      subroutine tfrst(word,cmd0,nvar,ntouch)
      use tfstk
      use ffs
      use ffs_pointer
      use tffitcode
      implicit none
      integer*8 j
      integer*4 l,i,lt,k,next
      integer*4 nvar,ntouch
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
        endif
      endif
      do i=1,nele
        l=klp(i)
        if(.not. cmd .or. tmatch(pname(ilist(2,latt(l))),word))then
          if(cmd .and. .not. exist1)then
            call cssetp(next)
            exist1=.true.
            exist=.true.
          endif
          lt=idtype(ilist(2,latt(l)))
          if(lt .eq. icMARK)then
            if(l .eq. 1)then
              call tmov(rlist(idval(ilist(2,latt(1)))+1),
     $             rlist(latt(1)+1),ntwissfun)
            endif
          else
            j=idval(ilist(2,latt(l)))
            if(all)then
              do k=1,kytbl(kwMAX,lt)-1
                rlist(latt(l)+k)=rlist(j+k)
              enddo
            endif
            do k=1,nvar
              if(ivarele(k) .eq. i)then
                valvar(k)=rlist(j+ivvar(k))
              endif
            enddo
            do k=1,ntouch
              if(itouchele(k) .eq. i)then
                rlist(latt(l)+itouchv(k))=rlist(j+itouchv(k))
              endif
            enddo
            if(ival(i) .gt. 0)then
              rlist(latt(l)+ival(i))=rlist(j+ival(i))*errk(1,l)
            endif
          endif
        endif
      enddo
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
 9000 call tffsadjust(ntouch)
      call tfinitvar(nvar)
      return
      end

      subroutine tffsresetall
      use tfstk
      use ffs
      use ffs_pointer
      use tffitcode
      implicit none
      integer*8 j,m
      integer*4 l,k,lt
      do l=1,nlat-1
        j=idval(ilist(2,latt(l)))
        lt=idtype(ilist(2,latt(l)))
        m=latt(l)
        do k=1,kytbl(kwMAX,lt)-1
          rlist(m+k)=rlist(j+k)
        enddo
      enddo
      return
      end
