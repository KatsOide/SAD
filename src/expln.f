      subroutine expln(idxl)
      use maccbk
      implicit none
      integer*4 idxl
      if(ilist(2,idval(idxl)) .le. 0)then
        call expnln(idxl)
      endif
      return
      end

      subroutine expnln(idxl)
      use tfstk, only:ilist,itastk,isp,idval,idtype,rlist
      use maccbk
      use sad_main
      use mackw
      implicit none
      type (sad_el), pointer :: el
      type (sad_comp), pointer :: cmp,cmps
      integer*4 idxl,isp0,i,iti,plen,orientation,
     $     idxerr,n,j,hsrchz,idxe,lpname,k
      integer*8 kp,idxpar,ia,idx0,idi
      real*8 frand
c      write(*,*)'expnln-0 ',idxl
      idx0=idval(idxl)
c      write(*,*)'expnln-0.1 ',pname(idxl)
      if(ilist(2,idx0) .gt. 0)then
        call tfreeln(ilist(2,idx0))
      endif
c      write(*,*)'expnln-0.2 ',ilist(2,idx0)
      ilist(2,idx0)=0
      isp0=isp
c     first stage
      call explnstk(idxl,1)
c      write(*,*)'expln: first stage done. ',isp,isp0
c     
c     Second stage : expand parameter list of each element
c
      n=isp-isp0
      idxe=hsrchz(pname(idxl)(:lpname(idxl))//"$EXPND")
      ia=kmelaloc(n,el)
      idval(idxe)=ia
c      write(*,*)'expnln-0.3 ',idxe,ia,n
      el%aux=0
      ilist(2,idx0)=idxe
      do j=1,n
        i=isp0+j
        idi=idval(itastk(1,i))
        iti=idtype(itastk(1,i))
        plen=ilist(1,idi)
        orientation=sign(1,itastk(2,i))
        kp=kmcompaloc(plen+kytbl(kwNPARAM,iti),cmp)
        el%comp(j)=kp
        cmp%id=itastk(1,i)
        cmp%nparam=kytbl(kwNPARAM,iti)
        cmp%orient=dble(orientation)
        call loc_comp(idi,cmps)
c     set Nominal values
        call tfvcopycmpall(cmps,cmp,plen)
c        do k=1,plen
c          call tfvcopycmp(cmps,cmp,k,1.d0)
c        enddo
c        cmp%value(1:plen)=rlist(idi+1:idi+plen)
c     and then add statistical error
        idxerr=ilist(2,idi)
        do while(idxerr .ne. 0)
          iti=idtype(itastk(1,i))
          if (kytbl(ilist(1,idxerr),iti) .ne. 0) then
            idxpar=kp+kytbl(ilist(1,idxerr),iti)
            if (ilist(2,idxerr+1) .gt. 0) then
              rlist(idxpar)=rlist(idxpar)+ frand(idxerr+1)
            else
              rlist(idxpar)=rlist(idxpar)
     &             +rlist(idxerr+2-ilist(2,idxerr+1))
              ilist(2,idxerr+1)=-mod(-ilist(2,idxerr+1)+1
     &             ,ilist(1,idxerr+1))
            endif
          endif
          idxerr=ilist(2,idxerr)
        enddo
      enddo
      isp=isp0
      return
      end

      recursive subroutine explnstk(idxl,direct)
      use maccbk
      use tfstk, only:ilist,itastk,isp,idval,idtype
      use maccode
      use mackw
      implicit none
      integer*8 idx,i1,i2,i
      integer*4 idxl
      integer*4 direct,idxi2,dir,j,lpname
      idx=idval(idxl)
      if(direct .gt. 0)then
        i1=idx+1
        i2=idx+ilist(1,idx)
      else
        i2=idx+1
        i1=idx+ilist(1,idx)
      endif
      do i=i1,i2,direct
        idxi2=ilist(2,i)
        dir=sign(1,direct*ilist(1,i))
        if (idtype(idxi2) .eq. icLINE) then
          if(idxl .eq. idxi2) then
            call errmsg('expnln',
     &           'definition contains a loop.',0,16)
          endif
          do j=1,abs(ilist(1,i))
            call explnstk(idxi2,dir)
          enddo
        else if (idtype(idxi2) .gt. icMXEL) then
          call errmsg('expnln',
     &         pname(idxi2)(:lpname(idxi2))
     $         //' is not an element ',0,16)
        else if (idtype(idxi2) .eq. icNULL) then
          call errmsg('expnln',
     &         pname(idxi2)(:lpname(idxi2))
     $         //' is not defined yet',0,0)
          call errmsg('expnln',
     &         'warnig:unable to expand line',0,0)
        else
          do j=1,abs(ilist(1,i))
            isp=isp+1
            itastk(1,isp)=idxi2
            itastk(2,isp)=dir
          enddo
        endif
      enddo
      return
      end

      subroutine tfreeln(idxe)
      use maccbk
      use mackw
      use tfstk, only:ilist,klist,tfree,tflocald
      use sad_main
      implicit none
      type (sad_el), pointer :: el
      type (sad_comp), pointer :: cmp
      integer*4 i,idxe,j
      integer*8 k,ip
      ip=idval(idxe)
      call loc_el(ip,el)
      idval(idxe)=0
      do i=1,el%nlat1-2
        k=el%comp(i)
        call loc_comp(k,cmp)
        do j=1,kytbl(kwMAX,idtype(cmp%id))-1
          call tflocald(cmp%dvalue(j))
        enddo
        call tfree(k-1)
      enddo
      call tfree(ip)
      return
      end
