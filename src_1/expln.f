      subroutine expln(idxl)
      use maccbk
      implicit none
      include 'inc/MACexpn.inc'
      integer*4 idxl
      if(ilist(2,idval(idxl)) .le. 0)then
        call expnln(idxl)
      endif
      return
      end

      subroutine expnln(idxl)
      use tfstk, only:ilist,klist,itastk,isp,idval,idtype,rlist
      use maccbk
      implicit none
      include 'inc/MACCODE.inc'
      include 'inc/MACKW.inc'
      include 'inc/MACexpn.inc'
      integer*4 idxl,iti,plen,orientation,
     $     idxerr,n,j,hsrchz,idxe,lenw
      integer*8 ktcaloc,kp,idxpar,ia,idx0,idi,isp0,i
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
      n=int(isp-isp0)
      idxe=hsrchz(pname(idxl)(1:lenw(pname(idxl)))//"$EXPND")
      ia=ktcaloc(n+1)
      idval(idxe)=ia
c      write(*,*)'expnln-0.3 ',idxe,ia,n
      ilist(1,ia)=n
      ilist(2,ia)=0
      ilist(2,idx0)=idxe
      do j=1,n
        i=isp0+j
c         write(*,*)'explnstk-i ',i,itastk(1,i)
        idi=idval(itastk(1,i))
        plen=ilist(1,idi)
        orientation=sign(1,itastk(2,i))
        kp=ktcaloc(plen+2+expnsize)
c        write(*,*)'expln-i ',j,kp,
c     $       itastk(1,i),idtype(itastk(1,i)),plen
        klist(ia+j)=kp
        ilist(1,kp)=plen+expnsize
        ilist(2,kp)=itastk(1,i)
        rlist(kp+plen+expnsize)=dble(orientation)
c     set Nominal values
        rlist(kp+1:kp+plen)=rlist(idi+1:idi+plen)
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
c      last=ktcaloc(4)
c      ilist(2,last-1)=0
c      klist(ia+n+1)=last
c         write(*,*)'explnstk-9 ',isp,isp0
      isp=isp0
      return
      end

      recursive subroutine explnstk(idxl,direct)
      use maccbk
      use tfstk, only:ilist,itastk,isp,idval,idtype
      implicit none
      integer*8 idx,i1,i2,i
      integer*4 idxl
      include 'inc/MACCODE.inc'
      include 'inc/MACKW.inc'
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
      use tfstk, only:ilist,klist
      implicit none
      integer*4 i,idxe
      integer*8 k,ip
      ip=idval(idxe)
      idval(idxe)=0
      do i=1,ilist(1,ip)
        k=klist(ip+i)
        if(klist(k+ilist(1,k)+1) .gt. 0)then
          call tfree(klist(k+ilist(1,k)+1))
        endif
        call tfree(k)
      enddo
      call tfree(ip)
      return
      end
