      subroutine expln(idxl)
      use maccbk
      implicit none
      include 'inc/MACexpn.inc'
      integer*4 idxl
      if(klist(idval(idxl)) .gt. 0) return
      call expnln(idxl)
      return
      end

      subroutine expnln(idxl)
      use tfstk, only:ilist,klist,itastk,isp,idval,idtype,rlist
      implicit none
      include 'inc/MACCODE.inc'
      include 'inc/MACKW.inc'
      include 'inc/MACexpn.inc'
      integer*4 isp0,i,idi,iti,plen,orientation,
     $     itcaloc,idxerr,n,j,idxl
      integer*8 ktcaloc,kp,idxpar,idx0,ia
      real*8 frand
c      write(*,*)'expnln-0 ',idxl
      idx0=idval(idxl)
c      write(*,*)'expnln-0.1 ',idx0
      if(klist(idx0) .gt. 0)then
        call tfreeln(klist(idx0))
      endif
c      write(*,*)'expnln-0.2 ',ilist(2,idx0)
      klist(idx0)=0
      isp0=isp
c     first stage
      call explnstk(idxl,1)
c      write(*,*)'expln: first stage done. ',isp,isp0
c     
c     Second stage : expand parameter list of each element
c     
      n=isp-isp0
      ia=ktcaloc(n+1)
      ilist(1,ia)=n
      ilist(idx0)=ia
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
      isp=isp0
      return
      end

      recursive subroutine explnstk(idxl,direct)
      use maccbk
      use tfstk, only:ilist,itastk,isp,idval,idtype
      implicit none
      include 'inc/MACCODE.inc'
      include 'inc/MACKW.inc'
      integer*8 idx,i,i1,i2,idxi2
      integer*4 direct,dir,j,lpname,idxl
      idx=idval(idxl)
      if(direct .gt. 0)then
        i1=idx+1
        i2=idx+ilist(1,idx-1)-2
      else
        i2=idx+1
        i1=idx+ilist(1,idx-1)-2
      endif
      do i=i1,i2,direct
        idxi2=klist(i)
        dir=sign(1,direct*ilist(2,i-1))
        if (idtype(idxi2) .eq. icLINE) then
          if(idxl .eq. idxi2) then
            call errmsg('expnln',
     &           'definition contains a loop.',0,16)
          endif
          do j=1,abs(ilist(2,i-1))
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
          do j=1,abs(ilist(2,i-1))
            isp=isp+1
            ktastk(isp)=idxi2
            itastk2(1,isp)=dir
          enddo
        endif
      enddo
      return
      end

      subroutine tfreeln(ip)
      use maccbk
      use tfstk, only:ilist,klist
      implicit none
      integer*4 i
      integer*8 k,ip
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
