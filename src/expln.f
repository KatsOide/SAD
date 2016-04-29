      subroutine expln(idxl)
      use maccbk
      use tfstk
      implicit none
      include 'inc/MACCODE.inc'
      include 'inc/MACKW.inc'
      include 'inc/MACMISC.inc'
      include 'inc/MACexpn.inc'
      integer*4 idxl,idx0,isp0,ia,i,idi,iti,plen,orientation,
     $     itcaloc,ip,idxerr,idxpar,n
      real*8 frand
      if(ilist(2,idval(idxl)) .gt. 0) return
     
      entry expnln(idxl)

      idx0=idval(idxl)
      if(ilist(2,idx0) .gt. 0)then
        call tfree(int8(ilist(2,idx0)))
      endif
      ilist(2,idx0)=0
      isp0=isp
c     first stage
      call explnstk(idxl,1)
c      write(*,*)'expln: first stage done. ',isp,isp0
c     
c     Second stage : expand parameter list of each element
c     
      n=isp-isp0
      ia=itcaloc(n+1)
      ilist(1,ia)=n
      ilist(2,ia)=0
      ilist(2,idx0)=ia
      do i=isp0+1,isp
        idi=idval(itastk(1,i))
        plen=ilist(1,idi)
        orientation=sign(1,itastk(2,i))
        ip=itcaloc(plen+1+expnsize)
c        write(*,*)'expln-i ',i-isp0,ip,
c     $       itastk(1,i),idtype(itastk(1,i)),plen
        ilist(1,ia+i-isp0)=itastk(1,i)
        ilist(2,ia+i-isp0)=ip
        ilist(1,ip)=plen+expnsize
        ilist(2,ip)=0
        rlist(ip+plen+expnsize)=dble(orientation)
c     set Nominal values
        rlist(ip+1:ip+plen)=rlist(idi+1:idi+plen)
c     and then add statistical error
        idxerr=ilist(2,idi)
        do while(idxerr .ne. 0)
          iti=idtype(itastk(1,i))
          if (kytbl(ilist(1,idxerr),iti) .ne. 0) then
            idxpar=ip+kytbl(ilist(1,idxerr),iti)
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
      use tfstk
      implicit none
      integer*4 idxl
      include 'inc/MACCODE.inc'
      include 'inc/MACKW.inc'
      include 'inc/MACMISC.inc'
      include 'inc/MACexpn.inc'
      integer*4 idx,direct,i1,i2,i,idxi2,dir,j,lpname
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
c        write(*,*)'explnstk-i ',i-i1+1,direct,ilist(1,i),dir
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
     $         //' is not a element ',0,16)
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
