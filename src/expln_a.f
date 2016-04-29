      recursive subroutine expln(idxl)
      use maccbk
      implicit none
      integer*4 idxl
      include 'inc/MACCODE.inc'
      include 'inc/MACKW.inc'
      include 'inc/MACMISC.inc'
      include 'inc/MACexpn.inc'
      integer*4 STACSZ
      parameter (STACSZ=100)
      integer*4 stack(3,STACSZ),pstack
      integer*4 llen,idx,plen,direct,orientation,idx0,idxi2
      integer*4 blksz,adr,newadr,i,i0,k,j,idxerr,idxpar,ip
      integer*4 lpname,mctaloc,mtaloc,idi
      real*8 frand
      direct=1

      adr(idx0)=ilist(2,idx0)
c     
      if(ilist(2,idval(idxl)) .gt. 0) return
c     
      entry expnln(idxl)
c     
      llen=0
c     first stage : expand element names into line
      idx=idval(idxl)
      idx0=idx
      blksz=max(pagesz,ilist(1,ilist(2,idx))+1)
      if(ilist(2,idx) .gt. 0) then
         call errmsg('expln',
     &        pname(idxl)(1:lpname(idxl))//
     $        ' is re-expanded again',0,0)
         do i=1,ilist(1,ilist(2,idx))
           call tfree(ilist(2,adr(idx0)+i))
c           call tfreem(idi,ilist(1,idi))
c    The original was below.  7/25/2011 K. Oide
c            call freeme(ilist(2,adr(idx0)+i),
c     &           ilist(1,ilist(2,adr(idx0)+i))+1)
         end do
         call tfree(ilist(2,idx))
c         call freeme(ilist(2,idx),ilist(1,ilist(2,idx))+1)   Is this +1 right???? 7/25/2011 K. Oide
      endif
      newadr=italoc(blksz)
c      newadr=mfalloc(blksz)
      if(newadr .le. 0) then
         call errmsg('expln',' no more free blocks',32,0)
         stop
      end if
      direct=1
      i0=0
      isp0=isp
      do while(i0 .le. ilist(1,idx))
        i0=i0+1
        i=((1-direct)*(ilist(1,idx)+1)+ 2*direct*i0)/2
        idxi2=ilist(2,idx+i)
        if (idtype(idxi2) .eq. icLINE) then
          if(idxl .eq. idxi2) then
            call errmsg('expnln',
     &           'definition contains a loop.',0,16)
          endif
          call expln(idxi2)


          isp=isp+1
          itastk(1,isp)=idx
          itastk(2,isp)=i0
          ivstk(1,isp)=direct
          if(ilist(1,idx+i).lt.0) then
            direct=-1*direct
          endif
          do k=1,abs(ilist(1,idx+i))-1
            isp=isp+1
            itastk(1,isp)=idval(idxi2)
            itastk(2,isp)=0
            ivstk(3,isp)=direct
          enddo
          idx=idval(idxi2)
          i0=0
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
          orientation=direct*isign(1,ilist(1,idx+i))
          do j=1,abs(ilist(1,idx+i))
            newadr=italoc(abs(ilist(1,idx+i))+1)



              ilist(2,newadr)=0
              do k=1,llen-1
                ilist(1,newadr+k)=ilist(1,adr(idx0)+k)
                ilist(2,newadr+k)=ilist(2,adr(idx0)+k)
              end do
              call tfreem(adr(idx0),blksz)
c     call freeme(adr(idx0),blksz)
              blksz=2*blksz
              ilist(2,idx0)=newadr
            end if
            ilist(1,adr(idx0)+llen)=idxi2
            ilist(2,adr(idx0)+llen)=orientation
          enddo
        endif
      enddo
      ilist(2,newadr)=0
      ilist(2,idx0)=newadr
c     
 1400 continue
      if ( pstack .gt. 0) then
         idx=stack(1,pstack)
         i0=stack(2,pstack)
         direct=stack(3,pstack)
         pstack=pstack-1
         go to 1100
      endif
      ilist(1,ilist(2,idx0))=llen
      ilist(2,ilist(2,idx0))=0
      if( llen  .gt. blksz-4 ) call errmsg('doexpn','expanded line is'
     &     //' too long.',0,32)
      call tfreem(ilist(2,idx0)+llen+1,blksz-llen-1)
c      call freeme(ilist(2,idx0)+llen+1,blksz-llen-1)
c     
c     Second stage : expand parameter list of each element
c     
      do i=ilist(2,idx0)+1,ilist(2,idx0)+llen
         plen=ilist(1, idval(ilist(1, i)))
         orientation=ilist(2,i)
         ip=mctaloc(plen + 1 + expnsize)
c    The original was below.  7/25/2011 K. Oide
c         ip=mcfallo(plen + 1 + expnsize )
c         write(*,*)'expln ',i,ip,orientation,plen,expnsize,
c     $        ' ',pname(ilist(1,i))
         ilist(2,i)=ip
         ilist(1,ip)=plen+expnsize
         ilist(2,ip)=0
         rlist(ip+plen+expnsize)=orientation
c     set Nominal values
         do  j=1,plen-1
           rlist(ip+j)=rlist(idval(ilist(1,i))+j)
         enddo
c     and then add statistical error
         idxerr=ilist(2,idval(ilist(1,i)))
 2200    continue
         if (idxerr .eq. 0) cycle
         if (kytbl(ilist(1,idxerr),idtype(ilist(1,i))) .ne. 0) then
            idxpar=ilist(2,i)
     &           +kytbl(ilist(1,idxerr),idtype(ilist(1,i)))
            if (ilist(2,idxerr+1) .gt. 0) then
               rlist(idxpar)=rlist(idxpar)
     &              + frand(idxerr+1)
            else
               rlist(idxpar)=rlist(idxpar)
     &              +rlist(idxerr+2-ilist(2,idxerr+1))
               ilist(2,idxerr+1)=-mod(-ilist(2,idxerr+1)+1
     &              ,ilist(1,idxerr+1))
            endif
         endif
         idxerr=ilist(2,idxerr)
         go to 2200
      enddo
      return
c     
      end
