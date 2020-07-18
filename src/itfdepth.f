      recursive integer*4 function itfdepth(k) result(id)
      use tfstk
      implicit none
      type (sad_descriptor) k
      type (sad_dlist), pointer ::kl
      integer*4 i
      if(ktflistq(k,kl))then
        if(ktfreallistq(kl))then
          id=2
        else
          id=2
          do i=1,kl%nl
            if(ktflistq(kl%dbody(i)))then
              id=max(id,itfdepth(kl%dbody(i))+1)
            endif
          enddo
        endif
      else
        id=1
      endif
      return
      end
      
      recursive function tflevelstk(k,kf,
     $     n1,n2,mode,ind,rind,ihead,ispmax,irtc)
     $     result(kx)
      use tfstk
      use efun
      implicit none
      type (sad_descriptor) kx
      type (sad_descriptor) ,intent(in):: k,kf
      type (sad_descriptor) kxi,ki
      type (sad_dlist), pointer :: klx,kl
      integer*4 ,intent(in):: mode,n1,n2,ispmax,ihead
      integer*4 ,intent(out):: irtc
      integer*4 i,isp1,m,idi,ispf,itfdepth,id,ind,ind1,
     $     itfpmatc
      real*8 ,intent(out)::rind(ind)
      logical*4 stack,map,indf
      integer*4 ,parameter::maxlevel=100000000
      logical*4 ,parameter:: stacktbl(0:7)=[
     $     .false.,.false., .true., .true., .true.,
     $     .false.,.false., .true.],
     $ indexf(0:7)=[
     $     .false.,.false.,.false.,.false., .true.,
     $      .true.,.false.,.false.],
     $ match(0:7)=[
     $     .false.,.false.,.false.,.false.,.false.,
     $      .true.,.true., .true.]
c            Level    Scan   Apply     Map MapIndx
c          Positio   Cases DeCases     
      integer*4 ,parameter ::ioff(0:7)=[0,1,1,1,1,0,0,0]
      irtc=0
      if(isp .ge. ispmax)then
        return
      endif
      stack=stacktbl(mode)
      indf=indexf(mode)
      ispf=isp+1
      isp=isp+ioff(mode)
      ind1=ind
c      call tfdebugprint(k,'tflevelstk',1)
c      write(*,*)mode,n1,n2,ihead,indf
      if(n2 .ne. 0 .and. ktflistq(k,kl))then
        isp1=isp
        m=kl%nl
        if(stack .and. ihead .eq. 1)then
          isp=isp+1
          dtastk(isp)=kl%head
        endif
        if(isp .ge. ispmax)then
          go to 4000
        endif
        if(n1 .ge. 0)then
          if(n2 .gt. 0)then
            if(ihead .eq. 0)then
              if(indf)then
                rind(ind+1)=0.d0
                kxi=tflevelstk(kl%head,kf,
     $               max(0,n1-1),n2-1,mode,ind+1,rind,
     $               0,ispmax,irtc)
              else
                kxi=tflevelstk(kl%head,kf,
     $               max(0,n1-1),n2-1,mode,0,rind,0,ispmax,irtc)
              endif
            endif
            if(isp .ge. ispmax)then
              go to 4000
            endif
            if(indf)then
              do i=1,m
                rind(ind+1)=i
                kxi=tflevelstk(kl%dbody(i),kf,
     $                 max(0,n1-1),n2-1,mode,ind+1,rind,
     $                 ihead,ispmax,irtc)
                if(irtc .ne. 0)then
                  go to 9000
                endif
                if(isp .ge. ispmax)then
                  go to 4000
                endif
              enddo
            else
              do i=1,m
                kxi=tflevelstk(kl%dbody(i),kf,
     $                 max(0,n1-1),n2-1,mode,0,rind,
     $                 ihead,ispmax,irtc)
                if(irtc .ne. 0)then
                  go to 9000
                endif
                if(isp .ge. ispmax)then
                  go to 4000
                endif
              enddo
            endif
          else
            ind1=ind+1
            if(ihead .eq. 0)then
              idi=itfdepth(kl%head)
              if(idi .ge. -n2)then
                if(indf)then
                  rind(ind1)=0.d0
                endif
                kxi=tflevelstk(kl%head,kf,
     $             max(0,n1-1),n2,mode,ind1,rind,ihead,
     $               ispmax,irtc)
                if(irtc .ne. 0)then
                  go to 9000
                endif
              elseif(stack)then
                isp=isp+1
                dtastk(isp)=kl%head
              endif
            endif
            if(isp .ge. ispmax)then
              go to 4000
            endif
            do i=1,m
              ki=kl%dbody(i)
              if(ktflistq(ki))then
                idi=itfdepth(ki)
              else
                idi=1
              endif
              if(idi .ge. -n2)then
                if(indf)then
                  rind(ind1)=i
                endif
                kxi=tflevelstk(ki,kf,
     $               max(0,n1-1),n2,mode,ind1,rind,
     $               ihead,ispmax,irtc)
                if(irtc .ne. 0)then
                  go to 9000
                endif
              elseif(stack)then
                isp=isp+1
                dtastk(isp)=ki
              endif
              if(isp .ge. ispmax)then
                go to 4000
              endif
            enddo
          endif
        else
          ind1=ind+1
          if(ihead .eq. 0)then
            ki=kl%head
            idi=itfdepth(ki)
            if(indf)then
              rind(ind1)=0.d0
            endif
            kxi=tflevelstk(ki,kf,
     $           max(0,idi+n1),min(n2,abs(n2)-1),
     $           mode,ind1,rind,ihead,ispmax,irtc)
            if(irtc .ne. 0)then
              go to 9000
            endif
            if(isp .ge. ispmax)then
              go to 4000
            endif
          endif
          do i=1,m
            ki=kl%dbody(i)
            if(ktflistq(ki))then
              idi=itfdepth(ki)
            else
              idi=1
            endif
            if(indf)then
              rind(ind1)=i
            endif
            kxi=tflevelstk(ki,kf,
     $           max(0,idi+n1),min(n2,abs(n2)-1),
     $           mode,ind1,rind,ihead,ispmax,irtc)
            if(irtc .ne. 0)then
              go to 9000
            endif
            if(isp .ge. ispmax)then
              exit
            endif
          enddo
        endif
 4000   if(stack)then
          kx=kxcompose(isp1+1)
          isp=isp1
        else
          kx=k
        endif
      else
        kx=k
      endif
      if(n1 .eq. 0)then
        map=.true.
        if(n2 .lt. 0)then
          id=itfdepth(k)
        else
          id=0
        endif
      elseif(n1 .lt. 0)then
        id=itfdepth(k)
        map=id+n1 .le. 0
      else
        map=.false.
      endif
      if(map)then
        if(n2 .ge. 0 .or. id .ge. -n2)then
          isp=isp+1
          dtastk(isp)=kx
          if(match(mode))then
            if(itfpmatc(kx%k,kf) .ge. 0)then
              if(mode .eq. 5)then
                if(ind .eq. 0)then
                  dtastk(isp)=dxnulll
                else
                  dtastk(isp)=kxm2l(rind,0,ind,1,.false.)
                endif
                return
              elseif(mode .eq. 6)then
                return
              endif
              kx=dxnull
            endif
            if(mode .ne. 7)then
              isp=isp-1
              return
            endif
          elseif(mode .eq. 3)then
            dtastk(ispf)=kf
            call tfefunrefc(ispf,kx,irtc)
          elseif(mode .eq. 2)then
            if(ktflistq(kx,klx))then
              isp=ispf
              call tfgetllstkall(klx)
              dtastk(ispf)=kf
              call tfefunrefc(ispf,kx,irtc)
            endif
          elseif(mode .eq. 4)then
            dtastk(ispf)=kf
            isp=isp+1
            if(ind .eq. 0)then
              dtastk(isp)=dxnulll
            else
              dtastk(isp)=kxm2l(rind,0,ind,1,.false.)
            endif
            call tfefunrefc(ispf,kx,irtc)
          elseif(mode .eq. 1)then
            dtastk(ispf)=kf
            call tfefunrefc(ispf,kx,irtc)
            isp=ispf-1
            if(irtc .eq. -2)then
              irtc=0
            endif
            return
          else
            return
          endif
        endif
      endif
      if(stack)then
        dtastk(ispf)=kx
        isp=ispf
      elseif(mode .eq. 1)then
        isp=ispf-1
      endif
      return
 9000 if(stack)then
        isp=ispf
      elseif(mode .eq. 1)then
        isp=ispf-1  
      endif
      return
      end

      subroutine tflevel(k,kl,kx,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) k,kl,kx
      type (sad_descriptor) tflevelstk
      integer*4 n1,n2,irtc,isp1
      real*8 rind(1)
      call tflevelspec(kl,n1,n2,irtc)
      if(irtc .ne. 0)then
        return
      endif
      if(n2 .ge. 0 .and. n1 .gt. n2 .or.
     $     n1 .lt. 0 .and. n1 .gt. n2)then
        kx=dxnulll
        return
      endif
c      write(*,*)'tflevel ',n1,n2
      isp1=isp
      kx=tflevelstk(k,sad_descr(ktfoper+mtfnull),n1,n2,
     $     0,0,rind,1,mstk,irtc)
      if(isp .le. isp1)then
        kx=dxnulll
      else
        kx=kxmakelist(isp1)
      endif
      isp=isp1
      irtc=0
      return
      end

      subroutine tflevelspec(k,n1,n2,irtc)
      use tfstk
      implicit none
      type (sad_rlist), pointer :: kl
      type (sad_descriptor) k
      integer*4 n1,n2,irtc,ivl,maxlevel,m,itfmessage
      parameter (maxlevel=2**30)
      real*8 vlmax,v
      parameter (vlmax=1.d8)
      if(ktfrealq(k,v))then
        ivl=int(max(-vlmax,min(vlmax,v)))
        n1=1
        n2=ivl
      elseif(.not. tfreallistq(k%k,kl))then
        n1=0
        n2=0
        irtc=itfmessage(9,'General::wrongtype',
     $       '"List of numbers for levelspec"')
        return
      else
        m=kl%nl
        if(m .eq. 1)then
          n1=int(max(-vlmax,min(vlmax,kl%rbody(1))))
          n2=n1
        elseif(m .eq. 2)then
          n1=int(max(-vlmax,min(vlmax,kl%rbody(1))))
          n2=int(max(-vlmax,min(vlmax,kl%rbody(2))))
        else
          n1=0
          n2=0
          irtc=itfmessage(9,'General::wrongval',
     $         '"n, {n}, or {n1, n2}","as level spec"')
          return
        endif
      endif
      irtc=0
      return
      end
