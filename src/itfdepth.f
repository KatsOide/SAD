      recursive integer*4 function itfdepth(k) result(id)
      use tfstk
      implicit none
      type (sad_descriptor) k
      type (sad_list), pointer ::kl
      integer*4 i
      if(ktflistqd(k,kl))then
        if(ktfreallistqo(kl))then
          id=2
        else
          id=2
          do i=1,kl%nl
            if(ktflistq(kl%body(i)))then
              id=max(id,itfdepth(kl%dbody(i))+1)
            endif
          enddo
        endif
      else
        id=1
      endif
      return
      end
      
      recursive subroutine tflevelstk(k,kf,kx,
     $     n1,n2,mode,ind,rind,ihead,ispmax,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) k,kx,kxi,ki,kf
      type (sad_list), pointer :: klx,kl
      integer*4 maxlevel,mode,m,n1,n2,irtc,i,isp1,
     $     idi,ispf,itfdepth,id,ind,ind1,
     $     itfpmatc,ihead,ispmax,ioff(0:7)
      real*8 rind(ind)
      logical*4 stack,map,stacktbl(0:7),indexf(0:7),indf,
     $     match(0:7)
      parameter (maxlevel=100000000)
      data stacktbl/
     $     .false.,.false., .true., .true., .true.,
     $     .false.,.false., .true./
      data indexf/
     $     .false.,.false.,.false.,.false., .true.,
     $      .true.,.false.,.false./
      data match/
     $     .false.,.false.,.false.,.false.,.false.,
     $      .true.,.true., .true./
c            Level    Scan   Apply     Map MapIndx
c          Positio   Cases DeCases     
      data ioff/0,1,1,1,1,0,0,0/
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
      if(n2 .ne. 0 .and. ktflistqd(k,kl))then
        isp1=isp
        m=kl%nl
        if(stack .and. ihead .eq. 1)then
          isp=isp+1
          ktastk(isp)=kl%head
        endif
        if(isp .ge. ispmax)then
          go to 4000
        endif
        if(n1 .ge. 0)then
          if(n2 .gt. 0)then
            if(ihead .eq. 0)then
              if(indf)then
                rind(ind+1)=0.d0
                call tflevelstk(kl%dbody(0),kf,kxi,
     $               max(0,n1-1),n2-1,mode,ind+1,rind,
     $               0,ispmax,irtc)
              else
                call tflevelstk(kl%dbody(0),kf,kxi,
     $               max(0,n1-1),n2-1,mode,0,rind,0,ispmax,irtc)
              endif
            endif
            if(isp .ge. ispmax)then
              go to 4000
            endif
            if(indf)then
              do i=1,m
                rind(ind+1)=i
                call tflevelstk(kl%dbody(i),kf,kxi,
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
                call tflevelstk(kl%dbody(i),kf,kxi,
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
              idi=itfdepth(kl%dbody(0))
              if(idi .ge. -n2)then
                if(indf)then
                  rind(ind1)=0.d0
                endif
                call tflevelstk(kl%dbody(0),kf,kxi,
     $             max(0,n1-1),n2,mode,ind1,rind,ihead,
     $               ispmax,irtc)
                if(irtc .ne. 0)then
                  go to 9000
                endif
              elseif(stack)then
                isp=isp+1
                ktastk(isp)=kl%head
              endif
            endif
            if(isp .ge. ispmax)then
              go to 4000
            endif
            do i=1,m
              ki=kl%dbody(i)
              if(ktflistqd(ki))then
                idi=itfdepth(ki)
              else
                idi=1
              endif
              if(idi .ge. -n2)then
                if(indf)then
                  rind(ind1)=i
                endif
                call tflevelstk(ki,kf,kxi,
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
            ki=kl%dbody(0)
            idi=itfdepth(ki)
            if(indf)then
              rind(ind1)=0.d0
            endif
            call tflevelstk(ki,kf,kxi,
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
            if(ktflistqd(ki))then
              idi=itfdepth(ki)
            else
              idi=1
            endif
            if(indf)then
              rind(ind1)=i
            endif
            call tflevelstk(ki,kf,kxi,
     $           max(0,idi+n1),min(n2,abs(n2)-1),
     $           mode,ind1,rind,ihead,ispmax,irtc)
            if(irtc .ne. 0)then
              go to 9000
            endif
            if(isp .ge. ispmax)then
              go to 4000
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
            if(itfpmatc(kx%k,kf) .lt. 0)then
              go to 7010
            endif
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
 7010       if(mode .ne. 7)then
              isp=isp-1
              return
            endif
          elseif(mode .eq. 3)then
            dtastk(ispf)=kf
            call tfefunrefc(ispf,kx%k,irtc)
          elseif(mode .eq. 2)then
            if(ktflistqd(kx,klx))then
              isp=ispf
              call tfgetllstkall(klx)
              dtastk(ispf)=kf
              call tfefunrefc(ispf,kx%k,irtc)
            endif
          elseif(mode .eq. 4)then
            dtastk(ispf)=kf
            isp=isp+1
            if(ind .eq. 0)then
              dtastk(isp)=dxnulll
            else
              dtastk(isp)=kxm2l(rind,0,ind,1,.false.)
            endif
            call tfefunrefc(ispf,kx%k,irtc)
          elseif(mode .eq. 1)then
            dtastk(ispf)=kf
            call tfefunrefc(ispf,kx,.true.,irtc)
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
      integer*4 n1,n2,irtc,isp1
      real*8 rind(1)
      call tflevelspec(kl,n1,n2,irtc)
      if(irtc .ne. 0)then
        return
      endif
      if(n2 .ge. 0 .and. n1 .gt. n2)then
        kx=dxnulll
        return
      endif
      if(n1 .lt. 0 .and. n1 .gt. n2)then
        kx=dxnulll
        return
      endif
c      write(*,*)'tflevel ',n1,n2
      isp1=isp
      call tflevelstk(k,ktfoper+mtfnull,kx,n1,n2,
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
      type (sad_list), pointer :: kl
      type (sad_descriptor) k
      integer*4 n1,n2,irtc,ivl,maxlevel,m,itfmessage
      parameter (maxlevel=2**30)
      real*8 vlmax
      parameter (vlmax=1.d8)
      if(ktfrealqd(k))then
        ivl=int(max(-vlmax,min(vlmax,rfromd(k))))
        n1=1
        n2=ivl
      elseif(ktflistqd(k,kl))then
        if(ktfnonreallistqo(kl))then
          irtc=itfmessage(9,'General::wrongtype',
     $         '"List of numbers for levelspec"')
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
            irtc=itfmessage(9,'General::wrongval',
     $           '"n, {n}, or {n1, n2}","as level spec"')
            return
          endif
        endif
      else
        irtc=itfmessage(9,'General::wrongtype',
     $       '"List of numbers for levelspec"')
        return
      endif
      irtc=0
      return
      end
