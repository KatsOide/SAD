      function tfmap(isp1,mode,ihead,irtc) result(kx)
      use tfstk
      implicit none
      type (sad_descriptor) kx,kf,kl
      type (sad_descriptor) tflevelstk,tfmap1
      integer*4 maxind
      parameter (maxind=4096)
      real*8 rind(maxind)
      integer*4 isp1,narg,irtc,n1,n2,ispf,
     $     isp0,ind,mode,ihead,itfmessage
      narg=isp-isp1
      if(narg .eq. 2)then
        if(mode .eq. 1 .or. mode .eq. 3)then
          kx=tfmap1(isp1,mode,irtc)
          return
        endif
        n1=1
        n2=1
        ispf=isp-1
      elseif(narg .eq. 3)then
        call tflevelspec(dtastk(isp),n1,n2,irtc)
        if(irtc .ne. 0)then
          return
        endif
        ispf=isp-2
      elseif(narg .eq. 1)then
        if(mode .ne. 0 .and. mode .ne. 4)then
          irtc=-1
        else
          irtc=itfmessage(9,'General::narg','"2 or 3"')
        endif
        return
      else
        if(mode .ne. 0 .and. mode .ne. 4)then
          irtc=itfmessage(9,'General::narg','"2 or 3"')
        else
          irtc=itfmessage(9,'General::narg','"1 or 2 or 3"')
        endif
        return
      endif
      kf=dtfcopy(dtastk(ispf))
      kl=dtfcopy(dtastk(ispf+1))
      ind=0
      isp0=isp
      kx=tflevelstk(kl,kf,
     $     n1,n2,mode,ind,rind,ihead,mstk,irtc)
      isp=isp0
      call tflocald(kf)
      call tflocald(kl)
      if(mode .eq. 1)then
        kx%k=ktfoper+mtfnull
        if(irtc .eq. -3 .or. irtc .eq. -2)then
          irtc=0
        endif
      endif
      return
      end

      function tfmap1(isp1,mode,irtc) result(kx)
      use tfstk
      use efun
      implicit none
      type (sad_descriptor) kx,k,ki,kf
      type (sad_dlist), pointer :: kl,kli
      integer*4 isp1,isp0,isp2,mode,irtc,l,i,itfdownlevel
      k=dtastk(isp)
      if(ktfnonlistq(k,kl))then
        if(mode .eq. 1)then
          kx%k=ktfoper+mtfnull
        else
          kx=k
        endif
        irtc=0
        return
      endif
      isp0=isp
      k=dtfcopy1(k)
      kf=dtfcopy(dtastk(isp1+1))
      if(mode .eq. 1)then
        isp2=isp+1
        do i=1,kl%nl
          dtastk(isp2)=kf
          isp=isp2+1
          dtastk(isp)=kl%dbody(i)
          levele=levele+1
          ki=tfefunref(isp2,.true.,irtc)
          l=itfdownlevel()
          if(irtc .ne. 0)then
            if(irtc .eq. -3)then
              go to 1000
            elseif(irtc .eq. -2)then
              irtc=0
            else
              go to 9000
            endif
          endif
        enddo
 1000   call tflocald(kf)
        call tflocal1(k%k)
        kx%k=ktfoper+mtfnull
        isp=isp2-1
        irtc=0
        return
      else
        do i=1,kl%nl
          isp=isp+1
          isp2=isp
          dtastk(isp)=kf
          isp=isp+1
          dtastk(isp)=kl%dbody(i)
          levele=levele+1
          ki=tfefunref(isp2,.true.,irtc)
          call tfconnect(ki,irtc)
          if(irtc .ne. 0)then
            go to 9000
          endif
          isp=isp2
          dtastk(isp)=ki
        enddo
        isp=isp+1
        dtastk(isp)=kl%head
        isp2=isp
        do i=isp0+1,isp2-1
          isp=isp+1
          ki=dtastk(i)
          if(ktflistq(ki,kli))then
            if(kli%head%k .eq. ktfoper+mtfnull)then
              isp=isp-1
              call tfgetllstkall(kli)
            else
              dtastk(isp)=ki
            endif
          else
            dtastk(isp)=ki
          endif
        enddo
        kx=tfefunref(isp2,.true.,irtc)
      endif
 9000 call tflocald(kf)
      call tflocal1(k%k)
      isp=isp0
      return
      end

      subroutine tfmapall(isp1,kx,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) kx,kf,kl
      type (sad_descriptor) tflevelstk
      integer*4 isp1,narg,irtc,ihead,itfmessage
      narg=isp-isp1
      ihead=1
      if(narg .eq. 3)then
        call tfgetoption('Heads',dtastk(isp),kx,irtc)
        if(irtc .eq. 0 .and. ktftrueq(kx%k))then
          ihead=0
        else
          irtc=itfmessage(9,'General::wrongopt',' ')
          return
        endif
      elseif(narg .ne. 2)then
        irtc=itfmessage(9,'General::narg','"2 (and options)"')
        return
      endif
      kf=dtfcopy(dtastk(isp1+1))
      kl=dtfcopy(dtastk(isp1+2))
      kx=tflevelstk(kl,kf,
     $     0,maxgeneration,3,0,(/0.d0/),ihead,mstk,irtc)
      call tflocald(kf)
      call tflocald(kl)
      return
      end

      subroutine tfapply(isp1,kx,irtc)
      use tfstk
      use efun
      implicit none
      type (sad_descriptor) kx,kf,kl
      type (sad_descriptor) tflevelstk
      type (sad_dlist), pointer :: klx
      integer*4 isp1,narg,irtc,n1,n2,ispf,itfmessage
      narg=isp-isp1
      if(narg .eq. 2)then
        kx=dtastk(isp)
        if(ktflistq(kx,klx))then
          isp=isp-1
          ispf=isp
          call tfgetllstkall(klx)
          levele=levele+1
          kx=tfefunref(ispf,.true.,irtc)
          call tfconnect(kx,irtc)
          isp=ispf
        else
          irtc=0
        endif
        return
      elseif(narg .eq. 3)then
        call tflevelspec(dtastk(isp),n1,n2,irtc)
        if(irtc .ne. 0)then
          return
        endif
        ispf=isp-2
      elseif(narg .eq. 1)then
        irtc=-1
        return
      else
        irtc=itfmessage(9,'General::narg','"2 or 3"')
        return
      endif
      kf=dtfcopy(dtastk(ispf))
      kl=dtfcopy(dtastk(ispf+1))
      kx=tflevelstk(kl,kf,n1,n2,2,0,[0.d0],1,mstk,irtc)
      call tflocald(kf)
      call tflocald(kl)
      return
      end

      subroutine tfcases(isp1,kx,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) kx,kr
      type (sad_dlist), pointer :: klx
      integer*4 isp1,irtc,narg,isp0,i,itfmessage,itr
      narg=isp-isp1
      if(narg .eq. 1)then
        irtc=-1
        return
      elseif(narg .lt. 2 .or. narg .gt. 4)then
        irtc=itfmessage(9,'General::narg','"2, 3, or 4"')
        return
      endif
      kr=dtastk(isp1+2)
      if(tfruleq(kr))then
        ktastk(isp1+2)=klist(ktfaddrd(kr)+1)
        itr=1
      else
        itr=-1
      endif
      call tfposition(isp1,kx,1,irtc)
      if(itr .lt. 0 .or. irtc .ne. 0)then
        return
      endif
      if(ktfnonlistq(kx,klx))then
        return
      endif
      isp0=isp
      do i=1,klx%nl
        isp=isp+1
        call tfreplace(klx%dbody(i),kr,ktastk(isp),
     $       .true.,.true.,.false.,irtc)
        if(irtc .ne. 0)then
          go to 100
        endif
      enddo
      kx=kxmakelist(isp0)
 100  isp=isp0
      return
      end

      subroutine tfposition(isp1,kx,icases,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) kf,kx,ki,kl
      type (sad_descriptor) tfgetoption1,tflevelstk
      type (sad_dlist), pointer :: kla,klx,kll
      type (sad_rlist), pointer :: klir
      integer*4 maxind
      parameter (maxind=4096)
      real*8 rind(maxind)
      integer*4 isp1,narg,irtc,n1,n2,ispf,ind,isp0,
     $     ihead,ispmax,icases,itfmessage,ispa,
     $     i,ii,nl,ispb,itfpmat,mstk0,iop
      logical*4 rep
      integer*4, parameter :: maxlevel=100000000
      type (sad_descriptor), save :: kxheads
      data kxheads%k /0/
      if(kxheads%k .eq. 0)then
        kxheads=kxsymbolz('Heads',5)
      endif
      narg=isp-isp1
      if(narg .eq. 1)then
        irtc=-1
        return
      endif
      if(icases .eq. 0)then
        ihead=0
      else
        ihead=1
      endif
      isp0=isp
      ispa=isp
 1    if(narg .eq. 2)then
        if(icases .eq. 0)then
          n1=0
          n2=maxlevel
        else
          n1=1
          n2=1
        endif
        ispmax=mstk
        ispf=ispa
      else
        if(tfruleq(dtastk(ispa),kla))then
          kx=tfgetoption1(kxheads,kla,rep)
          if(rep .and. ktfrealq(kx))then
            ispa=ispa-1
            if(kx%k .eq. 0)then
              ihead=1
            else
              ihead=0
            endif
            narg=narg-1
            if(narg .eq. 2)then
              go to 1
            endif
          endif
        endif
        if(narg .eq. 3 .or. narg .eq. 4)then
          ispf=ispa-narg+2
          call tflevelspec(dtastk(ispf+1),n1,n2,irtc)
          if(irtc .ne. 0)then
            return
          endif
          if(narg .eq. 3)then
            ispmax=mstk
          elseif(ktfrealq(ktastk(ispa)))then
            ispb=isp+int(rtastk(ispa))
            if(n1 .eq. 1 .and. n2 .eq. 1 .and. icases .ne. 2)then
              if(ispb .eq. isp+1)then
                kl=dtastk(ispf-1)
                if(ktflistq(kl,kll))then
                  kf=dtastk(ispf)
                  nl=kll%nl
                  if(.not. ktflistq(kf) .and. .not. ktfpatq(kf))then
                    if(ihead .eq. 0)then
                      ki=kll%head
                      if(tfsameq(ki,kf))then
                        ii=0
                        go to 200
                      endif
                    endif
                    if(ktfreallistq(kll))then
                      if(ktfrealq(kf))then
                        do i=1,nl
                          if(kll%dbody(i)%k .eq. kf%k)then
                            ii=i
                            go to 200
                          endif
                        enddo
                      endif
                    else
                      LOOP_I: do i=1,nl
                        if(tfsameq(kll%dbody(i),kf))then
                          ii=i
                          go to 200
                        endif
                      enddo LOOP_I
                    endif
                  else
                    kf=dtfcopy(kf)
                    mstk0=mstk
                    iop=iordless
                    iordless=0
                    call tfinitpat(isp0,kf%k)
                    if(isp .eq. isp0)then
                      do i=ihead,nl
                        if(itfpmat(kll%dbody(i),kf) .ge. 0)then
                          ii=i
                          go to 210
                        endif
                      enddo
                    else
                      isp=isp0
                      do i=ihead,nl
                        if(itfpmat(kll%dbody(i),kf) .ge. 0)then
                          ii=i
                          go to 220
                        endif
                      enddo
                      call tfresetpat(kf)
                    endif
                    mstk=mstk0
                    isp=isp0
                    iordless=iop
                    call tflocal(kf%k)
                  endif
                  kx=dxnulll
                  irtc=0
                  return
 220              call tfresetpat(kf)
 210              mstk=mstk0
                  isp=isp0
                  iordless=iop
                  call tflocal(kf%k)
 200              if(icases .eq. 0)then
                    ki=kxavaloc(0,1,klir)
                    klir%rbody(1)=dble(ii)
                  else
                    ki=dtfcopy(kll%dbody(ii))
                  endif
                  kx=kxadaloc(-1,1,klx)
                  klx%dbody(1)=ki
                  irtc=0
                  return
                endif   
              endif      
            endif
            ispmax=max(min(ispb,mstk),isp)
          else
            irtc=itfmessage(9,'General::wrongtype',
     $           '"Real number"')
            return
          endif
        else
          irtc=itfmessage(9,'General::narg',
     $         '"2, 3, or 4 (+ options)"')
          return
        endif
      endif
      kf=dtfcopy(dtastk(ispf))
      kl=dtfcopy(dtastk(ispf-1))
      ind=0
      kx=tflevelstk(kl,kf,
     $     n1,n2,5+icases,ind,rind,ihead,ispmax,irtc)
      call tflocald(kf)
      call tflocald(kl)
      if(irtc .ne. 0)then
        return
      endif
      if(icases .eq. 2)then
        call tfeevalref(kx,kx,irtc)
      else
        kx=kxmakelist(isp0)
      endif
      isp=isp0
      return
      end

      subroutine tfflatten(isp1,kx,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) kx,kh,kh0
      type (sad_dlist), pointer :: klx,kl
      integer*4 isp1,narg,level,irtc,itfmessage
      real*8 amaxl
      parameter (amaxl=1.d9)
      narg=isp-isp1
      irtc=-1
      if(narg .gt. 3)then
        irtc=itfmessage(9,'General::narg','"Less than 4"')
        return
      endif
      if(ktfnonlistq(ktastk(isp1+1),kl))then
        irtc=itfmessage(9,'General::wrongtype','"List or composition"')
        return
      endif
      kh=kl%head
      kh0=kh
      if(narg .eq. 1)then
        level=-1
      else
        if(ktfnonrealq(dtastk(isp1+2)))then
          irtc=itfmessage(9,'General::wrongtype','"Real number"')
          return
        endif
        level=int(min(amaxl,rtastk(isp1+2)))
        if(level .le. 0)then
          irtc=itfmessage(9,'General::wrongnum','"Positive"')
          return
        endif
        if(narg .eq. 3)then
          kh=dtastk(isp)
        endif
      endif
      call tfflattenstk(kl,level,kh,irtc)
      kx=kxmakelist(isp1+narg,klx)
      klx%head=dtfcopy(kh0)
      isp=isp1+narg
      return
      end

      subroutine tfflattenstk(kl,level0,kh,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) kh
      type (sad_dlist), intent(in):: kl
      type (sad_dlist), pointer :: list,kli
      integer*4 level,irtc,i,m,i0,level0,itfmessage,mstk0
c      write(*,*)'tfflattenstk ',sad_loc(kl%head)
      call dlist_dlist(kl,list)
      mstk0=mstk
      level=level0
      i0=1
      irtc=0
 1    m=list%nl
      if(isp+m-i0+1 .gt. mstk)then
        mstk=mstk0
        irtc=itfmessage(999,'General::stack',' ')
        return
      endif
      if(ktfreallistq(list))then
        dtastk(isp+i0:isp+m)=list%dbody(i0:m)
        isp=isp+m-i0+1
      else
        LOOP_I: do i=i0,m
          isp=isp+1
          dtastk(isp)=list%dbody(i)
          if(ktflistq(ktastk(isp),kli))then
            if(tfsameq(kli%head,kh))then
              if(level .ne. 0)then
                go to 100
              endif
              cycle LOOP_I
            endif
            if(kli%head%k .eq. ktfoper+mtfnull)then
              go to 100
            endif
          endif
        enddo LOOP_I
      endif
      if(mstk .ne. mstk0)then
        mstk=mstk+1
        call loc_dlist(ktastk(mstk),list)
        i0=itastk2(1,mstk)+1
        level=level+1
        go to 1
      endif
      return
 100  isp=isp-1
      ktastk(mstk)=ksad_loc(list%head%k)
      itastk2(1,mstk)=i
      mstk=mstk-1
      level=level-1
      list=>kli
      i0=1
      go to 1
      end

      subroutine tfiff(isp1,kx,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) kx
      integer*4 isp1,irtc,narg,itfmessage,j
      narg=isp-isp1
      if(narg .le. 1 .or. narg .gt. 4)then
        irtc=itfmessage(9,'General::narg','"2, 3, or 4"')
        return
      endif
      if(ktfnonrealq(dtastk(isp1+1)))then
        if(narg .ne. 4)then
          irtc=-1
          return
        else
          j=isp1+4
        endif
      else
        if(ktastk(isp1+1) .ne. 0)then
          j=isp1+2
        elseif(narg .lt. 3)then
          kx%k=ktfoper+mtfnull
          irtc=0
          return
        else
          j=isp1+3
        endif
      endif
      call tfeevalref(dtastk(j),kx,irtc)
      return
      end

      subroutine tfswitch(isp1,kx,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) kx,kxi
      integer*4 isp1,irtc,i,itfpmatc,itfmessage
      if(mod(isp-isp1,2) .eq. 0)then
        irtc=itfmessage(9,'General::narg','"odd number"')
        return
      endif
      do i=isp1+2,isp-1,2
        call tfeevalref(dtastk(i),kxi,irtc)
        if(irtc .ne. 0)then
          return
        endif
        if(itfpmatc(dtastk(isp1+1),kxi) .ge. 0)then
          call tfeevalref(dtastk(i+1),kx,irtc)
          return
        endif
      enddo
      irtc=-1
      return
      end

      subroutine tfwhich(isp1,kx,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) kx,kxi
      integer*4 isp1,irtc,i,itfmessage
      if(mod(isp-isp1,2) .ne. 0 .or. isp .le. isp1)then
        irtc=itfmessage(9,'General::narg','"even number"')
        return
      endif
      do i=isp1+1,isp-1,2
        call tfeevalref(dtastk(i),kxi,irtc)
        if(irtc .ne. 0)then
          return
        endif
        if(ktftrueq(kxi%k))then
          call tfeevalref(dtastk(i+1),kx,irtc)
          return
        endif
      enddo
      irtc=-1
      return
      end

      subroutine tfwhile(isp1,kx,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) kx,kr,kc,ke
      integer*4 isp1,itfdownlevel,l,itgetfpe,itfmessage,irtc
      logical*4 f
      if(isp .ne. isp1+2)then
        irtc=itfmessage(9,'Genearl::narg','"2"')
        return
      endif
      f=.true.
      kc=dtfcopy(dtastk(isp1+1))
      ke=dtfcopy(dtastk(isp))
      do while(f)
        levele=levele+1
        call tfeevalref(kc,kr,irtc)
        if(irtc .ne. 0)then
          go to 9000
        endif
        f=ktftrueq(kr%k)
        if(f)then
          call tfeevalref(ke,kx,irtc)
          if(irtc .ne. 0)then
            if(irtc .eq. -3)then
              irtc=0
              go to 9000
            elseif(irtc .eq. -2)then
              irtc=0
            else
              go to 9000
            endif
          endif
          if(itgetfpe() .ne. 0)then
            call tclrfpe
            irtc=itfmessage(9,'General::fpe','""')
            go to 9000
          endif
        endif
        l=itfdownlevel()
      enddo
      call tflocald(ke)
      call tflocald(kc)
      kx%k=ktfoper+mtfnull
      return
 9000 l=itfdownlevel()
      call tflocald(ke)
      call tflocald(kc)
      kx%k=ktfoper+mtfnull
      return
      end

      subroutine tfbreak(mode,narg,kx,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) kx
      integer*4 irtc,mode,narg,itfmessage
      if(narg .gt. 1)then
        irtc=itfmessage(9,'General::narg','"0"')
        return
      elseif(narg .eq. 1 .and.
     $       ktastk(isp) .ne. ktfoper+mtfnull)then
        irtc=itfmessage(9,'General::narg','"0"')
        return
      endif
      kx%k=ktfoper+mtfnull
      irtc=mode
      return
      end

      subroutine tfthrow(mode,k,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) k
      integer*4 irtc,mode,itfmessage
      if(modethrow .ne. -1)then
        irtc=itfmessage(999,'General::throwinthrow','""')
      else
        call tflocal(kerror)
        kerror=ktfcopy(k%k)
        irtc=mode
        modethrow=mode
      endif
      return
      end

      subroutine tfcatch(isp1,kx,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) kx
      integer*4 isp1,irtc,itfmessage
      if(isp-isp1 .ne. 1)then
        irtc=itfmessage(9,'General::narg','"1"')
        return
      endif
      call tfeevalref(dtastk(isp),kx,irtc)
      call tfcatchreturn(irtcthrow,kx,irtc)
      return
      end

      subroutine tfcatchreturn(mode,kx,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) kx
      integer*4 mode,irtc
      if(irtc .lt. -1)then
c        write(*,*)'tfcatchreturn ',mode,modethrow
        if(modethrow .eq. mode)then
          modethrow=-1
          kx%k=kerror
          call tflocal(kerror)
          kerror=0
          irtc=0
        endif
      endif
      return
      end

      subroutine tfselect(isp1,kx,irtc)
      use tfstk
      use efun
      implicit none
      type (sad_descriptor) kx,k,kf,ki,kxi
      type (sad_dlist), pointer ::kl
      integer*4 isp1,irtc,narg, m,nmax,n,
     $ isp0,isp2,i,isp00, itfmessage
      real*8 amaxl
      parameter (amaxl=1.d10)
      narg=isp-isp1
      if(narg .eq. 1)then
        irtc=-1
        return
      elseif(narg .eq. 0 .or. narg .gt. 3)then
        irtc=itfmessage(9,'General::narg','"1 or 2 or 3"')
        return
      endif
      k=dtastk(isp1+1)
      if(.not. ktflistq(k,kl))then
        irtc=itfmessage(9,'General::wrongtype','"List or composition"')
        return
      endif
      irtc=-1
      m=kl%nl
      isp00=isp
      if(narg .eq. 3)then
        if(ktfnonrealq(ktastk(isp)))then
          irtc=itfmessage(9,'General::wrongtype','"Real number"')
          return
        endif
        nmax=int(min(amaxl,rtastk(isp)))
        if(nmax .lt. 0)then
          return
        endif
      else
        nmax=m
        isp=isp+1
      endif
      irtc=0
      if(nmax .eq. 0)then
        kx=dxnulll
        isp=isp00
        return
      endif
      k=dtfcopy1(k)
      n=0
      kf=dtfcopy(dtastk(isp1+2))
      isp0=isp
      do i=1,m
        ki=kl%dbody(i)
        isp2=isp
        isp=isp+1
        dtastk(isp)=kf
        isp=isp+1
        dtastk(isp)=ki
        levele=levele+1
        kxi=tfefunref(isp2+1,.true.,irtc)
        call tfconnect(kxi,irtc)
        if(irtc .ne. 0)then
          go to 9000
        endif
        isp=isp2
        if(ktftrueq(kxi%k))then
          isp=isp+1
          dtastk(isp)=ki
          if(isp .ge. isp0+nmax)then
            go to 9100
          endif
        endif
      enddo
 9100 kx=kxmakelist(isp0)
 9000 call tflocald(kf)
      call tflocal1d(k)
      isp=isp00
      return
      end

      subroutine tfswitchcases(isp1,kx,mode,irtc)
      use tfstk
      use efun
      implicit none
      type (sad_descriptor) kx
      type (sad_dlist), pointer :: kl,kl2,klx
      integer*8 kak(4096)
      integer*4 isp1,irtc,mk(0:4096),km,mode,isp4,
     $     i,j,isp2,isp3,m,itfpmatc,kelm,isp0,itfmessage
      if(isp .eq. isp1+1)then
        irtc=-1
        return
      elseif(isp .ne. isp1+2)then
        irtc=itfmessage(9,'General::narg','"2"')
        return
      endif
      isp2=isp1+1
      if(ktfnonlistq(dtastk(isp2),kl2))then
        irtc=itfmessage(9,'General::wrongtype',
     $       '"List or composition for #1"')
        return
      endif
      if(tflistq(dtastk(isp),kl))then
        isp=isp2
        call tfgetllstkall(kl)
      else
        irtc=itfmessage(9,'General::wrongtype',
     $       '"List for #2"')
        return
      endif
      kelm=isp-isp2
      if(kelm .eq. 0)then
        kx=dxnulll
        irtc=0
        return
      endif
      mk(0:kelm)=0
      isp0=isp
      call tfgetllstkall(kl2)
      isp3=isp
      m=isp-isp0
      if(mode .eq. 0)then
        do i=isp0+1,isp
          do j=1,kelm
            if(itfpmatc(dtastk(i),dtastk(isp2+j)) .ge. 0)then
              km=j
              go to 10
            endif
          enddo
          km=0
 10       isp=isp+1
          itastk(1,isp)=km
          mk(km)=mk(km)+1
        enddo
      else
        do i=isp0+1,isp
          do j=1,kelm
            isp4=isp
            isp=isp+1
            ktastk(isp)=ktastk(isp2+j)
            isp=isp+1
            ktastk(isp)=ktastk(i)
            kx=tfefunref(isp4+1,.true.,irtc)
            if(irtc .ne. 0)then
              isp=isp2+1
              return
            endif
            isp=isp4
            if(ktftrueq(kx%k))then
              km=j
              go to 20
            endif
          enddo
          km=0
 20       isp=isp+1
          itastk(1,isp)=km
          mk(km)=mk(km)+1
        enddo
      endif
      kx=kxadaloc(-1,kelm,klx)
      do km=1,kelm
        if(mk(km) .eq. 0)then
          kak(km)=ktfaddr(ktfcopy1(kxnulll))
        else
          kak(km)=ktaaloc(0,mk(km))
        endif
        klx%dbody(km)%k=ktflist+kak(km)
        mk(km)=1
      enddo
      do i=1,m
        km=itastk(1,isp3+i)
        if(km .gt. 0)then
          call tfsetlist(ktfcopy(ktastk(isp0+i)),kak(km),mk(km))
          mk(km)=mk(km)+1
        endif
      enddo
      isp=isp2+1
      irtc=0
      return
      end

      subroutine tfnames(isp1,kx,irtc)
      use tfstk
      use tfcode
      use iso_c_binding
      implicit none
      type (sad_descriptor) kx
      type (sad_namtbl),pointer :: loc
      integer*8 i,j
      integer*4 isp1,irtc,isp0,npat,nv,l,itfmessage,n
      character*256 pat,tfgetstr,vname
      logical*4 tmatch
      isp0=isp
      do n=isp1+1,isp
        if(ktfnonstringq(ktastk(n)))then
          irtc=itfmessage(9,'General::wrongtype','"Character-string"')
          isp=isp0
          return
        endif
        pat=tfgetstr(dtastk(n),npat)
        if(pat .eq. '*')then
          do l=0,nsymhash
            j=itfcontext+l+1
            i=klist(j)
            do while(i .ne. j)
              call loc_namtbl(i,loc)
              nv=loc%str%nch
              vname=loc%str%str(1:nv)
              if(loc%symdef .ne. 0)then
                isp=isp+1
                dtastk(isp)=loc%str%alloc
              endif
              i=klist(i)
            enddo
          enddo
        else
          do l=0,nsymhash
            j=itfcontext+l+1
            i=klist(j)
            do while(i .ne. j)
              call loc_namtbl(i,loc)
              nv=loc%str%nch
              vname=loc%str%str(1:nv)
              if(loc%symdef .ne. 0 .and.
     $             tmatch(vname(1:nv),pat(1:npat)))then
                isp=isp+1
                dtastk(isp)=loc%str%alloc
              endif
              i=klist(i)
            enddo
          enddo
        endif
      enddo
      kx=kxmakelist(isp0)
      isp=isp0
      irtc=0
      return
      end
