      module casym
      use tfstk
      implicit none
      type (sad_descriptor) ,save :: kxcamonitor
      data kxcamonitor%k /0/
      integer*8, save::
     $     irl,icsconn,icarl,icacsconn,
     $     ics,ipos,icsl,irn,irnl,icscomm,iauto,istart,
     $     iv,isev,its,ivl,isevl,itsl,ivalcomm

      contains
        subroutine tfepicssyminit(irtc)
        use tfstk
        use eeval
        implicit none
          integer*4 irtc
        type (sad_descriptor) kax,kx
        kax=kxsymbolz('CaMonitor',9)
        kx=tfsyeval(kax,irtc)
        if(irtc .ne. 0 .or. ktfnonlistq(kx))then
          go to 9000
        endif
        kxcamonitor=dtfcopy(kx)
        irl=ktfsymbolz('rl',2)
        call tfclassmember(kxcamonitor,ktfsymbol+irl,kx%k,.false.,irtc)
        if(irtc .ne. 0 .or. ktfnonsymbolq(kx))then
          kxcamonitor%k=0
          go to 9000
        endif
        icarl=ktfaddr(kx)
        icsconn=ktfsymbolz('CS$Conn',7)
        call tfclassmember(kxcamonitor,ktfsymbol+icsconn,
     $       kx%k,.false.,irtc)
        if(irtc .ne. 0 .or. ktfnonsymbolq(kx))then
          kxcamonitor%k=0
          go to 9000
        endif
        icacsconn=ktfaddr(kx)
        ics=ktfsymbolz('cs',2)
        ipos=ktfsymbolz('pos',3)
        icsl=ktfsymbolz('csl',3)
        irn=ktfsymbolz('rn',2)
        irnl=ktfsymbolz('rnl',3)
        icscomm=ktfsymbolz('cscomm',6)
        iauto=ktfsymbolz('autostart',9)
        istart=ktfsymbolz('Start',5)
        iv=ktfsymbolz('v',1)
        isev=ktfsymbolz('sev',3)
        its=ktfsymbolz('ts',2)
        ivl=ktfsymbolz('vl',2)
        isevl=ktfsymbolz('sevl',4)
        itsl=ktfsymbolz('tsl',3)
        ivalcomm=ktfsymbolz('valcomm',7)
        return
 9000   if(irtc .gt. 0 .and. ierrorprint .ne. 0)then
          call tfreseterror
        endif
        return
          end subroutine
      end module

      subroutine tfepicsconstatcb(chid,istat)
      use tfstk
      use casym
      use eeval
      implicit none
      type (sad_descriptor) stat,kn,krnn,ki
      real*8 chid
      integer*8 kx, iastart,icarlch,kaa,
     $     ka,iacscomm,iarn
      integer*4 isp0,irtc,n,l,itfdownlevel,istat
      real*8 vn
      logical*4 ev
      levele=levele+1
      isp0=isp
      if(kxcamonitor%k .eq. 0)then
        call tfepicssyminit(irtc)
        if(irtc .ne. 0)then
          go to 9000
        endif
      endif
      call tfsetpendio
      isp=isp0+1
      ktastk(isp)=ktfsymbol+icarl
      isp=isp+1
      rtastk(isp)=chid
      call tfdeval(isp0+1,icarl,kx,1,.false.,ev,irtc)
      if(irtc .ne. 0 .or. ktfnonlistq(kx))then
        go to 9000
      endif
      icarlch=ktfaddr(kx)
      if(ilist(2,icarlch-1) .lt. 2)then
        go to 9000
      endif
      ka=klist(icarlch+1)
      kn=dlist(icarlch+2)
      if(ktfnonlistq(ka) .or. ktfnonrealq(kn,vn))then
        go to 9000
      endif
      kaa=ktfaddr(ka)
      stat%x(1)=istat
      if(vn .lt. 1.d0)then
        call tfcbsetsymbol(kaa,ics,stat,irtc)
        if(irtc .ne. 0)then
          go to 9000
        endif
      else
        call tfcbsetsymbol(kaa,ipos,kn,irtc)
        if(irtc .ne. 0)then
          go to 9000
        endif
        n=int(vn)
        call tfcbsetlist(kaa,ics,n,stat,irtc)
        if(irtc .ne. 0)then
          go to 9000
        endif
        call tfclassmember(dfromk(ktflist+kaa),ktfsymbol+irn,
     $       kx,.true.,irtc)
        if(irtc .ne. 0)then
          go to 9000
        endif
        iarn=ktfaddr(kx)
        if(ktfnonlistq(kx) .or. ilist(2,iarn-1) .lt. n)then
          go to 9000
        endif
        krnn=dlist(iarn+n)
        call tfcbsetsymbol(kaa,irnl,krnn,irtc)
        if(irtc .ne. 0)then
          go to 9000
        endif
        call tfcbsetsymbol(kaa,icsl,stat,irtc)
        if(irtc .ne. 0)then
          go to 9000
        endif
      endif
      call tfclassmember(dfromk(ktflist+kaa),ktfsymbol+icscomm,
     $     kx,.true.,irtc)
      if(irtc .ne. 0)then
        go to 9000
      endif
      iacscomm=ktfaddr(kx)
      if(ktflistq(kx))then
        if(klist(iacscomm) .eq. ktfoper+mtfhold)then
          ki=tfeevalref(dlist(iacscomm+1),irtc)
          if(irtc .gt. 0 .and. ierrorprint .ne. 0)then
            call tfreseterror
          endif
        endif
      endif
      call tfclassmember(dfromk(ktflist+kaa),ktfsymbol+iauto,
     $     kx,.true.,irtc)
      if(irtc .ne. 0)then
        go to 9000
      endif
      if(ktfrealq(kx) .and. kx .ne. 0)then
        if(rlist(icacsconn-4) .eq. stat%x(1))then
          call tfclassmember(dfromk(ktflist+kaa),ktfsymbol+istart,
     $         kx,.true.,irtc)
          if(irtc .ne. 0 .or. ktfnonsymbolq(kx))then
            go to 9000
          endif
          iastart=ktfaddr(kx)
          isp=isp0+1
          ktastk(isp)=kx
          if(vn .ge. 1.d0)then
            isp=isp+1
            rtastk(isp)=vn
          endif
          call tfdeval(isp0+1,iastart,kx,1,.false.,ev,irtc)
        endif
      endif
 9000 call tfresetpendio
      l=itfdownlevel()
      isp=isp0
      return
      end
      
      subroutine tfepicsvaluecb(chid,istat,jsev,t,itype,
     $     nc,karray)
      use tfstk
      use casym
      use eeval
      implicit none
      type (sad_descriptor) k,krnn,ka,kn,ki
      real*8 chid
      integer*8 karray(nc),kx,iarn,iavalcomm,kaa,icarlch
      integer*4 istat,jsev,itype,nc
      real*8 t,vn,stat
      integer*4 isp0,isp2,irtc,n,l,i,itfdownlevel
      logical*4 ev
      levele=levele+1
      isp0=isp
      if(kxcamonitor%k .eq. 0)then
        call tfepicssyminit(irtc)
        if(irtc .ne. 0)then
          go to 9000
        endif
      endif
      isp=isp0+1
      ktastk(isp)=ktfsymbol+icarl
      isp=isp+1
      rtastk(isp)=chid
      call tfdeval(isp0+1,icarl,kx,1,.false.,ev,irtc)
      if(irtc .ne. 0 .or. ktfnonlistq(kx))then
        go to 9000
      endif
      icarlch=ktfaddr(kx)
      if(ilist(2,icarlch-1) .lt. 2)then
        go to 9000
      endif
      ka=dlist(icarlch+1)
      kn=dlist(icarlch+2)
      if(ktfnonlistq(ka) .or. ktfnonrealq(kn))then
        go to 9000
      endif
      n=int(kn%x(1))
      stat=istat
      if(itype .eq. 0)then
        if(nc .eq. 1)then
          k%k=karray(1)
        else
          isp2=isp
          do i=1,nc
            isp=isp+1
            ktastk(isp)=karray(i)
          enddo
          k=kxmakelist(isp2)
          isp=isp2
        endif
      else
        if(nc .eq. 1) then
          k%k=karray(1)
        else
          k=kxm2l(rlist(ksad_loc(karray(1)):),0,nc,nc,.false.)
        endif
      endif
      kaa=ktfaddr(ka%k)
      if(n .lt. 1)then
        call tfcbsetsymbol(kaa,iv,k,irtc)
        if(irtc .ne. 0)then
          go to 9000
        endif
        call tfcbsetsymbol(kaa,isev,dfromr(dble(jsev)),irtc)
        if(irtc .ne. 0)then
          go to 9000
        endif
        call tfcbsetsymbol(kaa,its,dfromr(t),irtc)
        if(irtc .ne. 0)then
          go to 9000
        endif
      else
        call tfcbsetsymbol(kaa,ipos,dfromr(vn),irtc)
        if(irtc .ne. 0)then
          go to 9000
        endif
        call tfcbsetlist(kaa,iv,n,k,irtc)
        if(irtc .ne. 0)then
          go to 9000
        endif
        call tfcbsetlist(kaa,isev,n,dfromr(dble(jsev)),irtc)
        if(irtc .ne. 0)then
          go to 9000
        endif
        call tfcbsetlist(kaa,its,n,dfromr(t),irtc)
        if(irtc .ne. 0)then
          go to 9000
        endif
        call tfclassmember(dfromk(ktflist+kaa),ktfsymbol+irn,
     $       kx,.true.,irtc)
        if(irtc .ne. 0)then
          go to 9000
        endif
        iarn=ktfaddr(kx)
        if(ktfnonlistq(kx) .or. ilist(2,iarn-1) .lt. n)then
          go to 9000
        endif
        krnn=dlist(iarn+n)
        call tfcbsetsymbol(kaa,irnl,krnn,irtc)
        if(irtc .ne. 0)then
          go to 9000
        endif
        call tfcbsetsymbol(kaa,ivl,k,irtc)
        if(irtc .ne. 0)then
          go to 9000
        endif
        call tfcbsetsymbol(kaa,isevl,dfromr(dble(jsev)),irtc)
        if(irtc .ne. 0)then
          go to 9000
        endif
        call tfcbsetsymbol(kaa,itsl,dfromr(t),irtc)
        if(irtc .ne. 0)then
          go to 9000
        endif
      endif
      call tfclassmember(dfromk(ktflist+kaa),ktfsymbol+ivalcomm,
     $     kx,.true.,irtc)
      if(irtc .ne. 0)then
        go to 9000
      endif
      iavalcomm=ktfaddr(kx)
      if(ktflistq(kx))then
        if(klist(iavalcomm) .eq. ktfoper+mtfhold)then
          ki=tfeevalref(dlist(iavalcomm+1),irtc)
        endif
      endif
 9000 call tfresetpendio
      l=itfdownlevel()
      isp=isp0
      return
      end

      subroutine tfcbsetsymbol(ka,kv,k1,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) k1
      integer*8 ka,kv,kax,kx
      integer*4 irtc
      call tfclassmember(dfromk(ktflist+ka),ktfsymbol+kv,
     $     kx,.false.,irtc)
      if(irtc .ne. 0)then
        return
      elseif(ktfnonsymbolq(kx))then
        irtc=-1
        return
      endif
      kax=ktfaddr(kx)
      call tflocal(klist(kax-4))
      dlist(kax-4)=dtfcopy(k1)
      return
      end

      subroutine tfcbsetlist(ka,kv,n,k1,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) k1
      integer*8 ka,kv,kal,kax,kx
      integer*4 n,irtc,i
      call tfclassmember(dfromk(ktflist+ka),ktfsymbol+kv,
     $     kx,.false.,irtc)
      if(irtc .ne. 0)then
        return
      elseif(ktfnonsymbolq(kx))then
        go to 9000
      endif
      kax=ktfaddr(kx)
      if(ktfnonlistq(klist(kax-4)))then
        go to 9000
      endif
      kal=ktfaddr(klist(kax-4))
      if(ilist(2,kal-1) .lt. n)then
        go to 9000
      endif
      if(ktfreallistq(kal))then
        if(ktfnonrealq(k1))then
          ilist(2,kal-3)=lnonreallist
        endif
      else
        call tflocal(klist(kal+n))
        if(ktfrealq(k1))then
          if(ktfnonrealq(klist(kal+n)))then
            dlist(kal+n)=k1
            do i=1,ilist(2,kal-1)
              if(ktfnonrealq(klist(kal+i)))then
                return
              endif
            enddo
            ilist(2,kal-3)=ilist(2,kal-3)-lnonreallist
            return
          endif
        endif
      endif        
      dlist(kal+n)=dtfcopy(k1)
      return
 9000 irtc=-1
      return
      end
