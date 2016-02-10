      subroutine tffindroot(isp1,kx,irtc)
      use tfstk
      implicit none
      integer*4 nvmax,maxi0
      real*8 eps0
      parameter (nvmax=2048,maxi0=50,eps0=1.d-20)
      integer*8 kx,kav(nvmax),kav0(nvmax),ktfsymbolf,ke,kdl(nvmax),
     $     kax,ktfcopy,ktfmakelist
      integer*4 isp1,irtc,neq,nvar,itfmessage,isp2,i,maxi,ispv
      real*8 v0(nvmax),eps,vmin(nvmax),vmax(nvmax),d0,rfromk
      logical*4 trace,used
      integer*8 itfres
      data itfres /0/
      if(isp .lt. isp1+2)then
        irtc=itfmessage(9,'General::narg','"2 or more"')
        return
      endif
      if(itfres .eq. 0)then
        itfres=ktfsymbolf('Residual',8,.false.)
      endif
      maxi=maxi0
      eps=eps0
      trace=.false.
      used=.true.
      ispv=isp
      do i=isp,isp1+3,-1
        call tfgetoption('MaxIterations',ktastk(i),kx,irtc)
        if(irtc .eq. -1)then
          ispv=i
          go to 1
        endif
        if(irtc .ne. 0)then
          return
        endif
        if(ktfrealq(kx))then
          ispv=i-1
          maxi=int(rfromk(kx))
          cycle
        endif
        call tfgetoption('AccuracyGoal',ktastk(i),kx,irtc)
        if(irtc .ne. 0)then
          return
        endif
        if(ktfrealq(kx))then
          ispv=i-1
          eps=rfromk(kx)
          cycle
        endif
        call tfgetoption('Trace',ktastk(i),kx,irtc)
        if(irtc .ne. 0)then
          return
        endif
        if(ktfrealq(kx))then
          ispv=i-1
          trace=rfromk(kx) .ne. 0.d0
          cycle
        endif
        call tfgetoption('D',ktastk(i),kx,irtc)
        if(irtc .ne. 0)then
          return
        endif
        if(ktfrealq(kx))then
          ispv=i-1
          used=rfromk(kx) .ne. 0.d0
        endif
      enddo
 1    if(ispv .lt. isp1+2)then
        irtc=itfmessage(9,'General::narg','"2 or more"')
        return
      endif
      call tfsetupeqs(ktastk(isp1+1),ke,neq,irtc)
      if(irtc .ne. 0)then
        return
      endif
      call tfsetupvars(isp1+2,ispv,
     $     nvar,kav,kav0,v0,vmin,vmax,nvmax,irtc)
      if(irtc .ne. 0)then
        return
      endif
      ke=ktfcopy(ke)
      kdl(1:nvar)=ktfoper+mtfnull
      if(used)then
        call tfderiv(ke,nvar,kav,kdl,irtc)
        if(irtc .ne. 0)then
          call tflocal1(ke)
          go to 9000
        endif
      endif
      call tfnewton(ke,kav,v0,d0,kdl,
     $     vmin,vmax,neq,nvar,maxi,eps,trace,irtc)
      call tflocal1(ke)
      if(irtc .ne. 0)then
        go to 9000
      endif
      call tfassignrules(kav0,v0,nvar,kx)
      isp2=isp
      call tfgetllstkall(kx)
      call tfmakerulestk(ktfsymbol+itfres,int8(0))
      kax=ktfmakelist(isp2)
 9000 do i=1,nvar
        call tflocal(kdl(i))
        call tfdelete(kav(i)+4,.true.,.false.)
      enddo
      call tclrfpe
      return
      include 'inc/TFSF.inc'
      end

      subroutine tfassignrules(kav0,v0,nvar,kx)
      use tfstk
      implicit none
      integer*4 nvar,i
      integer*8 kav0(nvar),kx,kax,kai,ktadaloc,kfromr
      real*8 v0(nvar)
      kax=ktadaloc(-1,nvar)
      do i=1,nvar
        kai=ktadaloc(0,2)
        klist(kai)=ktfoper+mtfrule
        klist(kai+1)=ktfsymbol+ktfcopy1(kav0(i))
        klist(kai+2)=kfromr(v0(i))
        klist(kax+i)=ktflist+kai
      enddo
      kx=ktflist+kax
      return
      include 'inc/TFSF.inc'
      end

      subroutine tfnewton(ke,kav,v0,d0,kdl,vmin,vmax,
     $     neq,nvar,maxi,eps,trace,irtc)
      use tfstk
      implicit none
      integer*8 ke,kav(nvar),kdl(nvar),kax,kx
      integer*4 nvar,neq,maxi,irtc,i,j,iter
      real*8 v0(nvar),a(neq,nvar),f(neq),f0(neq),eps,dv(nvar),
     $     v(nvar),a0(neq,nvar),fact,fact1,fact2,d0,d,d1,d2,
     $     dg,am,sv,goal,s,tffsfmin,svi,vmin(nvar),vmax(nvar)
      logical*4 trace
      real*8 frac,factmin,svmin
      parameter (frac=1.d-7,factmin=1.d-4,svmin=1.d-7)
      iter=0
      call tfevalresidual(kav,v0,ke,f0,am,d0,nvar,neq,trace,irtc)
      goal=am*eps
      sv=0.d0
      do i=1,nvar
        sv=sv+abs(v0(i)*frac)+svmin
      enddo
      sv=sv/nvar
      fact=1.d0
      d1=d0
      fact1=0.d0
      if(irtc .ne. 0)then
        return
      endif
 1    if(d0 .lt. goal)then
        return
      endif
      iter=iter+1
      if(iter .gt. maxi)then
        return
      endif
      v(1:nvar)=v0(1:nvar)
      do i=1,nvar
c        call tfdebugprint(kdl(i),'tfnewton',2)
        call tfeevalref(kdl(i),kx,irtc)
        if(irtc .ne. 0)then
          return
        endif
        if(ktflistq(kx))then
          kax=ktfaddr(kx)
          if(ilist(2,kax-1) .eq. neq*2 .and.
     $         ktfreallistq(kax))then
            do j=1,neq
              a(j,i)=rlist(kax+j*2-1)-rlist(kax+j*2)
            enddo
            cycle
          endif
        elseif(ktfrealq(kx))then
          a(1:neq,i)=0.d0
          cycle
        endif
        svi=max(abs(v(i))*frac,sv)
        v(i)=v0(i)+svi
        call tfevalresidual(kav,v,ke,f,am,d1,nvar,neq,.false.,irtc)
        if(irtc .ne. 0)then
          return
        endif
        a(1:neq,i)=(f(1:neq)-f0(1:neq))/svi
        if(d1 .lt. d0)then
          v0(1:nvar)=v(1:nvar)
          f0(1:nvar)=f(1:nvar)
          d0=d1
        else
          v(i)=v0(i)
          rlist(kav(i))=v0(i)
        endif
      enddo
      a0=a
      f(1:nvar)=-f0(1:nvar)
      call tsolva(a,f,dv,neq,nvar,neq,1.d-8)
      dg=0.d0
      do i=1,neq
        s=0.d0
        do j=1,nvar
          s=s+a0(i,j)*dv(j)
        enddo
        dg=dg+f0(i)*s
      enddo
      dg=dg*2.d0
 2    v=min(vmax,max(vmin,v0+dv*fact))
c 2    do i=1,nvar
c        write(*,*)'newton ',i,v0(i),dv(i)
c        v(i)=min(vmax(i),max(vmin(i),v0(i)+dv(i)*fact))
c      enddo
      call tfevalresidual(kav,v,ke,f0,am,d,nvar,neq,trace,irtc)
      if(irtc .ne. 0)then
        return
      endif
      if(d .lt. d0)then
        v0=v
        d0=d
        fact=min(1.d0,fact*4.d0)
        fact1=0.d0
        d1=d
        go to 1
      else
        iter=iter+1
        if(iter .gt. maxi)then
          return
        endif
        d1=d
        d2=d1
        fact2=fact1
        fact1=fact
        fact=tffsfmin(fact1,fact2,d1,d2,d0,dg)
        if(fact .lt. factmin)then
          return
        endif
        go to 2
      endif
      include 'inc/TFSF.inc'
      end

      subroutine tfevalresidual(kav,v,ke,f,am,d,nvar,neq,trace,irtc)
      use tfstk
      implicit none
      integer*4 nvar,i,neq,l,itfuplevel,itfdownlevel,irtc,itfmessage
      integer*8 kav(nvar),ke,kx,kax
      real*8 v(nvar),f(neq),am,a,b,d
      logical*4 tflistqk,tfcomplexlistqk,trace
      do i=1,nvar
        if(trace)then
          write(*,*)'FindRoot Vars: ',i,v(i)
        endif
        rlist(kav(i))=v(i)
      enddo
      l=itfuplevel()
c      call tfdebugprint(ke,'evalresidual',2)
      call tfleval(ke,kx,.true.,irtc)
      if(irtc .ne. 0)then
        go to 9000
      endif
c      call tfdebugprint(kx,'==> ',2)
      if(.not. tflistqk(kx))then
        irtc=itfmessage(9,'General::wrongval','"Result of eqs","List"')
        go to 9000
      endif
      kax=ktfaddr(kx)
      if(ilist(2,kax-1) .ne. neq*2)then
        irtc=itfmessage(9,'General::wrongleng','"results","equations"')
        go to 9000
      endif
      if(ktfnonreallistq(kax))then
        if(tfcomplexlistqk(kx))then
          do i=1,neq
            if(ktfrealq(klist(kax+i*2-1)) .and.
     $           ktfrealq(klist(kax+i*2)))then
              f(i)=rlist(kax+i*2-1)-rlist(kax+i*2)
            else
              f(i)=1.d200
            endif
          enddo
          am=1.d300
          d=1.d300
        else
          irtc=itfmessage(9,'General::wrongval',
     $         '"results","numbers"')
        endif
        go to 9000
      endif
      am=0.d0
      d=0.d0
      do i=1,neq
        a=rlist(kax+i*2-1)
        b=rlist(kax+i*2)
        f(i)=a-b
        d=d+f(i)**2
        am=am+a**2+b**2
        if(trace)then
          write(*,*)'FindRoot eqs: ',i,a,b
        endif
      enddo
      if(trace)then
        write(*,*)'FindRoot Residual: ',d
      endif
 9000 l=itfdownlevel()
      return
      include 'inc/TFSF.inc'
      end

      subroutine tfsetupvars(isp1,isp2,
     $     nvar,kav,kav0,v0,vmin,vmax,nvmax,irtc)
      use tfstk
      implicit none
      integer*4 isp1,nvar,nvmax,i,j,isp0,isp2,irtc,ig,itfmessage
      integer*8 kav(nvmax),kav0(nvmax),kv,ka,k1,k2,k3,ka3,ka1,kavj,
     $     ki,ktnaloc1
      real*8 v0(nvmax),vmin(nvmax),vmax(nvmax),rfromk
      logical*4 tfnonlistqk
      nvar=isp2-isp1+1
      if(nvar .gt. nvmax)then
        irtc=itfmessage(9,'General::toomany','"variables"')
        return
      endif
      isp0=isp
      j=0
      do i=isp1,isp2
        ki=ktastk(i)
        if(tfnonlistqk(ki))then
          go to 8900
        endif
        ka=ktfaddr(ki)
        if(ilist(2,ka-1) .eq. 3)then
          k3=klist(ka+3)
          call tfeevalref(k3,k3,irtc)
          if(irtc .ne. 0)then
            go to 9000
          endif
          if(tfnonlistqk(k3))then
            go to 8900
          endif
          ka3=ktfaddr(k3)
          if(ktfnonreallistq(ka3) .or. ilist(2,ka3-1) .ne. 2)then
            go to 8900
          endif
          vmin(j+1)=rlist(ka3+1)
          vmax(j+1)=rlist(ka3+2)
        else
          if(ilist(2,ka-1) .ne. 2)then
            go to 8900
          endif
          vmin(j+1)=-1.d300
          vmax(j+1)=1.d300
        endif
        k1=klist(ka+1)
        if(ktfnonsymbolq(k1))then
          go to 8900
        endif
        ka1=ktfaddr(k1)
        j=j+1
        ig=max(0,ilist(2,ka1-1))
        kav0(j)=ka1
        kavj=ktnaloc1(ig,klist(ka1)-5)+4
        kav(j)=kavj
        call tflocal(klist(kavj))
        klist(kavj)=0
        isp=isp+2
        ktastk(isp-1)=kavj+4
        ktastk(isp)=ig
        k2=klist(ka+2)
        call tfeevalref(k2,kv,irtc)
        if(irtc .ne. 0)then
          go to 9000
        endif
        if(ktfnonrealq(kv))then
          irtc=itfmessage(9,'General::wrongval',
     $         '"initial value","Real"')
          go to 9000
        endif
        v0(j)=rfromk(kv)
      enddo
c      if(isp .gt. isp0)then
c        isp4=isp
c        call tfredefsymbol(isp0+1,isp4,
c     $       ite1,iae1,ve1,ite,iae,ve,rep)
c      endif
      isp=isp0
      return
 8900 irtc=itfmessage(9,'General::wrongtype',
     $     '"{var, ini} or {var, ini, {min, max}}"')
 9000 do i=1,j
        call tfdelete(kav(i)+4,.true.,.false.)
      enddo
      return
      include 'inc/TFSF.inc'
      end

      subroutine tfsetupeqs(kl,ke,neq,irtc)
      use tfstk
      implicit none
      integer*8 kl,kae,ke,kal,kei,kaei,ktfloadlstk,ktfstk2l,
     $     ktadaloc,ktfcopy
      integer*4 irtc,neq,j,i,itfmessage
      if(ktfnonlistq(kl))then
        go to 9000
      endif
      kal=ktfaddr(kl)
      if(klist(kal) .eq. ktfoper+mtfequal)then
        neq=1
        kae=ktfloadlstk(kal)
        klist(kae)=ktfoper+mtflist
        kae=ktfstk2l(kae)
        irtc=0
      elseif(klist(kal) .eq. ktfoper+mtflist)then
        neq=ilist(2,kal-1)
        kae=ktadaloc(-1,neq*2)
        do i=1,neq
          kei=klist(kal+i)
          if(ktflistq(kei))then
            kaei=ktfaddr(kei)
            if(klist(kaei) .eq. ktfoper+mtfequal)then
              klist(kae+i*2-1)=ktfcopy(klist(kaei+1))
              klist(kae+i*2  )=ktfcopy(klist(kaei+2))
              cycle
            endif
          endif
          do j=i*2-1,neq*2
            klist(kae+j)=ktfoper+mtfnull
          enddo
          go to 9000
        enddo
        irtc=0
      else
        go to 9000
      endif
      ke=ktflist+kae
      return
 9000 irtc=itfmessage(9,'General::wrongtype',
     $     '"f == g or List of them"')
      return
      include 'inc/TFSF.inc'
      end

      subroutine tffit(isp1,kx,irtc)
      use tfstk
      implicit none
      integer*4 nvmax,maxi0,irtc
      real*8 eps0
      parameter (nvmax=1024,maxi0=40,eps0=1.d-9)
      integer*4 isp1,nvar,i,maxi,ispv,isp2,n,m,
     $     iu,ig,itfmessage
      integer*8 kx,kav(nvmax),kdl(nvmax),kav0(nvmax),kas,kaiv,
     $     ke,kdp,ka1,kcv,kr,ktfsymbolf,ktfmaloc,ktaloc,
     $     ktfcopy,ktavaloc,ktfa2l,ktfmakelist,ktnaloc1
      real*8 v0(nvmax),r,gammaq,tinvgr,chin,rfromk,
     $     vmin(nvmax),vmax(nvmax),cut,cutoff,v0s(nvmax)
      logical*4 used
      integer*8 itfchisq,itfsigma,itfgood,itfconf,itfcov
      data itfchisq,itfsigma,itfgood,itfconf,itfcov/0,0,0,0,0/
      if(isp1 .gt. isp-4)then
        irtc=itfmessage(9,'General::narg','"4 or more"')
        return
      endif
      if(ktfnonsymbolq(ktastk(isp1+3)))then
        irtc=itfmessage(9,'General::wrongtype',
     $       '"Symbol for indep. var as #3"')
        return
      endif
      kas=ktfaddr(ktastk(isp1+3))
      if(itfchisq .eq. 0)then
        itfchisq=ktfsymbolf('ChiSquare',9,.true.)
        itfsigma=ktfsymbolf('StandardDeviation',17,.true.)
        itfgood=ktfsymbolf('GoodnessOfFit',13,.true.)
        itfconf=ktfsymbolf('ConfidenceInterval',18,.true.)
        itfcov=ktfsymbolf('CovarianceMatrix',16,.true.)
      endif
      maxi=maxi0
      used=.true.
      cutoff=0.d0
      ispv=isp
      do i=isp,isp1+4,-1
        call tfgetoption('MaxIterations',ktastk(i),kx,irtc)
        if(irtc .eq. -1)then
          ispv=i
          go to 1
        endif
        if(irtc .ne. 0)then
          return
        endif
        if(ktfrealq(kx))then
          ispv=i-1
          maxi=int(rfromk(kx))
          cycle
        endif
        call tfgetoption('D',ktastk(i),kx,irtc)
        if(irtc .ne. 0)then
          return
        endif
        if(ktfrealq(kx))then
          ispv=i-1
          used=rfromk(kx) .ne. 0.d0
          cycle
        endif
        call tfgetoption('Cutoff',ktastk(i),kx,irtc)
        if(irtc .ne. 0)then
          return
        endif
        if(ktfrealq(kx))then
          ispv=i-1
          cutoff=abs(rfromk(kx))
          cycle
        endif
      enddo
      irtc=itfmessage(9,'General::wrongopt',' ')
      return
 1    ig=max(0,ilist(2,kas-1))
      kaiv=ktnaloc1(ig,klist(kas)-5)+4
      call tfsetupvars(isp1+4,ispv,
     $     nvar,kav,kav0,v0,vmin,vmax,nvmax,irtc)
      ke=ktastk(isp1+2)
      if(irtc .ne. 0)then
        go to 9000
      endif
      kdp=ktfmaloc(ktastk(isp1+1),m,n,.false.,.true.,irtc)
      if(irtc .ne. 0)then
        if(irtc .eq. -1)then
          irtc=itfmessage(9,'General::wrongval',
     $         '"data","List of {x, y}, {x, y, dy}, or {x, y, dx, dy}"')
        endif
        go to 9100
      endif
      if(n .le. 1 .or. n .gt. 4)then
        irtc=itfmessage(9,'General::wrongval',
     $       '"data","List of {x, y}, {x, y, dy}, or {x, y, dx, dy}"')
        go to 9200
      endif
      kcv=ktaloc(nvar*m)
      iu=0
      ke=ktfcopy(ke)
      do i=1,nvar
        kdl(i)=ktfoper+mtfnull
      enddo
      if(used)then
        call tfderiv(ke,nvar,kav,kdl,irtc)
      endif
      call tflocal(klist(kaiv))
      klist(kaiv)=0
      if(irtc .eq. 0)then
        cut=cutoff
        if(n .le. 2)then
          if(cutoff .ne. 0.d0)then
            v0s=v0
c            call tmov(v0,v0s,nvar)
            call tffit1(rlist(kdp),n,m,ke,kaiv,nvar,kav,v0,
     $           kdl,vmin,vmax,r,rlist(kcv),eps0,maxi,0.d0,irtc)
            cut=sqrt(r/max(1,m-nvar))*cutoff
          endif
        endif
        call tffit1(rlist(kdp),n,m,ke,kaiv,nvar,kav,v0,
     $       kdl,vmin,vmax,r,rlist(kcv),eps0,maxi,cut,irtc)
      endif
      do i=1,nvar
        call tflocal(kdl(i))
      enddo
      call tflocal(ke)
      if(irtc .ne. 0)then
        call tfree(kcv)
        go to 9200
      endif
      call tfassignrules(kav0,v0,nvar,kr)
      isp2=isp
      call tfgetllstkall(kr)
      call tfmakerulestk(ktfsymbol+itfchisq,r)
      if(n .eq. 2)then
        call tfmakerulestk(ktfsymbol+itfgood,
     $       gammaq(dble(m-nvar)*.5d0,dble(m-nvar)*.5d0))
      else
        call tfmakerulestk(ktfsymbol+itfgood,
     $       gammaq(dble(m-nvar)*.5d0,r*.5d0))
      endif
      call tfcovmat(rlist(kcv),nvar,m)
      ka1=ktavaloc(-1,nvar)
      chin=tinvgr(dble(nvar))
      do i=0,nvar-1
        rlist(ka1+i+1)=sqrt(rlist(kcv+i*(m+1))*chin)
      enddo
      call tfmakerulestk(ktfsymbol+itfconf,ktflist+ka1)
      call tfmakerulestk(ktfsymbol+itfcov,ktflist+
     $     ktfa2l(kcv,m,nvar,m,.true.))
      kx=ktflist+ktfmakelist(isp2)
      isp=isp2-1
 9200 call tfree(kdp)
 9100 do i=1,nvar
        call tfdelete(kav(i)+4,.true.,.false.)
      enddo
 9000 call tfdelete(kaiv+4,.true.,.false.)
      call tclrfpe
      return
      include 'inc/TFSF.inc'
      end

      subroutine tfmakerulestk(kas,kx)
      use tfstk
      implicit none
      integer*8 kas,kx,ka1,ktadaloc,ktfcopy
      isp=isp+1
      ka1=ktadaloc(-1,2)
      klist(ka1)=ktfoper+mtfrule
      klist(ka1+1)=ktfcopy1(kas)
      klist(ka1+2)=ktfcopy(kx)
      ktastk(isp)=ktflist+ka1
      return
      include 'inc/TFSF.inc'
      end

      subroutine tfcovmat(c,n,ndim)
      implicit none
      integer*4 n,ndim,i,j,k
      real*8 c(ndim,n),ca(n,n),s
      do i=1,n
        do j=1,i
          s=0.d0
          do k=1,n
            s=s+c(k,i)*c(k,j)
          enddo
          ca(i,j)=s
        enddo
      enddo
      do i=1,n
        do j=1,i
          c(i,j)=ca(i,j)
          c(j,i)=c(i,j)
        enddo
      enddo
      return
      end

      subroutine tffit1(data,n,m,ke,kaiv,nvar,kav,v0,
     $     kdl,vmin,vmax,d0,a0,eps,maxi,cut,irtc)
      use tfstk
      implicit none
      integer*4 n,m,nvar,irtc,maxi,iter,i,j
      integer*8 ke,kaiv,kav(nvar),kdl(nvar),kaxvec,ktadaloc,
     $     ktavaloc
      real*8 data(n,m),v0(nvar),a0(m,nvar)
      real*8 a(m,nvar),abest(m,nvar),df(m),
     $     df0(m),d00,v00(nvar),w(nvar),
     $     vbest(nvar),dbest,dv(nvar),d2,eps,df2(m),svi,wi,db,
     $     fact,d0,d1,fact1,v(nvar),dg,s,d,fact2,tffsfmin,
     $     good,gammaq,ajump,vmin(nvar),vmax(nvar),sigma
      real*8 frac,factmin,svmin,factmin1,svdeps,cut
      parameter (frac=1.d-7,factmin=1.d-4,factmin1=1.d-4,
     $     svmin=1.d-7,svdeps=1.d-4)
      logical*4 newton,acalc
      v00=v0
      iter=0
      ajump=1.d0
      dbest=1.d300
      acalc=.false.
      kaxvec=ktadaloc(0,1)
      klist(kaxvec+1)=ktflist+ktavaloc(0,m)
      klist(kaxvec)=ktfcopy1(kxvect)
 21   call tfevalfit(df0,d0,data,n,m,ke,kaiv,nvar,kav,v0,
     $     kaxvec,.false.,cut,irtc)
      if(irtc .ne. 0)then
        call tflocal1(kaxvec)
        return
      endif
      if(d0 .lt. dbest)then
        dbest=d0
        vbest=v0
        if(acalc)then
          abest=a0
        endif
      endif
      newton=.true.
 11   fact=1.d0
      d1=d0
      fact1=0.d0
      d00=d0
 1    iter=iter+1
      if(iter .le. maxi)then
        v=v0
        do i=1,nvar
          call tfevalfit(df2,db,data,n,m,kdl(i),kaiv,nvar,kav,v,
     $         kaxvec,.true.,0.d0,irtc)
          if(irtc .eq. 0)then
            a(:,i)=df2
          elseif(irtc .eq. -1)then
            svi=min(max(svmin,abs(v0(i))*frac),
     $           (vmax(i)-vmin(i))*factmin)
            if(vmax(i)-v0(i) .ge. v0(i)-vmin(i))then
              svi=min(svi,(vmax(i)-v0(i))*.5d0)
            else
              svi=max(-svi,(vmin(i)-v0(i))*.5d0)
            endif
            v(i)=v0(i)+svi
            call tfevalfit(df2,db,data,n,m,ke,kaiv,nvar,kav,v,
     $           kaxvec,.false.,0.d0,irtc)
            if(irtc .ne. 0)then
              call tflocal1(kaxvec)
              return
            endif
            a(:,i)=(df2-df0)/svi
            if(db .lt. d0)then
              d0=db
              v0(i)=v0(i)+svi
              call tmov(df2,df0,m)
            else
              v(i)=v0(i)
            endif
          else
            call tflocal1(kaxvec)
            return
          endif
        enddo
        a0=a
        if(.not. acalc)then
          abest=a0
          acalc=.true.
        endif
        if(newton)then
          df=df0
          call tsolva(a,df,dv,m,nvar,m,1.d-6)
        else
          do i=1,nvar
            s=0.d0
            do j=1,m
              s=s+df0(j)*a0(j,i)
            enddo
            dv(i)=-s
          enddo
        endif
        dg=0.d0
        do i=1,m
          s=0.d0
          do j=1,nvar
            s=s+a0(i,j)*dv(j)
          enddo
          dg=dg+df0(i)*s
        enddo
        dg=dg*2.d0
        if(abs(dg) .lt. d0*eps)then
          good=gammaq(dble(m-nvar)*.5d0,d0*.5d0)
          if(good .gt. 0.001d0 .or. n .eq. 2)then
            iter=maxi
            go to 1
          else
            v0=max(vmin,min(vmax,v00+ajump*(v00-v0)))
            ajump=-ajump*2.d0
            go to 21
          endif          
        endif
 2      v=min(vmax,max(vmin,v0+dv*fact))
        call tfevalfit(df0,d,data,n,m,ke,kaiv,nvar,kav,v,
     $       kaxvec,.false.,cut,irtc)
        if(irtc .ne. 0)then
          call tflocal1(kaxvec)
          return
        endif
        if(d .lt. d0)then
          v0=v
          if(d .lt. dbest)then
            call tmov(v0,vbest,nvar)
            call tmov(a0,abest,m*nvar)
            dbest=d
          endif
          d0=d
          fact=min(1.d0,fact*4.d0)
          fact1=0.d0
          d1=d
          if(.not. newton)then
            newton=d .lt. 0.99d0*d00
            if(newton)then
              go to 11
            endif
          endif
          go to 1
        endif
        iter=iter+1
        if(iter .le. maxi)then
          d2=d1
          d1=d
          fact2=fact1
          fact1=fact
          fact=tffsfmin(fact1,fact2,d1,d2,d0,dg)
          if(newton)then
            if(fact .ge. factmin)then
              go to 2
            endif
            newton=.false.
            go to 11
          else
            if(fact .ge. factmin1)then
              go to 2
            else
              v0=max(vmin,min(vmax,2.d0*v00-v0))
              go to 21
            endif
          endif
        endif
      endif
      v0=vbest
      a0=abest
      d0=dbest
      if(cut .ne. 0.d0)then
        do j=1,m
          if(abs(df(j)) .gt. cut)then
            a0(j,:)=0.d0
          endif
        enddo
      endif
      call tsvdm(a0,0.d0,w,m,nvar,m,0,svdeps,.true.)
      if(n .eq. 2)then
        sigma=sqrt(d0/max(1,m-nvar))
        do i=1,nvar
          wi=w(i)*sigma
          a0(i,:)=a0(i,:)*wi
        enddo
      else
        do i=1,nvar
          wi=w(i)
          a0(i,:)=a0(i,:)*wi
        enddo
      endif
      call tflocal1(kaxvec)
      return
      include 'inc/TFSF.inc'
      end

      subroutine tfevalfit(df,r,a,n,m,ke,kaiv,nvar,kav,v,
     $     kaxvec,deriv,cutoff,irtc)
      use tfstk
      implicit none
      integer*4 n,m,nvar,irtc,i,itfuplevel,itfdownlevel,l,
     $     itfmessage
      integer*8 ke,kaiv,kav(nvar),kx,kavvec,k1,kaxvec,ka1,kax
      real*8 a(n,m),v(nvar),df(m),cutoff,r,vx,rfromk
      logical*4 tfcomplexqk,tfreallistq,deriv
      l=itfuplevel()
      do i=1,nvar
        rlist(kav(i))=v(i)
      enddo
      klist(kaiv)=ktflist+ktfcopy1(kaxvec)
      kavvec=ktfaddr(klist(kaxvec+1))
      do i=1,m
        rlist(kavvec+i)=a(1,i)
      enddo
      call tfeevalref(ke,kx,irtc)
      call tflocal1(kaxvec)
      klist(kaiv)=0
      if(irtc .ne. 0)then
        return
      endif
      if(ktflistq(kx))then
        kax=ktfaddr(kx)
        if(klist(kax) .eq. kxvect
     $       .and. ilist(2,kax-1) .eq. 1)then
          k1=klist(kax+1)
          ka1=ktfaddr(k1)
          if(tfreallistq(k1) .and. ilist(2,ka1-1) .eq. m)then
            if(deriv)then
              if(ktfnonreallistq(ka1))then
                do i=1,m
                  if(ktfrealq(klist(ka1+i)))then
                    df(i)=rlist(ka1+i)
                  else
                    df(i)=0.d0
                  endif
                enddo
              else
                df=rlist(ka1+1:ka1+m)
              endif
            else
              if(ktfnonreallistq(ka1))then
                do i=1,m
                  if(ktfrealq(klist(ka1+i)))then
                    df(i)=rlist(ka1+i)-a(2,i)
                  else
                    df(i)=1.d300
                  endif
                enddo
              else
                df=rlist(ka1+1:ka1+m)-a(2,:)
              endif
            endif
            go to 1000
          endif
        endif
      elseif(ktfrealq(kx))then
        if(deriv)then
          vx=rfromk(kx)
          df=vx
        else
          df=vx-a(2,:)
        endif
        go to 1000
      endif
      do i=1,m
        rlist(kaiv)=a(1,i)
        call tfeevalref(ke,kx,irtc)
        if(irtc .ne. 0)then
          go to 9000
        endif
c        call tfdebugprint(kx,'evalfit-i',2)
        if(ktfrealq(kx))then
          if(tfcomplexqk(kx))then
            if(deriv)then
              df(i)=0.d0
            else
              df(i)=1.d300
            endif
          else
            if(deriv)then
              irtc=-1
              go to 9000
            else
              irtc=itfmessage(9,'General::wrongval',
     $             '"result of equation","Real number"')
              go to 9000
            endif
          endif
        else
          vx=rfromk(kx)
          if(deriv)then
            df(i)=vx
          else
            df(i)=vx-a(2,i)
          endif
        endif
      enddo
 1000 r=0.d0
      if(deriv)then
        if(n .eq. 3)then
          do i=1,m
            df(i)=df(i)/a(3,i)
          enddo
        elseif(n .eq. 4)then
          do i=1,m
            df(i)=df(i)/a(4,i)
          enddo
        endif
      else
        if(cutoff .ne. 0.d0)then
          if(n .eq. 3)then
            do i=1,m
              df(i)=df(i)/a(3,i)
              r=r+min(cutoff,max(-cutoff,df(i)))**2
            enddo
          elseif(n .eq. 4)then
            do i=1,m
              df(i)=df(i)/a(4,i)
              r=r+min(cutoff,max(-cutoff,df(i)))**2
            enddo
          else
            do i=1,m
              r=r+min(cutoff,max(-cutoff,df(i)))**2
            enddo
          endif
        else
          if(n .eq. 3)then
            do i=1,m
              df(i)=df(i)/a(3,i)
              r=r+df(i)**2
            enddo
          elseif(n .eq. 4)then
            do i=1,m
              df(i)=df(i)/a(4,i)
              r=r+df(i)**2
            enddo
          else
            do i=1,m
              r=r+df(i)**2
            enddo
          endif
        endif
      endif
 9000 l=itfdownlevel()
      return
      include 'inc/TFSF.inc'
      end

      subroutine tfderiv(ke,nvar,kav,kdl,irtc)
      use tfstk
      implicit none
      integer*4 nvar,irtc,i,isp0
      integer*8 ke,kdl(nvar),kav(nvar),kr,kd,
     $     ks,ktfsymbolf,ktfcopy,ierr0
      logical*4 euv
      integer*8 iads,iader
      data iads,iader /0,0/
      if(iads .eq. 0)then
        iads=ktfsymbolf('System`D',8,.false.)
        call tfsyeval(iads,ks,irtc)
        if(irtc .ne. 0)then
          return
        endif
        iader=ktfsymbolf('CheckDerivative',15,.false.)
      endif
      do i=1,nvar
        klist(kav(i))=ktfsymbol+kav(i)+4
      enddo
      isp0=isp
      isp=isp0+1
      ktastk(isp)=ktfsymbol+iads
      do i=1,nvar
        isp=isp0+2
        ktastk(isp)=ke
        isp=isp+1
        ktastk(isp)=ktfsymbol+kav(i)+4
        ierr0=ierrorexp
        ierrorexp=1
c        call tfdebugprint(ke,'tfderiv D',2)
c        call tfdebugprint(ktastk(isp),'by ',2)
        call tfdeval(isp0+1,iads,kd,1,int8(0),.false.,euv,irtc)
c        call tfdebugprint(kd,'==>',2)
        ierrorexp=ierr0
        if(irtc .ne. 0)then
          go to 100
        endif
        isp=isp0+2
        ktastk(isp)=ktfsymbol+iader
        isp=isp+1
        ktastk(isp)=ktfcopy(kd)
        ierr0=ierrorexp
        ierrorexp=1
        call tfdeval(isp0+2,iader,kr,1,int8(0),.false.,euv,irtc)
        ierrorexp=ierr0
        if(irtc .eq. 0 .and. ktfrealq(kr))then
          kdl(i)=kd
        else
          call tflocal(kd)
          kdl(i)=ktfoper+mtfnull
          if(irtc .gt. 0. and. ierrorprint .ne. 0)then
            call tfreseterror
          endif
        endif
      enddo
 100  do i=1,nvar
        klist(kav(i))=0
      enddo
      irtc=0
      isp=isp0
      return
      include 'inc/TFSF.inc'
      end

      real*8 function tinvgr(an)
      implicit none
      real*8 an,c0,x0,gn,erfc,gamma,gammaq,df,dfdx,anh
      if(an .eq. 1.d0)then
        tinvgr=1.d0
        return
      endif
      c0=erfc(sqrt(0.5d0))
      anh=an*.5d0
      x0=anh
      gn=gamma(anh)
      df=gammaq(anh,x0)-c0
      do while(abs(df) .gt. 1.d-14)
        dfdx=-exp(-x0)*x0**(anh-1.d0)/gn
        x0=x0-df/dfdx
c        write(*,*)'tinvgr ',x0,dfdx,df
        df=gammaq(anh,x0)-c0
      enddo
      tinvgr=x0*2.d0
      return
      end
