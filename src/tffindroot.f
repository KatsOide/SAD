      module findr
      use tfstk
      real*8 ,parameter :: alarge=1.d300
      integer*8 ,save :: itfres=0
      type (sad_descriptor), save ::
     $     itfchisq,itfsigma,itfgood,itfconf,itfcov,itfdm
      data itfchisq%k,itfsigma%k,itfgood%k,itfconf%k,itfcov%k,
     $     itfdm%k
     $     /0,0,0,0,0,0/
      type(sad_descriptor) ,save::iads,iader
      data iads%k /0/

      type symv
      sequence
      type (sad_symdef), pointer :: p
      end type

      contains
      subroutine tfassignrules(sav0,v0,nvar,klx)
      implicit none
      type (sad_descriptor) kx
      type (sad_dlist), pointer ,intent(out) ::klx
      type (sad_dlist), pointer ::kli
      integer*4 ,intent(in)::nvar
      integer*4 i
      type (symv),intent(in):: sav0(nvar)
      real*8 ,intent(in):: v0(nvar)
      kx=kxadaloc(-1,nvar,klx)
      do i=1,nvar
        klx%dbody(i)=kxadaloc(0,2,kli)
        kli%head%k=ktfoper+mtfrule
        kli%dbody(1)=dtfcopy1(sad_descr(sav0(i)%p%sym))
        kli%rbody(2)=v0(i)
      enddo
      return
      end subroutine

      real*8 pure function fmincube(f1,f2,g1,g2,g0,dg) result(f)
      use macmath,only:m_euler
      implicit none
      real*8 ,intent(in):: f1,f2,g1,g2,g0,dg
      real*8 a,b,s
      if(f2 == 0.d0)then
        f=-.5d0*f1*dg/((g1-g0)/f1-dg)
      else
        a=((g1-g0)/f1**2-(g2-g0)/f2**2)/(f1-f2)+dg/f1/f2
        b=(-f2*(g1-g0)/f1**2+f1*(g2-g0)/f2**2)-dg*(f1+f2)/f1/f2
        s=sqrt(max(0.d0,b**2-3.d0*a*dg))
        if(b > 0.d0)then
          f=-dg/(s+b)
        else
          f=(s-b)/3.d0/a
        endif
c        f=merge(-dg/(s+b),(s-b)/3.d0/a,b > 0.d0)
      endif
      f=min(m_euler*f1,max(f1/16.d0,f))
      return
      end

      subroutine tffindroot(isp1,kx,irtc)
      use modul,only:tfdelete
      implicit none
      real*8 , parameter :: eps0=1.d-20, frac0=1.d-7
      integer*4 , parameter :: nvmax=2048, maxi0=50
      type (sad_descriptor) ,intent(out)::kx
      type (sad_descriptor) ke,kdl(nvmax)
      type (symv) , allocatable::sav(:),sav0(:)
      type (sad_dlist), pointer :: klx
      type (sad_rlist), pointer :: klo
      integer*4 ,intent(in):: isp1
      integer*4 ,intent(out):: irtc
      integer*4 neq,nvar,itfmessage,isp2,i,maxi,ispv
      real*8 eps,d0,frac
      real*8 , allocatable :: v0(:),vmin(:),vmax(:)
      logical*4 trace,used
      if(isp < isp1+2)then
        irtc=itfmessage(9,'General::narg','"2 or more"')
        return
      endif
      if(itfres == 0)then
        itfres=ktfsymbolz('Residual',8)
      endif
      maxi=maxi0
      eps=eps0
      trace=.false.
      used=.true.
      ispv=isp
      if(isp >= isp1+3)then
        if(tfreallistq(dtastk(isp),klo)
     $       .and. klo%nl >= 4)then
          ispv=isp-1
          maxi=int(klo%rbody(1))
          eps=klo%rbody(2)
          trace=klo%rbody(3) /= 0.d0
          used=klo%rbody(4) >= 1.d0
        endif          
      endif
      if(.not. used)then
        frac=merge(klo%rbody(4),frac0,klo%rbody(4) > 0.d0)
      endif
      if(ispv < isp1+2)then
        irtc=itfmessage(9,'General::narg','"2 or more"')
        return
      endif
      call tfsetupeqs(ktastk(isp1+1),ke%k,neq,irtc)
      if(irtc /= 0)then
        return
      endif
      allocate (sav(nvmax),sav0(nvmax),
     $     v0(nvmax),vmin(nvmax),vmax(nvmax))
      call tfsetupvars(isp1+2,ispv,
     $     nvar,sav,sav0,v0,vmin,vmax,nvmax,irtc)
      if(irtc /= 0)then
        deallocate (sav,sav0,v0,vmin,vmax)
        return
      endif
      ke=dtfcopy(ke)
      kdl(1:nvar)%k=ktfoper+mtfnull
c     write(*,*)'findroot-D ',used
      if(used)then
        call tfderiv(ke,nvar,sav,kdl,irtc)
        if(irtc /= 0)then
          call tflocal1(ke%k)
          go to 9000
        endif
c     call tfdebugprint(ke,'tffindroot-deriv',1)
      endif
      call tfnewton(ke%k,sav,v0,d0,kdl,
     $     vmin,vmax,neq,nvar,maxi,eps,trace,frac,irtc)
      call tflocal1(ke%k)
      if(irtc /= 0)then
        go to 9000
      endif
      call tfassignrules(sav0,v0,nvar,klx)
      isp2=isp
      call tfgetllstkall(klx)
      call tfmakerulestk(sad_descr(ktfsymbol+itfres),d0)
      kx=kxmakelist(isp2)
 9000 do i=1,nvar
        call tflocal(kdl(i))
        call tfdelete(sav(i)%p,.true.,.false.)
      enddo
      deallocate (sav,sav0,v0,vmin,vmax)
      call tclrfpe
      return
      end

      subroutine tfnewton(ke,sav,v0,d0,kdl,vmin,vmax,
     $     neq,nvar,maxi,eps,trace,frac,irtc)
      use eeval
      implicit none
      type (sad_dlist), pointer :: klx
      type (symv) sav(nvar)
      type (sad_descriptor) kdl(nvar),kx
      integer*8 ,intent(in):: ke
      integer*4 ,intent(in)::  nvar,neq,maxi
      integer*4 ,intent(out):: irtc
      integer*4 i,j,iter
      real*8 ,intent(inout):: v0(nvar),d0
      real*8 ,intent(in):: vmin(nvar),vmax(nvar),eps,frac
      real*8 ,allocatable,dimension(:):: f,f0,dv,v
      real*8 fact,fact1,fact2,d,d1,d2,dg,am,sv,goal,svi
      real*8 , allocatable :: a(:,:),a0(:,:)
      logical*4 ,intent(in):: trace
      real*8 , parameter :: factmin=1.d-4,svmin=1.d-7,svdtol=1.d-5
      allocate (f(neq),f0(neq),dv(nvar),v(nvar))
      iter=0
      call tfevalresidual(sav,v0,ke,f0,am,d0,nvar,neq,trace,irtc)
      goal=am*eps
      sv=svmin+sum(abs(v0))*frac
c     sv=0.d0
c     do i=1,nvar
c     sv=sv+abs(v0(i)*frac)+svmin
c     enddo
      sv=sv/nvar
      fact=1.d0
      d1=d0
      fact1=0.d0
      if(irtc /= 0)then
        return
      endif
      allocate (a(neq,nvar),a0(neq,nvar))
 1    if(d0 < goal)then
        deallocate (a0,a)
        return
      endif
      iter=iter+1
      if(iter > maxi)then
        deallocate (a0,a)
        return
      endif
      v=v0
      do i=1,nvar
        kx=tfeevalref(kdl(i),irtc)
        if(irtc /= 0)then
          deallocate (a0,a)
          return
        endif
        if(ktflistq(kx,klx))then
          if(klx%nl == neq*2 .and. ktfreallistq(klx))then
            do j=1,neq
              a(j,i)=klx%rbody(j*2-1)-klx%rbody(j*2)
            enddo
            cycle
          endif
        elseif(ktfrealq(kx))then
          a(1:neq,i)=0.d0
          cycle
        endif
        svi=max(abs(v(i))*frac,sv)
        v(i)=v0(i)+svi
        call tfevalresidual(sav,v,ke,f,am,d1,nvar,neq,.false.,irtc)
        if(irtc /= 0)then
          deallocate (a0,a)
          return
        endif
        a(1:neq,i)=(f(1:neq)-f0(1:neq))/svi
        if(d1 < d0)then
          v0=v
          f0=f
          d0=d1
        else
          v(i)=v0(i)
          sav(i)%p%value=dfromr(v0(i))
        endif
      enddo
      a0=a
      f=-f0
      call tsolva(a,f,dv,neq,nvar,neq,svdtol)
      dg=0.d0
      do i=1,neq
c     s=0.d0
c     do j=1,nvar
c     s=s+a0(i,j)*dv(j)
c     enddo
        dg=dg+f0(i)*sum(a0(i,1:nvar)*dv(1:nvar))
      enddo
      dg=dg*2.d0
 2    v=min(vmax,max(vmin,v0+dv*fact))
c     2    do i=1,nvar
c     write(*,*)'newton ',i,v0(i),dv(i)
c     v(i)=min(vmax(i),max(vmin(i),v0(i)+dv(i)*fact))
c     enddo
      call tfevalresidual(sav,v,ke,f0,am,d,nvar,neq,trace,irtc)
      if(irtc /= 0)then
        deallocate (a0,a)
        return
      endif
      if(d < d0)then
        v0=v
        d0=d
        fact=min(1.d0,fact*4.d0)
        fact1=0.d0
        d1=d
        go to 1
      else
        iter=iter+1
        if(iter > maxi)then
          deallocate (a0,a)
          return
        endif
        d2=d1
        d1=d
        fact2=fact1
        fact1=fact
        fact=fmincube(fact1,fact2,d1,d2,d0,dg)
        if(fact < factmin)then
          deallocate (a0,a)
          return
        endif
        go to 2
      endif
      end

      subroutine tfevalresidual(sav,v,ke,f,am,d,nvar,neq,trace,irtc)
      use sameq,only:tfcomplexlistqk
      use eeval
      implicit none
      integer*4 ,intent(in):: nvar,neq
      integer*4 ,intent(out):: irtc
      integer*4 i,l,itfmessage
      type (symv) ,intent(in):: sav(nvar)
      type (sad_descriptor) kx
      type (sad_dlist),pointer::kl
      integer*8 ke,kax
      real*8 ,intent(out):: f(neq),am,d
      real*8 ,intent(in):: v(nvar)
      real*8 a,b
      logical*4 ,intent(in):: trace
      do i=1,nvar
        if(trace)then
          write(*,*)'FindRoot Vars: ',i,v(i)
        endif
        sav(i)%p%value=dfromr(v(i))
      enddo
      l=itfuplevel()
      call loc_sad(ktfaddr(ke),kl)
      kx=tfleval(kl,.true.,irtc)
c     call tfdebugprint(kx,'evalres-2',1)
      if(irtc /= 0)then
        go to 9000
      endif
      if(.not. tflistq(kx))then
        irtc=itfmessage(9,'General::wrongval','"Result of eqs","List"')
        go to 9000
      endif
      kax=ktfaddr(kx%k)
      if(ilist(2,kax-1) /= neq*2)then
        irtc=itfmessage(9,'General::wrongleng','"results","equations"')
        go to 9000
      endif
      if(ktfnonreallistq(kax))then
        if(tfcomplexlistqk(kx))then
          do i=1,neq
            if(ktfrealq(klist(kax+i*2-1)) .and. ktfrealq(klist(kax+i*2)))then
              f(i)=rlist(kax+i*2-1)-rlist(kax+i*2)
            else
              f=alarge
            endif
c            f(i)=merge(rlist(kax+i*2-1)-rlist(kax+i*2),1.d200,
c     $           ktfrealq(klist(kax+i*2-1)) .and.
c     $           ktfrealq(klist(kax+i*2)))
          enddo
          am=alarge
          d=alarge
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
      end

      subroutine tfsetupvars(isp1,isp2,
     $     nvar,sav,sav0,v0,vmin,vmax,nvmax,irtc)
      use eeval
      use modul,only:tfdelete
      implicit none
      type (sad_descriptor) k1,ki,kv,k3
      type (sad_symbol), pointer :: sym1
      type (sad_dlist), pointer :: kli,kl3
      integer*4 ,intent(in):: isp1,isp2,nvmax
      integer*4 ,intent(out):: nvar,irtc
      integer*4 i,j,isp0,ig,itfmessage
      type (symv) ,intent(inout):: sav(nvmax),sav0(nvmax)
      real*8 ,intent(out):: v0(nvmax),vmin(nvmax),vmax(nvmax)
      nvar=isp2-isp1+1
      if(nvar > nvmax)then
        irtc=itfmessage(9,'General::toomany','"variables"')
        return
      endif
      isp0=isp
      j=0
      do i=isp1,isp2
        ki=dtastk(i)
        if(.not. tflistq(ki,kli))then
          go to 8900
        endif
        if(kli%nl == 3)then
          k3=tfeevalref(kli%dbody(3),irtc)
          if(irtc /= 0)then
            go to 9000
          endif
          if(.not. tflistq(k3,kl3))then
            go to 8900
          endif
          if(ktfnonreallistqo(kl3) .or. kl3%nl /= 2)then
            go to 8900
          endif
          vmin(j+1)=kl3%rbody(1)
          vmax(j+1)=kl3%rbody(2)
        else
          if(kli%nl /= 2)then
            go to 8900
          endif
          vmin(j+1)=-alarge
          vmax(j+1)=alarge
        endif
        k1=kli%dbody(1)
        if(ktfnonsymbolq(k1,sym1))then
          go to 8900
        endif
        j=j+1
        call sym_symdef(sym1,sav0(j)%p)
        ig=max(0,sym1%gen)
        call descr_symdef(kxnaloc1(ig,sym1%loc),sav(j)%p)
        call tflocald(sav(j)%p%value)
        sav(j)%p%value%k=0
c     isp=isp+2
c     ktastk(isp-1)=sad_loc(symd%sym%loc)
c     ktastk(isp)=ig
        kv=tfeevalref(kli%dbody(2),irtc)
        if(irtc /= 0)then
          go to 9000
        endif
        if(ktfnonrealq(kv,v0(j)))then
          irtc=itfmessage(9,'General::wrongval',
     $         '"initial value","Real"')
          go to 9000
        endif
        v0(j)=min(vmax(j),max(vmin(j),v0(j)))
      enddo
c     if(isp > isp0)then
c     isp4=isp
c     call tfredefsymbol(isp0+1,isp4,
c     $       ite1,iae1,ve1,ite,iae,ve,rep)
c     endif
      isp=isp0
      return
 8900 irtc=itfmessage(9,'General::wrongtype',
     $     '"{var, ini} or {var, ini, {min, max}}"')
 9000 do i=1,j
        call tfdelete(sav(i)%p,.true.,.false.)
      enddo
      return
      end

      subroutine tfsetupeqs(kl,ke,neq,irtc)
      use part,only:tfreplist
      implicit none
      type (sad_dlist), pointer :: list,listi,liste
      type (sad_descriptor) kei
      integer*8 ,intent(in):: kl
      integer*8 ,intent(out):: ke
      integer*8 kae
      integer*4 ,intent(out):: neq,irtc
      integer*4 i,itfmessage
      logical*4 eval
      if(ktfnonlistq(kl,list))then
        go to 9000
      endif
      if(list%head%k == ktfoper+mtfequal)then
        neq=1
        list=>tfduplist(list)
        call tfreplist(list,0,dfromk(ktfoper+mtflist),eval)
        kae=ksad_loc(list%head%k)
        irtc=0
      elseif(list%head%k == ktfoper+mtflist)then
        neq=list%nl
        kae=ktadaloc(-1,neq*2,liste)
        do i=1,neq
          kei=list%dbody(i)
          if(ktflistq(kei,listi))then
            if(listi%head%k == ktfoper+mtfequal)then
              liste%dbody(i*2-1)=dtfcopy(listi%dbody(1))
              liste%dbody(i*2  )=dtfcopy(listi%dbody(2))
              cycle
            endif
          endif
          liste%dbody(i*2-1:neq*2)%k=ktfoper+mtfnull
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
      end

      subroutine tffit(isp1,kx,irtc)
      use repl, only:tfgetoption
      use modul,only:tfdelete
      use maloc,only:ktfmaloc
      use iso_c_binding
      use gammaf
      implicit none
      type (sad_descriptor) ,intent(out)::kx
      integer*4 ,intent(in):: isp1
      integer*4 ,intent(out):: irtc
      type (sad_descriptor) ke
      type (sad_symbol), pointer :: sym
      type (sad_symdef), pointer :: symdv
      type (sad_dlist), pointer :: klr
      integer*4 ,parameter ::nvmax=1024,maxi0=40
      real*8 ,parameter:: eps0=1.d-9
      integer*4 nvar,i,maxi,ispv,isp2,n,m,iu,ig,itfmessage
      type (symv) , allocatable::sav(:),sav0(:)
      integer*8 kdp
      type (sad_descriptor) kdl(nvmax),kdm,kcv,kci
      real*8 , allocatable::v0(:),vmin(:),vmax(:),v0s(:)
      real*8 r,vx,cut,cutoff
      real*8 , pointer :: datap(:,:)
      logical*4 used
      if(isp1 > isp-4)then
        irtc=itfmessage(9,'General::narg','"4 or more"')
        return
      endif
      if(ktfnonsymbolq(ktastk(isp1+3),sym))then
        irtc=itfmessage(9,'General::wrongtype',
     $       '"Symbol for indep. var as #3"')
        return
      endif
      if(itfchisq%k == 0)then
        itfchisq=kxsymbolf('ChiSquare',9,.true.)
        itfsigma=kxsymbolf('StandardDeviation',17,.true.)
        itfgood=kxsymbolf('GoodnessOfFit',13,.true.)
        itfconf=kxsymbolf('ConfidenceInterval',18,.true.)
        itfcov=kxsymbolf('CovarianceMatrix',16,.true.)
        itfdm=kxsymbolf('DesignMatrix',12,.true.)
      endif
      maxi=maxi0
      used=.true.
      cutoff=0.d0
      ispv=isp
      do i=isp,isp1+4,-1
        call tfgetoption('MaxIterations',dtastk(i),kx,irtc)
        if(irtc == -1)then
          ispv=i
          allocate (sav(nvmax),sav0(nvmax),
     $         v0(nvmax),vmin(nvmax),vmax(nvmax),v0s(nvmax))
          go to 1
        endif
        if(irtc /= 0)then
          return
        endif
        if(ktfrealq(kx,maxi))then
          ispv=i-1
          cycle
        endif
        call tfgetoption('D',dtastk(i),kx,irtc)
        if(irtc /= 0)then
          return
        endif
        if(ktfrealq(kx,vx))then
          ispv=i-1
          used=vx /= 0.d0
          cycle
        endif
        call tfgetoption('Cutoff',dtastk(i),kx,irtc)
        if(irtc /= 0)then
          return
        endif
        if(ktfrealq(kx))then
          ispv=i-1
          cutoff=abs(kx%x(1))
          cycle
        endif
      enddo
      irtc=itfmessage(9,'General::wrongopt',' ')
      return
 1    ig=max(0,sym%gen)
      call descr_symdef(kxnaloc1(ig,sym%loc),symdv)
      call tfsetupvars(isp1+4,ispv,
     $     nvar,sav,sav0,v0,vmin,vmax,nvmax,irtc)
      if(irtc /= 0)then
        go to 9000
      endif
      ke=dtastk(isp1+2)
      kdp=ktfmaloc(dtastk(isp1+1),m,n,.false.,.true.,irtc)
      if(irtc /= 0)then
        if(irtc == -1)then
          irtc=itfmessage(9,'General::wrongval',
     $         '"data","List of {x, y}, {x, y, dy}, or {x, y, dx, dy}"')
        endif
        go to 9100
      endif
      if(n <= 1 .or. n > 4)then
        irtc=itfmessage(9,'General::wrongval',
     $       '"data","List of {x, y}, {x, y, dy}, or {x, y, dx, dy}"')
        go to 9200
      endif
      call c_f_pointer(c_loc(rlist(kdp)),datap,[n,m])
      iu=0
      ke=dtfcopy(ke)
      kdl(1:nvar)%k=ktfoper+mtfnull
      if(used)then
        call tfderiv(ke,nvar,sav,kdl,irtc)
      endif
      call tflocald(symdv%value)
      symdv%value%k=0
      if(irtc == 0)then
        cut=cutoff
        if(n <= 2)then
          if(cutoff /= 0.d0)then
            v0(1:nvar)=max(vmin(1:nvar),min(vmax(1:nvar),v0(1:nvar)))
            v0s(1:nvar)=v0(1:nvar)
            call tffit1(datap,n,m,ke,symdv,nvar,sav,v0,
     $           kdl,vmin,vmax,r,kdm,kcv,kci,eps0,maxi,0.d0,irtc)
            cut=sqrt(r/max(1,m-nvar))*cutoff
          endif
        endif
        v0(1:nvar)=max(vmin(1:nvar),min(vmax(1:nvar),v0(1:nvar)))
        call tffit1(datap,n,m,ke,symdv,nvar,sav,v0,
     $       kdl,vmin,vmax,r,kdm,kcv,kci,eps0,maxi,cut,irtc)
      endif
      v0(1:nvar)=max(vmin(1:nvar),min(vmax(1:nvar),v0(1:nvar)))
      do i=1,nvar
        call tflocal(kdl(i))
      enddo
      call tflocald(ke)
      if(irtc /= 0)then
        go to 9200
      endif
      call tfassignrules(sav0,v0,nvar,klr)
      isp2=isp
      call tfgetllstkall(klr)
      call tfmakerulestk(itfchisq,r)
      if(n == 2)then
        call tfmakerulestk(itfgood,
     $       sad_descr(gammaq(dble(m-nvar)*.5d0,dble(m-nvar)*.5d0)))
      else
        call tfmakerulestk(itfgood,
     $       sad_descr(gammaq(dble(m-nvar)*.5d0,r*.5d0)))
      endif
      call tfmakerulestk(itfconf,kci)
      call tfmakerulestk(itfcov,kcv)
      call tfmakerulestk(itfdm,kdm)
      kx=kxmakelist(isp2)
c     call tfdebugprint(kx,'tffit-8',3)
      isp=isp2-1
 9200 call tfree(kdp)
 9100 do i=1,nvar
        call tfdelete(sav(i)%p,.true.,.false.)
      enddo
 9000 call tfdelete(symdv,.true.,.false.)
      deallocate (sav,sav0,v0,vmin,vmax,v0s)
      call tclrfpe
      return
      end

      subroutine tfcovmat(c,ca,m,n,ndim)
      implicit none
      integer*4 ,intent(in):: n,ndim,m
      integer*4 i,j
      real*8 ,intent(in):: c(ndim,n)
      real*8 ,intent(out):: ca(n,n)
      real*8 s
      real*8, save:: rnan=0.d0
c     write(*,*)'covmat ',n,m,ndim
      if(rnan == 0.d0)then
        rnan=rtfnan
      endif
      do i=1,n
        do j=1,i
          s=sum(c(1:m,i)*c(1:m,j))
          ca(i,j)=merge(rnan,s,ktfenanq(s))
          ca(j,i)=ca(i,j)
        enddo
      enddo
      return
      end

      subroutine tffit1(data,n,m,ke,symdv,nvar,sav,v0,
     $     kdl,vmin,vmax,d0,kdm,kcv,kci,eps,maxi,cut,irtc)
      use gammaf
      implicit none
      type (sad_symdef) ,intent(inout):: symdv
      integer*4 ,intent(in):: n,m,nvar,maxi
      integer*4 ,intent(out):: irtc
      integer*4 iter,i,j
      type (symv) ,intent(in):: sav(nvar)
      type (sad_descriptor) ,intent(out):: kdm,kcv,kci
      type (sad_descriptor) ,intent(in):: ke,kdl(nvar)
      type (sad_rlist) , pointer :: klci
      integer*8 kaxvec
      real*8 ,intent(in):: data(n,m)
      real*8 ,intent(inout):: v0(nvar)
      real*8 , allocatable :: a0(:,:),a(:,:),abest(:,:),
     $     df(:),df0(:),v00(:),w(:),cv(:,:),vbest(:),dv(:),df2(:)
      real*8 ,intent(in):: eps,cut
      real*8 ,intent(out):: d0
      real*8 d00,d2,svi,db,dbest,
     $     fact,d1,fact1,v(nvar),dg,d,fact2,
     $     good,ajump,vmin(nvar),vmax(nvar),sigma
      real*8 chin,x
      real*8 ,parameter ::frac=1.d-7,factmin=1.d-4,factmin1=1.d-4,
     $     svmin=1.d-7,svdeps=1.d-4
      logical*4 newton,acalc
      allocate (a0(m,nvar),a(m,nvar),abest(m,nvar),df(m),
     $     df0(m),v00(nvar),w(nvar),cv(nvar,nvar),
     $     vbest(nvar),dv(nvar),df2(m))
      kdm%k=0
      kcv%k=0
      v0=max(vmin,min(vmax,v0))
      v00=v0
      iter=0
      ajump=1.d0
      dbest=alarge
      acalc=.false.
      kaxvec=ktadaloc(0,1)
      klist(kaxvec+1)=ktflist+ktavaloc(0,m)
      klist(kaxvec)=ktfcopy1(kxvect)
 21   v0=max(vmin,min(vmax,v0))
      call tfevalfit(df0,d0,data,n,m,ke,symdv,nvar,sav,v0,
     $     kaxvec,.false.,cut,irtc)
      if(irtc /= 0)then
        deallocate(a0,a,abest,df,df0,v00,w,cv,vbest,dv,df2)
        call tflocal1(kaxvec)
        return
      endif
      if(d0 < dbest)then
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
      if(iter <= maxi)then
        v=v0
        do i=1,nvar
          call tfevalfit(df2,db,data,n,m,kdl(i),symdv,nvar,sav,v,
     $         kaxvec,.true.,0.d0,irtc)
          if(irtc == 0)then
            a(:,i)=df2
          elseif(irtc == -1)then
            svi=min(max(svmin,abs(v0(i))*frac),
     $           (vmax(i)-vmin(i))*factmin)
            if(vmax(i)-v0(i) >= v0(i)-vmin(i))then
              svi=min(svi,(vmax(i)-v0(i))*.5d0)
            else
              svi=max(-svi,(vmin(i)-v0(i))*.5d0)
            endif
            v(i)=v0(i)+svi
            call tfevalfit(df2,db,data,n,m,ke,symdv,nvar,sav,v,
     $           kaxvec,.false.,0.d0,irtc)
            if(irtc /= 0)then
              deallocate(a0,a,abest,df,df0,v00,w,cv,vbest,dv,df2)
              call tflocal1(kaxvec)
              return
            endif
            a(:,i)=(df2-df0)/svi
            if(db < d0)then
              d0=db
              v0(i)=v0(i)+svi
              df0=df2
            else
              v(i)=v0(i)
            endif
          else
            deallocate(a0,a,abest,df,df0,v00,w,cv,vbest,dv,df2)
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
          df=-df0
          call tsolva(a,df,dv,m,nvar,m,svdeps)
        else
          do i=1,nvar
            dv(i)=-sum(df0(1:m)*a0(1:m,i))
          enddo
        endif
        dg=0.d0
        do i=1,m
          dg=dg+df0(i)*dot_product(a0(i,1:nvar),dv(1:nvar))
        enddo
        dg=dg*2.d0
        if(abs(dg) < d0*eps)then
          good=gammaq(dble(m-nvar)*.5d0,d0*.5d0)
          if(good > 0.001d0 .or. n == 2)then
            iter=maxi
            go to 1
          else
            do i=1,nvar
              v0(i)=max(vmin(i),min(vmax(i),
     $             v00(i)+ajump*(v00(i)-v0(i))))
            enddo
            ajump=-ajump*2.d0
            go to 21
          endif          
        endif
 2      v=min(vmax,max(vmin,v0+dv*fact))
        call tfevalfit(df0,d,data,n,m,ke,symdv,nvar,sav,v,
     $       kaxvec,.false.,cut,irtc)
        if(irtc /= 0)then
          deallocate(a0,a,abest,df,df0,v00,w,cv,vbest,dv,df2)
          call tflocal1(kaxvec)
          return
        endif
        if(d < d0)then
          v0=v
          if(d < dbest)then
            vbest=v0
            abest=a0
            dbest=d
          endif
          d0=d
          fact=min(1.d0,fact*4.d0)
          fact1=0.d0
          d1=d
          if(.not. newton)then
            newton=d < 0.99d0*d00
            if(newton)then
              go to 11
            endif
          endif
          go to 1
        endif
        iter=iter+1
        if(iter <= maxi)then
          d2=d1
          d1=d
          fact2=fact1
          fact1=fact
          fact=fmincube(fact1,fact2,d1,d2,d0,dg)
          if(newton)then
            if(fact >= factmin)then
              go to 2
            endif
            newton=.false.
            go to 11
          else
            if(fact >= factmin1)then
              go to 2
            else
              v0=max(vmin,min(vmax,2.d0*v00-v0))
c     do i=1,nvar
c     v0(i)=max(vmin(i),min(vmax(i),2.d0*v00(i)-v0(i)))
c     enddo
              go to 21
            endif
          endif
        endif
      endif
      v0=vbest
      a0=abest
      d0=dbest
      if(cut /= 0.d0)then
        do j=1,m
          if(abs(df(j)) > cut)then
            a0(j,:)=0.d0
          endif
        enddo
      endif
      call tsvdm(a0,0.d0,w,m,nvar,m,0,svdeps,.true.)
      if(n == 2)then
        sigma=sqrt(d0/max(1,m-nvar))
        do i=1,nvar
          a0(i,:)=a0(i,:)*w(i)*sigma
        enddo
      else
        do i=1,nvar
          a0(i,:)=a0(i,:)*w(i)
        enddo
      endif
      call tflocal1(kaxvec)
      kdm=kxm2l(abest,m,nvar,m,.false.)
      call tfcovmat(a0,cv,m,nvar,m)
      kcv=kxm2l(cv,nvar,nvar,nvar,.false.)
      kci=kxraaloc(-1,nvar,klci)
      chin=tinvgr(dble(nvar))
      do i=1,nvar
        x=sqrt(cv(i,i)*chin)
        klci%rbody(i)=merge(x,rtfnan,ktfrealq(transfer(x,i00)))
      enddo
      deallocate(a0,a,abest,df,df0,v00,w,cv,vbest,dv,df2)
      return
      end

      subroutine tfevalfit(df,r,a,n,m,ke,symdv,nvar,sav,v,
     $     kaxvec,deriv,cutoff,irtc)
      use eeval
      implicit none
      type (sad_symdef) ,intent(inout):: symdv
      type (sad_dlist), pointer :: klx,kl1
      type (sad_descriptor) k1,ke,kx
      integer*4 ,intent(out):: irtc
      integer*4 ,intent(in):: n,m,nvar
      integer*4 i,l,itfmessage
      type (symv) sav(nvar)
      integer*8 kavvec,kaxvec
      real*8 ,intent(out):: df(m),r
      real*8 ,intent(in):: a(n,m),v(nvar),cutoff
      real*8 vx
      logical*4 deriv
      l=itfuplevel()
      do i=1,nvar
        sav(i)%p%value=dfromr(v(i))
      enddo
      symdv%value%k=ktflist+ktfcopy1(kaxvec)
      kavvec=ktfaddr(klist(kaxvec+1))
      do i=1,m
        rlist(kavvec+i)=a(1,i)
      enddo
      kx=tfeevalref(ke,irtc)
      call tflocal1(kaxvec)
      symdv%value%k=0
      if(irtc /= 0)then
        return
      endif
      if(ktflistq(kx,klx))then
        if(klx%head%k == kxvect .and. klx%nl == 1)then
          k1=klx%dbody(1)
          if(tfcomplexnumlistqk(k1%k,kl1) .and. kl1%nl == m)then
            if(deriv)then
              if(ktfnonreallistqo(kl1))then
                do i=1,m
                  df(i)=merge(kl1%rbody(i),0.d0,ktfrealq(kl1%dbody(i)))
                enddo
              else
                do i=1,m
                  df(i)=kl1%rbody(i)
                enddo
              endif
            else
              if(ktfnonreallistqo(kl1))then
                do i=1,m
                  if(ktfrealq(kl1%dbody(i)))then
                    df(i)=kl1%rbody(i)-a(2,i)
                  else
                    df(i)=alarge
                  endif
c                  df(i)=merge(kl1%rbody(i)-a(2,i),alarge,
c     $                 ktfrealq(kl1%dbody(i)))
                enddo
              else
                do i=1,m
                  df(i)=kl1%rbody(i)-a(2,i)
                enddo
              endif
            endif
            go to 1000
          endif
        endif
      elseif(ktfrealq(kx,vx))then
        if(deriv)then
          df=vx
        else
          df=vx-a(2,:)
        endif
c        df=merge(vx,vx-a(2,:),deriv)
        go to 1000
      endif
      do i=1,m
        symdv%value=dfromr(a(1,i))
        kx=tfeevalref(ke,irtc)
        if(irtc /= 0)then
          go to 9000
        endif
        if(ktfnonrealq(kx))then
          if(tfcomplexq(kx))then
            df(i)=merge(0.d0,alarge,deriv)
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
          vx=kx%x(1)
          if(deriv)then
            df=vx
          else
            df=vx-a(2,i)
          endif
c          df(i)=merge(vx,vx-a(2,i),deriv)
        endif
      enddo
 1000 r=0.d0
      if(deriv)then
        if(n == 3)then
          df=df/a(3,:)
        elseif(n == 4)then
          df=df/a(4,:)
        endif
      else
        if(cutoff /= 0.d0)then
          if(n == 3)then
            do i=1,m
              df(i)=df(i)/a(3,i)
              r=r+min(cutoff,max(-cutoff,df(i)))**2
            enddo
          elseif(n == 4)then
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
          if(n == 3)then
            do i=1,m
              df(i)=df(i)/a(3,i)
              r=r+df(i)**2
            enddo
          elseif(n == 4)then
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
      end

      subroutine tfderiv(ke,nvar,sav,kdl,irtc)
      use eeval
      implicit none
      type (sad_descriptor) ,intent(in):: ke
      type (sad_descriptor) ,intent(out):: kdl(nvar)
      type (sad_descriptor) kd,kr,ks
      integer*4 ,intent(in):: nvar
      integer*4 ,intent(out):: irtc
      integer*4 i,isp0
      type (symv) ,intent(in):: sav(nvar)
      integer*8 ierr0
      logical*4 euv
      if(iads%k == 0)then
        iads=kxsymbolz('System`D',8)
        ks=tfsyeval(iads,irtc)
        if(irtc /= 0)then
          return
        endif
        iader=kxsymbolz('CheckDerivative',15)
      endif
      do i=1,nvar
        sav(i)%p%value=sad_descr(sav(i)%p%sym)
      enddo
      isp0=isp
      isp=isp0+1
      dtastk(isp)=iads
      do i=1,nvar
        isp=isp0+2
        dtastk(isp)=ke
        isp=isp+1
        dtastk(isp)=sad_descr(sav(i)%p%sym)
        ierr0=ierrorexp
        ierrorexp=1
c     call tfdebugprint(ke,'tfderiv D',2)
c     call tfdebugprint(dtastk(isp),'by ',2)
        call tfdeval(isp0+1,iads,kd,1,.false.,euv,irtc)
c     call tfdebugprint(kd,'==>',2)
        ierrorexp=ierr0
        if(irtc /= 0)then
          if(irtc > 0. and. ierrorprint /= 0)then
            call tfreseterror
          endif
          go to 100
        endif
        isp=isp0+2
        dtastk(isp)=iader
        isp=isp+1
        dtastk(isp)=dtfcopy(kd)
        ierr0=ierrorexp
        ierrorexp=1
        call tfdeval(isp0+2,iader,kr,1,.false.,euv,irtc)
        ierrorexp=ierr0
        if(irtc == 0 .and. ktfrealq(kr))then
          kdl(i)=kd
        else
          call tflocald(kd)
          kdl(i)%k=ktfoper+mtfnull
          if(irtc > 0. and. ierrorprint /= 0)then
            call tfreseterror
          endif
        endif
      enddo
 100  do i=1,nvar
        sav(i)%p%value%k=0
      enddo
      irtc=0
      isp=isp0
c     write(*,*)'tfderiv-end'
      return
      end

      real*8 function tinvgr(an)
      use gammaf
      implicit none
      real*8 ,intent(in)::an
      real*8 c0,x0,gn,erfc,df,dfdx,anh
      if(an == 1.d0)then
        tinvgr=1.d0
        return
      endif
      c0=erfc(sqrt(0.5d0))
      anh=an*.5d0
      x0=anh
      gn=gamma(anh)
      df=gammaq(anh,x0)-c0
      do while(abs(df) > 1.d-14)
        dfdx=-exp(-x0)*x0**(anh-1.d0)/gn
        x0=x0-df/dfdx
c     write(*,*)'tinvgr ',x0,dfdx,df
        df=gammaq(anh,x0)-c0
      enddo
      tinvgr=x0*2.d0
      return
      end

c      subroutine tfmakerulestk(k1,k2)
c      implicit none
c      type (sad_descriptor) ,intent(in)::k1,k2
c     call tfdebugprint(k1,'mkrs',1)
c     write(*,*)k2%k
c     call tfdebugprint(k2,'mkrs_f',1)
c      call mrs(k1,k2)
c      return
c      end

      end module
