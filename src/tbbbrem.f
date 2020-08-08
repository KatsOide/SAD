c
c  Subroutine for BBBrem, https://doi.org/10.1016/0010-4655(94)90085-X
c  R. Kleiss and H. Burkhardt, Comput.Phys.Commun. 81 (1994) 372-380   DOI 10.1016/0010-4655(94)90085-X
c  based on bbbrem.f sent by H. Burkhardt on Jul 1, 2019.
c
c  modified for SAD: K. Oide 20 Jul 2019.
c
      module bbbrem
      use macmath
      use macphys
      type labmom
      real*8 p2(4),q2(4),qk(4),weight
      end type
      real*8, parameter::alpha=finest,rme=elmass/1e9,
     $     tomb=(plankr*cveloc/elemch)**2/1.d18/1.d-31,
     $     twopi=m_2pi,rme2=rme**2
c     use SAD's tran()
      integer*4 :: nran=3

      contains
      subroutine bcube1(roots,rk0,mom,sigapp)
      implicit none
      type (labmom) mom
      real*8 roots,rk0
      real*8
c      double precision
c     ,p1,p2,q1,q2,qk,alpha,rme,tomb,pi,twopi,rme2,d1,d2,t,weight
     . p1(4),q1(4),d1,d2,t,
     $     s,rme2s,rls,z0,a1,a2,ac,sigapp,eb,pb,rin2pb,
     . z,y,q0,temp1,tmin,tmax,sy,w2,rlamx,b,rlam,eps,rl,vgam,
     . cgam,sgam,phi,phig,ql,qt,q(4),r0,w,rin2w,rinr0w,eta,phat1,
     . phat2,phat3,phatt,phatl,sfhat,cfhat,sthat,vthat,cthat,sfg,
     . cfg,temp2,veg,qkhat(4),c1,vnumx,vdenx,vnumn,vdenn,c2,rlabl,
     . rhat4,etar,rhat1,rhat2,rhat3,zz,s1,rind1,rind12,
     . rind2,rind22,aa0,aa1,aa2,aa3,aa4,rmex,rmap
c     always do init when called from SAD
c      integer*4, save:: init=0

c* experimental constants
c      common / experc / roots,rk0
c* momenta in the lab frame
c      common / labmom / p1(4),p2(4),q1(4),q2(4),qk(4),d1,d2,t,weight

c      data init/0/
* start of initialization
c      if(init.eq.0) then
c        init = 1
c        write(*,901)
c  901   format(' ******************************************'/,
c     .         ' *** bbbrem : beam-beam bremsstrahlung  ***'/,
c     .         ' *** authors r. kleiss and h. burkhardt ***'/,
c     .         ' ******************************************')

* physical constants
c        alpha = 1d0/137.036d0
c        rme = 0.51099906d-03
c        tomb = 3.8937966d05 /1d6

* mathematical constants
c        pi = 4*datan(1d0)
c        twopi = 2*pi

* derived constants
c        write(*,902) roots,rk0
c  902   format(' total energy              ',f15.6,' gev'/,
c     .         ' minimum photon energy     ',f15.6,' * beam energy')
c        if(rk0.le.0d0.or.rk0.ge.1d0) then
c          write(*,*) 'wrong value for rk0'
c          stop
c        endif
        s = roots**2
c        rme2 = rme**2
        rme2s = rme2/s
        rls = -dlog(rme2s)
        z0 = rk0/(1-rk0)

* approximate total cross section
        a1 = dlog((1+z0)/z0)
        a2 = (dlog(1+z0))/z0
        ac = a1/(a1+a2)
        sigapp = 8*alpha**3/rme2*(-dlog(rme2s))*(a1+a2)*tomb
c        write(*,903) sigapp
c  903   format(' approximate cross section ',d15.6,' millibarn')

* the initial-state momenta
        eb = roots*0.5d0
        pb = dsqrt(eb*eb-rme2)
        rin2pb = 0.5d0/pb
        p1(1) = 0
        p1(2) = 0
        p1(3) = -pb
        p1(4) = eb
        q1(1) = 0
        q1(2) = 0
        q1(3) = pb
        q1(4) = eb

* end of initialization
c      endif

* generate z
      if(random(1.d0).lt.ac) then
         temp1=random(2.d0)
         z = 1d0/(temp1*(dexp(a1*random(3.d0))-1))
      else
        z = z0/random(4.d0)
      endif

* bounds on t
      y = rme2s*z
      q0 = eb*y
      temp1 = pb*pb-eb*q0
c    correction by hbu sent by email on Jul. 19, 2019
c      temp2 = temp1*temp1-rme2-q0*q0
      temp2 = temp1*temp1-rme2*q0*q0
* exit if temp2<0 (very very rare): put weight to 0
* the `else' clause extends to the end of the routine
      if(temp2.lt.0d0) then
        write(*,904) temp2
  904   format(' y too large: delta_t^2 = ',d15.6)
        mom%weight = 0d0
      else
        tmin = -2*(temp1+dsqrt(temp2))
        tmax = rme2*s*y*y/tmin

* generate t
        sy = s*y
        w2 = sy+rme2
        temp1 = sy+tmax
        rlamx = dsqrt(temp1*temp1-4*w2*tmax)
        if(temp1.le.0d0) then
          temp1 = rlamx-temp1
        else
          temp1 = -4*w2*tmax/(rlamx+temp1)
        endif
    1   continue
          b = dexp(random(5.d0)*dlog(1+2*sy/temp1))
          t = -b*z*z*rme2/((b-1)*(b*z+b-1))
          if(t.lt.tmin) then
            write(*,905) t,tmin
  905       format(' t = ',d15.6,'   < t_min =',d15.6)
            goto 1
          endif
        continue

* generate cgam
        rlam = dsqrt((sy-t)*(sy-t)-4*rme2*t)
        eps = 4*rme2*w2/(rlam*(rlam+w2+rme2-t))
        rl = dlog((2+eps)/eps)
        vgam = eps*(dexp(random(6.d0)*rl)-1)
        cgam = 1-vgam
        sgam = dsqrt(vgam*(2-vgam))

* generate azimuthal angles
        phi = twopi*random(7.d0)
        phig = twopi*random(8.d0)
c centralize phi : 3 Aug 2019 K. Oide
c        phi=twopi*(random(7.d0)-0.5d0)
c        phig=twopi*(random(8.d0)-0.5d0)

* construct momentum transfer q(mu)
        ql = (2*eb*q0-t)*rin2pb
        qt = dsqrt((tmax-t)*(t-tmin))*rin2pb
        q(1) = qt*dsin(phi)
        q(2) = qt*dcos(phi)
        q(3) = ql
        q(4) = q0

* construct momentum of outgoing positron in lab frame
        mom%q2(1) = q1(1)-q(1)
        mom%q2(2) = q1(2)-q(2)
        mom%q2(3) = q1(3)-q(3)
        mom%q2(4) = q1(4)-q(4)

* find euler angles of p1(mu) in cm frame
        r0 = eb+q0
        w = dsqrt(w2)
        rin2w = 0.5d0/w
        rinr0w = 1d0/(r0+w)
        eta = -(sy+2*w*q0+t)*rin2w*rinr0w
        phat1 = -q(1)*(1+eta)
        phat2 = -q(2)*(1+eta)
        phat3 = pb*eta-ql*(1+eta)
        phatl = rlam*rin2w
        phatt = dsqrt(phat1*phat1+phat2*phat2)
        sfhat = phat1/phatt
        cfhat = phat2/phatt
        sthat = phatt/phatl
        if(phat3.gt.0d0) then
          vthat = sthat*sthat/(1-dsqrt(1-sthat*sthat))
        else
          vthat = sthat*sthat/(1+dsqrt(1-sthat*sthat))
        endif
        cthat = vthat-1

* rotate using these euler angles to get the qk direction in the cm
        sfg = dsin(phig)
        cfg = dcos(phig)
        temp1 = sgam*sfg
        temp2 = cthat*sgam*cfg+sthat*cgam
        veg = vthat+vgam-vthat*vgam-sthat*sgam*cfg
        qkhat(4) = sy*rin2w
        qkhat(1) = qkhat(4)*( cfhat*temp1+sfhat*temp2)
        qkhat(2) = qkhat(4)*(-sfhat*temp1+cfhat*temp2)
        qkhat(3) = qkhat(4)*(veg-1)

* boost the photon momentum to the lab frame
        temp1 = pb*qkhat(3)
        if(temp1.gt.0d0) then
          temp2 = (rme2*qkhat(4)*qkhat(4)
     .            +pb*pb*(qkhat(1)*qkhat(1)+qkhat(2)*qkhat(2)))
     .            /(eb*qkhat(4)+temp1)
        else
          temp2 = eb*qkhat(4)-temp1
        endif
        mom%qk(4) = (temp2+qkhat(4)*q(4)+qkhat(1)*q(1)
     .                +qkhat(2)*q(2)+qkhat(3)*q(3))/w
        temp1 = (mom%qk(4)+qkhat(4))*rinr0w
        mom%qk(1) = qkhat(1)+temp1*q(1)
        mom%qk(2) = qkhat(2)+temp1*q(2)
        mom%qk(3) = qkhat(3)+temp1*(-pb+ql)

* construct p2 by momentum conservation
        mom%p2(1) = -mom%q2(1)-mom%qk(1)
        mom%p2(2) = -mom%q2(2)-mom%qk(2)
        mom%p2(3) = -mom%q2(3)-mom%qk(3)
        mom%p2(4) = -mom%q2(4)-mom%qk(4)+roots

* impose cut on the photon energy: qk(4)>eb*rk0
        if(mom%qk(4).lt.eb*rk0) then
          mom%weight = 0d0
        else

* the event is now accepted: compute matrix element and weight
* compute fudge factor c1
          c1 = dlog(1+z)/dlog((2+eps)/eps)

* compute fudge factor c2
          temp1 = sy-tmax
          vnumx = dsqrt(temp1*temp1-4*rme2*tmax)+temp1
          temp1 = sy+tmax
          if(temp1.lt.0d0) then
            vdenx = dsqrt(temp1*temp1-4*w2*tmax)-temp1
          else
            vdenx = -4*w2*tmax/(dsqrt(temp1*temp1-4*w2*tmax)+temp1)
          endif
          temp1 = sy-tmin
          vnumn = dsqrt(temp1*temp1-4*rme2*tmin)+temp1
          temp1 = sy+tmin
          if(temp1.lt.0d0) then
            vdenn = dsqrt(temp1*temp1-4*w2*tmin)-temp1
          else
            vdenn = -4*w2*tmin/(dsqrt(temp1*temp1-4*w2*tmin)+temp1)
          endif
          c2 = 2*rls/dlog((vnumx*vdenn)/(vdenx*vnumn))

* compute vector (small) r in cm frame, and (big) z
          rlabl = (t-2*rme2*y)*rin2pb
          rhat4 = -(2*rme2*y+(1-y)*t)*rin2w
          etar = rhat4*rinr0w
          rhat1 = -q(1)*(1+etar)
          rhat2 = -q(2)*(1+etar)
          rhat3 = rlabl+(pb-ql)*etar
          zz = s*(rhat4*qkhat(4)-rhat1*qkhat(1)
     .         -rhat2*qkhat(2)-rhat3*qkhat(3))

* the other invariants
          s1 = 4*eb*(eb-mom%qk(4))
          d1 = sy*rlam*(eps+vgam)*rin2w*rin2w
          d2 = 0.5d0*sy

* the exact matrix element
          rind1 = 1d0/d1
          rind12 = rind1*rind1
          rind2 = 1d0/d2
          rind22 = rind2*rind2
c* kleiss-burkhardt cross section multiplied by t**2
          temp1 = s+t-2*d2
          temp2 = s1+t+2*d1
          aa0 = (s*s+s1*s1+temp1*temp1+temp2*temp2)*rind1*rind2*(-t)
          aa1 = -4*rme2*zz*zz*rind12*rind22
          aa2 = -8*rme2*(d1*d1+d2*d2)*rind1*rind2
          aa3 = 16*rme2*rme2*(d1-d2)*zz*rind12*rind22
          aa4 = -16*rme2*rme2*rme2*(d1-d2)**2*rind12*rind22
          rmex = aa0+aa1+aa2+aa3+aa4

* the approximate matrix element without c1,2, multiplied by t**2
          rmap = 4*s*s*rind1*rind2*(-t)*c1*c2

* the weight
          mom%weight = rmex/rmap*sigapp
CHBU      IF(T.GT.-3.56D-21) weight=0.D0  ! CHBU cutoff

* the weight is now defined for both accepted and rejected events
        endif
      endif
c 1001 format(3(a10,d15.6))
      return
      end subroutine

      real*8 function random(r)
      implicit none
      real*8 r,tran
c      common / ranchn / nran
c      save
c      data init/0/
c      if(init.eq.0) then
c        init = 1
c        if(nran.eq.1) then
c          write(*,*) ' random: very crude random numbers chosen'
c        elseif(nran.eq.2) then
c          write(*,*) ' random: quasi-random numbers chosen'
c        else
c          write(*,*) ' random: error: nran = ',nran,' not allowed'
c          stop
c        endif
c      endif
      if(nran.eq.1) then
        random = rando1(r)
      elseif(nran.eq.2) then
        random = rando2(r)
      elseif(nran .eq. 3)then
c   use SAD's tran()
        random=tran()
      endif
      return
      end function

      real*8 function rando1(r)
      implicit none
      save
      real*8 r,rm
      integer*4 :: init=0
      integer*4 m,ia,ic,iseed
      if(init.eq.0) then
        init = 1
        m = 2**20
        rm = 1d0/(1d0*m)
        ia = 2**9+5
        ic = 1
        iseed = 100461
      endif
    1 iseed = mod(ia*iseed+ic,m)
      if(iseed.eq.0) goto 1
      rando1 = iseed*rm
      return
      end function

      real*8 function rando2(r)
      use macmath, nprime=>smallprime
      implicit none
      save
      integer*4 , parameter::n = 3
      integer*4  init,k,j
      real*8 s(10),ss,r
      data init/0/
      if(init.eq.0) then
        init = 1
        do 1 k = 1,n
          s(k) = dsqrt(nprime(k)*1d0)
    1   continue
      endif
      do 3 k = 1,n
        ss = 0
        do 2 j=1,n
          ss = ss+s(j)
    2   continue
        s(k) = dmod(ss,1d0)
    3 continue
      rando2 = s(1)
      return
      end function

      end module

      subroutine tbcube1(isp1,kx,irtc)
c  interface to function BBBrem1
      use tfstk
      use bbbrem
      implicit none
      integer*4 isp1,irtc,itfmessage
      type (sad_descriptor) kx
      type (labmom) mom
      type (sad_dlist), pointer:: kl
      type (sad_rlist), pointer:: klp
      real*8 ars,drs,sigapp
      real*8 , parameter :: gev=1.d9,mbarn=1.d-31
      if(isp .ne. isp1+2)then
        irtc=itfmessage(9,'General::narg','"2"')
        return
      endif
      if(.not. ktfrealq(dtastk(isp1+1),ars))then
        irtc=itfmessage(9,'General::wrongtype','"Sqrt[s] (eV) for #1"')
        return
      endif
      if(.not. ktfrealq(dtastk(isp1+2),drs))then
        irtc=itfmessage(9,'General::wrongtype','"Min. ds/s for #2"')
        return
      endif
      if(drs .le. 0.d0 .or. drs .ge. 1.d0)then
        irtc=itfmessage(9,'General::wrongval','"0 < #2 <1","ds/s"')
        return
      endif
      call bcube1(ars/gev,drs,mom,sigapp)
      kx=kxaaloc(-1,2,kl)
      kl%dbody(1)=kxraaloc(0,4,klp)
      klp%rbody(1:4)=(/-1.d0,-1.d0,-1.d0,1.d0/)*mom%p2*gev
      kl%dbody(2)=dfromr(sigapp*mbarn)
      irtc=0
      return
      end
