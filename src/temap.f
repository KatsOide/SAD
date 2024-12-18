      subroutine temap(np,np0,x,px,y,py,z,g,dv,sx,sy,sz,l,nt,kptbl)
      use tfstk
      use efun
      use temw, only:tmulbs
      use ffs_flag, only:calpol
      use tfcsi,only:icslfnm,lfnm
      use mathfun
      implicit none
      type alist
        type (sad_rlist), pointer :: p
      end type
      type (alist) kav(9)
      type (sad_descriptor) kx
      type (sad_dlist), pointer :: klx,kl
      type (sad_rlist), pointer :: klrk
      integer, parameter :: nkptbl = 6
      integer*4 ,intent(in):: np0,l,nt
      integer*4 ,intent(inout):: kptbl(np0,nkptbl),np
      integer*4 i,j,k,isp0,irtc,m,kptmp(nkptbl),nc
      real*8 ,intent(inout):: x(np0),px(np0),y(np0),py(np0),
     $     z(np0),g(np0),dv(np0),sx(np0),sy(np0),sz(np0)
      character*2 ord
      real*8 xa(9),pr,sr,phir
      integer*8,save:: iem=0,ifv
      logical*4 dodrop,doinject
      if(itfcontext .le. 0)then
        return
      endif
      if(iem == 0)then
        iem=ktfsymbolz('ExternalMap',11)
        ifv=ktsalocb(0,'TRACK',5)
      endif
      nc=merge(9,7,calpol)
      call tclrfpe
      levele=levele+1
      isp0=isp
      isp=isp+1
      ktastk(isp)=ktfsymbol+iem
      isp=isp+1
      ktastk(isp)=ktfstring+ifv
      isp=isp+1
      rtastk(isp)=l
      isp=isp+1
      rtastk(isp)=nt
      isp=isp+1
      dtastk(isp)=kxadaloc(-1,nc,kl)
      do k=1,nc
        kl%dbody(k)=kxavaloc(0,np0,kav(k)%p)
      enddo
      do i=1,np0
        j=kptbl(i,1)
        if((j .gt. np) .or. (kptbl(j,4) /= 0))then
          do concurrent (k=1:nc)
            kav(k)%p%rbody(i)=0.d0
          enddo
        else
          pr=g(j)
          kav(1)%p%rbody(i)=x(j)
          kav(2)%p%rbody(i)=px(j)*(1.d0+pr)
          kav(3)%p%rbody(i)=y(j)
          kav(4)%p%rbody(i)=py(j)*(1.d0+pr)
          kav(5)%p%rbody(i)=z(j)
          kav(6)%p%rbody(i)=pr
          kav(nc)%p%rbody(i)=1.d0
          if(calpol)then
            kav(7)%p%rbody(i)=sy(j)
            kav(8)%p%rbody(i)=atan(sz(j),sx(j))
          endif            
        endif
      enddo
      kx=tfefunref(isp0+1,.false.,irtc)
      if(irtc /= 0)then
        levele=itfdownlevel()
        isp=isp0
        if(ierrorprint /= 0)then
          call tfaddmessage(' ',2,icslfnm())
        endif
        write(lfnm,*)' Error in ExternalMap of ',l,ord(l),' element at ',
     $       nt,ord(nt),' turn.'
        return
      elseif(.not. tflistq(kx,klx) .or. klx%nl /= 7)then
        go to 9000
      endif
      do k=1,nc
        if(.not. tfreallistq(klx%dbody(k),klrk))then
          go to 9000
        endif
        if(klrk%nl /= np0)then
          go to 9000
        endif
        kav(k)%p=>klrk
      enddo
      dodrop=.false.
      doinject=.false.
      do i=1,np0
        j=kptbl(i,1)
        if(kav(nc)%p%rbody(i) /= 0.d0)then
          if(.not. (j .le. np)) then
c     Case: dropped before MAP element
            doinject=.true.
          endif
c     Re-activate particle slot[j] if kptbl(j,4) != 0
c     Note: kptbl(j,4) MUST be `0' for alive particles
          kptbl(j,4)=0
          kptbl(j,5)=nt+1
c     Copy-in ExternalMap[] result for alive/reinject case
          xa(1)=kav(1)%p%rbody(i)
          xa(2)=kav(2)%p%rbody(i)
          xa(3)=kav(3)%p%rbody(i)
          xa(4)=kav(4)%p%rbody(i)
          xa(5)=kav(5)%p%rbody(i)
          xa(6)=kav(6)%p%rbody(i)
          call tconv(xa,xa,-1)
          x (j)=xa(1)
          px(j)=xa(2)
          y (j)=xa(3)
          py(j)=xa(4)
          z (j)=xa(5)
          g (j)=xa(6)
          dv(j)=xa(7)
          if(calpol)then
            sy(j)=kav(7)%p%rbody(i)
            sr=1.d0+sqrt1(-sy(j)**2)
            phir=kav(8)%p%rbody(i)
            sx(i)=sr*cos(phir)
            sz(i)=sr*sin(phir)
          endif
        elseif((j .le. np) .and. (kptbl(j,4) == 0))then
c     Lose particle slot[j] at current MAP element[l]
          dodrop=.true.
          kptbl(j,4)=l
          kptbl(j,5)=nt
        endif              
      enddo

      if(doinject)then
        i=np+1
        m=np0
        do while(i .le. m)
          if(kptbl(i,4) /= 0)then
c     Search alive paricle from tail: (i, m]
            do while((i .lt. m) .and. (kptbl(m,4) /= 0))
              m=m-1
            enddo
            if(kptbl(m,4) == 0)then
c     Swap dead particle slot[i] with tail alive particle slot[m]
              j=kptbl(m,2)
              k=kptbl(i,2)
c     - Update maps between particle ID and array index
              kptbl(k,1)=m
              kptbl(j,1)=i
              kptbl(m,2)=k
              kptbl(i,2)=j
c     - Swap kptbl except forward/backward[kptbl(*,1)/kptbl(*,2)]
              kptmp(  3:nkptbl) = kptbl(m,3:nkptbl)
              kptbl(m,3:nkptbl) = kptbl(i,3:nkptbl)
              kptbl(i,3:nkptbl) = kptmp(  3:nkptbl)
c     - Swap particle coordinates
              xa(1)=x (i)
              xa(2)=px(i)
              xa(3)=y (i)
              xa(4)=py(i)
              xa(5)=z(i)
              xa(6)=g (i)
              xa(7)=dv(i)

              x (i)=x (m)
              px(i)=px(m)
              y (i)=y (m)
              py(i)=py(m)
              z (i)=z (m)
              g (i)=g (m)
              dv(i)=dv(m)

              x (m)=xa(1)
              px(m)=xa(2)
              y (m)=xa(3)
              py(m)=xa(4)
              z (m)=xa(5)
              g (m)=xa(6)
              dv(m)=xa(7)

              if(calpol)then
                xa(1)=sx(i)
                xa(2)=sy(i)
                xa(3)=sz(i)
                
                sx(i)=sx(m)
                sy(i)=sy(m)
                sz(i)=sz(m)

                sx(m)=xa(1)
                sy(m)=xa(2)
                sz(m)=xa(3)
              endif
            endif
            m=m-1
          endif
          i=i+1
        enddo
        np=m
      endif

      if(dodrop)then
      else
c        call tfdebugprint(kx,'temap',1)
      endif
 9000 levele=itfdownlevel()
      isp=isp0
      return
      end

      subroutine temapp(isp1,kx,irtc)
      use tfstk
      use tftr
      use maloc,only:tfmsize
      implicit none
      type (sad_descriptor) ,intent(out):: kx
      integer*4 ,intent(in):: isp1
      integer*4 ,intent(out):: irtc
      type (sad_dlist), pointer :: kzl,kxl
      type (sad_rlist), pointer :: kzi,kpl,kpri
      real*8 xa
      integer*8 ka,ka0
      integer*4 n,np,nprm,i,itfmessage
      if(isp /= isp1+3)then
        irtc=itfmessage(9,'General::narg','"3"')
        return
      endif
      if(.not. ktfrealq(dtastk(isp1+1),xa))then
        irtc=itfmessage(9,'General::wrongtype','"address for #1"')
        return
      endif
      call tfmsize(dtastk(isp1+2),n,np,irtc)
      if(irtc /= 0)then
        return
      endif
      if(n /= 7 .and. n /= 9)then
        irtc=itfmessage(9,'General::wrongtype','"particles for #2"')
        return
      endif
      if(.not. tfreallistq(dtastk(isp1+3),kpl))then
        irtc=itfmessage(9,'General::wrongtype','"parameters for #3"')
        return
      endif
      call descr_dlist(dtastk(isp1+2),kzl)
      nprm=kpl%nl
      ka0=int8(xa)/8
      if(ka0 /= 0)then
        klist(ka0)=npz
        klist(ka0+1)=nprm
        ka=ka0+ipn+2
        do i=1,n-1
          call descr_rlist(kzl%dbody(i),kzi)
          rlist(ka+(i-1)*npz:ka+(i-1)*npz+np-1)=kzi%rbody(1:np)
        enddo
        call descr_rlist(kzl%dbody(n),kzi)
        rlist(ka+8*npz:ka+8*npz+np-1)=kzi%rbody(1:np)
        ka=ka0+9*npz+2+npri*nprm
        rlist(ka:ka+nprm-1)=kpl%rbody(1:nprm)
      endif
      kx=kxadaloc(-1,3,kxl)
      kxl%rbody(1)=dble(ipn)
      kxl%rbody(2)=dble(npri)
      if(iprid == 0)then
        kxl%rbody(3)=dble(dble(getpid()))
      else
        kxl%dbody(3)=kxavaloc(0,npr,kpri)
        kpri%rbody(1:npr)=dble(ipr(1:npr))
      endif
      return
      end

      function teunmapp(isp1,irtc) result(kx)
      use tfstk
      use tftr
      use ffs_flag, only:calpol
      implicit none
      type (sad_descriptor) kx
      integer*4 ,intent(in):: isp1
      integer*4 ,intent(out):: irtc
      type (sad_dlist), pointer :: kzl,kl
      type (sad_rlist), pointer :: kzi,kpl
      real*8 xa
      integer*8 ka,ka0
      integer*4 nprm,i,itfmessage,n,np
      kx=dxnullo
      if(isp /= isp1+2)then
        irtc=itfmessage(9,'General::narg','"2"')
        return
      endif
      if(.not. ktfrealq(dtastk(isp1+1),xa))then
        irtc=itfmessage(9,'General::wrongtype','"address for #1"')
        return
      endif
      if(.not. ktfrealq(dtastk(isp),np))then
        irtc=itfmessage(9,'General::wrongtype',
     $       '"number of particles for #2"')
        return
      endif
      ka0=int8(xa)/8
      nprm=int(klist(ka0+1))
      kx=kxadaloc(-1,2,kl)
      n=merge(8,6,calpol)
      kl%dbody(1)=kxadaloc(0,n+1,kzl)
      kl%dbody(2)=kxavaloc(0,nprm,kpl)
      ka=ka0+ipn+2
      do i=1,n
        kzl%dbody(i)=kxavaloc(0,np,kzi)
        kzi%rbody(1:np)=rlist(ka+(i-1)*npz:ka+(i-1)*npz+np-1)
      enddo
      kzl%dbody(n+1)=kxavaloc(0,np,kzi)
      kzi%rbody(1:np)=rlist(ka+8*npz:ka+8*npz+np-1)
c      write(*,*)'unmap ',ipn,n,np,kzi%rbody(4)
      ka=ka0+9*npz+2+npr*nprm
      kpl%rbody(1:nprm)=rlist(ka:ka+nprm-1)
      irtc=0
      return
      end

      subroutine temape(trans,cod,beam,l)
      use tfstk
      use tmacro
      use sad_main, only:iaidx
      use temw,only:tmulbs
      use efun
      use maloc,only:ktfmalocp
      use sad_basics
      use iso_c_binding
      use tfcsi,only:icslfnm,lfnm
      implicit none
      type (sad_descriptor) kx,k1,k2,k3,k4
      type (sad_dlist) ,pointer ::klx,kl1
      integer*8 kat1,kbm,krt
      integer*4 l,isp0,n,m,irtc,i,j
      real*8 trans(6,6),cod(6),beam(42)
      real*8 ,pointer::trat1(:,:)
      character*2 ord
      integer*8 , save :: ifv=0,iem=0
c      iaidx(m,n)=((m+n+abs(m-n))**2+2*(m+n)-6*abs(m-n))/8
      if(iem == 0)then
        iem=ktfsymbolz('`ExternalMap',12)
        ifv=ktsalocb(0,'EMIT',4)
      endif
      levele=levele+1
      isp0=isp
      isp=isp+1
      ktastk(isp)=ktfsymbol+iem
      isp=isp+1
      ktastk(isp)=ktfstring+ifv
      isp=isp+1
      rtastk(isp)=l
      isp=isp+1
      dtastk(isp)=kxm2l(cod,0,6,1,.false.)
      kx=tfefunref(isp0+1,.false.,irtc)
      if(irtc /= 0)then
        levele=itfdownlevel()
        if(ierrorprint /= 0)then
          call tfaddmessage(' ',2,icslfnm())
        endif
        write(lfnm,*)' Error in ExternalMap(EMIT) of ',l,ord(l),' element.'
        isp=isp0
        return
      endif
      if(tfnonlistq(kx,klx))then
        go to 9000
      endif
      if(klx%nl == 3 .and. klx%body(1) == ktfsymbol+iem)then
        go to 9000
      endif
      k1=klx%dbody(1)
      if(ktfnonlistq(k1,kl1))then
        go to 9100
      endif
      if(kl1%nl /= 6 .or. ktfnonreallistqo(kl1))then
        go to 9100
      endif
      cod=kl1%rbody(1:6)
      k2=klx%dbody(2)
      kat1=ktfmalocp(k2,n,m,.false.,.false.,
     $     .false.,.false.,irtc)
      if(irtc /= 0)then
        go to 9100
      endif
      if(n /= 6 .or. m /= 6)then
        call tfree(kat1)
        go to 9100
      endif
      call c_f_pointer(c_loc(rlist(kat1)),trat1,[6,6])
      kbm=0
      krt=0
      if(klx%nl /= 2)then
        if(klx%nl == 4)then
          k3=klx%dbody(3)
          krt=ktfmalocp(k3,n,m,.false.,.false.,
     $         .false.,.false.,irtc)
          if(irtc /= 0)then
            go to 9110
          elseif(n /= 6 .or. m /= 6)then
            go to 9120
          endif
          k4=klx%dbody(4)
          kbm=ktfmalocp(k4,n,m,.false.,.false.,
     $         .false.,.false.,irtc)
          if(irtc /= 0)then
            go to 9120
          elseif(n /= 6 .or. m /= 6)then
            call tfree(kbm)
            go to 9120
          endif
        else
          go to 9100
        endif
      endif
c      write(*,'(1p6g12.4)')(trat1(i,1:6),i=1,6)
      call tmultr(trans,trat1,irad)
      if(kbm /= 0)then
        do i=0,35
          rlist(kat1+i)=rlist(kat1+i)+rlist(krt+i)
        enddo
        call tmulbs(beam,trat1,.true.)
        do i=1,6
          do j=i,6
            beam(iaidx(i,j))=beam(iaidx(i,j))+rlist(kbm+(j-1)*6+i-1)
          enddo
        enddo
        call tmuld6(trans,rlist(krt))
        call tfree(kbm)
        call tfree(krt)
      else
        call tmulbs(beam,trat1,.true.)
      endif
      call tfree(kat1)
 9000 levele=itfdownlevel()
      isp=isp0
      return
 9120 call tfree(krt)
 9110 call tfree(kat1)
 9100 levele=itfdownlevel()
      isp=isp0
      write(lfnm,*)'ExternalMap(EMIT) of ',l,ord(l),
     $     ' element did not return ',
     $     '{cod(6), trans(6,6)} or ',
     $     '{cod(6), trans(6,6), dtrans(6,6), dbeam(6,6)}.'
      return
      end

      subroutine qemap(trans,cod,l,coup,err)
      use tfstk
      use tmacro
      use efun
      use maloc,only:tfl2m,tfmsize
      use sad_basics
      use tfcsi,only:icslfnm,lfnm
      implicit none
      type (sad_descriptor) :: kx
      type (sad_dlist), pointer :: kxl, k2l
      type (sad_rlist), pointer :: k1l
      integer*4 ,intent(in):: l
      integer*4 isp0,n,m,irtc
      real*8 ,intent(inout)::  trans(6,6),cod(6)
      character*2 ord
      logical*4 ,intent(out):: err,coup
      integer*8, save:: ifv=0,iem=0
      if(iem == 0)then
        iem=ktfsymbolz('ExternalMap',11)
        ifv=ktsalocb(0,'OPTICS',6)
      endif
      err=.true.
      levele=levele+1
      isp0=isp
      isp=isp+1
      ktastk(isp)=ktfsymbol+iem
      isp=isp+1
      ktastk(isp)=ktfstring+ifv
      isp=isp+1
      rtastk(isp)=l
      isp=isp+1
      dtastk(isp)=kxm2l(cod,0,6,1,.false.)
      kx=tfefunref(isp0+1,.false.,irtc)
      if(irtc /= 0)then
        levele=itfdownlevel()
        if(ierrorprint /= 0)then
          call tfaddmessage(' ',2,icslfnm())
        endif
        write(lfnm,*)' Error in ExternalMap(OPTICS) of ',l,ord(l),
     $       ' element.'
        isp=isp0
        return
      endif
      if(ktfnonlistq(kx,kxl) .or. kxl%nl == 3 .and.
     $     kxl%body(0) == ktfsymbol+iem)then
        go to 9200
      endif
      if(kxl%nl /= 2)then
        go to 9100
      endif
      if(tfnonreallistq(kxl%dbody(1),k1l) .or. k1l%nl /= 6)then
        go to 9100
      endif
      cod=k1l%rbody(1:6)
      call tfmsize(kxl%dbody(2),n,m,irtc)
      if(irtc /= 0 .or. n /= 6 .or. m /= 6)then
        go to 9100
      endif
      call descr_sad(kxl%dbody(2),k2l)
      call tfl2m(k2l,trans,6,6,.false.)
      coup=trans(1,3) /= 0.d0 .or. trans(1,4) /= 0.d0 .or.
     $     trans(2,3) /= 0.d0 .or. trans(2,4) /= 0.d0
      err=.false.
 9000 levele=itfdownlevel()
      isp=isp0
      return
 9100 write(lfnm,*)'ExternalMap(OPTICS) of ',l,ord(l),
     $     ' element did not return ',
     $     '{cod(6), trans(6,6)}.'
      go to 9000
 9200 call tinitr(trans)
      go to 9000
      end

      subroutine tgmap(l)
      use tfstk
      use ffs_pointer
      use tmacro
      use geolib
      use efun
      use tfcsi,only:icslfnm,lfnm
      implicit none
      type (sad_descriptor) kx
      integer*8 ktfgeol,kax,k1,k2,k11,k12,ka1,ka11,ka12,kdb
      integer*4 l,isp0,irtc
      character*2 ord
      logical*4 err
      integer*8 ifv,iem
      data ifv,iem/0,0/
      if(iem == 0)then
        iem=ktfsymbolz('ExternalMap',11)
        ifv=ktsalocb(0,'GEO',3)
      endif
      err=.true.
c      iat=itfm2l(cod,0,6,1,.false.)
      levele=levele+1
      isp0=isp
      isp=isp+1
      ktastk(isp)=ktfsymbol+iem
      isp=isp+1
      ktastk(isp)=ktfstring+ifv
      isp=isp+1
      rtastk(isp)=l
      isp=isp+1
      kdb=ktfgeol(geo(1,1,l))
      ktastk(isp)=ktflist+kdb
c      ktastk(isp)=ktflist+ktfgeol(geo(1,1,l))
      isp=isp+1
      rtastk(isp)=pos(l)
      kx=tfefunref(isp0+1,.false.,irtc)
      if(irtc /= 0)then
        levele=itfdownlevel()
        if(ierrorprint /= 0)then
          call tfaddmessage(' ',2,icslfnm())
        endif
        write(lfnm,*)' Error in ExternalMap(GEO) of ',l,ord(l),' element.'
        isp=isp0
        return
      endif
      if(ktfnonlistq(kx))then
        go to 9000
      endif
      kax=ktfaddr(kx)
      if(ilist(2,kax-1) == 4 .and. klist(kax) == ktfsymbol+iem)then
        go to 9000
      endif
      if(ilist(2,kax-1) /= 2)then
        go to 9100
      endif
      k1=klist(kax+1)
      if(ktfnonlistq(k1))then
        go to 9100
      endif
      ka1=ktfaddr(k1)
      if(ilist(2,ka1-1) /= 2)then
        go to 9100
      endif
      k11=klist(ka1+1)
      if(ktfnonlistq(k11))then
        go to 9100
      endif
      ka11=ktfaddr(k11)
      if(ilist(2,ka11-1) /= 3 .or. ktfnonreallistq(ka11))then
        go to 9100
      endif
      k12=klist(ka1+1)
      if(ktfnonlistq(k12))then
        go to 9100
      endif
      ka12=ktfaddr(k12)
      if(ilist(2,ka12-1) /= 3 .or. ktfnonreallistq(ka12))then
        go to 9100
      endif
      k2=klist(kax+2)
      if(ktfnonrealq(k2,pos(l+1)))then
        go to 9100
      endif
      call tmov(rlist(ka11+1),geo(1,4,l+1),3)
      geo(:,:,l+1)=tfchitogeo(rlist(ka12+1:ka12+3))
      levele=itfdownlevel()
      isp=isp0
      return
 9000 geo(:,:,l+1)=geo(:,:,l)
c     call tmov(geo(1,1,l),geo(1,1,l+1),12)
      pos(l+1)=pos(l)
      levele=itfdownlevel()
      isp=isp0
      return
 9100 write(lfnm,*)'ExternalMap(GEO) of ',l,ord(l),
     $     ' element did not return ',
     $     '{{{GX, GY, GX}, {CHI1, CHI2, CHI3}}, S}.'
      go to 9000
      end
