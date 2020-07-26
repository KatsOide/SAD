      subroutine temap(np,np0,x,px,y,py,z,g,dv,l,nt,kptbl)
      use tfstk
      use efun
      use temw, only:tmulbs
      implicit none
      type alist
        type (sad_rlist), pointer :: p
      end type
      type (alist), dimension(7) :: kav
      type (sad_descriptor) kx
      type (sad_dlist), pointer :: klx,kl
      type (sad_rlist), pointer :: klrk
      integer, parameter :: nkptbl = 6
      integer*4 l,nt,np,np0,kptbl(np0,nkptbl),itfdownlevel
      real*8 x(np0),px(np0),y(np0),py(np0),z(np0),g(np0),dv(np0)
      character*2 ord
      real*8 xa(7), pr
      integer*4 i,j,k,isp0,irtc,m
      integer*8,save:: iem=0,ifv
      logical*4 dodrop,doinject
      integer*4 kptmp(nkptbl)
      if(itfcontext .le. 0)then
        return
      endif
      if(iem .eq. 0)then
        iem=ktfsymbolz('ExternalMap',11)
        ifv=ktsalocb(0,'TRACK',5)
      endif
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
      dtastk(isp)=kxadaloc(-1,7,kl)
      do k=1,7
        kl%dbody(k)=kxavaloc(0,np0,kav(k)%p)
      enddo
      do i=1,np0
        j=kptbl(i,1)
        if((j .gt. np) .or. (kptbl(j,4) .ne. 0))then
          kav(1)%p%rbody(i)=0.d0
          kav(2)%p%rbody(i)=0.d0
          kav(3)%p%rbody(i)=0.d0
          kav(4)%p%rbody(i)=0.d0
          kav(5)%p%rbody(i)=0.d0
          kav(6)%p%rbody(i)=0.d0
          kav(7)%p%rbody(i)=0.d0
        else
          pr=g(j)
          kav(1)%p%rbody(i)=x(j)
          kav(2)%p%rbody(i)=px(j)*(1.d0+pr)
          kav(3)%p%rbody(i)=y(j)
          kav(4)%p%rbody(i)=py(j)*(1.d0+pr)
          kav(5)%p%rbody(i)=z(j)
          kav(6)%p%rbody(i)=pr
          kav(7)%p%rbody(i)=1.d0
        endif
      enddo
      kx=tfefunref(isp0+1,.false.,irtc)
      if(irtc .ne. 0)then
        levele=itfdownlevel()
        isp=isp0
        if(ierrorprint .ne. 0)then
          call tfaddmessage(' ',2,6)
        endif
        write(*,*)' Error in ExternalMap of ',l,ord(l),' element at ',
     $       nt,ord(nt),' turn.'
        return
      elseif(ktflistq(kx,klx))then
        if(klx%head%k .ne. ktfoper+mtflist .or.
     $       klx%nl .ne. 7 .or. ktfreallistq(klx))then
          go to 9000
        endif
        do k=1,7
          if(.not. tfreallistq(klx%dbody(k),klrk))then
            go to 9000
          endif
          if(klrk%nl .ne. np0)then
            go to 9000
          endif
          kav(k)%p=>klrk
        enddo
        dodrop=.false.
        doinject=.false.
        do i=1,np0
          j=kptbl(i,1)
          if(kav(7)%p%rbody(i) .ne. 0.d0)then
            if(.not. (j .le. np)) then
c           Case: dropped before MAP element
              doinject=.true.
            endif
c           Re-activate particle slot[j] if kptbl(j,4) != 0
c           Note: kptbl(j,4) MUST be `0' for alive particles
            kptbl(j,4)=0
            kptbl(j,5)=nt+1
c           Copy-in ExternalMap[] result for alive/reinject case
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
          elseif((j .le. np) .and. (kptbl(j,4) .eq. 0))then
c           Lose particle slot[j] at current MAP element[l]
            dodrop=.true.
            kptbl(j,4)=l
            kptbl(j,5)=nt
          endif              
        enddo

        if(doinject)then
          i=np+1
          m=np0
          do while(i .le. m)
             if(kptbl(i,4) .ne. 0)then
c               Search alive paricle from tail: (i, m]
                do while((i .lt. m) .and. (kptbl(m,4) .ne. 0))
                   m=m-1
                enddo
                if(kptbl(m,4) .eq. 0)then
c                  Swap dead particle slot[i] with tail alive particle slot[m]
                   j=kptbl(m,2)
                   k=kptbl(i,2)
c                  - Update maps between particle ID and array index
                   kptbl(k,1)=m
                   kptbl(j,1)=i
                   kptbl(m,2)=k
                   kptbl(i,2)=j
c                  - Swap kptbl except forward/backward[kptbl(*,1)/kptbl(*,2)]
                   kptmp(  3:nkptbl) = kptbl(m,3:nkptbl)
                   kptbl(m,3:nkptbl) = kptbl(i,3:nkptbl)
                   kptbl(i,3:nkptbl) = kptmp(  3:nkptbl)
c                  - Swap particle coordinates
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
                endif
                m=m-1
             endif
             i=i+1
          enddo
          np=m
        endif

        if(dodrop)then
        endif
      else
c        call tfdebugprint(kx,'temap',1)
      endif
 9000 levele=itfdownlevel()
      isp=isp0
      return
      end

      subroutine temape(trans,cod,beam,l)
      use tfstk
      use tmacro
      use sad_main, only:iaidx
      use temw,only:tmulbs
      use efun
      use iso_c_binding
      implicit none
      type (sad_descriptor) kx
      integer*8 k1,k2,k3,k4,kax,
     $     ktfmalocp,ka1,kat1,kbm,krt
      integer*4 l,isp0,itfdownlevel,n,m,irtc,i,j
      real*8 trans(6,6),cod(6),beam(42)
      real*8 ,pointer::trat1(:,:)
      character*2 ord
      integer*8 , save :: ifv=0,iem=0
c      iaidx(m,n)=((m+n+abs(m-n))**2+2*(m+n)-6*abs(m-n))/8
      if(iem .eq. 0)then
        iem=ktfsymbolz('ExternalMap',11)
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
      if(irtc .ne. 0)then
        levele=itfdownlevel()
        if(ierrorprint .ne. 0)then
          call tfaddmessage(' ',2,6)
        endif
        write(*,*)' Error in ExternalMap(EMIT) of ',l,ord(l),
     $       ' element.'
        isp=isp0
        return
      endif
      if(ktfnonlistq(kx))then
        go to 9000
      endif
      kax=ktfaddr(kx)
      if(ilist(2,kax-1) .eq. 3 .and.
     $     klist(kax) .eq. ktfsymbol+iem)then
        go to 9000
      endif
      k1=klist(kax+1)
      if(ktfnonlistq(k1))then
        go to 9100
      endif
      ka1=ktfaddr(k1)
      if(ilist(2,ka1-1) .ne. 6 .or. ktfnonreallistq(ka1))then
        go to 9100
      endif
      cod=rlist(ka1+1:ka1+6)
      k2=klist(kax+2)
      kat1=ktfmalocp(k2,n,m,.false.,.false.,
     $     .false.,.false.,irtc)
      if(irtc .ne. 0)then
        go to 9100
      endif
      if(n .ne. 6 .or. m .ne. 6)then
        call tfree(kat1)
        go to 9100
      endif
      call c_f_pointer(c_loc(rlist(kat1)),trat1,[6,6])
      kbm=0
      krt=0
      if(ilist(2,kax-1) .ne. 2)then
        if(ilist(2,kax-1) .eq. 4)then
          k3=klist(kax+3)
          krt=ktfmalocp(k3,n,m,.false.,.false.,
     $         .false.,.false.,irtc)
          if(irtc .ne. 0)then
            go to 9110
          elseif(n .ne. 6 .or. m .ne. 6)then
            go to 9120
          endif
          k4=klist(kax+4)
          kbm=ktfmalocp(k4,n,m,.false.,.false.,
     $         .false.,.false.,irtc)
          if(irtc .ne. 0)then
            go to 9120
          elseif(n .ne. 6 .or. m .ne. 6)then
            call tfree(kbm)
            go to 9120
          endif
        else
          go to 9100
        endif
      endif
      call tmultr(trans,trat1,irad)
      if(kbm .ne. 0)then
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
      write(*,*)'ExternalMap(EMIT) of ',l,ord(l),
     $     ' element did not return ',
     $     '{cod(6), trans(6,6)} or ',
     $     '{cod(6), trans(6,6), dtrans(6,6), dbeam(6,6)}.'
      return
      end

      subroutine qemap(trans,cod,l,coup,err)
      use tfstk
      use tmacro
      use efun
      implicit none
      type (sad_descriptor) :: kx
      type (sad_dlist), pointer :: kxl, k2l
      type (sad_rlist), pointer :: k1l
      integer*4 l,isp0,itfdownlevel,n,m,irtc
      real*8 trans(6,6),cod(6)
      character*2 ord
      logical*4 err,coup
      integer*8, save:: ifv=0,iem=0
      if(iem .eq. 0)then
        iem=ktfsymbolz('ExternalMap',11)
        ifv=ktsalocb(0,'OPTICS',6)
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
      dtastk(isp)=kxm2l(cod,0,6,1,.false.)
c      itastk(2,isp)=iat
      kx=tfefunref(isp0+1,.false.,irtc)
      if(irtc .ne. 0)then
        levele=itfdownlevel()
        if(ierrorprint .ne. 0)then
          call tfaddmessage(' ',2,6)
        endif
        write(*,*)' Error in ExternalMap(OPTICS) of ',l,ord(l),
     $       ' element.'
        isp=isp0
        return
      endif
      if(ktfnonlistq(kx,kxl))then
        go to 9200
      endif
      if(kxl%nl .eq. 3 .and.
     $     kxl%body(0) .eq. ktfsymbol+iem)then
        go to 9200
      endif
      if(kxl%nl .ne. 2)then
        go to 9100
      endif
      if(tfnonreallistq(kxl%dbody(1),k1l) .or. k1l%nl .ne. 6)then
        go to 9100
      endif
      cod=k1l%rbody(1:6)
      call tfmsize(kxl%dbody(2),n,m,irtc)
      if(irtc .ne. 0 .or. n .ne. 6 .or. m .ne. 6)then
        go to 9100
      endif
      call descr_sad(kxl%dbody(2),k2l)
      call tfl2m(k2l,trans,6,6,.false.)
      coup=trans(1,3) .ne. 0.d0 .or. trans(1,4) .ne. 0.d0 .or.
     $     trans(2,3) .ne. 0.d0 .or. trans(2,4) .ne. 0.d0
      err=.false.
 9000 levele=itfdownlevel()
      isp=isp0
      return
 9100 write(*,*)'ExternalMap(OPTICS) of ',l,ord(l),
     $     ' element did not return ',
     $     '{cod(6), trans(6,6)}.'
      go to 9000
 9200 call tinitr(trans,6)
      go to 9000
      end

      subroutine tgmap(l)
      use tfstk
      use ffs_pointer
      use tmacro
      use geolib
      use efun
      implicit none
      type (sad_descriptor) kx
      integer*8 ktfgeol,kax,k1,k2,k11,k12,ka1,ka11,ka12,kdb
      integer*4 l,isp0,irtc,itfdownlevel
      real*8 rfromk
      character*2 ord
      logical*4 err
      integer*8 ifv,iem
      data ifv,iem/0,0/
      if(iem .eq. 0)then
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
      if(irtc .ne. 0)then
        levele=itfdownlevel()
        if(ierrorprint .ne. 0)then
          call tfaddmessage(' ',2,6)
        endif
        write(*,*)' Error in ExternalMap(GEO) of ',l,ord(l),
     $       ' element.'
        isp=isp0
        return
      endif
      if(ktfnonlistq(kx))then
        go to 9000
      endif
      kax=ktfaddr(kx)
      if(ilist(2,kax-1) .eq. 4 .and.
     $     klist(kax) .eq. ktfsymbol+iem)then
        go to 9000
      endif
      if(ilist(2,kax-1) .ne. 2)then
        go to 9100
      endif
      k1=klist(kax+1)
      if(ktfnonlistq(k1))then
        go to 9100
      endif
      ka1=ktfaddr(k1)
      if(ilist(2,ka1-1) .ne. 2)then
        go to 9100
      endif
      k11=klist(ka1+1)
      if(ktfnonlistq(k11))then
        go to 9100
      endif
      ka11=ktfaddr(k11)
      if(ilist(2,ka11-1) .ne. 3 .or. ktfnonreallistq(ka11))then
        go to 9100
      endif
      k12=klist(ka1+1)
      if(ktfnonlistq(k12))then
        go to 9100
      endif
      ka12=ktfaddr(k12)
      if(ilist(2,ka12-1) .ne. 3 .or. ktfnonreallistq(ka12))then
        go to 9100
      endif
      k2=klist(kax+2)
      if(ktfnonrealq(k2))then
        go to 9100
      endif
      call tmov(rlist(ka11+1),geo(1,4,l+1),3)
      geo(:,:,l+1)=tfchitogeo(rlist(ka12+1:ka12+3))
      pos(l+1)=rfromk(k2)
      levele=itfdownlevel()
      isp=isp0
      return
 9000 geo(:,:,l+1)=geo(:,:,l)
c     call tmov(geo(1,1,l),geo(1,1,l+1),12)
      pos(l+1)=pos(l)
      levele=itfdownlevel()
      isp=isp0
      return
 9100 write(*,*)'ExternalMap(GEO) of ',l,ord(l),
     $     ' element did not return ',
     $     '{{{GX, GY, GX}, {CHI1, CHI2, CHI3}}, S}.'
      go to 9000
      end
