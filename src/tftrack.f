      subroutine tftrack(isp1,kx,irtc)
      use tfstk
      use ffs
      use ffs_pointer
      use tffitcode
      use tfshare
      use ffs_wake
      use tfcsi
      implicit none
      type (sad_descriptor) kx,kx1,kx2,ks,kp
      type (sad_list), pointer :: klx,kl
      integer, parameter :: nkptbl = 6
      integer*4, parameter :: npparamin=9,npnlatmin=30000
      integer*8 kz,kzp,kpz,kzf,kaxl,ktfmalocp,ktfresetparticles
      integer*4 isp1,irtc,narg,itfloc,outfl0,ld,ls,mc,npz,npa,np00,
     $     ipr(100),npr,np1,fork_worker,iprid, ne,nend,
     $     npp,ipn,m,itfmessage,nt,mt,kseed,j
      integer*8 ikptblw,ikptblm
      real*8 trf00,p00,vcphic0,vcalpha0

      logical*4 dapert0,normal
      narg=isp-isp1
      if(narg .gt. 4)then
        irtc=itfmessage(9,'General::narg','"1, 2, 3, or 4"')
        return
      endif
      nt=1
      nend=nt
      mt=1
      if(narg .ge. 3)then
        if(ktastk(isp1+3) .eq. ktfoper+mtfnull)then
        elseif(ktfrealq(ktastk(isp1+3)))then
          nt=int(rtastk(isp1+3))
          nend=nt
        else
          go to 9100
        endif
        if(narg .eq. 4)then
          if(ktastk(isp) .eq. ktfoper+mtfnull)then
          elseif(ktfrealq(ktastk(isp)))then
            nend=int(rtastk(isp))
            if(nend .lt. nt)then
              irtc=itfmessage(9,'General::wrongval',
     $             '"#4 >= #3"')
              return
            endif
            mt=nend-nt+1
          else
            go to 9100
          endif
        endif
      endif

      ld=nlat
      if(narg .ge. 2)then
        if(ktastk(isp1+2) .eq. ktfoper+mtfnull)then
        else
          ld=itfloc(ktastk(isp1+2),irtc)
          if(irtc .ne. 0)then
            return
          endif
        endif
      endif
      if(.not. tflistqd(dtastk(isp1+1),kl))then
        go to 9000
      endif
      ks=kl%dbody(1)
      ls=itfloc(ks%k,irtc)
      if(irtc .ne. 0)then
        return
      endif
      if(ls .eq. nlat)then
        ls=1
      endif
      if(ld .le. ls)then
        mt=mt+1
      endif
      kp=kl%dbody(2)
      if(.not. tflistqd(kp))then
        go to 9000
      endif
      call tfmsize(kp,mc,npz,irtc)
      if(irtc .ne. 0)then
        return
      endif
      if(.not. ((mc .eq. 7) .or. (mc .eq. 8)))then
        go to 9000
      endif
      if(npz .lt. 0)then
        go to 9000
      endif
      call tftclupdate(7)
      call tfsetparam
      wake=(twake .or. lwake) .and. trpt
      kwakep=0
      kwakeelm=0
      nwakep=0
      if(wake)then
        call tffssetupwake(icslfno(),irtc)
        if(irtc .ne. 0)then
          return
        endif
        if(nwakep .eq. 0)then
          wake=.false.
        else
          nparallel=1
        endif
      endif
      ipn=0
      npr=0
      npp=npz
      kz=ktfmalocp(kp%k,mc,npz,.false.,.true.,.true.,.true.,irtc)
      if(irtc .ne. 0)then
        go to 8900
      endif
      if(photons)then
        call tphotoninit()
      endif
      if(radlight)then
        rad=.true.
        trpt=.true.
        call tlinit(npz,h0,rlist(ifgeo+12*(ls-1)))
      endif
      if(nparallel .gt. 1)then
        kseed=0
        ne=ld-ls
        if(ne .le. 0)then
          ne=ne+nlat
        endif
        if(npz .gt. nparallel*npparamin
     $       .and. npz*ne .gt. npnlatmin)then
          irtc=1
          ikptblm=ktfallocshared(npz*((nkptbl+1)/2))
          npr=nparallel-1
          np1=npz/nparallel+1
          ipn=0
          npr=0
          do while(ipn+np1 .lt. npz)
            kseed=kseed+2
            npr=npr+1
            iprid=fork_worker()
            if(iprid .eq. 0)then
              call tfaddseed(kseed,irtc)
              if(irtc .ne. 0)then
                return
              endif
              npr=-1
              npp=np1
              go to 10
            endif
            ipr(npr)=iprid
            ipn=ipn+np1
          enddo
          npp=npz-ipn
        endif
      endif
 10   trf00=trf0
      vcphic0=vcphic
      vcalpha0=vcalpha
      if(trpt)then
        trf0=0.d0
        vcphic=0.d0
        vcalpha=1.d0
      endif
      kpz=ktaloc(npp)
      ikptblw=ktaloc(npp*((nkptbl+1)/2))
      kzp=kz+ipn
      kzf=kzp+npz*6
      p00=pgev
c      pgev=rgetgl1('MOMENTUM')
      pgev=rlist(ifgamm+ls-1)*amass
c      call tclrparaall
      call tphyzp
      call tsetdvfs
      call tfsetparticles(rlist(kzp),rlist(kzf),rlist(kpz),
     $     ilist(1,ikptblw),npp,npa,npz,mc,nlat,nt)
      if(npa .gt. 0)then
        outfl0=outfl
        outfl=0
        dapert0=dapert
        dapert=.false.
        np00=np0
        np0=npp
        if(mt .gt. 1)then
          call tturn0(npa,latt,ls,nlat,
     $         rlist(kzp),      rlist(kzp+npz),  rlist(kzp+npz*2),
     $         rlist(kzp+npz*3),rlist(kzp+npz*4),rlist(kzp+npz*5),
     $         rlist(kzp+npz*6),rlist(kpz),ilist(1,ikptblw),nt,normal)
          nt=nt+1
          mt=mt-1
          ls=1
          pgev=rlist(ifgamm)*amass
          call tphyzp
        endif
        do m=1,mt-1
          if(npr .ge. 0)then
            call tftclupdate(7)
          endif
          if(npa .le. 0)then
            nt=nt+mt-1
            mt=0
            exit
          endif
          call tturn0(npa,latt,1,nlat,
     $         rlist(kzp),      rlist(kzp+npz),  rlist(kzp+npz*2),
     $         rlist(kzp+npz*3),rlist(kzp+npz*4),rlist(kzp+npz*5),
     $         rlist(kzp+npz*6),rlist(kpz),ilist(1,ikptblw),nt,normal)
          nt=nt+1
          mt=mt-1
        enddo
        if(ld .le. ls)then
          normal=.false.
        elseif(mt .ge. 1 .and. npa .gt. 0)then
          call tturn0(npa,latt,ls,ld,
     $         rlist(kzp),      rlist(kzp+npz),  rlist(kzp+npz*2),
     $         rlist(kzp+npz*3),rlist(kzp+npz*4),rlist(kzp+npz*5),
     $         rlist(kzp+npz*6),rlist(kpz),ilist(1,ikptblw),nt,normal)
        else
          normal=.true.
        endif
        np0=np00
        outfl=outfl0
        dapert=dapert0
        rlist(kzf:kzf+npa-1)=1.d0
        rlist(kzf+npa:kzf+npp-1)=0.d0
      endif
      trf0=trf00
      vcphic=vcphic0
      vcalpha=vcalpha0
      irtc=0
      if(npr .ne. 0)then
c       - Copy particle ID to array index map[kptbl(#,1)] with ipn offset
        ilist(ipn+1:ipn+npp,ikptblm)=ilist(1:npp,ikptblw)+ipn
c       - Don't copy kptbl(#,2) becuase reversed map is not used at post process
c       - Copy iptbl(*,3:nkptbl)
        do j=3,nkptbl
          ilist( ipn+(j-1)*npz+1:ipn+(j-1)*npz+npp,ikptblm)=
     $         ilist((j-1)*npp+1:    (j-1)*npp+npp,ikptblw)
        enddo
        call tffswait(iprid,npr+1,ipr,int8(0),'tftrack',irtc)
        kaxl=ktfresetparticles(rlist(kz),ilist(1,ikptblm),npz,nlat,nend)
        call tfreeshared(ikptblm)
c        if(mapfree(iptbl(ikptblm+1)) .ne. 0)then
c          write(*,*)'???tftrack-munmap error.'
c        endif
      else
        kaxl=ktfresetparticles(rlist(kz),ilist(1,ikptblw),npz,nlat,nend)
      endif
      call tmunmapp(kz)
      if(photons)then
        call tphotonlist()
      endif
      pgev=p00
      call tphyzp
      call tfree(ikptblw)
      call tfree(kpz)
      kx=kxadaloc(-1,2,klx)
      if(normal)then
        klx%rbody(1)=dble(ld)
      else
        klx%rbody(1)=dble(ls)
      endif
      klx%body(2)=ktflist+ktfcopy1(kaxl)
      if(radlight)then
        kx1=kx
        kx=kxadaloc(-1,2,klx)
        klx%dbody(1)=dtfcopy1(kx1)
        call tlresult(kx2)
        klx%dbody(2)=dtfcopy1(kx2)
      endif
      call tclrfpe
 8900 if(wake)then
        call tffsclearwake
      endif
      return
 9000 irtc=itfmessage(9,'General::wrongtype',
     $     '"{loc,{{x1,..},{px1,..},{y1,..},{py1,..},'//
     $     '{z1,..},{dp1,..},{flg1,..}[,{lost1,..}]}} for #2"')
      return
 9100 irtc=itfmessage(9,'General::wrongval',
     $     '"[nbegin [, nend]] for #3 and #4"')
      return
      end

      subroutine tfsurvivedparticles(isp1,kx,irtc)
      use tfstk
      implicit none
      type llist
        type (sad_list), pointer :: kl
      end type
      type (llist) klxi(7), kli(7)
      type (sad_descriptor) kx
      type (sad_list), pointer :: kll,klf,klx
      integer*4 isp1,irtc,npx,i,j,itfmessage,np,ii
      if(isp .ne. isp1+1)then
        irtc=itfmessage(9,'General::narg','"1"')
        return
      endif
      if(.not. tflistqd(dtastk(isp),kll))then
        go to 9000
      endif
      if(kll%nl .ne. 7)then
        go to 9000
      endif
      if(ktfreallistqo(kll))then
        go to 9000
      endif
      if(.not. tfreallistq(kll%body(7),klf))then
        go to 9000
      endif
      np=klf%nl
      npx=0
      do i=1,np
        if(klf%rbody(i) .ne. 0.d0)then
          npx=npx+1
        endif
      enddo
      irtc=0
      if(npx .eq. np)then
        kx=dtastk(isp)
        return
      endif
      kx=kxadalocnull(-1,7,klx)
      if(npx .eq. 0)then
        do i=1,7
          klx%dbody(i)=dtfcopy1(dxnulll)
        enddo
        return
      endif
      do i=1,7
        if(.not. tfreallistqd(kll%dbody(i),kli(i)%kl))then
          go to 8900
        endif
        if(kli(i)%kl%nl .ne. np)then
          go to 8900
        endif
        klx%dbody(i)=kxavaloc(0,npx,klxi(i)%kl)
      enddo
      j=0
      do i=1,np
        if(klf%rbody(i) .ne. 0.d0)then
          j=j+1
          klxi(1)%kl%rbody(j)=kli(1)%kl%rbody(i)
          klxi(2)%kl%rbody(j)=kli(2)%kl%rbody(i)
          klxi(3)%kl%rbody(j)=kli(3)%kl%rbody(i)
          klxi(4)%kl%rbody(j)=kli(4)%kl%rbody(i)
          klxi(5)%kl%rbody(j)=kli(5)%kl%rbody(i)
          klxi(6)%kl%rbody(j)=kli(6)%kl%rbody(i)
          klxi(7)%kl%rbody(j)=1.d0
        endif
      enddo
      return
 8900 do ii=1,i-1
        klxi(ii)%kl%rbody(1:npx)=0.d0
      enddo
 9000 irtc=itfmessage(9,'General::wrongtype',
     $     '"{loc,{{x1,..},{px1,..},{y1,..},{py1,..},'//
     $     '{z1,..},{dp1,..},{flg1,..}}}"')
      return
      end
      
      integer*8 function ktfresetparticles(zx,iptbl,np,nlati,tend)
      use tfstk
      use ffs
      implicit none
      integer*8 ka,kaj(9)
      integer*4 np,iptbl(np,6),i,j,k,nlati,tend,nv
      real*8 zx(np,7)
      nv=7
      if(lossmap)then
        nv=9
      endif
      call tconvm(np,zx(1,2),zx(1,4),zx(1,6),0.d0,1)
      ka=ktadaloc(-1,nv)
      do j=1,nv
        kaj(j)=ktavaloc(0,np)
        klist(ka+j)=ktflist+kaj(j)
      enddo

      do i=1,np
        k=iptbl(i,1)
        rlist(kaj(1)+i)=zx(k,1)
        rlist(kaj(2)+i)=zx(k,2)
        rlist(kaj(3)+i)=zx(k,3)
        rlist(kaj(4)+i)=zx(k,4)
        rlist(kaj(5)+i)=zx(k,5)
        rlist(kaj(6)+i)=zx(k,6)
        if(iptbl(k,4) .eq. 0)then
          rlist(kaj(7)+i)=1.d0
        else
          rlist(kaj(7)+i)=0.d0
        endif
      enddo
      if(lossmap)then
        do i=1,np
          k=iptbl(iptbl(i,1),4)
          if((-nlati .le. k) .and. (k .le. nlati))then
            rlist(kaj(8)+i)=dble(k)
          else
            rlist(kaj(8)+i)=0.d0
          endif
          if(k .ne. 0)then
            rlist(kaj(9)+i)=dble(iptbl(iptbl(i,1),5))
          else
            rlist(kaj(9)+i)=dble(tend+1)
          endif
        enddo
      endif
      ktfresetparticles=ka
      return
      end

      subroutine tfsetparticles(zx,dv,pz,iptbl,np,npa,npc,mc,nlat,
     $     tbegin)
      implicit none
      integer, parameter :: nkptbl = 6
      real(8),    intent(inout) :: zx(npc,mc), dv(np), pz(np)
      integer(4), intent(out)   :: iptbl(np,nkptbl)
      integer(4), intent(in)    :: np
      integer(4), intent(out)   :: npa
      integer(4), intent(in)    :: npc, mc, nlat, tbegin
      integer(4) :: i, p, iptmp(nkptbl)
      real(8) :: x1, px1, y1, py1, z1, g1
c     zx(npc,mc): Initial particle distribution
c            npc: Number of particles
c             mc: Number of rows
c                 7: (x, x', y, y', z, delta, alive)
c                 8: (x, x', y, y', z, delta, alive, lost-position)
c                 9: (x, x', y, y', z, delta, alive, lost-position, lost-turns)

c     Initialize iptbl
      iptbl(1:np,1:nkptbl)=0
      do i=1,np
c       Initialize map between particle ID and array index
        iptbl(i,1)=i
        iptbl(i,2)=i
c       Initialize particle lost position[0:alive, -nlat-1:initial-dead]
        if(zx(i,7) .eq. 0.d0)then
          iptbl(i,4)=-nlat-1
          iptbl(i,5)=tbegin-1
        else
          iptbl(i,4)=0
          iptbl(i,5)=0
        endif
      enddo

c     Load particle lost position/turn from zx(#,8/9)
      if(mc .ge. 8)then
        do i=1,np
          if(iptbl(i,4) .ne. 0)then
            p=abs(nint(zx(i,8)))
            if(p .gt. 0 .and. nlat .ge. p)then
              iptbl(i,4)=-p
            endif
          endif
        enddo
      endif

c     Load particle lost turn from zx(#,9)
      if(mc .ge. 9)then
        do i=1,np
          if(iptbl(i,4) .ne. 0)then
            iptbl(i,5)=nint(zx(i,9))
          endif
        enddo
      endif

c     Compact dead particles into list-tail
      i=1
      npa=np
      do while(i .le. npa)
         if(iptbl(i,4) .ne. 0)then
c           Search avlive particle from tail: (i, npa]
            do while((i .lt. npa) .and. (iptbl(npa,4) .ne. 0))
               npa=npa-1
            enddo
            if(iptbl(npa,4) .eq. 0)then
c              Swap dead particle slot[i] with tail alive particle slot[npa]
c              - Update maps between partice ID and array index
c                Note: maps is initialized and slot swapping MUST be once
               iptbl(i,  1)=npa
               iptbl(npa,1)=i
               iptbl(npa,2)=i
               iptbl(i,  2)=npa
c              - Swap iptbl except forward/backward[iptbl(*,1)/iptbl(*,2)]
               iptmp(    3:nkptbl) = iptbl(npa,3:nkptbl)
               iptbl(npa,3:nkptbl) = iptbl(i,  3:nkptbl)
               iptbl(i,  3:nkptbl) = iptmp(    3:nkptbl)
c              - Swap particle coordinates
               x1 =zx(i,1)
               px1=zx(i,2)
               y1 =zx(i,3)
               py1=zx(i,4)
               z1 =zx(i,5)
               g1 =zx(i,6)
 
               zx(i,1)=zx(npa,1)
               zx(i,2)=zx(npa,2)
               zx(i,3)=zx(npa,3)
               zx(i,4)=zx(npa,4)
               zx(i,5)=zx(npa,5)
               zx(i,6)=zx(npa,6)
 
               zx(npa,1)=x1
               zx(npa,2)=px1
               zx(npa,3)=y1
               zx(npa,4)=py1
               zx(npa,5)=z1
               zx(npa,6)=g1
            endif
            npa=npa-1
         endif
         i=i+1
      enddo
      if(npa .gt. 0)then
        call tconvm(npa,zx(1,2),zx(1,4),zx(1,6),dv,-1)
      endif
      return
      end

      subroutine tfwriteunformatted(lfn,a,n)
      implicit none
      integer*4 lfn,n
      real*8 a(n)
      write(lfn)a
      return
      end

      subroutine tfreadunformatted(lfn,a,n,irtc)
      implicit none
      integer*4 lfn,n,irtc
      real*8 a(n)
      read(lfn,ERR=900)a
      irtc=0
      return
 900  irtc=997
      return
      end

      subroutine tfaddseed(kseed,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) kx,kn
      type (sad_list), pointer :: klx
      integer*4 kseed,isp1,irtc
      real*8 vn
      isp1=isp
      isp=isp1+1
      dtastk(isp)=dxnull
      call SeedRandom(isp1,kx,irtc)
      if(irtc .ne. 0)then
        return
      elseif(ktfnonlistqd(kx,klx))then
        irtc=-1
        return
      endif
      if(klx%nl .lt. 2)then
        irtc=-1
        return
      endif
c      call tfdebugprint(kx,'tfaddseed',1)
      kn=klx%dbody(2)
      if(ktfnonrealqd(kn,vn))then
        irtc=-1
        return
      endif
      isp=isp1+1
      rtastk(isp)=vn+kseed
      call SeedRandom(isp1,kx,irtc)
c      call tfdebugprint(kx,'=> ',1)
      isp=isp1
      irtc=0
      return
      end
