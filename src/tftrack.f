      module tftr
      use tfstk
      integer*4 ,parameter :: nprmax=256, nkptbl = 6, nkptbl1=3,
     $     npparamin=9,npnlatmin=3000
      integer*4 :: npz=0,ipn=0,npr=0,npri=0,iprid=0,ipr(nprmax)

      contains
      subroutine tftrack(isp1,kx,irtc)
      use ffs
      use ffs_pointer
      use tffitcode
      use tfshare
      use ffs_wake
      use tparastat
      use photontable,only:tphotoninit,tphotonlist
      use trexc
      use tfcsi
      use maloc,only:ktfmalocp,tfmsize
      use iso_c_binding
      implicit none
      integer*4 ,intent(in):: isp1
      integer*4 ,intent(out):: irtc
      type (sad_descriptor) ,intent(out):: kx
      type (sad_descriptor) kx1,kx2,ks,kp
      type (sad_dlist), pointer :: klx,kl
      integer*8 kz,kzp,kzf,kaxl,kdv,kpsx,kpsy,kpsz
      integer*4 narg,itfloc,outfl0,ld,ls,mc,npa,np00,
     $     np1,ne,nend,npara,npp,m,itfmessage,nt,mt,kseed,mcf
      integer*8 ikptblw,ikptblm
      real*8 trf00,p00,vcalpha0
      real*8 , pointer::zx(:,:),zx0(:,:)
      integer*4 , pointer::iptbl(:,:),jptbl(:,:)
      logical*4 :: dapert0,normal=.true.
      character*32 autos
      narg=isp-isp1
      if(narg > 4)then
        irtc=itfmessage(9,'General::narg','"1, 2, 3, or 4"')
        return
      endif
      nt=1
      nend=1
      mt=1
      if(narg >= 3)then
        if(ktastk(isp1+3) == ktfoper+mtfnull)then
        elseif(ktfrealq(ktastk(isp1+3),nt))then
          nend=nt
        else
          go to 9100
        endif
        if(narg == 4)then
          if(ktastk(isp) == ktfoper+mtfnull)then
          elseif(ktfrealq(ktastk(isp),nend))then
            if(nend < nt)then
              irtc=itfmessage(9,'General::wrongval','"#4 >= #3"')
              return
            endif
            mt=nend-nt+1
          else
            go to 9100
          endif
        endif
      endif

      ld=nlat
      if(narg >= 2)then
        if(ktastk(isp1+2) == ktfoper+mtfnull)then
        else
          ld=itfloc(ktastk(isp1+2),irtc)
          if(irtc /= 0)then
            return
          endif
        endif
      endif
      if(.not. tflistq(dtastk(isp1+1),kl))then
        go to 9000
      endif
      ks=kl%dbody(1)
      ls=itfloc(ks%k,irtc)
      if(irtc /= 0)then
        return
      endif
      if(ls > 0)then
        ls=mod(ls-1,nlat-1)+1
      else
        ls=nlat-mod(1-ls,nlat-1)
      endif
      if(ld > 1)then
        ld=mod(ld-2,nlat-1)+2
      else
        ld=nlat+1-mod(2-ld,nlat-1)
      endif
      if(ld <= ls)then
        mt=mt+1
      endif
      kp=kl%dbody(2)
      if(.not. tflistq(kp))then
        go to 9000
      endif
      call tfmsize(kp,mc,npz,irtc)
      if(irtc /= 0)then
        return
      endif
      if(npz < 0)then
        go to 9000
      endif
      if(calpol)then
        if(.not. (mc >= 9 .and. mc <= 11))then
          go to 9000
        endif
        mcf=9
      else
        if(.not. (mc >= 7 .and. mc <= 9))then
          go to 9000
        endif
        mcf=7
      endif
      call tfsetupex
      call tftclupdate(7)
      call tfsetparam
      wake=(twake .or. lwake) .and. trpt
      kwakep=0
      kwakeelm=0
      nwakep=0
      npara=min(nprmax,max(nparallel,1))
      if(wake)then
        call tffssetupwake(icslfnm(),irtc)
        if(irtc /= 0)then
          return
        endif
        if(nwakep == 0)then
          wake=.false.
        else
          npara=1
        endif
      endif
      ipn=0
      npr=0
      npp=npz
      kz=ktfmalocp(kp,mc,npz,.false.,.true.,.true.,.true.,irtc)
      if(irtc /= 0)then
        go to 8900
      endif
      call c_f_pointer(c_loc(rlist(kz)),zx0,[npz,mc])
      if(photons)then
        npara=1
        if(rad .and. rfluct)then
          call tphotoninit()
        else
          irtc=itfmessage(9,'Track::photons','""')
          return
        endif
      endif
      if(radlight)then
        rad=.true.
        trpt=.true.
        call tlinit(npz,h0,rlist(ifgeo+12*(ls-1)))
      endif
      call tfevals('`ExtMap$@InitMap['//autos(dble(npz))//']',kx,irtc)
      if(npara > 1)then
        kseed=0
        ne=ld-ls
        if(ne <= 0)then
          ne=ne+nlat
        endif
        npara=min(npara,npz/npparamin+1,ne*npz/npnlatmin+1)
        if(npara > 1)then
          irtc=1
          ikptblm=ktfallocshared(npz*nkptbl1)
          npr=npara-1
          np1=npz/npara
          ipn=0
          npr=0
          do while(ipn+np1 < npz .and. npr < npara)
            kseed=kseed+2
            npri=npr
            npr=npr+1
            iprid=itffork()
c            write(*,'(a,6i7)')'fork ',iprid,ipn
            if(iprid == 0)then
              call tfaddseed(kseed,irtc)
              if(irtc /= 0)then
                write(*,*)'addseed-error ',irtc
                call exit_without_hooks(0)
              endif
              npr=-1
              npp=np1
              call tsetintm(-1.d0)
              go to 10
            endif
            ipr(npr)=iprid
            ipn=ipn+np1
          enddo
          npri=npr
          npp=npz-ipn
        endif
      endif
 10   trf00=trf0
      vcalpha0=vcalpha
      if(trpt)then
        trf0=0.d0
        vcalpha=1.d0
      endif
      kdv=ktaloc(npz)
      ikptblw=ktaloc(npp*nkptbl1)
      kzp=kz+ipn
      kzf=kzp+npz*(mcf-1)
      if(calpol)then
        kpsx=kzp+npz*6
        kpsy=kzp+npz*7
        kpsz=kzp+npz*8
      else
        kpsx=kdv
        kpsy=kdv
        kpsz=kdv
      endif
      p00=pgev
      pgev=rlist(ifgamm+ls-1)*amass
      call tphyzp
      call tsetdvfs
      call c_f_pointer(c_loc(rlist(kzp)),zx,[npz,mc])
      call c_f_pointer(c_loc(ilist(1,ikptblw)),iptbl,[npp,nkptbl])
      call c_f_pointer(c_loc(ilist(1,ikptblm)),jptbl,[npz,nkptbl])
      call tfsetparticles(zx,rlist(kdv:kdv+npp-1),
     $     iptbl,npp,npa,npz,mc,nlat,nt,mcf)
c      jptbl(ipn+1:ipn+npp,1)=iptbl(1:npp,1)+ipn
c      write(*,'(a,6i7,1p8g15.7)')'after-setparticles',iprid,npz,npa,ipn,m,jptbl(398,1)
c     $     zx0(iptbl(398,1),1),zx0(iptbl(398,1),7)
      if(npa > 0)then
        outfl0=outfl
        outfl=0
        dapert0=dapert
        dapert=.false.
        np00=np0
        np0=npp
c        jptbl(ipn+1:ipn+npp,1)=iptbl(1:npp,1)+ipn
c        do m=1,npp
c          if(iptbl(m,1)+ipn == jptbl(398,1))then
c            write(*,'(a,6i7,1p8g15.7)')'tftrack-setparticles',iprid,npz,npa,ipn,m,jptbl(398,1),
c     $           zx0(jptbl(398,1),1),zx0(jptbl(398,1),7)
c            exit
c          endif
c        enddo
        if(mt > 1)then
          call tturn(npa,ls,nlat,
     $         zx(1:npa,1),zx(1:npa,2),zx(1:npa,3),zx(1:npa,4),
     $         zx(1:npa,5),zx(1:npa,6),
     $         rlist(kdv),rlist(kpsx),rlist(kpsy),rlist(kpsz),
     $         iptbl,nt,normal)
          nt=nt+1
          mt=mt-1
          ls=1
          pgev=rlist(ifgamm)*amass
          call tphyzp
        endif
        do m=1,mt-1
          if(npr >= 0)then
            call tftclupdate(7)
          endif
          if(npa <= 0)then
            nt=nt+mt-1
            mt=0
            exit
          endif
          call tturn(npa,1,nlat,
     $         zx(1:npa,1),zx(1:npa,2),zx(1:npa,3),zx(1:npa,4),
     $         zx(1:npa,5),zx(1:npa,6),
     $         rlist(kdv),rlist(kpsx),rlist(kpsy),rlist(kpsz),
     $         iptbl,nt,normal)
          nt=nt+1
          mt=mt-1
        enddo
        normal=.true.
        if(ld <= ls)then
          normal=.false.
        elseif(mt >= 1 .and. npa > 0)then
          call tturn(npa,ls,ld,
     $         zx(1:npa,1),zx(1:npa,2),zx(1:npa,3),zx(1:npa,4),
     $         zx(1:npa,5),zx(1:npa,6),
     $         rlist(kdv),rlist(kpsx),rlist(kpsy),rlist(kpsz),
     $         iptbl,nt,normal)
        endif
        np0=np00
        outfl=outfl0
        dapert=dapert0
      endif
      trf0=trf00
      vcalpha=vcalpha0
      irtc=0
      if(npr /= 0)then
c       - Copy particle ID to array index map[kptbl(#,1)] with ipn offset
        jptbl(ipn+1:ipn+npp,1)=iptbl(1:npp,1)+ipn
c        ilist(ipn+1:ipn+npp,ikptblm)=ilist(1:npp,ikptblw)+ipn
c       - Don't copy kptbl(#,2) becuase reversed map is not used at post process
c       - Copy iptbl(*,3:nkptbl)
        jptbl(ipn+1:ipn+npp,3:nkptbl)=iptbl(1:npp,3:nkptbl)
c        do m=1,npp
c          if(iptbl(m,1)+ipn == jptbl(398,1))then
c            write(*,'(a,6i7,1p8g15.7)')'tftrack-wait      ',iprid,npz,npa,ipn,m,jptbl(398,1),
c     $           zx0(jptbl(398,1),1),zx0(jptbl(398,1),7)
c            exit
c          endif
c        enddo
        call tffswait(iprid,npr+1,ipr,i00,100000,'tftrack',irtc)
c        write(*,'(a,6i7,1p8g15.7)')'tftrack-afterwait ',iprid,npz,npa,ipn,0,jptbl(398,1),
c     $         zx0(jptbl(398,1),1),zx0(jptbl(398,1),7)
        kaxl=ktfresetparticles(zx0,jptbl,npz,nlat,nend,mc)
        call tfreeshared(ikptblm)
      else
        kaxl=ktfresetparticles(zx0,iptbl,npz,nlat,nend,mc)
      endif
      call tmunmapp(kz)
      if(.not. ktfrealq(kx))then
        call tfevals('`ExtMap$@ResetMap[]',kx,irtc)
      endif
      if(photons)then
        call tphotonlist()
      endif
      pgev=p00
      call tphyzp
      call tfree(ikptblw)
      call tfree(kdv)
c      if(.not. calpol)then
c        call tfree(kpsx)
c        call tfree(kpsy)
c        call tfree(kpsz)
c      endif
      kx=kxadaloc(-1,2,klx)
      klx%rbody(1)=merge(dble(ld),dble(ls),normal)
      klx%dbody(2)%k=ktflist+ktfcopy1(kaxl)
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
     $     '{z1,..},{dp1,..}[,{sz1...},{psis1...}]'//
     $     ',{flg1,..}[,{lost1,..}]}} for #2"')
      return
 9100 irtc=itfmessage(9,'General::wrongval',
     $     '"[nbegin [, nend]] for #3 and #4"')
      return
c 9200 irtc=itfmessage(9,'FFS:nospin','""')
c      return
      end

      subroutine tfsurvivedparticles(isp1,kx,irtc)
      use ffs_flag, only:calpol
      implicit none
      type llist
        type (sad_rlist), pointer :: kl
      end type
      type (llist) klxi(10), kli(10)
      type (sad_descriptor) kx
      type (sad_dlist), pointer :: kll,klx
      type (sad_rlist), pointer :: klf
      integer*4 isp1,irtc,npx,i,j,k,itfmessage,np,ii,kcf
      if(isp /= isp1+1)then
        irtc=itfmessage(9,'General::narg','"1"')
        return
      endif
      if(.not. tflistq(dtastk(isp),kll))then
        go to 9000
      endif
      if(calpol)then
        if(kll%nl /= 9 .and. kll%nl /= 11)then
          go to 9000
        endif
        kcf=9
      else
        if(kll%nl /= 7 .and. kll%nl /= 9)then
          go to 9000
        endif
        kcf=7
      endif
      if(ktfreallistq(kll))then
        go to 9000
      endif
      if(.not. tfreallistq(kll%dbody(kcf),klf))then
        go to 9000
      endif
      np=klf%nl
      npx=0
      do i=1,np
        if(klf%rbody(i) /= 0.d0)then
          npx=npx+1
        endif
      enddo
      irtc=0
      if(npx == np)then
        kx=dtastk(isp)
        return
      endif
      kx=kxadalocnull(-1,kll%nl,klx)
      if(npx == 0)then
        do i=1,kll%nl
          klx%dbody(i)=dtfcopy1(dxnulll)
        enddo
        return
      endif
      do i=1,kll%nl
        if(.not. tfreallistq(kll%dbody(i),kli(i)%kl))then
          go to 8900
        endif
        if(kli(i)%kl%nl /= np)then
          go to 8900
        endif
        klx%dbody(i)=kxavaloc(0,npx,klxi(i)%kl)
      enddo
      j=0
      do i=1,np
        if(klf%rbody(i) /= 0.d0)then
          j=j+1
          do k=1,kll%nl
            klxi(k)%kl%rbody(j)=kli(k)%kl%rbody(i)
          enddo
          klxi(kcf)%kl%rbody(j)=1.d0
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
      
      integer*8 function ktfresetparticles(zx,iptbl,np,nlati,tend,mc)
      use ffs
      use tspin
      use ffs_flag, only:calpol
      use trexc
      implicit none
      integer*4 ,intent(in):: np,iptbl(np,6),nlati,mc,tend
      real*8 ,intent(inout):: zx(np,mc)
      integer*8 ka,kaj(11)
      integer*4 i,j,k,mcf,nv
      mcf=merge(9,7,calpol)
      nv=merge(mcf+2,mcf,lossmap)
      ka=ktadaloc(-1,nv)
      do j=1,nv
        kaj(j)=ktavaloc(0,np)
        klist(ka+j)=ktflist+kaj(j)
      enddo
      call tconvm(np,zx(:,2),zx(:,4),zx(:,6))
      do i=1,np
        k=iptbl(i,1)
        rlist(kaj(1:6)+i)=zx(k,1:6)
        rlist(kaj(mcf)+i)=merge(1.d0,0.d0,iptbl(k,4) == 0)
      enddo
      if(calpol)then
        do i=1,np
          k=iptbl(i,1)
          rlist(kaj(8)+i)=zx(k,8)/norm2(zx(k,7:9))
          rlist(kaj(7)+i)=atan(zx(k,9),zx(k,7))
        enddo
      endif
      if(lossmap)then
        do i=1,np
          k=iptbl(iptbl(i,1),4)
          rlist(kaj(mcf+1)+i)=merge(dble(k),0.d0,
     $         (-nlati <= k) .and. (k <= nlati))
          rlist(kaj(mcf+2)+i)=merge(dble(iptbl(iptbl(i,1),5)),
     $         dble(tend+1),k /= 0)
        enddo
      endif
      ktfresetparticles=ka
      return
      end

      subroutine tfsetparticles(zx,dv,iptbl,np,npa,npc,mc,nlat,tbegin,mcf)
      use tspin
      use ffs_flag,only:calpol
      use trexc
      use mathfun, only: sqrt1
      implicit none
      real(8),    intent(inout) :: zx(npc,mc), dv(np)
      integer(4), intent(out)   :: iptbl(np,nkptbl)
      integer(4), intent(in)    :: np
      integer(4), intent(out)   :: npa
      integer(4), intent(in)    :: npc, mc, nlat, tbegin, mcf
      integer(4) :: i, p, iptmp(nkptbl)
      real(8) :: x1, px1, y1, py1, z1, g1, st, sv(2)
c     zx(npc,mc): Initial particle distribution
c            npc: Number of particles
c             mc: Number of rows
c             mcf: Row for flag
c             if(nocalpol)
c                 7: (x, x', y, y', z, delta, alive)
c                 8: (x, x', y, y', z, delta, alive, lost-position)
c                 9: (x, x', y, y', z, delta, alive, lost-position, lost-turns)
c             if(calpol)
c                 9:  (x, x', y, y', z, delta, phis, sz, alive)
c                 10: (x, x', y, y', z, delta, phis, sz, alive, lost-position)
c                 11: (x, x', y, y', z, delta, phis, sz, alive, lost-position, lost-turns)

c     Initialize iptbl
c      iptbl(1:np,1:nkptbl)=0
      iptbl=0
      do i=1,np
c       Initialize map between particle ID and array index
        iptbl(i,1)=i
        iptbl(i,2)=i
c       Initialize particle lost position[0:alive, -nlat-1:initial-dead]
        if(zx(i,mcf) == 0.d0)then
          iptbl(i,4)=-nlat-1
          iptbl(i,5)=tbegin-1
        else
          iptbl(i,4)=0
          iptbl(i,5)=0
        endif
      enddo

c     Load particle lost position/turn from zx(#,8/9)
      if(mc >= mcf+1)then
        do i=1,np
          if(iptbl(i,4) /= 0)then
            p=abs(nint(zx(i,mcf+1)))
            if(p > 0 .and. nlat >= p)then
              iptbl(i,4)=-p
            endif
          endif
        enddo
      endif

c     Load particle lost turn from zx(#,9)
      if(mc >= mcf+2)then
        do i=1,np
          if(iptbl(i,4) /= 0)then
            iptbl(i,5)=nint(zx(i,mcf+2))
          endif
        enddo
      endif

c     Compact dead particles into list-tail
      i=1
      npa=np
      do while(i <= npa)
         if(iptbl(i,4) /= 0)then
c           Search avlive particle from tail: (i, npa]
            do while((i < npa) .and. (iptbl(npa,4) /= 0))
               npa=npa-1
            enddo
            if(iptbl(npa,4) == 0)then
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
               
               if(calpol)then
                 sv=zx(i,7:8)
                 zx(i,7:8)=zx(npa,7:8)
                 zx(npa,7:8)=sv
               endif
            endif
            npa=npa-1
         endif
         i=i+1
      enddo
      if(npa > 0)then
        call tconvm(npa,zx(:,2),zx(:,4),zx(:,6),dv)
        if(calpol)then
          do i=1,npa
            st=1.d0+sqrt1(-zx(i,8)**2)
            zx(i,9)=st*sin(zx(i,7))
            zx(i,7)=st*cos(zx(i,7))
          enddo
        endif
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
      implicit none
      integer*4 ,intent(in):: kseed
      integer*4 ,intent(out):: irtc
      type (sad_descriptor) kx,kn
      type (sad_dlist), pointer :: klx
      integer*4 isp1
      real*8 vn
      isp1=isp
      isp=isp1+1
      ktastk(isp)=ktfoper+mtfnull
      call SeedRandom(isp1,kx,irtc)
      if(irtc /= 0)then
        return
      elseif(ktfnonlistq(kx,klx))then
        irtc=-1
        return
      endif
      if(klx%nl < 2)then
        irtc=-1
        return
      endif
      kn=klx%dbody(2)
      if(ktfnonrealq(kn,vn))then
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

      end module
