      module dfun

      contains
      subroutine tdfun(iqcol,lfp,nqcola,nqcola1,kdp,df1,error)
      use tfstk
      use ffs
      use ffs_pointer
      use ffs_fit
      use tffitcode
      use gfun
      implicit none
      real*8 ,parameter :: abmax=1.d-8
      real*8 ,parameter::factor=0.97d0,dmax=1.d10
      integer*4 ,parameter::npeak=10
      integer*4 ,intent(out):: nqcola,nqcola1,
     $     iqcol(*),lfp(2,maxcond),kdp(*)
      real*8 ,intent(out):: df1(*)
      real*8 vpeak(npeak),vf,vb,ve,v,vf1,ve0
      logical*4 ,intent(out):: error
      integer*4 ipeak(npeak),j,i,ka,kf,kp,kp1,mp,idp,kpb,kpe,
     $     k,ip,irtc,m
      logical*4 maxfit,ttrans(-nfam:nfam),tftype1fit
c      do j=1,nfcol
c        if(flv%kfit(flv%kfitp(j)) == mfitnx .or.
c     $       flv%kfit(flv%kfitp(j)) == mfitny)then
c          do m=nfam1,nfam
c            detr0=utwiss(mfitr1,m,1)*utwiss(mfitr4,m,1)
c     $           -utwiss(mfitr2,m,1)*utwiss(mfitr3,m,1)
c            ttrans(m)=.false.
c            if(detr0 < 1.d0)then
c              do i=2,nut
c                if(utwiss(mfitr1,m,i)*utwiss(mfitr4,m,i)
c     $               -utwiss(mfitr2,m,i)*utwiss(mfitr3,m,i)
c     1               >= 1.d0)then
c                  ttrans(m)=.true.
c                  exit
c                endif
c              enddo
c            else
c              do i=2,nut
c                if(utwiss(mfitr1,m,i)*utwiss(mfitr4,m,i)
c     $               -utwiss(mfitr2,m,i)*utwiss(mfitr3,m,i)
c     1               < 1.d0)then
c                  ttrans(m)=.true.
c                  exit
c                endif
c              enddo
c            endif
c          enddo
c          exit
c        endif
      do m=-nfam,nfam
        ttrans(m)=.false.
      enddo
      error=.false.
      i=1
      do j=1,nfcol
        ka=flv%kfitp(j)
        kf=flv%kfit(ka)
        if(kf > mfitchi3)then
          cycle
        endif
        kp=flv%ifitp(ka)
        kp1=flv%ifitp1(ka)
        if(flv%mfitp(ka) /= 0)then
          mp=(abs(flv%mfitp(ka))-1)/2
          do20:     do idp=nfam1,nfam
            vf=flv%fitval(ka)
            if((kfam(idp) == 0 .and. idp >= -mp .and. idp <= mp)
     $           .or. jfam(idp) >= -mp .and. jfam(idp) <= mp)then
              if(kf >= mfitdx .and. kf <= mfitdpy .and.
     $             (idp < -nfr .or. idp > nfr ))then
c     $             .or. inicond .and. idp /= 0))then
                cycle
              elseif(idp /= 0 .and. kf > mfittry)then
                cycle
              endif
              if(idp /= 0 .or. mp*2+1 /= abs(flv%mfitp(ka)))then
                if(kp == nlat)then
                  if(idp >= -1 .and. idp <= 1)then
                    if(kf == mfitnx)then
                      vf=vf+(dp(idp)-dp0)*xixf
                    endif
                    if(kf == mfitny)then
                      vf=vf+(dp(idp)-dp0)*xiyf
                    endif
                  endif
                endif
                maxfit=flv%mfitp(ka) < 0
                if(kp /= kp1)then
                  maxfit=maxfit .and. .not. tftype1fit(kf)
                  kpb=min(kp,kp1)
                  kpe=max(kp,kp1)
                  if(maxfit)then
                    if(kf <= mfitdpy
     $                   .or. kf >= mfitpex .and.
     $                   kf <= mfitgmz)then
                      call tfpeak(idp,kf,kpb,kpe,ipeak,vpeak,npeak)
                      do k=1,npeak
                        ip=ipeak(k)
                        if(ip <= 0)then
                          cycle do20
                        endif
                        df1(i)=tdfun1(vf,vpeak(k),kf,maxfit,idp,ttrans(idp))
                        if(df1(i) /= 0.d0)then
                          iqcol(i)=j
                          lfp(1,i)=ip
                          lfp(2,i)=0
                          kdp(i)=idp
                          i=i+1
                          if(i > maxcond)then
                            error=.true.
                            return
                          endif
                        endif
                      enddo
                    endif
                  else
                    if(kf <= mfitbmagz .or.
     $                   (kf >= mfitleng .and. kf <= mfitchi3))then
                      maxfit=flv%mfitp(ka) < 0
                      select case (kf)
                      case (mfitbmagx,mfitbmagy,mfitbmagz)
                        ve=tgetbmagu2(kpb,kpe,idp,kf)
                        vb=vf
                      case default
                        vb=tgfun(kf,kpb,idp)
                        ve=tgfun(kf,kpe,idp)
                      end select
                      vf1=vb
                      call tfgetfitval(nlist(kf),kpb,kpe,dp(idp),
     $                     iuid(idp),kfam(idp),vb,ve,vf,vf1,irtc)
                      if(irtc == -1)then
                        cycle
                      endif
                      if(maxfit)then
                        ve=ve-vf1
                        vf1=vf
                      endif
                      df1(i)=tdfun1(vf1,ve,kf,maxfit,idp,ttrans(idp))
c                      write(*,'(a,3i5,l2,1p10g12.4)')'tdfun ',ka,j,kf,maxfit,vf1,ve,df1(i)
                      if(cell .and. ka > nfc0 .and. ka <= nfc0+2 .and. abs(df1(i)) < abmax)then
                        cycle
                      endif
                      if(.not. maxfit)then
                        select case (kf)
                        case (mfitbx,mfitby,mfitbz,mfitgmx,mfitgmy,mfitgmz)
                          df1(i)=df1(i)+log(vf)
                        case (mfitbmagx,mfitbmagy,mfitbmagz)
                        case default
                          df1(i)=df1(i)+vf
                        end select
                      endif
                      iqcol(i)=j
                      lfp(1,i)=kpe
                      lfp(2,i)=kpb
                      kdp(i)=idp
                      i=i+1
                      if(i > maxcond)then
                        error=.true.
                        return
                      endif
                    endif
                  endif
                else
                  v=tgfun(kf,kp,idp)
                  vf1=vf
                  call tfgetfitval(nlist(kf),kp,0,dp(idp),
     $                 iuid(idp),kfam(idp),vf,v,vf,vf1,irtc)
                  if(irtc == -1)then
                    cycle
                  endif
                  df1(i)=tdfun1(vf1,v,kf,maxfit,idp,ttrans(idp))
                  if(maxfit .and. df1(i) == 0.d0)then
                    cycle
                  endif
                  iqcol(i)=j
                  lfp(1,i)=kp
                  lfp(2,i)=0
                  kdp(i)=idp
                  i=i+1
                  if(i > maxcond)then
                    error=.true.
                    return
                  endif
                endif
              endif
            endif
          enddo do20
        endif
      enddo
      nqcola=i-1
      nqcola1=nqcola
      call tffsfitfun(nqcola,df1,iqcol,kdp,maxcond,error)
      return
      end

      subroutine tfgetfitval(funname,kp,kp1,dp,iuid,kfam,vf,v,vf0,vf1,irtc)
      use tfstk
      use ffs
      use ffs_pointer, only:mult
      use ffs_fit, only:inicond
      use tffitcode
      use eeval
      implicit none
      integer*4 ,intent(in):: kp,kp1,iuid,kfam
      integer*4 ,intent(out):: irtc
      real*8 ,intent(in):: dp,vf,v,vf0
      character*8 ,intent(in):: funname
      type (sad_descriptor) kx
      type (sad_dlist),pointer,save::klv,klv1
      type (sad_rlist),pointer,save::klid
      integer*4 lenw,itfuplevel,itfdownlevel
      real*8 vf1
      logical retry,retry1
      character*(MAXPNAME+8) name,name1
      integer*4 level,ln,ln1
      type (sad_descriptor) ,save::kfv,kfv1,kfid
      integer*8 ifvloc,ifvfun,ifvloc1
      save ifvloc,ifvfun,ifvloc1
      data kfv%k /0/
      if(kfv%k == 0)then
        kfid=kxavaloc(0,2,klid)
        kfv =kxadaloc(0,5,klv)
        klv%head=dtfcopy(kxsymbolz('`FitValue',9))
        ifvloc=ktsalocb(0,'                ',MAXPNAME+8)
        ifvloc1=ktsalocb(0,'                ',MAXPNAME+8)
        ifvfun=ktsalocb(0,'        ',MAXPNAME)
        klv%body(1)=ktfstring+ifvloc
        klv%body(2)=ktfstring+ifvfun
        klv%dbody(3)=kfid
        kfv1=kxadaloc(0,7,klv1)
        klv1%head=dtfcopy(klv%head)
        klv1%body(1)=ktfstring+ktfcopy1(ifvloc)
        klv1%body(2)=ktfstring+ktfcopy1(ifvloc1)
        klv1%body(3)=ktfstring+ktfcopy1(ifvfun)
        klv1%dbody(4)=dtfcopy(kfid)
      endif
      irtc=0
      vf1=vf
      call tfpadstr(funname,ifvfun+1,len_trim(funname))
      ilist(1,ifvfun)=len_trim(funname)
c      call tfdebugprint(kfid,'gfv-1',1)
      klid%rbody(1)=dble(merge(iuid,kfam,inicond))
      klid%rbody(2)=dp
      call elname(kp,name)
      ln=lenw(name)
      if(kp1 == 0)then
        retry1=.false.
        klv%rbody(4)=vf
        klv%rbody(5)=v
      else
        call elname(kp1,name1)
        ln1=lenw(name1)
        retry1=kp1 /= nlat
        klv1%rbody(5)=vf0
        klv1%rbody(6)=vf
        klv1%rbody(7)=v
      endif
      retry=kp /= nlat
 100  call tfpadstr(name,ifvloc+1,ln)
      ilist(1,ifvloc)=ln
      if(kp1 /= 0)then
        call tfpadstr(name1,ifvloc1+1,ln1)
        ilist(1,ifvloc1)=ln1
      endif
      call tclrfpe
      levele=itfuplevel()
      if(kp1 == 0)then
c        call tfdebugprint(kfv,'FitValue-1',1)
        kx=tfleval(klv,.true.,irtc)
      else
c        call tfdebugprint(kfv1,'FitValue-2',1)
        kx=tfleval(klv1,.true.,irtc)
      endif
c      call tfdebugprint(kx,'==> ',1)
      level=itfdownlevel()
 110  if(irtc /= 0)then
        if(ierrorprint /= 0)then
          call tfaddmessage(' ',0,6)
        endif
        call termes(6,
     $       'Error in FitValue '//
     $       funname//' at '//name,' ')
      elseif(ktfrealq(kx,vf1))then
      elseif(kx%k == ktfoper+mtfnull)then
        irtc=-1
      elseif(retry)then
        retry=.false.
        if(mult(kp) == 0)then
c     Generate singlet element name with suffix number(.###)
          call elnameK(kp,name)
          ln=lenw(name)
          go to 100
        elseif(nelvx(ilist(kp,ifele1))%klp == kp)then
c     Remove suffix number(.###) if head of multiple elements
c     Note: index(name,'.') > 0 if mult(kp) != 0
          ln=index(name,'.')-1
          name(ln+1:) = ' '
          go to 100
        elseif(retry1)then
c     No more candidate for 1st argument(name),
c     however, we need scan candidates for 2nd argument(name1)
          go to 110
        endif
        vf1=vf
      elseif(retry1)then
        retry1=.false.
c     Reset `kp'-element name in name(1:ln)
        retry=kp /= nlat
        if(retry)then
          call elname(kp,name)
          ln=lenw(name)
        endif
        if(mult(kp1) == 0)then
c     Generate singlet element name with suffix number(.###)
          call elnameK(kp1,name1)
          ln1=lenw(name1)
          go to 100
        elseif(nelvx(ilist(kp1,ifele1))%klp == kp1)then
c     Remove suffix number(.###) if singlet or head of multiple elements
c     Note: index(name1,'.') > 0 if kp1 != 0
          ln1=index(name1,'.')-1
          name1(ln1+1:) = ' '
          go to 100
        endif
      endif
      return
      end

      subroutine tffsfitfun(nqcol,df,iqcol,kdp,maxcond,error)
      use tfstk
      use tfcsi, only:icslfno
      use eeval
      implicit none
      integer*4 ,intent(in):: maxcond
      integer*4 ,intent(inout):: nqcol
      integer*4 ,intent(out):: iqcol(maxcond),kdp(maxcond)
      real*8 ,intent(out):: df(maxcond)
      logical*4 ,intent(out):: error
      type (sad_rlist), pointer :: klx
      type (sad_descriptor) kx
      type (sad_descriptor) ,save :: kff
      data kff%k /0/
      integer*4 itfuplevel,itfdownlevel,i,m,level,irtc
      if(kff%k == 0)then
        kff=kxsymbolz('`FitFunction',12)
      endif
      level=itfuplevel()
      kx%k=0
      kx=tfsyeval(kff,irtc)
c      call tfdebugprint(kx,'fitfun',3)
      if(irtc /= 0)then
        level=itfdownlevel()
        call tfaddmessage(' ',2,icslfno())
        call termes(6,'Error in FitFunction ',' ')
        error=.true.
        return
      elseif(ktfrealq(kx,df(nqcol+1)))then
        if(nqcol >= maxcond)then
          error=.true.
          level=itfdownlevel()
          return
        endif
        nqcol=nqcol+1
        iqcol(nqcol)=-1
        kdp(nqcol)=0
      elseif(tfreallistq(kx%k,klx))then
        m=klx%nl
        if(m > 0)then
          if(m+nqcol > maxcond)then
            error=.true.
            level=itfdownlevel()
            return
          endif
          df(nqcol+1:nqcol+m)=klx%rbody(1:m)
          kdp(nqcol+1:nqcol+m)=0
          do i=1,m
            iqcol(nqcol+i)=-i
          enddo
          nqcol=nqcol+m
        endif
      endif
c      call tfmemcheckprint('FitFunction-end',.true.,irtc)
      level=itfdownlevel()
      error=.false.
      return
      end

      real*8 function tdfun1(vf,v,kf,maxfit,kdp,ttrans)
      use macmath
      use tffitcode
      implicit none
      real*8 ,parameter::factor=1.d0
      integer*4 ,intent(in):: kf,kdp
      real*8 ,intent(in):: vf,v
      logical*4 ,intent(in):: maxfit,ttrans
      real*8 vfa
      select case (kf)
      case (mfitbx,mfitby,mfitbz,mfitgmx,mfitgmy,mfitgmz)
        if(maxfit)then
          if(v > vf)then
            tdfun1=log(vf*factor/v)
          else
            tdfun1=0.d0
          endif
        else
          tdfun1=log(vf/v)
        endif

      case (mfitax,mfitay,mfitaz)
        if(maxfit)then
          vfa=abs(vf)
          if(v > vfa)then
            tdfun1=atan(vfa*factor)-atan(v)
          elseif(v < -vfa)then
            tdfun1=-atan(vfa*factor)-atan(v)
          else
            tdfun1=0.d0
          endif
        else
          tdfun1=atan(vf)-atan(v)
        endif

      case (mfitex,mfitey)
        if(maxfit)then
          vfa=abs(vf)
          if(v > vfa)then
            tdfun1=vfa-v
          else
            tdfun1=max(-vfa-v,0.d0)
          endif
        else
          tdfun1=merge(-vf,vf,kdp < 0)-v
        endif
        return

      case (mfitchi1,mfitchi2,mfitchi3)
        if(maxfit)then
          vfa=abs(vf)
          if(v > vfa)then
            tdfun1=vfa-v
          else
            tdfun1=max(-vfa-v,0.d0)
          endif
        else
          tdfun1=vf-v
        endif
        do while(tdfun1 < -pi)
          tdfun1=tdfun1+pi2
        enddo
        do while(tdfun1 > pi)
          tdfun1=tdfun1-pi2
        enddo
        
      case (mfitnx,mfitny,mfitnz)
c     vf1=pi2*(anint(vf/pi2)+sign(.5d0*sin(.5d0*vf)**2,sin(vf)))
c     v1=pi2*(anint(v/pi2)+sign(.5d0*sin(.5d0*v)**2,sin(v)))
        if(maxfit)then
          tdfun1=min(0.d0,vf-v)
        else
          tdfun1=vf-v
        endif
        if(ttrans)then
          tdfun1=tdfun1-anint(tdfun1/pi2)*pi2
        endif

c      case (mfitbmagx,mfitbmagy,mfitbmagz)
c        tdfun1=sqrt(vf-1.d0)-sqrt(v-1.d0)
c
      case default
        if(maxfit)then
          vfa=abs(vf)
          if(v > vfa)then
            tdfun1=vfa-v
          else
            tdfun1=max(-vfa-v,0.d0)
          endif
        else
          tdfun1=vf-v
        endif
c        write(*,'(a,i5,l2,1p10g12.4)')'tdfun1 ',kf,maxfit,vf,v,tdfun1
      end select
      return

      end

      subroutine tfpeak(idp,kf,ibegin,iend,ipeak,vpeak,npeak)
      use ffs_pointer
      use tffitcode
      use gfun
      implicit none
      integer*4 ,intent(in)::ibegin,iend,kf,npeak,idp
      integer*4 ,intent(out)::ipeak(npeak)
      real*8 ,intent(out)::vpeak(npeak)
      integer*4 i,j,k
      real*8 va,va0,va1
      vpeak=0.d0
      ipeak=0
      va0=0.d0
      va=abs(tgfun(kf,ibegin,idp))
      do i=ibegin,iend
        va1=abs(tgfun(kf,min(i+1,iend),idp))
c        if(kf == mfitgmy)then
c          write(*,*)'tfpeak ',i,va,va0,va1
c        endif
        if(va > va0 .and. va >= va1)then
          do j=1,npeak
            if(va > abs(vpeak(j)))then
              do k=npeak,j+1,-1
                vpeak(k)=vpeak(k-1)
                ipeak(k)=ipeak(k-1)
              enddo
              vpeak(j)=abs(tgfun(kf,i,idp))
              ipeak(j)=i
              exit
            endif
          enddo
        endif
        va0=va
        va=va1
      enddo
      return
      end

      subroutine qdcell(dtwiss,kk,ll,idp,iv,ctrans,iclast,nfam,nut,disp,nzcod)
      use tfstk
      use ffs
      use ffs_pointer
      use tffitcode
      implicit none
      integer*4 ,intent(in):: nfam,nut,kk,ll,idp,iv
      integer*4 ,intent(inout):: iclast(-nfam:nfam)
      real*8 ,intent(out):: dtwiss(mfittry)
      real*8 ,intent(inout):: ctrans(27,-nfam:nfam)
      integer*4 k,l
      real*8 dpsix,dpsiy,cosmux,cosmuy,sinmux,sinmuy,
     $     bxr,byr,trx,try,
     $     x11,x12,x21,x22,y11,y12,y21,y22,
     $     dx11,dx12,dx21,dx22,dy11,dy12,dy21,dy22,
     $     detr,r11,r12,r21,r22,c1,ddetr,ddetr0,
     $     dr11,dr12,dr21,dr22,
     $     s11,s12,s13,s14,s21,s22,s23,s24,
     $     s31,s32,s33,s34,s41,s42,s43,s44,
     $     sr11,sr12,sr13,sr14,sr21,sr22,sr23,sr24,
     $     sr31,sr32,sr33,sr34,sr41,sr42,sr43,sr44,
     $     tm11,tm12,tm13,tm14,tm21,tm22,tm23,tm24,
     $     tm31,tm32,tm33,tm34,tm41,tm42,tm43,tm44,
     $     dalphx,dbetax,dalphy,dbetay,
     $     dsr11,dsr12,dsr13,dsr14,dsr21,dsr22,dsr23,dsr24,
     $     dsr31,dsr32,dsr33,dsr34,dsr41,dsr42,dsr43,dsr44,
     $     dtm11,dtm12,dtm13,dtm14,dtm21,dtm22,dtm23,dtm24,
     $     dtm31,dtm32,dtm33,dtm34,dtm41,dtm42,dtm43,dtm44,
     $     ds11,ds12,ds13,ds14,ds21,ds22,ds23,ds24,
     $     ds31,ds32,ds33,ds34,ds41,ds42,ds43,ds44,
     $     da11,da12,da21,da22,da13,da23,
     $     db11,db12,db21,db22,db13,db23,
     $     dmux,dmuy,sinl,cosl,
     $     dr1,dr2,dr3,dr4,
     $     cod1,cod2,cod3,cod4,
     $     dcod1,dcod2,dcod3,dcod4,
     $     ddcod1,ddcod2,ddcod3,ddcod4,
     $     dea,depa,deb,depb,dc,dc1,detimx,detimy,
     $     ts11,ts21,ts31,ts41,ts12,ts22,ts32,ts42,
     $     dts11,dts21,dts31,dts41,dts12,dts22,dts32,dts42,
     $     tds11,tds21,tds31,tds41,tds12,tds22,tds32,tds42,
     $     ts13,ts23,ts33,ts43,ts14,ts24,ts34,ts44,
     $     dts13,dts23,dts33,dts43,dts14,dts24,dts34,dts44,
     $     tds13,tds23,tds33,tds43,tds14,tds24,tds34,tds44,
     $     ex,epx,ey,epy,pex,pepx,pey,pepy,dex,depx,dey,depy,
     $     thx,thy,dcosx,dcosy
c      logical*4 cell0,disp,nzcod,htrx,htry,normal
      logical*4 cell0,disp,nzcod,normal,xstab,ystab
c-deb
c      print *,'--------- enter qdcell --k,l=',k,l,'-------'
c          ------------------------------------
      k=itwissp(kk)
      l=itwissp(ll)
c      htrx=.false.
c      htry=.false.
      cell0=cell
      if(cell)then
        dpsix = utwiss(mfitnx,idp,nut) - utwiss(mfitnx,idp,1)
        dpsiy = utwiss(mfitny,idp,nut) - utwiss(mfitny,idp,1)
        thx=tan(0.5d0*dpsix)
        sinmux=2.d0*thx/(1.d0+thx**2)
        dcosx=thx*sinmux
        cosmux=1.d0-dcosx
        detimx=2.d0*dcosx
        thy=tan(0.5d0*dpsiy)
        sinmuy=2.d0*thy/(1.d0+thy**2)
        dcosy=thy*sinmuy
        cosmuy=1.d0-dcosy
        detimy=2.d0*dcosy
c        cosmux = cos(dpsix)
c        sinmux = sin(dpsix)
c        cosmuy = cos(dpsiy)
c        sinmuy = sin(dpsiy)
c        detimx = 4.d0*sin(.5d0*dpsix)**2
c        detimy = 4.d0*sin(.5d0*dpsiy)**2
        bxr=sqrt(utwiss(mfitbx,idp,nut)/utwiss(mfitbx,idp,1))
        trx=(bxr+1.d0/bxr)*cosmux
     1      +(bxr*utwiss(mfitax,idp,1)-utwiss(mfitax,idp,nut)/bxr)
     $       *sinmux
        byr=sqrt(utwiss(mfitby,idp,nut)/utwiss(mfitby,idp,1))
        try=(byr+1.d0/byr)*cosmuy
     1      +(byr*utwiss(mfitay,idp,1)-utwiss(mfitay,idp,nut)/byr)
     $       *sinmuy
        xstab=abs(trx) < 2.d0
        ystab=abs(try) < 2.d0
        cell=xstab .and. ystab
c$$$        htrx=abs(trx) >= 2.d0
c$$$        htry=abs(try) >= 2.d0
c$$$        if(htrx)then
c$$$          akx=(utwiss(mfitax,idp,nut)-utwiss(mfitax,idp,1))
c$$$     $         /utwiss(mfitbx,idp,1)
c$$$        endif
c$$$        if(htry)then
c$$$          aky=(utwiss(mfitay,idp,nut)-utwiss(mfitay,idp,1))
c$$$     $         /utwiss(mfitby,idp,1)
c$$$        endif
c$$$        cell=.not. (htrx .or. htry)
      else
        sinmux=0.d0
        cosmux=1.d0
        sinmuy=0.d0
        cosmuy=1.d0
        dpsix=0.d0
        dpsiy=0.d0
        detimx=0.d0
        detimy=0.d0
      endif
      call qdtwis(dtwiss,ctrans,iclast,
     $     kk,ll,idp,iv,nfam,nut,disp,nzcod)
c          ------------------------------------
      cell=cell0
      if(cell)then
C======Get Transfer Matrix tm**  =========================
        x11=cosmux+utwiss(mfitax,idp,l)*sinmux
        x22=cosmux-utwiss(mfitax,idp,l)*sinmux
        x12=utwiss(mfitbx,idp,l)*sinmux
        x21=-(1.d0+utwiss(mfitax,idp,l)**2)*sinmux/utwiss(mfitbx,idp,l)
c     x21=-((utwiss(mfitax,idp,l)-utwiss(mfitax,idp,1))*cosmux
c     $     +(1.d0+utwiss(mfitax,idp,1)*utwiss(mfitax,idp,l))*sinmux)
c     $     /utwiss(mfitbx,idp,l)
        y11=cosmuy+utwiss(mfitay,idp,l)*sinmuy
        y22=cosmuy-utwiss(mfitay,idp,l)*sinmuy
        y12=utwiss(mfitby,idp,l)*sinmuy
c     y21=-((utwiss(mfitay,idp,l)-utwiss(mfitay,idp,1))*cosmuy
c     $     +(1.d0+utwiss(mfitay,idp,1)*utwiss(mfitay,idp,l))*sinmuy)
c     $     /utwiss(mfitby,idp,l)
        y21=-(1.d0+utwiss(mfitay,idp,l)**2)*sinmuy/utwiss(mfitby,idp,l)
        dalphx = (1.d0+utwiss(mfitax,idp,l)**2)*dtwiss(mfitax)
        dbetax = 0.5d0*dtwiss(mfitbx)
        dalphy = (1.d0+utwiss(mfitay,idp,l)**2)*dtwiss(mfitay)
        dbetay = 0.5d0*dtwiss(mfitby)
C     
C     ( transformation matrix at point l )
        r11 = utwiss(mfitr1,idp,l)
        r12 = utwiss(mfitr2,idp,l)
        r21 = utwiss(mfitr3,idp,l)
        r22 = utwiss(mfitr4,idp,l)
        detr=r11*r22-r12*r21
        c1 = sqrt(1.d0 - detr)
        s11 = c1
        s12 = 0.d0
        s13 = -r22
        s14 =  r12
        s21 = 0.d0
        s22 = c1
        s23 =  r21
        s24 = -r11
        s31 =  r11
        s32 =  r12
        s33 = c1
        s34 = 0.d0
        s41 =  r21
        s42 =  r22
        s43 = 0.d0
        s44 = c1
        sr11 = c1
        sr12 = 0.d0
        sr13 =  r22
        sr14 = -r12
        sr21 = 0.d0
        sr22 = c1
        sr23 = -r21
        sr24 =  r11
        sr31 = -r11
        sr32 = -r12
        sr33 = c1
        sr34 = 0.d0
        sr41 = -r21
        sr42 = -r22
        sr43 = 0.d0
        sr44 = c1
        tm11 =sr11*(x11*s11)
     1       +sr13*(y11*s31 + y12*s41)
     1       +sr14*(y21*s31 + y22*s41)
        tm12 =sr11*(x12*s22)
     1       +sr13*(y11*s32 + y12*s42)
     1       +sr14*(y21*s32 + y22*s42)
        tm13 =sr11*(x11*s13 + x12*s23)
     1       +sr13*(y11*s33)
     $       +sr14*(y21*s33)
        tm14 =sr11*(x11*s14 + x12*s24)
     1       +sr13*(y12*s44)
     $       +sr14*(y22*s44)
        tm21 =sr22*(x21*s11)
     1       +sr23*(y11*s31 + y12*s41)
     1       +sr24*(y21*s31 + y22*s41)
        tm22 =sr22*(x22*s22)
     1       +sr23*(y11*s32 + y12*s42)
     1       +sr24*(y21*s32 + y22*s42)
        tm23 =sr22*(x21*s13 + x22*s23)
     1       +sr23*(y11*s33)
     $       +sr24*(y21*s33)
        tm24 =sr22*(x21*s14 + x22*s24)
     1       +sr23*(y12*s44)
     $       +sr24*(y22*s44)
        tm31 =sr31*(x11*s11)
     $       +sr32*(x21*s11)
     1       +sr33*(y11*s31 + y12*s41)
        tm32 =sr31*(x12*s22)
     $       +sr32*(x22*s22)
     1       +sr33*(y11*s32 + y12*s42)
        tm33 =sr31*(x11*s13 + x12*s23)
     1       +sr32*(x21*s13 + x22*s23)
     1       +sr33*(y11*s33)
        tm34 =sr31*(x11*s14 + x12*s24)
     1       +sr32*(x21*s14 + x22*s24)
     1       +sr33*(y12*s44)
        tm41 =sr41*(x11*s11)
     $       +sr42*(x21*s11)
     1       +sr44*(y21*s31 + y22*s41)
        tm42 =sr41*(x12*s22)
     $       +sr42*(x22*s22)
     1       +sr44*(y21*s32 + y22*s42)
        tm43 =sr41*(x11*s13 + x12*s23)
     1       +sr42*(x21*s13 + x22*s23)
     1       +sr44*(y21*s33)
        tm44 =sr41*(x11*s14 + x12*s24)
     1       +sr42*(x21*s14 + x22*s24)
     1       +sr44*(y22*s44)
C     ( 4*4 transfer matrix( assuming correct COD found ))
C     
C     M  = s-inv. Md .s
C     
C======Get derivative of Transfer Matrix dtm**  ==========
C     
C     dM  = dsf-inv. Md .s + sf-inv. dMd .s
C     
        dx11=dbetax*x11
     1       +(-sinmux+utwiss(mfitax,idp,l)*cosmux)*dtwiss(mfitnx)
        dx12=dbetax*x12
     1       +utwiss(mfitbx,idp,l)*cosmux*dtwiss(mfitnx)
        dx21=-x21*dbetax+(
     1       -(utwiss(mfitax,idp,l)*sinmux+cosmux)*dalphx
     1       -(1.d0+utwiss(mfitax,idp,l)**2)*cosmux*dtwiss(mfitnx))
     1       /utwiss(mfitbx,idp,l)
        dx22=-dbetax*x22
     1       -sinmux*dalphx
     1       -(sinmux+utwiss(mfitax,idp,l)*cosmux)*dtwiss(mfitnx)
c---  
        dy11=dbetay*y11
     1       +(-sinmuy+utwiss(mfitay,idp,l)*cosmuy)*dtwiss(mfitny)
        dy12=dbetay*y12
     1       +utwiss(mfitby,idp,l)*cosmuy*dtwiss(mfitny)
        dy21=-y21*dbetay+(
     1       -(utwiss(mfitay,idp,l)*sinmuy+cosmuy)*dalphy
     1       -(1.d0+utwiss(mfitay,idp,l)**2)*cosmuy*dtwiss(mfitny))
     1       /utwiss(mfitby,idp,l)
        dy22=-dbetay*y22
     1       -sinmuy*dalphy
     1       -(sinmuy+utwiss(mfitay,idp,l)*cosmuy)*dtwiss(mfitny)
c---  
c     c^2+Det(R)=1   ---> dc = - d(Det(R))/(2c)
        dr11 = dtwiss(mfitr1)
        dr12 = dtwiss(mfitr2)
        dr21 = dtwiss(mfitr3)
        dr22 = dtwiss(mfitr4)
        ddetr = dr11*r22 + dr22*r11
     1       -dr12*r21 - dr21*r12
        dc1   = -0.5d0*ddetr/c1
        dsr11 = dc1
        dsr12 = 0.d0
        dsr13 =  dr22
        dsr14 = -dr12
        dsr21 = 0.d0
        dsr22 = dc1
        dsr23 = -dr21
        dsr24 =  dr11
        dsr31 = -dr11
        dsr32 = -dr12
        dsr33 = dc1
        dsr34 = 0.d0
        dsr41 = -dr21
        dsr42 = -dr22
        dsr43 = 0.d0
        dsr44 = dc1
        dtm11 =  dsr11*(x11*s11)+sr11*(dx11*s11)
     1       +dsr13*(y11*s31 + y12*s41)+sr13*(dy11*s31 + dy12*s41)
     1       +dsr14*(y21*s31 + y22*s41)+sr14*(dy21*s31 + dy22*s41)
        dtm12 =  dsr11*(x12*s22)+sr11*(dx12*s22)
     1       +dsr13*(y11*s32 + y12*s42)+sr13*(dy11*s32 + dy12*s42)
     1       +dsr14*(y21*s32 + y22*s42)+sr14*(dy21*s32 + dy22*s42)
        dtm13 =  dsr11*(x11*s13 + x12*s23)+sr11*(dx11*s13 + dx12*s23)
     1       +dsr13*(y11*s33)+sr13*(dy11*s33)
     1       +dsr14*(y21*s33)+sr14*(dy21*s33)
        dtm14 =  dsr11*(x11*s14 + x12*s24)+sr11*(dx11*s14 + dx12*s24)
     1       +dsr13*(y12*s44)+sr13*(dy12*s44)
     1       +dsr14*(y22*s44)+sr14*(dy22*s44)
        dtm21 =  
     1       +dsr22*(x21*s11)+sr22*(dx21*s11)
     1       +dsr23*(y11*s31 + y12*s41)+sr23*(dy11*s31 + dy12*s41)
     1       +dsr24*(y21*s31 + y22*s41)+sr24*(dy21*s31 + dy22*s41)
        dtm22 =  
     1       +dsr22*(x22*s22)+sr22*(dx22*s22)
     1       +dsr23*(y11*s32 + y12*s42)+sr23*(dy11*s32 + dy12*s42)
     1       +dsr24*(y21*s32 + y22*s42)+sr24*(dy21*s32 + dy22*s42)
        dtm23 =  
     1       +dsr22*(x21*s13 + x22*s23)+sr22*(dx21*s13 + dx22*s23)
     1       +dsr23*(y11*s33)+sr23*(dy11*s33 )
     1       +dsr24*(y21*s33)+sr24*(dy21*s33 )
        dtm24 =  
     1       +dsr22*(x21*s14 + x22*s24)+sr22*(dx21*s14 + dx22*s24)
     1       +dsr23*(y12*s44)+sr23*(dy12*s44)
     1       +dsr24*(y22*s44)+sr24*(dy22*s44)
        dtm31 =  dsr31*(x11*s11)+sr31*(dx11*s11)
     1       +dsr32*(x21*s11)+sr32*(dx21*s11)
     1       +dsr33*(y11*s31 + y12*s41)+sr33*(dy11*s31 + dy12*s41)
        dtm32 =  dsr31*(x12*s22)+sr31*(dx12*s22)
     1       +dsr32*(x22*s22)+sr32*(dx22*s22)
     1       +dsr33*(y11*s32 + y12*s42)+sr33*(dy11*s32 + dy12*s42)
        dtm33 =  dsr31*(x11*s13 + x12*s23)+sr31*(dx11*s13 + dx12*s23)
     1       +dsr32*(x21*s13 + x22*s23)+sr32*(dx21*s13 + dx22*s23)
     1       +dsr33*(y11*s33 )+sr33*(dy11*s33 )
        dtm34 =  dsr31*(x11*s14 + x12*s24)+sr31*(dx11*s14 + dx12*s24)
     1       +dsr32*(x21*s14 + x22*s24)+sr32*(dx21*s14 + dx22*s24)
     1       +dsr33*(y12*s44)+sr33*(dy12*s44)
        dtm41 =  dsr41*(x11*s11)+sr41*(dx11*s11)
     1       +dsr42*(x21*s11)+sr42*(dx21*s11)
     1       +dsr44*(y21*s31 + y22*s41)+sr44*(dy21*s31 + dy22*s41)
        dtm42 =  dsr41*(x12*s22)+sr41*(dx12*s22)
     1       +dsr42*(x22*s22)+sr42*(dx22*s22)
     1       +dsr44*(y21*s32 + y22*s42)+sr44*(dy21*s32 + dy22*s42)
        dtm43 =  dsr41*(x11*s13 + x12*s23)+sr41*(dx11*s13 + dx12*s23)
     1       +dsr42*(x21*s13 + x22*s23)+sr42*(dx21*s13 + dx22*s23)
     1       +dsr44*(y21*s33 )+sr44*(dy21*s33 )
        dtm44 =  dsr41*(x11*s14 + x12*s24)+sr41*(dx11*s14 + dx12*s24)
     1       +dsr42*(x21*s14 + x22*s24)+sr42*(dx21*s14 + dx22*s24)
     1       +dsr44*(y22*s44)+sr44*(dy22*s44)
C======Get derivative of transformation matrix dr*  ======
c     with periodic boundary condition
c     write(*,'(4(1p4g12.4/))')
c     $     tm11,tm12,tm13,tm14,
c     1     tm21,tm22,tm23,tm24,
c     1     tm31,tm32,tm33,tm34,
c     1     tm41,tm42,tm43,tm44,
c     1     dtm11,dtm12,dtm13,dtm14,
c     1     dtm21,dtm22,dtm23,dtm24,
c     1     dtm31,dtm32,dtm33,dtm34,
c     1     dtm41,dtm42,dtm43,dtm44
c     -------------------------------
        call qdmdia(tm11,tm12,tm13,tm14,
     1       tm21,tm22,tm23,tm24,
     1       tm31,tm32,tm33,tm34,
     1       tm41,tm42,tm43,tm44,
     1       dtm11,dtm12,dtm13,dtm14,
     1       dtm21,dtm22,dtm23,dtm24,
     1       dtm31,dtm32,dtm33,dtm34,
     1       dtm41,dtm42,dtm43,dtm44,
     1       dr1,dr2,dr3,dr4,r11)
c     -------------------------------
        dtwiss(mfitr1) = dr1
        dtwiss(mfitr2) = dr2
        dtwiss(mfitr3) = dr3
        dtwiss(mfitr4) = dr4
        ddetr0=dr1*r22+r11*dr4-dr2*r21-r12*dr3
        normal=utwiss(mfitdetr,idp,l) < 1.d0
        if(normal)then
          dtwiss(mfitdetr)=ddetr0
        else
          dtwiss(mfitdetr)=-ddetr0
        endif
C======Get derivative of transformation matrix ds**  =====
c     with periodic boundary condition
        dc   = -0.5d0*(ddetr0)/c1
c     write(*,'(a,1p5g15.7)')'qdcell-r ',dr1,dr2,dr3,dr4,ddetr0
        ds11 = dc
        ds12 = 0.d0
        ds13 = -dr4
        ds14 = dr2
        ds21 = 0.d0
        ds22 = dc
        ds23 = dr3
        ds24 = -dr1
        ds31 = dr1
        ds32 = dr2
        ds33 = dc
        ds34 = 0.d0
        ds41 = dr3
        ds42 = dr4
        ds43 = 0.d0
        ds44 = dc
C------dA,dB:derivative of diagonalized transfer matrix ------
        ts11=tm11*s11-tm13*s31-tm14*s41
        ts21=tm21*s11-tm23*s31-tm24*s41
        ts31=tm31*s11-tm33*s31-tm34*s41
        ts41=tm41*s11-tm43*s31-tm44*s41
        ts12=tm12*s22-tm13*s32-tm14*s42
        ts22=tm22*s22-tm23*s32-tm24*s42
        ts32=tm32*s22-tm33*s32-tm34*s42
        ts42=tm42*s22-tm43*s32-tm44*s42
        dts11=dtm11*s11-dtm13*s31-dtm14*s41
        dts21=dtm21*s11-dtm23*s31-dtm24*s41
        dts31=dtm31*s11-dtm33*s31-dtm34*s41
        dts41=dtm41*s11-dtm43*s31-dtm44*s41
        dts12=dtm12*s22-dtm13*s32-dtm14*s42
        dts22=dtm22*s22-dtm23*s32-dtm24*s42
        dts32=dtm32*s22-dtm33*s32-dtm34*s42
        dts42=dtm42*s22-dtm43*s32-dtm44*s42
        tds11=tm11*ds11-tm13*ds31-tm14*ds41
        tds21=tm21*ds11-tm23*ds31-tm24*ds41
        tds31=tm31*ds11-tm33*ds31-tm34*ds41
        tds41=tm41*ds11-tm43*ds31-tm44*ds41
        tds12=tm12*ds22-tm13*ds32-tm14*ds42
        tds22=tm22*ds22-tm23*ds32-tm24*ds42
        tds32=tm32*ds22-tm33*ds32-tm34*ds42
        tds42=tm42*ds22-tm43*ds32-tm44*ds42
        ts13=-tm11*s13-tm12*s23+tm13*s33
        ts23=-tm21*s13-tm22*s23+tm23*s33
        ts33=-tm31*s13-tm32*s23+tm33*s33
        ts43=-tm41*s13-tm42*s23+tm43*s33
        ts14=-tm11*s14-tm12*s24+tm14*s44
        ts24=-tm21*s14-tm22*s24+tm24*s44
        ts34=-tm31*s14-tm32*s24+tm34*s44
        ts44=-tm41*s14-tm42*s24+tm44*s44
        dts13=-dtm11*s13-dtm12*s23+dtm13*s33
        dts23=-dtm21*s13-dtm22*s23+dtm23*s33
        dts33=-dtm31*s13-dtm32*s23+dtm33*s33
        dts43=-dtm41*s13-dtm42*s23+dtm43*s33
        dts14=-dtm11*s14-dtm12*s24+dtm14*s44
        dts24=-dtm21*s14-dtm22*s24+dtm24*s44
        dts34=-dtm31*s14-dtm32*s24+dtm34*s44
        dts44=-dtm41*s14-dtm42*s24+dtm44*s44
        tds13=-tm11*ds13-tm12*ds23+tm13*ds33
        tds23=-tm21*ds13-tm22*ds23+tm23*ds33
        tds33=-tm31*ds13-tm32*ds23+tm33*ds33
        tds43=-tm41*ds13-tm42*ds23+tm43*ds33
        tds14=-tm11*ds14-tm12*ds24+tm14*ds44
        tds24=-tm21*ds14-tm22*ds24+tm24*ds44
        tds34=-tm31*ds14-tm32*ds24+tm34*ds44
        tds44=-tm41*ds14-tm42*ds24+tm44*ds44
        da11 = ds11*ts11+ ds13*ts31+ ds14*ts41
     1       + s11*dts11+ s13*dts31+ s14*dts41
     1       + s11*tds11+ s13*tds31+ s14*tds41
        da12 = ds11*ts12+ ds13*ts32+ ds14*ts42
     1       + s11*dts12+ s13*dts32+ s14*dts42
     1       + s11*tds12+ s13*tds32+ s14*tds42
        da21 = ds22*ts21+ ds23*ts31+ ds24*ts41
     1       + s22*dts21+ s23*dts31+ s24*dts41
     1       + s22*tds21+ s23*tds31+ s24*tds41
        da22 = ds22*ts22+ ds23*ts32+ ds24*ts42
     1       + s22*dts22+ s23*dts32+ s24*dts42
     1       + s22*tds22+ s23*tds32+ s24*tds42
        db11 = ds31*ts13+ ds32*ts23+ ds33*ts33
     $       + s31*dts13+ s32*dts23+ s33*dts33
     $       + s31*tds13+ s32*tds23+ s33*tds33
        db12 = ds31*ts14+ ds32*ts24+ ds33*ts34
     $       + s31*dts14+ s32*dts24+ s33*dts34
     $       + s31*tds14+ s32*tds24+ s33*tds34
        db21 = ds41*ts13+ ds42*ts23+ ds44*ts43
     $       + s41*dts13+ s42*dts23+ s44*dts43
     $       + s41*tds13+ s42*tds23+ s44*tds43
        db22 = ds41*ts14+ ds42*ts24+ ds44*ts44
     $       + s41*dts14+ s42*dts24+ s44*dts44
     $       + s41*tds14+ s42*tds24+ s44*tds44
c     -deb
c     print *,'da11,da12=',da11,da12
c     print *,'da21,da22=',da21,da22
C======Get derivative of twiss parameters  ===============
C     
C     Changed by K. O.  6/5/1992
C     Dtwiss(mfitnx) and dtwiss(mfitny) below are still not perfect.
C     The change of R at l is not included yet.
C     
        if(sinmux/=0.d0) then
          dmux=-0.5d0*(da11+da22)/sinmux
          dtwiss(mfitax)=(.5d0*(da11-da22)
     $         -utwiss(mfitax,idp,l)*cosmux*dmux)/sinmux
          dtwiss(mfitbx)=(da12-utwiss(mfitbx,idp,l)*cosmux*dmux)/x12
        else
          dmux=0.d0
        endif
        if(k < l)then
          sinl=sin(dpsix-utwiss(mfitnx,idp,l))
          cosl=cos(dpsix-utwiss(mfitnx,idp,l))
          dtwiss(mfitnx)=dmux
     $         +sinl*((cosl+utwiss(mfitax,idp,l)*sinl)*dtwiss(mfitbx)
     1         -sinl*dtwiss(mfitax))
        else
          sinl=sin(utwiss(mfitnx,idp,l))
          cosl=cos(utwiss(mfitnx,idp,l))
          dtwiss(mfitnx)=
     $         -sinl*((cosl-utwiss(mfitax,idp,l)*sinl)*dtwiss(mfitbx)
     1         +sinl*dtwiss(mfitax))
        endif
        dtwiss(mfitax)=dtwiss(mfitax)/(1.d0+utwiss(mfitax,idp,l)**2)
        if(sinmuy/=0.d0) then
          dmuy=-0.5d0*(db11+db22)/sinmuy
          dtwiss(mfitay)=(.5d0*(db11-db22)
     $         -utwiss(mfitay,idp,l)*cosmuy*dmuy)/sinmuy
          dtwiss(mfitby)=(db12-utwiss(mfitby,idp,l)*cosmuy*dmuy)/y12
        else
          dmuy=0.d0
        end if
        if(k < l)then
C     
C     A bug was found below on 7/1/1994. sin(dpsix- ....)
C     
          sinl=sin(dpsiy-utwiss(mfitny,idp,l))
          cosl=cos(dpsiy-utwiss(mfitny,idp,l))
          dtwiss(mfitny)=dmuy+sinl
     1         *((cosl+utwiss(mfitay,idp,l)*sinl)*dtwiss(mfitby)
     1         -sinl*dtwiss(mfitay))
        else
          sinl=sin(utwiss(mfitny,idp,l))
          cosl=cos(utwiss(mfitny,idp,l))
          dtwiss(mfitny)=-sinl
     1         *((cosl-utwiss(mfitay,idp,l)*sinl)*dtwiss(mfitby)
     1         +sinl*dtwiss(mfitay))
        endif
        dtwiss(mfitay)=dtwiss(mfitay)/(1.d0+utwiss(mfitay,idp,l)**2)
c-------dtwiss(mfitax) = dalpha/(1+alpha^2)
c     dtwiss(mfitbx) = dbeta/beta
C======Get derivative of orbit  ==========================
c..........going to 2*2 world
        if(normal)then
          cod1 = s11*utwiss(mfitdx,idp,l)
     1         + s13*utwiss(mfitdy,idp,l)+s14*utwiss(mfitdpy,idp,l)
          cod2 = s22*utwiss(mfitdpx,idp,l)
     1         + s23*utwiss(mfitdy,idp,l)+s24*utwiss(mfitdpy,idp,l)
          cod3 = s31*utwiss(mfitdx,idp,l)+s32*utwiss(mfitdpx,idp,l)
     1         + s33*utwiss(mfitdy,idp,l)
          cod4 = s41*utwiss(mfitdx,idp,l)+s42*utwiss(mfitdpx,idp,l)
     1         + s44*utwiss(mfitdpy,idp,l)
          ddcod1 = s11*dtwiss(mfitdx)
     1         + s13*dtwiss(mfitdy)+s14*dtwiss(mfitdpy)
          ddcod2 = s22*dtwiss(mfitdpx)
     1         + s23*dtwiss(mfitdy)+s24*dtwiss(mfitdpy)
          ddcod3 = s31*dtwiss(mfitdx)+s32*dtwiss(mfitdpx)
     1         + s33*dtwiss(mfitdy)
          ddcod4 = s41*dtwiss(mfitdx)+s42*dtwiss(mfitdpx)
     1         + s44*dtwiss(mfitdpy)
        else
          cod3 = s11*utwiss(mfitdx,idp,l)
     1         + s13*utwiss(mfitdy,idp,l)+s14*utwiss(mfitdpy,idp,l)
          cod4 = s22*utwiss(mfitdpx,idp,l)
     1         + s23*utwiss(mfitdy,idp,l)+s24*utwiss(mfitdpy,idp,l)
          cod1 = s31*utwiss(mfitdx,idp,l)+s32*utwiss(mfitdpx,idp,l)
     1         + s33*utwiss(mfitdy,idp,l)
          cod2 = s41*utwiss(mfitdx,idp,l)+s42*utwiss(mfitdpx,idp,l)
     1         + s44*utwiss(mfitdpy,idp,l)
          ddcod3 = s11*dtwiss(mfitdx)
     1         + s13*dtwiss(mfitdy)+s14*dtwiss(mfitdpy)
          ddcod4 = s22*dtwiss(mfitdpx)
     1         + s23*dtwiss(mfitdy)+s24*dtwiss(mfitdpy)
          ddcod1 = s31*dtwiss(mfitdx)+s32*dtwiss(mfitdpx)
     1         + s33*dtwiss(mfitdy)
          ddcod2 = s41*dtwiss(mfitdx)+s42*dtwiss(mfitdpx)
     1         + s44*dtwiss(mfitdpy)
        endif
c..........dx = (I-M)^-1.(dM.x+delx)
        if(detimx/=0.d0) then
          dcod1 = ((1.d0-x22)*ddcod1+x12*ddcod2)/detimx
          dcod2 = (x21*ddcod1+(1.d0-x11)*ddcod2)/detimx
        else
c     write(*,*)'qdcell @ src/qdcell.f: ',
c     $        '2x2 horizontal matrix (I - A) is degenerated!',
c     $        '(FIXME)'
          dcod1=0.d0
          dcod2=0.d0
        end if
        if(detimy/=0.d0) then
          dcod3 = ((1.d0-y22)*ddcod3+y12*ddcod4)/detimy
          dcod4 = (y21*ddcod3+(1.d0-y11)*ddcod4)/detimy
        else
c     write(*,*)'qdcell @ src/qdcell.f: ',
c     $        '2x2 vetical matrix (I - B) is degenerated!',
c     $        '(FIXME)'
          dcod3=0.d0
          dcod4=0.d0
        end if
c..........going back to 4*4 world
        if(normal)then
          dtwiss(mfitdx)  = sr11*dcod1+sr13*dcod3+sr14*dcod4
          dtwiss(mfitdpx) = sr22*dcod2+sr23*dcod3+sr24*dcod4
          dtwiss(mfitdy)  = sr31*dcod1+sr32*dcod2+sr33*dcod3
          dtwiss(mfitdpy) = sr41*dcod1+sr42*dcod2+sr44*dcod4
        else
          dtwiss(mfitdy)  = sr11*dcod1+sr13*dcod3+sr14*dcod4
          dtwiss(mfitdpy) = sr22*dcod2+sr23*dcod3+sr24*dcod4
          dtwiss(mfitdx)  = sr31*dcod1+sr32*dcod2+sr33*dcod3
          dtwiss(mfitdpx) = sr41*dcod1+sr42*dcod2+sr44*dcod4
        endif
        if(normal)then
          da13=s11*dtwiss(mfitpex)
     $         +s13*dtwiss(mfitpey)+s14*dtwiss(mfitpepy)
          da23=s22*dtwiss(mfitpepx)
     $         +s23*dtwiss(mfitpey)+s24*dtwiss(mfitpepy)
          db13=s33*dtwiss(mfitpey)
     $         +s31*dtwiss(mfitpex)+s32*dtwiss(mfitpepx)
          db23=s44*dtwiss(mfitpepy)
     $         +s41*dtwiss(mfitpex)+s42*dtwiss(mfitpepx)
        else
          db13=s11*dtwiss(mfitpex)
     $         +s13*dtwiss(mfitpey)+s14*dtwiss(mfitpepy)
          db23=s22*dtwiss(mfitpepx)
     $         +s23*dtwiss(mfitpey)+s24*dtwiss(mfitpepy)
          da13=s33*dtwiss(mfitpey)
     $         +s31*dtwiss(mfitpex)+s32*dtwiss(mfitpepx)
          da23=s44*dtwiss(mfitpepy)
     $         +s41*dtwiss(mfitpex)+s42*dtwiss(mfitpepx)
        endif
        if( detimx/=0.d0 ) then
          dea =(da13+x12*da23-x22*da13)/detimx
          depa=(da23+x21*da13-x11*da23)/detimx
        else
          dea=0.d0
          depa=0.d0
        end if
        if( detimy/=0.d0 ) then
          deb =(db13+y12*db23-y22*db13)/detimy
          depb=(db23+y21*db13-y11*db23)/detimy
        else
          deb=0.d0
          depb=0.d0
        end if
        ex =utwiss(mfitex, idp,l)
        epx=utwiss(mfitepx,idp,l)
        ey =utwiss(mfitey, idp,l)
        epy=utwiss(mfitepy,idp,l)
        pex =sr11*ex  +sr13*ey  +sr14*epy
        pepx=sr22*epx +sr23*ey  +sr24*epy
        pey =sr31*ex  +sr32*epx +sr33*ey
        pepy=sr41*ex  +sr42*epx +sr44*epy
        dex =ds11*pex +ds13*pey +ds14*pepy
        depx=ds22*pepx+ds23*pey +ds24*pepy
        dey =ds31*pex +ds32*pepx+ds33*pey
        depy=ds41*pex +ds42*pepx+ds44*pepy
        dtwiss(mfitex) =dea+dex
        dtwiss(mfitepx)=depa+depx
        dtwiss(mfitey) =deb+dey
        dtwiss(mfitepy)=depb+depy
        if(normal)then
          dtwiss(mfitpex) =sr11*dea +sr13*deb +sr14*depb
          dtwiss(mfitpepx)=sr22*depa+sr23*deb +sr24*depb
          dtwiss(mfitpey) =sr31*dea +sr32*depa+sr33*deb
          dtwiss(mfitpepy)=sr41*dea +sr42*depa+sr44*depb
        else
          dtwiss(mfitpey) =sr11*dea +sr13*deb +sr14*depb
          dtwiss(mfitpepy)=sr22*depa+sr23*deb +sr24*depb
          dtwiss(mfitpex) =sr31*dea +sr32*depa+sr33*deb
          dtwiss(mfitpepx)=sr41*dea +sr42*depa+sr44*depb
        endif
      endif
      call dgammag(dtwiss(mfitgmx),dtwiss(mfitbmagx),
     $     dtwiss(mfitax),dtwiss(mfitbx),
     $     utwiss(mfitax,idp,l),utwiss(mfitbx,idp,l),
     $     twiss(ll,-1,mfitax),twiss(ll,-1,mfitbx))
      call dgammag(dtwiss(mfitgmy),dtwiss(mfitbmagy),
     $     dtwiss(mfitay),dtwiss(mfitby),
     $     utwiss(mfitay,idp,l),utwiss(mfitby,idp,l),
     $     twiss(ll,-1,mfitay),twiss(ll,-1,mfitby))
      return
      end subroutine

      real*8 pure function dbmag2(da,db,a,b,a1,b1) result(v)
      implicit none
      real*8 ,intent(in):: da,db,a,b,a1,b1
      real*8 dg,g,g1
      dg=2.d0*a*da-db
      g =(1.d0+ a**2)/b
      g1=(1.d0+a1**2)/b1
      v=.5d0*(b1*g*dg+g1*b*db)-a1*g*b*da
      return
      end function

      subroutine dgammag(dg,dbmag,da,db,a,b,a1,b1)
      implicit none
      real*8 ,intent(out):: dg,dbmag
      real*8 ,intent(in):: da,db,a,b,a1,b1
      real*8 g,g1
      dg=2.d0*a*da-db
      g =(1.d0+ a**2)/b
      g1=(1.d0+a1**2)/b1
      dbmag=.5d0*(b1*g*dg+g1*b*db)-a1*g*b*da
      return
      end

      subroutine tffsqu(nqcol,nqcol1,nvar,nqumax,ifquw,ifqu,
     $     free,nlat,nele,nfam,nfam1,nut,
     $     nparallel,cell,lfno,irtc)
      use tfstk
      use ffs, only:flv,ffs_bound,nvevx,nelvx,tsetintm
      use ffs_pointer
      use tffitcode
      use tfshare
      use mackw
      implicit none
      type (ffs_bound) fbound
      integer*8 ,intent(inout):: ifquw,ifqu
      integer*4 ,intent(in):: nqcol,nqcol1,nvar,nlat,nele,
     $     lfno,nut,nfam,nfam1,nparallel
      integer*4 ,intent(out):: irtc,nqumax
      logical*4 ,intent(in):: free(nele),cell
      logical*4 , allocatable,dimension(:,:)::col
      integer*8 itmmapp,kcm,kkqu,kqu,ic,iec
      integer*4 npp,nqu,k,kk,i,kq,j,kf,lp,idp,iv,kkf,kkq,
     $     ii,ltyp,jj,kc,ik1,iclast(-nfam:nfam),ik,nk,kk1,
     $     ip,ipr,istep,npr(nparallel),ll,ll1
      real*8 s,dtwiss(mfittry),coup,posk,wk,ctrans(4,7,-nfam:nfam)
      logical*4 disp,nzcod
      integer*4 , parameter :: minnqu=512
      call tffscoupmatrix(kcm,lfno)
      irtc=0
      allocate(col(2,nqcol))
      nqu=max(minnqu,nqcol*nvar)
      do9000: do
        if(nqu > nqumax)then
          if(nqumax > 0)then
            call tfree(ifquw)
            call tmunmapp(ifqu)
          endif
          ifqu=itmmapp(nqu)
          if(ifqu <= 0)then
            exit
          endif
          ifquw=ktaloc(nqu)
          if(ifquw <= 0)then
            call tmunmapp(ifqu)
            exit
          endif
          nqumax=nqu
        endif
        call tffsbound(fbound)
        rlist(ifqu:ifqu+nqu-1)=0.d0
        npp=min(nvar,nparallel)
        ipr=-1
        if(npp > 1)then
          istep=npp
          ip=0
          do while(ipr /= 0 .and. ip < npp-1)
            ip=ip+1
            ipr=itffork()
            npr(ip)=ipr
          enddo
          if(ipr /= 0)then
            ip=ip+1
          else
            call tsetintm(-1.d0)
          endif
        else
          ip=1
          istep=1
        endif  
        do k=1,nlat-1
          kc=icomp(k)
          kk=iele1(kc)
          kk1=iele1(k)
          iec=kele2(k)
          nk=merge(0,ilist(1,iec),iec == 0)
          if(kk > 0 .and. free(kk) .or. iec /= 0)then
            posk=pos(k)
            wk=1.d0
            if(k == fbound%lb)then
              wk=1.d0-fbound%fb
            endif
            if(k == fbound%le)then
              wk=fbound%fe
            endif
            ltyp=idtypec(k)
            do ii=ip,nvar,istep
              if(kk <= 0 .or.
     $             (.not. free(kk) .and. .not. free(kk1)))then
                ik1=1
              else
                ik1=0
                do ik=1,nk
                  if(nvevx(ii)%ivvar == ilist(2,iec+ik))then
                    ik1=1
                    exit
                  endif
                enddo
              endif
              do11: do ik=ik1,nk
                if(ik /= 0)then
                  if(kcm == 0)then
                    cycle
                  else
                    ic=kcm+(ilist(1,iec+ik)-1)*nvar-1
                    coup=rlist(ic+ii)*wk
                    if(coup /= 0.d0)then
                      iv=ilist(2,iec+ik)
                    else
                      cycle
                    endif
                  endif
                else
                  iv=nvevx(ii)%ivvar
                  if(iv == nelvx(kk)%ival .and.
     $                 nvevx(ii)%ivarele == kk .and.
     $                 (nvevx(ii)%ivcomp == 0 .or.
     $                 nvevx(ii)%ivcomp == kc))then
                    coup=couple(k)*wk
                  elseif(iv /= nelvx(kk)%ival .and.
     $                   nvevx(ii)%ivarele == kk1 .and.
     $                   (nvevx(ii)%ivcomp == 0 .or.
     $                   nvevx(ii)%ivcomp == k))then
                    coup=wk
                  else
                    cycle
                  endif
                endif
                iclast(nfam1:nfam)=0
                do concurrent (kq=1:nqcol1)
                  col(1,kq)=.true.
                  col(2,kq)=flv%lfp(2,kq) > 0
                enddo
                do kq=1,nqcol1
                  idp=flv%kdp(kq)
                  kqu=(ii-1)*nqcol+kq+ifqu-1
                  s=coup
                  do i=1,2
                    if(col(i,kq))then
                      j=flv%iqcol(kq)
                      kf=flv%kfit(flv%kfitp(j))
                      lp=flv%lfp(i,kq)
                      select case (kf)
                      case (mfitleng)
                        if(kytbl(kwL,ltyp) == iv)then
                          if(posk < pos(lp) .and. posk >= 0.d0)then
                            rlist(kqu)=rlist(kqu)+s
                          elseif(posk >= pos(lp) .and. posk < 0.d0)then
                            rlist(kqu)=rlist(kqu)-s
                          endif
                        endif
                      case (mfitbmagx,mfitbmagy,mfitbmagz)
c                        write(*,'(a,l2,5i5,1p10g12.4)')'ffsqu-bmag ',cell,i,kq,kf,k,lp,flv%lfp(3-i,kq)
                        call qdcell(dtwiss,k,lp,idp,iv,ctrans,iclast,
     $                       nfam,nut,.false.,.false.)
                        if(flv%lfp(3-i,kq) == 0)then
                          select case (kf)
                          case (mfitbmagx)
                            rlist(kqu)=coup*dtwiss(mfitbmagx)
                          case (mfitbmagy)
                            rlist(kqu)=coup*dtwiss(mfitbmagy)
                          case (mfitbmagz)
                            rlist(kqu)=coup*dtwiss(mfitbmagz)
                          end select
                        elseif(cell .or. k < lp)then
                          ll=itwissp(lp)
                          ll1=itwissp(flv%lfp(3-i,kq))
                          select case (kf)
                          case (mfitbmagx)
                            rlist(kqu)=rlist(kqu)+coup*
     $                           dbmag2(dtwiss(mfitax),dtwiss(mfitbx),
     $                           utwiss(mfitax,idp,ll),utwiss(mfitbx,idp,ll),
     $                           utwiss(mfitax,idp,ll1),utwiss(mfitbx,idp,ll1))
                          case (mfitbmagy)
                            rlist(kqu)=rlist(kqu)+coup*
     $                           dbmag2(dtwiss(mfitay),dtwiss(mfitby),
     $                           utwiss(mfitay,idp,ll),utwiss(mfitby,idp,ll),
     $                           utwiss(mfitay,idp,ll1),utwiss(mfitby,idp,ll1))
                          case (mfitbmagz)
                            rlist(kqu)=rlist(kqu)+coup*
     $                           dbmag2(dtwiss(mfitaz),dtwiss(mfitbz),
     $                           utwiss(mfitaz,idp,ll),utwiss(mfitbz,idp,ll),
     $                           utwiss(mfitaz,idp,ll1),utwiss(mfitbz,idp,ll1))
                          end select
                        endif
                      case default
                        if(kf <= mfittry)then
                          if(k < lp .or. posk < pos(lp) .or. cell)then
                            nzcod= kf == mfitdz.or. kf == mfitddp
                            disp=kf >= mfitex .and. kf <= mfitepy .or. nzcod
                            if(.not. disp)then
                              do kkq=kq+1,nqcol1
                                jj=flv%iqcol(kkq)
                                kkf=flv%kfit(flv%kfitp(jj))
                                disp=kkf >= mfitex .and. kkf <= mfitepy
                                if(disp)then
                                  exit
                                endif
                              enddo
                            endif
                            call qdcell(dtwiss,
     $                           k,lp,idp,iv,ctrans,iclast,nfam,nut,disp,nzcod)
                            rlist(kqu)=rlist(kqu)+s*dtwiss(kf)
                            do kkq=kq+1,nqcol1
                              jj=flv%iqcol(kkq)
                              kkf=flv%kfit(flv%kfitp(jj))
                              if(kkf <= mfito .and. idp == flv%kdp(kkq))then
                                kkqu=(ii-1)*nqcol+kkq+ifqu-1
                                if(col(1,kkq) .and. flv%lfp(1,kkq) == lp)then
                                  rlist(kkqu)=rlist(kkqu)+dtwiss(kkf)*coup
                                  col(1,kkq)=.false.
                                elseif(col(2,kkq) .and. flv%lfp(2,kkq) == lp)then
                                  rlist(kkqu)=rlist(kkqu)-dtwiss(kkf)*coup
                                  col(2,kkq)=.false.
                                endif
                              endif
                            enddo
                          endif
                        else
                          call tdgeo(s,rlist(kqu),kf,lp,k,ltyp,iv,nut,nfam)
                        endif
                      end select
                    endif
                    s=-s
                  enddo
                enddo
              enddo do11
            enddo
          endif
        enddo
        call tffswait(ipr,npp,npr,i00,0,'tffsqu',irtc)
        if(kcm /= 0)then
          call tfree(kcm)
        endif
        return
      enddo do9000
      call termes(lfno,'?Too many conditions*variables.',' ')
      if(nqumax > 0)then
        ifqu=itmmapp(nqumax)
        ifquw=ktaloc(nqumax)
      endif
      irtc=20002
      if(kcm /= 0)then
        call tfree(kcm)
      endif
      return
      end subroutine

      end module
