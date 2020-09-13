      subroutine tdfun(iqcol,lfp,nqcola,nqcola1,kdp,df1,error)
      use tfstk
      use ffs
      use ffs_pointer
      use ffs_fit
      use tffitcode
      implicit none
      real*8 ,parameter :: abmax=1.d-8
      real*8 ,parameter::factor=0.97d0,dmax=1.d10
      integer*4 ,parameter::npeak=10
      integer*4 ,intent(out):: nqcola,nqcola1,
     $     iqcol(*),lfp(2,maxcond),kdp(*)
      real*8 ,intent(out):: df1(*)
      real*8 vpeak(npeak),vf,tdfun1,tgfun,vb,ve,v,vf1
      logical*4 ,intent(out):: error
      integer*4 ipeak(npeak),j,i,ka,kf,kp,kp1,mp,idp,kpb,kpe,
     $     k,ip,irtc,m
      logical*4 maxfit,ttrans(-nfam:nfam),tftype1fit
c      do j=1,nfcol
c        if(flv%kfit(flv%kfitp(j)) .eq. mfitnx .or.
c     $       flv%kfit(flv%kfitp(j)) .eq. mfitny)then
c          do m=nfam1,nfam
c            detr0=utwiss(mfitr1,m,1)*utwiss(mfitr4,m,1)
c     $           -utwiss(mfitr2,m,1)*utwiss(mfitr3,m,1)
c            ttrans(m)=.false.
c            if(detr0 .lt. 1.d0)then
c              do i=2,nut
c                if(utwiss(mfitr1,m,i)*utwiss(mfitr4,m,i)
c     $               -utwiss(mfitr2,m,i)*utwiss(mfitr3,m,i)
c     1               .ge. 1.d0)then
c                  ttrans(m)=.true.
c                  exit
c                endif
c              enddo
c            else
c              do i=2,nut
c                if(utwiss(mfitr1,m,i)*utwiss(mfitr4,m,i)
c     $               -utwiss(mfitr2,m,i)*utwiss(mfitr3,m,i)
c     1               .lt. 1.d0)then
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
c        write(*,*)'TDFUN ',j,ka,kf,mfitgx
        if(kf .gt. mfitchi3)then
          cycle
        endif
        kp=flv%ifitp(ka)
        kp1=flv%ifitp1(ka)
        if(flv%mfitp(ka) .ne. 0)then
          mp=(abs(flv%mfitp(ka))-1)/2
          do20:     do idp=nfam1,nfam
            vf=flv%fitval(ka)
            if((kfam(idp) .eq. 0 .and. idp .ge. -mp .and. idp .le. mp)
     $           .or. jfam(idp) .ge. -mp .and. jfam(idp) .le. mp)then
              if(kf .ge. mfitdx .and. kf .le. mfitdpy .and.
     $             (idp .lt. -nfr .or. idp .gt. nfr ))then
c     $             .or. inicond .and. idp .ne. 0))then
                cycle
              elseif(idp .ne. 0 .and. kf .gt. mfittry)then
                cycle
              endif
              if(idp .ne. 0 .or. mp*2+1 .ne. abs(flv%mfitp(ka)))then
                if(kp .eq. nlat)then
                  if(idp .ge. -1 .and. idp .le. 1)then
                    if(kf .eq. mfitnx)then
                      vf=vf+(dp(idp)-dp0)*xixf
                    endif
                    if(kf .eq. mfitny)then
                      vf=vf+(dp(idp)-dp0)*xiyf
                    endif
                  endif
                endif
                maxfit=flv%mfitp(ka) .lt. 0
                if(kp .ne. kp1)then
                  maxfit=maxfit .and. .not. tftype1fit(kf)
                  kpb=min(kp,kp1)
                  kpe=max(kp,kp1)
                  if(maxfit)then
                    if(kf .le. mfitdpy
     $                   .or. kf .ge. mfitpex .and.
     $                   kf .le. mfitgmz)then
                      call tfpeak(idp,kf,kpb,kpe,ipeak,vpeak,npeak)
                      do k=1,npeak
                        ip=ipeak(k)
                        if(ip .le. 0)then
                          cycle do20
                        endif
                        df1(i)=tdfun1(vf,vpeak(k),
     $                       kf,maxfit,idp,ttrans(idp))
                        if(df1(i) .ne. 0.d0)then
                          iqcol(i)=j
                          lfp(1,i)=ip
                          lfp(2,i)=0
                          kdp(i)=idp
                          i=i+1
                          if(i .gt. maxcond)then
                            error=.true.
                            return
                          endif
                        endif
                      enddo
                    endif
                  else
                    if(kf .le. mfitgmz .or.
     $                   (kf .ge. mfitleng .and. kf .le. mfitchi3))then
                      maxfit=flv%mfitp(ka) .lt. 0
                      vb=tgfun(kf,kpb,idp)
                      ve=tgfun(kf,kpe,idp)
                      vf1=vb
                      call tfgetfitval(nlist(kf),
     $                     kpb,kpe,dp(idp),
     $                     iuid(idp),kfam(idp),
     $                     vb,ve,vf,vf1,irtc)
                      if(irtc .eq. -1)then
                        cycle
                      endif
                      if(maxfit)then
                        ve=ve-vf1
                        vf1=vf
                      endif
                      df1(i)=tdfun1(vf1,ve,kf,maxfit,idp,ttrans(idp))
                      if(cell .and. ka .gt. nfc0 .and. ka .le. nfc0+4
     $                   .and. abs(df1(i)) .lt. abmax)then
                        cycle
                      endif
                      if(.not. maxfit)then
                        df1(i)=df1(i)+merge(log(vf),vf,
     $                       kf .eq. mfitbx .or. kf .eq. mfitby
     $                       .or. kf .eq. mfitbz)
                      endif
                      iqcol(i)=j
                      lfp(1,i)=kpe
                      lfp(2,i)=kpb
                      kdp(i)=idp
                      i=i+1
                      if(i .gt. maxcond)then
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
                  if(irtc .eq. -1)then
                    cycle
                  endif
                  df1(i)=tdfun1(vf1,v,kf,maxfit,idp,ttrans(idp))
                  if(maxfit .and. df1(i) .eq. 0.d0)then
                    cycle
                  endif
                  iqcol(i)=j
                  lfp(1,i)=kp
                  lfp(2,i)=0
                  kdp(i)=idp
                  i=i+1
                  if(i .gt. maxcond)then
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

      subroutine tfgetfitval(funname,kp,kp1,dp,
     $     iuid,kfam,vf,v,vf0,vf1,irtc)
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
      if(kfv%k .eq. 0)then
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
      if(kp1 .eq. 0)then
        retry1=.false.
        klv%rbody(4)=vf
        klv%rbody(5)=v
      else
        call elname(kp1,name1)
        ln1=lenw(name1)
        retry1=kp1 .ne. nlat
        klv1%rbody(5)=vf0
        klv1%rbody(6)=vf
        klv1%rbody(7)=v
      endif
      retry=kp .ne. nlat
 100  call tfpadstr(name,ifvloc+1,ln)
      ilist(1,ifvloc)=ln
      if(kp1 .ne. 0)then
        call tfpadstr(name1,ifvloc1+1,ln1)
        ilist(1,ifvloc1)=ln1
      endif
      call tclrfpe
      levele=itfuplevel()
      if(kp1 .eq. 0)then
c        call tfdebugprint(kfv,'FitValue-1',1)
        kx=tfleval(klv,.true.,irtc)
      else
c        call tfdebugprint(kfv1,'FitValue-2',1)
        kx=tfleval(klv1,.true.,irtc)
      endif
c      call tfdebugprint(kx,'==> ',1)
      level=itfdownlevel()
 110  if(irtc .ne. 0)then
        if(ierrorprint .ne. 0)then
          call tfaddmessage(' ',0,6)
        endif
        call termes(6,
     $       'Error in FitValue '//
     $       funname//' at '//name,' ')
      elseif(ktfrealq(kx,vf1))then
      elseif(kx%k .eq. ktfoper+mtfnull)then
        irtc=-1
      elseif(retry)then
        retry=.false.
        if(mult(kp) .eq. 0)then
c     Generate singlet element name with suffix number(.###)
          call elnameK(kp,name)
          ln=lenw(name)
          go to 100
        elseif(nelvx(ilist(kp,ifele1))%klp .eq. kp)then
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
        retry=kp .ne. nlat
        if(retry)then
          call elname(kp,name)
          ln=lenw(name)
        endif
        if(mult(kp1) .eq. 0)then
c     Generate singlet element name with suffix number(.###)
          call elnameK(kp1,name1)
          ln1=lenw(name1)
          go to 100
        elseif(nelvx(ilist(kp1,ifele1))%klp .eq. kp1)then
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
      if(kff%k .eq. 0)then
        kff=kxsymbolz('`FitFunction',12)
      endif
      level=itfuplevel()
      kx%k=0
      kx=tfsyeval(kff,irtc)
c      call tfdebugprint(kx,'fitfun',3)
      if(irtc .ne. 0)then
        level=itfdownlevel()
        call tfaddmessage(' ',2,icslfno())
        call termes(6,'Error in FitFunction ',' ')
        error=.true.
        return
      elseif(ktfrealq(kx,df(nqcol+1)))then
        if(nqcol .ge. maxcond)then
          error=.true.
          level=itfdownlevel()
          return
        endif
        nqcol=nqcol+1
        iqcol(nqcol)=-1
        kdp(nqcol)=0
      elseif(tfreallistq(kx%k,klx))then
        m=klx%nl
        if(m .gt. 0)then
          if(m+nqcol .gt. maxcond)then
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
          if(v .gt. vf)then
            tdfun1=log(vf*factor/v)
          else
            tdfun1=0.d0
          endif
          return
        else
          tdfun1=log(vf/v)
        endif
        return

      case (mfitax,mfitay,mfitaz)
        if(maxfit)then
          vfa=abs(vf)
          if(v .gt. vfa)then
            tdfun1=atan(vfa*factor)-atan(v)
          elseif(v .lt. -vfa)then
            tdfun1=-atan(vfa*factor)-atan(v)
          else
            tdfun1=0.d0
          endif
          return
        else
          tdfun1=atan(vf)-atan(v)
        endif
        return

      case (mfitex,mfitey)
        if(maxfit)then
          vfa=abs(vf)
          if(v .gt. vfa)then
            tdfun1=vfa-v
          else
            tdfun1=max(-vfa-v,0.d0)
          endif
          return
        else
          tdfun1=merge(-vf,vf,kdp .lt. 0)-v
        endif
        return

      case (mfitchi1,mfitchi2,mfitchi3)
        if(maxfit)then
          vfa=abs(vf)
          if(v .gt. vfa)then
            tdfun1=vfa-v
          else
            tdfun1=max(-vfa-v,0.d0)
          endif
          return
        else
          tdfun1=vf-v
        endif
        do while(tdfun1 .lt. -pi)
          tdfun1=tdfun1+pi2
        enddo
        do while(tdfun1 .gt. pi)
          tdfun1=tdfun1-pi2
        enddo
        return
        
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
        return

      case default
        if(maxfit)then
          vfa=abs(vf)
          if(v .gt. vfa)then
            tdfun1=vfa-v
          else
            tdfun1=max(-vfa-v,0.d0)
          endif
          return
        else
          tdfun1=vf-v
        endif
c        write(*,*)'tdfun1 ',kf,maxfit,vf,v,tdfun1
        return
      end select
      return

      end

      subroutine tfpeak(idp,kf,ibegin,iend,ipeak,vpeak,npeak)
      use ffs_pointer
      use tffitcode
      implicit none
      integer*4 ,intent(in)::ibegin,iend,kf,npeak,idp
      integer*4 ,intent(out)::ipeak(npeak)
      real*8 ,intent(out)::vpeak(npeak)
      integer*4 i,j,k
      real*8 va,va0,va1,tgfun
      vpeak=0.d0
      ipeak=0
      va0=0.d0
      va=abs(tgfun(kf,ibegin,idp))
      do i=ibegin,iend
        va1=abs(tgfun(kf,min(i+1,iend),idp))
c        if(kf .eq. mfitgmy)then
c          write(*,*)'tfpeak ',i,va,va0,va1
c        endif
        if(va .gt. va0 .and. va .ge. va1)then
          do j=1,npeak
            if(va .gt. abs(vpeak(j)))then
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
