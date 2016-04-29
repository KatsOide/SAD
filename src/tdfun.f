      subroutine tdfun(flv,
     1     nfcol,iqcol,lfp,nqcol,nqcol1,
     1     kdp,nfr,nfam,nfam1,jfam,kfam,
     $     inicond,iuid,
     $     df,utwiss,itwissp,nut,
     1     tracex,tracey,pos,geo,
     $     latt,mult,dp,maxcond,
     $     nlist,error)
      use tfstk
      use ffs
      use ffslocal, only: ffslocalv
      use tffitcode
      implicit none
      type (ffslocalv) flv
      real*8 factor,dmax
      integer*4 npeak
      parameter (factor=0.97d0,dmax=1.d10)
      parameter (npeak=10)
      integer*4 nqcol,nut,nfr,itwissp(nlat),maxcond,
     $     nfam,nqcol1,nfam1,j,nfcol,i,ka,kf,kp,kp1,mp,
     $     idp,kpb,kpe,k,ip,irtc,m
      integer*4 
     1     kdp(*),iqcol(*),
     1     lfp(2,maxcond),mult(*),latt(2,nlat),ipeak(npeak),
     $     iuid(-nfam:nfam)
      real*8 utwiss(ntwissfun,-nfam:nfam,nut),
     $     df(*),dp(-nfam:nfam),pos(nlat),geo(3,4,nlat),
     $     tracex(-nfam:nfam),tracey(-nfam:nfam),
     $     vpeak(npeak),vf,tdfun1,tgfun,vb,ve,v,vf1
      logical*4 error,maxfit,ttrans(-nfam:nfam),tftype1fit,inicond
      character*8 nlist(*)
      integer*4 jfam(-nfam:nfam),kfam(-nfam:nfam)
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
      do 10 j=1,nfcol
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
     $             (idp .lt. -nfr .or. idp .gt. nfr .or.
     $             inicond .and. idp .ne. 0))then
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
c           write(*,*)'TDFUN ',j,i,idp,kf,kp,kp1,maxcond,flv%mfitp(ka)
                if(kp .ne. kp1)then
                  kpb=min(kp,kp1)
                  kpe=max(kp,kp1)
 3010             if(maxfit)then
                    if(tftype1fit(kf))then
                      maxfit=.false.
                      goto 3010
                    elseif(kf .le. mfitdpy
     $                     .or. kf .ge. mfitpex .and.
     $                     kf .le. mfitpepx)then
                      call tfpeak(utwiss,itwissp,idp,kf,kpb,kpe,
     $                     ipeak,vpeak,npeak,nlat,nfam,ntwissfun)
                      do 110 k=1,npeak
                        ip=ipeak(k)
                        if(ip .le. 0)then
                          cycle do20
                        endif
                        df(i)=tdfun1(vf,vpeak(k),
     $                       kf,maxfit,idp,ttrans(idp))
c                        write(*,*)'tdfun ',k,ip,kf,maxfit,vpeak(k),df(i)
                        if(df(i) .ne. 0.d0)then
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
 110                  continue
                    endif
                  else
                    if(kf .le. mfitpepy .or.
     $                   (kf .ge. mfitleng .and. kf .le. mfitgz))then
                      maxfit=flv%mfitp(ka) .lt. 0
                      vb=tgfun(kf,kpb,idp,utwiss,itwissp,
     $                     tracex,tracey,pos,geo,nfam)
                      ve=tgfun(kf,kpe,idp,utwiss,itwissp,
     $                     tracex,tracey,pos,geo,nfam)
                      vf1=vb
                      call tfgetfitval(latt,mult,nlist(kf),
     $                     kpb,kpe,dp(idp),
     $                     iuid(idp),kfam(idp),inicond,
     $                     vb,ve,vf,vf1,irtc)
                      if(irtc .eq. -1)then
c                        write(*,*)kf,vb,ve,vf1
                        cycle
                      endif
                      if(maxfit)then
                        ve=ve-vf1
                        vf1=vf
                      endif
                      df(i)=tdfun1(vf1,ve,kf,maxfit,idp,ttrans(idp))
                      if(tftype1fit(kf))then
                        if(.not. maxfit)then
                          df(i)=vf+df(i)
                        endif
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
                  v=tgfun(kf,kp,idp,utwiss,itwissp,
     $                 tracex,tracey,pos,geo,nfam)
                  vf1=vf
                  call tfgetfitval(latt,mult,nlist(kf),kp,0,dp(idp),
     $                 iuid(idp),kfam(idp),inicond,vf,v,vf,vf1,irtc)
                  if(irtc .eq. -1)then
                    cycle
                  endif
                  df(i)=tdfun1(vf1,v,kf,maxfit,idp,ttrans(idp))
c                  write(*,*)'tdfun ',vf1,v,df(i)
                  if(maxfit .and. df(i) .eq. 0.d0)then
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
 10   continue
      nqcol=i-1
      nqcol1=nqcol
      call tffsfitfun(nqcol,df,iqcol,kdp,maxcond,error)
      return
      end

      subroutine tfgetfitval(latt,mult,nlist,kp,kp1,dp,
     $     iuid,kfam,inicond,vf,v,vf0,vf1,irtc)
      use tfstk
      use ffs
      use tffitcode
      implicit none
      type (sad_descriptor) kx
      integer*4 latt(2,*),mult(*)
      character*8 nlist
      integer*4 kp,kp1
      real*8 dp
      integer*4 iuid,kfam
      logical*4 inicond
      real*8 vf,v,vf0,vf1
      integer*4 irtc
c     
      integer*4 lenw,itfuplevel,itfdownlevel
c
      logical retry,retry1
      character*(MAXPNAME+8) name,name1
      integer*4 level,ln,ln1
      integer*8 ifv,ifvh,ifvloc,ifvfun,ifid,ifv1,ifvloc1
      save ifv,ifvh,ifvloc,ifvfun,ifid,ifv1,ifvloc1
      data ifv /0/
      if(ifv .eq. 0)then
        ifv =ktadaloc(0,5)
        ifv1=ktadaloc(0,7)
        ifid=ktavaloc(0,2)
        ifvh=ktfsymbolz('`FitValue',9)
        ifvloc=ktsalocb(0,'                ',MAXPNAME+8)
        ifvloc1=ktsalocb(0,'                ',MAXPNAME+8)
        ifvfun=ktsalocb(0,'        ',MAXPNAME)
        klist(ifv)=ktfsymbol+ktfcopy1(ifvh)
        klist(ifv+1)=ktfstring+ifvloc
        klist(ifv+2)=ktfstring+ifvfun
        klist(ifv+3)=ktflist+ifid
        klist(ifv1)=ktfsymbol+ktfcopy1(ifvh)
        klist(ifv1+1)=ktfstring+ktfcopy1(ifvloc)
        klist(ifv1+2)=ktfstring+ktfcopy1(ifvloc1)
        klist(ifv1+3)=ktfstring+ktfcopy1(ifvfun)
        klist(ifv1+4)=ktflist+ktfcopy1(ifid)
      endif
      irtc=0
      call tfpadstr(nlist,ifvfun+1,len_trim(nlist))
      ilist(1,ifvfun)=len_trim(nlist)
      if(inicond)then
        rlist(ifid+1)=dble(iuid)
      else
        rlist(ifid+1)=dble(kfam)
      endif
      rlist(ifid+2)=dp
      call elname1(latt,kp,mult,name,.false.)
      ln=lenw(name)
      if(kp1 .eq. 0)then
        retry1=.false.
        rlist(ifv+4)=vf
        rlist(ifv+5)=v
      else
        call elname1(latt,kp1,mult,name1,.false.)
        ln1=lenw(name1)
        retry1=kp1 .ne. nlat
        rlist(ifv1+5)=vf0
        rlist(ifv1+6)=vf
        rlist(ifv1+7)=v
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
c        call tfdebugprint(ktflist+ifv,'FitValue-1',3)
        call tfleval(klist(ifv-3),kx,.true.,irtc)
      else
c        call tfdebugprint(ktflist+ifv1,'FitValue-2',3)
        call tfleval(klist(ifv1-3),kx,.true.,irtc)
      endif
c      call tfdebugprint(kx,'==> ',3)
c      write(*,*)'kp: ',kp,'kp1: ',kp1
      level=itfdownlevel()
 110  if(irtc .ne. 0)then
        if(ierrorprint .ne. 0)then
          call tfaddmessage(' ',0,6)
        endif
        call termes(6,
     $       'Error in FitValue '//
     $       nlist//' at '//name,' ')
      elseif(ktfrealqd(kx,vf1))then
      elseif(kx%k .eq. ktfoper+mtfnull)then
        irtc=-1
      elseif(retry)then
        retry=.false.
        if(mult(kp) .eq. 0)then
c     Generate singlet element name with suffix number(.###)
          call elname1(latt,kp,mult,name,.true.)
          ln=lenw(name)
          go to 100
        elseif(ilist(ilist(kp,ifele1),ifklp) .eq. kp)then
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
      elseif(retry1)then
        retry1=.false.
c     Reset `kp'-element name in name(1:ln)
        retry=kp .ne. nlat
        if(retry)then
          call elname1(latt,kp,mult,name,.false.)
          ln=lenw(name)
        endif
c
        if(mult(kp1) .eq. 0)then
c     Generate singlet element name with suffix number(.###)
          call elname1(latt,kp1,mult,name1,.true.)
          ln1=lenw(name1)
          go to 100
        elseif(ilist(ilist(kp1,ifele1),ifklp) .eq. kp1)then
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
      implicit none
      type (sad_list), pointer :: klx
      type (sad_descriptor) kx
      integer*8 kff
      integer*4 maxcond,nqcol,iqcol(maxcond),kdp(maxcond)
      real*8 df(maxcond)
      integer*4 l,itfuplevel,itfdownlevel,i,m,level,irtc,icslfno
      logical*4 error
      save kff
      data kff/0/
      if(kff .eq. 0)then
        kff=ktfcopy(ktfsymbolz('`FitFunction',12))
      endif
      l=itfuplevel()
      kx%k=0
      call tfsyeval(kff,kx,irtc)
c      call tfdebugprint(kx,'fitfun',3)
c      write(*,*)'with ',irtc
      if(irtc .ne. 0)then
        level=itfdownlevel()
        call tfaddmessage(' ',2,icslfno())
        call termes(6,'Error in FitFunction ',' ')
        error=.true.
        return
      elseif(ktfrealqd(kx,df(nqcol+1)))then
        if(nqcol .ge. maxcond)then
          error=.true.
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
      level=itfdownlevel()
      error=.false.
      return
      end

      real*8 function tdfun1(vf,v,kf,maxfit,kdp,ttrans)
      implicit none
      include 'inc/MACMATH.inc'
      real*8 factor
      parameter (factor=1.d0)
      integer*4 kf,kdp
      real*8 vf,v,vfa
      logical*4 maxfit,ttrans
      go to (
     $     1130,1110,1210,1130,1110,1210,1170,1190,1170,1190,
     $     1190,1190,1190,1190,1190,1190,1190,1190,1190,1190,
     $     1190,1190,1190,1190,1190,1190,1190,1190,1190,1190,
     $     1190,1410,1410,1410
     $     ),kf
      go to 1190
c     $  'AX   ','BX   ','NX   ','AY   ','BY   ','NY   ','EX   ','EPX  ','EY   ','EPY  ',
c     $  'R1   ','R2   ','R3   ','R4   ','DETR ','DX   ','DPX  ','DY   ','DPY  ','DZ   ',
c     $  'DDP  ','PEX  ','PEPX ','PEY  ','PEPY ','TRX  ','TRY  ','LENG ','GX   ','GY   ',
c     $  'GZ   ','CHI1 ','CHI2 ','CHI3 ','DEX  ','DEPX ','DEY  ','DEPY ','DDX  ','DDPX ',
c     $  'DDY  ','DDPY ','PDEX ','PDEPX','PDEY ','PDEPY'/
1110  if(maxfit)then
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
1130  if(maxfit)then
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
1170  if(maxfit)then
        vfa=abs(vf)
        if(v .gt. vfa)then
          tdfun1=vfa-v
        elseif(v .lt. -vfa)then
          tdfun1=-vfa-v
        else
          tdfun1=0.d0
        endif
        return
      else
        if(kdp .lt. 0)then
          tdfun1=-vf-v
        else
          tdfun1=vf-v
        endif
      endif
      return
1190  if(maxfit)then
        vfa=abs(vf)
        if(v .gt. vfa)then
          tdfun1=vfa-v
        elseif(v .lt. -vfa)then
          tdfun1=-vfa-v
        else
          tdfun1=0.d0
        endif
        return
      else
        tdfun1=vf-v
      endif
      return
 1210 continue
c      vf1=pi2*(anint(vf/pi2)+sign(.5d0*sin(.5d0*vf)**2,sin(vf)))
c      v1=pi2*(anint(v/pi2)+sign(.5d0*sin(.5d0*v)**2,sin(v)))
      if(maxfit)then
        if(v .gt. vf)then
          tdfun1=vf-v
        else
          tdfun1=0.d0
        endif
      else
        tdfun1=vf-v
      endif
      if(ttrans)then
        tdfun1=tdfun1-anint(tdfun1/pi2)*pi2
      endif
      return
 1410 if(maxfit)then
        vfa=abs(vf)
        if(v .gt. vfa)then
          tdfun1=vfa-v
        elseif(v .lt. -vfa)then
          tdfun1=-vfa-v
        else
          tdfun1=0.d0
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
      end

      subroutine tfpeak(utwiss,itwissp,idp,kf,
     $ibegin,iend,ipeak,vpeak,npeak,nlat,nfam,ntwissfun)
      implicit none
      integer*4 nlat,nfam,ibegin,iend,kf,npeak,idp,ipeak(npeak),
     $     itwissp(nlat),ntwissfun
      real*8 vpeak(npeak),utwiss(ntwissfun,-nfam:nfam,*)
      integer*4 i,j,k
      real*8 va,va0,va1
      do 10 i=1,npeak
        vpeak(i)=0.d0
        ipeak(i)=0
10    continue
      va0=0.d0
      va=abs(utwiss(kf,idp,itwissp(ibegin)))
      do 110 i=ibegin,iend
        va1=abs(utwiss(kf,idp,itwissp(min(i+1,iend))))
        if(va .gt. va0 .and. va .ge. va1)then
          do 120 j=1,npeak
            if(va .gt. abs(vpeak(j)))then
              do 130 k=npeak,j+1,-1
                vpeak(k)=vpeak(k-1)
                ipeak(k)=ipeak(k-1)
 130          continue
              vpeak(j)=utwiss(kf,idp,itwissp(i))
              ipeak(j)=i
c              write(*,*)'tfpeak ',j,i,kf,idp,itwissp(i),vpeak(j)
              go to 111
            endif
 120      continue
        endif
 111    va0=va
        va=va1
 110  continue
      return
      end
