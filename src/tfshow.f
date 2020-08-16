      subroutine tfshow(stab,df,mfpnta,mfpnta1,
     $     kx,irtc,ret,lfno)
      use tfstk
      use ffs
      use ffs_pointer
      use ffs_fit
      use tffitcode
      implicit none
      integer*8 kx,kax,kax1,kax2,kax3,kax4,kax3i,kaxi,kaxi4
      integer*4 namel,mfpnta,mfpnta1,
     $     lfno,lw,lf,lfs,nn,mf,m,i,k,l,jshowi,
     $     ncc,iq,kpa,kpb,jm,lfs1,ln,
     $     mm,lb,ip,j,nf,nl
      parameter (namel=11)
      integer*4 icalc1(3,ndim1)
      real*8 df(nqcol),fm,tgfun
      integer*4 jshow(64),itfgetrecl,irtc
      character*5 fun,form,forms
      character*31 name
      character*15 name1,namea
      character*10 autofg
      character*11 vout
      character*256 buf2
      character*256 buf0,buf1
      character*3 famlabel
      logical*4 err,stab,trx,try,ret,tftype1fit
      external trim
      lw=max(79,min(256,itfgetrecl()))
      lf=max(6,min(10,(lw-20)/(2*nfam+2)))
      if(lf .gt. 6)then
        lfs=max(6,min(8,lf,lw-lf*(2*nfam+1)-22))
        nn=(lw-22-lfs)/lf
      else
        lfs=max(6,min(8,lf,lw-lf*(2*nfam+1)-20))
        nn=(lw-20-lfs)/lf
      endif        
      write(form,'(i2,''.'',i1,1x)')lfs,min(6,lfs-3)
      call trim(form)
      forms='S'//form(1:len_trim(form))
      write(form,'(i2,''.'',i1,1x)')lf,min(6,lf-3)
      call trim(form)
      mf=nfam1
      if(nfam .gt. 0)then
        nn=min(nn,nfam-mf+1)
        m=0
        do i=mf,nfam
          do k=1,m
            if(residual(jshow(k)) .lt. residual(i))then
              do l=m,k,-1
                jshow(l+1)=jshow(l)
              enddo
              jshow(k)=i
              go to 301
            endif
          enddo
          jshow(m+1)=i
 301      m=min(m+1,nn)
        enddo
        do i=1,nn
          if(jshow(i) .eq. 0)then
            go to 302
          endif
        enddo
        jshow(nn)=0
 302    do i=2,nn
          jshowi=jshow(i)
          do j=1,i-1
            if(dp(jshow(j)) .gt. dp(jshow(i)))then
              do k=i-1,j,-1
                jshow(k+1)=jshow(k)
              enddo
              jshow(j)=jshowi
              go to 303
            endif
          enddo
 303      continue
        enddo
        do m=1,nn
          buf0((m-1)*lf+1:m*lf)=
     $         autofg(dp(jshow(m))+dp0,form)
        enddo
        vout(1:lfs+3)=' DP '
        if(lfs .gt. 6)then
          buf1='           '//vout(1:lfs+3)
     $         //'       '//buf0(1:nn*lf)
        else
          buf1='          '//vout(1:lfs+3)
     $         //'      '//buf0(1:nn*lf)
        endif
        write(lfno,'(A)')buf1(1:len_trim(buf1))
        if(nfam .gt. nfr)then
          do m=1,nn
            buf0((m-1)*lf+1:m*lf)='   '//famlabel(kfam(jshow(m)))
          enddo
          vout(1:lfs+3)=' '
          if(lfs .gt. 6)then
            buf1='                '//vout(1:lfs+3)
     $           //'  '//buf0(1:nn*lf)
          else
            buf1='               '//vout(1:lfs+3)
     $           //' '//buf0(1:nn*lf)
          endif
          write(lfno,'(A)')buf1(1:len_trim(buf1))
        elseif(inicond)then
          do m=1,nn
            write(buf0((m-1)*lf+1:m*lf),'(i6)')iuid(jshow(m))
          enddo
          vout(1:lfs+3)=' Orbit ID '
          if(lfs .gt. 6)then
            buf1='           '//vout(1:lfs+3)
     $           //'       '//buf0(1:nn*lf)
          else
            buf1='          '//vout(1:lfs+3)
     $           //'      '//buf0(1:nn*lf)
          endif
          write(lfno,'(A)')buf1(1:len_trim(buf1))
        endif
        do m=1,nn
          buf0((m-1)*lf+1:m*lf)=
     $         autofg(residual(jshow(m)),form)
        enddo
        vout(1:lfs+3)=' Res.'
        if(lfs .gt. 6)then
          buf1='           '//vout(1:lfs+3)
     $         //'       '//buf0(1:nn*lf)
        else
          buf1='          '//vout(1:lfs+3)
     $         //'      '//buf0(1:nn*lf)
        endif
        write(lfno,'(A)')buf1(1:len_trim(buf1))
        mm=nn
      else
        do j=mf,nfam
          jshow(j-mf+1)=j
        enddo
        mm=nfam-mf+1
      endif
      ncc=flv%ncalc
      trx=.false.
      try=.false.
      do 3000 i=1,ncc
        icalc1(1,i)=flv%icalc(1,i)
        icalc1(2,i)=flv%icalc(2,i)
        icalc1(3,i)=flv%icalc(3,i)
        if(flv%icalc(3,i) .eq. mfittrx)then
          trx=.true.
        endif
        if(flv%icalc(3,i) .eq. mfittry)then
          try=.true.
        endif
3000  continue
      if(.not. stab)then
        if(.not. trx)then
          ncc=ncc+1
          icalc1(1,ncc)=nlat
          icalc1(2,ncc)=nlat
          icalc1(3,ncc)=mfittrx
        endif
        if(.not. try)then
          ncc=ncc+1
          icalc1(1,ncc)=nlat
          icalc1(2,ncc)=nlat
          icalc1(3,ncc)=mfittry
        endif
      endif
      do i=nqcol1+1,nqcol
        ncc=ncc+1
        icalc1(1,ncc)=nlat
        icalc1(2,ncc)=nlat
        icalc1(3,ncc)=nqcol1-i
      enddo
      do 3010 i=1,nqcol1
        iq=flv%kfitp(flv%iqcol(i))
        if(flv%ifitp(iq) .ne. flv%ifitp1(iq))then
          if(flv%lfp(2,i) .eq. 0)then
            call txcalc(icalc1,ncc,flv%lfp(1,i),flv%lfp(1,i),
     $           flv%kfit(iq),.true.,err)
          endif
        elseif(flv%mfitp(iq) .lt. 0)then
          call txcalc(icalc1,ncc,flv%lfp(1,i),flv%lfp(1,i),flv%kfit(iq),
     $         .true.,err)
        endif
3010  continue
      kax4=0
      kaxi4=0
      if(ret)then
        nf=nfam-mf+1
        irtc=0
        kax =ktadaloc(-1,4)
        kax1=ktavaloc(0,nf)
        kax2=ktavaloc(0,nf)
        kax3=ktadaloc(0,nf)
        kax4=ktadaloc(0,ncc)
        do i=mf,nfam
          rlist(kax1+i-mf+1)=dp(i)+dp0
          rlist(kax2+i-mf+1)=dble(merge(iuid(i),kfam(i),inicond))
          kax3i=ktavaloc(0,3)
          klist(kax3+i-mf+1)=ktflist+kax3i
          rlist(kax3i+1)=residual(i)
          rlist(kax3i+2)=merge(1.d0,0.d0,optstat(i)%stabx)
          rlist(kax3i+3)=merge(1.d0,0.d0,optstat(i)%staby)
        enddo
        klist(kax+1)=ktflist+kax1
        klist(kax+2)=ktflist+kax2
        klist(kax+3)=ktflist+kax3
        klist(kax+4)=ktflist+kax4
        kx=ktflist+kax
      endif
      ip=0
      do 10 i=1,ncc
        kpa=icalc1(1,i)
        kpb=icalc1(2,i)
        k=icalc1(3,i)
        jm=0
        fm=1.d100
        if(k .gt. 0)then
          do 20 j=1,flv%nfc
            if(flv%mfitp(j) .gt. 0)then
              if(kpa .eq. flv%ifitp(j) .and. kpb .eq. flv%ifitp1(j))then
                if(flv%kfit(j) .eq. k)then
                  jm=j
                  go to 21
                endif
              endif
            elseif(flv%mfitp(j) .lt. 0)then
              if(tftype1fit(k))then
                if(kpa .eq. flv%ifitp(j) .and.
     $               kpb .eq. flv%ifitp1(j))then
                  if(flv%kfit(j) .eq. k)then
                    jm=j
                    go to 21
                  endif
                endif
              else
                if(kpa .ge. flv%ifitp(j) .and.
     $               kpb .le. flv%ifitp1(j))then
                  if(flv%kfit(j) .eq. k)then
                    if(abs(flv%fitval(j)) .lt. fm)then
                      jm=j
                      fm=flv%fitval(j)
                    endif
                  endif
                endif
              endif
            endif
 20       continue
 21       if(jm .eq. 0)then
            vout=' #######'
            vout(lfs+1:lfs+3)='  #'
            nl=len_trim(nlist(k))
            fun=nlist(k)(1:nl)
          else
            lfs1=lfs+1
            write(vout(lfs1:lfs1+2),9001)abs(flv%mfitp(jm))-1
 9001       format(I3)
c     write(*,*)jm,flv%mfitp(jm),lfs,vout(lfs+1:lfs+2)
            if(flv%mfitp(jm) .gt. 0)then
              if(flv%ifitp(jm) .ne. flv%ifitp1(jm))then
                if(k .eq. mfitbx .or. k .eq. mfitby
     $               .or. k .eq. mfitbz)then
                  vout(1:lfs)=autofg(flv%fitval(jm),forms)
                else
                  vout(1:lfs)=autofg(flv%fitval(jm)/scale(k),forms)
c                  x=tgfun(k,flv%ifitp(jm),0)
c                  vout(1:lfs)=autofg(x/scale(k),forms)
                endif
              else
                vout(1:lfs)=autofg(flv%fitval(jm)/scale(k),forms)
              endif
              nl=len_trim(nlist(k))
              fun=nlist(k)(1:nl)
            else
              vout(1:lfs)=autofg(flv%fitval(jm)/scale(k),forms)
              nl=len_trim(nlist(k))
              fun=nlist(k)(1:nl)//'M'
            endif
          endif
        else
          vout=' 0.0     '
          vout(lfs+1:lfs+3)='  1'
          fun(4:)=autofg(dble(-k),'S10.7')
          call trim(fun(4:))
          fun(1:3)='FUN'
        endif
        call elname(kpa,namea)
        if(kpa .ne. kpb)then
          call elname(kpb,name1)
          ln=len_trim(namea)
c          name=namea
c          name(ln+1:namel)='/'//name1
          name=namea(1:ln)//'/'//name1
        else
          name1=' '
          name=namea
        endif
        if(kpa .eq. mfpnta .and. kpb .eq. mfpnta1)then
          name(namel:namel)='f'
        endif
        if(ret)then
          kaxi=ktadaloc(0,4)
          dlist(kaxi+1)=kxsalocb(0,namea,len_trim(namea))
          dlist(kaxi+2)=kxsalocb(0,name1,len_trim(name1))
          dlist(kaxi+3)=kxsalocb(0,fun,len_trim(fun))
          kaxi4=ktavaloc(0,nf)
          klist(kaxi+4)=ktflist+kaxi4
          klist(kax4+i)=ktflist+kaxi
        endif
        if(kpa .ne. kpb)then
          if(k .eq. mfitbx .or. k .eq. mfitby
     $         .or. k .eq. mfitbz)then
            do m=1,mm
              j=jshow(m)
              buf0((m-1)*lf+1:m*lf)=
     1             autofg((tgfun(k,kpb,j)/tgfun(k,kpa,j)),form)
            enddo
            buf0(mm*lf+1:)=' '
            if(ret)then
              do j=mf,nfam
                rlist(kaxi4+j-mf+1)=
     $               tgfun(k,kpb,j)/tgfun(k,kpa,j)
              enddo
            endif
          else
            do m=1,mm
              j=jshow(m)
              buf0((m-1)*lf+1:m*lf)=
     1             autofg((tgfun(k,kpb,j)-tgfun(k,kpa,j))/scale(k),form)
            enddo
            buf0(mm*lf+1:)=' '
            if(ret)then
              do j=mf,nfam
                rlist(kaxi4+j-mf+1)=
     $               tgfun(k,kpb,j)-tgfun(k,kpa,j)
              enddo
            endif
          endif
        elseif(k .le. mfittry .and. k .gt. 0)then
          do 1010 m=1,mm
            j=jshow(m)
c            if(k .le. mfitddp)then
c              write(*,*)'tfshow ',m,j,k,kpb,nfam,
c     $             utwiss(1,0,1),
c     $             utwiss(k,j,itwissp(kpb))
c            endif
            buf0((m-1)*lf+1:m*lf)=
     1           autofg(tgfun(k,kpb,j)/scale(k),form)
1010      continue
          buf0(mm*lf+1:)=' '
          if(ret)then
            do j=mf,nfam
              rlist(kaxi4+j-mf+1)=tgfun(k,kpb,j)
            enddo
          endif
        elseif(k .gt. 0)then
          buf0=autofg(tgfun(k,kpb,0)/scale(k),'10.6')
          if(ret)then
            do j=mf,nfam
              rlist(kaxi4+j-mf+1)=tgfun(k,kpb,j)
            enddo
          endif
        else
          buf0=autofg(df(nqcol1-k),'10.6')
          if(ret)then
            do j=mf,nfam
              rlist(kaxi4+j-mf+1)=df(nqcol1-k)
            enddo
          endif
        endif
        if(lfs .gt. 6)then
          buf1=name(1:namel)//' '//fun//vout(1:lfs+3)//' '//buf0
        else
          buf1=name(1:namel)//     fun//vout(1:lfs+3)//     buf0
        endif
        lb=len_trim(buf1)
        if(ip .ne. 0 .and. ip+lb+1 .gt. lw)then
          write(lfno,'(a)')buf2(1:ip)
          ip=0
        endif
        if(ip .eq. 0)then
          buf2(1:lb)=buf1(1:lb)
          ip=lb
        else
          buf2(ip+1:ip+lb+1)=' '//buf1(1:lb)
          ip=ip+lb+1
        endif
10    continue
      if(ip .gt. 0)then
        write(lfno,'(A)')buf2(1:ip)
      endif
      return
      end

      logical*4 function tftype1fit(k)
      use tffitcode
      implicit none
      integer*4 ,intent(in):: k
      tftype1fit=k .eq. mfitnx .or. k .eq. mfitny .or.
     $     (k .ge. mfitleng .and. k .le. mfitgz)
      return
      end

      character*3 function famlabel(k)
      implicit none
      integer*4 ,intent(in):: k
      if(k .eq. 0)then
        famlabel='   '
      elseif(k .gt. 9)then
        write(famlabel,'(a,i2)')'x',k
      elseif(k .gt. 0)then
        write(famlabel,'(a,i1,a)')'x',k,' '
      elseif(k .gt. -10)then
        write(famlabel,'(a,i1,a)')'y',-k,' '
      else
        write(famlabel,'(a,i2)')'y',-k
      endif
      return
      end
