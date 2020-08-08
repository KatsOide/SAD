      subroutine tracka(latt,kptbl,x,px,y,py,z,g,dv,pz,
     $                  sjx,sjxjx,sjy,sjyjy,sjz,sjzjz,zbuf,lsp,zbufa)
      use tfstk
      use ffs_flag
      use tmacro
      use ffs_pointer, only:idelc,idtypec
      implicit none
      type (sad_symdef), pointer :: symd
      integer*8 latt(nlat)
      integer*4 kptbl(np0,6)
      real*8 x(np0),px(np0),y(np0),py(np0),z(np0),g(np0),dv(np0),pz(np0)
      real*8 sjx(np0),sjxjx(np0),sjy(np0),sjyjy(np0),sjz(np0),sjzjz(np0)
      complex*16 zbuf(lsp,3,np0)
      integer*4 lsp
      logical zbufa

      logical cmplot1,fourie1
      integer*8 kal,kaf,
     $     kv,kaj,ke,kax,ix0pl,iy0pl,ixjpl,ixnpl,ixppl,
     $     izwrk,izexp,iajx,iaord,ib,iaexp,iiord,izwork,lp0
      integer*4 np,jpcm,level,nsmear,
     $     n,j,i,ip,jp,irtc,nsm,lpa,k,
     $     kp,lpl,kpl,lfnplt,maxord,nres,ndim,lcoeff,lv
      real*8 emx,rgetgl1,emz,es,ajx,ajy,ajyjy,ajz,ajxjx,ajzjz
      real*8 xa(8),sa(6),ss(6,6)

      character*32 vname
      data vname/'RESULTSOFTRACKING'/
      character*2 ord
      integer*4 itfuplevel,itfdownlevel
      lp0=latt(1)+kytbl(kwmax,idtypec(1))+1
      np=np0
      emx=sqrt(abs(rgetgl1('EMITX')+rgetgl1('EMITY')))
      if(rfsw)then
        emz=sqrt(abs(rgetgl1('EMITZ')))
        if(emz .le. 0.d0)then
          emz=abs(rgetgl1('SIGE'))
        endif
      else
        emz=abs(rgetgl1('SIGE'))
      endif
      call tsetdvfs
      if(trpt)then
        call ttinit(latt,x,px,y,py,z,g,dv,pz)
        call ttstat(np0,x,px,y,py,z,g,dv,pz,
     1              'Entrance',sa,ss,es,
     $              .false.,.true.,outfl)
      else
        call tinip(np,x,px,y,py,z,g,dv,emx,emz,codin,dvfs,.false.)
      endif
      if(.not. trpt)then
        call tclr(sjx,np0)
        call tclr(sjxjx,np0)
        call tclr(sjy,np0)
        call tclr(sjyjy,np0)
        call tclr(sjz,np0)
        call tclr(sjzjz,np0)
      endif
      jpcm=abs(kptbl(1,3))
      level= itfuplevel()
      kaf=ktfsymbolz(vname,len_trim(vname),symd)
      kax=0
      if(symd%downval .eq. 0)then
        kv=0
      else
        kv=1
        kal=ktadaloc(0,2)
        klist(kal)=ktfsymbol+ktfcopy1(kaf)
      endif
      cmplot1=jpcm .ne. 0
      fourie1=fourie .and. lsp .gt. 1
      nsmear=max(100,nturn/50)/100*100
      do n=1,nturn+1
        if(.not. trpt)then
          if(kv .gt. 0)then
            level= itfuplevel()
            if(kax .gt. 0)then
              call tflocal1(kax)
            endif
            kax=ktadaloc(0,7)
            klist(kal+2)=ktflist+kax
            do j=1,7
              kaj=ktavaloc(0,np0)
              klist(kax+j)=ktflist+kaj
            enddo
          endif
          do 7020 i=1,np0
            ip=kptbl(i,1)
            if(ip .le. np)then
              xa(1)=x(ip)
              xa(2)=px(ip)
              xa(3)=y(ip)
              xa(4)=py(ip)
              xa(5)=z(ip)
              xa(6)=g(ip)
              call tconv(xa,xa,1)
              if(kv .gt. 0)then
                xa(7)=1.d0
                do j=1,7
                  rlist(ktfaddr(klist(kax+j))+i)=xa(j)
                enddo
              endif
              call tsub(xa,codin,xa,6)
              jp=abs(kptbl(i,3))
              if(jp .ne. 0)then
                if(kptbl(i,3) .lt. 0)then
                  call liemap(-2,xa,rlist(jp+6),1023)
                  call tmap(xa,xa,1)
                else
                  call tmap(xa,xa,1)
                endif
                call tmov(xa,rlist(jp),6)
              else
                call tmap(xa,xa,1)
              endif
              if(zbufa)then
                zbuf(n,1,i)=dcmplx(xa(1),xa(2))
                zbuf(n,2,i)=dcmplx(xa(3),xa(4))
                zbuf(n,3,i)=dcmplx(xa(5),xa(6))
              endif
              ajx=xa(1)**2+xa(2)**2
              ajxjx=ajx**2
              ajy=xa(3)**2+xa(4)**2
              ajyjy=ajy**2
              ajz=xa(5)**2+xa(6)**2
              ajzjz=ajz**2
              sjxjx(i)=sjxjx(i)+ajxjx
              sjx(i)=sjx(i)+ajx
              sjyjy(i)=sjyjy(i)+ajyjy
              sjy(i)=sjy(i)+ajy
              sjzjz(i)=sjzjz(i)+ajzjz
              sjz(i)=sjz(i)+ajz
            elseif(kv .gt. 0)then
              call tclr(xa,7)
              do j=1,7
                rlist(ktfaddr(klist(kax+j))+i)=xa(j)
              enddo
            endif
7020      continue
          if(kv .gt. 0)then
            rlist(kal+1)=dble(n-1)
            call tfleval(klist(kal-3),ke,.true.,irtc)
            level= itfdownlevel()
            if(irtc .ne. 0)then
              write(outfl,*)
     $             'Function execution error at ',
     $             n,ord(n),' turn: code =',irtc
            endif
          endif
          if(smearp .and.
     $         (mod(n,nsmear) .eq. 0 .or. n .eq. nturn+1))then
            nsm=mod(n-1,nsmear)+1
            if(nsm .gt. 1)then
              call tsmear(n,nsmear,np,kptbl,
     1             sjx,sjxjx,sjy,sjyjy,sjz,sjzjz,outfl)
            endif
            call tclr(sjx,np0)
            call tclr(sjxjx,np0)
            call tclr(sjy,np0)
            call tclr(sjyjy,np0)
            call tclr(sjz,np0)
            call tclr(sjzjz,np0)
          endif
        endif
        do 7030 i=1,nspect
          lpa=ilist(2,lspect+i)
          jp=ilist(1,lpa+1)
          if(kptbl(jp,1).le. np)then
            k=ilist(1,lpa)
            if(k .gt. 0)then
              kp=abs(kptbl(jp,3))+mod(k-1,3)*2
            else
              kp=abs(kptbl(jp,3))+mod(-k-1,3)*2+6
            endif
c         write(*,*)jp,k,kptbl(jp,3),rlist(kp),rlist(kp+1)
            call tspect(i,n,rlist(kp))
          endif
7030    continue
        do 7110 i=1,nplot
          lpl=ilist(2,lplot+i)
          jp=ilist(1,lpl+1)
          if(kptbl(jp,1) .le. np)then
            if(n .ge. ilist(1,lpl+2))then
              kpl=ilist(1,lplot+i)+(n-ilist(1,lpl+2)+1)*2-1
              k=ilist(1,lpl)
              if(k .gt. 0)then
                kp=abs(kptbl(jp,3))+mod(k-1,3)*2
              else
                kp=abs(kptbl(jp,3))+mod(-k-1,3)*2+6
              endif
              rlist(kpl  )=rlist(kp  )
              rlist(kpl+1)=rlist(kp+1)
            endif
          endif
7110    continue
        if(n .gt. nturn)then
          go to 1010
        endif
        call tturn(np,x,px,y,py,z,g,dv,pz,kptbl,n)
        if(np .le. 0)then
          go to 1010
        endif
        continue
      enddo
1010  if(trpt .and. outfl .gt. 0)then
        call ttstat(np,x,px,y,py,z,g,dv,pz,
     1              'Exit',sa,ss,es,
     $              .false.,.true.,outfl)
      else
        lfnplt=nint(rgetgl1('PHSPLOTS'))
        if(lfnplt .gt. 0)then
          ix0pl=ktaloc(np0*3)
          iy0pl=ktaloc(np0*3)
          ixjpl=ktaloc(np0*3)
          ixnpl=ktaloc(np0*3)
          ixppl=ktaloc(np0*3)
          izwrk=ktaloc(lsp*2)
          if(zbufa)then
            call tphplt(np,zbuf,lsp,
     1              rlist(ix0pl),rlist(iy0pl),rlist(ixjpl),rlist(ixnpl),
     1              rlist(ixppl),
     1              rlist(izwrk),
     1              nturn+1,np0,kptbl,lfnplt)
          else
            write(outfl,'(a)')'Not enough memory for PHSPLOTS.'
          endif
          call tfree(izwrk)
          call tfree(ixppl)
          call tfree(ixnpl)
          call tfree(ixjpl)
          call tfree(iy0pl)
          call tfree(ix0pl)
        endif
        if(fourie1)then
          maxord=nint(rgetgl1('MAXORDER'))
          nres=nint(rgetgl1('ADDTERMS'))
          if(rfsw)then
            ndim=3
            lcoeff=((maxord+1)*(3+maxord*(1+maxord*2)))/3+2+nres
            lv=(lcoeff-3)*2+3
          else
            ndim=2
            lcoeff=maxord*(maxord+1)+2+nres
            lv=(lcoeff-2)*2+2
          endif
          izexp=ktaloc(lsp*6)
          iajx =ktaloc(lsp*3)
          izwork=ktaloc(lsp*2)
          iaord=ktaloc(lcoeff*3)
          ib=ktaloc(lv)
          iaexp=ktaloc(lsp*2)
          iiord=ktaloc((lcoeff*3)/2+1)
        else
          ib=1
          iiord=1
          maxord=0
          iaord=1
          iajx=1
          izexp=1
          iaexp=1
          izwork=1
        endif
        if(zbufa)then
          call wjfit(np,zbuf,lsp,np0,kptbl,maxord,lcoeff,
     1             nres,
     1             rlist(ib),
     1             rlist(iiord),rlist(iaord),rlist(iajx),
     1             rlist(izexp),rlist(iaexp),rlist(izwork),
     1             rfsw,outfl)
        elseif(smearp)then
          write(outfl,'(a)')'Not enough memory for FOURIER.'
        endif
        if(fourie1)then
          call tfree(iiord)
          call tfree(iaexp)
          call tfree(ib)
          call tfree(iaord)
          call tfree(izwork)
          call tfree(iajx)
          call tfree(izexp)
        endif
      endif
      call tltrm(kptbl)
      if(kv .gt. 0)then
        call tflocal1(kal)
      endif
      level= itfdownlevel()
      return
      end
