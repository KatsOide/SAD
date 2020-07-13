      subroutine wjfit(np,zbuf,lsp,np0,kptbl,maxord,lcoeff,
     1                  nres,
     1                  b,iord,aord,ajx,zexp,aexp,zwork,rfsw,lfno)
      use tfstk
      implicit none
      integer*4 nbig,itmax,nfre
      real*8 conv
      parameter (nbig=16,conv=1.d-10,itmax=30,nfre=8)
      integer*4 np,lsp,lfno,i,j,k,l,n,np0,kptbl(np0,6),nres
      integer*4 maxord,lcoeff,lv,ll,kk,id,iter,ibigb(nbig),ndim
      real*8 aord(3,*),ajx(lsp,3),sj(3),aj,s1,s2,s3,s,smax(3),srms(3)
      real*8 b(lcoeff*2),bigb(nbig)
      real*8 u,r,u1,u2,p1,p2,p3,f,r1,r2,r3
      real*8 f1(3,nfre),fa(3,nfre),fx,fy,fz
      complex*16 zbuf(lsp,3,np),zexp(lsp,3),aexp(lsp),zs,zave(3)
      complex*16 zwork(lsp)
      integer*4 iord(3,*)
      character*6 tag
      character*14 label(2)
      character*15 vname
      logical*4 rfsw
      integer*8 iv,ial,iai,iax,iaf
      data vname/'TunesByTracking'/
      data label/'  |mx|+|my|   ','|mx|+|my|+|mz|'/
      if(maxord .gt. 0)then
        if(rfsw)then
          ndim=3
        else
          ndim=2
        endif
        call wiord(maxord,iord,aord,rfsw)
        lv=(lcoeff-ndim)*2+ndim
        write(lfno,*)
        write(lfno,9201)label(ndim-1),maxord,lv-nres*2,
     1                  maxord+1,label(ndim-1),nres*2
9201    format(' <Fourier analysis>      ',a,'<=',i3,
     1         '    fitting variables:',i5,/,
     1         '                 ',i5,' <=',a,'     ',
     1         '    fitting variables:',i5)
      else
        if(itfcontext .gt. 0)then
          iv=ktfsymbolz(vname,len(vname))-4
          call tflocal(klist(iv))
          ial=ktadalocnull(0,np0)
          klist(iv)=ktflist+ial
        else
          iv=0
        endif
        do i=1,np0
          if(kptbl(i,1) .gt. np)then
            tag=',lost:'
          else
            tag=':     '
          endif
          call wtunem(zbuf(1,1,i),zbuf(1,1,i),lsp,f1,fa,nfre)
          write(lfno,9301)i,tag,
     1         (k,(f1(j,k),fa(j,k),j=1,3),k=1,nfre)
9301      format(' P.',i3,a,
     1    ' nux       ampl        nuy       ampl        nuz       ampl',
     $         /,(i6,2x,3(0pf10.6,1p,g12.4)))
          if(iv .gt. 0)then
            iai=ktadaloc(0,4)
            klist(ial+i)=ktflist+iai
            if(kptbl(i,1) .gt. np)then
              rlist(iai+1)=0.d0
            else
              rlist(iai+1)=1.d0
            endif
            do j=1,3
              iax=ktadaloc(0,nfre)
              klist(iai+j+1)=ktflist+iax
              do k=1,nfre
                iaf=ktavaloc(0,2)
                klist(iax+k)=ktflist+iaf
                rlist(iaf+1)=f1(j,k)
                rlist(iaf+2)=fa(j,k)
              enddo
            enddo
          endif
        enddo
        return
      endif
      do i=1,np0
        if(kptbl(i,1) .gt. np)then
          tag=',lost:'
        else
          tag=':     '
        endif
        write(lfno,9203)i,tag
9203    format('   Particle #',i3,a)
        call wfres(zbuf(1,1,i),zwork,lsp,
     1             maxord+1,lcoeff,nres,
     1             iord,aord,
     1             ndim,fx,fy,fz)
        write(lfno,9204)fx,fy,fz,(iord(k,lcoeff),k=1,3)
9204    format('          Tune = ',3f15.7,/,
     1         ' Last variable = (',3i4,' )')
        call tclr(b,lv)
        do k=1,ndim
          zs=(0.d0,0.d0)
          do n=1,lsp
            zs=zs+zbuf(n,k,i)
          enddo
          zs=zs/lsp
          zave(k)=zs
          do n=1,lsp
            aj=abs(zbuf(n,k,i)-zs)
            if(aj .eq. 0.d0)then
              zexp(n,k)=(1.d0,0.d0)
            else
              zexp(n,k)=(zbuf(n,k,i)-zs)/aj
            endif
            ajx(n,k)=aj**2
          enddo
        enddo
        if(.not. rfsw)then
          call tclr(ajx(1,3),lsp)
          zs=(0.d0,0.d0)
          do n=1,lsp
            zs=zs+zbuf(n,3,i)
          enddo
          zave(3)=zs/lsp
        endif
        iter=0
101     r=0.d0
        s=0.d0
        do l=1,ndim
          u=0.d0
          do n=1,lsp
            u=u+ajx(n,l)
          enddo
          u=u/lsp
          b(l)=b(l)+u
          r=r+u**2
          s=s+b(l)**2
          do n=1,lsp
            ajx(n,l)=ajx(n,l)-u
          enddo
        enddo
        do n=1,lsp
          aexp(n)=(1.d0,0.d0)
        enddo
        do l=ndim+1,lcoeff
          ll=2*l-ndim-1
          do k=1,ndim
            id=iord(k,l)-iord(k,l-1)
            if(id .eq. 0)then
            elseif(id .eq. 1)then
              do n=1,lsp
                aexp(n)=aexp(n)*zexp(n,k)
              enddo
            elseif(id .eq. -1)then
              do n=1,lsp
                aexp(n)=aexp(n)*conjg(zexp(n,k))
              enddo
            elseif(id .gt. 0)then
              do n=1,lsp
                aexp(n)=aexp(n)*zexp(n,k)**id
              enddo
            else
              do n=1,lsp
                aexp(n)=aexp(n)*conjg(zexp(n,k))**(-id)
              enddo
            endif
          enddo
          u1=0.d0
          u2=0.d0
          p1=0.d0
          p2=0.d0
          p3=0.d0
          do n=1,lsp
            f=aord(1,l)*ajx(n,1)+aord(2,l)*ajx(n,2)+aord(3,l)*ajx(n,3)
            p1=p1+imag(aexp(n))**2
            u1=u1-f*imag(aexp(n))
            p2=p2+dble (aexp(n))**2
            u2=u2-f*dble (aexp(n))
            p3=p3+dble (aexp(n))*imag(aexp(n))
          enddo
          f=(aord(1,l)**2+aord(2,l)**2+aord(3,l)**2)*(p1*p2-p3**2)
          u =(p2*u1-p3*u2)/f
          u2=(p1*u2-p3*u1)/f
          u1=u
          b(ll  )=b(ll  )+u1
          b(ll+1)=b(ll+1)+u2
          r=r+u1**2+u2**2
          s=s+b(ll)**2+b(ll+1)**2
          do n=1,lsp
            f=u1*imag(aexp(n))+u2*dble(aexp(n))
            ajx(n,1)=ajx(n,1)+aord(1,l)*f
            ajx(n,2)=ajx(n,2)+aord(2,l)*f
            ajx(n,3)=ajx(n,3)+aord(3,l)*f
          enddo
        enddo
        if(sqrt(r/s) .lt. conv)then
        else
          iter=iter+1
          if(iter .lt. itmax)then
            go to 101
          endif
        endif
        call tclr(bigb,nbig)
        do k=1,nbig
          ibigb(k)=0
        enddo
        s1=b(1)
        s2=b(2)
        if(rfsw)then
          s3=b(3)
        else
          s3=0.d0
        endif
        r1=0.d0
        r2=0.d0
        r3=0.d0
        LOOP_L: do l=ndim+1,lcoeff
          ll=2*l-ndim-1
          s=sqrt(b(ll)**2+b(ll+1)**2)
          s1=s1+abs(aord(1,l)*s)
          s2=s2+abs(aord(2,l)*s)
          s3=s3+abs(aord(3,l)*s)
          r1=r1+(aord(1,l)*s)**2
          r2=r2+(aord(2,l)*s)**2
          r3=r3+(aord(3,l)*s)**2
          s=s*(aord(1,l)**2+aord(2,l)**2+aord(3,l)**2)
          do k=1,nbig
            if(s .ge. bigb(k))then
              do  kk=nbig,k+1,-1
                bigb(kk)=bigb(kk-1)
                ibigb(kk)=ibigb(kk-1)
              enddo
              bigb(k)=s
              ibigb(k)=l
              cycle LOOP_L
            endif
          enddo
        enddo LOOP_L
        smax(1)=s1
        smax(2)=s2
        smax(3)=s3
        srms(1)=sqrt(r1*.5d0)
        srms(2)=sqrt(r2*.5d0)
        srms(3)=sqrt(r3*.5d0)
        sj(1)=0.d0
        sj(2)=0.d0
        sj(3)=0.d0
        do n=1,lsp
          sj(1)=sj(1)+ajx(n,1)**2
          sj(2)=sj(2)+ajx(n,2)**2
          sj(3)=sj(3)+ajx(n,3)**2
        enddo
        sj(1)=sqrt(sj(1)/lsp)
        sj(2)=sqrt(sj(2)/lsp)
        sj(3)=sqrt(sj(3)/lsp)
        write(lfno,9101)(zave(k),k=1,3),
     1  '         2J        rms deviation     rms d2J      maximum 2J'
     1  ,(b(k),sj(k),srms(k),smax(k),k=1,ndim)
9101    format(
     1  '      <X>        <PX>         <Y>        <PY>   ',
     1  '      <Z>         <P>',/,1p,6g12.4/,a,/,
     1         ' x',1p,4g15.7/,
     1         ' y',1p,4g15.7/,:,
     1         ' z',1p,4g15.7)
        write(lfno,*)'  Big 16 coefficients:'
        write(lfno,*)(' mx  my  mz       coefficients      ',k=1,2)
        write(lfno,9102)((iord(j,ibigb(k)),j=1,3),
     1                 b(2*ibigb(k)-ndim-1),b(2*ibigb(k)-ndim),k=1,nbig)
9102    format(i4,2i4,1p,2g12.4,3i4,2g12.4)
c        call wjdraw(ajx,lsp,ndim,i,maxord)
      enddo
      return
      end
