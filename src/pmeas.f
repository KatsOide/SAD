      real*8 function pmeas(latt,twiss,gammab,size,observ,iobs,datas,
     $     sto,synchb,lfno)
      use tfstk
      use ffs
      use tffitcode
      implicit real*8 (a-h,o-z)
c----- Measure luminosity and store data (if sto). ------
      parameter (npara=59,item=20,itemn=item/2)
      integer*4 irtc
      logical stab,lumi,sto,synchb
      character*(*) observ
      integer*8 latt(nlat)
      dimension twiss(nlat,-ndim:ndim,ntwissfun),gammab(*),
     $     size(21,nlat)
      dimension iobs(*),datas(itemn,2,*)
      dimension params(npara),paramss(npara),x(42),u(42)
c
      save
c
c     begin initialize for preventing compiler warning
      pmeas=0.d0
c     end   initialize for preventing compiler warning
      nobs=iobs(1)
      if(observ.eq.'LUMI') then
        ilum=italoc(5*nobs)
        lumi=.true.
      endif
      itr=italoc(72)
      ibeam=italoc(42)
      ictrb=italoc(441)
      chargea=charge
      if(observ.eq.'LUMI') then
        charge=-charge
        call tphyzp
      endif
      emx0=0d0
 110  continue
      call tclr(codin,6)
      call tclr(beamin,21)
      call temit(rlist(itr),codin,rlist(ibeam),rlist(ictrb),
     $     .true.,i00,i00,i00,i00,
     $     .true.,params,stab,15,0)
      if(.not.stab) then
        call permes(' ','pmeas(temit) --> Unstable orbit.',' ',lfno)
      endif
      if(charge.gt.0) then
        ic=1
      else
        ic=2
      endif
      if(synchb) then
        amus0=abs(params(9))*pi2
        amus1=amus0
        amusstep=0d0
        mphi2=7
        call temits(ndim,ntwissfun,mphi2,
     $       amus0,amus1,amusstep,
     $       emxe,emye,rese,paramss,
     $       15,0,irtc)
        fey=max(1d0,emye/paramss(23))
        fex=emxe/paramss(22)
        if(fey.eq.1d0) fex=1d0
        emx0=emx0+emxe
        write(lfno,'(4(a,1pg15.8))')' Nuz=',abs(params(9)),' fex=',fex,
     $       ' fey=',fey,' demi/emi=',rese
        do i=1,nobs
          lp=iobs(i+1)
c         print *,'*** i=',i,' ***'
          call tmov(size(1,lp),x,21)
c         print *,'x:'
c         write(*,'(1p,6d11.4)') (x(j),j=1,21)
          call pbmconv(twiss,lp,x,u,.true.)
c         print *,'u:'
c         write(*,'(1p,6d11.4)') (u(j),j=1,21)
          call pbmemi(u,emxn,emyn,ax,bx,ay,by,.true.)          
c         print *,emxn,emyn,ax,bx,ay,by
          call pbmemi(u,emxn*fex,emyn*fey,ax,bx,ay,by,.false.)
          call pbmconv(twiss,lp,u,x,.false.)
c         print *,'x:'
c         write(*,'(1p,6d11.4)') (x(j),j=1,21)
          call tmov(x,size(1,lp),21)
        enddo
      else
        fex=1d0
        fey=1d0
        emx0=emx0+params(22)
      endif
      if(observ.eq.'EMIY') then
        pmeas=params(23)/params(22)*fey/fex
        if(sto) then
c         ... Accumulate data to the history buffer ....
          do i=1,nobs
            j=iobs(i+1)
            call pmeas1(twiss,j,size(1,j),datas(1,ic,i))
          enddo
        endif
      elseif(observ.eq.'SIGY') then
        s=0d0
        do i=1,nobs
          j=iobs(i+1)
          s=s + (size(6,j)
     $         /twiss(j,ndim,5))/(size(1,j)/twiss(j,ndim,2))
          if(sto) call pmeas1(twiss,j,size(1,j),datas(1,ic,i))
        enddo
        pmeas=s/dble(nobs)
      elseif(observ.eq.'LUMI') then
        if(lumi) then
c         write(*,'(1p,21d10.3)')(size(i,iobs(2)),i=1,21)
          do 112 i=1,nobs
            j=iobs(i+1)
            k=ilum+5*(i-1)
            rlist(k)=size(1,j)
            rlist(k+1)=size(4,j)
            rlist(k+2)=size(6,j)
            rlist(k+3)=twiss(j,0,mfitdx)
            rlist(k+4)=twiss(j,0,mfitdy)
            if(sto) then
c         ... Accumulate data to the history buffer ....
              call pmeas1(twiss,j,size(1,j),datas(1,ic,i))
            endif
 112      continue
          lumi=.not.lumi
          charge=-charge
          call tphyzp
          go to 110 
        endif
c       .. evaluate luminosity.
        s=0d0
        q=0d0
c       write(*,'(10i5)') (iobs(i),i=2,nobs+1)
c       write(*,'(1p,5d11.3)') ((rlist(ilum+5*(i-1)+k),k=0,4),i=1,nobs)
c       print *,charge
c       write(*,'(1p,5d11.3)') (size(1,iobs(i+1)),size(4,iobs(i+1)),
c    $       size(6,iobs(i+1)),twiss(iobs(i+1),0,15),
c    $       twiss(iobs(i+1),0,17),i=1,nobs)
        do 113 i=1,nobs
          j=iobs(i+1)
          k=ilum+5*(i-1)
          sa1=size(1,j)*size(6,j)-size(4,j)**2
          a1=0.5d0/sa1*size(6,j)
          b1=-0.5d0/sa1*size(4,j)
          c1=0.5d0/sa1*size(1,j)
          sa2=rlist(k)*rlist(k+2)-rlist(k+1)**2
          a2=0.5d0/sa2*rlist(k+2)
          b2=-0.5d0/sa2*rlist(k+1)
          c2=0.5d0/sa2*rlist(k)
          x0=twiss(j,0,mfitdx)-rlist(k+3)
          y0=twiss(j,0,mfitdy)-rlist(k+4)
          t=(a1+a2)*(c1+c2)-(b1+b2)**2
          f=x0**2*( (a1*a2*(c1+c2)-a2*b1**2-a1*b2**2)/t )
          g=2d0*x0*y0*( (b2*a1*c1+b1*a2*c2-b1*b2*(b1+b2))/t )
          h=y0**2*( (c1*c2*(a1+a2)-c2*b1**2-c1*b2**2)/t )
          s=s + exp(-(f+g+h))/pi* sqrt((a1*c1-b1**2)*(a2*c2-b2**2)/t)
c         s=s + 1d0/pi* sqrt((a1*c1-b1**2)*(a2*c2-b2**2)/t)
          q=q+1d0/sqrt(twiss(j,ndim,2)*twiss(j,ndim,5))
c         write(*,'(a,1p,10d11.4)')'a1 b1 c1 a2 b2 c2 t x0 y0:',a1,b1,
c    $         c1,a2,b2,c2,t,x0,y0
c         write(*,'(a,1p,10d11.4)')'f g h s:',f,g,h,s
          if(sto) then
c         ... Accumulate data for history. ....
            call pmeas1(twiss,j,size(1,j),datas(1,ic,i))
          endif
 113    continue
c       ... s/nobs is nearly equal to 1/(4pi*sigx*sigy) ..
        if(observ.eq.'LUMI') emx0=emx0*0.5d0
c       ... pmeas=emx/emy ...
        pmeas=(q/(2d0*pi2*s*emx0))**2
        datas(10,1,1)=pmeas
c       print *,'s emx q pmeas:',sngl(s),sngl(emx0),sngl(q),sngl(pmeas)
      endif
      charge=chargea
      call tphyzp
      call tfree(int8(ictrb))
      call tfree(int8(ibeam))
      call tfree(int8(itr))
      if(observ.eq.'LUMI') call tfree(int8(ilum))
      return
      end

      subroutine pmeas1(twiss,j,beam,datas)
c     ... Accumulate data to the history buffer ....
      use tfstk
      use ffs
      use tffitcode
      implicit real*8 (a-h,o-z)
      dimension twiss(nlat,-ndim:ndim,ntwissfun),beam(42),ue(4),xe(4)
      dimension datas(9)
c
      save
c
c
      ue(1)=twiss(j,0,7)
      ue(2)=twiss(j,0,8)
      ue(3)=twiss(j,0,9)
      ue(4)=twiss(j,0,10)
      call mc2to4(twiss,0,j,ue,xe)
      datas(1)=sqrt(beam(1))
      datas(2)=abs(xe(1))*sqrt(beam(21))
      datas(3)=sqrt(beam(6))
      datas(4)=abs(xe(3))*sqrt(beam(21))
      datas(5)=atan2(-2d0*beam(4),beam(1)-beam(6))
     $     /2d0/sqrt(beam(6)/beam(1))
      datas(6)=twiss(j,0,mfitdx)
      datas(7)=twiss(j,0,mfitdy)
      datas(8)=-beam(2)/beam(3)
      datas(9)=-beam(9)/beam(10)
c     write(*,'(i3,1x,1p,9d9.2)')charge,(datas(l),l=1,9)
      return
      end

      subroutine pbmemi(beam,emxn,emyn,ax,bx,ay,by,emival)
      use tfstk
      use ffs
      use tffitcode
      implicit real*8 (a-h,o-z)
c-----Calc normal_mode emittance from normal_mode beam_matrix.      
      logical emival
c
      dimension beam(42)
c
      ia(i,j)=((i+j+abs(i-j))**2+2*(i+j)-6*abs(i-j))/8
c
      if(emival) then
        emxn=(beam(ia(1,1))-beam(ia(1,6))**2)*
     $       (beam(ia(2,2))-beam(ia(2,6))**2)-
     $       (beam(ia(1,2))-beam(ia(1,6))*beam(ia(2,6)))**2
        emxn=sqrt(emxn)
        emyn=(beam(ia(3,3))-beam(ia(3,6))**2)*
     $       (beam(ia(4,4))-beam(ia(4,6))**2)-
     $       (beam(ia(3,4))-beam(ia(3,6))*beam(ia(4,6)))**2
        emyn=sqrt(emyn)
        bx=(beam(ia(1,1))-beam(ia(1,6))**2)/emxn
        ax=-(beam(ia(1,2))-beam(ia(1,6))*beam(ia(2,6)))/emxn
        by=(beam(ia(3,3))-beam(ia(3,6))**2)/emyn
        ay=-(beam(ia(3,4))-beam(ia(3,6))*beam(ia(4,6)))/emyn
      else
        beam(ia(1,1))=bx*emxn+beam(ia(1,6))**2
        beam(ia(1,2))=-ax*emxn+beam(ia(1,6))*beam(ia(2,6))
        beam(ia(2,2))=(1d0+ax**2)/bx*emxn+beam(ia(2,6))**2
        beam(ia(3,3))=by*emyn+beam(ia(3,6))**2
        beam(ia(3,4))=-ay*emyn+beam(ia(3,6))*beam(ia(4,6))
        beam(ia(4,4))=(1d0+ay**2)/by*emyn+beam(ia(4,6))**2
      endif 
      return
      end

      subroutine pbmconv(twiss,lp,x,u,tonorm)
      use tfstk
      use ffs
      use tffitcode
      implicit real*8 (a-h,o-z)
      logical tonorm
      dimension twiss(nlat,-ndim:ndim,ntwissfun),x(21),u(21)
c
      save
c
c-----if tonorm=.true. then real beam_matrix x is converted to 
c                      normal_mode u.
c     if tonorm=.false. then normal_mode beam_matrix x is converted to 
c                      real_space beam_matrix u.
c
c     ia(i,j)=((i+j+abs(i-j))**2+2*(i+j)-6*abs(i-j))/8
c
      r11=twiss(lp,0,11)
      r12=twiss(lp,0,12)
      r21=twiss(lp,0,13)
      r22=twiss(lp,0,14)
      det=r11*r22-r12*r21
      if(tonorm) then
        if(det.gt.1d0) then
          qr=sqrt((det-1d0)/det)
          a=sqrt(det)
          s11=-r12*qr
          s12=r22*qr
          s21=r11*qr
          s22=-r21*qr
          s33=r21*qr
          s34=r22*qr
          s43=r11*qr
          s44=r12*qr
          s13=a
          s14=0d0
          s23=0d0
          s24=a
          s31=a
          s32=0d0
          s41=0d0
          s42=a
          s16=twiss(lp,0,7)
          s26=twiss(lp,0,8)
          s36=twiss(lp,0,9)
          s46=twiss(lp,0,10)
        else
          a=sqrt(1d0-det)
          s11=a
          s12=0d0
          s21=0d0
          s22=a
          s33=a
          s34=0d0
          s43=0d0
          s44=a
          s13=-r22
          s14=r12
          s23=r21
          s24=-r11
          s31=r11
          s32=r12
          s41=r21
          s42=r22
          s16=twiss(lp,0,7)
          s26=twiss(lp,0,8)
          s36=twiss(lp,0,9)
          s46=twiss(lp,0,10)
        endif
      else
        if(det.gt.1d0) then
          qr=sqrt((det-1d0)/det)
          a=sqrt(det)
          s11=-r21*qr
          s12=-r22*qr
          s21=-r11*qr
          s22=-r12*qr
          s33=r12*qr
          s34=-r22*qr
          s43=-r11*qr
          s44=r21*qr
          s13=a
          s14=0d0
          s23=0d0
          s24=a
          s31=a
          s32=0d0
          s41=0d0
          s42=a
          s16=-s11*twiss(lp,0,7) - s12*twiss(lp,0,8) -
     $         s13*twiss(lp,0,9) - s14*twiss(lp,0,10)
          s26=-s21*twiss(lp,0,7) - s22*twiss(lp,0,8) -
     $         s23*twiss(lp,0,9) - s24*twiss(lp,0,10)
          s36=-s31*twiss(lp,0,7) - s32*twiss(lp,0,8) -
     $         s33*twiss(lp,0,9) - s34*twiss(lp,0,10)
          s46=-s41*twiss(lp,0,7) - s42*twiss(lp,0,8) -
     $         s43*twiss(lp,0,9) - s44*twiss(lp,0,10)
        else
          a=sqrt(1d0-det)
          s11=a
          s12=0d0
          s21=0d0
          s22=a
          s33=a
          s34=0d0
          s43=0d0
          s44=a
          s13=r22
          s14=-r12
          s23=-r21
          s24=r11
          s31=-r11
          s32=-r12
          s41=-r21
          s42=-r22
          s16=-s11*twiss(lp,0,7) - s12*twiss(lp,0,8) -
     $         s13*twiss(lp,0,9) - s14*twiss(lp,0,10)
          s26=-s21*twiss(lp,0,7) - s22*twiss(lp,0,8) -
     $         s23*twiss(lp,0,9) - s24*twiss(lp,0,10)
          s36=-s31*twiss(lp,0,7) - s32*twiss(lp,0,8) -
     $         s33*twiss(lp,0,9) - s34*twiss(lp,0,10)
          s46=-s41*twiss(lp,0,7) - s42*twiss(lp,0,8) -
     $         s43*twiss(lp,0,9) - s44*twiss(lp,0,10)
        endif
      endif
      if(detr.gt.1d0) then
        u(1)=s11*(s11*x(1) + s12*x(2) + s13*x(4) + s16*x(16)) + 
     -       s12*(s11*x(2) + s12*x(3) + s13*x(5) + s16*x(17)) + 
     -       s13*(s11*x(4) + s12*x(5) + s13*x(6) + s16*x(18)) + 
     -       s16*(s11*x(16) + s12*x(17) + s13*x(18) + s16*x(21))
        u(2)=s11*(s21*x(1) + s22*x(2) + s24*x(7) + s26*x(16)) + 
     -       s12*(s21*x(2) + s22*x(3) + s24*x(8) + s26*x(17)) + 
     -       s13*(s21*x(4) + s22*x(5) + s24*x(9) + s26*x(18)) + 
     -       s16*(s21*x(16) + s22*x(17) + s24*x(19) + s26*x(21))
        u(3)=s21*(s21*x(1) + s22*x(2) + s24*x(7) + s26*x(16)) + 
     -       s22*(s21*x(2) + s22*x(3) + s24*x(8) + s26*x(17)) + 
     -       s24*(s21*x(7) + s22*x(8) + s24*x(10) + s26*x(19)) + 
     -       s26*(s21*x(16) + s22*x(17) + s24*x(19) + s26*x(21))
        u(4)=s11*(s31*x(1) + s33*x(4) + s34*x(7) + s36*x(16)) + 
     -       s12*(s31*x(2) + s33*x(5) + s34*x(8) + s36*x(17)) + 
     -       s13*(s31*x(4) + s33*x(6) + s34*x(9) + s36*x(18)) + 
     -       s16*(s31*x(16) + s33*x(18) + s34*x(19) + s36*x(21))
        u(5)=s21*(s31*x(1) + s33*x(4) + s34*x(7) + s36*x(16)) + 
     -       s22*(s31*x(2) + s33*x(5) + s34*x(8) + s36*x(17)) + 
     -       s24*(s31*x(7) + s33*x(9) + s34*x(10) + s36*x(19)) + 
     -       s26*(s31*x(16) + s33*x(18) + s34*x(19) + s36*x(21))
        u(6)=s31*(s31*x(1) + s33*x(4) + s34*x(7) + s36*x(16)) + 
     -       s33*(s31*x(4) + s33*x(6) + s34*x(9) + s36*x(18)) + 
     -       s34*(s31*x(7) + s33*x(9) + s34*x(10) + s36*x(19)) + 
     -       s36*(s31*x(16) + s33*x(18) + s34*x(19) + s36*x(21))
        u(7)=s11*(s42*x(2) + s43*x(4) + s44*x(7) + s46*x(16)) + 
     -       s12*(s42*x(3) + s43*x(5) + s44*x(8) + s46*x(17)) + 
     -       s13*(s42*x(5) + s43*x(6) + s44*x(9) + s46*x(18)) + 
     -       s16*(s42*x(17) + s43*x(18) + s44*x(19) + s46*x(21))
        u(8)=s21*(s42*x(2) + s43*x(4) + s44*x(7) + s46*x(16)) + 
     -       s22*(s42*x(3) + s43*x(5) + s44*x(8) + s46*x(17)) + 
     -       s24*(s42*x(8) + s43*x(9) + s44*x(10) + s46*x(19)) + 
     -       s26*(s42*x(17) + s43*x(18) + s44*x(19) + s46*x(21))
        u(9)=s31*(s42*x(2) + s43*x(4) + s44*x(7) + s46*x(16)) + 
     -       s33*(s42*x(5) + s43*x(6) + s44*x(9) + s46*x(18)) + 
     -       s34*(s42*x(8) + s43*x(9) + s44*x(10) + s46*x(19)) + 
     -       s36*(s42*x(17) + s43*x(18) + s44*x(19) + s46*x(21))
        u(10)=s42*(s42*x(3) + s43*x(5) + s44*x(8) + s46*x(17)) + 
     -       s43*(s42*x(5) + s43*x(6) + s44*x(9) + s46*x(18)) + 
     -       s44*(s42*x(8) + s43*x(9) + s44*x(10) + s46*x(19)) + 
     -       s46*(s42*x(17) + s43*x(18) + s44*x(19) + s46*x(21))
        u(11)=s11*x(11) + s12*x(12) + s13*x(13) + s16*x(20)
        u(12)=s21*x(11) + s22*x(12) + s24*x(14) + s26*x(20)
        u(13)=s31*x(11) + s33*x(13) + s34*x(14) + s36*x(20)
        u(14)=s42*x(12) + s43*x(13) + s44*x(14) + s46*x(20)
        u(15)=x(15)
        u(16)=s11*x(16) + s12*x(17) + s13*x(18) + s16*x(21)
        u(17)=s21*x(16) + s22*x(17) + s24*x(19) + s26*x(21)
        u(18)=s31*x(16) + s33*x(18) + s34*x(19) + s36*x(21)
        u(19)=s42*x(17) + s43*x(18) + s44*x(19) + s46*x(21)
        u(20)=x(20)
        u(21)=x(21)
      else
        u(1)=s11*(s11*x(1) + s13*x(4) + s14*x(7) + s16*x(16)) + 
     -       s13*(s11*x(4) + s13*x(6) + s14*x(9) + s16*x(18)) + 
     -       s14*(s11*x(7) + s13*x(9) + s14*x(10) + s16*x(19)) + 
     -       s16*(s11*x(16) + s13*x(18) + s14*x(19) + s16*x(21))
        u(2)=s11*(s22*x(2) + s23*x(4) + s24*x(7) + s26*x(16)) + 
     -       s13*(s22*x(5) + s23*x(6) + s24*x(9) + s26*x(18)) + 
     -       s14*(s22*x(8) + s23*x(9) + s24*x(10) + s26*x(19)) + 
     -       s16*(s22*x(17) + s23*x(18) + s24*x(19) + s26*x(21))
        u(3)=s22*(s22*x(3) + s23*x(5) + s24*x(8) + s26*x(17)) + 
     -       s23*(s22*x(5) + s23*x(6) + s24*x(9) + s26*x(18)) + 
     -       s24*(s22*x(8) + s23*x(9) + s24*x(10) + s26*x(19)) + 
     -       s26*(s22*x(17) + s23*x(18) + s24*x(19) + s26*x(21))
        u(4)=s11*(s31*x(1) + s32*x(2) + s33*x(4) + s36*x(16)) + 
     -       s13*(s31*x(4) + s32*x(5) + s33*x(6) + s36*x(18)) + 
     -       s14*(s31*x(7) + s32*x(8) + s33*x(9) + s36*x(19)) + 
     -       s16*(s31*x(16) + s32*x(17) + s33*x(18) + s36*x(21))
        u(5)=s22*(s31*x(2) + s32*x(3) + s33*x(5) + s36*x(17)) + 
     -       s23*(s31*x(4) + s32*x(5) + s33*x(6) + s36*x(18)) + 
     -       s24*(s31*x(7) + s32*x(8) + s33*x(9) + s36*x(19)) + 
     -       s26*(s31*x(16) + s32*x(17) + s33*x(18) + s36*x(21))
        u(6)=s31*(s31*x(1) + s32*x(2) + s33*x(4) + s36*x(16)) + 
     -       s32*(s31*x(2) + s32*x(3) + s33*x(5) + s36*x(17)) + 
     -       s33*(s31*x(4) + s32*x(5) + s33*x(6) + s36*x(18)) + 
     -       s36*(s31*x(16) + s32*x(17) + s33*x(18) + s36*x(21))
        u(7)=s11*(s41*x(1) + s42*x(2) + s44*x(7) + s46*x(16)) + 
     -       s13*(s41*x(4) + s42*x(5) + s44*x(9) + s46*x(18)) + 
     -       s14*(s41*x(7) + s42*x(8) + s44*x(10) + s46*x(19)) + 
     -       s16*(s41*x(16) + s42*x(17) + s44*x(19) + s46*x(21))
        u(8)=s22*(s41*x(2) + s42*x(3) + s44*x(8) + s46*x(17)) + 
     -       s23*(s41*x(4) + s42*x(5) + s44*x(9) + s46*x(18)) + 
     -       s24*(s41*x(7) + s42*x(8) + s44*x(10) + s46*x(19)) + 
     -       s26*(s41*x(16) + s42*x(17) + s44*x(19) + s46*x(21))
        u(9)=s31*(s41*x(1) + s42*x(2) + s44*x(7) + s46*x(16)) + 
     -       s32*(s41*x(2) + s42*x(3) + s44*x(8) + s46*x(17)) + 
     -       s33*(s41*x(4) + s42*x(5) + s44*x(9) + s46*x(18)) + 
     -       s36*(s41*x(16) + s42*x(17) + s44*x(19) + s46*x(21))
        u(10)=s41*(s41*x(1) + s42*x(2) + s44*x(7) + s46*x(16)) + 
     -       s42*(s41*x(2) + s42*x(3) + s44*x(8) + s46*x(17)) + 
     -       s44*(s41*x(7) + s42*x(8) + s44*x(10) + s46*x(19)) + 
     -       s46*(s41*x(16) + s42*x(17) + s44*x(19) + s46*x(21))
        u(11)=s11*x(11) + s13*x(13) + s14*x(14) + s16*x(20)
        u(12)=s22*x(12) + s23*x(13) + s24*x(14) + s26*x(20)
        u(13)=s31*x(11) + s32*x(12) + s33*x(13) + s36*x(20)
        u(14)=s41*x(11) + s42*x(12) + s44*x(14) + s46*x(20)
        u(15)=x(15)
        u(16)=s11*x(16) + s13*x(18) + s14*x(19) + s16*x(21)
        u(17)=s22*x(17) + s23*x(18) + s24*x(19) + s26*x(21)
        u(18)=s31*x(16) + s32*x(17) + s33*x(18) + s36*x(21)
        u(19)=s41*x(16) + s42*x(17) + s44*x(19) + s46*x(21)
        u(20)=x(20)
        u(21)=x(21)
      endif
      return
      end
