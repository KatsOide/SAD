      subroutine tsol(np,x,px,y,py,z,g,dv,pz,
     $     latt,k,kstop,ke,sol,kptbl,la,n,
     $     nwak,nextwake,out)
      use kyparam
      use tfstk
      use ffs
      use ffs_wake
      use ffs_pointer,
     $     only:idelc,direlc,elatt,idtypec,idvalc,pnamec,lpnamec
      use sad_main
      implicit none
      real*8 conv
      parameter (conv=3.d-16)
      type (sad_comp), pointer::cmp
      integer*4 la1,la
      parameter (la1=15)
      integer*4 k,kbz,np
      real*8 x(np0),px(np0),y(np0),py(np0),z(np0),g(np0),dv(np0),pz(np0)
      real*8 tfbzs,fw,bzs,rho,al,theta,phi,phix,phiy,rhoe,
     $     psi1,psi2,bz1,rho1,dx,dy,rot,fb1,fb2,chi1m,chi2m,rtaper,
     $     harm,w,ph
      integer*8 latt(nlat),l1,lp
      integer*4 kptbl(np0,6),nwak,nextwake,n,
     $     i,ke,l,lt,mfr,itab(np),izs(np),
     $     kdx,kdy,krot,kstop,kb,lwl,lwt
      integer*8 iwpl,iwpt,itp
      logical*4 sol,enarad,dir,out,fringe,autophi
      l1=latt(k)
      if(sol)then
        kb=k
      else
        kb=k+1
      endif
      do 10 i=kb,nlat
        if(idtypec(i) .eq. icSOL)then
          if(rlist(idvalc(i)+ky_BND_SOL)
     $         .ne. 0.d0)then
            ke=i
            go to 20
          endif
        endif
 10   continue
      write(*,*)' ???-TRACK-?Missing end of solenoid ',
     $     pname(idelc(k))(1:lpnamec(k))
      ke=nlat
 20   bzs=tfbzs(k,kbz)
      if(bzs .eq. 0.d0)then
        rho=1.d50
      else
        rho=1.d0/bzs
      endif
      if(.not. sol)then
        call trots(np,x,px,y,py,z,g,dv,
     $       rlist(l1+ky_CHI1_SOL),
     $       rlist(l1+ky_CHI2_SOL),
     $       rlist(l1+ky_CHI3_SOL),
     $       rlist(l1+ky_DX_SOL),
     $       rlist(l1+ky_DY_SOL),
     $       rlist(l1+ky_DZ_SOL),
     $       .true.)
        fringe=rlist(l1+ky_FRIN_SOL) .eq. 0.d0      
        if(fringe)then
          call tsfrin(np,x,px,y,py,z,g,bzs)
        endif
        if(rad .and. rlist(l1+ky_RAD_SOL) .eq. 0.d0)then
          call tserad(np,x,px,y,py,g,dv,l1,rho)
        endif
        sol=.true.
      endif
      iwpt=0
      iwpl=0
      do l=kb,min(ke,kstop)
        if(la .le. 0)then
          call tapert(l,latt,x,px,y,py,z,g,dv,pz,kptbl,np,n,
     $         0.d0,0.d0,0.d0,0.d0,
     $         -alost,-alost,alost,alost,0.d0,0.d0,0.d0,0.d0)
          if(np .le. 0)then
            return
          endif
          la=la1
        else
          la=la-1
        endif
        lt=idtypec(l)
        lp=elatt%comp(l)
        call loc_comp(lp,cmp)
        if(l .eq. nextwake)then
          iwpl=abs(kwaketbl(1,nwak))
          if(iwpl .ne. 0)then
            lwl=(ilist(1,iwpl-1)-2)/2
          else
            lwl=0
          endif
          iwpt=abs(kwaketbl(2,nwak))
          if(iwpt .ne. 0)then
            lwt=(ilist(1,iwpt-1)-2)/2
          else
            lwt=0
          endif
          fw=(abs(charge)*e*pbunch*anbunch/amass)/np0*.5d0
          kdx=kytbl(kwDX,lt)
          if(kdx .ne. 0)then
            dx=cmp%value(kdx)
          else
            dx=0.d0
          endif
          kdy=kytbl(kwDY,lt)
          if(kdy .ne. 0)then
            dy=cmp%value(kdy)
          else
            dy=0.d0
          endif
          krot=kytbl(kwROT,lt)
          if(krot .ne. 0)then
            rot=cmp%value(krot)
          else
            rot=0.d0
          endif
          call txwake(np,x,px,y,py,z,g,dv,
     $         dx,dy,rot,
     $         int(anbunch),
     $         fw,lwl,rlist(iwpl),lwt,rlist(iwpt),
     $         p0,h0,itab,izs,.true.)
        endif
        if(lt .eq. icDRFT)then
          al=cmp%value(ky_L_DRFT)
          if(spac)then
            call spdrift_solenoid(np,x,px,y,py,z,g,dv,pz,al,bzs,
     $           cmp%value(ky_RADI_DRFT),n,l,latt,kptbl)
          elseif(rad .and. cmp%value(ky_RAD_DRFT) .eq. 0.d0)then
            call tsdrad(np,x,px,y,py,z,g,dv,al,rho)
          else
            call tdrift_solenoid(np,x,px,y,py,z,g,dv,pz,al,bzs)
          endif
        elseif(lt .eq. icBEND)then
          al=cmp%value(ky_L_BEND)
          theta=cmp%value(ky_ROT_BEND)
     $         +cmp%value(ky_DROT_BEND)
          phi=cmp%value(ky_ANGL_BEND)+cmp%value(ky_K0_BEND)
          phiy= phi*cos(theta)
          phix= phi*sin(theta)
          enarad=rad .and. al .ne. 0.d0
     $         .and. cmp%value(ky_RAD_BEND) .eq. 0.d0
          call tdrift(np,x,px,y,py,z,g,dv,pz,al,bzs,phiy,phix)
        elseif(lt .eq. icQUAD)then
          al=cmp%value(ky_L_QUAD)
          dir=direlc(l) .gt. 0.d0
          if(dir)then
            mfr=nint(cmp%value(ky_FRMD_QUAD))
          else
            mfr=nint(cmp%value(ky_FRMD_QUAD))
            mfr=mfr*(11+mfr*(2*mfr-9))/2
          endif
          rtaper=1.d0
          if(rad .and. radcod .and. radtaper)then
            rtaper=1.d0+(gettwiss(mfitddp,l)+gettwiss(mfitddp,l+1))*.5d0
          endif
          itp=cmp%param
          call tquads(np,x,px,y,py,z,g,dv,pz,l,al,
     $         cmp%value(ky_K1_QUAD)*rtaper,bzs,
     $         cmp%value(ky_DX_QUAD),cmp%value(ky_DY_QUAD),
     1         rlist(itp+4),rlist(itp+2),rlist(itp+3),
     1         cmp%value(ky_RAD_QUAD),cmp%value(ky_FRIN_QUAD) .eq. 0.d0,
     $         rlist(itp+6)*rtaper,rlist(itp+7)*rtaper,
     $         rlist(itp+8)*rtaper,rlist(itp+9)*rtaper,
     $         mfr,cmp%value(ky_EPS_QUAD),l,dir)
        elseif(lt .eq. icMULT)then
          al=cmp%value(ky_L_MULT)
          dir=direlc(l) .gt. 0.d0
          if(dir)then
            mfr=nint(cmp%value(ky_FRMD_MULT))
            psi1=cmp%value(ky_ANGL_MULT)*
     $           cmp%value(ky_E1_MULT)
            psi2=cmp%value(ky_ANGL_MULT)*
     $           cmp%value(ky_E2_MULT)
            fb1=cmp%value(ky_FB1_MULT)
            fb2=cmp%value(ky_FB2_MULT)
            chi1m=cmp%value(ky_CHI1_MULT)
            chi2m=cmp%value(ky_CHI2_MULT)
          else
            mfr=nint(cmp%value(ky_FRMD_MULT))
            mfr=mfr*(11+mfr*(2*mfr-9))/2
            psi1=cmp%value(ky_ANGL_MULT)*
     $           cmp%value(ky_E2_MULT)
            psi2=cmp%value(ky_ANGL_MULT)*
     $           cmp%value(ky_E1_MULT)
            fb2=cmp%value(ky_FB1_MULT)
            fb1=cmp%value(ky_FB2_MULT)
            chi1m=-cmp%value(ky_CHI1_MULT)
            chi2m=-cmp%value(ky_CHI2_MULT)
          endif
          harm=cmp%value(ky_HARM_MULT)
          if(harm .eq. 0.d0)then
            w=pi2*cmp%value(ky_FREQ_MULT)/c
          else
            w=omega0*harm/c
          endif
          autophi=cmp%value(ky_APHI_MULT) .ne. 0.d0
          ph=cmp%value(ky_DPHI_MULT)
          if(autophi)then
            ph=ph+gettwiss(mfitdz,l)*w
          endif
          rtaper=1.d0
          if(rad .and. radcod .and. radtaper)then
            rtaper=1.d0+(gettwiss(mfitddp,l)+gettwiss(mfitddp,l+1))*.5d0
          endif
          itp=cmp%param
          call tmulti(np,x,px,y,py,z,g,dv,pz,al,
     $         cmp%value(ky_K0_MULT),bzs,
     $         cmp%value(ky_ANGL_MULT),psi1,psi2,
     1         cmp%value(ky_DX_MULT),cmp%value(ky_DY_MULT),
     $         cmp%value(ky_DZ_MULT),
     $         chi1m,chi2m,cmp%value(ky_ROT_MULT),
     $         cmp%value(ky_DROT_MULT),
     $         cmp%value(ky_EPS_MULT),cmp%value(ky_RAD_MULT) .eq. 0.d0,
     $         cmp%value(ky_FRIN_MULT) .eq. 0.d0,
     $         rlist(itp+1)*rtaper,rlist(itp+2)*rtaper,
     $         rlist(itp+3)*rtaper,rlist(itp+4)*rtaper,
     $         mfr,fb1,fb2,
     $         cmp%value(ky_VOLT_MULT),w,cmp%value(ky_PHI_MULT),ph,
     $         cmp%value(ky_RADI_MULT),rtaper,autophi,
     $         n,l,latt,kptbl)
        elseif(lt .eq. icSOL)then
          if(l .eq. ke)then
            fringe=cmp%value(ky_FRIN_SOL) .eq. 0.d0      
            if(fringe)then
              call tsfrin(np,x,px,y,py,z,g,-bzs)
            endif
            if(rad .and. cmp%value(ky_RAD_SOL) .eq. 0.d0)then
              call tserad(np,x,px,y,py,g,dv,lp,rho)
            endif
            call trots(np,x,px,y,py,z,g,dv,
     $           cmp%value(ky_CHI1_SOL),
     $           cmp%value(ky_CHI2_SOL),
     $           cmp%value(ky_CHI3_SOL),
     $           cmp%value(ky_DX_SOL),
     $           cmp%value(ky_DY_SOL),
     $           cmp%value(ky_DZ_SOL),
     $           .false.)
          else
            bz1=tfbzs(l,kbz)
            fringe=cmp%value(ky_FRIN_SOL) .eq. 0.d0      
            if(fringe)then
              call tsfrin(np,x,px,y,py,z,g,bz1-bzs)
            endif
            if(bz1 .eq. bzs)then
              rhoe=1.d50
            else
              rhoe=1/(bz1-bzs)
            endif
            if(bz1 .eq. 0.d0)then
              rho1=1.d50
            else
              rho1=1.d0/bz1
            endif
            rho=rho1
            bzs=bz1
            if(rad .and. cmp%value(ky_RAD_SOL) .eq. 0.d0)then
              call tserad(np,x,px,y,py,g,dv,lp,rhoe)
            endif
          endif
        elseif(lt .eq. icMAP)then
          call temap(np,np0,x,px,y,py,z,g,dv,l,n,kptbl)
        elseif(lt .eq. icAprt)then
          call tapert1(l,latt,x,px,y,py,z,g,dv,pz,
     1         kptbl,np,n)
          if(np .le. 0)then
            return
          endif
        elseif(out .and. lt .ne. icMARK .and. lt .ne. icMONI)then
          out=.false.
          write(*,*)
     $'Only DRIFT, BEND, QUAD, SOL, MULT, MAP, APERT, MARK, MON '//
     $         'are supported in SOL: ',
     $         lt
        endif
        if(l .eq. nextwake .and. l .ne. ke)then
          call txwake(np,x,px,y,py,z,g,dv,
     $         dx,dy,rot,
     $         int(anbunch),
     $         fw,lwl,rlist(iwpl),lwt,rlist(iwpt),
     $         p0,h0,itab,izs,.false.)
          nwak=nwak+1
          if(nwak .gt. nwakep)then
            nextwake=0
          else
            nextwake=iwakeelm(nwak)
          endif
        endif
      enddo
      return
      end
c     
      subroutine trads(x,px,y,py,g,dv,brad,al)
      use tfstk
      use ffs_flag
      use tmacro
      implicit none 
      real*8 x,px,y,py,g,dv,brad,al,
     $     alc,er,pr,p,hh,dp,de,h1,tran
      alc=al*crad
      if(rfluct)then
        er=c/amass*erad
        dv=(tran()-.5d0)*3.46410161513775461d0
        pr=1.d0+g
        p=p0*pr
        hh=1.d0+p**2
        dp=-hh*brad*alc
        de=er*sqrt(hh*brad)/p*dp*hh
        g=g+dp+sqrt(abs(de))*dv
      else
        pr=1.d0+g
        hh=1.d0+(p0*pr)**2
        dp=-hh*brad*alc
        g=g+dp
      endif
      pr=1.d0+g
      h1=p2h(p0*pr)
c      h1=p0*pr*sqrt(1.d0+1.d0/(p0*pr)**2)
      dv=-g*(1.d0+pr)/h1/(h1+pr*h0)
      return
      end

      subroutine trots(np,x,px,y,py,z,g,dv,
     $     chi1,chi2,chi3,dx,dy,dz,ent)
      use tfstk
      implicit none
      integer*4 np,i
      real*8 x(np),px(np),y(np),py(np),z(np),g(np),dv(np),
     $     chi1,chi2,chi3,
     $     cchi1,schi1,cchi2,schi2,cchi3,schi3,
     $     r11,r12,r13,r21,r22,r23,r31,r32,r33,
     $     pr,pxi,pyi,a,pzi,xi,yi,xf,yf,zf,pxf,pyf,pzf,
     $     dx,dy,dz
      logical*4 ent
      cchi1=cos(chi1)
      schi1=sin(chi1)
      cchi2=cos(chi2)
      schi2=sin(chi2)
      cchi3=cos(chi3)
      schi3=sin(chi3)
      r11= cchi1*cchi3+schi1*schi2*schi3
      r12=-cchi2*schi3
      r13= schi1*cchi3-cchi1*schi2*schi3
      r21=-schi1*schi2*cchi3+cchi1*schi3
      r22= cchi2*cchi3
      r23= cchi1*schi2*cchi3+schi1*schi3
      r31=-schi1*cchi2
      r32=-schi2
      r33= cchi1*cchi2
      if(ent)then
        do i=1,np
          pr=1.d0+g(i)
          pxi=px(i)
          pyi=py(i)
          a=min(.9999d0,pxi**2+pyi**2)
          pzi=1.d0+sqrt1(-a)
c          pzi=1.d0-a/(sqrt(1.d0-a)+1.d0)
          xi=x(i)
          yi=y(i)
          xf =r11*xi +r12*yi
          yf =r21*xi +r22*yi
          zf =r31*xi +r32*yi
          pxf=r11*pxi+r12*pyi+r13*pzi
          pyf=r21*pxi+r22*pyi+r23*pzi
          pzf=r31*pxi+r32*pyi+r33*pzi
          px(i)=pxf
          py(i)=pyf
          x(i)=xf-pxf/pzf*zf+dx
          y(i)=yf-pyf/pzf*zf+dy
          z(i)=z(i)+zf/pzf+dz
        enddo
      else
        do i=1,np
          pr=(1.d0+g(i))
          pxi=px(i)
          pyi=py(i)
          a=min(.9999d0,pxi**2+pyi**2)
          pzi=1.d0+sqrt1(-a)
c          pzi=1.d0-a/(sqrt(1.d0-a)+1.d0)
          xi=x(i)-dx
          yi=y(i)-dy
          xf =r11*xi +r12*yi
          yf =r21*xi +r22*yi
          zf =r31*xi +r32*yi
          pxf=r11*pxi+r12*pyi+r13*pzi
          pyf=r21*pxi+r22*pyi+r23*pzi
          pzf=r31*pxi+r32*pyi+r33*pzi
          x(i)=xf-pxf/pzf*zf
          y(i)=yf-pyf/pzf*zf
          z(i)=z(i)+zf/pzf+dz-dv(i)*dz
          px(i)=pxf
          py(i)=pyf
        enddo
      endif
      return
      end

      subroutine tsfrin(np,x,px,y,py,z,g,dbz)
      implicit none
      integer*4 np
      real*8 x(np),px(np),y(np),py(np),z(np),g(np),dbz
      integer*4 i
      real*8 x0,y0,px0,py0,z0,p,bq,bp,pr,pphi,phi0,w,c,s
c     Apply identity map if dbz == 0
      if(dbz .eq. 0.d0) then
        return
      endif
      do i=1,np
c       Map definition
c         [position @ enter]
c          x0  :=  x(i)
c         px0  := px(i)  [Px / p]
c          y0  :=  y(i)
c         py0  := py(i)  [Py / p]
c          z0  :=  z(i)
c          p   := 1 + g(i)
c
c         r    := hypot(x0, y0)
c         b    := .5 * dbz / p
c         bq   := -.25 * b
c         bp   :=  .5  * b
c         pr   := x0 * px0 + y0 * py0
c         pphi := x0 * py0 - y0 * px0
c         r1   := r * exp(bq * pphi)
c         phi  := atan(y0, x0)
c         phi1 := phi + bq * pr
c
c         [position @ exit]
c          x1  := x(i)  <- r1 * cos(phi1)
c          y1  := y(i)  <- r1 * sin(phi1)
c         px1  := px(i) <- (x1 * pr - y1 * pphi) / r1**2
c         py1  := py(i) <- (y1 * pr + x1 * pphi) / r1**2
c          z1  := z(i)  <- z(i) + bp * pr * pphi
c
c       Map expansion
c         phi = atan2(y0, x0) -> x0 = r * cos(phi), y0 = r * sin(phi)
c         phi0 := bq * pr
c         cos(phi1) = cos(phi + phi0)
c                   = cos(phi) * cos(phi0) - sin(phi) * sin(phi0)
c                   = x0/r * cos(phi0) - y0/r * sin(phi0)
c                   = (x0 * cos(phi0) - y0 * sin(phi0)) / r
c         sin(phi1) = sin(phi + phi0)
c                   = sin(phi) * cos(phi0) + cos(phi) * sin(phi0)
c                   = y0/r * cos(phi0) + x0/r * sin(phi0)
c                   = (y0 * cos(phi0) + x0 * sin(phi0)) / r
c
c          x1 = r1 * cos(phi1)
c             = (cos(phi0) *  x0 - sin(phi0) *  y0) * (r1 / r)
c          y1 = r1 * sin(phi1)
c             = (cos(phi0) *  y0 + sin(phi0) *  x0) * (r1 / r)
c         px1 = (cos(phi1) *  pr - sin(phi1) * pphi) / r1
c             = (cos(phi0) * px0 - sin(phi0) * py0) / (r1 / r)
c         py1 = (sin(phi1) *  pr + cos(phi1) * pphi) / r1
c             = (cos(phi0) * py0 + sin(phi0) * px0) / (r1 / r)
c          z1 = z0 + bp * pr * pphi
c
c         w := (r1 / r) = exp(bq * pphi)
c         c := cos(phi0)
c         s := sin(phi0)
c
c       Optimized map code
        x0  = x(i)
        y0  = y(i)
        px0 = px(i)
        py0 = py(i)
        z0  = z(i)
        p   = 1.d0 + g(i)

c       b     = .5d0 * dbz / p
c       bq    = -.25d0 * b
c       bp    =  .50d0 * b
        bq    = -.125d0 * dbz / p
        bp    = -2.d0 * bq
        pr    = x0 * px0 + y0 * py0
        pphi  = x0 * py0 - y0 * px0
        phi0  = bq * pr
        w     = exp(bq * pphi)
        c     = cos(phi0)
        s     = sin(phi0)

        x(i)  = (c *  x0 - s *  y0) * w
        y(i)  = (c *  y0 + s *  x0) * w
        px(i) = (c * px0 - s * py0) / w
        py(i) = (c * py0 + s * px0) / w
        z(i)  = z0 + bp * pr * pphi
      enddo
      return
      end

      subroutine tsfrie(trans,cod,dbz)
      implicit none
      real*8 trans(6,6),cod(6),dbz
      real*8 x0,px0,y0,py0,z0,p,bq,bp,pr,pphi,phi0,w,c,s,
     $     x1,px1,y1,py1,z1

c     Apply identity map if dbz == 0
c      if(dbz .eq. 0.d0)then
c        call tinitr(trans)
c        return
c      endif

c     Map definition
c       [position @ enter]
c        x0  := cod(1)
c       px0  := cod(2)
c        y0  := cod(3)
c       py0  := cod(4)
c        z0  := cod(5)
c        p   := 1 + cod(6)
c
c       r    := hypot(x0, y0)
c       b    := .5 * dbz / p**2
c       bq   := -.25 * b
c       bp   :=  .5  * b / p
c       pr   := x0 * px0 + y0 * py0
c       pphi := x0 * py0 - y0 * px0
c       r1   := r * exp(bq * pphi)
c       phi  := atan(y0, x0)
c       phi1 := phi + bq * pr
c
c       [position @ exit]
c        x1  := cod(1) <- r1 * cos(phi1)
c        y1  := cod(3) <- r1 * sin(phi1)
c       px1  := cod(2) <- (x1 * pr - y1 * pphi) / r1**2
c       py1  := cod(4) <- (y1 * pr + x1 * pphi) / r1**2
c        z1  := cod(5) <- z0 + bp * pr * pphi
c
c     Map expansion
c       phi = atan2(y0, x0) -> x0 = r * cos(phi), y0 = r * sin(phi)
c       phi0 := bq * pr
c       cos(phi1) = cos(phi + phi0)
c                 = cos(phi) * cos(phi0) - sin(phi) * sin(phi0)
c                 = x0/r * cos(phi0) - y0/r * sin(phi0)
c                 = (x0 * cos(phi0) - y0 * sin(phi0)) / r
c       sin(phi1) = sin(phi + phi0)
c                 = sin(phi) * cos(phi0) + cos(phi) * sin(phi0)
c                 = y0/r * cos(phi0) + x0/r * sin(phi0)
c                 = (y0 * cos(phi0) + x0 * sin(phi0)) / r
c
c        x1 = r1 * cos(phi1)
c           = (cos(phi0) *  x0 - sin(phi0) *  y0) * (r1 / r)
c        y1 = r1 * sin(phi1)
c           = (cos(phi0) *  y0 + sin(phi0) *  x0) * (r1 / r)
c       px1 = (cos(phi1) *  pr - sin(phi1) * pphi) / r1
c           = (cos(phi0) * px0 - sin(phi0) * py0) / (r1 / r)
c       py1 = (sin(phi1) *  pr + cos(phi1) * pphi) / r1
c           = (cos(phi0) * py0 + sin(phi0) * px0) / (r1 / r)
c        z1 = z0 + bp * pr * pphi
c
c       w := (r1 / r) = exp(bq * pphi)
c       c := cos(phi0)
c       s := sin(phi0)
c
c     Derivatives
c       dbq    =  bq * (-2 / p) * dp = bp * dp
c       dbp    =  bp * (-3 / p) * dp
c       dpr    =  px0 * dx + x0 * dpx + py0 * dy + y0 * dpy
c       dpphi  =  py0 * dx - y0 * dpx - px0 * dx + x0 * dpy
c       dphi0  =  bq * dpr + pr * dbq
c              =  bq * (px0 * dx + x0 * dpx + py0 * dy + y0 * dpy)
c               + bp * pr * dp
c       dw     =  w * (bq * dpphi + pphi * bp * dp)
c       d(1/w) =  (1 / w) * (-1 / w) * dw
c              = -(1 / w) * (bq * dpphi + pphi * bp * dp)
c       dc     = -s * dphi0
c       ds     =  c * dphi0
c
c        dx1 =  ( c *  dx - s *  dy) * w + (-s *  x0 - c *  y0) * w * dphi0
c             + ( x1 / w) * dw
c            =  (c *  dx - s *  dy) * w
c             -  y1 * dphi0 +  x1 * (bq * dpphi + bp * pphi * dp)
c            =  c * w *  dx + bq * (  x1 * ( py0) -   y1 * px0) *  dx
c                           + bq * (  x1 * ( -y0) -   y1 *  x0) * dpx
c             - s * w *  dy + bq * (  x1 * (-px0) -   y1 * py0) *  dy
c                           + bq * (  x1 * (  x0) -   y1 *  y0) * dpy
c                           + bp * (  x1 *  pphi  -   y1 *  pr) *  dp
c
c       dpx1 =  ( c * dpx - s * dpy) / w + (-s * px0 - c * py0) / w * dphi0
c             + (px1 * w) * d(1/w)
c            =  (c * dpx - s * dpy) / w
c             - py1 * dphi0 - px1 * (bq * dpphi + pphi * bp * dp)
c            =                bq * (-px1 * ( py0) -  py1 * px0) *  dx
c             + c / w * dpx + bq * (-px1 * ( -y0) -  py1 *  x0) * dpx
c                           + bq * (-px1 * (-px0) -  py1 * py0) *  dy
c             - s / w * dpy + bq * (-px1 * (  x0) -  py1 *  y0) * dpy
c                           + bp * (-px1 *  pphi  -  py1 *  pr) *  dp
c
c        dy1 =  ( c *  dy + s *  dx) * w + (-s *  y0 + c *  x0) * w * dphi0
c             + ( y1 / w) * dw
c            =  ( s *  dx + c *  dy) * w
c             +  x1 * dphi0 +  y1 * (bq * dpphi + bp * pphi * dp)
c            =  s * w *  dx + bq * (  y1 * ( py0) +   x1 * px0) *  dx
c                           + bq * (  y1 * ( -y0) +   x1 *  x0) * dpx
c             + c * w *  dy + bq * (  y1 * (-px0) +   x1 * py0) *  dy
c                           + bq * (  y1 * (  x0) +   x1 *  y0) * dpy
c                           + bp * (  y1 *  pphi  +   x1 *  pr) *  dp
c
c       dpy1 =  ( c * dpy + s * dpx) / w + (-s * py0 + c * px0) / w * dphi0
c             + (py1 * w) * d(1/w)
c            =  (c * dpx - s * dpy) / w
c             + px1 * dphi0 - py1 * (bq * dpphi + pphi * bp * dp)
c            =                bq * (-py1 * ( py0) +  px1 * px0) *  dx
c             + s / w * dpx + bq * (-py1 * ( -y0) +  px1 *  x0) * dpx
c                           + bq * (-py1 * (-px0) +  px1 * py0) *  dy
c             + c / w * dpy + bq * (-py1 * (  x0) +  px1 *  y0) * dpy
c                           + bp * (-py1 *  pphi  +  px1 *  pr) *  dp
c
c       dz1 = dz + pr * pphi * dbp + bp * pr * dpphi + bp + pphi * dpr
c           =                 bp * (  pr * ( py0) + pphi * px0) *  dx
c                           + bp * (  pr * ( -y0) + pphi *  x0) * dpx
c                           + bp * (  pr * (-px0) + pphi * py0) *  dy
c                           + bp * (  pr * (  x0) + pphi *  y0) * dpy
c             +          dz
c                           - 3 * bp * pr * pphi / p            *  dp
c
c     Optimized map code
      x0  = cod(1)
      px0 = cod(2)
      y0  = cod(3)
      py0 = cod(4)
      z0  = cod(5)
      p   = 1.d0 + cod(6)

c     b     = .5d0 * dbz / p**2
c     bq    = -.25d0 * b
c     bp    =  .50d0 * b / p
      bq    = -.125d0 * dbz / p**2
      bp    = -2.d0 * bq / p
      pr    = x0 * px0 + y0 * py0
      pphi  = x0 * py0 - y0 * px0
      phi0  = bq * pr
      w     = exp(bq * pphi)
      c     = cos(phi0)
      s     = sin(phi0)

      x1  = (c *  x0 - s *  y0) * w
      y1  = (c *  y0 + s *  x0) * w
      px1 = (c * px0 - s * py0) / w
      py1 = (c * py0 + s * px0) / w
      z1  = z0 + bp * pr * pphi

      cod(1) =  x1
      cod(2) = px1
      cod(3) =  y1
      cod(4) = py1
      cod(5) =  z1

c     Transfer Matrix: trans(i,j) := dcod(i)@exit / dcod(j)@enter
      trans(1,1) =  c * w + bq * (  x1 * ( py0) -   y1 * px0)
      trans(1,2) =          bq * (  x1 * ( -y0) -   y1 *  x0)
      trans(1,3) = -s * w + bq * (  x1 * (-px0) -   y1 * py0)
      trans(1,4) =          bq * (  x1 * (  x0) -   y1 *  y0)
      trans(1,5) =  0.d0
      trans(1,6) =          bp * (  x1 *   pphi -   y1 *  pr)

      trans(2,1) =          bq * (-px1 * ( py0) -  py1 * px0)
      trans(2,2) =  c / w + bq * (-px1 * ( -y0) -  py1 *  x0)
      trans(2,3) =          bq * (-px1 * (-px0) -  py1 * py0)
      trans(2,4) = -s / w + bq * (-px1 * (  x0) -  py1 *  y0)
      trans(2,5) =  0.d0
      trans(2,6) =          bp * (-px1 *   pphi -  py1 *  pr)

      trans(3,1) =  s * w + bq * (  y1 * ( py0) +   x1 * px0)
      trans(3,2) =          bq * (  y1 * ( -y0) +   x1 *  x0)
      trans(3,3) =  c * w + bq * (  y1 * (-px0) +   x1 * py0)
      trans(3,4) =          bq * (  y1 * (  x0) +   x1 *  y0)
      trans(3,5) =  0.d0
      trans(3,6) =          bp * (  y1 *   pphi +   x1 *  pr)

      trans(4,1) =          bq * (-py1 * ( py0) +  px1 * px0)
      trans(4,2) =  s / w + bq * (-py1 * ( -y0) +  px1 *  x0)
      trans(4,3) =          bq * (-py1 * (-px0) +  px1 * py0)
      trans(4,4) =  c / w + bq * (-py1 * (  x0) +  px1 *  y0)
      trans(4,5) =  0.d0
      trans(4,6) =          bp * (-py1 *   pphi +  px1 *  pr)

      trans(5,1) =          bp * (  pr * ( py0) + pphi * px0)
      trans(5,2) =          bp * (  pr * ( -y0) + pphi *  x0)
      trans(5,3) =          bp * (  pr * (-px0) + pphi * py0)
      trans(5,4) =          bp * (  pr * (  x0) + pphi *  y0)
      trans(5,5) =  1.d0
      trans(5,6) = -3.d0 * bp * pr * pphi / p

      trans(6,1) =  0.d0
      trans(6,2) =  0.d0
      trans(6,3) =  0.d0
      trans(6,4) =  0.d0
      trans(6,5) =  0.d0
      trans(6,6) =  1.d0

      return
      end
