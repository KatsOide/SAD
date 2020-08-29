      subroutine qtwiss(twiss,idp,la,lb,over)
      use tfstk
      use ffs
      use tffitcode
      implicit none
      integer*4 idp,la,lb
      real*8 twiss(nlat*(2*ndim+1),ntwissfun)
      real*8 trans(4,5),cod(6)
      logical*4 over
      call qtwiss1(twiss,idp,la,lb,trans,cod,.false.,over)
      return
      end

      subroutine qtwiss1(twiss,idp,la,lb,tr,cod,mat,over)
      use kyparam
      use tfstk
      use ffs
      use ffs_pointer, only:elatt,idelc,direlc,idtypec,idvalc,
     $     tsetfringepe
      use tffitcode
      use sad_main
      use ffs_seg
      use tfcsi, only:icslfno
      implicit none
      type (sad_comp), pointer :: cmp
      type (sad_dlist), pointer :: lsegp
      integer*4 idp,la,lb,ip0,l1,i,l,ip1,ip,ltyp,
     $     j,mfr,itgetfpe,k,ibb,ibg,ntfun,ipa,irtc
      integer*8 le,lp,ld
      real*8 twiss(nlat*(2*ndim+1),ntwissfun),epschop
      parameter (epschop=1.d-30)
      real*8 trans(4,5),cod(6),tr(4,5),rxy(4,5),trans6(6,6),
     $     r1,r2,r3,r4,detr,rr,sqrdet,trtr,bx0,by0,
     $     ax0,ay0,al,pxi,pyi,pxisq,pyisq,pzi,ale,alz,psi1,psi2,
     $     theta0,x,px,y,dpsix,dpsiy,bz,
     $     pr,a,dpz,trf00,dtheta,
     $     apsi1,apsi2,sspc0,sspc,vcalpha0,fb1,fb2,
     $     ak1,ftable(4),dir
      logical*4 over,coup,normal,mat,calpol0,insmat,err,seg,wspaccheck
      real*8 a11,a12,a13,a14,a21,a22,a23,a24,a31,a32,a33,
     $     a34,a41,a42,a43,a44,a15,a25,a35,a45
      real*8 u11,u12,u13,u14,u21,u22,u23,u24,u31,u32,u33,
     $     u34,u41,u42,u43,u44,u15,u25,u35,u45,cod1(6)
      equivalence (a11,trans(1,1)),(a12,trans(1,2)),(a13,trans(1,3)),
     1            (a14,trans(1,4)),(a21,trans(2,1)),(a22,trans(2,2)),
     1            (a23,trans(2,3)),(a24,trans(2,4)),
     1            (a31,trans(3,1)),(a32,trans(3,2)),(a33,trans(3,3)),
     1            (a34,trans(3,4)),(a41,trans(4,1)),(a42,trans(4,2)),
     1            (a43,trans(4,3)),(a44,trans(4,4)),
     1            (a15,trans(1,5)),(a25,trans(2,5)),(a35,trans(3,5)),
     1            (a45,trans(4,5))
      equivalence (u11,rxy(1,1)),(u12,rxy(1,2)),(u13,rxy(1,3)),
     1            (u14,rxy(1,4)),(u21,rxy(2,1)),(u22,rxy(2,2)),
     1            (u23,rxy(2,3)),(u24,rxy(2,4)),
     1            (u31,rxy(3,1)),(u32,rxy(3,2)),(u33,rxy(3,3)),
     1            (u34,rxy(3,4)),(u41,rxy(4,1)),(u42,rxy(4,2)),
     1            (u43,rxy(4,3)),(u44,rxy(4,4)),
     1            (u15,rxy(1,5)),(u25,rxy(2,5)),(u35,rxy(3,5)),
     1            (u45,rxy(4,5))
      real*8, parameter :: almostone=1.d0-1.d-16
c     begin initialize for preventing compiler warning
      normal=.true.
      sqrdet=0.d0
      dpsix=0.d0
      dpsiy=0.d0
      ax0=0.d0
      ay0=0.d0
      bx0=0.d0
      by0=0.d0
      rr=0.d0
c     end   initialize for preventing compiler warning
      if(wspaccheck())then
        over=.true.
        return
      endif
      call tclrfpe
      irad=6
      trf00=trf0
      vcalpha0=vcalpha
      if(trpt)then
        trf0=0.d0
        vcalpha=1.d0
      endif
      calpol0=calpol
      calpol=.false.
      over=.false.
      insmat=.false.
      ip0=(ndim+idp)*nlat
      ipa=ip0+la
      if(mat)then
        tr=0.d0
        tr(1,1)=1.d0
        tr(3,3)=1.d0
        tr(2,2)=1.d0
        tr(4,4)=1.d0
      else
        do i=mfitex,ntwissfun
          if(abs(twiss(ipa,i)) .lt. epschop)then
            twiss(ipa,i)=0.d0
          endif
        enddo
        cod=twiss(ipa,mfitdx:mfitddp)
        r1=twiss(ipa,mfitr1)
        r2=twiss(ipa,mfitr2)
        r3=twiss(ipa,mfitr3)
        r4=twiss(ipa,mfitr4)
        detr=r1*r4-r2*r3
        sqrdet=sqrt(1.d0-detr)
        normal=twiss(ipa,mfitdetr) .lt. 1.d0
      endif
      dvfs=0.d0
      call tesetdv(cod(6))
      sspc0=rlist(ifpos+la-1)
      ntfun=merge(ntwissfun,mfitdetr,orbitcal)
      do l=la+1,lb
c        call tfmemcheckprint1('qtwiss',l,.false.)
c        if(mod(l,1) .eq. 0)then
c          write(*,*)'qtwiss1 ',l,la,lb
c        endif
        l1=l-1
        ip1=ip0+l1
        ip=ip1+1
        ltyp=idtypec(l1)
        if(.not. orbitcal)then
          cod=twiss(ip1,mfitdx:mfitddp)
          call tesetdv(cod(6))
        endif
c        if(l .gt. 20200 .and. l .lt. 20300)then
c        if(l .gt. 20200 .and. mod(l,100) .eq. 0)then
c          write(*,'(a,2i5,1p6g15.7)')'qtwiss1 ',l,ltyp,cod
c        endif
        if(ltyp .gt. icMARK)then
          if(.not. mat)then
            twiss(ip,1:ntfun)=twiss(ip1,1:ntfun)
          endif
        else
          if(mat)then
            if(itgetfpe() .ne. 0)then
              call tclrfpe
              over=.true.
              go to 9000
            endif
            trtr=tr(1,1)+tr(2,2)+tr(3,3)+tr(4,4)
            if(trtr .gt. 0.d0)then
            elseif(trtr .le. 0.d0)then
            else
              over=.true.
              go to 9000
            endif
            if(insmat)then
              if(ltyp .eq. icINS)then
                insmat=.false.
              endif
              go to 1010
            endif
          else
            bx0=twiss(ip1,mfitbx)
            by0=twiss(ip1,mfitby)
            if(bx0 .gt. 0.d0 .and. by0 .gt. 0.d0 .and.
     $           itgetfpe() .eq. 0)then
            else
              do j=ip1-1,ip0+la,-1
                if(twiss(j,mfitbx) .gt. 0.d0
     $               .and. twiss(j,mfitby) .gt. 0.d0)then
                  ip1=j
                  go to 12
                endif
              enddo
              ip1=ip0+la
 12           continue
              do concurrent (k=1:ntfun)
                twiss(ip1+1:ip0+lb,k)=twiss(ip1,k)
              enddo
              over=.true.
              go to 9000
            endif
            ax0=twiss(ip1,mfitax)
            ay0=twiss(ip1,mfitay)
          endif
          call tfbndsol(l1,ibg,ibb)
          if(ibg .gt. 0)then
            call qsol(trans,cod,l1,coup)
            go to 20
          endif
          le=elatt%comp(l1)
          ld=idvalc(l1)
          lp=merge(ld,le,ideal)
          call loc_comp(lp,cmp)
          seg=tcheckseg(cmp,ltyp,al,lsegp,irtc)
c          write(*,*)'qtwiss ',seg,ltyp,al,irtc
c          if(seg)then
c            write(*,*)'qtwiss-seg ',cmp%value(1),cmp%value(ky_K1_MULT)
c          endif
          if(irtc .ne. 0)then
            call tffserrorhandle(l,irtc)
            go to 1010
          endif
          dir=direlc(l1)
          select case (ltyp)

          case (icDRFT)
            pr=1.d0+cod(6)
            pxi=cod(2)
            pyi=cod(4)
            pxisq=pxi**2
            pyisq=pyi**2
            a=pxisq+pyisq
            dpz=-a/pr/(1.d0+sqrt(max(0.d0,1.d0-a/pr**2)))
            pzi=pr+dpz
            ale=al/pzi
            alz=ale/pzi**2
            a12=ale+pxisq*alz
            a14=pxi*pyi*alz
            a15=-pr*pxi*alz
            a32=a14
            a34=ale+pyisq*alz
            a35=-pr*pyi*alz
            cod(1)=cod(1)+pxi*ale
            cod(3)=cod(3)+pyi*ale
            cod(5)=cod(5)+dpz*ale-dvemit*al
            if(mat .and. .not. wspac)then
              tr(1,1)=tr(1,1)+a12*tr(2,1)+a14*tr(4,1)
              tr(1,2)=tr(1,2)+a12*tr(2,2)+a14*tr(4,2)
              tr(1,3)=tr(1,3)+a12*tr(2,3)+a14*tr(4,3)
              tr(1,4)=tr(1,4)+a12*tr(2,4)+a14*tr(4,4)
              tr(1,5)=tr(1,5)+a12*tr(2,5)+a14*tr(4,5)+a15
              tr(3,1)=tr(3,1)+a32*tr(2,1)+a34*tr(4,1)
              tr(3,2)=tr(3,2)+a32*tr(2,2)+a34*tr(4,2)
              tr(3,3)=tr(3,3)+a32*tr(2,3)+a34*tr(4,3)
              tr(3,4)=tr(3,4)+a32*tr(2,4)+a34*tr(4,4)
              tr(3,5)=tr(3,5)+a32*tr(2,5)+a34*tr(4,5)+a35
              go to 10
            else
              coup=a14 .ne. 0.d0
              a11=1.d0
              a13=0.d0
              a21=0.d0
              a22=1.d0
              a23=0.d0
              a24=0.d0
              a25=0.d0
              a31=0.d0
              a33=1.d0
              a41=0.d0
              a42=0.d0
              a43=0.d0
              a44=1.d0
              a45=0.d0
              go to 20
            endif

          case (icBEND)
            if(dir .gt. 0.d0)then
              psi1=cmp%value(ky_E1_BEND)
              psi2=cmp%value(ky_E2_BEND)
              apsi1=cmp%value(ky_AE1_BEND)
              apsi2=cmp%value(ky_AE2_BEND)
              fb1=cmp%value(ky_F1_BEND)
     $             +cmp%value(ky_FB1_BEND)
              fb2=cmp%value(ky_F1_BEND)
     $             +cmp%value(ky_FB2_BEND)
            else
              psi1=cmp%value(ky_E2_BEND)
              psi2=cmp%value(ky_E1_BEND)
              apsi1=cmp%value(ky_AE2_BEND)
              apsi2=cmp%value(ky_AE1_BEND)
              fb2=cmp%value(ky_F1_BEND)
     $             +cmp%value(ky_FB1_BEND)
              fb1=cmp%value(ky_F1_BEND)
     $             +cmp%value(ky_FB2_BEND)
            endif
            dtheta=cmp%value(ky_DROT_BEND)
            theta0=cmp%value(ky_ROT_BEND)
            cod1=cod
            call qbend(trans,cod,al,
     $           cmp%value(ky_ANGL_BEND)+cmp%value(ky_K0_BEND),
     1           cmp%value(ky_ANGL_BEND),psi1,psi2,apsi1,apsi2,
     $           cmp%value(ky_K1_BEND),
     1           cmp%value(ky_DX_BEND),
     $           cmp%value(ky_DY_BEND),
     $           theta0,dtheta,
     $           fb1,fb2,
     $           nint(cmp%value(ky_FRMD_BEND)),
     $           cmp%value(ky_FRIN_BEND) .eq. 0.d0,
     $           cmp%value(ky_EPS_BEND),
     1           coup)
            go to 20

          case (icQUAD)
            mfr=nint(cmp%value(ky_FRMD_QUAD))
            if(dir .lt. 0.d0)then
              mfr=mfr*(11+mfr*(2*mfr-9))/2
            endif
            ak1=cmp%value(ky_K1_QUAD)
            call tsetfringepe(cmp,icQUAD,ftable)
            call qquad(trans,cod,al,
     1           ak1,cmp%value(ky_DX_QUAD),cmp%value(ky_DY_QUAD),
     1           cmp%value(ky_ROT_QUAD),
     $           cmp%value(ky_FRIN_QUAD) .eq. 0.d0,
     $           ftable(1),ftable(2),ftable(3),ftable(4),
     $           mfr,cmp%value(ky_EPS_QUAD),
     $           cmp%value(ky_KIN_QUAD) .eq. 0.d0,
     $           cmp%value(ky_CHRO_QUAD) .ne. 0.d0,
     $           coup)
            go to 20

          case (icSEXT, icOCTU, icDECA, icDODECA)
            call qthin(trans,cod,ltyp,al,cmp%value(ky_K_THIN),
     1           cmp%value(ky_DX_THIN),cmp%value(ky_DY_THIN),
     $           cmp%value(ky_ROT_THIN),coup)
          go to 20

          case (icSOL)
            write(*,*)'Qtwiss: implementation error of solenoid ',l1
          go to 1010

          case (icST)
            write(*,*)'Use BEND with ANGLE=0 for STEER.'
          call abort

          case (icMULT)
            bz=0.d0
            if(seg)then
c          call tfevals('Print["PROF-QT0: ",LINE["PROFILE","Q1"]]',
c     $             kxx,irtc)
              call qmultseg(trans,cod,l1,cmp,lsegp,bz,coup)
c          call tfevals('Print["PROF-QT1: ",LINE["PROFILE","Q1"]]',
c     $             kxx,irtc)
            else
              call qmult1(trans,cod,l1,cmp,bz,.true.,coup)
            endif
            go to 20

          case (icTEST)
            call qtest(trans,cod,al,cmp%value(ky_ANGL_TEST),coup)
          go to 20

          case (icCAVI)
            mfr=nint(cmp%value(ky_FRMD_CAVI))
            if(direlc(l1) .ge. 0.d0)then
            else
              mfr=mfr*(11+mfr*(2*mfr-9))/2
            endif
            call qcav(trans,cod,l1,al,
     $           cmp%value(ky_VOLT_CAVI)+cmp%value(ky_DVOLT_CAVI),
     $           cmp%value(ky_HARM_CAVI),
     $           cmp%value(ky_PHI_CAVI),cmp%value(ky_FREQ_CAVI),
     $           cmp%value(ky_DX_CAVI),cmp%value(ky_DY_CAVI),
     $           cmp%value(ky_ROT_CAVI),
     $           cmp%value(ky_V1_CAVI),cmp%value(ky_V20_CAVI),
     $           cmp%value(ky_V11_CAVI),cmp%value(ky_V02_CAVI),
     $           cmp%value(ky_FRIN_CAVI) .eq. 0.d0,mfr,
     $           cmp%value(ky_APHI_CAVI) .ne. 0.d0,
     $           coup)
            go to 20

          case (icTCAV)
            call qtcav(trans,cod,
     $         al,cmp%value(ky_K0_TCAV),cmp%value(ky_HARM_TCAV),
     $         cmp%value(ky_PHI_TCAV),cmp%value(ky_FREQ_TCAV),
     $         cmp%value(ky_DX_TCAV),cmp%value(ky_DY_TCAV),
     $         cmp%value(ky_ROT_TCAV),coup)
            go to 20

          case (icMAP)
            call qemap(trans6,cod,l1,coup,err)
            if(err)then
              go to 1010
            endif
            call qcopymat(trans,trans6,.false.)
            go to 20

          case (icINS)
            call qins(trans,cod,l1,idp,
     $           cmp%value(ky_DIR_INS) .ge. 0.d0,
     $           cmp%value(1),cmp%value(ky_DIR_INS+1),coup,
     $           mat,insmat)
            if(insmat)then
              go to 1010
            endif
            if(coup)then
              go to 20
            else
              go to 10
            endif

          case (icCOORD)
            call qcoord(trans,cod,
     1           cmp%value(ky_DX_COORD),cmp%value(ky_DY_COORD),
     $           cmp%value(ky_DZ_COORD),cmp%value(ky_CHI1_COORD),
     $           cmp%value(ky_CHI2_COORD),cmp%value(ky_CHI3_COORD),
     1           cmp%value(ky_DIR_COORD) .eq. 0.d0,coup)
            go to 20

          case default
          end select
 1010     if(wspac)then
            if(l .eq. lb)then
              trans=0.d0
              trans(1,1)=1.d0
              trans(2,2)=1.d0
              trans(3,3)=1.d0
              trans(4,4)=1.d0
              coup=.false.
              go to 20
            endif
          endif
          if(.not. mat)then
            twiss(ip,1:ntfun)=twiss(ip1,1:ntfun)
          endif
          go to 10

 20       if(wspac)then
            sspc=(rlist(ifpos+l1)+rlist(ifpos+l1-1))*.5d0
            call qwspac(trans,cod,sspc-sspc0,
     $           rlist(ifsize+(l1-1)*21),coup,l1)
c            write(*,*)'qtwiss-qwsapc ',l1,ifsize,rlist(ifsize+(l1-1)*21)
            sspc0=sspc
            coup=.true.
          endif
          if(mat)then
            if(coup)then
              x =tr(1,1)
              px=tr(2,1)
              y =tr(3,1)
              tr(1,1)=a11*x+a12*px+a13*y+a14*tr(4,1)
              tr(2,1)=a21*x+a22*px+a23*y+a24*tr(4,1)
              tr(3,1)=a31*x+a32*px+a33*y+a34*tr(4,1)
              tr(4,1)=a41*x+a42*px+a43*y+a44*tr(4,1)
              x =tr(1,2)
              px=tr(2,2)
              y =tr(3,2)
              tr(1,2)=a11*x+a12*px+a13*y+a14*tr(4,2)
              tr(2,2)=a21*x+a22*px+a23*y+a24*tr(4,2)
              tr(3,2)=a31*x+a32*px+a33*y+a34*tr(4,2)
              tr(4,2)=a41*x+a42*px+a43*y+a44*tr(4,2)
              x =tr(1,3)
              px=tr(2,3)
              y =tr(3,3)
              tr(1,3)=a11*x+a12*px+a13*y+a14*tr(4,3)
              tr(2,3)=a21*x+a22*px+a23*y+a24*tr(4,3)
              tr(3,3)=a31*x+a32*px+a33*y+a34*tr(4,3)
              tr(4,3)=a41*x+a42*px+a43*y+a44*tr(4,3)
              x =tr(1,4)
              px=tr(2,4)
              y =tr(3,4)
              tr(1,4)=a11*x+a12*px+a13*y+a14*tr(4,4)
              tr(2,4)=a21*x+a22*px+a23*y+a24*tr(4,4)
              tr(3,4)=a31*x+a32*px+a33*y+a34*tr(4,4)
              tr(4,4)=a41*x+a42*px+a43*y+a44*tr(4,4)
              x =tr(1,5)
              px=tr(2,5)
              y =tr(3,5)
              tr(1,5)=a11*x+a12*px+a13*y+a14*tr(4,5)+a15
              tr(2,5)=a21*x+a22*px+a23*y+a24*tr(4,5)+a25
              tr(3,5)=a31*x+a32*px+a33*y+a34*tr(4,5)+a35
              tr(4,5)=a41*x+a42*px+a43*y+a44*tr(4,5)+a45
            else
              x =tr(1,1)
              y =tr(3,1)
              tr(1,1)=a11*x+a12*tr(2,1)
              tr(2,1)=a21*x+a22*tr(2,1)
              tr(3,1)=a33*y+a34*tr(4,1)
              tr(4,1)=a43*y+a44*tr(4,1)
              x =tr(1,2)
              y =tr(3,2)
              tr(1,2)=a11*x+a12*tr(2,2)
              tr(2,2)=a21*x+a22*tr(2,2)
              tr(3,2)=a33*y+a34*tr(4,2)
              tr(4,2)=a43*y+a44*tr(4,2)
              x =tr(1,3)
              y =tr(3,3)
              tr(1,3)=a11*x+a12*tr(2,3)
              tr(2,3)=a21*x+a22*tr(2,3)
              tr(3,3)=a33*y+a34*tr(4,3)
              tr(4,3)=a43*y+a44*tr(4,3)
              x =tr(1,4)
              y =tr(3,4)
              tr(1,4)=a11*x+a12*tr(2,4)
              tr(2,4)=a21*x+a22*tr(2,4)
              tr(3,4)=a33*y+a34*tr(4,4)
              tr(4,4)=a43*y+a44*tr(4,4)
              x =tr(1,5)
              y =tr(3,5)
              tr(1,5)=a11*x+a12*tr(2,5)+a15
              tr(2,5)=a21*x+a22*tr(2,5)+a25
              tr(3,5)=a33*y+a34*tr(4,5)+a35
              tr(4,5)=a43*y+a44*tr(4,5)+a45
            endif
          else
            call qmat2twiss(trans,ip,l,twiss,dpsix,dpsiy,coup,normal)
          endif
          if(.not. mat)then
            if(orbitcal)then
              twiss(ip,mfitdx:mfitddp)=cod
            endif
            twiss(ip,mfitnx)=twiss(ip1,mfitnx)+dpsix
            twiss(ip,mfitny)=twiss(ip1,mfitny)+dpsiy
          endif
        endif
 10     continue
        call limitcod(cod)
      enddo
 9000 calpol=calpol0
      trf0=trf00
      vcalpha=vcalpha0
      return
      end

      subroutine qmat2twiss(trans,ip,l,
     $     twiss,dpsix,dpsiy,coup,normal)
      use tfstk
      use ffs
      use tffitcode
      use ffs_pointer, only:gammab
      implicit none
      integer*4 ip1,ip,l1,l
      real*8 trans(4,5),twiss(nlat*(2*ndim+1),ntwissfun),dpsix,dpsiy,
     $     ax0,bx0,ay0,by0
      real*8 r1,r2,r3,r4,sqrdet,detp,epx0,epy0,ex0,ey0
      real*8 u11,u12,u13,u14,u21,u22,u23,u24,u31,u32,u33,
     $     u34,u41,u42,u43,u44,u15,u25,u35,u45,
     $     r15,r25,r35,r45,aa,bb,cc,dd,detm,gr,
     $     t1,t2,t3,t4
      logical*4 coup,normal
      ip1=ip-1
      l1=l-1
      twiss(ip,mfitaz:mfitzpy)=twiss(ip1,mfitaz:mfitzpy)
      ax0=twiss(ip1,mfitax)
      bx0=twiss(ip1,mfitbx)
      ay0=twiss(ip1,mfitay)
      by0=twiss(ip1,mfitby)
      if(coup)then
        r1=twiss(ip1,mfitr1)
        r2=twiss(ip1,mfitr2)
        r3=twiss(ip1,mfitr3)
        r4=twiss(ip1,mfitr4)
        sqrdet=sqrt(1.d0-r1*r4+r2*r3)
        if(normal)then
          r15=sqrdet*twiss(ip1,mfitex)
     1         +r4*twiss(ip1,mfitey)-r2*twiss(ip1,mfitepy)
          r25=sqrdet*twiss(ip1,mfitepx)
     1         -r3*twiss(ip1,mfitey)+r1*twiss(ip1,mfitepy)
          r35=sqrdet*twiss(ip1,mfitey)
     1         -r1*twiss(ip1,mfitex)-r2*twiss(ip1,mfitepx)
          r45=sqrdet*twiss(ip1,mfitepy)
     1         -r3*twiss(ip1,mfitex)-r4*twiss(ip1,mfitepx)
          u11= trans(1,1)*sqrdet-trans(1,3)*r1-trans(1,4)*r3
          u12= trans(1,2)*sqrdet-trans(1,3)*r2-trans(1,4)*r4
          u13= trans(1,1)*r4-trans(1,2)*r3+trans(1,3)*sqrdet
          u14=-trans(1,1)*r2+trans(1,2)*r1+trans(1,4)*sqrdet
          u15= trans(1,1)*r15+trans(1,2)*r25+trans(1,3)*r35
     $         +trans(1,4)*r45+trans(1,5)
          u21= trans(2,1)*sqrdet-trans(2,3)*r1-trans(2,4)*r3
          u22= trans(2,2)*sqrdet-trans(2,3)*r2-trans(2,4)*r4
          u23= trans(2,1)*r4-trans(2,2)*r3+trans(2,3)*sqrdet
          u24=-trans(2,1)*r2+trans(2,2)*r1+trans(2,4)*sqrdet
          u25= trans(2,1)*r15+trans(2,2)*r25+trans(2,3)*r35
     $         +trans(2,4)*r45+trans(2,5)
          u31= trans(3,1)*sqrdet-trans(3,3)*r1-trans(3,4)*r3
          u32= trans(3,2)*sqrdet-trans(3,3)*r2-trans(3,4)*r4
          u33= trans(3,1)*r4-trans(3,2)*r3+trans(3,3)*sqrdet
          u34=-trans(3,1)*r2+trans(3,2)*r1+trans(3,4)*sqrdet
          u35= trans(3,1)*r15+trans(3,2)*r25+trans(3,3)*r35
     $         +trans(3,4)*r45+trans(3,5)
          u41= trans(4,1)*sqrdet-trans(4,3)*r1-trans(4,4)*r3
          u42= trans(4,2)*sqrdet-trans(4,3)*r2-trans(4,4)*r4
          u43= trans(4,1)*r4-trans(4,2)*r3+trans(4,3)*sqrdet
          u44=-trans(4,1)*r2+trans(4,2)*r1+trans(4,4)*sqrdet
          u45= trans(4,1)*r15+trans(4,2)*r25+trans(4,3)*r35
     $         +trans(4,4)*r45+trans(4,5)
        else
          r15=-r1*twiss(ip1,mfitex)-r2*twiss(ip1,mfitepx)
     1         +sqrdet*twiss(ip1,mfitey)
          r25=-r3*twiss(ip1,mfitex)-r4*twiss(ip1,mfitepx)
     1         +sqrdet*twiss(ip1,mfitepy)
          r35=sqrdet*twiss(ip1,mfitex)
     1         +r4*twiss(ip1,mfitey)-r2*twiss(ip1,mfitepy)
          r45=sqrdet*twiss(ip1,mfitepx)
     1         -r3*twiss(ip1,mfitey)+r1*twiss(ip1,mfitepy)
          u11=-trans(1,1)*r1-trans(1,2)*r3+trans(1,3)*sqrdet
          u12=-trans(1,1)*r2-trans(1,2)*r4+trans(1,4)*sqrdet
          u13= trans(1,1)*sqrdet+trans(1,3)*r4-trans(1,4)*r3
          u14= trans(1,2)*sqrdet-trans(1,3)*r2+trans(1,4)*r1
          u15= trans(1,1)*r15+trans(1,2)*r25+trans(1,3)*r35
     $         +trans(1,4)*r45+trans(1,5)
          u21=-trans(2,1)*r1-trans(2,2)*r3+trans(2,3)*sqrdet
          u22=-trans(2,1)*r2-trans(2,2)*r4+trans(2,4)*sqrdet
          u23= trans(2,1)*sqrdet+trans(2,3)*r4-trans(2,4)*r3
          u24= trans(2,2)*sqrdet-trans(2,3)*r2+trans(2,4)*r1
          u25= trans(2,1)*r15+trans(2,2)*r25+trans(2,3)*r35
     $         +trans(2,4)*r45+trans(2,5)
          u31=-trans(3,1)*r1-trans(3,2)*r3+trans(3,3)*sqrdet
          u32=-trans(3,1)*r2-trans(3,2)*r4+trans(3,4)*sqrdet
          u33= trans(3,1)*sqrdet+trans(3,3)*r4-trans(3,4)*r3
          u34= trans(3,2)*sqrdet-trans(3,3)*r2+trans(3,4)*r1
          u35= trans(3,1)*r15+trans(3,2)*r25+trans(3,3)*r35
     $         +trans(3,4)*r45+trans(3,5)
          u41=-trans(4,1)*r1-trans(4,2)*r3+trans(4,3)*sqrdet
          u42=-trans(4,1)*r2-trans(4,2)*r4+trans(4,4)*sqrdet
          u43= trans(4,1)*sqrdet+trans(4,3)*r4-trans(4,4)*r3
          u44= trans(4,2)*sqrdet-trans(4,3)*r2+trans(4,4)*r1
          u45= trans(4,1)*r15+trans(4,2)*r25+trans(4,3)*r35
     $         +trans(4,4)*r45+trans(4,5)
        endif
        detp=(u11*u22-u21*u12+u33*u44-u43*u34)*.5d0
        normal=detp .gt. xyth
        if(normal)then
          sqrdet=sqrt(detp)
          u11=u11/sqrdet
          u12=u12/sqrdet
          u21=u21/sqrdet
          u22=u22/sqrdet
          u33=u33/sqrdet
          u34=u34/sqrdet
          u43=u43/sqrdet
          u44=u44/sqrdet
          r1=-u31*u22+u32*u21
          r2= u31*u12-u32*u11
          r3=-u41*u22+u42*u21
          r4= u41*u12-u42*u11
          twiss(ip,mfitr1)=r1
          twiss(ip,mfitr2)=r2
          twiss(ip,mfitr3)=r3
          twiss(ip,mfitr4)=r4
          twiss(ip,mfitdetr)=r1*r4-r2*r3
          sqrdet=sqrt(1.d0-twiss(ip,mfitdetr))
c     write(*,'(a,1p6g15.7)')
c     $               'qtwiss-coup-n  ',r1,r2,r3,r4,r1*r4-r2*r3,detp
          twiss(ip,mfitex)=u15*sqrdet
     1         -twiss(ip,mfitr4)*u35+twiss(ip,mfitr2)*u45
          twiss(ip,mfitepx)=u25*sqrdet
     1         +twiss(ip,mfitr3)*u35-twiss(ip,mfitr1)*u45
          twiss(ip,mfitey)=u35*sqrdet+
     1         twiss(ip,mfitr1)*u15+twiss(ip,mfitr2)*u25
          twiss(ip,mfitepy)=u45*sqrdet+
     1         twiss(ip,mfitr3)*u15+twiss(ip,mfitr4)*u25
          aa=bx0*u11-ax0*u12
          cc=(u12-ax0*aa)/bx0
          twiss(ip,mfitax) =-(u21*aa+u22*cc)
          twiss(ip,mfitbx) =u11*aa+u12*cc
          dpsix=atan2(u12,aa)
          bb=by0*u33-ay0*u34
          dd=(u34-ay0*bb)/by0
          twiss(ip,mfitay) =-(u43*bb+u44*dd)
          twiss(ip,mfitby) =u33*bb+u34*dd
          dpsiy=atan2(u34,bb)
        else
          detm=(u31*u42-u32*u41+u13*u24-u14*u23)*.5d0
          sqrdet=sqrt(abs(detm))
          u31=u31/sqrdet
          u32=u32/sqrdet
          u41=u41/sqrdet
          u42=u42/sqrdet
          u13=u13/sqrdet
          u14=u14/sqrdet
          u23=u23/sqrdet
          u24=u24/sqrdet
          r1=-u11*u42+u12*u41
          r2= u11*u32-u12*u31
          r3=-u21*u42+u22*u41
          r4= u21*u32-u22*u31
          sqrdet=sqrt(1.d0-r1*r4+r2*r3)
          twiss(ip,mfitex) =u35*sqrdet-r4*u15+r2*u25
          twiss(ip,mfitepx)=u45*sqrdet+r3*u15-r1*u25
          twiss(ip,mfitey) =u15*sqrdet+r1*u35+r2*u45
          twiss(ip,mfitepy)=u25*sqrdet+r3*u35+r4*u45
          twiss(ip,mfitr1)=r1
          twiss(ip,mfitr2)=r2
          twiss(ip,mfitr3)=r3
          twiss(ip,mfitr4)=r4
          twiss(ip,mfitdetr)=1.d0+xyth-r1*r4+r2*r3
          aa=bx0*u31-ax0*u32
          cc=(u32-ax0*aa)/bx0
          twiss(ip,mfitax) =-(u41*aa+u42*cc)
          twiss(ip,mfitbx) =u31*aa+u32*cc
          dpsix=atan2(u32,aa)
          bb=by0*u13-ay0*u14
          dd=(u14-ay0*bb)/by0
          twiss(ip,mfitay) =-(u23*bb+u24*dd)
          twiss(ip,mfitby) =u13*bb+u14*dd
          dpsiy=atan2(u14,bb)
        endif
      else
        if(trpt)then
          gr=gammab(l)/gammab(l1)
          ex0 =twiss(ip1,mfitex)*gr
          epx0=twiss(ip1,mfitepx)*gr
          ey0 =twiss(ip1,mfitey)*gr
          epy0=twiss(ip1,mfitepy)*gr
        else
          ex0 =twiss(ip1,mfitex)
          epx0=twiss(ip1,mfitepx)
          ey0 =twiss(ip1,mfitey)
          epy0=twiss(ip1,mfitepy)
        endif
        if(normal)then
          t1= twiss(ip1,mfitr1)*trans(2,2)
     $         -twiss(ip1,mfitr2)*trans(2,1)
          t2=-twiss(ip1,mfitr1)*trans(1,2)
     $         +twiss(ip1,mfitr2)*trans(1,1)
          t3= twiss(ip1,mfitr3)*trans(2,2)
     $         -twiss(ip1,mfitr4)*trans(2,1)
          t4=-twiss(ip1,mfitr3)*trans(1,2)
     $         +twiss(ip1,mfitr4)*trans(1,1)
          r1= trans(3,3)*t1+trans(3,4)*t3
          r2= trans(3,3)*t2+trans(3,4)*t4
          r3= trans(4,3)*t1+trans(4,4)*t3
          r4= trans(4,3)*t2+trans(4,4)*t4
          twiss(ip,mfitr1)=r1
          twiss(ip,mfitr2)=r2
          twiss(ip,mfitr3)=r3
          twiss(ip,mfitr4)=r4
          twiss(ip,mfitdetr)=r1*r4-r2*r3
          sqrdet=sqrt(1.d0-twiss(ip,mfitdetr))
          twiss(ip,mfitex) =
     1         trans(1,1)*ex0+trans(1,2)*epx0+sqrdet*trans(1,5)
     $         -r4*trans(3,5)+r2*trans(4,5)
          twiss(ip,mfitepx) =
     1         trans(2,1)*ex0+trans(2,2)*epx0+sqrdet*trans(2,5)
     $         +r3*trans(3,5)-r1*trans(4,5)
          twiss(ip,mfitey) =
     1         trans(3,3)*ey0+trans(3,4)*epy0+sqrdet*trans(3,5)
     $         +r1*trans(1,5)+r2*trans(2,5)
          twiss(ip,mfitepy)=
     1         trans(4,3)*ey0+trans(4,4)*epy0+sqrdet*trans(4,5)
     $         +r3*trans(1,5)+r4*trans(2,5)
          aa=bx0*trans(1,1)-ax0*trans(1,2)
          cc=(trans(1,2)-ax0*aa)/bx0
          twiss(ip,mfitax) =-(trans(2,1)*aa+trans(2,2)*cc)
          twiss(ip,mfitbx) =trans(1,1)*aa+trans(1,2)*cc
          dpsix=atan2(trans(1,2),aa)
          bb=by0*trans(3,3)-ay0*trans(3,4)
          dd=(trans(3,4)-ay0*bb)/by0
          twiss(ip,mfitay) =-(trans(4,3)*bb+trans(4,4)*dd)
          twiss(ip,mfitby) =trans(3,3)*bb+trans(3,4)*dd
          dpsiy=atan2(trans(3,4),bb)
        else
          r1=twiss(ip1,mfitr1)
          r2=twiss(ip1,mfitr2)
          r3=twiss(ip1,mfitr3)
          r4=twiss(ip1,mfitr4)
          sqrdet=sqrt(1.d0-r1*r4+r2*r3)
          t1= r1*trans(4,4)-r2*trans(4,3)
          t2=-r1*trans(3,4)+r2*trans(3,3)
          t3= r3*trans(4,4)-r4*trans(4,3)
          t4=-r3*trans(3,4)+r4*trans(3,3)
          r1= trans(1,1)*t1+trans(1,2)*t3
          r2= trans(1,1)*t2+trans(1,2)*t4
          r3= trans(2,1)*t1+trans(2,2)*t3
          r4= trans(2,1)*t2+trans(2,2)*t4
          twiss(ip,mfitr1)=r1
          twiss(ip,mfitr2)=r2
          twiss(ip,mfitr3)=r3
          twiss(ip,mfitr4)=r4
          twiss(ip,mfitdetr)=1.d0+xyth-r1*r4+r2*r3
          twiss(ip,mfitex) =
     1         trans(3,3)*ex0+trans(3,4)*epx0+sqrdet*trans(3,5)
     $         +(-r4*trans(1,5)+r2*trans(2,5))
          twiss(ip,mfitepx) =
     1         trans(4,3)*ex0+trans(4,4)*epx0+sqrdet*trans(4,5)
     $         +( r3*trans(1,5)-r1*trans(2,5))
          twiss(ip,mfitey) =
     1         trans(1,1)*ey0+trans(1,2)*epy0+sqrdet*trans(1,5)
     $         +( r1*trans(3,5)+r2*trans(4,5))
          twiss(ip,mfitepy)=
     1         trans(2,1)*ey0+trans(2,2)*epy0+sqrdet*trans(2,5)
     $         +( r3*trans(3,5)+r4*trans(4,5))
          aa=bx0*trans(3,3)-ax0*trans(3,4)
          cc=(trans(3,4)-ax0*aa)/bx0
          twiss(ip,mfitax) =-(trans(4,3)*aa+trans(4,4)*cc)
          twiss(ip,mfitbx) =trans(3,3)*aa+trans(3,4)*cc
          dpsix=atan2(trans(3,4),aa)
          bb=by0*trans(1,1)-ay0*trans(1,2)
          dd=(trans(1,2)-ay0*bb)/by0
          twiss(ip,mfitay) =-(trans(2,1)*bb+trans(2,2)*dd)
          twiss(ip,mfitby) =trans(1,1)*bb+trans(1,2)*dd
          dpsiy=atan2(trans(1,2),bb)
        endif
      endif
      return
      end

      subroutine qtrans(la,lb,trans,cod,over)
      use tfstk
      use ffs
      use ffs_pointer
      use tffitcode
      implicit none
      type (ffs_bound) fbound
      integer*4 la,lb
      real*8 trans(4,5),cod(6)
      logical*4 over
      call tffsbound1(la,lb,fbound)
c      write(*,*)'qtrans ',la,lb,la1,lb1,fra,frb
      cod(5)=0.d0
      call qcod(1,fbound,trans,cod,.true.,over)
      return
      end

      subroutine qcod(idp,fbound,trans,cod0,codfnd,over)
      use tfstk
      use ffs
      use ffs_pointer
      use tffitcode
      use iso_c_binding
      implicit none
      type (ffs_bound) fbound
      real*8 conv,cx,sx,ax,bx,cy,sy,ay,by,r0,dcod(6)
      integer*4 , parameter :: itmax=30
      real*8 , parameter :: conv0=1.d-19,conv1=1.d-10,
     $     factmin=1.d-3
      integer*4 idp,it
      real*8 r,fact
      real*8 , pointer :: ptwiss(:,:)
      real*8 trans(4,5),cod(6),cod0(6),trans1(4,5),transb(4,5),
     $     transe(4,5),ftwiss(ntwissfun),trans2(4,5),cod00(6)
      logical*4 over,codfnd,stab
      it=0
      r0=1.d100
      cod00=cod0
      fact=.5d0
      conv=min(conv1,conv0*(1.d0+(cod0(6)/0.001d0)**2))
      stab=.false.
      call c_f_pointer(c_loc(rlist(iftwis)),
     $     ptwiss,[nlat*(2*ndim+1),ntwissfun])
      do while(it .le. itmax)
        cod=cod0
        if(fbound%fb .gt. 0.d0)then
          call qtwissfrac1(ftwiss,transb,cod,idp,
     $         fbound%lb,fbound%fb,1.d0,.true.,.true.,over)
          call qtwiss1(ptwiss,idp,fbound%lb+1,fbound%le,
     $         trans1,cod,.true.,over)
c          do i=1,5
            trans2(1,1:5)=
     $         trans1(1,1)*transb(1,1:5)+trans1(1,2)*transb(2,1:5)
     $        +trans1(1,3)*transb(3,1:5)+trans1(1,4)*transb(4,1:5)
            trans2(2,1:5)=
     $         trans1(2,1)*transb(1,1:5)+trans1(2,2)*transb(2,1:5)
     $        +trans1(2,3)*transb(3,1:5)+trans1(2,4)*transb(4,1:5)
            trans2(3,1:5)=
     $         trans1(3,1)*transb(1,1:5)+trans1(3,2)*transb(2,1:5)
     $        +trans1(3,3)*transb(3,1:5)+trans1(3,4)*transb(4,1:5)
            trans2(4,1:5)=
     $         trans1(4,1)*transb(1,1:5)+trans1(4,2)*transb(2,1:5)
     $        +trans1(4,3)*transb(3,1:5)+trans1(4,4)*transb(4,1:5)
c          enddo
          trans2(:,5)=trans2(:,5)+trans1(:,5)
        else
          call qtwiss1(ptwiss,idp,fbound%lb,fbound%le,
     $         trans2,cod,.true.,over)
        endif
        if(fbound%fe .gt. 0.d0)then
          call qtwissfrac1(ftwiss,transe,cod,idp,
     $         fbound%le,0.d0,fbound%fe,.true.,.true.,over)
c          do i=1,5
            trans(1,1:5)=
     $           transe(1,1)*trans2(1,1:5)+transe(1,2)*trans2(2,1:5)
     $          +transe(1,3)*trans2(3,1:5)+transe(1,4)*trans2(4,1:5)
            trans(2,1:5)=
     $           transe(2,1)*trans2(1,1:5)+transe(2,2)*trans2(2,1:5)
     $          +transe(2,3)*trans2(3,1:5)+transe(2,4)*trans2(4,1:5)
            trans(3,1:5)=
     $           transe(3,1)*trans2(1,1:5)+transe(3,2)*trans2(2,1:5)
     $          +transe(3,3)*trans2(3,1:5)+transe(3,4)*trans2(4,1:5)
            trans(4,1:5)=
     $           transe(4,1)*trans2(1,1:5)+transe(4,2)*trans2(2,1:5)
     $          +transe(4,3)*trans2(3,1:5)+transe(4,4)*trans2(4,1:5)
c          enddo
          trans(:,5)=trans(:,5)+transe(:,5)
        else
          trans=trans2
        endif
        call resetnan(cod,1.d300)
        if(.not. orbitcal)then
          codfnd=.true.
        endif
c        write(*,'(a,2l3,1p6g15.7)')'qcod ',codfnd,over,cod
        if(codfnd)then
          cod0=cod
          return
        endif
        if(over)then
          codfnd=.false.
          return
        endif
        cx=.5d0*(trans(1,1)+trans(2,2))
        cy=.5d0*(trans(3,3)+trans(4,4))
        if(abs(cx) .gt. 1.d0)then
c          if(stab)then
c            it=it+1
c            cod0=(2.d0*cod0+cod00)/3.d0
c            cycle
c          endif
          cx=1.d0/cx
        endif
        if(abs(cy) .gt. 1.d0)then
c          if(stab)then
c            it=it+1
c            cod0=(2.d0*cod0+cod00)/3.d0
c            cycle
c          endif
          cy=1.d0/cy
        else
          stab=abs(cx) .le. 1.d0
        endif
        sx=sign(sqrt(max(1.d-6,1.d0-cx**2)),trans(1,2))
        ax=.5d0*(trans(1,1)-trans(2,2))/sx
        bx=trans(1,2)/sx
        sy=sign(sqrt(max(1.d-6,1.d0-cy**2)),trans(3,4))
        ay=.5d0*(trans(3,3)-trans(4,4))/sy
        by=trans(3,4)/sy
        dcod=cod-cod0
        r=dcod(1)**2/bx+bx*(dcod(2)+ax/bx*dcod(1))**2
     $       +dcod(3)**2/by+by*(dcod(4)+ay/by*dcod(3))**2
        if(ktfenanq(r))then
          r=1.d300
        endif
        if(r .le. conv)then
          codfnd=.true.
          return
        endif
c        write(*,'(a,i5,1p7g14.6)')'qcod ',it,r,r0,fact,cod0(1:4)
        it=it+1
        if(r .gt. r0)then
          if(fact .lt. factmin)then
            fact=fact*16.d0
            cod0=(1.d0+fact)*cod00-fact*cod0
          else
            cod0=(1.d0-fact)*cod00+fact*cod0
            fact=abs(fact)*.5d0
          endif
        else
          fact=min(0.5d0,fact*2.d0)
          r0=r
          cod00=cod0
          cod(1)=-cod(1)+trans(1,1)*cod0(1)+trans(1,2)*cod0(2)
     $         +trans(1,3)*cod0(3)+trans(1,4)*cod0(4)
          cod(2)=-cod(2)+trans(2,1)*cod0(1)+trans(2,2)*cod0(2)
     $         +trans(2,3)*cod0(3)+trans(2,4)*cod0(4)
          cod(3)=-cod(3)+trans(3,1)*cod0(1)+trans(3,2)*cod0(2)
     $         +trans(3,3)*cod0(3)+trans(3,4)*cod0(4)
          cod(4)=-cod(4)+trans(4,1)*cod0(1)+trans(4,2)*cod0(2)
     $         +trans(4,3)*cod0(3)+trans(4,4)*cod0(4)
          trans1=trans
          trans1(1,1)=trans1(1,1)-1.d0
          trans1(2,2)=trans1(2,2)-1.d0
          trans1(3,3)=trans1(3,3)-1.d0
          trans1(4,4)=trans1(4,4)-1.d0
          call tsolvg(trans1,cod,cod0,4,4,4)
        endif
      enddo
      cod=cod0
      codfnd=.false.
      return
      end

      subroutine qtwissfrac(ftwiss,l,fr,over)
      use ffs
      implicit none
      integer*4 , intent(in)::l
      real*8 , intent(out)::ftwiss(ntwissfun)
      real*8 , intent(in)::fr
      real*8 gv(3,4)
      logical*4 , intent(out)::over
      call qtwissfracgeo(ftwiss,gv,l,fr,.false.,over)
      return
      end

      subroutine qtwissfracgeo(ftwiss,gv,l,fr,cgeo,over)
      use tfstk
      use ffs
      use ffs_pointer
      use tffitcode
      use temw, only:etwiss2ri,tfetwiss,tinv6,iaez
      use geolib
      implicit none
      type (sad_descriptor) dsave(kwMAX)
      type (sad_comp) , pointer :: cmp
      integer*4 , intent(in)::l
      integer*4 nvar,le,itfdownlevel,irtc
      real*8 , intent(in)::fr
      real*8 , intent(out)::ftwiss(ntwissfun),gv(3,4)
      real*8 trans(6,6),cod(6),gr,sgr,sgr2,gr1,
     $     tw1(ntwissfun),beam(21),srot(3,9)
      logical*4 , intent(in)::cgeo
      logical*4 , intent(out)::over
      logical*4 sol,rt,chg,cp0,normal
      if(calc6d)then
        cp0=codplt
        codplt=.false.
        rt=radtaper
        sol=.false.
        irad=6
        levele=levele+1
        call qfracsave(l,dsave,nvar,.true.)
        call compelc(l,cmp)
        call qfracseg(cmp,cmp,0.d0,fr,chg,irtc)
        if(irtc .ne. 0)then
          call tffserrorhandle(l,irtc)
        else
          tw1=twiss(l,0,1:ntwissfun)
          cod=tw1(mfitdx:mfitddp)
          trans=tinv6(etwiss2ri(tw1,normal))
          call tturne1(trans,cod,beam,srot,
     $         iaez,0,
     $         .false.,sol,rt,.true.,l,l)
        endif
        if(cgeo)then
          gv=tfgeo1s(l)
        endif
        if(chg)then
          call qfracsave(l,dsave,nvar,.false.)
        endif
        le=itfdownlevel()
        if(trpt)then
          gr=gammab(l+1)/gammab(l)
          sgr2=1.d0+(gr-1.d0)*fr
          sgr=sqrt(sgr2)
          gr1=gr/sgr
          trans(1,:)=trans(1,:)*sgr
          trans(3,:)=trans(3,:)*sgr
          trans(5,:)=trans(5,:)*sgr
          trans(2,:)=trans(2,:)*gr1
          trans(4,:)=trans(4,:)*gr1
          trans(6,:)=trans(6,:)*gr1
          cod(2)=cod(2)*gr/sgr2
          cod(4)=cod(4)*gr/sgr2
          cod(6)=(1.d0+cod(6))*gr/sgr2-1.d0
c          write(*,'(a,1p8g15.7)')'qtwissfrac ',fr,sgr2,cod
        endif
c        write(*,'(1p6g15.7)')(trans(i,1:6),i=1,6)
c        call tinv6(trans,ri)
        ftwiss=tfetwiss(tinv6(trans),cod,normal)
        ftwiss(mfitnx)=ftwiss(mfitnx)+twiss(l,0,mfitnx)
        ftwiss(mfitny)=ftwiss(mfitny)+twiss(l,0,mfitny)
        ftwiss(mfitnz)=ftwiss(mfitnz)+twiss(l,0,mfitnz)
c        write(*,'(a,i5,1p8g14.6)')'qtwissfrac ',l,fr,gr,ftwiss(1:mfitny)
        over=.false.
        codplt=cp0
      else
        call qtwissfrac1geo(ftwiss,gv,trans,cod,
     $       0,l,0.d0,fr,.false.,.true.,.false.,over)
      endif
      return
      end

      subroutine qtwissfrac1(ftwiss,
     $     trans,cod,idp,l,fr1,fr2,mat,force,over)
      use ffs
      implicit none
      integer*4 idp,l
      real*8 , intent(out)::ftwiss(ntwissfun)
      real*8 gv(3,4),trans(4,5),cod(6),fr1,fr2
      logical*4 over,mat,force
      call qtwissfrac1geo(ftwiss,gv,
     $     trans,cod,idp,l,fr1,fr2,mat,.false.,force,over)
      return
      end

      subroutine qtwissfrac1geo(ftwiss,gv,
     $     trans,cod,idp,l,fr1,fr2,mat,cgeo,force,over)
      use tfstk
      use ffs
      use ffs_pointer
      use tffitcode
      use geolib, only:tfgeo1s
      implicit none
      type (sad_descriptor) dsave(kwMAX)
      type (sad_comp) ,pointer :: cmp,cmp0
      integer*4 ,intent(in):: idp,l
      integer*4 nvar,le,itfdownlevel,irtc
      real*8 , intent(out)::ftwiss(ntwissfun),gv(3,4),
     $     trans(4,5),cod(6)
      real*8 ,intent(in):: fr1,fr2
      real*8 twisss(ntwissfun),gb0,gb1,dgb
      logical*4 , intent(in):: mat,force,cgeo
      logical*4 , intent(out):: over
      logical*4 chg
      levele=levele+1
      call qfracsave(l,dsave,nvar,.true.)
      call compelc(l,cmp)
      cmp0=>cmp
      if(ideal)then
        call loc_comp(idval(cmp%id),cmp)
      endif
      call qfracseg(cmp,cmp0,fr1,fr2,chg,irtc)
      gb0=gammab(l)
      gb1=gammab(l+1)
      dgb=gb1-gb0
      if(chg .or. force)then
        gammab(l)=gb0+fr1*dgb
        gammab(l+1)=gb0+fr2*dgb
        pgev=gammab(l)*amass
        call tphyzp
        if(mat)then
          call qtwiss1(twiss,idp,l,l+1,
     $         trans,cod,.true.,over)
        elseif(force)then
          call qtwiss(twiss,idp,l,l+1,over)
        else
          twisss(1:ntwissfun)=twiss(l+1,idp,1:ntwissfun)
          call qtwiss(twiss,idp,l,l+1,over)
          ftwiss(1:ntwissfun)=twiss(l+1,idp,1:ntwissfun)
          twiss(l+1,idp,1:ntwissfun)=twisss(1:ntwissfun)
        endif
        gammab(l)=gb0
        gammab(l+1)=gb1
        if(cgeo)then
          gv=tfgeo1s(l)
        endif
        if(chg)then
          call qfracsave(l,dsave,nvar,.false.)
        endif
      else
        ftwiss(1:ntwissfun)=twiss(l+1,idp,1:ntwissfun)
      endif
      le=itfdownlevel()
      return
      end

      subroutine qfracsave(l,dsave,nvar,save)
      use tfstk
      use mackw
      use ffs_flag
      use ffs_pointer, only:idelc,elatt,idtypec,idvalc
      use kyparam
      implicit none
      type (sad_descriptor) dsave(*)
      integer*4 nvar,l,lt
      integer*8 i
      logical*4 save
      i=merge(idvalc(l),elatt%comp(l),ideal)
      lt=idtypec(l)
      if(save)then
        nvar=kytbl(kwMAX,lt)-1
        dsave(1:nvar)=dlist(i+1:i+nvar)
        if(lt .eq. icMULT)then
          dsave(nvar+1)=dlist(i+p_PROF_MULT)
        endif
      else
        dlist(i+1:nvar)=dsave(1:nvar)
        if(lt .eq. icMULT)then
          dlist(i+p_PROF_MULT)=dsave(nvar+1)
        endif
      endif
      return
      end

      subroutine qfraccomp(cmp,rx1,rx2,ideal,chg)
      use kyparam
      use tfstk
      use sad_main
      use ffs_pointer, only:idelc,direlc,idtypec,idvalc
      use mackw
      implicit none
      type (sad_comp) :: cmp
      integer*4 lt
      real*8 r,dl,rx1,rx2,f1,f2,fr0
      logical*4 ideal,chg
      chg=.false.
      r=rx2-rx1
      if(r .eq. 1.d0)then
        return
      endif
      f1=merge(1.d0,0.d0,rx1 .eq. 0.d0)
      f2=merge(1.d0,0.d0,rx2 .eq. 0.d0)
      lt=idtype(cmp%id)
      select case (lt)

      case (icDRFT)
        go to 9000
        
      case (icBEND)
        if(cmp%value(ky_FRMD_BEND) .eq. 0)then
          cmp%value(ky_F1_BEND)=0.d0
        endif
        cmp%value(ky_FRMD_BEND)=-f1-2.d0*f2
        if(r .ne. 0.d0)then
c          if(cmp%orient .gt. 0.d0)then
          if(cmp%ori)then
            cmp%value(ky_E1_BEND)=
     $           cmp%value(ky_E1_BEND)*f1/r
            cmp%value(ky_E2_BEND)=
     $           cmp%value(ky_E2_BEND)*f2/r
          else
            cmp%value(ky_E1_BEND)=
     $           cmp%value(ky_E1_BEND)*f2/r
            cmp%value(ky_E2_BEND)=
     $           cmp%value(ky_E2_BEND)*f1/r
          endif
        endif
        cmp%value(ky_ANGL_BEND)=
     $       cmp%value(ky_ANGL_BEND)*r
        cmp%value(ky_K1_BEND)=cmp%value(ky_K1_BEND)*r
        cmp%value(ky_K0_BEND)=cmp%value(ky_K0_BEND)*r
c     write(*,*)'qfraccomp ',r,
c     $     cmp%value(ky_ANGL_BEND),
c     $     cmp%value(ky_K1_BEND),
c     $     cmp%value(ky_K0_BEND)
        go to 9000

      case (icQUAD)
        cmp%value(ky_K1_QUAD)=cmp%value(ky_K1_QUAD)*r

      case (icSEXT,icOCTU,icDECA,icDODECA)
        cmp%value(ky_K_THIN)=cmp%value(ky_K_THIN)*r
        go to 9000

      case (icMULT)
        cmp%value(ky_K0_MULT:ky_SK21_MULT-1)=
     $       cmp%value(ky_K0_MULT:ky_SK21_MULT-1)*r
        dl=(1.d0-r)*cmp%value(ky_L_MULT)*.5d0
        cmp%value(ky_DX_MULT)=
     $       cmp%value(ky_DX_MULT)-dl*sin(cmp%value(ky_CHI1_MULT))
        cmp%value(ky_DY_MULT)=
     $       cmp%value(ky_DY_MULT)-dl*sin(cmp%value(ky_CHI2_MULT))
        cmp%value(ky_DZ_MULT)=cmp%value(ky_DZ_MULT)
     $       +dl*(2.d0*sin(.5d0*cmp%value(ky_CHI1_MULT))**2*
     $       cos(cmp%value(ky_CHI2_MULT)))
        cmp%value(ky_VOLT_MULT)=cmp%value(ky_VOLT_MULT)*r
        cmp%value(ky_DVOLT_MULT)=cmp%value(ky_DVOLT_MULT)*r
        cmp%value(ky_W1_MULT)=cmp%value(ky_W1_MULT)*r
        if(cmp%value(ky_ANGL_MULT) .ne. 0.d0)then
          if(cmp%value(ky_FRMD_MULT) .eq. 0.d0)then
            cmp%value(ky_FB1_MULT)=0.d0
            cmp%value(ky_FB2_MULT)=0.d0
          endif
          cmp%value(ky_FRMD_BEND)=-f1-2.d0*f2
          if(cmp%ori)then
            cmp%value(ky_E1_MULT)=
     $           cmp%value(ky_E1_MULT)*f1/r
            cmp%value(ky_E2_MULT)=
     $           cmp%value(ky_E2_MULT)*f2/r
          else
            cmp%value(ky_E1_MULT)=
     $           cmp%value(ky_E1_MULT)*f2/r
            cmp%value(ky_E2_MULT)=
     $           cmp%value(ky_E2_MULT)*f1/r
          endif
          cmp%value(ky_ANGL_MULT)=
     $         cmp%value(ky_ANGL_MULT)*r
          go to 9000
        endif

      case (icCAVI)
        cmp%value(ky_VOLT_CAVI)=cmp%value(ky_VOLT_CAVI)*r
        cmp%value(ky_DVOLT_CAVI)=cmp%value(ky_DVOLT_CAVI)*r
        cmp%value(ky_RANV_CAVI)=cmp%value(ky_RANV_CAVI)*r
        cmp%value(ky_V1_CAVI)=cmp%value(ky_V1_CAVI)*r

      case(icTCAV)
        cmp%value(ky_K0_TCAV)=cmp%value(ky_K0_TCAV)*r
        cmp%value(ky_RANK_TCAV)=cmp%value(ky_RANK_TCAV)*r
        go to 9000

      case default
        return
      end select
      fr0=cmp%value(kytbl(kwFRMD,lt))
      if(fr0 .eq. 0.d0 .or. fr0 .eq. 2)then
        if(kytbl(kwFB1,lt) .ne. 0)then
          cmp%value(kytbl(kwFB1,lt))=0.d0
        endif
      endif
      if(fr0 .eq. 0.d0 .or. fr0 .eq. 1)then
        if(kytbl(kwFB2,lt) .ne. 0)then
          cmp%value(kytbl(kwFB2,lt))=0.d0
        endif
      endif
      if(fr0 .eq. 0)then
        fr0=3.d0
      endif
      cmp%value(kytbl(kwFRMD,lt))=0.d0
      if(f1 .ne. 0.d0)then
        if(cmp%ori)then
          if(fr0 .eq. 3.d0 .or. fr0 .eq. 1.d0)then
            cmp%value(kytbl(kwFRMD,lt))=1.d0
          endif
        else
          if(fr0 .eq. 3.d0 .or. fr0 .eq. 2.d0)then
            cmp%value(kytbl(kwFRMD,lt))=2.d0
          endif
        endif
      endif
      if(f2 .ne. 0.d0)then
        if(cmp%ori)then
          if(fr0 .eq. 3.d0 .or. fr0 .eq. 2.d0)then
            cmp%value(kytbl(kwFRMD,lt))=cmp%value(kytbl(kwFRMD,lt))+2.d0
          endif
        else
          if(fr0 .eq. 3.d0 .or. fr0 .eq. 1.d0)then
            cmp%value(kytbl(kwFRMD,lt))=cmp%value(kytbl(kwFRMD,lt))+1.d0
          endif
        endif
      endif
      if(cmp%value(kytbl(kwFRMD,lt)) .eq. 0.d0)then
        cmp%value(kytbl(kwFRMD,lt))=-4.d0
      endif
 9000 if(kytbl(kwL,lt) .ne. 0)then
        cmp%value(kytbl(kwL,lt))=cmp%value(kytbl(kwL,lt))*r
      endif
      chg=.true.
      if(.not. ideal)then
        cmp%update=cmp%nparam .le. 0
      endif
      return
      end

      subroutine qfracseg(cmp,cmp0,fr1,fr2,chg,irtc)
      use ffs
      use sad_main
      use tfstk
      use kyparam
      use ffs_seg
      implicit none
      type (sad_comp) cmp,cmp0
      type (sad_dlist), pointer :: lprof,lsegp,ls1
      type (sad_rlist), pointer :: kl
      real*8 fr1,fr2,df1,df2,al1,al2,s,al,al0
      integer*4 i1,i2,k,i,n,irtc,lt,j1,j2,js
      logical*4 chg,chg1
      chg=.false.
      irtc=0
      i1=0
      i2=0
      df1=fr1
      df2=fr2
      lt=idtype(cmp%id)
      k=kytbl(kwPROF,lt)
      if(k .eq. 0 .or. .not. tflistq(cmp%dvalue(k),lprof))then
        go to 100
      endif
      al0=cmp%value(kytbl(kwL,lt))
      if(al0 .eq. 0.d0)then
        go to 100
      endif
      if(ideal)then
        chg=.true.
        cmp0%dvalue(p_PROF_MULT)%k=ktfoper+mtfnull
        call tsetupseg(cmp0,lprof,lsegp,irtc)
        if(irtc .ne. 0)then
          return
        endif
      else
        call descr_sad(cmp%dvalue(p_PROF_MULT),lsegp)
      endif
      call descr_sad(lsegp%dbody(1),ls1)
      call descr_sad(ls1%dbody(2),kl)
      n=kl%nl
      if(n .le. 1)then
        i1=1
        i2=1
        go to 100
      endif
      al=sum(kl%rbody(1:n))*al0
      al1=al*fr1
      al2=al*fr2
      s=0.d0
      if(cmp%ori)then
        j1=1
        j2=n
        js=1
      else
        j1=n
        j2=1
        js=-1
      endif
      do i=j1,j2,js
        s=s+kl%rbody(i)*al0
        if(i1 .eq. 0 .and. s .ge. al1)then
          i1=i
          df1=1.d0-(s-al1)/(kl%rbody(i)*al0)
        endif
        if(s .ge. al2)then
          i2=i
          df2=1.d0-(s-al2)/(kl%rbody(i)*al0)
          if(i1 .eq. 0)then
            i1=i2
            df1=df2
          endif
          go to 100
        endif
      enddo
      if(i1 .eq. 0)then
        i1=j2
        df1=1.d0
      endif
      i2=j2
      df2=1.d0
 100  if(i1 .eq. 0)then
        call qfraccomp(cmp,fr1,fr2,ideal,chg1)
      else
        call qfraccompseg(cmp,cmp0,i1,df1,i2,df2,lsegp,ideal,chg1)
      endif
      chg=chg .or. chg1
      return
      end

      subroutine qfraccompseg(cmp,cmp0,i1,rx1,i2,rx2,lsegp,ideal,chg)
      use kyparam
      use tfstk
      use sad_main
      use ffs_pointer, only:idelc,direlc,idtypec,idvalc
      use mackw
      implicit none
      type (sad_comp) :: cmp,cmp0
      type (sad_dlist) lsegp
      type (sad_dlist), pointer:: lsegp1,lk,lk0
      integer*4 i1,i2,k,nseg,i,nk,is,j
      real*8 rx1,rx2,f1,f2,r1,r2
      logical*4 ideal,chg
      nseg=abs(i2-i1)+1
      nk=lsegp%nl
      cmp0%dvalue(p_PROF_MULT)=kxadaloc(-1,nk,lsegp1)
      do k=1,nk
        call descr_sad(lsegp%dbody(k),lk0)
        lsegp1%dbody(k)=kxadaloc(0,2,lk)
        lk%dbody(1)=lk0%dbody(1)
        lk%dbody(2)=kxavaloc(0,nseg)
      enddo
      is=merge(1,-1,i2 .ge. i1)
      do i=i1,i2,is
        if(i1 .eq. i2)then
          f1=merge(1.d0,0.d0,rx1 .eq. 0.d0)
          f2=merge(1.d0,0.d0,rx2 .eq. 0.d0)
          r1=rx1
          r2=rx2
        elseif(i .eq. i1)then
          r1=rx2
          f1=merge(1.d0,0.d0,rx1 .eq. 0.d0)
          r2=1.d0
        elseif(i .eq. i2)then
          r1=0.d0
          r2=rx2
          f2=merge(1.d0,0.d0,rx2 .eq. 0.d0)
        else
          r1=0.d0
          r2=1.d0
        endif
        j=i+1-merge(i1,i2,is .gt. 0)
        call qputfracseg(lsegp1,j,r2-r1,lsegp,i)
c fringes are not taken into account yet...
      enddo
      chg=.true.
      if(.not. ideal)then
        cmp%update=.false.
        cmp%updateseg=.true.
      endif
      return
      end

      subroutine qputfracseg(larg,i1,r,larg0,i)
      use ffs
      use sad_main
      use tfstk
      implicit none
      type (sad_dlist) :: larg,larg0
      type (sad_dlist) ,pointer :: lk,lk0
      type (sad_rlist), pointer :: lkv,lkv0
      real*8 r
      integer*4 k,i,i1
      do k=1,larg%nl
        call descr_sad(larg%dbody(k),lk)
        call descr_sad(lk%dbody(2),lkv)
        call descr_sad(larg0%dbody(k),lk0)
        call descr_sad(lk0%dbody(2),lkv0)
        lkv%rbody(i1)=r*lkv0%rbody(i)
c        write(*,'(a,3i5,1p2g15.7)')'qputfracseg ',k,i1,i,r,lkv0%rbody(i)
      enddo
      return
      end

      subroutine qmultseg(trans,cod,l1,cmp,lsegp,bz,coup)
      use kyparam
      use tfstk
      use ffs
      use ffs_pointer, only:tsetfringepe
      use tffitcode
      use sad_main
      implicit none
      type (sad_comp) :: cmp
      type (sad_dlist) :: lsegp
      type (sad_dlist) , pointer :: lk
      type (sad_rlist), pointer :: lak,lal,lkv
      real*8 :: rsave(cmp%ncomp2)
      real*8 trans(4,5),cod(6),bz
      integer*8 kk
      integer*4 l1,i,nseg,i1,i2,istep,k,nk,k1,k2
      logical*4 coup,coup1
      integer*4 , parameter :: nc=ky_PROF_MULT-1
      coup=.false.
      rsave(1:nc)=cmp%value(1:nc)
      nk=lsegp%nl
      call descr_sad(lsegp%dbody(1),lal)
      call descr_sad(lal%dbody(2),lak)
      nseg=lak%nl
      if(cmp%ori)then
        i1=1
        i2=nseg
        istep=1
      else
        i1=nseg
        i2=1
        istep=-1
      endif
      do i=i1,i2,istep
        do k=1,nc
          if(integv(k,icMULT))then
            cmp%value(k)=rsave(k)*lak%rbody(i)
          else
            cmp%value(k)=rsave(k)
          endif
        enddo
        do k=1,nk
          call descr_sad(lsegp%dbody(k),lk)
          call descr_sad(lk%dbody(2),lkv)
          kk=ktfaddr(lsegp%dbody(k))
          k1=ilist(1,kk+1)
          k2=ilist(2,kk+1)
          if(k1 .eq. k2)then
            cmp%value(k1)=0.d0
          endif
          cmp%value(k1)=cmp%value(k1)+rsave(k2)*lkv%rbody(i)
        enddo
        call qmult1(trans,cod,l1,cmp,bz,i .eq. i1,coup1)
        coup=coup .or. coup1
      enddo
      cmp%value(1:nc)=rsave(1:nc)
      return
      end

      subroutine qmult1(trans,cod,l1,cmp,bz,ini,coup)
      use kyparam
      use tfstk
      use ffs
      use ffs_pointer, only:tsetfringepe
      use tffitcode
      use sad_main
      implicit none
      type (sad_comp) :: cmp
      real*8 trans(4,5),cod(6),ftable(4),psi1,psi2,apsi1,apsi2,
     $     fb1,fb2,chi1m,chi2m,al,bz,phi
      integer*4 mfr,l1
      logical*4 coup,ini
      al=cmp%value(ky_L_MULT)
      phi=cmp%value(ky_ANGL_MULT)
      mfr=nint(cmp%value(ky_FRMD_MULT))
      if(cmp%ori)then
        psi1=cmp%value(ky_E1_MULT)
        psi2=cmp%value(ky_E2_MULT)
        apsi1=cmp%value(ky_AE1_MULT)
        apsi2=cmp%value(ky_AE2_MULT)
        fb1=cmp%value(ky_FB1_MULT)
        fb2=cmp%value(ky_FB2_MULT)
        chi1m=cmp%value(ky_CHI1_MULT)
        chi2m=cmp%value(ky_CHI2_MULT)
      else
        mfr=mfr*(11+mfr*(2*mfr-9))/2
        psi1=cmp%value(ky_E2_MULT)
        psi2=cmp%value(ky_E1_MULT)
        apsi1=cmp%value(ky_AE2_MULT)
        apsi2=cmp%value(ky_AE1_MULT)
        fb2=cmp%value(ky_FB1_MULT)
        fb1=cmp%value(ky_FB2_MULT)
        chi1m=-cmp%value(ky_CHI1_MULT)
        chi2m=-cmp%value(ky_CHI2_MULT)
      endif
      call tsetfringepe(cmp,icMULT,ftable)
      call qmult(trans,cod,l1,al,
     $     cmp%value(ky_K0_MULT),bz,
     $     phi,psi1,psi2,apsi1,apsi2,
     1     cmp%value(ky_DX_MULT),cmp%value(ky_DY_MULT),
     $     cmp%value(ky_DZ_MULT),
     $     chi1m,chi2m,cmp%value(ky_ROT_MULT),
     $     cmp%value(ky_DROT_MULT),
     $     cmp%value(ky_EPS_MULT),
     $     cmp%value(ky_FRIN_MULT) .eq. 0.d0,
     $     ftable(1),ftable(2),ftable(3),ftable(4),
     $     mfr,fb1,fb2,
     $     cmp%value(ky_K0FR_MULT) .eq. 0.d0,
     $     cmp%value(ky_VOLT_MULT)+cmp%value(ky_DVOLT_MULT),
     $     cmp%value(ky_HARM_MULT),
     $     cmp%value(ky_PHI_MULT),cmp%value(ky_FREQ_MULT),
     $     cmp%value(ky_W1_MULT),
     $     cmp%value(ky_APHI_MULT) .ne. 0.d0,ini,
     $     coup)
      return
      end

      subroutine qthin(trans,cod,nord,al,ak,
     1                 dx,dy,theta,coup)
      implicit none
      integer*4 nord
      real*8 trans(4,5),cod(6),transe(6,12),beam(42),srot(3,9),
     $     dx,dy,theta,ak,al
      logical*4 coup
      call tinitr(transe)
      call tthine(transe,cod,beam,srot,nord,al,ak,
     1     dx,dy,theta,.false.,1)
      call qcopymat(trans,transe,.false.)
      coup=trans(1,3) .ne. 0.d0 .or. trans(1,4) .ne. 0.d0 .or.
     $     trans(2,3) .ne. 0.d0 .or. trans(2,4) .ne. 0.d0
      return
      end

      subroutine qquad(trans,cod,al,ak,
     1dx,dy,theta,fringe,f1in,f2in,f1out,f2out,mfring,eps0,
     $     kin,achro,coup)
      implicit none
      integer*4 mfring
      real*8 trans(4,5),cod(6),transe(6,12),beam(42),srot(3,9),
     $     dx,dy,theta,ak,eps0,al,f1in,f2in,f1out,f2out
      logical*4 fringe,coup,kin,achro
      call tinitr(transe)
      call tquade(transe,cod,beam,srot,al,ak,0.d0,
     1     dx,dy,theta,.false.,fringe,f1in,f2in,f1out,f2out,mfring,eps0,
     $     kin,achro,.false.)
      call qcopymat(trans,transe,.false.)
      coup=trans(1,3) .ne. 0.d0 .or. trans(1,4) .ne. 0.d0 .or.
     $     trans(2,3) .ne. 0.d0 .or. trans(2,4) .ne. 0.d0
      return
      end
