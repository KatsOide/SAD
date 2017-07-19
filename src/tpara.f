      module tparastat
        use sad_main
        use tmacro
        real*8, save :: rstat(2,icMARK)=0.d0
        logical*4, save :: lstat(2,icMARK)=.true.

        contains
        logical*4 function tparacheck(l,cmp)
        use ffs_flag
        implicit none
        type (sad_comp) :: cmp
        integer*4 l

        select case (l)
        case (icMULT,icCAVI)
          tparacheck=rstat(1,l) .ne. amass .or.
     $         rstat(2,l) .ne. charge .or.
     $         lstat(1,l) .neqv. trpt
          if(tparacheck)then
            rstat(1,l)=amass
            rstat(2,l)=charge
            lstat(1,l)=trpt
          else
            tparacheck=cmp%update .eq. 0
          endif

        case default
          tparacheck=cmp%update .eq. 0

        end select
        return
        end function

        subroutine tpara(cmp)
        use kyparam
        use tfstk
        use ffs_pointer, only:idelc,idvalc,compelc,tsetfringep
        use ffs_flag, only:trpt
        implicit none
        type (sad_comp) :: cmp
        integer*4 ltyp
        real*8 phi,al,psi1,psi2,theta,dtheta,w,akk,sk1,
     $       fb1,fb2,harm,vnominal,frmd

        cmp%update=1
        ltyp=idtype(cmp%id)
        if(kytbl(kwNPARAM,ltyp) .eq. 0)then
          return
        endif
        select case (ltyp)
        case (icBEND)
          phi=cmp%value(ky_ANGL_BEND)
          if(cmp%orient .gt. 0.d0)then
            psi1=cmp%value(ky_E1_BEND)*phi+cmp%value(ky_AE1_BEND)
            psi2=cmp%value(ky_E2_BEND)*phi+cmp%value(ky_AE2_BEND)
            fb1=cmp%value(ky_F1_BEND)+cmp%value(ky_FB1_BEND)
            fb2=cmp%value(ky_F1_BEND)+cmp%value(ky_FB2_BEND)
          else
            psi1=cmp%value(ky_E2_BEND)*phi+cmp%value(ky_AE2_BEND)
            psi2=cmp%value(ky_E1_BEND)*phi+cmp%value(ky_AE1_BEND)
            fb2=cmp%value(ky_F1_BEND)+cmp%value(ky_FB1_BEND)
            fb1=cmp%value(ky_F1_BEND)+cmp%value(ky_FB2_BEND)
          endif
          if(cmp%value(ky_FRMD_BEND) .eq. 0.d0)then
            fb1=0.d0
            fb2=0.d0
          endif
          al=cmp%value(ky_L_BEND)
          w=phi-psi1-psi2
          if((fb1 .ne. 0.d0 .or. fb2 .ne. 0.d0) .and.
     1         al .ne. 0.d0 .and. phi .ne. 0.d0)then
            al=al-((phi*fb1)**2+(phi*fb2)**2)/al/48.d0
     1           *sin(.5d0*w)/sin(.5d0*phi)
          endif
          cmp%value(p_L_BEND)=al
          dtheta=cmp%value(ky_DROT_BEND)
          theta=cmp%value(ky_ROT_BEND)+dtheta
          cmp%value(p_COSPSI1_BEND)=cos(psi1)
          cmp%value(p_SINPSI1_BEND)=sin(psi1)
          cmp%value(p_COSPSI2_BEND)=cos(psi2)
          cmp%value(p_SINPSI2_BEND)=sin(psi2)
          cmp%value(p_COSTHETA_BEND)=cos(theta)
          cmp%value(p_SINTHETA_BEND)=sin(theta)
          cmp%value(p_COSW_BEND)=cos(w)
          cmp%value(p_SINW_BEND)=sin(w)
          if(cmp%value(p_COSW_BEND) .ge. 0.d0)then
            cmp%value(p_SQWH_BEND)=cmp%value(p_SINW_BEND)**2
     $           /(1.d0+cmp%value(p_COSW_BEND))
          else
            cmp%value(p_SQWH_BEND)=1.d0-cmp%value(p_COSW_BEND)
          endif
          cmp%value(p_SINWP1_BEND)=sin(phi-psi2)
          cmp%value(p_DPHIX_BEND)=phi*sin(.5d0*dtheta)**2
          cmp%value(p_DPHIY_BEND)=.5d0*phi*sin(dtheta)
          cmp%value(p_THETA_BEND)=theta
          cmp%value(p_FB1_BEND)=fb1
          cmp%value(p_FB2_BEND)=fb2

        case (icQUAD)
          al=cmp%value(ky_L_QUAD)
          if(al .ne. 0.d0)then
            akk=cmp%value(ky_K1_QUAD)/al
            cmp%value(p_SQRTK_QUAD)  =sqrt(abs(akk))
            call tsetfringep(cmp,icQUAD,cmp%orient,akk,
     $           cmp%value(p_AKF1F_QUAD:p_AKF2B_QUAD))
          else
            cmp%value(p_SQRTK_QUAD)=1.d100
            cmp%value(p_AKF1F_QUAD:p_AKF2B_QUAD)=0.d0
          endif
          theta=cmp%value(ky_ROT_QUAD)
          cmp%value(p_COSTHETA_QUAD)=cos(theta)
          cmp%value(p_SINTHETA_QUAD)=sin(theta)
          cmp%value(p_THETA_QUAD)=theta
          frmd=cmp%value(ky_FRMD_QUAD)
          if(cmp%orient .lt. 0.d0)then
            frmd=frmd*(11.d0+frmd*(2.d0*frmd-9.d0))/2.d0
          endif
          cmp%value(p_FRMD_QUAD)=frmd

        case (icSEXT,icOCTU,icDECA,icDODECA)
          theta=cmp%value(ky_ROT_THIN)
          cmp%value(p_COSTHETA_THIN)=cos(theta)
          cmp%value(p_SINTHETA_THIN)=sin(theta)
          cmp%value(p_THETA_THIN)=theta

        case (icUND)
          call undinit(cmp%value(1),cmp%value(p_PARAM_UND))

        case (icWIG)
          call twigp()

        case (icMULT)
          al=cmp%value(ky_L_MULT)
          phi=cmp%value(ky_ANGL_MULT)
          if(cmp%orient .gt. 0.d0)then
            psi1=cmp%value(ky_E1_MULT)*phi+cmp%value(ky_AE1_MULT)
            psi2=cmp%value(ky_E2_MULT)*phi+cmp%value(ky_AE2_MULT)
            fb1=cmp%value(ky_FB1_MULT)
            fb2=cmp%value(ky_FB2_MULT)
          else
            psi1=cmp%value(ky_E2_MULT)*phi+cmp%value(ky_AE2_MULT)
            psi2=cmp%value(ky_E1_MULT)*phi+cmp%value(ky_AE1_MULT)
            fb2=cmp%value(ky_FB1_MULT)
            fb1=cmp%value(ky_FB2_MULT)
          endif
          frmd=cmp%value(ky_FRMD_MULT)
          if(cmp%orient .lt. 0.d0)then
            frmd=frmd*(11.d0+frmd*(2.d0*frmd-9.d0))/2.d0
          endif
          cmp%value(p_FRMD_MULT)=frmd
          if(frmd .ne. 3.d0 .and. frmd .ne. 1.d0)then
            fb1=0.d0
          endif
          if(frmd .ne. 3.d0 .and. frmd .ne. 2.d0)then
            fb2=0.d0
          endif
          w=phi-psi1-psi2
          if((fb1 .ne. 0.d0 .or. fb2 .ne. 0.d0) .and.
     1         al .ne. 0.d0 .and. phi .ne. 0.d0)then
            al=al-((phi*fb1)**2+(phi*fb2)**2)/al/48.d0
     $           *sin(.5d0*w)/sin(.5d0*phi)
          endif
          cmp%value(p_L_MULT)=al
          cmp%value(p_PSI1_MULT)=psi1
          cmp%value(p_PSI2_MULT)=psi2
          cmp%value(p_FB1_MULT)=fb1
          cmp%value(p_FB2_MULT)=fb2
          if(al .ne. 0.d0)then
            sk1=cmp%value(ky_SK1_MULT)
            if(sk1 .eq. 0.d0)then
              akk=cmp%value(ky_K1_MULT)/al
            else
              akk=sqrt(cmp%value(ky_K1_MULT)**2+sk1**2)/al
            endif
            call tsetfringep(cmp,icMULT,cmp%orient,akk,
     $           cmp%value(p_AKF1F_MULT:p_AKF2B_MULT))
          else
            cmp%value(p_AKF1F_MULT:p_AKF2B_MULT)=0.d0
          endif
          harm=cmp%value(ky_HARM_MULT)
          if(harm .eq. 0.d0)then
            w=pi2*cmp%value(ky_FREQ_MULT)/c
          else
            w=omega0*harm/c
          endif
          if(trpt)then
            vnominal=cmp%value(ky_VOLT_MULT)/amass*abs(charge)
     $           *sin(-cmp%value(ky_PHI_MULT)*sign(1.d0,charge))
          else
            vnominal=0.d0
          endif
          cmp%value(p_W_MULT)=w
          cmp%value(p_VNOMINAL_MULT)=vnominal

        case (icCAVI)
          frmd=cmp%value(ky_FRMD_CAVI)
          if(cmp%orient .lt. 0.d0)then
            frmd=frmd*(11.d0+frmd*(2.d0*frmd-9.d0))/2.d0
          endif
          cmp%value(p_FRMD_CAVI)=frmd
          harm=cmp%value(ky_HARM_CAVI)
          if(harm .eq. 0.d0)then
            w=pi2*cmp%value(ky_FREQ_CAVI)/c
          else
            w=omega0*harm/c
          endif
          if(trpt)then
            vnominal=cmp%value(ky_VOLT_CAVI)/amass*abs(charge)
     $           *sin(-cmp%value(ky_PHI_CAVI)*sign(1.d0,charge))
          else
            vnominal=0.d0
          endif
          cmp%value(p_W_CAVI)=w
          cmp%value(p_VNOMINAL_CAVI)=vnominal

        case (icBEAM)
          call bbinit(cmp%value(1),cmp%value(p_PARAM_BEAM))

        case (icPROT)
          call phsinit(cmp%value(1),cmp%value(p_PARAM_Prot))

        case default
        end select
        return
        end subroutine
      end module
