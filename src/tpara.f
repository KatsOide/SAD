      subroutine tpara()
      use tfstk
      use tmacro
      use sad_main
      use ffs_pointer, only:idelc,direlc,elatt,idvalc,compelc,
     $     tsetfringep
      use wsbb, only:nblist
      implicit none
      type (sad_comp), pointer :: cmp
      integer*4 lele,l
      integer*8 ib,id
      real*8 dir,phi,f1,al,psi1,psi2,theta,dtheta,w,akk,sk1
      if(tparaed)then
        return
      endif
      call tclrpara(elatt,nlat-1)
      lele=0
      do l=1,nlat-1
        call compelc(l,cmp)
        lele=idtype(cmp%id)
        dir=direlc(l)
        al=cmp%value(kytbl(kwL,lele))
c     go to (1010,1200,1010,1400,1010,1600,1010,1600,1010,1600,
c     1         1010,1600,1010,1010,1010,1010,1010,1800,1900,1010,
c     1         2100,2200,1010,1010,1010,1010,1010,1010,1010,1010,
c     1         1010,1010,1010,1010,1010,3600,3700,1010,1010,1010),lele
c     go to 5000
        select case (lele)

        case (icBEND)
          ib=ktaloc(14)
          id=idvalc(l)
          phi=cmp%value(kytbl(kwANGL,icBEND))
          cmp%param=ib
          f1=cmp%value(15)
          if(f1 .ne. 0.d0 .and. nint(cmp%value(16)) .ne. 0 .and.
     1         al .ne. 0.d0 .and. phi .ne. 0.d0)then
            al=al
     1           -(phi*f1)**2/cmp%value(1)/24.d0
     1           *sin(.5d0*phi*(1.d0-cmp%value(3)-cmp%value(4)))
     $           /sin(.5d0*phi)
          endif
          rlist(ib  )=al
          if(dir .gt. 0.d0)then
            psi1=cmp%value(kytbl(kwE1,icBEND))*phi
     $           +cmp%value(kytbl(kwAE1,icBEND))
            psi2=cmp%value(kytbl(kwE2,icBEND))*phi
     $           +cmp%value(kytbl(kwAE2,icBEND))
          else
            psi1=cmp%value(kytbl(kwE2,icBEND))*phi
     $           +cmp%value(kytbl(kwAE2,icBEND))
            psi2=cmp%value(kytbl(kwE1,icBEND))*phi
     $           +cmp%value(kytbl(kwAE1,icBEND))
          endif
          dtheta=cmp%value(kytbl(kwDROT,icBEND))
          theta=cmp%value(kytbl(kwROT,icBEND))+dtheta
          rlist(ib+1)=cos(psi1)
          rlist(ib+2)=sin(psi1)
          rlist(ib+3)=cos(psi2)
          rlist(ib+4)=sin(psi2)
          rlist(ib+5)=cos(theta)
          rlist(ib+6)=sin(theta)
          w=phi-psi1-psi2
          rlist(ib+7)=cos(w)
          rlist(ib+8)=sin(w)
          if(rlist(ib+7) .ge. 0.d0)then
            rlist(ib+9)=rlist(ib+8)**2/(1.d0+rlist(ib+7))
          else
            rlist(ib+9)=1.d0-rlist(ib+7)
          endif
          rlist(ib+10)=sin(phi-psi2)
          rlist(ib+11)=phi*sin(.5d0*dtheta)**2
          rlist(ib+12)=.5d0*phi*sin(dtheta)
          rlist(ib+13)=theta

        case (icQUAD)
          ib=ktaloc(10)
          cmp%param=ib
          if(al .ne. 0.d0)then
            akk=cmp%value(kytbl(kwK1,icQUAD))/al
            rlist(ib)  =sqrt(abs(akk))
            call tsetfringep(cmp,icQUAD,dir,akk,rlist(ib+6:ib+9))
          else
            rlist(ib)=1.d100
            rlist(ib+6:ib+9)=0.d0
          endif
          theta=cmp%value(4)
          rlist(ib+2)=cos(theta)
          rlist(ib+3)=sin(theta)
          rlist(ib+4)=theta
          rlist(ib+5)=1.d0

        case (icSEXT,icOCTU,icDECA,icDODECA)
          ib=ktaloc(5)
          cmp%param=ib
          theta=cmp%value(4)
          rlist(ib+2)=cos(theta)
          rlist(ib+3)=sin(theta)
          rlist(ib+4)=theta

        case (icUND)
          ib=ktaloc(20)
          cmp%param=ib
          call undinit(cmp%value(1),rlist(ib))

        case (icWIG)
          call twigp()

        case (icST)
          ib=ktaloc(5)
          cmp%param=ib
          theta=cmp%value(3)
          rlist(ib+2)=cos(theta)
          rlist(ib+3)=sin(theta)
          rlist(ib+4)=theta

        case (icMULT)
          ib=ktaloc(5)
          cmp%param=ib
          if(al .ne. 0.d0)then
            sk1=cmp%value(kytbl(kwSK1,icMULT))
            if(sk1 .eq. 0.d0)then
              akk=cmp%value(kytbl(kwK1,icMULT))/al
            else
              akk=sqrt(cmp%value(kytbl(kwK1,icMULT))**2+sk1**2)/al
            endif
            call tsetfringep(cmp,icMULT,dir,akk,rlist(ib+1:ib+4))
          else
            rlist(ib+1:ib+4)=0.d0
          endif

        case (icBEAM)
          ib=ktaloc(nblist)
c     This number should be 2*nslicemax+200 and equal to nblist.        
          cmp%param=ib
          call bbinit(cmp%value(1),rlist(ib))

        case (icPROT)
          ib=ktaloc(240)
          cmp%param=ib
          call phsinit(cmp%value(1),rlist(ib))

        case default
        end select
      enddo
      tparaed=.true.
      return
      end
