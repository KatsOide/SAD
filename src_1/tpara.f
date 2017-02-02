      subroutine tpara(latt)
      use tfstk
      use tmacro
      implicit none
      integer*4 lele,l,lele1
      integer*8 latt(nlat),lp,lp1,ib,ktaloc,id
      real*8 dir,phi,f1,al,psi1,psi2,theta,dtheta,w,akk,sk1
      if(tparaed)then
        return
      endif
      call tclrpara(latt,nlat-1)
      lele=0
      lp=0
      do 1010 l=1,nlat-1
        lele1=lele
        lele=idtype(ilist(2,latt(l)))
        lp1=lp
        lp=latt(l)
        dir=rlist(lp+ilist(1,lp))
        go to (1010,1200,1010,1400,1010,1600,1010,1600,1010,1600,
     1         1010,1600,1010,1010,1010,1010,1010,1800,1900,1010,
     1         2100,2200,1010,1010,1010,1010,1010,1010,1010,1010,
     1         1010,1010,1010,1010,1010,3600,3700,1010,1010,1010),lele
        go to 5000
1200    ib=ktaloc(14)
        id=idval(ilist(2,latt(l)))
        phi=rlist(lp+kytbl(kwANGL,icBEND))
        klist(lp+ilist(1,lp)+1)=ib
        f1=rlist(lp+15)
        al=rlist(lp+kytbl(kwL,icBEND))
        if(f1 .ne. 0.d0 .and. nint(rlist(lp+16)) .ne. 0 .and.
     1     al .ne. 0.d0 .and. phi .ne. 0.d0)then
          al=al
     1      -(phi*f1)**2/rlist(lp+1)/24.d0
     1       *sin(.5d0*phi*(1.d0-rlist(lp+3)-rlist(lp+4)))/sin(.5d0*phi)
        endif
        rlist(ib  )=al
        if(dir .gt. 0.d0)then
          psi1=rlist(lp+kytbl(kwE1,icBEND))*phi
     $         +rlist(lp+kytbl(kwAE1,icBEND))
          psi2=rlist(lp+kytbl(kwE2,icBEND))*phi
     $         +rlist(lp+kytbl(kwAE2,icBEND))
        else
          psi1=rlist(lp+kytbl(kwE2,icBEND))*phi
     $         +rlist(lp+kytbl(kwAE2,icBEND))
          psi2=rlist(lp+kytbl(kwE1,icBEND))*phi
     $         +rlist(lp+kytbl(kwAE1,icBEND))
        endif
        dtheta=rlist(lp+kytbl(kwDROT,icBEND))
        theta=rlist(lp+kytbl(kwROT,icBEND))+dtheta
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
        go to 1010
 1400   ib=ktaloc(8)
        klist(lp+ilist(1,lp)+1)=ib
        if(rlist(lp+1) .ne. 0.d0)then
          akk=rlist(lp+2)/rlist(lp+1)
          rlist(ib)  =sqrt(abs(akk))
          rlist(ib+5)=-akk*rlist(lp+10)*abs(rlist(lp+10))/24.d0
          rlist(ib+6)=akk*rlist(lp+11)
        else
          rlist(ib)=1.d35
          rlist(ib+5)=0.d0
          rlist(ib+6)=0.d0
        endif
        theta=rlist(lp+4)
        rlist(ib+2)=cos(theta)
        rlist(ib+3)=sin(theta)
        rlist(ib+7)=theta
        rlist(ib+4)=1.d0
        go to 1010
 1600   ib=ktaloc(5)
        klist(lp+ilist(1,lp)+1)=ib
        theta=rlist(lp+4)
        rlist(ib+2)=cos(theta)
        rlist(ib+3)=sin(theta)
        rlist(ib+4)=theta
        go to 1010
 1800   ib=ktaloc(20)
        klist(lp+ilist(1,lp)+1)=ib
        call undinit(rlist(lp+1),rlist(ib))
        goto 1010
 1900   call twigp(lp)
        go to 1010
 2100   ib=ktaloc(5)
        klist(lp+ilist(1,lp)+1)=ib
        theta=rlist(lp+3)
        rlist(ib+2)=cos(theta)
        rlist(ib+3)=sin(theta)
        rlist(ib+4)=theta
        go to 1010
 2200   ib=ktaloc(4)
        klist(lp+ilist(1,lp)+1)=ib
        al=rlist(lp+kytbl(kwL,icMULT))
        if(al .ne. 0.d0)then
          if(rlist(lp+kytbl(kwANGL,icMULT)) .eq. 0.d0)then
            sk1=rlist(lp+kytbl(kwSK1,icMULT))
            if(sk1 .eq. 0.d0)then
              akk=rlist(lp+kytbl(kwK1,icMULT))/al
            else
              akk=sqrt(rlist(lp+kytbl(kwK1,icMULT))**2+sk1**2)/al
            endif
            rlist(ib+1)=-akk*rlist(lp+12)*abs(rlist(lp+12))/24.d0
            rlist(ib+2)=akk*rlist(lp+13)
          else
            rlist(ib+1)=0.d0
            rlist(ib+2)=0.d0
          endif
        else
          rlist(ib+1)=0.d0
          rlist(ib+2)=0.d0
        endif
        go to 1010
 3600   ib=ktaloc(1200)
c      This number should be 2*nslicemax+200 and equal to nblist.        
        klist(lp+ilist(1,lp)+1)=ib
        call bbinit(rlist(lp+1),rlist(ib))
        go to 1010
 3700   ib=ktaloc(240)
        klist(lp+ilist(1,lp)+1)=ib
        call phsinit(rlist(lp+1),rlist(ib))
        go to 1010
5000    continue
1010  continue
      tparaed=.true.
      return
      end

      subroutine tclrpara(latt,nl)
      use tfstk
      use tmacro
      implicit none
      integer*4 nl,i,iwpl,iwpt
      integer*8 latt(nl),lp
      do i=1,nl
        lp=latt(i)
        if(lp .gt. 0)then
          if(klist(lp+ilist(1,lp)+1) .gt. 0)then
            call tfree(klist(lp+ilist(1,lp)+1))
            klist(lp+ilist(1,lp)+1)=0
          endif
        endif
        if(idtype(ilist(2,lp)) .eq. 31)then
          iwpl=ilist(1,lp+kytbl(kwLWAK,icCAVI))
          if(iwpl .gt. 0)then
            if(ilist(2,iwpl) .gt. 0)then
              call tfree(int8(ilist(2,iwpl)))
              ilist(2,iwpl)=0
            endif
          endif
          iwpt=ilist(1,lp+kytbl(kwTWAK,icCAVI))
          if(iwpt .gt. 0)then
            if(ilist(2,iwpt) .gt. 0)then
              call tfree(int8(ilist(2,iwpt)))
              ilist(2,iwpt)=0
            endif
          endif
        endif
      enddo
      tparaed=.false.
      return
      end
