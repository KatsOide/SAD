      subroutine tpara(latt)
      use tfstk
      implicit none
      include 'inc/TMACRO1.inc'
      integer*4 lele,lp,l,lele1,lp1,ib,id,italoc
      integer*4 latt(2,nlat)
      real*8 dir,phi,f1,al,psi1,psi2,theta,dtheta,w,akk,sk1
      data tparaed/.false./
      if(tparaed)then
        return
      endif
      call tclrpara(latt,nlat-1)
      lele=0
      lp=0
      do 1010 l=1,nlat-1
        lele1=lele
        lele=idtype(latt(1,l))
        lp1=lp
        lp=latt(2,l)
        dir=rlist(lp+ilist(1,lp))
        ilist(2,lp)=0
        go to (1010,1200,1010,1400,1010,1600,1010,1600,1010,1600,
     1         1010,1600,1010,1010,1010,1010,1010,1800,1900,1010,
     1         2100,2200,1010,1010,1010,1010,1010,1010,1010,1010,
     1         1010,1010,1010,1010,1010,3600,3700,1010,1010,1010),lele
        go to 5000
1200    ib=italoc(14)
        id=idval(latt(1,l))
        phi=rlist(lp+kytbl(kwANGL,icBEND))
        ilist(2,lp)=ib
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
 1400   ib=italoc(8)
        ilist(2,lp)=ib
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
 1600   ib=italoc(5)
        ilist(2,lp)=ib
        theta=rlist(lp+4)
        rlist(ib+2)=cos(theta)
        rlist(ib+3)=sin(theta)
        rlist(ib+4)=theta
        go to 1010
 1800   ib=italoc(20)
        ilist(2,lp)=ib
        call undinit(rlist(lp+1),rlist(ib))
        goto 1010
1900      call twigp(lp)
        go to 1010
 2100   ib=italoc(5)
        ilist(2,lp)=ib
        theta=rlist(lp+3)
        rlist(ib+2)=cos(theta)
        rlist(ib+3)=sin(theta)
        rlist(ib+4)=theta
        go to 1010
 2200   ib=italoc(4)
        ilist(2,lp)=ib
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
 3600   ib=italoc(1200)
c      This number should be 2*nslicemax+200 and equal to nblist.        
        ilist(2,lp)=ib
        call bbinit(rlist(lp+1),rlist(ib))
        go to 1010
 3700   ib=italoc(240)
        ilist(2,lp)=ib
        call phsinit(rlist(lp+1),rlist(ib))
        go to 1010
5000    continue
1010  continue
      tparaed=.true.
      return
      end

      subroutine tclrpara(latt,nl)
      use tfstk
      implicit none
      include 'inc/TMACRO1.inc'
      integer*4 nl,latt(2,nl),i,lp,iwpl,iwpt
      do 100 i=1,nl
        lp=latt(2,i)
        if(lp .gt. 0)then
          if(ilist(2,lp) .gt. 0)then
            call tfree(int8(ilist(2,lp)))
            ilist(2,lp)=0
          endif
        endif
        if(idtype(latt(1,i)) .eq. 31)then
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
 100  continue
      tparaed=.false.
      return
      end
