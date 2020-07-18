      subroutine tfbessel(isp1,kx,mode,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) kx
      integer*4 isp1,irtc,mode,itfmessage
      complex*16 cbesselj,cbessely,cbesseli,cbesselk
      external cbesselj,cbessely,cbesseli,cbesselk
      if(isp .ne. isp1+2)then
        irtc=itfmessage(9,'General::narg','"2"')
        return
      endif
      if(mode .eq. 0)then
        call tfbesself(dtastk(isp1+1),dtastk(isp),kx,cbesselj,irtc)
      elseif(mode .eq. 1)then
        call tfbesself(dtastk(isp1+1),dtastk(isp),kx,cbessely,irtc)
      elseif(mode .eq. 2)then
        call tfbesself(dtastk(isp1+1),dtastk(isp),kx,cbesseli,irtc)
      elseif(mode .eq. 3)then
        call tfbesself(dtastk(isp1+1),dtastk(isp),kx,cbesselk,irtc)
      endif
      return
      end

      recursive subroutine tfbesself(k1,k2,kx,cfun,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) k1,k2,kx
      type (sad_dlist), pointer :: kl1,kl2
      integer*4 irtc,itfmessage,n1,n2,isp0,isp2,i
      real*8 x1,x2
      complex*16 cfun,cx,c1,c2
      external cfun
      if(tfcomplexnumlistqk(k1%k,kl1))then
        n1=kl1%nl
        if(tfcomplexnumlistqk(k2%k,kl2))then
          if(n1 .ne. kl2%nl)then
            irtc=itfmessage(9,'General::equalleng','"#1 and #2"')
            return
          endif
          isp0=isp
          call tfgetllstkall(kl1)
          call tfgetllstkall(kl2)
          isp2=isp
          do i=1,n1
            isp=isp+1
            call tfbesself(dtastk(isp0+i),dtastk(isp0+n1+i),
     $           dtastk(isp),cfun,irtc)
            if(irtc .ne. 0)then
              isp=isp0
              return
            endif
          enddo
          kx=kxmakelist(isp2)
          isp=isp0
        elseif(tfnumberq(k2))then
          isp0=isp
          call tfgetllstkall(kl1)
          isp2=isp
          do i=1,n1
            isp=isp+1
            call tfbesself(dtastk(isp0+i),k2,dtastk(isp),cfun,irtc)
            if(irtc .ne. 0)then
              isp=isp0
              return
            endif
          enddo
          kx=kxmakelist(isp2)
          isp=isp0
        else
          irtc=-1
          return
        endif
      elseif(tfcomplexnumlistqk(k2%k,kl2))then
        n2=kl2%nl
        isp0=isp
        call tfgetllstkall(kl2)
        isp2=isp
        do i=1,n2
          isp=isp+1
          call tfbesself(k1,dtastk(isp0+i),dtastk(isp),cfun,irtc)
          if(irtc .ne. 0)then
            isp=isp0
            return
          endif
        enddo
        kx=kxmakelist(isp2)
        isp=isp0
      elseif(ktfrealq(k1,x1))then
        if(ktfrealq(k2,x2))then
          kx=dfromr(dble(cfun(dcmplx(x1,0.d0),
     $         dcmplx(x2,0.d0))))
        elseif(tfcomplexq(k2,c2))then
          cx=cfun(dcmplx(x1,0.d0),c2)
          go to 10
        else
          irtc=-1
          return
        endif
      elseif(tfcomplexq(k1,c1))then
        if(ktfrealq(k2,x2))then
          cx=cfun(c1,dcmplx(x2,0.d0))
        elseif(tfcomplexq(k2,c2))then
          cx=cfun(c1,c2)
        else
          irtc=-1
          return
        endif
        go to 10
      else
        irtc=-1
        return
      endif
      irtc=0
      return
 10   if(imag(cx) .ne. 0.d0)then
        kx=kxcalocv(-1,dble(cx),imag(cx))
      else
        kx=dfromr(dble(cx))
      endif
      irtc=0
      return
      end

      complex*16 function cbesseli(cn,z)
      implicit none
      complex*16 cn,z,cbesselj
      if(imag(z) .eq. 0.d0 .and.
     $     imag(cn) .eq. 0.d0)then
        cbesseli=dble((0.d0,-1.d0)**cn*
     $       cbesselj(cn,dcmplx(0.d0,dble(z))))
      else
        cbesseli=(0.d0,-1.d0)**cn*
     $       cbesselj(cn,dcmplx(-imag(z),dble(z)))
      endif
      return
      end

      complex*16 function cbesselk(cn,z)
      use macmath
      implicit none
c     Including euler(Euler-Mascheroni constant)/pi symbol
      complex*16 cn,z,cbessely,cbesselj
      if(imag(z) .eq. 0.d0 .and.
     $     imag(cn) .eq. 0.d0)then
        cbesselk=m_pi_2*dble((0.d0,1.d0)**(cn+1.d0)*(
     $       cbesselj(cn,dcmplx(0.d0,dble(z)))+
     $       (0.d0,1.d0)*cbessely(cn,dcmplx(0.d0,dble(z)))))
      else
        cbesselk=m_pi_2*(0.d0,1.d0)**(cn+1.d0)*(
     $       cbesselj(cn,dcmplx(-imag(z),dble(z)))+
     $       (0.d0,1.d0)*cbessely(cn,dcmplx(-imag(z),dble(z))))
      endif
      return
      end

      recursive complex*16 function cbessely(cn,z)
     $     result(cb)
      use macmath
      use mathfun
      use gammaf,only:cgamma
      implicit none
c     Including m_pi_2 symbol
      integer*4 i,itmax,m
      parameter (itmax=26)
      complex*16 cn,z,cbesselj,cj,cg1,cg2,cgamm1,cgamm2,ca3,
     $     cnf,cp,cq,cr,cs0,cs1,zp,cx,cxp,csigma,
     $     ca1,ca2,cf,cg,cf1,cf2,clogz,zhi,cnfpi,cs
      real*8 az,xcn,an0,anf,acnf,ak
      az=abs(z)
      if(az .gt. 12.2d0)then
        call cbesasym(cn,z,cj,cb)
      else
        xcn=dble(cn)
        if(xcn .lt. 0.d0)then
          cb=-sin(pi*cn)*cbesselj(-cn,z)+
     $         cos(pi*cn)*cbessely(-cn,z)
        else
          an0=aint(xcn)
          anf=xcn-an0
          if(anf .gt. .5d0)then
            anf=anf-1.d0
            an0=an0+1.d0
          endif
          cnf=dcmplx(anf,imag(cn))
          acnf=abs(cnf)
          zhi=2.d0/z
          clogz=log(zhi)
          if(acnf .eq. 0.d0)then
            zp=1.d0
            cp=1.d0/pi
            cq=1.d0/pi
            cg1=-euler
            cg2=1.d0
            ca1=1.d0
            ca2=1.d0
            ca3=1.d0
            cr=0.d0
          else
            zp=exp(cnf*clogz)
            cgamm1=cgamma(1.d0+cnf)
            cgamm2=cgamma(1.d0-cnf)
            cp=cgamm1*zp/pi
            cq=cgamm2/zp/pi
            cnfpi=cnf*pi
            cr=2.d0/cnf*sin(cnfpi*.5d0)**2
            if(acnf .lt. 1.d-3)then
              cg1=(-euler+cnf**2*(0.042002635034095236d0+
     $             0.042197734555544337d0*cnf**2))
            else
              cg1=(1.d0/cgamm2-1.d0/cgamm1)*.5d0/cnf
            endif
            ca1=cnfpi/sin(cnfpi)
            csigma=cnf*clogz
            ca2=tcsinh(csigma)/csigma
            ca3=tccosh(csigma)
            cg2=.5d0*(1.d0/cgamm2+1.d0/cgamm1)
          endif
          cf=2.d0/pi*ca1*(ca3*cg1+ca2*clogz*cg2)
          cs0=cf+cr*cq
          cs1=cp
          ak=0.d0
          cx=-(z*.5d0)**2
          cxp=1.d0
          m=min(itmax,max(int(az*2.3d0),10))
          do i=1,m
            ak=ak+1.d0
            cxp=cxp*cx/ak
            cf1=1.d0/(ak-cnf)
            cf2=1.d0/(ak+cnf)
            cf=(ak*cf+cp+cq)*cf1*cf2
            cp=cp*cf1
            cq=cq*cf2
            cg=cf+cr*cq
            cs0=cs0+cxp*cg
            cs1=cs1+cxp*(cp-ak*cg)
          enddo
          if(an0 .eq. 0.d0)then
            cb=-cs0
          elseif(an0 .eq. 1.d0)then
            cb=-cs1*zhi
          else
            cs0=-cs0
            cs1=-cs1*zhi
            do i=1,int(an0)-1
              cnf=cnf+1.d0
              cs=cnf*zhi*cs1-cs0
              cs0=cs1
              cs1=cs
            enddo
            cb=cs1
          endif
        endif
      endif
      return
      end

      recursive complex*16 function cbesselj(cn,z)
     $     result(cb)
      use gammaf
      implicit none
      integer*4 nmax,i
      complex*16 cn,z,cbessely,cj,cj1,cj2,cs,zi,
     $     cnk,cn2k,cg
      real*8 az,acn,an,ak,xcn
c     Initialization to avoid compiler warning
      nmax=0
c
      az=abs(z)
      acn=abs(cn)
      if(az .gt. acn*.25d0+20.d0)then
        call cbesasym(cn,z,cb,cbessely)
      else
        if(az .ne. 0.d0)then
          nmax=int(max(1.d0,20.d0/max(log10(2.d0/az),1.d0),
     $         2.d0*az-dble(cn)))
        endif
        if(cn .eq. (0.d0,0.d0))then
          if(az .eq. 0.d0)then
            cb=1.d0
          else
            cs=1.d-100
            cj2=0.d0
            cj1=1.d-100
            an=nmax*2
            zi=2.d0/z
            do i=2,nmax
              cj2=an*cj1*zi-cj2
              cj1=(an-1.d0)*zi*cj2-cj1
              an=an-2.d0
              cs=cs+cj1
            enddo
            cj=zi*(2.d0*zi*cj1-cj2)-cj1
            cb=cj/(cs*2.d0+cj)
          endif
        else
          xcn=dble(cn)
          if(xcn .lt. 0.d0 .and. cn .eq. aint(xcn))then
            if(mod(-xcn,2.d0) .eq. 0.d0)then
              cb=cbesselj(-cn,z)
            else
              cb=-cbesselj(-cn,z)
            endif
          elseif(az .eq. 0.d0)then
            cb=0.d0
          else
            cj2=0.d0
            cj1=1.d-100
            ak=dble(nmax)
            cnk=cn+ak
            cn2k=cnk+ak
            cg=cgamma(cnk)/factorial(ak)
            cs=cn2k*cg*cj1
            zi=2.d0/z
            do i=1,nmax
              cj2=cn2k*cj1*zi-cj2
              cj1=(cn2k-1.d0)*zi*cj2-cj1
              cn2k=cn2k-2.d0
              cnk=cnk-1.d0
              cg=ak*cg/cnk
              ak=ak-1.d0
              cs=cs+cj1*cn2k*cg
            enddo
            cb=cj1*(.5d0*z)**cn/cs
          endif
        endif
      endif
      return
      end

      subroutine cbesasym(cn,z,cbesj,cbesy)
      use macmath
      implicit none
c     Including m_pi_2/m_pi_4/m_2_pi symbol
      complex*16 cn,z,cp,cq,cmu,cbesj,cbesy,z1,chi,csqrtz,
     $     ccoschi,csinchi
      cmu=4.d0*cn**2
      z1=1.d0/64.d0/z**2
      cp=1.d0-(cmu-1.d0)*(cmu-9.d0)*z1/2.d0*
     $     (1.d0-(cmu-25.d0)*(cmu-49.d0)*z1/12.d0*
     $     (1.d0-(cmu-81.d0)*(cmu-121.d0)*z1/30.d0*
     $     (1.d0-(cmu-169.d0)*(cmu-225.d0)*z1/56.d0*
     $     (1.d0-(cmu-289.d0)*(cmu-361.d0)*z1/90.d0*
     $     (1.d0-(cmu-441.d0)*(cmu-529.d0)*z1/132.d0*
     $     (1.d0-(cmu-25.d0**2)*(cmu-27.d0**2)*z1/(13.d0*14.d0)*
     $     (1.d0-(cmu-29.d0**2)*(cmu-31.d0**2)*z1/(15.d0*16.d0)*
     $     (1.d0-(cmu-33.d0**2)*(cmu-35.d0**2)*z1/(17.d0*18.d0)*
     $     (1.d0-(cmu-37.d0**2)*(cmu-39.d0**2)*z1/(19.d0*20.d0)*
     $     (1.d0-(cmu-41.d0**2)*(cmu-43.d0**2)*z1/(21.d0*22.d0)*
     $     (1.d0-(cmu-45.d0**2)*(cmu-47.d0**2)*z1/(23.d0*24.d0)
     $     )))))))))))
      cq=(cmu-1.d0)/8.d0/z*
     $     (1.d0-(cmu-9.d0)*(cmu-25.d0)*z1/6.d0*
     $     (1.d0-(cmu-49.d0)*(cmu-81.d0)*z1/20.d0*
     $     (1.d0-(cmu-121.d0)*(cmu-169.d0)*z1/42.d0*
     $     (1.d0-(cmu-225.d0)*(cmu-289.d0)*z1/72.d0*
     $     (1.d0-(cmu-361.d0)*(cmu-441.d0)*z1/110.d0*
     $     (1.d0-(cmu-23.d0**2)*(cmu-25.d0**2)*z1/(12.d0*13.d0)*
     $     (1.d0-(cmu-27.d0**2)*(cmu-29.d0**2)*z1/(14.d0*15.d0)*
     $     (1.d0-(cmu-31.d0**2)*(cmu-33.d0**2)*z1/(16.d0*17.d0)*
     $     (1.d0-(cmu-35.d0**2)*(cmu-37.d0**2)*z1/(18.d0*19.d0)*
     $     (1.d0-(cmu-39.d0**2)*(cmu-41.d0**2)*z1/(20.d0*21.d0)*
     $     (1.d0-(cmu-43.d0**2)*(cmu-45.d0**2)*z1/(22.d0*23.d0)
     $     )))))))))))
      chi=z-cn*m_pi_2-m_pi_4
      ccoschi=cos(chi)
      csinchi=sin(chi)
      csqrtz=sqrt(m_2_pi/z)
      cbesj=csqrtz*(cp*ccoschi-cq*csinchi)
      cbesy=csqrtz*(cp*csinchi+cq*ccoschi)
      return
      end
