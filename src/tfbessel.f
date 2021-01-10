      module bes

      contains
      complex*16 pure function cbesseli(cn,z)
      implicit none
      complex*16 ,intent(in):: cn,z
      cbesseli=(0.d0,1.d0)**cn*cbesselj(cn,dcmplx(imag(z),-dble(z)))
      return
      end

      complex*16 pure function cbesselk(cn,z)
      use macmath
      implicit none
c     Including euler(Euler-Mascheroni constant)/pi symbol
      complex*16 ,intent(in):: cn,z
      cbesselk=-m_pi_2*(0.d0,1.d0)**(1.d0-cn)*(
     $     cbesselj(cn,dcmplx(imag(z),-dble(z)))-
     $     (0.d0,1.d0)*cbessely(cn,dcmplx(imag(z),-dble(z))))
      return
      end

      recursive complex*16 pure function cbessely(cn,z)
     $     result(cb)
      use macmath
      use mathfun
      use gammaf,only:cgamma
      implicit none
c     Including m_pi_2 symbol
      integer*4 i,itmax,m
      parameter (itmax=26)
      complex*16 ,intent(in):: cn,z
      complex*16 cg1,cg2,cgamm1,cgamm2,ca3,
     $     cnf,cp,cq,cr,cs0,cs1,zp,cx,cxp,csigma,
     $     ca1,ca2,cf,cg,cf1,cf2,clogz,zhi,cnfpi,cs,cb2(2)
      real*8 az,xcn,an0,anf,acnf,ak,by,c,s,bj
      if(imag(cn) .eq. 0.d0 .and. imag(z) .eq. 0.d0)then
        if(dble(z) .ge. 0.d0)then
          cb=dcmplx(rrbessely(dble(cn),dble(z)),0.d0)
        else
          by=rrbessely(dble(cn),-dble(z))
          bj=rrbesselj(dble(cn),-dble(z))
          c=cos(m_pi*dble(cn))
          s=sin(m_pi*dble(cn))
          cb=dcmplx(c*by,2.d0*c*bj-s*by)
        endif
        return
      endif
      az=abs(z)
      if(az .gt. 12.2d0)then
        cb2=cbesasym(cn,z)
        cb=cb2(2)
      else
        xcn=dble(cn)
        if(xcn .lt. 0.d0)then
          cb=-sin(m_pi*cn)*cbesselj(-cn,z)+
     $         cos(m_pi*cn)*cbessely(-cn,z)
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
            cp=1.d0/m_pi
            cq=1.d0/m_pi
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
            cp=cgamm1*zp/m_pi
            cq=cgamm2/zp/m_pi
            cnfpi=cnf*m_pi
            cr=2.d0/cnf*sin(cnfpi*.5d0)**2
            if(acnf .lt. 1.d-3)then
              cg1=(-euler+cnf**2*(0.042002635034095236d0+
     $             0.042197734555544337d0*cnf**2))
            else
              cg1=(1.d0/cgamm2-1.d0/cgamm1)*.5d0/cnf
            endif
            ca1=cnfpi/sin(cnfpi)
            csigma=cnf*clogz
            if(csigma .eq. (0.d0,0.d0))then
              ca2=(1.d0,0.d0)
            else
              ca2=tcsinh(csigma)/csigma
            endif
            ca3=tccosh(csigma)
            cg2=.5d0*(1.d0/cgamm2+1.d0/cgamm1)
          endif
          cf=2.d0/m_pi*ca1*(ca3*cg1+ca2*clogz*cg2)
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

      recursive real*8 pure function rrbessely(cn,z)
     $     result(cb)
      use macmath
      use mathfun
      implicit none
c     Including m_pi_2 symbol
      integer*4 i,itmax,m
      parameter (itmax=26)
      real*8 ,intent(in):: cn,z
      real*8 cg1,cg2,cgamm1,cgamm2,ca3,
     $     cp,cq,cr,cs0,cs1,zp,cx,cxp,csigma,cx2(2),
     $     ca1,ca2,cf,cg,cf1,cf2,clogz,zhi,anfpi,cs
      real*8 az,an0,anf,aanf,ak
      az=abs(z)
      if(az .gt. 12.2d0)then
        cx2=rrbesasym(cn,z)
        cb=cx2(2)
      else
        if(cn .lt. 0.d0)then
          cb=-sin(m_pi*cn)*rrbesselj(-cn,z)+
     $         cos(m_pi*cn)*rrbessely(-cn,z)
        else
          an0=aint(cn)
          anf=cn-an0
          if(anf .gt. .5d0)then
            anf=anf-1.d0
            an0=an0+1.d0
          endif
          aanf=abs(anf)
          zhi=2.d0/z
          clogz=log(zhi)
          if(aanf .eq. 0.d0)then
            zp=1.d0
            cp=1.d0/m_pi
            cq=1.d0/m_pi
            cg1=-euler
            cg2=1.d0
            ca1=1.d0
            ca2=1.d0
            ca3=1.d0
            cr=0.d0
          else
            zp=exp(anf*clogz)
            cgamm1=gamma(1.d0+anf)
            cgamm2=gamma(1.d0-anf)
            cp=cgamm1*zp/m_pi
            cq=cgamm2/zp/m_pi
            anfpi=anf*m_pi
            cr=2.d0/anf*sin(anfpi*.5d0)**2
            if(aanf .lt. 1.d-3)then
              cg1=(-euler+anf**2*(0.042002635034095236d0+
     $             0.042197734555544337d0*anf**2))
            else
              cg1=(1.d0/cgamm2-1.d0/cgamm1)*.5d0/anf
            endif
            ca1=anfpi/sin(anfpi)
            csigma=anf*clogz
            if(csigma .eq. 0.d0)then
              ca2=1.d0
            else
              ca2=sinh(csigma)/csigma
            endif
            ca3=cosh(csigma)
            cg2=.5d0*(1.d0/cgamm2+1.d0/cgamm1)
          endif
          cf=2.d0/m_pi*ca1*(ca3*cg1+ca2*clogz*cg2)
          cs0=cf+cr*cq
          cs1=cp
          ak=0.d0
          cx=-(z*.5d0)**2
          cxp=1.d0
          m=min(itmax,max(int(az*2.3d0),10))
          do i=1,m
            ak=ak+1.d0
            cxp=cxp*cx/ak
            cf1=1.d0/(ak-anf)
            cf2=1.d0/(ak+anf)
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
              anf=anf+1.d0
              cs=anf*zhi*cs1-cs0
              cs0=cs1
              cs1=cs
            enddo
            cb=cs1
          endif
        endif
      endif
      return
      end

      recursive complex*16 pure function cbesselj(cn,z)
     $     result(cb)
      use gammaf
      implicit none
      integer*4 nmax,i
      complex*16 ,intent(in):: cn,z
      complex*16 cj,cj1,cj2,cs,zi,cnk,cn2k,cg,cb2(2)
      real*8 az,acn,an,ak,xcn
      if(imag(cn) .eq. 0.d0 .and. imag(z) .eq. 0.d0)then
        if(dble(z) .ge. 0.d0)then
          cb=dcmplx(rrbesselj(dble(cn),dble(z)),0.d0)
        else
          cb=rrbesselj(dble(cn),-dble(z))*dcmplx(-1.d0,0.d0)**dble(cn)
        endif
        return
      endif
c     Initialization to avoid compiler warning      
      nmax=0
c     
      az=abs(z)
      acn=abs(cn)
      if(az .gt. acn*.25d0+20.d0)then
        cb2=cbesasym(cn,z)
        cb=cb2(1)
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
          if(xcn .lt. 0.d0 .and. cn .eq. anint(xcn))then
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
            cg=exp(cloggamma1(cnk)-aloggamma1(ak))/cnk
c            cg=cgamma(cnk)/factorial(ak)
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

      complex*16 pure function cbesasym(cn,z) result(cx)
      use macmath
      implicit none
c     Including m_pi_2/m_pi_4/m_2_pi symbol
      dimension cx(2)
      complex*16 ,intent(in):: cn,z
      complex*16 cp,cq,cmu,z1,chi,csqrtz,ccoschi,csinchi
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
      cx=csqrtz*[cp*ccoschi-cq*csinchi,cp*csinchi+cq*ccoschi]
      return
      end

      recursive real*8 pure function rrbesselj(cn,z)
     $     result(cb)
      use gammaf
      implicit none
      integer*4 nmax,i
      real*8 ,intent(in):: cn,z
      real*8 cj,cj1,cj2,cs,zi,cnk,cn2k,cg,az,acn,an,ak,cx(2)
c     Initialization to avoid compiler warning
      nmax=0
c     
      az=abs(z)
      acn=abs(cn)
      if(az .gt. acn*.25d0+20.d0)then
        cx=rrbesasym(cn,z)
        cb=cx(1)
      else
        if(az .ne. 0.d0)then
          nmax=int(max(1.d0,20.d0/max(log10(2.d0/az),1.d0),
     $         2.d0*az-dble(cn)))
        endif
        if(cn .eq. 0.d0)then
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
          if(cn .lt. 0.d0 .and. cn .eq. anint(cn))then
            if(mod(-cn,2.d0) .eq. 0.d0)then
              cb=rrbesselj(-cn,z)
            else
              cb=-rrbesselj(-cn,z)
            endif
          elseif(az .eq. 0.d0)then
            cb=0.d0
          else
            cj2=0.d0
            cj1=1.d-100
            ak=dble(nmax)
            cnk=cn+ak
            cn2k=cnk+ak
            cg=exp(aloggamma1(cnk)-aloggamma1(ak))/cnk
c            cg=cgamma(cnk)/factorial(ak)
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

      real*8 pure function rrbesasym(cn,z) result(cx)
      use macmath
      implicit none
c     Including m_pi_2/m_pi_4/m_2_pi symbol
      dimension cx(2)
      real*8 ,intent(in):: cn,z
      real*8 cp,cq,z1,chi,csqrtz,ccoschi,csinchi,cmu
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
      cx=csqrtz*[cp*ccoschi-cq*csinchi,cp*csinchi+cq*ccoschi]
      return
      end

      end module
