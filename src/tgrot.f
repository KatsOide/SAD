      module geolib
      real*8 ,parameter :: geoini(3,3)=
     $     reshape(
     $     (/0.d0,-1.d0, 0.d0,
     $       0.d0, 0.d0,-1.d0,
     $       1.d0, 0.d0, 0.d0/) , (/3,3/))

      contains
      real*8 function tfchi(geo,i)
      implicit none
      real*8 , intent(in)::geo(3,3)
      integer*4 ,intent(in):: i
      if(i .eq. 2)then
        tfchi=asin(geo(3,3))
      elseif(i .eq. 1)then
        if(geo(2,3) .eq. 0.d0)then
          tfchi=0.d0
        else
          tfchi=atan(geo(2,3),geo(1,3))
        endif
      else
        if(geo(3,1) .eq. 0.d0)then
          tfchi=0.d0
        else
          tfchi=atan(geo(3,1),-geo(3,2))
        endif
      endif
      return
      end function

      real*8 function tfchitogeo(chi) result(geo)
      implicit none
      real*8 , intent(in)::chi(3)
      dimension geo(3,3)
      real*8 schi1,cchi1,schi2,cchi2,schi3,cchi3
      cchi1=cos(chi(1))
      cchi2=cos(chi(2))
      cchi3=cos(chi(3))
      schi1=sin(chi(1))
      schi2=sin(chi(2))
      schi3=sin(chi(3))
      geo(1,1)= schi1*cchi3-cchi1*schi2*schi3
      geo(2,1)=-cchi1*cchi3-schi1*schi2*schi3
      geo(3,1)= cchi2*schi3
      geo(1,2)= schi1*schi3+cchi1*schi2*cchi3
      geo(2,2)=-cchi1*schi3+schi1*schi2*cchi3
      geo(3,2)=-cchi2*cchi3
      geo(1,3)= cchi1*cchi2
      geo(2,3)= schi1*cchi2
      geo(3,3)= schi2
      return
      end function

      real*8 function tsetg(chi) result(geo)
      implicit none
      dimension geo(3,3)
      real*8 ,intent(in):: chi(3)
      real*8 c1,c2,c3,s1,s2,s3
      c1=cos(chi(1))
      s1=sin(chi(1))
      c2=cos(chi(2))
      s2=sin(chi(2))
      c3=cos(chi(3))
      s3=sin(chi(3))
      geo(1,1)=c1*c2
      geo(2,1)=s1*c2
      geo(3,1)=s2
      geo(1,2)=-c3*s1+c1*s2*s3
      geo(2,2)=c1*c3+s1*s2*s3
      geo(3,2)=-c2*s3
      geo(1,3)=-c1*s2*c3-s1*s3
      geo(2,3)=-s1*s2*c3-c1*s3
      geo(3,3)= c2*c3
      return
      end function

      real*8 function tfrotgeo(geo,chi) result(geo1)
      implicit none
      real*8, intent(in)::geo(3,3),chi(3)
      real*8 cschi1,snchi1,cschi2,snchi2,cschi3,snchi3,
     $     g1(3),g2(3)
      dimension geo1(3,3)
      cschi1=cos(chi(1))
      snchi1=sin(chi(1))
      cschi2=cos(chi(2))
      snchi2=sin(chi(2))
      cschi3=cos(chi(3))
      snchi3=sin(chi(3))
      g1       = geo(:,1) *cschi3+geo(:,2)*snchi3
      geo1(:,2)=-geo(:,1) *snchi3+geo(:,2)*cschi3
      g2       = geo1(:,2)*snchi2+geo(:,3)*cschi2
      geo1(:,2)= geo1(:,2)*cschi2-geo(:,3)*snchi2
      geo1(:,3)=        g2*cschi1+      g1*snchi1
      geo1(:,1)=       -g2*snchi1+      g1*cschi1
      return
      end function

      real*8 function tgrot(geo0,geo1) result(a)
      implicit none
      real*8 ,intent(in)::geo0(3,3),geo1(3,3)
      dimension a(3)
      real*8 a13,a21
      a(2)=asin(
     1     -geo0(1,2)*geo1(1,3)-geo0(2,2)*geo1(2,3)-geo0(3,2)*geo1(3,3))
      a13= geo0(1,1)*geo1(1,3)+geo0(2,1)*geo1(2,3)+geo0(3,1)*geo1(3,3)
      if(a13 .eq. 0.d0)then
        a(1)=0.d0
      else
        a(1)=-atan(a13,
     $       geo0(1,3)*geo1(1,3)+geo0(2,3)*geo1(2,3)
     $       +geo0(3,3)*geo1(3,3))
      endif
      a21= geo0(1,2)*geo1(1,1)+geo0(2,2)*geo1(2,1)+geo0(3,2)*geo1(3,1)
      if(a21 .eq. 0.d0)then
        a(3)=0.d0
      else
        a(3)=atan(-a21,
     $       geo0(1,2)*geo1(1,2)+geo0(2,2)*geo1(2,2)
     $       +geo0(3,2)*geo1(3,2))
      endif
      return
      end function

      real*8 function tfgeofrac(lxp) result(gv)
      use ffs_pointer
      implicit none
      integer*4 , intent(in)::lxp
      dimension gv(3,4)
      real*8 geo1(3,4),pos0,gam0
      geo1=geo(:,:,lxp+1)
      pos0=pos(lxp+1)
      gam0=gammab(lxp+1)
      call tfgeo1(lxp,lxp+1,.true.,.false.)
      gv=geo(:,:,lxp+1)
      geo(:,:,lxp+1)=geo1
      gammab(lxp+1)=gam0
      pos(lxp+1)=pos0
      return
      end function

      real*8 function  tforbitgeo(geo,cod) result(og)
      implicit none
      dimension og(3,4)
      real*8 , intent(in):: geo(3,4),cod(4)
      real*8 chi2i,cchi2i,schi2i,chi1i,cchi1i,schi1i,g1(3)
      chi2i =-asin(min(1.d0,max(-1.d0,cod(4))))
      cchi2i=cos(chi2i)
      schi2i=sin(chi2i)
      chi1i =-asin(min(1.d0,max(-1.d0,cod(2)/cchi2i)))
      cchi1i=cos(chi1i)
      schi1i=sin(chi1i)
      og(:,4)=geo(:,4)+geo(:,1)*cod(1)+geo(:,2)*cod(3)
      og(:,1)= cchi1i*geo(:,1)+schi1i*geo(:,3)
      g1     =-schi1i*geo(:,1)+cchi1i*geo(:,3)
      og(:,3)= cchi2i*g1-schi2i*geo(:,2)
      og(:,2)= schi2i*g1+cchi2i*geo(:,2)
      return
      end function

      end module
