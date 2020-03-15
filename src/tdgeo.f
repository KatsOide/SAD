      subroutine tdgeo(couple1,qu,kf,ne,k,ltyp,iv,nut,nfam)
      use kyparam
      use tfstk
      use ffs
      use ffs_pointer
      use tffitcode
      use mathfun
      implicit none
      integer*4 kf,ne,k,ltyp,iv,nut,nfam,ip,j,j1,j2,
     $     k1
      integer*4 kbe,kbg
      real*8 dtrans(36),dcod(6),dcod0(6),
     $     trans1(4,5),s,couple1,qu,theta,al,drj1,drj2,c1,c2,
     $     cost,sint,yr1,yr2,yr3,xr,dz1,dz2,dx3,dy3,dx,dy,
     $     dr1,dr2,dr3,v
      if(ne .eq. iorgr)then
        return
      endif
      if(ne .gt. iorgr .and. (k .ge. ne .or. k .lt. iorgr) .or.
     1   ne .lt. iorgr .and. (k .lt. ne .or. k .ge. iorgr))then
        return
      endif
      s=couple1*sign(1.d0,dble(ne-iorgr))
      if(ltyp .eq. 1)then
        if(kf .ge. mfitgx .and. kf .le. mfitgz)then
          qu=qu+s*geo(kf-mfitgx+1,3,k)
        endif
      elseif(ltyp .eq. icBEND)then
        if(iv .eq. ky_ANGL_BEND)then
          ip=idelc(k)
          theta=rlist(idval(ip)+5)
          if(kf .ge. mfitgx .and. kf .le. mfitgz)then
            al=rlist(idval(ip)+1)
            j=kf-mfitgx+1
            j1=mod(j,3)+1
            j2=6-j-j1
            drj1=geo(j1,4,ne)-geo(j1,4,k+1)
            drj2=geo(j2,4,ne)-geo(j2,4,k+1)
            v=rlist(latt(k)+2)
            if(v .eq. 0.d0)then
              c1=0.d0
              c2=-.5d0
            else
              c1=sinc(v)/v**2
              c2=(sin(v*.5d0)-v*cos(v*.5d0))*2.d0*sin(v*.5d0)/v**2
            endif
            if(theta .ne. 0.d0)then
              cost=cos(theta)
              sint=sin(theta)
              yr1=-cost*geo(j1,2,k)-sint*geo(j1,1,k)
              yr2=-cost*geo(j2,2,k)-sint*geo(j2,1,k)
              xr =cost*geo(j,1,k)-sint*geo(j,2,k)
            else
              yr1=-geo(j1,2,k)
              yr2=-geo(j2,2,k)
              xr =geo(j,1,k)
            endif
            qu=qu+(yr1*drj2-yr2*drj1
     1           +al*(c1*geo(j,3,k)+c2*xr))*s
          elseif(kf .ge. mfitchi1 .and. kf .le. mfitchi3)then
            if(theta .ne. 0.d0)then
              cost=cos(theta)
              sint=sin(theta)
              yr1=-cost*geo(1,2,k)-sint*geo(1,1,k)
              yr2=-cost*geo(2,2,k)-sint*geo(2,1,k)
              yr3=-cost*geo(3,2,k)-sint*geo(3,1,k)
            else
              yr1=-geo(1,2,k)
              yr2=-geo(2,2,k)
              yr3=-geo(3,2,k)
            endif
            select case (kf)
            case (mfitchi2)
              qu=qu+s*(yr1*geo(2,3,ne)-yr2*geo(1,3,ne))
     $             /hypot(geo(1,3,ne),geo(2,3,ne))
c     $             /sqrt(geo(1,3,ne)**2+geo(2,3,ne)**2)
            case (mfitchi1)
              dz1=-geo(2,3,ne)*yr3+geo(3,3,ne)*yr2
              dz2=-geo(3,3,ne)*yr1+geo(1,3,ne)*yr3
              qu=qu+(geo(1,3,ne)*dz2-geo(2,3,ne)*dz1)
     $             /(geo(1,3,ne)**2+geo(2,3,ne)**2)*s
            case default
              dx3=-geo(1,1,ne)*yr2+geo(2,1,ne)*yr1
              dy3=-geo(1,2,ne)*yr2+geo(2,2,ne)*yr1
              qu=qu+(-geo(3,2,ne)*dx3+geo(3,1,ne)*dy3)
     $             /(geo(3,1,ne)**2+geo(3,2,ne)**2)*s
            end select
          endif
        elseif(iv .eq. ky_L_BEND)then
          if(kf .ge. mfitgx .and. kf .le. mfitgz)then
            qu=qu+s*geo(kf-mfitgx+1,3,k+1)
          endif
        endif
      elseif(ltyp .le. icDODECA .or. ltyp .eq. icMULT)then
        call tfbndsol(k,kbg,kbe)
        if(kbe .ne. 0)then
          ip=idelc(k)
          dx=rlist(idval(ip)+kytbl(kwDX,ltyp))
          dy=rlist(idval(ip)+kytbl(kwDY,ltyp))
          k1=k+1
          call qdtrans(ltyp,itwissp(k),k,k1,
     $         iv,dtrans,dcod0,0,nfam,nut)
c      write(*,'(a,4i5,1p8g13.5)')'qdtrans ',k1,kbe,
c     $         itwissp(k),itwissp(kbe),
c     $         utwiss(mfitr1:mfitr4,0,itwissp(k1)),
c     $         utwiss(mfitr1:mfitr4,0,itwissp(kbe))
          call tftmatu(
     $         utwiss(1,0,itwissp(k1)),utwiss(1,0,itwissp(kbe)),
     $         0.d0,0.d0,
     $         trans1,k1,kbe,.false.,.false.)
          dcod(1:4)=trans1(1:4,1)*dcod0(1)+trans1(1:4,2)*dcod0(2)
     $             +trans1(1:4,3)*dcod0(3)+trans1(1:4,4)*dcod0(4)
          yr1=-dcod(4)*geo(1,1,kbe)+dcod(2)*geo(1,2,kbe)
          yr2=-dcod(4)*geo(2,1,kbe)+dcod(2)*geo(2,2,kbe)
          yr3=-dcod(4)*geo(3,1,kbe)+dcod(2)*geo(3,2,kbe)
          dr1=geo(1,4,ne)-geo(1,4,kbe)
          dr2=geo(2,4,ne)-geo(2,4,kbe)
          dr3=geo(3,4,ne)-geo(3,4,kbe)
          select case (kf)
          case (mfitgx)
            qu=qu+s*(-dr2*yr3+dr3*yr2
     $           +dcod(1)*geo(1,1,kbe)+dcod(3)*geo(1,2,kbe))
          case (mfitgy)
            qu=qu+s*(-dr3*yr1+dr1*yr3
     $           +dcod(1)*geo(2,1,kbe)+dcod(3)*geo(2,2,kbe))
          case (mfitgz)
            qu=qu+s*(-dr1*yr2+dr2*yr1
     $           +dcod(1)*geo(3,1,kbe)+dcod(3)*geo(3,2,kbe))
          case (mfitchi2)
            qu=qu+s*(yr1*geo(2,3,ne)-yr2*geo(1,3,ne))
     $           /hypot(geo(1,3,ne),geo(2,3,ne))
c     $           /sqrt(geo(1,3,ne)**2+geo(2,3,ne)**2)
          case (mfitchi1)
            dz1=-geo(2,3,ne)*yr3+geo(3,3,ne)*yr2
            dz2=-geo(3,3,ne)*yr1+geo(1,3,ne)*yr3
            qu=qu+(geo(1,3,ne)*dz2-geo(2,3,ne)*dz1)
     $           /(geo(1,3,ne)**2+geo(2,3,ne)**2)*s
          case default
            dx3=-geo(1,1,ne)*yr2+geo(2,1,ne)*yr1
            dy3=-geo(1,2,ne)*yr2+geo(2,2,ne)*yr1
            qu=qu+(-geo(3,2,ne)*dx3+geo(3,1,ne)*dy3)
     $           /(geo(3,1,ne)**2+geo(3,2,ne)**2)*s
          end select
        endif           
      endif
      return
      end
