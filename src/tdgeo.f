      subroutine tdgeo(latt,utwiss,itwissp,gammab,
     $couple,qu,kf,ne,k,geo,ltyp,iv,nut,nfam)
      use tfstk
      use ffs
      use tffitcode
      implicit none
      integer*4 kf,ne,k,ltyp,iv,nut,nfam,ip,j,j1,j2,
     $     k1
      integer*4 latt(2,nlat),itwissp(nlat),kbe,kbg
      real*8 geo(3,4,nlat),utwiss(ntwissfun,-nfam:nfam,nut),
     $     gammab(nlat),dtrans(36),dcod(6),dcod0(6),
     $     trans1(4,5),s,couple,qu,theta,al,drj1,drj2,c1,c2,
     $     cost,sint,yr1,yr2,yr3,xr,dz1,dz2,dx3,dy3,dx,dy,
     $     dr1,dr2,dr3,v,sinc
      if(ne .eq. iorgr)then
        return
      endif
      if(ne .gt. iorgr .and. (k .ge. ne .or. k .lt. iorgr) .or.
     1   ne .lt. iorgr .and. (k .lt. ne .or. k .ge. iorgr))then
        return
      endif
      s=couple*sign(1.d0,dble(ne-iorgr))
      if(ltyp .eq. 1)then
        if(kf .ge. mfitgx .and. kf .le. mfitgz)then
          qu=qu+s*geo(kf-mfitgx+1,3,k)
        endif
      elseif(ltyp .eq. icBEND)then
        if(iv .eq. kytbl(kwANGL,icBEND))then
          ip=latt(1,k)
          theta=rlist(idval(ip)+5)
          if(kf .ge. mfitgx .and. kf .le. mfitgz)then
            al=rlist(idval(ip)+1)
            j=kf-mfitgx+1
            j1=mod(j,3)+1
            j2=6-j-j1
            drj1=geo(j1,4,ne)-geo(j1,4,k+1)
            drj2=geo(j2,4,ne)-geo(j2,4,k+1)
            v=rlist(latt(2,k)+2)
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
            j=kf-mfitchi1+1
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
            if(j .eq. 2)then
              qu=qu+s*(yr1*geo(2,3,ne)-yr2*geo(1,3,ne))
     $             /sqrt(geo(1,3,ne)**2+geo(2,3,ne)**2)
            elseif(j .eq. 1)then
              dz1=-geo(2,3,ne)*yr3+geo(3,3,ne)*yr2
              dz2=-geo(3,3,ne)*yr1+geo(1,3,ne)*yr3
              qu=qu+(geo(1,3,ne)*dz2-geo(2,3,ne)*dz1)
     $             /(geo(1,3,ne)**2+geo(2,3,ne)**2)*s
            else
              dx3=-geo(1,1,ne)*yr2+geo(2,1,ne)*yr1
              dy3=-geo(1,2,ne)*yr2+geo(2,2,ne)*yr1
              qu=qu+(-geo(3,2,ne)*dx3+geo(3,1,ne)*dy3)
     $             /(geo(3,1,ne)**2+geo(3,2,ne)**2)*s
            endif
          endif
        elseif(iv .eq. kytbl(kwL,icBEND))then
          if(kf .ge. mfitgx .and. kf .le. mfitgz)then
            qu=qu+s*geo(kf-mfitgx+1,3,k+1)
          endif
        endif
      elseif(ltyp .le. icDODECA .or. ltyp .eq. icMULT)then
        call tfbndsol(k,kbg,kbe)
        if(kbe .ne. 0)then
          ip=latt(1,k)
          dx=rlist(idval(ip)+kytbl(kwDX,ltyp))
          dy=rlist(idval(ip)+kytbl(kwDY,ltyp))
          k1=k+1
          call qdtrans(latt,nlat,ltyp,itwissp(k),k,k1,
     $         iv,dtrans,dcod0,0,
     $         utwiss,gammab,nfam,nut)
          call tftmatu(
     $         utwiss(1,0,itwissp(k1)),utwiss(1,0,itwissp(kbe)),
     $         0.d0,0.d0,
     $         gammab,trans1,k1,kbe,.false.,.false.)
          dcod(1)=trans1(1,1)*dcod0(1)+trans1(1,2)*dcod0(2)
     $           +trans1(1,3)*dcod0(3)+trans1(1,4)*dcod0(4)
          dcod(2)=trans1(2,1)*dcod0(1)+trans1(2,2)*dcod0(2)
     $           +trans1(2,3)*dcod0(3)+trans1(2,4)*dcod0(4)
          dcod(3)=trans1(3,1)*dcod0(1)+trans1(3,2)*dcod0(2)
     $           +trans1(3,3)*dcod0(3)+trans1(3,4)*dcod0(4)
          dcod(4)=trans1(4,1)*dcod0(1)+trans1(4,2)*dcod0(2)
     $           +trans1(4,3)*dcod0(3)+trans1(4,4)*dcod0(4)
          if(kf .ge. mfitgx .and. kf .le. mfitgz)then
            dr1=geo(1,4,ne)-geo(1,4,kbe)
            dr2=geo(2,4,ne)-geo(2,4,kbe)
            dr3=geo(3,4,ne)-geo(3,4,kbe)
            yr1=-dcod(4)*geo(1,1,kbe)+dcod(2)*geo(1,2,kbe)
            yr2=-dcod(4)*geo(2,1,kbe)+dcod(2)*geo(2,2,kbe)
            yr3=-dcod(4)*geo(3,1,kbe)+dcod(2)*geo(3,2,kbe)
            j=kf-mfitgx+1
            if(j .eq. 1)then
              qu=qu+s*(-dr2*yr3+dr3*yr2
     $             +dcod(1)*geo(1,1,kbe)+dcod(3)*geo(1,2,kbe))
            elseif(j .eq. 2)then
              qu=qu+s*(-dr3*yr1+dr1*yr3
     $             +dcod(1)*geo(2,1,kbe)+dcod(3)*geo(2,2,kbe))
            else
              qu=qu+s*(-dr1*yr2+dr2*yr1
     $             +dcod(1)*geo(3,1,kbe)+dcod(3)*geo(3,2,kbe))
            endif
          elseif(kf .ge. mfitchi1 .and. kf .le. mfitchi3)then
            yr1=-dcod(4)*geo(1,1,kbe)+dcod(2)*geo(1,2,kbe)
            yr2=-dcod(4)*geo(2,1,kbe)+dcod(2)*geo(2,2,kbe)
            yr3=-dcod(4)*geo(3,1,kbe)+dcod(2)*geo(3,2,kbe)
            j=kf-mfitchi1+1
            if(j .eq. 2)then
              qu=qu+s*(yr1*geo(2,3,ne)-yr2*geo(1,3,ne))
     $             /sqrt(geo(1,3,ne)**2+geo(2,3,ne)**2)
            elseif(j .eq. 1)then
              dz1=-geo(2,3,ne)*yr3+geo(3,3,ne)*yr2
              dz2=-geo(3,3,ne)*yr1+geo(1,3,ne)*yr3
              qu=qu+(geo(1,3,ne)*dz2-geo(2,3,ne)*dz1)
     $             /(geo(1,3,ne)**2+geo(2,3,ne)**2)*s
            else
              dx3=-geo(1,1,ne)*yr2+geo(2,1,ne)*yr1
              dy3=-geo(1,2,ne)*yr2+geo(2,2,ne)*yr1
              qu=qu+(-geo(3,2,ne)*dx3+geo(3,1,ne)*dy3)
     $             /(geo(3,1,ne)**2+geo(3,2,ne)**2)*s
            endif
          endif
        endif           
      endif
      return
      end
