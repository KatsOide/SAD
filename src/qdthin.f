      subroutine qdthin(dtrans,dcod,nord,al,ak,
     $     k1,idp,dx,dy,theta,iv,nut)
      use ffs
      use ffs_pointer
      use tffitcode
      use multa, only:fact
      implicit none
      integer*4 ,intent(in):: nut,k1,idp,iv,nord
      integer*4 kord
      real*8 ,intent(inout):: dtrans(4,5),dcod(6)
      real*8 ,intent(in):: al,ak,dx,dy
      real*8 trans1(4,5),cod(6),
     $     theta,ala,alb,pr,aki,daki,
     $     a21,a23,a41,a43,a26,a46,
     $     da21,da23,da41,da43,da26,da46
      complex*16 cx,cx1,cx0
      if(iv .eq. 4)then
        call qdrotate(dtrans,dcod,k1,itwissp(k1),idp,dx,dy,nut)
      else
        call qtentu(trans1,cod,utwiss(1,idp,itwissp(k1)),.true.)
        call qchg(trans1,cod,-dx,-dy,theta,.true.)
c     begin initialize for preventing compiler warning
        ala=0.d0
        alb=0.d0
c     end   initialize for preventing compiler warning
        if(al .ne. 0.d0)then
          ala=al/6.d0
          alb=al/1.5d0
          cod(1)=cod(1)+cod(2)*ala
          cod(3)=cod(3)+cod(4)*ala
          trans1(1,:)=trans1(1,:)+trans1(2,:)*ala
          trans1(3,:)=trans1(3,:)+trans1(4,:)*ala
        endif
        pr=1.d0+cod(6)
        kord=nord/2-1
        daki=merge(1.d0,0.5d0,al .eq. 0.d0)/fact(kord)/pr
        aki=ak*daki
        cx0=dcmplx(cod(1),-cod(3))
        cx=merge((1.d0,0.d0),cx0**kord,kord .le. 0)
        cod(2)=cod(2)-aki*dble(cx)
        cod(4)=cod(4)-aki*imag(cx)
        dcod(1)=0.d0
        dcod(2)=-daki*dble(cx)
        dcod(3)=0.d0
        dcod(4)=-daki*imag(cx)
        if(kord .gt. 0)then
          cx1=merge(kord*cx0**(kord-1),(1.d0,0.d0),
     $         kord .gt. 1)
          a21=-aki*dble(cx1)
          a23=-aki*imag(cx1)
          a41=-aki*imag(cx1)
          a43= aki*dble(cx1)
          a26=aki/pr*dble(cx)
          a46=aki/pr*imag(cx)
          da21=-daki*dble(cx1)
          da23=-daki*imag(cx1)
          da41=-daki*imag(cx1)
          da43= daki*dble(cx1)
          da26=daki/pr*dble(cx)
          da46=daki/pr*imag(cx)
          dtrans(1,:)=0.d0
          dtrans(2,:)=da21*trans1(1,:)+da23*trans1(3,:)
          dtrans(3,:)=0.d0
          dtrans(4,:)=da41*trans1(1,:)+da43*trans1(3,:)
          trans1(2,:)=trans1(2,:)+ a21*trans1(1,:)+ a23*trans1(3,:)
          trans1(4,:)=trans1(4,:)+ a41*trans1(1,:)+ a43*trans1(3,:)
          dtrans(2,5)=dtrans(2,5)+da26
          dtrans(4,5)=dtrans(4,5)+da46
          trans1(2,5)=trans1(2,5)+a26
          trans1(4,5)=trans1(4,5)+a46
        endif
        if(al .eq. 0.d0)then
          go to 3000
        endif
        cod(1)=cod(1)+cod(2)*alb
        cod(3)=cod(3)+cod(4)*alb
        dcod(1)=dcod(1)+dcod(2)*alb
        dcod(3)=dcod(3)+dcod(4)*alb
        trans1(1,:)=trans1(1,:)+trans1(2,:)*alb
        trans1(3,:)=trans1(3,:)+trans1(4,:)*alb
        dtrans(1,:)=dtrans(1,:)+dtrans(2,:)*alb
        dtrans(3,:)=dtrans(3,:)+dtrans(4,:)*alb
        cx0=dcmplx(cod(1),-cod(3))
        cx=merge((1.d0,0.d0),cx0**kord,kord .le. 0)
        dcod(2)=dcod(2)-daki*dble(cx)
        dcod(4)=dcod(4)-daki*imag(cx)
        if(kord .gt. 0)then
          cx1=merge(kord*cx0**(kord-1),(1.d0,0.d0),
     $         kord .gt. 1)
          a21=-aki*dble(cx1)
          a23=-aki*imag(cx1)
          a41=-aki*imag(cx1)
          a43= aki*dble(cx1)
          da21=-daki*dble(cx1)
          da23=-daki*imag(cx1)
          da41=-daki*imag(cx1)
          da43= daki*dble(cx1)
          da26=daki/pr*dble(cx)
          da46=daki/pr*imag(cx)
          dtrans(2,:)=dtrans(2,:)+ a21*dtrans(1,:)+ a23*dtrans(3,:)
     1         +da21*trans1(1,:)+da23*trans1(3,:)
          dtrans(4,:)=dtrans(4,:)+ a41*dtrans(1,:)+ a43*dtrans(3,:)
     1         +da41*trans1(1,:)+da43*trans1(3,:)
          dtrans(2,5)=dtrans(2,5)+da26
          dtrans(4,5)=dtrans(4,5)+da46
          dcod(2)=dcod(2)+a21*dcod(1)+a23*dcod(3)
          dcod(4)=dcod(4)+a41*dcod(1)+a43*dcod(3)
        endif
        dcod(1)=dcod(1)+dcod(2)*ala
        dcod(3)=dcod(3)+dcod(4)*ala
        dtrans(1,:)=dtrans(1,:)+dtrans(2,:)*ala
        dtrans(3,:)=dtrans(3,:)+dtrans(4,:)*ala
 3000   continue
        call qchg(dtrans,dcod,0.d0,0.d0,-theta,.false.)
      endif
      return
      end
