      subroutine qdthin(dtrans,dcod,nord,al,ak,
     $     k1,idp,dx,dy,theta,iv,nut)
      use ffs
      use ffs_pointer
      use tffitcode
      implicit none
      integer*4 nut,k1,idp,iv,nord,i,kord
      real*8 trans1(4,5),dtrans(4,5),dcod(6),cod(6),al,ak,
     $     dx,dy,theta,ala,alb,pr,aki,daki,
     $     a21,a23,a41,a43,a26,a46,
     $     da21,da23,da41,da43,da26,da46
      complex*16 cx,cx1,cx0
      real*8 fact(0:10)
      data fact / 1.d0,  1.d0,   2.d0,   6.d0,   24.d0,   120.d0,
     1          720.d0,5040.d0,40320.d0,362880.d0,3628800.d0 /
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
          do 1110 i=1,5
            trans1(1,i)=trans1(1,i)+trans1(2,i)*ala
            trans1(3,i)=trans1(3,i)+trans1(4,i)*ala
 1110     continue
        endif
        pr=1.d0+cod(6)
        kord=nord/2-1
        if(al .eq. 0.d0)then
          daki=1.d0/fact(kord)/pr
        else
          daki=1.d0/fact(kord)/pr*.5d0
        endif
        aki=ak*daki
        cx0=dcmplx(cod(1),-cod(3))
        if(kord .le. 0)then
          cx=(1.d0,0.d0)
        else
          cx=cx0**kord
        endif
        cod(2)=cod(2)-aki*dble(cx)
        cod(4)=cod(4)-aki*imag(cx)
        dcod(1)=0.d0
        dcod(2)=-daki*dble(cx)
        dcod(3)=0.d0
        dcod(4)=-daki*imag(cx)
        if(kord .gt. 0)then
          if(kord .gt. 1)then
            cx1=kord*cx0**(kord-1)
          else
            cx1=(1.d0,0.d0)
          endif
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
          do 10 i=1,5
            dtrans(1,i)=0.d0
            dtrans(2,i)=da21*trans1(1,i)+da23*trans1(3,i)
            dtrans(3,i)=0.d0
            dtrans(4,i)=da41*trans1(1,i)+da43*trans1(3,i)
            trans1 (2,i)=trans1(2,i)+ a21*trans1(1,i)+ a23*trans1(3,i)
            trans1 (4,i)=trans1(4,i)+ a41*trans1(1,i)+ a43*trans1(3,i)
 10       continue
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
        do 1120 i=1,5
          trans1(1,i)=trans1(1,i)+trans1(2,i)*alb
          trans1(3,i)=trans1(3,i)+trans1(4,i)*alb
          dtrans(1,i)=dtrans(1,i)+dtrans(2,i)*alb
          dtrans(3,i)=dtrans(3,i)+dtrans(4,i)*alb
 1120   continue
        cx0=dcmplx(cod(1),-cod(3))
        if(kord .le. 0)then
          cx=(1.d0,0.d0)
        else
          cx=cx0**kord
        endif
        dcod(2)=dcod(2)-daki*dble(cx)
        dcod(4)=dcod(4)-daki*imag(cx)
        if(kord .gt. 0)then
          if(kord .gt. 1)then
            cx1=kord*cx0**(kord-1)
          else
            cx1=(1.d0,0.d0)
          endif
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
          do 20 i=1,5
            dtrans(2,i)=dtrans(2,i)+ a21*dtrans(1,i)+ a23*dtrans(3,i)
     1           +da21*trans1(1,i)+da23*trans1(3,i)
            dtrans(4,i)=dtrans(4,i)+ a41*dtrans(1,i)+ a43*dtrans(3,i)
     1           +da41*trans1(1,i)+da43*trans1(3,i)
 20       continue
          dtrans(2,5)=dtrans(2,5)+da26
          dtrans(4,5)=dtrans(4,5)+da46
          dcod(2)=dcod(2)+a21*dcod(1)+a23*dcod(3)
          dcod(4)=dcod(4)+a41*dcod(1)+a43*dcod(3)
        endif
        dcod(1)=dcod(1)+dcod(2)*ala
        dcod(3)=dcod(3)+dcod(4)*ala
        do 1130 i=1,5
          dtrans(1,i)=dtrans(1,i)+dtrans(2,i)*ala
          dtrans(3,i)=dtrans(3,i)+dtrans(4,i)*ala
 1130   continue
 3000   continue
        call qchg(dtrans,dcod,0.d0,0.d0,-theta,.false.)
      endif
      return
      end
