      subroutine TWTRANS(R0,R0I)
      use tmacro
      implicit real*8(a-h,o-z)
      integer*4 ntwissfun
      parameter (ntwissfun=20)
      COMMON /TWISSP/TWPIN(ntwissfun),TWPOT(ntwissfun),
     +   AX,BX,PX,AY,BY,PY,EX,EPX,EY,EPY,R1,R2,R3,R4,DX(4)
      COMMON /TMATR/TM(4,4),VM(4),TO(4,4),VO(4),TU(4,4),VU(4)
      DIMENSION TTX(3,3),TTY(3,3),R0(4,4),R0I(4,4),DX1(4)
      TTX(1,1)=TU(1,1)*TU(1,1)
      TTX(1,2)=-2.*TU(1,1)*TU(1,2)
      TTX(1,3)=TU(1,2)*TU(1,2)
      TTX(2,1)=-TU(1,1)*TU(2,1)
      TTX(2,2)=1.+2.*TU(1,2)*TU(2,1)
      TTX(2,3)=-TU(1,2)*TU(2,2)
      TTX(3,1)=TU(2,1)*TU(2,1)
      TTX(3,2)=-2.*TU(2,2)*TU(2,1)
      TTX(3,3)=TU(2,2)*TU(2,2)
      TTY(1,1)=TU(3,3)*TU(3,3)
      TTY(1,2)=-2.*TU(3,3)*TU(3,4)
      TTY(1,3)=TU(3,4)*TU(3,4)
      TTY(2,1)=-TU(3,3)*TU(4,3)
      TTY(2,2)=1.+2.*TU(3,4)*TU(4,3)
      TTY(2,3)=-TU(3,4)*TU(4,4)
      TTY(3,1)=TU(4,3)*TU(4,3)
      TTY(3,2)=-2.*TU(4,4)*TU(4,3)
      TTY(3,3)=TU(4,4)*TU(4,4)
      GX=(1+AX*AX)/BX
      GY=(1+AY*AY)/BY
      BX1=TTX(1,1)*BX+TTX(1,2)*AX+TTX(1,3)*GX
      AX1=TTX(2,1)*BX+TTX(2,2)*AX+TTX(2,3)*GX
*     GX1=TTX(3,1)*BX+TTX(3,2)*AX+TTX(3,3)*GX
*     GX1=(1+AX1*AX1)/BX1
      BY1=TTY(1,1)*BY+TTY(1,2)*AY+TTY(1,3)*GY
      AY1=TTY(2,1)*BY+TTY(2,2)*AY+TTY(2,3)*GY
*     GY1=TTY(3,1)*BY+TTY(3,2)*AY+TTY(3,3)*GY
*     GY1=(1+AY1*AY1)/BY1
      EX1=TU(1,1)*EX+TU(1,2)*EPX+VU(1)
      EPX1=TU(2,1)*EX+TU(2,2)*EPX+VU(2)
      EY1=TU(3,3)*EY+TU(3,4)*EPY+VU(3)
      EPY1=TU(4,3)*EY+TU(4,4)*EPY+VU(4)
      AX=AX1
      AY=AY1
      BX=BX1
      BY=BY1
      EX=EX1
      EPX=EPX1
      EY=EY1
      EPY=EPY1
      CALL VMUL4(TO,DX,DX1)
      DX(1)=DX1(1)+VO(1)*DP0
*    +      DP0*(R0I(1,1)*EX+R0I(1,2)*EPX+R0I(1,3)*EY+R0I(1,4)*EPY)
      DX(2)=DX1(2)+VO(2)*DP0
*    +      DP0*(R0I(2,1)*EX+R0I(2,2)*EPX+R0I(2,3)*EY+R0I(2,4)*EPY)
      DX(3)=DX1(3)+VO(3)*DP0
*    +      DP0*(R0I(3,1)*EX+R0I(3,2)*EPX+R0I(3,3)*EY+R0I(3,4)*EPY)
      DX(4)=DX1(4)+VO(4)*DP0
*    +      DP0*(R0I(4,1)*EX+R0I(4,2)*EPX+R0I(4,3)*EY+R0I(4,4)*EPY)
      END
