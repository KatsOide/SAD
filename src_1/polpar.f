      subroutine polpar(id0,ld,al,ak00,ak0,ak1,az,cod)
      use tfstk
      use ffs
      use tffitcode
      implicit real*8(a-h,o-z)
      real*8 cod(6)
      if(irad .eq. 6)then
        return
      endif
      ilist(1,ipolid+ipelm-1)=id0
      ilist(2,ipolid+ipelm-1)=ld
      rlist(ipoll+ipelm-1)=al
      icod=ipolo+(ipelm-1)*6
      call tmov(cod,rlist(icod),6)
      x=cod(1)
      y=cod(3)
      id=id0/10
      if(id .gt. 6 .or. id .eq. 0)then
        return
      endif
      ib=ipolb+(ipelm-1)*21
      if(id .eq. 1)then
        bx  =0.d0
        bxx =0.d0
        bxy =0.d0
        bxxy=0.d0
        by0 =0.d0
        by  =0.d0
        byx =0.d0
        byy =0.d0
        byxx=0.d0
        byyy=0.d0
        bz  =az*brho
      elseif(id .eq. 2)then
        b   =ak1*brho
        bx  = b*y
        bxx =0.d0
        bxy = b
        bxxy=0.d0
        by0 = ak00*brho
        by  = b*x+ak0*brho
        byx = b
        byy =0.d0
        byxx=0.d0
        byyy=0.d0
        bz  =az*brho
      elseif(id .eq. 4)then
        b=ak1*brho
        bx  = b*y
        bxx =0.d0
        bxy = b
        bxxy=0.d0
        by0 =0.d0
        by  = b*x
        byx = b
        byy =0.d0
        byxx=0.d0
        byyy=0.d0
        bz  =az*brho
      elseif(id .eq. 6)then
        b=ak1*brho
        bx  = b*x*y
        bxx = b*y
        bxy = b*x
        bxxy= b
        by0 = 0.d0
        by  = b*(x-y)*(x+y)*.5d0
        byx = b*x
        byy =-b*y
        byxx= b
        byyy=-b
        bz  =0.d0
      else
c     begin initialize for preventing compiler warning
        bx  =0.d0
        bxx =0.d0
        bxy =0.d0
        bxxy=0.d0
        by0 =0.d0
        by  =0.d0
        byx =0.d0
        byy =0.d0
        byxx=0.d0
        byyy=0.d0
        bz  =0.d0
c     end   initialize for preventing compiler warning
      endif
      rlist(ib   )=0.d0
      rlist(ib+1 )=bx
      rlist(ib+2 )=bxx
      rlist(ib+3 )=bxy
      rlist(ib+4 )=0.d0
      rlist(ib+5 )=bxxy
      rlist(ib+6 )=0.d0
      rlist(ib+7 )=by0
      rlist(ib+8 )=by
      rlist(ib+9 )=byx
      rlist(ib+10)=byy
      rlist(ib+11)=byxx
      rlist(ib+12)=0.d0
      rlist(ib+13)=byyy
      rlist(ib+14)=bz
      rlist(ib+15)=bz
      rlist(ib+16)=0.d0
      rlist(ib+17)=0.d0
      rlist(ib+18)=0.d0
      rlist(ib+19)=0.d0
      rlist(ib+20)=0.d0
      return
      end
