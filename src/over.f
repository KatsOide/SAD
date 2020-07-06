      logical function over(cs,csa,dcs,scale)
c      judge if two character strings overlap or not
      implicit real*8 (a-h,o-z)
      dimension cs(7),csa(7)
      over=.False.
      cl =1.5* 0.042857*cs(7)*scale
      cla=1.5* 0.042857*csa(7)*scale
      clu=1.2*  0.038462*cs(7)*scale
      clua=1.2* 0.038462*csa(7)*scale
      cld =1.2* 0.049444*cs(7)*scale
      clda=1.2* 0.049444*csa(7)*scale
      if(cs(4).ne.csa(4)) then
        dcs=0
        over=.False.
      elseif(cs(4).eq.-1..and.cs(6).ge.0.
     &                   .or.cs(4).eq.1. .and.cs(6).lt.0.) then
        dcs=sqrt((cs(1)-csa(1))**2+(cs(2)-csa(2))**2)-(clua+cld)
        if( dcs.ge.0) then
          over=.False.
        else
          dcs=sqrt(-dcs)/scale/0.049444/1.2
          over=.True.
        endif
      elseif(cs(4).eq.-1..and.cs(6).lt.0.
     &                   .or.cs(4).eq.1. .and.cs(6).ge.0.) then
        ln=int(min(cs(5),csa(5)))
        x=cs(1)+cl *2.*ln*cos(cs(3))
        y=cs(2)+cl *2.*ln*sin(cs(3))
        xa=csa(1)+cla*2.*ln*cos(csa(3))
        ya=csa(2)+cla*2.*ln*sin(csa(3))
        dcs=sqrt( (x-xa)**2+(y-ya)**2 )-(clua+cld)
        if(dcs.ge.0) then
          over=.False.
        else
          dcs=sqrt(-dcs)/scale/0.049444/1.2
          over=.True.
        endif
      endif
      return
      end
