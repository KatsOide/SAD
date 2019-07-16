      subroutine geodrw(geo,word,lfnd,title,case)
      use tfstk
      use ffs
      use tffitcode
      use sad_main
      use ffs_pointer, only:elatt
      use tfcsi
      implicit real*8 (a-h,o-z)
c       idval(latt(1,i)+5) ---> rotation
c       s.t.       x(beam)=x cos(t) - y sin(t)
      parameter (nkey=2)
      dimension geo(3,4,nlat)
      dimension window(4), pn1(3), pn2(3), pn3(3), ds(3)
      dimension cs(7),csa(7)
      character*(MAXPNAME) name
      character*(8) keywrd(nkey)
      character*(*) word,title,case
      character*30 dat
      logical tmatch, over, exist, rght, other
      data keywrd /'RIGHT   ','LEFT    '/
      data window/2.,10.,1.,9./ parity /-1./
      goto 9000
 9000 continue
      call termes(lfno,
     $     'A message from SAD/FFS: replace GEO command '//
     $      'with GeometryPlot[].',
     $  ' ')
      return
c$$$      
      alpha=getva(exist)
      if(.not.exist) then
        alpha=0.
      endif
      beta=getva(exist)
      if(.not.exist) then
        beta=0.
      endif
      gamma=getva(exist)
      if(.not.exist) then
        gamma=0.
      endif
      alpha=alpha*pi/180.
      beta=beta*pi/180.
      gamma=gamma*pi/180.
      ymag=getva(exist)
      if(.not.exist) then
        ymag=1.
      endif
      rght=.true.
      other=.true.
      j=0
    1 j=j+1
      call peekwd(word,next)
      do 10 i=1,2
        if(word.eq.keywrd(i)) then
          call cssetp(next)
          if(j.eq.1) then
            if(i.eq.2) then
              rght=.false.
            endif
          elseif(j.eq.2) then
            if(i.eq.1) then
              if(rght) then
                other=.false.
              endif
            else
              if(.not.rght) then
                other=.false.
              endif
            endif
          endif
          goto 1
        elseif(j.eq.2) then
          other=.false.
        endif
   10 continue
      if(rght) then
        parity=-1.
      else
        parity=1.
      endif
c ... write title
      write(lfnd,*)'NEWFRAME;SET FONT DUPLEX;SET TITLE SIZE -2'
      write(lfnd,1001) (window(i),i=1,4)
 1001 format(' SET WINDOW X ',F5.2,1X,F5.2,' Y ',F5.2,1X,F5.2)
      write(lfnd,*)'TITLE TOP ''',title(1:min(59,lene(title))),''''
      write(lfnd,*)'CASE      ''',case(1:min(59,lene(case))) ,''''
      call fdate1(dat)
      write(lfnd,*)'SET TITLE SIZE -1.6'
      write(lfnd,*)'TITLE 6.9 9.15 ''',dat,  ''''
      write(lfnd,*)'SET TITLE SIZE -3'
c ... calc frame vector of projection plane
      pn1(1)=cos(alpha)*cos(beta)
      pn1(2)=sin(alpha)*cos(beta)
      pn1(3)=sin(beta)
      pn2(1)=-sin(alpha)*cos(gamma)+cos(alpha)*sin(beta)*sin(gamma)
      pn2(2)=cos(alpha)*cos(gamma)+sin(alpha)*sin(beta)*sin(gamma)
      pn2(3)=cos(beta)*sin(gamma)
      pn3(1)=sin(alpha)*sin(gamma)+cos(alpha)*sin(beta)*cos(gamma)
      pn3(2)=-cos(alpha)*sin(gamma)+sin(alpha)*sin(beta)*cos(gamma)
      pn3(3)=cos(beta)*cos(gamma)
c ... find limits in projection plane
      xmin=1.d74
      xmax=-1.d74
      ymin=1.d74
      ymax=-1.d74
      do 20 i=1,nlat-1
        x=geo(1,4,i)*pn1(1)+ geo(2,4,i)*pn1(2)+ geo(3,4,i)*pn1(3)
        y=geo(1,4,i)*pn2(1)+ geo(2,4,i)*pn2(2)+ geo(3,4,i)*pn2(3)
        xmax=max(x,xmax)
        xmin=min(x,xmin)
        ymax=max(y,ymax)
        ymin=min(y,ymin)
   20 continue
      span=1.10*max(xmax-xmin,ymax-ymin)/2d0
      xc=(xmin+xmax)/2d0
      yc=(ymin+ymax)/2d0
      xmin=xc-span
      xmax=xc+span
      ymin=yc-span
      ymax=yc+span
      if(abs(xmax) .lt. 1.d-50) xmax=0d0
      if(abs(xmin) .lt. 1.d-50) xmin=0d0
      if(abs(ymax) .lt. 1.d-50) ymax=0d0
      if(abs(ymin) .lt. 1.d-50) ymin=0d0
      scale=span*2./(window(2)-window(1))
      write(lfnd,*)
     1  'SET LIMIT X ',sngl(xmin),' ',sngl(xmax)
      yc=(ymin+ymax)*.5d0
      write(lfnd,*)
     2  'SET LIMIT Y ',sngl(yc+(ymin-yc)/ymag),' ',
     $     sngl(yc+(ymax-yc)/ymag)
      write(lfnd,*)'SET OUTLINE ALL OFF'
      write(lfnd,*)'SET LABELS ALL OFF;SET TICKS ALL OFF'
c ... write frame
      call tdinit(lfnd,'JOIN 1','SOLID   ')
      call tdput(xmin,ymin)
      xp=xmin+pn1(1)*span*0.1
      yp=ymin+pn1(2)*span*0.1
      call tdput(xp,yp)
      call tdterm
      call tdinit(lfnd,'JOIN 1','DOTDASH ')
      call tdput(xmin,ymin)
      xp=xmin+pn2(1)*span*0.1
      yp=ymin+pn2(2)*span*0.1
      call tdput(xp,yp)
      call tdterm
      call tdinit(lfnd,'JOIN 1','DOTS    ')
      call tdput(xmin,ymin)
      xp=xmin+pn3(1)*span*0.1
      yp=ymin+pn3(2)*span*0.1
      call tdput(xp,yp)
      call tdterm
c ... write beam line
      call tdinit(lfnd,'JOIN 1','SOLID   ')
      xa=geo(1,4,1)*pn1(1)+ geo(2,4,1)*pn1(2)+ geo(3,4,1)*pn1(3)
      ya=geo(1,4,1)*pn2(1)+ geo(2,4,1)*pn2(2)+ geo(3,4,1)*pn2(3)
      call tdput(xa,ya)
      do 30 i=2,nlat-1
        x=geo(1,4,i)*pn1(1)+ geo(2,4,i)*pn1(2)+ geo(3,4,i)*pn1(3)
        y=geo(1,4,i)*pn2(1)+ geo(2,4,i)*pn2(2)+ geo(3,4,i)*pn2(3)
        if(abs(x-xa)+abs(y-ya).ne.0d0) then
          call tdput(x,y)
          xa=x
          ya=y
        endif
   30 continue
      call tdterm
c ... write element box
      bh=0.015*span
      qh=0.03*span
      sh=0.0225*span
      fh=0.003*span
      call tdinit(lfnd,'JOIN 1','SOLID   ')
      x=geo(1,4,1)*pn1(1)+ geo(2,4,1)*pn1(2)+ geo(3,4,1)*pn1(3)
      y=geo(1,4,1)*pn2(1)+ geo(2,4,1)*pn2(2)+ geo(3,4,1)*pn2(3)
      call tdput(x,y)
      i=1
   41 ds(1)=geo(1,4,i+1)-geo(1,4,i)
      ds(2)=geo(2,4,i+1)-geo(2,4,i)
      ds(3)=geo(3,4,i+1)-geo(3,4,i)
      dx=ds(1)*pn1(1)+ ds(2)*pn1(2)+ ds(3)*pn1(3)
      dy=ds(1)*pn2(1)+ ds(2)*pn2(2)+ ds(3)*pn2(3)
      if(dx*dx+dy*dy.eq.0)then
        i=i+1
        goto 41
      else
        angl0=atan2(dy,dx)
      endif
      do 40 i=1,nlat-1
        x=geo(1,4,i)*pn1(1)+ geo(2,4,i)*pn1(2)+ geo(3,4,i)*pn1(3)
        y=geo(1,4,i)*pn2(1)+ geo(2,4,i)*pn2(2)+ geo(3,4,i)*pn2(3)
        call tdput(x,y)
        ds(1)=geo(1,4,i+1)-geo(1,4,i)
        ds(2)=geo(2,4,i+1)-geo(2,4,i)
        ds(3)=geo(3,4,i+1)-geo(3,4,i)
        dx=ds(1)*pn1(1)+ ds(2)*pn1(2)+ ds(3)*pn1(3)
        dy=ds(1)*pn2(1)+ ds(2)*pn2(2)+ ds(3)*pn2(3)
        dss=abs(dx)+abs(dy)
        if(dss.ne.0) then
          angl=atan2(dy,dx)
        else
          angl=angl0
        endif
        n=idcomp(elatt,i)
        id=idtype(n)
        if(id.eq.1) then
          bw=fh
        elseif(id.eq.2) then
          bw=bh
        elseif(id.eq.4) then
          bw=qh
        elseif(id.le. 32) then
          bw=sh*.5d0
        elseif(id .eq. icMARK)then
          bw=fh
        else
          bw=0
        endif
        if(bw.ne.0.) then
          bx=x  -bw*sin(angl)
          by=y  +bw*cos(angl)
          call tdput(bx,by)
          if(dss.ne.0d0) then
            bx=bx+ dx
            by=by+ dy
            call tdput(bx,by)
          endif
          bx=bx + 2.*bw*sin(angl)
          by=by - 2.*bw*cos(angl)
          call tdput(bx,by)
          if(dss.ne.0d0) then
            bx=bx- dx
            by=by- dy
            call tdput(bx,by)
          endif
          call tdput(x,y)
        endif
        if(dss.ne.0d0) call tdput(x,y)
        angl0=angl
   40 continue
      call tdterm
c ... write element name
      spname=max(bh,qh,sh)+0.03*span
      csize=0.03*span/(.038462+.049444)/scale
      csmin=0.02*span/(.038462+.049444)/scale
      x=geo(1,4,1)*pn1(1)+ geo(2,4,1)*pn1(2)+ geo(3,4,1)*pn1(3)
      y=geo(1,4,1)*pn2(1)+ geo(2,4,1)*pn2(2)+ geo(3,4,1)*pn2(3)
      i=1
   51 ds(1)=geo(1,4,i+1)-geo(1,4,i)
      ds(2)=geo(2,4,i+1)-geo(2,4,i)
      ds(3)=geo(3,4,i+1)-geo(3,4,i)
      dx=ds(1)*pn1(1)+ ds(2)*pn1(2)+ ds(3)*pn1(3)
      dy=ds(1)*pn2(1)+ ds(2)*pn2(2)+ ds(3)*pn2(3)
      if(dx*dx+dy*dy.eq.0)then
        i=i+1
        goto 51
      else
        angl0=atan2(dy,dx)
      endif
      do 50 i=1,nlat-1
        ds(1)=geo(1,4,i+1)-geo(1,4,i)
        ds(2)=geo(2,4,i+1)-geo(2,4,i)
        ds(3)=geo(3,4,i+1)-geo(3,4,i)
        dx=ds(1)*pn1(1)+ ds(2)*pn1(2)+ ds(3)*pn1(3)
        dy=ds(1)*pn2(1)+ ds(2)*pn2(2)+ ds(3)*pn2(3)
        dss=sqrt(dx*dx+dy*dy)
        if(dss.ne.0) then
          angl=atan2(dy,dx)
        else
          angl=angl0
        endif
        n=idcomp(elatt,i)
        id=idtype(n)
        if( tmatch(pname(n),word) ) then
          call cssetp(next)
          if( pname(n).ne.pname(idcomp(elatt,i+1)) ) then
            cs(1)=x+dx/2.-spname*dy/dss*parity
            cs(2)=y+dy/2.+spname*dx/dss*parity
            cs(3)=atan2(dy,dx)+parity*pi/2.
            cs(4)=parity
            cs(5)=lene(pname(n))
            cs(6)=angl-angl0
            cs(7)=csize
            if( i.ne.1.and. over(cs,csa,dcs,scale) ) then
              if(csize-dcs .lt. csmin) then
                cs(4)=-parity
                cs(1)=x+dx/2.+spname*dy/dss*parity
                cs(2)=y+dy/2.-spname*dx/dss*parity
              else
                cs(7)=csize-dcs
              endif
            endif
            xp=cs(1)
            yp=cs(2)
            name=' '
            if(cs(4).eq.-1.) then
              name=pname(n)
              ap=atan2(dy,dx)*180./pi-90.
              if(rght .or. .not.rght.and.other) then
                write(lfnd,1002) xp,yp,cs(7),ap,'''',name(1:8),''''
              endif
            else
              name(9-int(cs(5)):8)=pname(n)
              ap=atan2(dy,dx)*180./pi-90.
              xp=cs(1)-(1.5*8*0.042857*cs(7)*scale )*dy/dss
              yp=cs(2)+(1.5*8*0.042857*cs(7)*scale )*dx/dss
              if(.not.rght .or. rght.and.other) then
                write(lfnd,1002) xp,yp,cs(7),ap,'''',name(1:8),''''
              endif
            endif
 1002       format(' TITLE ',2(F9.3,1X),'DATA SIZE ',F4.1,
     &                                                ' ANGLE ',F6.1,3A)
            csa=cs
          endif
        endif
        x=x+dx
        y=y+dy
        angl0=angl
   50 continue
      return
      end
