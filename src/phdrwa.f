      subroutine phdrwa(ep,zp,npart,ax,title,case,ip)
      implicit none
      include 'inc/MACMATH.inc'
      integer nbin,npara
      parameter (nbin=100,npara=3)

      integer*4 npart,ip
      real*8 ep(npart),zp(npart)
      integer*4 ee(nbin),zz(nbin),ax(2)
      character*(*) title,case

      character*30 dat
      character name(6,2)*4,tag(6)*2
      real*8 x,ydata,p,fwk,wei,covar,sige0,emean0,sigz0,zmean0,
     $     xlimit,ylimit,xmin,ymin,xmax,ymax,
     $     dbin,zmean1,sigz1,sum,chisq,yy,emean1,sige1
      integer*4 lisp,i,nexcl,lenw,k,l,ibin,nex,iex1,iex2
      dimension x(nbin),ydata(nbin),p(npara),fwk(npara,npara),wei(nbin),
     1                 lisp(npara),covar(npara,npara)

      external gaus3

      data name /'X   ','P0X1','Y   ','P0Y1','Z   ','P0Z1',
     1           'L   ','LXLX','L   ','LXLX','L   ','LXLX'/
      data tag /'x ','px','y ','py','z ','pz'/

      sige0=0.
      emean0=0.
      sigz0=0.
      zmean0=0.
      do 5 i=1,npart
        emean0=emean0+ep(i)
        sige0=sige0+ep(i)*ep(i)
        zmean0=zmean0+zp(i)
    5   sigz0=sigz0+zp(i)*zp(i)
      emean0=emean0/npart
      zmean0=zmean0/npart
      sige0=sqrt(sige0/npart-emean0**2)
      sigz0=sqrt(sigz0/npart-zmean0**2)
      nexcl=2*((npart*0.05)/2)
      xlimit=4.*sigz0
      ylimit=4.*sige0
      xmin=zmean0-xlimit
      xmax=zmean0+xlimit
      ymin=emean0-ylimit
      ymax=emean0+ylimit
      write(ip,*)'NEWFRAME;SET FONT DUPLEX;SET TITLE SIZE -2'
      write(ip,*) ' SET WINDOW X 2.2 TO 9.2 Y 2 TO 9'
      write(ip,*)'TITLE TOP ''',title(1:min(59,lenw(title))),''''
      write(ip,*)'CASE      ''',case(1:min(59,lenw(case))) ,''''
      call fdate1(dat)
      write(ip,*)'SET TITLE SIZE -1.3'
      write(ip,*)'TITLE 7.35 9.15 ''',dat,  ''''
      write(ip,*)'SET TITLE SIZE -3'
      write(ip,*) ' SET LIMIT X ',sngl(xmin),' TO ',sngl(xmax)
      write(ip,*) ' SET LIMIT Y ',sngl(ymin),' TO ',sngl(ymax)
      write(ip,*) ' SET LABELS RIGHT OFF BOTTOM OFF'
      k=ax(2)
      l=4-3*mod(k,2)
      write(ip,*) ' TITLE LEFT ''',name(k,1)(:l),''''
      write(ip,*) ' CASE       ''',name(k,2)(:l),''''
      write(ip,*) ' SET SYMBOL .M'
      call tdinit(ip,'PLOT',' ')
      do 10 i=1,npart
        call tdput(zp(i),ep(i))
   10 continue
      call tdterm
c ---  projection plot
      do 110 i=1,nbin
        ee(i)=0
  110   zz(i)=0
      do 100 i=1,npart
        if(ep(i).gt.ymin.and.ep(i).lt.ymax) then
          ibin=int( (ep(i)-ymin)/(2.d0*ylimit)*nbin )+1
          ee(ibin)=ee(ibin)+1
        endif
        if(zp(i).gt.xmin.and.zp(i).lt.xmax) then
          ibin=int( (zp(i)-xmin)/(2.d0*xlimit)*nbin )+1
          zz(ibin)=zz(ibin)+1
        endif
  100 continue
c     write(*,'(10I5)') (zz(i),i=1,nbin)
c  --- plot z_projection
      write(ip,*) ' SET WINDOW X 2.2 TO 9.2 Y 1 TO 2'
      write(ip,*) ' SET LABELS OFF; SET TICKS ALL OFF'
      write(ip,*) ' SET LIMIT X ',sngl(xmin),' TO ',sngl(xmax)
      k=ax(1)
      l=4-3*mod(k,2)
      write(ip,*) ' TITLE BOTTOM ''',name(k,1)(:l),''''
      write(ip,*) ' CASE         ''',name(k,2)(:l),''''
      dbin=2.*xlimit/nbin
      call tdinit(ip,'JOIN 1','10')
      do 120 i=1,nbin
        call tdput(xmin+(dble(i)-1.d0)*dbin,dble(zz(i)))
        call tdput(xmin+ dble(   i)   *dbin,dble(zz(i)))
  120 continue
      call tdput(xmax,0.d0)
      call tdterm
c
c  ... data fit to gaussian
c  ..... fit z_data
      nex=0
      iex1=0
      iex2=0
      do 140 i=1,nbin
        nex=nex+zz(i)
        if(nex.ge.nexcl/2) then
          iex1=i
          goto 141
        endif
  140 continue
  141 do 142 i=nbin,1,-1
        nex=nex+zz(i)
        if(nex.ge.nexcl) then
          iex2=i
          goto 143
        endif
  142 continue
  143 zmean1=0.
      sigz1=0.d0
      sum=0.d0
      do 144 i=iex1+1 ,iex2-1
        sum=sum+zz(i)
        sigz1=sigz1+zz(i)*i*i
  144   zmean1=zmean1+zz(i)*i
      zmean1=zmean1/(npart-nex)
      sigz1=sqrt(sigz1/(npart-nex)-zmean1**2)
      write(6,*) ' Statistics of projected distributions'
      write(6,*) '  rms(',tag(ax(1)),')',sngl(sigz1*dbin),
     1           '    (',1.e2*real(npart-nex)/real(npart),
     1                            ' % particles )'
      p(1)=sum/sqrt(pi2)/sigz1
      p(2)=sigz1
      p(3)=zmean1
      do 145 i=1,nbin
c        if(zz(i).ne.0) then
c          wei(i)=1./zz(i)
c        else
          wei(i)=1.
c        endif
        x(i)=i
  145   ydata(i)=zz(i)
      do 146 i=1,npara
        lisp(i)=i
  146 continue
      call nlfit(x,ydata,wei,nbin,p,npara,lisp,covar,fwk,chisq,gaus3)
      write(6,*) '  rms(',tag(ax(1)),')',sngl(p(2)*dbin),
     1           '    ( Gaussian fit using all particles )'
      call tdinit(ip,'JOIN 1','10')
      do 150 i=1,nbin
        call gaus3(x(i),p,yy,fwk,npara)
        call tdput(xmin+(dble(i)-.5d0)*dbin,yy)
  150 continue
      call tdterm
      write(ip,*) ' SET WINDOW X 2.2 TO 9.2 Y 1 TO 2'
      write(ip,*) ' SET LIMIT X ',sngl(xmin),' TO ',sngl(xmax)
      write(ip,*) ' SET TICKS OFF BOTTOM ON;SET LABELS OFF BOTTOM ON;'
      write(ip,*) ' PLOT AXIS'
c  ... write sigz,sige, etc., on screen
      k=ax(1)
      l=4-3*mod(k,2)
      write(ip,2100) name(k,1)(:l),name(k,2)(:l),sigz0,
     z               name(k,1)(:l),name(k,2)(:l),p(2)*dbin
 2100 format(     ' TITLE ''rms(',a,')'' 11 8 SIZE 1.8'/
     1           '  CASE  ''    ',a,''''/
     1           '  TITLE ''   ',1PE9.2,''' 11 7.6 SIZE 1.8'/
     1           '  TITLE ''S(',a,')0Gauss1 '' 11 7.2 SIZE 1.8'/
     1           '  CASE  ''G ',a,' X     X'''/
     1           '  TITLE ''   ',1PE9.2,''' 11 6.8 SIZE 1.8')
c  --- plot e_projection
      write(ip,*) ' SET WINDOW X 9.2 TO 10.2 Y 2 TO 9'
      write(ip,*) ' SET LABELS OFF; SET TICKS ALL OFF'
      write(ip,*) ' SET LIMIT Y ',sngl(ymin),' TO ',sngl(ymax)
      dbin=2.*ylimit/nbin
      call tdinit(ip,'JOIN 1','10')
      do 130 i=1,nbin
        call tdput(dble(-ee(i)),ymin+(dble(i)-1.d0)*dbin)
        call tdput(dble(-ee(i)),ymin+dble(  i)*dbin)
  130 continue
      call tdput(0d0,ymax)
      call tdterm
c  ..... fit e_data
      nex=0
      iex1=0
      iex2=0
      do 160 i=1,nbin
        nex=nex+ee(i)
        if(nex.ge.nexcl/2) then
          iex1=i
          goto 161
        endif
  160 continue
  161 do 162 i=nbin,1,-1
        nex=nex+ee(i)
        if(nex.ge.nexcl) then
          iex2=i
          goto 163
        endif
  162 continue
  163 emean1=0.d0
      sige1=0.d0
      sum=0.d0
      do 164 i=iex1+1 ,iex2-1
        sum=sum+ee(i)
        sige1=sige1+ee(i)*i*i
 164    emean1=emean1+ee(i)*i
      emean1=emean1/(npart-nex)
      sige1=sqrt(sige1/(npart-nex)-emean1**2)
      write(6,*) '  rms(',tag(ax(2)),')',sngl(sige1*dbin),
     1           '    (',1.e2*real(npart-nex)/real(npart),
     1                            ' % particles )'
      p(1)=sum/sqrt(pi2)/sige1
      p(2)=sige1
      p(3)=emean1
      do 165 i=1,nbin
c        if(ee(i).ne.0) then
c          wei(i)=1./ee(i)
c        else
          wei(i)=1.
c        endif
  165   ydata(i)=ee(i)
      call nlfit(x,ydata,wei,nbin,p,npara,lisp,covar,fwk,chisq,gaus3)
      write(6,*) '  rms(',tag(ax(2)),')',sngl(p(2)*dbin),
     1           '    ( Gaussian fit using all particles )'
      call tdinit(ip,'JOIN 1','10')
      do 170 i=1,nbin
        call gaus3(x(i),p,yy,fwk,npara)
        call tdput(-yy,ymin+(dble(i)-.5d0)*dbin)
  170 continue
      call tdterm
      write(ip,*) ' SET WINDOW X 9.2 TO 10.2 Y 2 TO 9'
      write(ip,*) ' SET LIMIT Y ',sngl(ymin),' TO ',sngl(ymax)
      write(ip,*) ' SET TICKS OFF; SET LABELS OFF'
      write(ip,*) ' PLOT AXIS'
c  ... write sigz,sige, etc., on screen
      k=ax(2)
      l=4-3*mod(k,2)
      write(ip,2110) name(k,1)(:l),name(k,2)(:l),sige0,
     z               name(k,1)(:l),name(k,2)(:l),p(2)*dbin
 2110 format(     ' TITLE ''rms(',a,')'' 11 6 SIZE 1.8'/
     1           '  CASE  ''    ',a,''''/
     1           '  TITLE ''   ',1PE9.2,''' 11 5.6 SIZE 1.8'/
     1           '  TITLE ''S(',a,')0Gauss1 '' 11 5.2 SIZE 1.8'/
     1           '  CASE  ''G ',a,' X     X'''/
     1           '  TITLE ''   ',1PE9.2,''' 11 4.8 SIZE 1.8')
      return
      end




