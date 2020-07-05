      subroutine pstati3(data,x,ex,n,ghisto,title,case,lfnd,lfno)
      use tfstk
      use ffs
      use tffitcode
      include 'inc/TFMACRO.inc'
      parameter (ndi=8)
      logical ghisto
      character*(*) title,case
      character*11 autofg,vout(7,ndi)
      character*80 name,namf
c     character*1 prefix,prefixc
      character*3 sufix
      dimension data(ndi,*),x(ndi),ex(ndi),ceil(ndi),floor(ndi)
      call pclr(x,ndi)
      call pclr(ex,ndi)
      do 10 i=1,ndi
        ceil(i)=-1d31
        floor(i)=1d31
        do 20 k=1,n
          ceil(i)=max(ceil(i),data(i,k))
          floor(i)=min(floor(i),data(i,k))
 20     continue
        do 22 k=1,n
          x(i)=x(i)+data(i,k)
 22     continue
        x(i)=x(i)/n
        do 26 k=1,n
          ex(i)=ex(i)+(data(i,k)-x(i))**2
 26     continue
        ex(i)=sqrt(ex(i)/max(n-1,1))
   10 continue
      do 30 i=1,ndi
        vout(1,i)=autofg(x(i),'11.8')
        vout(2,i)=autofg(ex(i),'11.8')
        vout(3,i)=autofg(floor(i),'11.8')
        vout(4,i)=autofg(ceil(i),'11.8')
 30   continue
      it=italoc(n)
      do i=1,ndi
        do k=1,n
          rlist(it-1+k)=data(i,k)
        enddo
        if(i.eq.6) then
          icl=int(n*0.68268)
          cl=abs(pselect(n-icl,rlist(it),n))
          cl=max(cl,abs(pselect(icl,rlist(it),n)))
          vout(5,i)=autofg(cl,'11.8')
          icl=int(n*0.86638)
          cl=abs(pselect(n-icl,rlist(it),n))
          cl=max(cl,abs(pselect(icl,rlist(it),n)))
          vout(6,i)=autofg(cl,'11.8')
          icl=int(n*0.9545)
          cl=abs(pselect(n-icl,rlist(it),n))
          cl=max(cl,abs(pselect(icl,rlist(it),n)))
          vout(7,i)=autofg(cl,'11.8')
        else
          icl=int(n*0.68268)
          cl=pselect(icl,rlist(it),n)
          vout(5,i)=autofg(cl,'11.8')
          icl=int(n*0.86638)
          cl=pselect(icl,rlist(it),n)
          vout(6,i)=autofg(cl,'11.8')
          icl=int(n*0.9545)
          cl=pselect(icl,rlist(it),n)
          vout(7,i)=autofg(cl,'11.8')
        endif
      enddo
      call tfree(int8(it))
      write(lfno,9102) ((vout(j,i),j=1,7),i=1,ndi)
9102  format(
     $     1x,'Emittance X  ',a,t27,a,t43,a,t59,a,t75,a,t91,a,t107,a/
     $     1x,'Emittance Y  ',a,t27,a,t43,a,t59,a,t75,a,t91,a,t107,a/
     $     1x,'Emittance Z  ',a,t27,a,t43,a,t59,a,t75,a,t91,a,t107,a/
     $     1x,'Energy sprd  ',a,t27,a,t43,a,t59,a,t75,a,t91,a,t107,a/
     $     1x,'Bunch Lngth  ',a,t27,a,t43,a,t59,a,t75,a,t91,a,t107,a/
     $     1x,'Beam tilt    ',a,t27,a,t43,a,t59,a,t75,a,t91,a,t107,a/
     $     1X,'Beam size x  ',a,t27,a,t43,a,t59,a,t75,a,t91,a,t107,a/
     $     1X,'Beam size y  ',a,t27,a,t43,a,t59,a,t75,a,t91,a,t107,a)
c----- HistoPlot -------------
      if(ghisto .and. n.gt.1) then
        nbin=10*max(log10(dble(n)),1d0)
        ig=italoc(nbin)
        call pstatihist(rlist(ig:ig+nbin-1),xmin,xmax,nbin,title,case,name,namf,
     $       0,0,1,lfnd)
        do i=1,ndi
          do j=1,nbin
            rlist(ig-1+j)=0d0
          enddo
c         xmin=min(floor(i),0d0)
          xmin=floor(i)
c         xmax=max(0d0,ceil(i))
          xmax=ceil(i)
          wbin=(xmax-xmin)/dble(nbin)
          do k=1,n
            j=min(1+int((data(i,k)-xmin)/wbin),nbin)
            rlist(ig-1+j)=rlist(ig-1+j)+1d0
          enddo
          iscale=int(log10(max(ceil(i),-floor(i))))
c         iscale=(iscale/3+isign(1,iscale))*3
          if(iscale.lt.0) iscale=iscale-1
          xmin=xmin*10.**(-iscale)
          xmax=xmax*10.**(-iscale)
c         prefixc=' '
c         if(iscale.eq.0) then
c           prefix =' '
c         elseif(iscale.eq.3) then
c           prefix ='K'
c           prefixc='L'
c         elseif(iscale.eq.6) then
c           prefix ='M'
c         elseif(iscale.eq.9) then
c           prefix ='G'
c         elseif(iscale.eq.12) then
c           prefix ='T'
c         elseif(iscale.eq.-3) then
c           prefix ='M'
c           prefixc='L'
c         elseif(iscale.eq.-6) then
c           prefix ='M'
c           prefixc='G'
c         elseif(iscale.eq.-9) then
c           prefix ='N'
c           prefixc='L'
c         elseif(iscale.eq.-12) then
c           prefix ='P'
c           prefixc='L'
c         elseif(iscale.eq.-15) then
c           prefix ='F'
c           prefixc='L'
c         elseif(iscale.eq.-18) then
c           prefix ='A'
c           prefixc='L'
c         endif
          name=' '
          namf=' '
          if(iscale.ne.0) then
            if(iscale.le.-10) then
              write(sufix,'(i3)')iscale
            elseif(iscale.lt.0 .or. iscale.ge.10) then
              write(sufix,'(i2)')iscale
            else
              write(sufix,'(i1)')iscale
            endif
            name=' (102'//sufix(1:lene(sufix))//'3'
            namf=' '
            namf='    X'//namf(1:lene(sufix))//'X'
          endif
          if(i.eq.1) then
c           name='E0X1 ('//prefix//'M)'
c           namf='GXLX  '//prefixc//'L '
            name='E0X1'//name(1:lene(name))//'M)'
            namf='GXLX'//namf(1:lene(namf))//'L'
          elseif(i.eq.2) then
c           name='E0Y1 ('//prefix//'M)'
c           namf='GXLX  '//prefixc//'L '
            name='E0Y1'//name(1:lene(name))//'M)'
            namf='GXLX'//namf(1:lene(namf))//'L'
          elseif(i.eq.3) then
c           name='E0Z1 ('//prefix//'M)'
c           namf='GXLX  '//prefixc//'L '
            name='E0Z1'//name(1:lene(name))//'M)'
            namf='GXLX'//namf(1:lene(namf))//'L'
          elseif(i.eq.4) then
            name='DE/E'//name(1:lene(name))//')'
            namf='F   '//namf(1:lene(namf))
          elseif(i.eq.5) then
c           name='S0Z1 ('//prefix//'M)'
c           namf='GXLX  '//prefixc//'L '
            name='S0Z1'//name(1:lene(name))//'M)'
            namf='GXLX'//namf(1:lene(namf))//'L'
          elseif(i.eq.6) then
c           name='Q ('//prefix//'RAD)'
c           namf='G  '//prefixc//'LLL '
            name='Q'//name(1:lene(name))//'RAD)'
            namf='G'//namf(1:lene(namf))//'LLL'
          elseif(i.eq.7) then
c           name='S0X1 ('//prefix//'M)'
c           namf='GXLX  '//prefixc//'L '
            name='S0X1'//name(1:lene(name))//'M)'
            namf='GXLX'//namf(1:lene(namf))//'L'
          elseif(i.eq.8) then
c           name='S0Y1 ('//prefix//'M)'
c           namf='GXLX  '//prefixc//'L '
            name='S0Y1'//name(1:lene(name))//'M)'
            namf='GXLX'//namf(1:lene(namf))//'L'
          endif
          call pstatihist(rlist(ig:ig+nbin-1),xmin,xmax,nbin,
     $         title,case,name,namf,
     $         3,3,i,lfnd)
        enddo
        call tfree(int8(ig))
      endif
      return
      end
c
      subroutine pstatihist(data,xmin,xmax,nbin,title,case,label,labelc,
     $     iwh,iwv,iw,lfnd)
      implicit real*8 (a-h,o-z)
      parameter (xlmargin=.75d0,ydmargin=1.5d0,xofset=.75d0,
     $     yofset=0.5d0,wh0=13.d0-xofset,wv0=8.3d0-yofset)
      character*(*) title,case,label,labelc
      character*30 dat
      dimension data(nbin)
c
      if(iwh.eq.0 .or. iwv.eq.0) then
        call fdate1(dat)
        write(lfnd,*)'NEWFRAME;SET FONT DUPLEX;SET TITLE SIZE -2.5'
        write(lfnd,*)'SET WINDOW X 3.3 12.1 Y 2.5  8.35'
        write(lfnd,*)
     1              'TITLE 8.5 9.1 CENTER '' '''
        write(lfnd,*)'MORE ''',title(1:min(59,lene(title))),''''
        write(lfnd,*)'CASE ''',case(1:min(59,lene(case))) ,''''
        write(lfnd,*)'SET TITLE SIZE -1.6'
cslac   write(lfnd,*)'TITLE 7 8.5 ''',dat,''''
        write(lfnd,*)'TITLE SIZE 1 10.4 8.5 ''',dat,  ''''
        return
      endif
      wh=wh0/dble(iwh)
      wv=wv0/dble(iwv)
      ih=min(mod((iw-1),iwh)+1,iwh)
      iv=min((iw-1)/iwh+1,iwv)
c     who=(ih-1)*wh+0.5d0*wh
      who=(ih-1)*wh+0.5d0*wh + xofset
      wvo=(iwv-iv)*wv+0.5d0*wv + yofset
      wx1=who-0.5d0*min(wh,wv)
      wx2=who+0.5d0*min(wh,wv)
      wy1=wvo-0.5d0*min(wh,wv)
      wy2=wvo+0.5d0*min(wh,wv)
      slmargin=xlmargin/sqrt(dble(max(iwh,iwv)))
      sdmargin=ydmargin/sqrt(dble(max(iwh,iwv)))
c
      write(lfnd,*)'SET WINDOW X ',sngl(wx1+slmargin),sngl(wx2)
      write(lfnd,*)'SET WINDOW Y ',sngl(wy1+sdmargin),sngl(wy2)
      write(lfnd,*)'SET LABELS SIZE=',
     $     sngl(-2d0/sqrt(dble(max(iwh,iwv))))
      write(lfnd,*)'SET TICKS SIZE=',
     $     sngl(-.1d0/sqrt(dble(max(iwh,iwv))))
      write(lfnd,*)'TITLE BOTTOM ''',label(1:lene(label)),''''
      write(lfnd,*)'CASE         ''',labelc(1:lene(labelc)),''''
      write(lfnd,*)'TITLE LEFT ''NEVENTS/BIN'''
      write(lfnd,*)'CASE       '' LLLLLL LLL'''
      if(xmin.eq.0d0) then
        write(lfnd,*)'SET LIMIT X ',sngl(xmin),' TO ',
     $       sngl(1.1*xmax-0.1*xmin)
      elseif(xmax.eq.0d0) then
        write(lfnd,*)'SET LIMIT X ',sngl(1.1*xmin-0.1*xmax),' TO ',
     $       sngl(xmax)
      else
        write(lfnd,*)'SET LIMIT X ',sngl(1.1*xmin-0.1*xmax),' TO ',
     $       sngl(1.1*xmax-0.1*xmin)
      endif
      call tdinit(lfnd,'JOIN 1','10')
      call tdput(xmin,0.d0)
      bin=(xmax-xmin)/dble(nbin)
      do i=1,nbin
        call tdput(xmin+bin*(i-1),data(i))
        call tdput(xmin+bin*i,data(i))
      enddo
      call tdput(xmax,0.d0)
      call tdterm
c
      return
      end
