      subroutine pmbdrw(datas,iobs,nobs,ndata,title,case,lfno)
c----- Draw accumulated data of PMEAS. -------      
      implicit real*8 (a-h,o-z)
      parameter (item=20,itemn=item/2)
      character*(*) title,case
c     character line*131,lepton(2)*4,autofg*8,name*8
      dimension iobs(nobs+1),datas(itemn,2,nobs,ndata)
c     data lepton/'(e+)','(e-)'/
c
c..... plot min and max of emiy/emix ....     
      write(lfno,*)'NEWFRAME;SET FONT DUPLEX;SET TITLE SIZE -3'
      write(lfno,*)'SET WINDOW X 3.3 12.1 Y 2.5  8.35'
      write(lfno,*)'TITLE 8.5 9.1 CENTER '' '''
      write(lfno,*)'MORE  ''',title(1:min(59,lene(title))),''''
      write(lfno,*)'CASE '' ',case(1:min(59,lene(case))) ,''''
      write(lfno,*)'SET TITLE SIZE -2.6'
      write(lfno,*)'TITLE BOTTOM ''NUMBER OF STEP'''
      write(lfno,*)'CASE         '' LLLLL LL  LLL'''
      write(lfno,*)'TITLE LEFT ''E0Y1/E0X1'' SIZE -3'
      write(lfno,*)'CASE       ''GXLX GXLX'''
      write(lfno,*)'SET SYMBOL ''. '' SIZE 0.02'
      call tdinit(lfno,'PLOT;JOIN 1','10 ')
      ymax=abs(datas(10,1,1,1))
      ymin=ymax
      do 8 i=1,ndata
        y=abs(datas(10,1,1,i))
        if(y.lt.ymin) ymin=y
        ymax=max(y,datas(10,2,1,i),ymax)
 8    continue
      write(lfno,'(a,i4,a,1p,2g13.5)')'SET LIMIT X 0 ',ndata,' Y ',
     $     1.05*ymin-0.05*ymax,1.05*ymax-0.05*ymin
      yp=abs(datas(10,1,1,1))
      do 10 i=1,ndata
        y=datas(10,1,1,i)
        if(y.lt.0d0) then 
          yp=abs(y)
        endif
        call tdput(dble(i),yp)
 10   continue
      call tdterm
      yp=datas(10,2,1,1)
      call tdinit(lfno,'PLOT;JOIN 1','.2 .04 ')
      do 12 i=1,ndata
        y=datas(10,2,1,i)
        if(y.ne.0d0) yp=y
        call tdput(dble(i),yp)
 12   continue
      call tdterm
      return
      end
