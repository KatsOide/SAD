      Subroutine openP$(name,x1,x2,y1,y2)
      real*8 x1,x2,y1,y2
      character*(*) name
      include 'inc/MACFILE.inc'
C
      character*8 cdate1,tim
      character cdate*8,ctime*12
      data ctime /'00:00:00****'/
      call datetime(cdate1,tim)
      cdate=cdate1(7:8)//"/"//cdate1(5:6)//"/"//cdate1(3:4)
      ctime=tim(1:2)//';'//tim(3:4)//':'//tim(5:6)
      write(pltfl,'(''('',T5,A8,T14,A8,T25,''MEMBER NAME '',A8,'//
     &             'T61,''TOPDRAW'',T70,'')'')')cdate,ctime(:8),name
      write(pltfl,*) 'newp;new frame  '
      write(pltfl,*) 'set font duplex'
c
      Entry openw$(x1,x2,y1,y2)
c
      if(x1 .ne. x2)
     &  write(pltfl,*) 'set window x from',min(x1,x2),'to',max(x1,x2)
      if(y1 .ne. y2)
     &  write(pltfl,*) 'set window y from',min(y1,y2),'to',max(y1,y2)
      return
c
      Entry limit$(x1,x2,y1,y2)
c
      if(x1 .ne. x2)
     &  write(pltfl,*) 'set limit x from',min(x1,x2),'to',max(x1,x2)
      if(y1 .ne. y2)
     &  write(pltfl,*) 'set limit y from',min(y1,y2),'to',max(y1,y2)
      return
c
      end
