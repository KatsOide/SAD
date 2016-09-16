      Subroutine tdjin$(lvar,span,center,bottom,left,top,right)
      use maccbk
      use macfile
      implicit real*8 (a-h,o-z)
      integer lvar
      real*8 span,center,t(3),dt
      character*(*) bottom,left,top,right
      character cdate*9,cdate1*8, ctime*12,fmt*80
      DATA CTIME /'00:00:00****'/
      character*6 time1
      integer*4 lene
      CHARACTER*1 SYMB(10)
      DIMENSION PLIST(10),PT(10)
      DATA NPEAK/10/
      DATA NCHECK/5/
      DATA SYMB/'1','2','3','4','5','6','7','8','9','A'/
c
      call datetime(cdate1,time1)
      cdate=cdate1(7:8)//"/"//cdate1(5:6)//"/"//cdate1(3:4)
      ctime=time1(1:2)//':'//time1(3:4)//':'//time1(5:6)
c.....for debug
c     call ptrace('tdjin',1)
c.....end debug
      len=ilist(1,lvar)-1
      t(1)=center-0.5d0*span
      dt=span/dble(len)
      T0=T(1)
c.....for debug
c     print *,len,center,span
c.....end debug
c
      call openP$('$tdjoin',1.8d0,10.1d0,0.7d0,9d0)
      write(pltfl,*) ' set scale y logarithmic'
      write(pltfl,'(A20,T25,1PG11.4,A4,1PG11.4)')
     &            ' set limits x from ',t(1),' to ',t(1)+span
      fmt=top(:lene(top))//'  '//cdate//'/'//ctime(:8)
      call title$(fmt,'TOP')
      call title$(bottom,'BOTTOM')
      call title$(left,'LEFT')
      call title$(right,'RIGHT')
      do 500 i=lvar+1,lvar+len
        if(abs(rlist(i)) .le. 1.d-38) rlist(i)=0.0d0
  500 continue
      do 1000 i=0,len-1-mod(len,3),3
        t(2)=t(1)+dt
        t(3)=t(2)+dt
        write(pltfl,'(1H ,3(2(1PG11.4),'';''))')
     &                   (t(j),rlist(lvar+i+j),j=1,3)
        t(1)=t(3)+dt
        if((i .gt. 1) .and. (mod(i-1,1008) .eq. 0)) then
          write(pltfl,*)'join 1'
        endif
 1000 continue
      do 1100 i=mod(len,3)-1,0,-1
        write(pltfl,'(1H ,(2(1PG11.4),'';''))')t(1),rlist(lvar+len-i)
        t(1)=t(1)+dt
 1100 continue
      write(pltfl,*)'JOIN 1'
C
C  SEARCH FOR PEAKS OF A SPECTRUM
C
      DO 2000 I=1,NPEAK
         PLIST(I)=0.
         PT(I)   =0.
 2000 CONTINUE
C
      DO 2100 I=LVAR+1,LVAR+LEN-1
         IF(RLIST(I).LE.RLIST(I+1)) GOTO 2100
         PDATA=RLIST(I)
C        WRITE(6,*) PDATA,I,'CAN 1'
         DO 2010 II=1,NCHECK
            I1=MOD(I+II-LVAR-1,LEN)+LVAR+1
            I2=MOD(I-II-LVAR-1,LEN)+LVAR+1
            IF(RLIST(I+II).GE.PDATA) GOTO 2100
            IF(RLIST(I-II).GE.PDATA) GOTO 2100
 2010    CONTINUE
C        WRITE(6,*) PDATA,I,'CAN 2'
         DO 2020 II=1,NPEAK
            IF(PLIST(II).GE.PDATA) GOTO 2020
            DO 2015 III=NPEAK-1,II,-1
               PLIST(III+1)=PLIST(III)
               PT(III+1)   =PT(III)
 2015       CONTINUE
C           WRITE(6,*) PDATA,I,II,'SET'
            PLIST(II)=PDATA
            PT(II)   =(I-LVAR-1)*DT+T0
            GOTO 2100
 2020    CONTINUE
C
 2100 CONTINUE
C
      DO 2200 I=1,NPEAK
         PLIST0=PLIST(I)*1.5
         WRITE(PLTFL,*) ' SET SYMBOL ',SYMB(I),';',PT(I),PLIST0
 2200 CONTINUE
      WRITE(PLTFL,*) ' PLOT; SET TITLE SIZE 1.8'
      WRITE(PLTFL,*) ' TITLE 10.6 8 " TUNE  AMPLITUDE"'
      DO 2210 I=1,NPEAK
         WRITE(PLTFL,600) 8-0.5*I,I,PT(I),PLIST(I)
  600    FORMAT(1H ,'TITLE 10.6',F5.1,'"',I2,1X,F5.4,1PE9.2,'"')
 2210 CONTINUE
C
c.....for debug
c     call ptrace('tdjin',-1)
c.....end debug
      return
      end
