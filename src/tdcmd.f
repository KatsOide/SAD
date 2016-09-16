      subroutine tdcmd$(cmd,lvar,bottom,left,top,right,limit)
      use maccbk
      use macfile
      implicit real*8 (a-h,o-z)
      integer lvar
      character*(*) cmd,bottom,left,top,right
      logical limit
c
      integer*4 lene
      character*8 cdate1
      character*6 time1
      character cdate*9,ctime*12,fmt*80
      data ctime/'00:00:00****'/
      real*8 ymax,ymin
c
      call datetime(cdate1,time1)
      cdate=cdate1(7:8)//"/"//cdate1(5:6)//"/"//cdate1(3:4)
      ctime=time1(1:2)//':'//time1(3:4)//':'//time1(5:6)
c.....for debug
c     call ptrace('tdcmd',1)
c.....end debug
      if (lvar .gt. 0) then
      len=ilist(1,lvar)-1
      call openP$(cmd,3.d0,11.3d0,0.7d0,9d0)
      fmt=top(:lene(top))//'  '//cdate//' '//ctime(:8)
      call title$(fmt,'TOP')
      call title$(bottom,'BOTTOM')
      call title$(left,'LEFT')
      call title$(right,'RIGHT')
c.....Top drawer cannot hadle the number smaller than 1e-38
      ymax=rlist(lvar+2)
      ymin=ymax
      do 500 i=lvar+1,lvar+len-1,2
        if(abs(rlist(i)) .le. 1.d-38) rlist(i)=0.0d0
        if(abs(rlist(i)) .gt. 1.d+37) rlist(i)=1.d+37
        if(abs(rlist(i+1)) .le. 1.d-38) rlist(i+1)=0.0d0
        if(abs(rlist(i+1)) .gt. 1.d+37) rlist(i+1)=1.d+37
        ymax=max(ymax,rlist(i+1))
        ymin=min(ymin,rlist(i+1))
  500 continue
      if( limit) then
        call limit$(rlist(lvar+1),rlist(lvar+len-1),ymin,ymax)
      endif
      do 1000 i=1,len-mod(len,6),6
        write(pltfl,'(1H ,3(2(1PG11.4),'';''))')
     &              (rlist(lvar+i+J),J=0,5)
        if((i .gt. 1) .and. (mod(i-1,1008) .eq. 0)) then
          write(pltfl,*)cmd
c         write(errfl,*)i,cmd
          write(pltfl,'(1H ,(2(1PG11.4),'';''))')
     &                (rlist(lvar+i+J),J=4,5)
        endif
 1000 continue
      do 1100 i=mod(len,6)-1,0,-2
        j=len-i
        write(pltfl,'(1H ,(2(1PG11.4),'';''))')
     &               rlist(lvar+j),rlist(lvar+j+1)
 1100 continue
      endif
      write(pltfl,*)cmd
c.....for debug
c     call ptrace('tdcmd',-1)
c.....end debug
      return
      end
