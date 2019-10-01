      integer*4 function hsrch(id)
      use maccbk
      implicit none
      character*(*) id,idw*(MAXPNAME+1)
      integer*4 hsrchz,slen
c for debug
c      print *,'hsrch>',lenw(id),id(:lenw(id))
c end debug
      if (id(1:1) .eq. '$') then
        hsrch=hsrchz(id)
      else
        slen=min(len_trim(id),MAXPNAME)
        idw(1:1)='$'
        idw(2:)=id(:slen)
        hsrch=hsrchz(idw(:slen+1))
      endif
c for debug
c      print *,'hsrch<',slen,id(:slen),idw(:lenw(idw)),hsrch
c end debug
      return
      end

      integer*4 function hsrchz(token)
      use maccbk
      use mackw
      implicit none
      character*(*) token
      character wtoken*(MAXPNAME)
c
      integer*4 loopc,slen,slenw,iw,isum,i
c
      slen=len_trim(token)
      wtoken=token(:slen)
      if(slen .le. 0) then
         print *,'hsrchz warning: NULL token as input'
         hsrchz=0
         rlist(7)=0.d0
      endif
      isum=0
      do i=1,slen
        iw=ichar(wtoken(i:i))
c         read(wtoken(i:i),'(A1)')iw
         isum=iw+2*abs(isum)
      end do
c      isum=isum*isum
      hsrchz=mod(abs(isum), HTMAX)+1
c       write(*,*)'@hsrchz ',token(:slen),slen," -> ", hsrchz
      do loopc=1,HTMAX
         slenw=lpname(hsrchz)
         if(slenw .le. 0) then
            pname(hsrchz)=wtoken(:slen)
            lpname(hsrchz)=slen
c            print *,"@hsrchz insert ",pname(hsrchz)
c     $           ,hsrchz, slen, slenw,loopc
            return
         endif
         if(slenw .eq. slen)then
           if(pname(hsrchz)(:slenw) .eq. wtoken(:slen)) then
c             print *,'@hsrchz found ',pname(hsrchz)(:slenw),hsrchz
             return
           endif
         endif
c            print *,'@hsrchz next ',pname(hsrchz)(:slenw),hsrchz
         hsrchz=mod(hsrchz,HTMAX)+1
      end do
c     
      continue
      call errmsg('hsrchz',
     &           'No more free space in hash table',
     &           0,16)
      end

      integer function hsrchz1(token)
      use maccbk
      use maccode
      implicit none
      character*(*) token
      character wtoken*(MAXPNAME)
c     
      integer*4 loopc,slen,slenw,iw,isum,i
c     
      slen=len_trim(token)
      wtoken=token(:slen)
      if(slen .le. 0) then
        print *,'hsrchz warning: NULL token as input'
        hsrchz1=0
         rlist(7)=0.d0
      endif
      isum=0
      do i=1,slen
        iw=ichar(wtoken(i:i))
        isum=iw+2*abs(isum)
      end do
c      isum=isum*isum
      hsrchz1=mod(abs(isum), HTMAX)+1
      do loopc=1,HTMAX
        slenw=lpname(hsrchz1)
        if(slenw .le. 0) then
          hsrchz1=0
          return
        endif
        if(slenw .eq. slen)then
          if(pname(hsrchz1)(:slenw) .eq. wtoken(:slen)) then
            return
          endif
        endif
        hsrchz1=mod(hsrchz1,HTMAX)+1
      end do
c     
      continue
      call errmsg('hsrchz1',
     &     'No more free space in hash table',
     &     0,16)
      end
