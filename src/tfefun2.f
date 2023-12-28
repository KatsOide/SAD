c     
      subroutine tfdeffun2
      implicit none
      integer*4 i,itfunaloc,map(32),ieval(32)
c     Initialize map/ieval array
      do i=1,32
         map(i)=0
         ieval(i)=0
      enddo
c
c     itfunaloc(fname, f_id,narg,map,ieval)
c     f_id = function id
c     narg = number of arguments
c     map(i)= 
c     ieval(i)=
c     map(narg+1)=ieval(narg+1)=0
c
      i=itfunaloc('CaOpen1',2001,1,map,ieval,0)
      i=itfunaloc('CaClose1',2002,1,map,ieval,0)
      i=itfunaloc('CaRead1',2003,1,map,ieval,0)
      i=itfunaloc('CaWrite1',2004,2,map,ieval,0)
      i=itfunaloc('CaPendIo',2005,1,map,ieval,0)
      i=itfunaloc('CaChName',2006,1,map,ieval,0)
c
      i=itfunaloc('CaOpen2',2007,1,map,ieval,0)
      i=itfunaloc('CaRead2',2008,1,map,ieval,0)
      i=itfunaloc('CaWrite2',2009,2,map,ieval,0)
      i=itfunaloc('CaConStatus1',2010,2,map,ieval,0)
c      i=itfunaloc('CaStatus2',2011,2,map,ieval,0)

      return
      end

c     evaluate function with id larger than 2000
      subroutine tfefun2(isp1,id0,k,kx,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) kd
      type (sad_descriptor) ,intent(in) ::k
      type (sad_descriptor) ,intent(out) ::kx
      integer*4 ,intent(in)::isp1
      integer*4 ,intent(out)::irtc
      integer*8 kv,kav,kax
      integer*4 isp0,id0,id,narg,i,j
      real*8 wtime,v
      integer*4 istatus,itfmessage
      real*8 nfCaOpen,nfCaClose,nfCaWrite,nfCaPendIO
      integer*4 maxrpt
      parameter(maxrpt=3)
      logical*4 wait(256)
c     
      irtc=-1
      id=id0-2000
      narg=isp-isp1
c
c      write(*,*)'tfefun2: ',id
      go to (2001,2002,2003,2004,2005,2000,
     $     2007,2008,2009,2010),id
c     CaOpen1,CaClose,CaRead1,CaWrite1,CaPendIO,CaName
c     CaOpen2,CaRead2,CaWrite2,CaConStatus1

 2000 continue
      write(*,*)'Wrong implementation of function.  ID = ',id0
      return
c     CaOpen1[ch_name]
 2001 continue
      if(narg /= 1)then
         irtc=itfmessage(9,'General::narg','"1"')
      else
         kx=dfromr(dble(nfCaOpen(k,irtc)))
      endif
      return
c     CaOpen2[ch_name,...]
 2007 continue
      if(narg == 0) then
         irtc=itfmessage(9,'General::narg','"1 or more"')
         return
      endif
      isp0=isp
c 12/9/1998 K. Oide
      wtime=min(4.d0,max(0.5d0,narg*0.2d0))
      isp=isp0+narg
      do j=1,maxrpt
        do i=1,narg
          if(ktfstringq(ktastk(isp1+i))) then
            kav=ktfaddr(ktastk(isp1+i))
            call CaSearch(ilist(1,kav+1),
     $           ilist(1,kav),rtastk(isp0+i),irtc)
            if(irtc /= 0)then
              go to 3007
            endif
          else
            go to 3007
          endif
        end do
        call CaPendIo(wtime, irtc)
        if(irtc == 0)then
          kx%k=ktflist+ktfmakelist(isp0)
          isp=isp0
          return
        endif
        wtime=wtime*1.414
      end do
 3007 irtc=itfmessage(9,'CA::open','""')
      isp=isp0
      return
c
 2002 continue
      if(narg /= 1)then
         irtc=itfmessage(9,'General::narg','"1"')
      else
         kx=dfromr(dble(nfCaClose(dtastk(isp1+1),irtc)))
      endif
      return
c     CaRead1[ch]
 2003 continue
      if (narg /= 1) then
         irtc=itfmessage(9,'General::narg','"1"')
      else
         call nfCaRead(isp1,kx,irtc)
      endif
      return
c     CaRead2[ch,...]
 2008 continue
      if(narg == 0) then
         irtc=itfmessage(9,'General::narg','"1 or more"')
         return
      endif
      isp0=isp
c 
c     6/10/1998 K. Oide
      wtime=min(4.d0,max(0.5d0,narg*0.2d0))
      istatus=0
      irtc=0
      do i=1,narg
        if(ktfrealq(ktastk(isp1+i)))then
          call caConStatus(rtastk(isp1+i),istatus)
          if(istatus /= 0)then
            irtc=-1
          endif
        else
          irtc=-1
        endif
        if (irtc /= 0) then
          irtc=itfmessage(9,'CA::unconnected','"channel in the list"')
          return
        endif
      end do
      irtc=0
c
      kax=ktadalocnull(-1,narg)
      do j=1,maxrpt
        do i=1,narg
          isp=isp+1
          if (ktfrealq(ktastk(i+isp1),v)) then
            kv=ktaaloc(0,4)
            call nfCaReadNoWait(kv,v,wait(i),irtc)
            if(irtc == 0) then
              klist(kax+i)=ktflist+kv
            else
              call tflocal1(kv)
            endif
          endif
        end do
c     wait  for CA Operatio to finish
        call CaPendIo(wtime, irtc)
        if(irtc == 0) go to 3008 
        wtime=wtime*1.414
      end do
c
      print *,"CaPendIO expired"
      irtc=itfmessage(9,'CA::read','""')
      return
c
 3008 continue
      do i = 1,narg
        kd=dlist(kax+i)
        if(kd%k /= ktfoper+mtfnull .and. wait(i))then
          call nfCaReadFinish(kd,irtc)
          if(irtc /= 0) then
            print *," failed to receive value at ",i,"-th element code=",irtc
          end if
        endif
      end do
      kx%k=ktflist+kax
      irtc=0
      return

c     CaWrite1
 2004 continue
      if(narg /= 2)then
         irtc=itfmessage(9,'General::narg','"2"')
      else 
         kx=dfromr(dble(nfCaWrite(isp1,kx,irtc)))
      endif
      return
c     CaWrite2[{ch,val},....]
 2009 continue
      if(narg == 0) then
         irtc=itfmessage(9,'General::narg','"1 or more"')
         return
      endif
      kx%k=ktfoper+mtfnull
c 
      do i=1,narg
         if (ktflistq(ktastk(isp1+1))) then
            call nfCaWriteNoWait(dtastk(isp1+i),irtc)
         endif
      end do
c
      irtc=0
      return
c     CaPendIO
 2005 continue
      if(narg /= 1)then
         irtc=itfmessage(9,'General::narg','"1"')
      else
         v=rtastk(isp1+1)
         nfCaPendIO=-1
         if(ktfnonrealq(k))then
           irtc=itfmessage(9,'General::wrongtype','"Real number"')
           return
         endif
         call CaPendIO(rfromd(k),irtc)
         if (irtc /= 0) then
           irtc=itfmessage(9,'CA::pendio','""')
         end if
      endif
      return
c  CaConStatus
 2010 continue
      if(narg /= 1)then
         irtc=itfmessage(9,'General::narg','"1"')
      else
         call CaConStatus(rtastk(isp1+1),irtc)
         kx=dfromr(dble(irtc))
         irtc=0
      endif
      return
      end
c     
      real*8 function nfCaOpen(kx, irtc)
      use tfstk
      implicit none
      type (sad_descriptor) ,intent(in):: kx
      integer*8 kax
      integer*4 ,intent(out):: irtc
      integer*4 itfmessage

      nfCaOpen=-1
      if(ktfnonstringq(kx))then
         nfCaOpen=0
         irtc=itfmessage(9,'General::wrongtype','"Character-string"')
         return
      endif
      kax=ktfaddr(kx)
      call CaSearch(ilist(1,kax+1),ilist(1,kax),nfCaOpen,irtc)
      if (irtc /= 0) then
         irtc=itfmessage(9,'CA::search','""')
         return
      endif
      call CaPendIo(10.d0, irtc)
      if (irtc /= 0) then
         irtc=itfmessage(9,'CA::Channel','"open time out"')
      endif
      return
      end
c
      real*8 function nfCaClose(kx,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) ,intent(in)::  kx
      integer*4 irtc,itfmessage

      nfCaClose=-1
      if(ktfnonrealq(kx))then
         nfCaClose=0
         irtc=itfmessage(9,'General::wrongtype','"Real number"')
         return
      endif
      call CaClose(rfromd(kx),irtc)
      if (irtc == 0) then
         nfCaClose=0
      else
        irtc=itfmessage(9,'CA::close','""')
        nfCaClose=-1
      endif
      return
      end
c
      real*8 function nfCaPendIo(kx,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) ,intent(in):: kx
      integer*4 irtc,itfmessage

      nfCaPendIO=-1
      if(ktfnonrealq(kx))then
         nfCaPendIO=0
         irtc=itfmessage(9,'General::wrongtype','"Real number"')
         return
      endif
      call CaPendIO(rfromd(kx),irtc)
      if (irtc == 0) then
         nfCaPendIO=0
      else
        irtc=itfmessage(9,'CA::pendio','""')
         nfCaPendIO=-1
      endif
      return
      end
c     
c/* database field types */
c     #define   DBF_STRING	0
c     #define	DBF_INT		1
c     #define	DBF_SHORT	1
c     #define	DBF_FLOAT	2
c     #define	DBF_ENUM	3
c     #define	DBF_CHAR	4
c     #define	DBF_LONG	5
c     #define	DBF_DOUBLE	6
c     #define   DBF_NO_ACCESS	7
c     #define	LAST_TYPE	DBF_DOUBLE
c
      subroutine nfCaRead(isp1,kx,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) ,intent(out):: kx
      integer*8 kax
      integer*4 ,intent(in)::  isp1
      integer*4 ,intent(out):: irtc
      real*8 chid
      integer*4 istatus,etype,ecount,itfmessage

      chid=rtastk(isp1+1)
      call caConStatus(chid,istatus)
      if (istatus /= 0) then 
        irtc=itfmessage(9,'CA::Channel','"Not connected"')
        return
      endif
      call CaGetType(chid, etype, ecount)
      kax=ktadaloc(-1,4)
      irtc=-1
      if (etype == 0 .or. etype == 4 .or. etype == 3)then
c     STRING or CHAR or ENUM
         if(ecount == 1) then
            call CaGetString(chid, kax, irtc)
         else
            call CaGetStringArray(chid,kax,irtc)
         endif
      else
         if(ecount == 1) then
            call CaGetReal(chid, kax, irtc)
         else if (ecount .gt. 1) then
            call CaGetRealArray(chid,kax,irtc)
         endif
      endif
      if(irtc /= 0)then
        irtc=itfmessage(9,'CA::read','""')
      endif
      kx%k=ktflist+kax
      return
      end
c
      subroutine nfCaReadNoWait(ka,chid,wait,irtc)
      use tfstk
      implicit none
      integer*8 ,intent(out):: ka
      integer*4 ,intent(out):: irtc
      real*8 ,intent(in):: chid
      integer*4 etype,ecount,itfmessage
      logical*4 ,intent(out):: wait

      wait=.true.
      call CaGetType(chid, etype, ecount)
      irtc=-1
      if (etype == 0 .or. etype == 4 .or. etype == 3)then
c     STRING or CHAR or ENUM
        if(ecount == 1)then
          call CaGetStringNoWait(chid,ka,irtc)
        else
          call CaGetStringArray(chid,ka,irtc) ! temporarily using wait version
          wait=.false.
        endif
        ilist(2,ka-3)=lnonreallist 
      else
        call CaGetRealNoWait(chid, ka, irtc)
      endif
      if(irtc /= 0)then
        irtc=itfmessage(9,'CA::read','""')
      endif
      return
      end
c
      subroutine nfCaReadFinish(k,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) ,intent(in) ::k
      type (sad_dlist) ,pointer ::kl
      integer*8 ka,kx
      integer*4 ,intent(out)::  irtc
      integer*4 itfmessage
      if(ktfnonlistq(k)) then
         irtc=itfmessage(9,'General::wrongtype','"List"')
         return
      endif
      ka=ktfaddr(k)
      kx=klist(ka+1)
c      if(ktfnonrealq(kx)) then
c         irtc=itfmessage(9,'General::wrongtype','"Real number"')
c         return
c      endif
c      etype=int(klist(ka+4))  ! suspicious
      if(ktfrealq(k) .or. tfreallistq(k)) then
         call CaGetRealFinish(ka,irtc)
         if(irtc /= 0)then
           irtc=itfmessage(9,'CA::read','""')
         endif
      else if(ktfstringq(k)) then
         call CaGetStringFinish(ka,irtc)
         if(irtc /= 0)then
           irtc=itfmessage(9,'CA::read','""')
         endif
      else
         irtc=itfmessage(9,'General::wrongtype','"Real or RealList or String"')
      endif
      return
      end
c     
      real*8 function nfCaWrite(isp1,kx,irtc)
      use tfstk
      implicit none
      type (sad_dlist) ,pointer:: kl
      type (sad_descriptor) ,intent(out):: kx
      type (sad_descriptor) ::k
      integer*8 kax
      integer*4 ,intent(in):: isp1
      integer*4 ,intent(out):: irtc
      
      real*8 val,chid
      integer*4 etype,ecount,nc
      integer*4 itfmessage
      character*256 str,tfgetstr

      nfCaWrite=0.0d0
      if(ktfnonrealq(ktastk(isp1+1)))then
         irtc=itfmessage(9,'General::wrongtype','"Real number"')
         return
      endif
      chid=rtastk(isp1+1)
      irtc=0
      nfCaWrite=0.0d0
      if(ktfrealq(ktastk(isp1+2),val))then
         call CaPutReal(chid, val, irtc)
         kx=dfromr(val)
         nfCaWrite=val
      else if(ktfstringq(ktastk(isp1+2))) then
         str=tfgetstr(dtastk(isp1+2),nc)
         call CaPutString(chid,str(:nc),nc,irtc)
         kx%k=ktfstring + ktsalocb(-1, str, nc)
         nfCaWrite=0.0d0
      else if(ktflistq(dtastk(isp1+2),kl)) then
         call CaGetType(chid,etype,ecount)
         ecount=min(ecount,kl%nl)
         if(ecount .le. 0) then
           irtc=-1
         else
           k=kl%dbody(1)
           if(ktfrealq(k)) then
             call CaPutRealArray(chid, ecount, kax, irtc)
             kx%k=ktflist+kax
           else if(ktfstringq(k)) then
             call CaPutStringArray(chid, ecount, kax, irtc)
             kx%k=ktflist+kax
           else
             irtc=itfmessage(9,'General::wrongtype',
     $            '"Real ot String"')
             return
           endif
         endif
      else
         irtc=itfmessage(9,'General::wrongtype',
     $       '"Real, String, or List of them"')
         return
      endif
      if(irtc /= 0)then
        irtc=itfmessage(9,'CA::write','""')
      endif
      return
      end
c
      subroutine nfCaWriteNoWait(kx,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) ,intent(in):: kx
      type (sad_descriptor) kv,k
      type (sad_dlist) ,pointer ::kl,kl1
      integer*8 kav
      integer*4 ,intent(out):: irtc
      real*8 chid
      integer*4 il,istatus,etype,ecount,nc,itfmessage
      character*256 str,tfgetstr
      call descr_sad(kx,kl)
      il=kl%nl
      if (il .le. 0) then
        irtc=-1
      else if (il == 1) then
        print *,"nfCaWriteNoWait: No value to write"
      else if (il /= 2) then
        print *,"nfCaWriteNoWait: too many arguments"
      endif
      if(ktfnonrealq(kl%dbody(1))) then
        print *,"Invalid channel ID"
        irtc=itfmessage(9,'General::wrongtype','"Real number"')
        return
      endif
      chid=kl%rbody(1)
      kv=kl%dbody(2)
      call caConStatus(chid,istatus)
      if (istatus /= 0) then 
        irtc=itfmessage(9,'CA::Channel','""')
        return
      endif
      irtc=0
      if(ktfrealq(kv)) then
        call CaPutRealNowait(chid, kv, irtc)
      else if(ktfstringq(kv)) then
        str=tfgetstr(kv,nc)
        nc=len_trim(str)
        call CaPutStringNowait(chid,str(:nc),nc,irtc)
      else if(ktflistq(kv,kl1)) then
        kav=ktfaddr(kv)
        call CaGetType(chid,etype,ecount)
        ecount=min(ecount,kl1%nl)
        if(ecount .le. 0) then
          irtc=-1
        else
          k=kl%dbody(1)
          if(ktfrealq(k)) then
            call CaPutRealArrayNowait(chid, ecount, kav, irtc)
          else if(ktfstringq(k)) then
            call CaPutStringArrayNowait(chid, ecount, kav, irtc)
          else
            irtc=itfmessage(9,'General::wrongtype',
     $           '"Real or String"')
            return
          endif
        endif
      else
        irtc=itfmessage(9,'General::wrongtype',
     $       '"Real, String, or List of them"')
        return
      endif
      if(irtc /= 0)then
        irtc=itfmessage(9,'CA::write','""')
      endif
      return
      end
