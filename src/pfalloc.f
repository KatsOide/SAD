c     SAD internal memory chunk allocator
c
c     Structure of Memory Chunk List
c
c  /- ilist(2, ilistroot or ilistroot+1)
c  |
c  |    ____________      _________            ________
c  \->| size | ptr |--->|size| ptr|--> ....->|size| 0 |
c     |~~~~~~~~~~~~|    |         |          |    |   |
c     |~~~~~~~~~~~~|    |         |           ~~~~~~~~
c     |~~~~~~~~~~~~|     ~~~~~~~~~
c      ~~~~~~~~~~~~
c

c     Public API
c     Initialize memory block list root node
      subroutine inimem
      use tfstk
      use iso_c_binding
      use macfile
      use tfmem
      implicit none
c      integer*4 memsize

c      call talocinit
c      call lminit(klist(0), 8)

c      if(MAXMEM0 .gt. 2) then
c        memsize = MAXMEM0
c        ilistroot = RBASE
c      else
c        memsize = MAXMEM
c        ilistroot = lmalloc(memsize,irtc)
c        if(irtc .ne. 0) then
c          stop 'inimem is failed to allocate root memory node'
c        endif
c      endif
c$$$      memsize=4
c$$$c      memsize=MAXMEM0
c$$$      ilistroot=ktaloc(memsize)
c$$$
c$$$      ilist(1,ilistroot)=0
c$$$      ilist(2,ilistroot)=ilistroot+2
c$$$      if(memsize .gt. inipage*pagesz) then
c$$$         ilist(1,ilistroot+2)=inipage*pagesz-2
c$$$         ilist(2,ilistroot+2)=0
c$$$         ilist(1,ilistroot+1)=0
c$$$         ilist(2,ilistroot+1)=ilistroot+inipage*pagesz
c$$$         ilist(1,ilist(2,ilistroot+1))=(memsize-inipage*pagesz)
c$$$         ilist(2,ilist(2,ilistroot+1))=0
c$$$      else
c$$$         ilist(1,ilistroot+2)=memsize-2
c$$$         ilist(2,ilistroot+2)=0
c$$$         ilist(1,ilistroot+1)=0
c$$$         ilist(2,ilistroot+1)=0
c$$$      endif
c      write(*,*)'@ inimem'
c      call pr_mem_map()
      return
      end

c     Allocate small memory block initialized by Real[0]
      integer*4 function mcfallo(n)
      use maccbk
      use tfstk
      use macfile
      implicit none
      integer n
c      integer*4 i,mfalloc
      
      write(*,*)'mcfallo is obsolete: use mctaloc. ',n
      call abort
c      mcfallo=mfalloc(n)
c      if (n .ge. 1) then
c         do i=0,n-1
c            rlist(mcfallo+i)=0.0d0
c         enddo
c      endif
c      return
      end

c     Allocate/Free small memory block
      integer*4 function mfalloc(n)
      use maccbk
      use tfstk
      use macfile
      implicit none
      integer*4 n
c      integer*4 pfalloc
c      character*(MAXLLEN) msgstr
c      integer*4 currp,next,msize
c      integer*4 lene

      write(*,*)'mfalloc is obsolete: use mtaloc. ',n
      call abort
c      call chkflst(ilistroot,'mfalloc in')
c$$$      if(n .le. 0) then
c$$$         mfalloc=0
c$$$         if(n .eq. -1) then
c$$$            currp=ilistroot
c$$$ 1300       continue
c$$$            mfalloc=max(mfalloc,ilist(1,currp))
c$$$            if(ilist(2,currp) .ne. 0) then
c$$$               currp=ilist(2,currp)
c$$$               go to 1300
c$$$            endif
c$$$         elseif(n .ne. 0) then
c$$$            call errmsg('mfalloc','invalid argument for mfalloc',0,0)
c$$$            write(msgstr,*)'argument is ',n
c$$$            call errmsg('freeme',msgstr(:lene(msgstr)),0,0)
c$$$         endif
c$$$         return
c$$$      endif
c$$$
c$$$      currp=ilistroot
c$$$ 1000 continue
c$$$      next=ilist(2,currp)
c$$$      if(next .le. 0) then
c$$$         if(next .lt. 0) then
c$$$            call errmsg('mfalloc',' negative pointer',0,0)
c$$$         endif
c$$$         if(n .ge. pagesz) then
c$$$            msize=n
c$$$         else
c$$$            msize=pagesz
c$$$         endif
c$$$         next=pfalloc(msize)
c$$$         if(next .le. 0) then
c$$$            mfalloc=0
c$$$            return
c$$$         endif
c$$$         ilist(1,next)=msize
c$$$         ilist(2,next)=0
c$$$         call insflst(ilistroot,next,currp)
c$$$c     call chkflst(ilistroot,'mfalloc insert new block')
c$$$         if(ilist(1,next) .lt. msize) then
c$$$            call errmsg('mfalloc', '???',0,0)
c$$$            stop
c$$$         endif
c$$$      endif
c$$$      msize=ilist(1,next)
c$$$      if(msize .lt. 0) then
c$$$         call errmsg('mfalloc',' Abnormal memory block',32,0)
c$$$         stop
c$$$      endif
c$$$c     call chkflst(ilistroot,'mfalloc itermediate')
c$$$      if(msize .lt. n) then
c$$$         currp=next
c$$$         go to 1000
c$$$      endif
c$$$      if(msize .eq. n) then
c$$$         ilist(2,currp)=ilist(2,next)
c$$$         mfalloc=next
c$$$      else
c$$$         ilist(1,next)=msize-n
c$$$         mfalloc=next+ilist(1,next)
c$$$         ilist(1,mfalloc)=n
c$$$      endif
c$$$      call chkflst(ilistroot,'mfalloc out')
c$$$c      write(*,*)'@ mfalloc'
c$$$c      call pr_mem_map()
c$$$      return
      end

      subroutine freeme(pt,size)
      use maccbk
      use tfstk
      use macfile
      implicit none
      integer*4 pt,size
c      character*(MAXLLEN) msgstr
c      integer*4 currp,next
c      integer*4 lene

      write(*,*)'freeme is obsolete. use tfreem. ',pt,size
      call abort
c$$$      if(pt .eq. 0 .or. size .eq. 0) then
c$$$         return
c$$$      else if(pt .le. 0 .or. size .le. 0)then
c$$$         write(msgstr,*)' invalid parameter(s) freeme(',pt,size,')'
c$$$         call errmsg('freeme',msgstr(:lene(msgstr)),0,0)
c$$$         call pr_mem_map
c$$$         return
c$$$      endif
c$$$
c$$$      call chkflst(ilistroot,'mfree in')
c$$$      next=pt
c$$$      ilist(1,next)=size
c$$$      ilist(2,next)=0
c$$$      call insflst(ilistroot,next,currp)
c$$$      call chkflst(ilistroot,'mfree insert')
c$$$      if(ilist(1,next) .ge. pagesz) then
c$$$         ilist(2,currp)=ilist(2,next)
c$$$         call pfree(next,ilist(1,next))
c$$$      endif
c$$$      call chkflst(ilistroot,'mfree out')
c$$$c      write(*,*)'@ freeme'
c$$$c      call pr_mem_map()
c$$$      return
      end

c     Print small memory/page block list
      subroutine pr_mem_map()
      use maccbk
      use macfile
      implicit none

      write(errfl,*) ' Small memory block list'
      call pr_mem_list(ilistroot)
      write(errfl,*) ' Page block list'
      call pr_mem_list(ilistroot+1)
      return  
      end

c     Internal APIs
c     Allocate/Free large memory page
      integer*4 function pfalloc(size)
      use maccbk
      use macfile
      implicit none
      integer*4 size
      character*(MAXLLEN) msgstr
      integer*4 currp,next,psize,irtc,lene,lmalloc

      if(size .le. 0)then
         pfalloc = 0
         write(msgstr,*)'invalid argument ',size,' is given to pfalloc.'
         call errmsg('pfalloc',msgstr(:lene(msgstr)),0,0)
         return
      end if

      call chkflst(ilistroot+1,'pfalloc in')
      currp=ilistroot+1
 1000 continue
      next=ilist(2,currp)
c      write(*,*)'pfalloc-1 ',currp,next,ilist(1,currp),ilist(1,next)
      if(next .le. 0) then
         if(next .lt. 0) then
            write(msgstr,*)' negative pointer(next)',next
            call errmsg('pfalloc',msgstr(:lene(msgstr)),0,0)
            pfalloc=-1
            return
         endif
c        Insert new allocated page fragment
         next=lmalloc(size,irtc)
         if(irtc .ne. 0) then
            pfalloc=0
            return
         endif
         ilist(1,next)=size
         ilist(2,next)=0
         ilist(2,currp)=next
      endif
      psize=ilist(1,next)
      if(psize .le. 0) then
         write(*,*)'pfalloc ',size,next,ilistroot+1,currp,ilist(1,next)
         write(msgstr,*)'invalid data in page list. argument is ',size
     $        , ' (ilistroot+1,currp, next,mem left)=('
     $        ,ilistroot+1,currp,next,ilist(1,next),')'
         call errmsg('pfalloc',msgstr(:lene(msgstr)),0,0)
c         call pr_mem_map
         call chkflst(ilistroot+1,'pfalloc')
         stop
      endif

      if(psize .lt. size) then
c        Next page fragment size is too small. Check another page fragment!
         currp=next
         go to 1000
      endif

      if(psize .eq. size) then
c        Unlink next page fragment from page fragment list
         ilist(2,currp)=ilist(2,next)
         pfalloc=next
      else
c        Cut `size' page fragment from tail of next page fragment
         ilist(1,next)=psize-size
         pfalloc=next+ilist(1,next)
         if(ilist(1,next) .lt. pagesz) then
c           Move next page fragment into small page list
            ilist(2,currp)=ilist(2,next)
            call insflst(ilistroot,next,currp)
c            call chkflst(ilistroot,'pfalloc move small block')
c            call pr_mem_map()
         endif
      endif
      ilist(1,pfalloc)=size
      ilist(2,pfalloc)=0 
      call chkflst(ilistroot+1,'pfalloc out')
      return
      end

      subroutine pfree(pt,size)
      use maccbk
      use macfile
      implicit none
      integer*4 pt,size
      character*(MAXLLEN) msgstr
      integer*4 currp,next

      if(size .le. 0 .or.  pt .le. 0) then
         write(msgstr,*)' invalid parameter(s) pfree(',pt,size,')'
         call errmsg('pfree',msgstr,0,0)
         if(pt .eq. 0 .or. size .eq. 0) then
            call chkflst(ilistroot+1,'pfree zero pointer/size')
            call pr_mem_map
         else
            call chkflst(ilistroot+1,'pfree negative pointer/size')
         endif
         return
      endif

      call chkflst(ilistroot+1,'pfree enter')
      next=pt
      ilist(1,next)=size
      ilist(2,next)=0
      call insflst(ilistroot+1,next,currp)
      call chkflst(ilistroot+1,'pfree exit')
c      call pr_mem_map
      return
      end

c     Insert memory block to allocated block list
      subroutine insflst(pt,insb,prev)
      use maccbk
      use macfile
      implicit none
      integer*4 pt,insb,prev
      integer*4 parent,currp,next,size

      parent=pt
      currp=pt
      next=ilist(2,pt)
 1000 continue
      if(next .lt. 0) then
         call errmsg('insflst',' negative pointer',0,0)
         stop
      endif
      if((next .eq. 0) .or. ((next .gt. insb) .and.
     $     ((insb .gt. currp) .or. (currp .eq. pt)))) then
         ilist(2,insb)=next
         ilist(2,currp)=insb
      else
         parent=currp
         currp=next
         next=ilist(2,currp)
         go to 1000
      endif

      next=ilist(2,currp)
      if(ilist(2,next) .ne. 0) then
         size=ilist(1,next)
         if(size .le. 0) then
            call errmsg('insflst', ' bloken flist',0,0)
            stop
         endif
         if(next+size .eq. ilist(2,next)) then
c           Merge tailing page fragment into next page(insb)
            ilist(2,next)=ilist(2,next+size)
            ilist(1,next)=size+ilist(1,next+size)
         endif
      endif
      if(parent .ne. pt) then
         size=ilist(1,currp)
         if(size .le. 0) then
            call errmsg('insflst', ' bloken flist',0,0)
            stop
         endif
         if(currp+size .eq. next)then
c           Merge next page fragment(insb) into currp page
            ilist(2,currp)=ilist(2,next)
            ilist(1,currp)=size+ilist(1,next)
            next=currp
            currp=parent
         endif
      endif
      insb=next
      prev=currp
      return
      end

c     Check allocated block list structure
      subroutine chkflst(pt,caller0)
      use maccbk
      use macfile
      implicit none
      integer*4 pt
      character*(*) caller0
      integer*4 currp,next,size,np,l
      integer clen
      parameter(clen=64)
      character*(clen) caller

      l=min(clen,len(caller0))
      caller=caller0(:l)

      np=ilist(2,pt)
      currp=ilist(2,pt)
      next=ilist(2,pt)
 1000 continue
      if(next .eq. 0) then
         return
      endif
      if(next .lt. 0) then
         write(errfl,*)'chkflst(',pt,')'
         call errmsg('chkflst',' negative pointer '//caller(:l),0,0)
         call pr_mem_map()
         stop
      endif
      if(next .lt. currp) then
         write(errfl,*)'chkflst(',pt,')'
         call errmsg('chkflst',' inconsistent '//caller(:l),0,0)
         call pr_mem_map()
         stop
      endif
      size=ilist(1,next)
      if(size .le. 0) then
         write(errfl,*)'chkflst(',pt,')'
         if(size .eq. 0) then
            call errmsg('chkflst',' block size=0 '//caller(:l),0,0)
         else
            call errmsg('chkflst',' block size<0 '//caller(:l),0,0)
         endif
         call pr_mem_map()
         stop
      endif
      currp=next
      next=ilist(2,currp)
      if(next .eq. np) then
         write(errfl,*)'chkflst(',pt,')'
         call errmsg('chkflst',' looped flist '//caller(:l),0,0)
         stop
      endif
      go to 1000
      end

c     Print allocated block list
      subroutine pr_mem_list(flist)
      use maccbk
      use macfile
      implicit none
      integer*4 flist
      character*80 line
      integer*4 currp,cnt,memrst,memmax

      memmax=0
      memrst=0
      currp=ilist(2,flist)
      cnt=1
 1000 continue
      write(line(cnt:cnt+15),'(I7,A1,I7,A1)')
     &     currp,' ',ilist(1,currp),'!'
      if(currp .ne. 0) then
         memrst=memrst+ilist(1,currp)
         memmax=max(memmax,ilist(1,currp))
         cnt=cnt+16
         if(cnt .gt. 80 - 16 ) then
            write(errfl,*)line(:cnt-1)
            cnt=1
         endif
         if(ilist(2,currp) .ne. 0) then
            if(ilist(2,currp) .lt. currp) then
               write(errfl,*)' invalid memory map',currp,ilist(2,currp),
     $              ilist(1,ilist(2,currp))
               return
            endif
            currp=ilist(2,currp)
            go to 1000
         endif
      endif
      write(errfl,*)line(:cnt-1)
      write(errfl,*)'Total free area = ',memrst,
     $     ' Maximum block size =',memmax
      return
      end

c     End of File
