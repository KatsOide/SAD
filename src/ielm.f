      integer*4 function ielm(word,exist)
      implicit none
      integer*4 ielmf
      character*(*) word
      real*8 frac
      logical*4 exist
      ielm=ielmf(word,frac,exist,0)
      return
      end

      integer*4 function ielmex(word,exist,lfn) result(iv)
      use tfstk
      use tmacro, only:nlat
      implicit none
      integer*4 ielme,lfn,irtc,lw
      integer*8 iep
      character*(*) word
      real*8 v
      logical*4 exist
      type (sad_descriptor) kx
      exist=.false.
      iv=0
      iep=ierrorprint
      ierrorprint=0
      irtc=0
      lw=len_trim(word)
      if(lw .gt. 0)then
        if(word(1:1) .ne. "^")then
          call tfevalb(word,kx,irtc)
          ierrorprint=iep
          if(irtc .eq. 0 .and. ktfrealq(kx,v))then
            if(v .ge. 0.d0)then
              iv=int(v+0.499)
            else
              iv=int(nlat+1+v+0.5d0)
            endif
            iv=max(1,min(nlat,iv))
            exist=.true.
            return
          endif
          if(irtc .ne. 0)then
            call tfreseterror
          endif
        endif
        iv=ielme(word,exist,lfn)
      endif
      return
      end function

      integer*4 function ielme(word,exist,lfn)
      implicit none
      integer*4 ielmf,lfn
      character*(*) word
      real*8 frac
      logical*4 exist
      ielme=ielmf(word,frac,exist,lfn)
      return
      end

      integer*4 function ielmf(word,frac,exist,lfn)
      use tfstk
      use ffs
      use tffitcode
      use ffs_pointer,only:ielma
      implicit none
      type (sad_descriptor) kx
      integer*4 lw,iord,ln,i,ip,im,lfn
      character*(*) word
      character*64 ordw
      character*(MAXPNAME) name
      real*8 frac
      integer*4 ioff,m,ipm, irtc,idot,ielmh
      logical*4 exist
      lw=len_trim(word)
      idot=index(word(1:lw),'.')
      if(idot .gt. 0)then
        ln=idot-1
        ip=index(word(1:lw),'+')
        im=index(word(1:lw),'-')
        if(ip .le. 0)then
          ipm=im
        else
          if(im .le. 0)then
            ipm=ip
          else
            ipm=min(ip,im)
          endif
        endif
        if(ipm .gt. 0 .and. ipm .lt. idot)then
          ln=ipm-1
          iord=0
        else
          if(ipm .gt. 0)then
            ordw=word(idot+1:ipm-1)
          else
            ordw=word(idot+1:lw)
          endif
          call tfeval(ordw,1,m,kx,.false.,irtc)
          iord=int(rfromd(kx))
          if(irtc .ne. 0 .or. ktfnonrealq(kx))then
            if(irtc .gt. 0 .and. ierrorprint .ne. 0)then
              call tfreseterror
            endif
            ielmf=1
            exist=.false.
            return
          endif
        endif
      else
        iord=0
        ip=index(word(1:lw),'+')
        im=index(word(1:lw),'-')
        if(ip .le. 0)then
          ipm=im
        else
          if(im .le. 0)then
            ipm=ip
          else
            ipm=min(ip,im)
          endif
        endif
        if(ipm .gt. 0)then
          ln=ipm-1
        else
          ln=lw
        endif
      endif
      name=word(1:ln)
      ielmf=0
      if(ipm .ne. 0)then
        call tfevals(word(ipm:lw),kx%k,irtc)
        if(irtc .ne. 0 .or. ktfnonrealq(kx))then
          if(irtc .gt. 0 .and. ierrorprint .ne. 0)then
            call tfreseterror
          endif
          ielmf=1
          exist=.false.
          return
        endif
        ioff=int(rfromd(kx))
        frac=rfromd(kx)-ioff
      else
        ioff=0
        frac=0.d0
      endif
      if(name .eq. '$$$' .or. name .eq. '***')then
        i=nlat
      elseif(name .eq. '^^^')then
        i=1
      else
        i=ielmh(name,iord)
        if(i .eq. 0)then
          if(lfn .ne. 0 .and. idot .gt. 0)then
            call termes(lfn,'?Undefined location ',
     $           word(1:lw))
          endif
          ielmf=nlat
          exist=.false.
          return
        endif
      endif
      ielmf=ielma(i+ioff)
      exist=.true.
      return
      end

      integer*4 function igelm(word,exist)
      use tfstk
      use ffs
      use ffs_pointer
      use tffitcode
      implicit none
      integer*4 ielm
      character*(*) word
      character*256 wordp
      logical exist
      call getwdl2(word,wordp)
      igelm=ielm(wordp,exist)
      return
      end

      integer*4 function igelme(word,exist,lfn)
      use tfstk
      use ffs
      use ffs_pointer
      use tffitcode
      implicit none
      integer*4 ielme,lfn
      character*(*) word
      character*256 wordp
      logical exist
      call getwdl2(word,wordp)
      igelme=ielme(wordp,exist,lfn)
      return
      end

c     Search beamline element position
c     * name     : element name
c     * iord    0: singlet element or head of multiple elements
c     *        >0: number suffix of multiple elements
c     *
c     * result  0: not found
c     *        >0: found
      integer function ielmh(name,iord)
      use tfstk
      use ffs
      use tffitcode
      use ffs_pointer, only:idelc,pnamec
      implicit none
      character*(*) name
      integer*8 j,i
      integer*4 iord
      integer*4 itehash,k
      character*(MAXPNAME) name1
      name1=name
      k=itehash(name1,MAXPNAME)*2
      j=klist(ielmhash+k+2)
      if(j .ne. 0)then
        do i=j,j+ilist(1,ielmhash+k+1)-1
          ielmh=ilist(1,i)
          if(name1 .eq. pnamec(ielmh))then
            if(ilist(ielmh,ifmult) .eq. iord)return
            if(iord .ne. 0)cycle
            if(nelvx(ilist(ielmh,ifele1))%klp .eq. ielmh)return
c     Note: ilist(i, ifmult) == 0 if ilist(ilist(i, ifele1), ifklp) == i
c     *     by tfinit(), tfinimult() initialization
          endif
        enddo
      endif
      ielmh=0
      return
      end

      subroutine tfhashelement
      use tfstk
      use ffs
      use ffs_pointer
      use tffitcode
      implicit none
      integer*8 k
      integer*4 itehash,nelm(0:nelmhash),j,n,i
      nelm=0
      do i=1,nlat-1
        j=itehash(pnamec(i),MAXPNAME)
        nelm(j)=nelm(j)+1
      enddo
      k=ktaloc((nelmhash+1)*2+1)
      klist(k)=ielmhash
      do j=0,nelmhash
        n=nelm(j)
        ilist(1,k+j*2+1)=0
        if(n .gt. 0)then
          klist(k+j*2+2)=ktaloc(n)
        else
          klist(k+j*2+2)=0
        endif
      enddo
      do i=1,nlat-1
        j=itehash(pnamec(i),MAXPNAME)*2
        ilist(1,klist(k+j+2)+ilist(1,k+j+1))=i
        ilist(1,k+j+1)=ilist(1,k+j+1)+1
      enddo
      ielmhash=k
      return
      end

      integer*4 function itehash(name,nc)
      use tfstk
      use ffs
      use tffitcode
      implicit none
      integer*4 nc,i,ih,nh
      parameter (nh=nelmhash)
      character name(nc)
      ih=0
      do i=1,nc
        if(name(i) .eq. ' ')then
          exit
        endif          
        ih=ih+ichar(name(i))
      enddo
      itehash=iand(ih,nh)
      return
      end

      subroutine tfresethash
      use tfstk
      use ffs
      use tffitcode
      implicit none
      integer*8 j
      integer*4 i
      j=ielmhash
      do i=0,nelmhash*2,2
        if(ilist(1,j+i+1) .ne. 0)then
          call tfree(klist(j+i+2))
        endif
      enddo
      ielmhash=klist(j)
      call tfree(j)
      return
      end
