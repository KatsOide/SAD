      integer*4 function ielm(word,exist)
      implicit none
      integer*4 ielmf
      character*(*) ,intent(in):: word
      real*8 frac
      logical*4 ,intent(out):: exist
      ielm=ielmf(word,frac,exist,0)
      return
      end

      integer*4 function ielmex(word,exist,lfn) result(iv)
      use tfstk
      use tmacro, only:nlat
      implicit none
      integer*4 ,intent(in):: lfn
      integer*4 ielme,irtc,lw
      integer*8 iep
      character*(*) ,intent(in):: word
      real*8 v
      logical*4 ,intent(out):: exist
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
            iv=int(merge(v+0.499,nlat+1+v+0.5d0,v .ge. 0.d0))
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
      integer*4 ,intent(in):: lfn
      integer*4 ielmf
      character*(*) ,intent(in):: word
      real*8 frac
      logical*4 ,intent(out):: exist
      ielme=ielmf(word,frac,exist,lfn)
      return
      end

      integer*4 function ielmf(word,frac,exist,lfn)
      use tfstk
      use ffs
      use tffitcode
      use ffs_pointer,only:ielma
      implicit none
      type (sad_descriptor) kx,tfeval
      integer*4 ,intent(in):: lfn
      integer*4 lw,iord,ln,i,ip,im
      character*(*) ,intent(in):: word
      character*64 ordw
      character*(MAXPNAME) name
      real*8 ,intent(out):: frac
      integer*4 ioff,m,ipm, irtc,idot,ielmh
      logical*4 ,intent(out):: exist
      lw=len_trim(word)
      idot=index(word(1:lw),'.')
      if(idot .gt. 0)then
        ln=idot-1
        ip=index(word(1:lw),'+')
        im=index(word(1:lw),'-')
        if(ip .le. 0)then
          ipm=im
        else
          ipm=merge(ip,min(ip,im),im .le. 0)
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
          kx=tfeval(ordw,1,m,.false.,irtc)
          iord=int(kx%x(1))
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
          ipm=merge(ip,min(ip,im),im .le. 0)
        endif
        ln=merge(ipm-1,lw,ipm .gt. 0)
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
        ioff=int(kx%x(1))
        frac=kx%x(1)-ioff
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
            call termes(lfn,'?Undefined location ',word(1:lw))
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
      character*(*) ,intent(in):: word
      character*256 wordp
      logical*4 ,intent(out):: exist
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
      integer*4 ,intent(in):: lfn
      integer*4 ielme
      character*(*) ,intent(in):: word
      character*256 wordp
      logical*4 ,intent(out):: exist
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
      character*(*) ,intent(in):: name
      integer*8 j,i
      integer*4 ,intent(in):: iord
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
        klist(k+j*2+2)=merge(ktaloc(n),i00,n .gt. 0)
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
      integer*4 ,intent(in):: nc
      integer*4 i,ih,nh
      parameter (nh=nelmhash)
      character ,intent(in):: name(nc)
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
