      integer*4 function ielm(latt,word,istart,mult,exist)
      implicit none
      integer*4 latt(2,1),istart,mult(1),ielmf
      character*(*) word
      real*8 frac
      logical*4 exist
      ielm=ielmf(word,frac,exist)
      return
      end

      integer*4 function ielmf(word,frac,exist)
      use tfstk
      use ffs
      use tffitcode
      implicit none
      type (sad_descriptor) kx
      integer*4 lw,iord,ln,i,j,lenw,ip,im
      character*(*) word
      character*64 ordw
      character*(MAXPNAME) name
      real*8 frac
      integer*4 ioff,m,ipm, irtc,idot,ielmh,ist1
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
          ist1=1
          call tfeval(ordw,lenw(ordw),ist1,m,kx,.false.,irtc)
          iord=int(rfromd(kx))
          if(irtc .ne. 0 .or. ktfnonrealqd(kx))then
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
        if(irtc .ne. 0 .or. ktfnonrealqd(kx))then
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
          ielmf=nlat
          exist=.false.
          return
        endif
      endif
      if(trpt)then
        ielmf=min(max(i+ioff,1),nlat)
      else
        j=i+ioff
        do while(j .le. 0)
          j=j+nlat
        enddo
        ielmf=mod(j-1,nlat)+1
      endif
      exist=.true.
      return
      end

      integer function igelm(latt,word,mult,exist)
      use tfstk
      use ffs
      use tffitcode
      implicit none
      integer*4 latt(2,nlat),mult(nlat),ielm
      character*(*) word
      character*255 wordp
      logical exist
      call getwdl2(word,wordp)
      igelm=ielm(latt,wordp,1,mult,exist)
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
      implicit none
      character*(*) name
      integer*4 iord
      integer*4 itehash,k,j,i
      character*(MAXPNAME) name1
      name1=name
      k=itehash(name1,MAXPNAME)
      j=ilist(2,ielmhash+k)
      if(j .ne. 0)then
        do i=j,j+ilist(1,ielmhash+k)-1
          ielmh=ilist(1,i)
          if(name1 .eq. pname(ilist(1,ilattp+ielmh)))then
            if(ilist(ielmh,ifmult) .eq. iord)return
            if(iord .ne. 0)cycle
            if(ilist(ilist(ielmh,ifele1),ifklp) .eq. ielmh)return
c     Note: ilist(i, ifmult) == 0 if ilist(ilist(i, ifele1), ifklp) == i
c     *     by tfinit(), tfinimult() initialization
          endif
        enddo
      endif
      ielmh=0
      return
      end

      subroutine tfhashelement(latt)
      use tfstk
      use ffs
      use tffitcode
      implicit none
      integer*4 latt(2,nlat),itehash,nelm(0:nelmhash),
     $     i,j,k,l,n,italoc
      do j=0,nelmhash
        nelm(j)=0
      enddo
      do i=1,nlat-1
        j=itehash(pname(latt(1,i)),MAXPNAME)
        nelm(j)=nelm(j)+1
      enddo
      k=italoc(nelmhash+1)
      ilist(2,k-1)=ielmhash
      do j=0,nelmhash
        n=nelm(j)
        ilist(1,k+j)=0
        if(n .gt. 0)then
          ilist(2,k+j)=italoc(n)
        else
          ilist(2,k+j)=0
        endif
      enddo
      do i=1,nlat-1
        j=itehash(pname(latt(1,i)),MAXPNAME)
        l=ilist(2,k+j)+ilist(1,k+j)
        ilist(1,l)=i
        ilist(1,k+j)=ilist(1,k+j)+1
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
      integer*1 name(nc)
      ih=0
      do i=1,nc
        if(name(i) .eq. ichar(' '))then
          itehash=iand(ih,nh)
          return
        endif          
        ih=ih+name(i)
      enddo
      itehash=iand(ih,nh)
      return
      end

      subroutine tfresethash
      use tfstk
      use ffs
      use tffitcode
      implicit none
      integer*4 i,j
      j=ielmhash
      do i=0,nelmhash
        if(ilist(1,j+i) .ne. 0)then
          call tfree(int8(ilist(2,j+i)))
        endif
      enddo
      ielmhash=ilist(2,j-1)
      call tfree(int8(j))
      return
      end
