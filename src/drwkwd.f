      subroutine drwkwd(word,wordp,kwd,idim,ipw,sc,nkey,lfno,exist,err)
      use maccbk
      use tfcsi,only:icslfnm
      implicit none
C*DEC added by Y.Tange 10-Jan-1995
      integer*4 nkey
C*DEC End
      character*(*) word,wordp,kwd(nkey)
      character*(8) word1
      real*8 v,sc(2,nkey),dflt
C*DEC changed by Y.Tange 10-Jan-1995
      integer*4 idim(nkey),ipw(2,nkey),lfno
C*DEC End
      integer*4 ip,id,i,j,ia,icat1,ipage,iwin,ifany,lw,lenw
      logical*4 exist,err
      external trim
C*DEC changed
      parameter (dflt=-1d20)
C*DEC
c     begin initialize for preventing compiler warning
      ia=0
c     end   initialize for preventing compiler warning
      err=.false.
      do 3 i=1,nkey
        ipw(1,i)=0
        ipw(2,i)=0
        sc(1,i)=dflt
        sc(2,i)=-dflt
3     continue
      ipage=1
      iwin=1
      go to 4
2     call getwdl2(word,wordp)
4     if(word .eq. ' ')then
        exist=.true.
        return
      endif
      lw=lenw(word)
      ip=1
1     id=ifany(word,'&/',ip)
      if(id .le. 0)then
        word1=word(ip:)
        ip=lw+1
      else
        if(id .eq. ip)then
          word1=word(ip:ip)
          if(word1 .eq. '&')then
            iwin=iwin+1
          else
            ipage=ipage+1
            iwin=1
          endif
          ip=ip+1
          go to 1
        else
          word1=word(ip:id-1)
          ip=id
        endif
      endif
      if(word1 .eq. ' ')then
        go to 2
      endif
      exist=.false.
      do 10 i=1,nkey
        if(kwd(i) .eq. word1)then
          ia=i
          exist=.true.
          if(ipw(1,i) .ne. 0)then
            call termes('?Duplicate graph ',word1)
            err=.true.
            return
          endif
          icat1=-9999
          do 20 j=1,nkey
            if(ipw(1,j) .eq. ipage .and. ipw(2,j) .eq. iwin)then
              if(idim(j) .eq. idim(i))then
                go to 21
              else
                if(icat1 .eq. -9999)then
                  icat1=idim(j)
                elseif(icat1 .ne. idim(j))then
                  call termes('?Too many scales in a window ',word1)
                  err=.true.
                  return
                endif
              endif
            endif
20        continue
21        ipw(1,i)=ipage
          ipw(2,i)=iwin
          go to 1
        endif
10    continue
      if(word1.eq.'@') then
        if(sc(1,ia).eq.dflt) then
          sc(1,ia)=2*dflt
        else
          sc(2,ia)=2*dflt
        endif
        goto 1
      endif
      call trim(word1)
      if(index('-+.01234567890',word1(1:1)) .gt. 0)then
        read(word1,*,err=11) v
        if(sc(1,ia).eq.dflt) then
          if(ia.eq.2 .or. ia.eq.5) then
            sc(1,ia)=sqrt(max(0d0,v))
          else
            sc(1,ia)=v
          endif
        else
          sc(2,ia)=sc(1,ia)
          if(ia.eq.2 .or. ia.eq.5) then
            sc(1,ia)=sqrt(max(0d0,v))
          else
            sc(1,ia)=v
          endif
        endif
        goto 1
      endif
11    exist=.false.
      return
      end
