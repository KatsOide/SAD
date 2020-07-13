      Integer*8 function mkplst(idxl)
      use maccbk
      use mackw
      use tfmem, only:ktaloc,tfree
      implicit real*8 (a-h,o-z)
      integer*8 idxl,i
      integer*4 STKSIZ,allmem
      parameter (STKSIZ =1024, allmem=32768)
      integer mlen,isp
      integer*8 pstack,mtop
c
      pstack=ktaloc(STKSIZ)
      isp=0
c      allmem=mtaloc(-1)
      mtop=ktaloc(allmem)
      mlen=0
      call push(int8(0),0,mlen,ilist(1,mtop),allmem)
c
      level=0
 1000 continue
      llen=ilist(1,idxl)
      i=idxl+1
 1100 continue
      if(llen .eq. 0) then
        if(level .eq. 0) go to 2000
c
        level=level-1
        call pop(i,llen,isp,ilist(1,pstack),STKSIZ)
      else
        if(lstchk(ilist(2,i),ilist(1,mtop),mlen,allmem) .eq. 0) then
          call push(int8(0),ilist(2,i),mlen,ilist(1,mtop),allmem)
        endif
        if (idtype(ilist(2,i)) .eq. icLINE) then
          call push(i,llen,isp,ilist(1,pstack),STKSIZ)
          idxl=idval(ilist(2,i))
          level=level+1
          go to 1000
        endif
      endif
      llen=llen-1
      i=i+1
      go to 1100
c
 2000 continue
      ilist(1,mtop)=mlen
      ilist(2,mtop)=0
      if(allmem .gt. mlen+3)then
        ilist(1,mtop+mlen)=allmem-mlen
        call tfree(mtop+mlen+1)
      endif
      mkplst=mtop
c     print *,'mklst',mtop,mlen
      return
      end
