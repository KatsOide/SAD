      subroutine phdrw(ix,np,word,title,case,exist,lfno)
      use tfstk
      use tmacro
      implicit real*8(a-h,o-z)
      integer*4 ax(2)
      character*(*) word,title,case
      character*8 coord(6)
      logical*4 abbrev,exist,fconv,lod,append
      real*8 s(8)
      data coord  /'X       ','PX      ','Y       ',
     z             'PY      ','Z       ','PZ      '/
      fconv=.false.
      lfn=lfno
c .... parsing word
      if(.not.exist) goto 2
    1 call getwdl(word)
    2 continue
      ax(1)=0
      ax(2)=0
      do 11 j=1,3
        if(word(j:j).eq.'-') then
          do 12 k=1,6
            if(word(:j-1).eq.coord(k)(:j-1)) then
              ax(1)=k
            endif
            if(word(j+1:j+2).eq.coord(k)(:2)) then
              ax(2)=k
            endif
   12     continue
          goto 19
        endif
   11 continue
      do 10 i=1,6
        if(word.eq.coord(i)) then
          ax(1)=i
          ax(2)=i+2*(i-2*(i/2))-1
          goto 19
        endif
   10 continue
      if(word.eq.'OUT'.or.abbrev(word,'APP_END','_')) then
        append=abbrev(word,'APP_END','_')
        lfn=int(getva(exist))
        if(exist)then
          if(lfn .le. 0)then
            lfn=6
          endif
          goto 1
        else
          lfn=98
          inquire(unit=lfn,opened=lod,iostat=ios)
          if(lod) close(lfn)
          call getwrd(word)
          call texpfn(word)
          if(append) then
            open(lfn,file=word,status='UNKNOWN',
     1           access='APPEND',ERR=6101)
          else
            open(lfn,file=word,status='UNKNOWN',
     1           ERR=6101)
          endif
          goto 1
        endif
6101    lfn=6
        call termes(lfno,'?File open error ',word)
        goto 999
      endif
   19 continue
c
      if(ax(1)*ax(2).eq.0) then
        exist=.false.
        goto 999
      endif
c .... convert to canonical variables
      if(.not.fconv) then
        ix1=italoc(np0*8)
        fconv=.true.
      endif
      do 100 i=1,np
        do 101 j=1,8
  101     s(j)=rlist(ix+np0*(j-1)+i-1)
        call tconv(s,s,1)
        do 102 j=1,8
  102     rlist(ix1+np0*(j-1)+i-1)=s(j)
  100 continue
      call phdrwa(rlist(ix1+np0*(ax(2)-1)),rlist(ix1+np0*(ax(1)-1)),
     z            np,ax,title,case,lfn)
      goto 1
c
  999 continue
      if(fconv) then
        call tfree(int8(ix1))
      endif
      if(lfn .ne. lfno) close(lfn)
      return
      end
