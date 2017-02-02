      subroutine tffile(word,lfnstk,lfopen,lfret,lfrecl,lflinep,
     $     maxlfn,lfni,lfno,lfnb,init,exist)
      use tfstk
      use ffs
      use tffitcode
      implicit none
      type (sad_descriptor) kx
      type (sad_string), pointer :: str
      real*8 vx
      integer*4 maxlfn,lfno,i,lfni0,lfni1,nc,next,icsmrk,icslinep,
     $     icslrecl,itype,lfno1,j
      integer*4 lfnstk(maxlfn),lfret(maxlfn),lfrecl(maxlfn),
     $     itfpeeko,lfnb,lfni,lflinep(maxlfn)
      logical*4 lfopen(maxlfn),init,exist,termin,rew,app,abbrev,
     $     clo,ret
      character*(*) word
      character*255 tfconvstr
      exist=.true.
      init=.false.
      clo=.false.
      if(word .eq. 'END')then
        init=lfnp .gt. lfnb
        lfni0=lfni
        call tfclose(lfnb,int(lfnp),lfnstk,lfopen,lfret,lfrecl,
     $       lflinep,maxlfn,lfni,lfnb)
        lfnp=lfnb
        lfno=6
        outfl=lfno
        if(lfnb .eq. 1)then
          close(98)
        endif
      elseif(abbrev(word,'TERM_INATE','_') .or.
     $       abbrev(word,'CLO_SE','_'))then
        clo=abbrev(word,'CLO_SE','_')
        call peekwd(word,next)
        if(abbrev(word,'IN_PUT','_'))then
          termin=.true.
          call cssetp(next)
        elseif(abbrev(word,'OUT_PUT','_'))then
          termin=.false.
          call cssetp(next)
          if(clo)then
            close(lfno)
          else
            close(98)
          endif
        else
          termin=.true.
        endif
        if(termin)then
          call tfclose(int(lfnp),int(lfnp),lfnstk,lfopen,lfret,lfrecl,
     $         lflinep,maxlfn,lfni,lfnb)
        else
          lfno=6
          outfl=lfno
        endif
        init=termin
      elseif(abbrev(word,'IN_PUT','_')  .or. word .eq. 'GET'
     1       .or. word .eq. 'READ')then
        rew=word .eq. 'READ'
        itype=itfpeeko(kx,next)
        if(ktfrealqd(kx,vx))then
          lfni1=int(vx+.5d0)
          call cssetp(next)
          write(word,'(''ftn'',i2.2)')lfni1
          lfnp=lfnp+1
          lfopen(lfnp)=.false.
        elseif(ktfstringqd(kx))then
          call cssetp(next)
          word=tfconvstr(kx,nc,'*')
          do 8020 j=51,97
            do 8030 i=lfnb,int(lfnp)
              if(j .eq. lfnstk(i))then
                go to 8020
              endif
 8030       continue
            lfni1=j
            go to 8021
 8020     continue
          goto 8101
 8021     if(word .eq. ' ')then
            call termes(lfno,'?Missing filename for IN_PUT',' ')
            init=.true.
            return
          endif
          call texpfn(word)
          open(lfni1,file=word,status='OLD',ERR=8101)
          lfnp=lfnp+1
          lfopen(lfnp)=.true.
        else
          call termes(lfno,'?Missing filename for IN_PUT',' ')
          init=.true.
          return
        endif
        if(lfnp .ge. maxlfn)then
          call termes(lfno,'?Number of input files exceeds limit',' ')
          init=.true.
          return
        endif
        lfnstk(lfnp)=lfni1
        lfret(lfnp)=icsmrk()
        lfrecl(lfnp)=icslrecl()
        lflinep(lfnp)=icslinep()
        call cssetp(lfrecl(lfnp))
        call cssetlinep(lfrecl(lfnp))
        lfni=lfni1
        if(rew)then
          rewind(lfni)
        endif
        init=.true.
        return
 8101   call termes(lfno,'?File open error ',word)
        call cssets(9998)
        init=.true.
      elseif(abbrev(word,'EXE_CUTE','_'))then
        ret=.true.
        itype=itfpeeko(kx,next)
        if(.not. ktfstringqd(kx,str))then
          call termes(lfno,'?Missing string for EXE_CUTE.',' ')
          init=.true.
          return
        endif
        call cssetp(next)
        if(lfnp .ge. maxlfn)then
          call termes(lfno,'?Number of input files exceeds limit',' ')
          init=.true.
          return
        endif
        lfnp=lfnp+1
        lfni=0
        lfopen(lfnp)=.false.
        lfnstk(lfnp)=0
        lfret(lfnp)=icsmrk()
        lfrecl(lfnp)=icslrecl()
        lflinep(lfnp)=icslinep()
        call cssetp(lfrecl(lfnp))
        call cssetlinep(lfrecl(lfnp))
        call setbuf(str%str,str%nch)
        init=.false.
        return
      elseif(abbrev(word,'OUT_PUT','_') .or. word .eq. 'PUT'
     1       .or. abbrev(word,'APP_END','_'))then
        app=abbrev(word,'APP_END','_')
        itype=itfpeeko(kx,next)
c        write(*,*)'tffile-0 ',itype,ia
        if(ktfrealqd(kx,vx))then
          lfno1=int(vx+.5d0)
          call cssetp(next)
        elseif(ktfstringqd(kx))then
          lfno1=98
          close(lfno1)
          call cssetp(next)
          word=tfconvstr(kx,nc,'*')
c          write(*,*)'tffile-1 ',word(1:len_trim(word))
          call texpfn(word)
c          write(*,*)'tffile-2 ',word(1:len_trim(word))
          if(app)then
            open(lfno1,file=word,status='UNKNOWN',
     1           access='APPEND',ERR=6101)
          else
            open(lfno1,file=word,status='UNKNOWN',
     1           ERR=6101)
          endif
        else
          call termes(lfno,'?Missing filename for OUT_PUT',' ')
          init=.true.
          return
        endif
        lfno=lfno1
        outfl=lfno
        exist=.true.
        return
 6101   lfno=6
        call termes(lfno,'?File open error ',word)
        init=.true.
        call cssets(9998)
      else
        exist=.false.
      endif
      return
      end

      subroutine texpfn(file)
      implicit none
      character*(*) file
      integer*4 l
c
      l=len_trim(file)
      if(file(1:1) .eq. '''' .or. file(1:1) .eq. '"')then
        file=file(2:l-1)
      endif
      call cfexptilde(file)
      return
      end

      subroutine tfclose(lfnp1,lfnp,lfnstk,lfopen,lfret,lfrecl,
     $     lflinep,maxlfn,lfni,lfnb)
      implicit none
      integer*4 lfnp1,lfnp,maxlfn,lfnstk(maxlfn),lfni0,
     $     lfni,i,lfret(maxlfn),lfrecl(maxlfn),lfnp0,lfnb,
     $     lflinep(maxlfn)
      logical*4 lfopen(maxlfn)
      lfni0=lfni
      do 10 i=lfnp1,lfnp
        if(lfopen(i))then
          if(lfnstk(i) .gt. 0)then
            close(lfnstk(i))
          endif
        endif
10    continue
      lfnp0=max(1,lfnb-1,lfnp1-1)
      lfni=lfnstk(lfnp0)
      if(lfret(lfnp1) .gt. 0)then
        call cssetp(lfret(lfnp1))
        call cssetl(lfrecl(lfnp1))
        call cssetlinep(lflinep(lfnp1))
      elseif(lfni0 .ne. lfni)then
        call skipline
      endif
      lfnp=lfnp0
      return
      end
