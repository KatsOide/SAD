      subroutine tffile(word,lfnb,init,exist)
      use tfstk
      use ffs
      use tffitcode
      use readbuf, only:trbopenmap
      use tfrbuf
      use tfcsi
      use ffsfile
      implicit none
      type (sad_descriptor) kx
      type (sad_string), pointer :: str
      real*8 vx
      integer*4 , intent(in)::lfnb
      integer*4 i,lfni1,nc,next,itype,lfno1,isp0,irtc,
     $     itfdownlevel,itfpeeko
      integer*4 ,save:: lfni0=0
      logical*4 ,intent(out)::init,exist
      logical*4 termin,rew,app,abbrev,clo,ret
      character*(*) word
      character*256 tfconvstr
      exist=.true.
      init=.false.
      clo=.false.
      if(abbrev(word,'SUSP_END','_') .or. word .eq. 'END')then
        if(suspend)then
          init=lfnp .gt. lfnb
          lfni0=lfni
          call tfclose(lfnb,lfnb)
          lfnp=lfnb
          lfno=6
          outfl=lfno
          call trbassign(lfni)
        endif
        return
      elseif(abbrev(word,'TERM_INATE','_') .or.
     $       abbrev(word,'CLO_SE','_'))then
        clo=abbrev(word,'CLO_SE','_')
        call peekwd(word,next)
        if(abbrev(word,'IN_PUT','_'))then
          termin=.true.
          ipoint=next
        elseif(abbrev(word,'OUT_PUT','_'))then
          termin=.false.
          ipoint=next
          if(clo)then
            close(lfno)
          else
            close(98)
          endif
        else
          termin=.true.
        endif
        if(termin)then
          call tfclose(int(lfnp),lfnb)
        else
          lfno=6
          outfl=lfno
        endif
        init=termin
        return
      elseif(abbrev(word,'EXE_CUTE','_'))then
        ret=.true.
        itype=itfpeeko(kx,next)
        if(.not. ktfstringq(kx,str))then
          call termes(lfno,'?Missing string for EXE_CUTE.',' ')
          init=.true.
          return
        endif
        ipoint=next
        isp0=isp
        isp=isp+1
        dtastk(isp)=kx
        isp=isp+1
        rtastk(isp)=dble(outfl)
        levele=levele+1
        call tfffs(isp0,kx,irtc)
        i=itfdownlevel()
        isp=isp0
        return
      elseif(abbrev(word,'OUT_PUT','_') .or. word .eq. 'PUT'
     1       .or. abbrev(word,'APP_END','_'))then
        app=abbrev(word,'APP_END','_')
        itype=itfpeeko(kx,next)
        if(ktfrealq(kx,vx))then
          lfno1=int(vx+.5d0)
          ipoint=next
        elseif(ktfstringq(kx))then
          lfno1=98
          close(lfno1)
          ipoint=next
          word=tfconvstr(kx,nc,'*')
          call texpfn(word)
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
        ios=9998
        return
      elseif(abbrev(word,'RES_UME','_'))then
        lfni1=lfni0
        if(lfni1 .eq. 0)then
          return
        endif
        lfni0=0
        call trbassign(lfni1)
        lfnp=lfnp+1
        rew=.false.
      elseif(abbrev(word,'IN_PUT','_')  .or. word .eq. 'GET'
     1       .or. word .eq. 'READ')then
        rew=word .eq. 'READ'
        itype=itfpeeko(kx,next)
c        call tfdebugprint(kx,'IN',1)
        if(ktfrealq(kx,vx))then
          lfni1=int(vx+.5d0)
        elseif(ktfstringq(kx,str))then
          call trbopenmap(str%str(1:str%nch),kx,irtc)
          if(Irtc .ne. 0)then
            call termes(lfno,'?File open error for IN_PUT',
     $           str%str(1:str%nch))
            return
          endif
          lfni1=int(rfromd(kx))
        else
          call termes(lfno,'?Missing filename for IN_PUT',' ')
          init=.true.
          return
        endif
        ipoint=next
        if(lfnp .ge. maxlfn)then
          call termes(lfno,'?Number of input files exceeds limit',' ')
          init=.true.
          return
        endif
        if(lfni1 .le. 0)then
          return
        endif
c     write(word,'(''ftn'',i2.2)')lfni1
        call trbassign(lfni1)
        lfnp=lfnp+1
      else
        exist=.false.
        return
      endif
      lfnstk(lfnp)=lfni1
      lfni=lfni1
      if(rew)then
        rewind(lfni)
      endif
      init=.true.
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

      subroutine tfclose(lfnp1,lfnb)
      use tfcsi
      use tfrbuf
      use ffsfile
      implicit none
      integer*4 lfnp1,lfni0,lfnp0,lfnb
      lfni0=lfni
      lfnp0=max(1,lfnb-1,lfnp1-1)
      lfni=lfnstk(lfnp0)
      if(lfni .ne. lfni0)then
        call trbassign(lfni)
      endif
      if(lfni0 .ne. lfni)then
        call skipline
      endif
      lfnp=lfnp0
      return
      end
