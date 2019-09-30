      subroutine tfdbun(word,bunch,nbunch,lfno,err)
      use ffs
      use tffitcode
      real*8 bunch(nbunch),prn(nbunch)
      character*(*) word
      logical*4 err,abbrev,exist,add
      add=.false.
1     call peekwd(word,next)
      if(abbrev(word,'R_ANDOM','_'))then
        ipoint=next
        r=getva(exist)
        if(.not. exist)then
          call termes(lfno,'?Missing value for DBUNCH R_ANDOM.',' ')
          err=.true.
          return
        endif
        call tgauss_array(prn(1),nbunch)
        if(add)then
          do i=1,nbunch
            bunch(i)=bunch(i)+r*prn(i)
          enddo
        else
          do i=1,nbunch
            bunch(i)=r*prn(i)
          enddo
        endif
      elseif(abbrev(word,'A_DD','_'))then
        ipoint=next
        add=.true.
        go to 1
      else
        n=int(getva(exist))
        if(.not. exist)then
          call termes(lfno,
     1          '?Missing bunch number for DBUNCH.',' ')
          err=.true.
          return
        endif
        if(n .le. 0 .or. n .gt. nbunch)then
          call termes(lfno,
     1          '?Bunch number out of range.',' ')
          err=.true.
          return
        endif
        r=getva(exist)
        if(.not. exist)then
          call termes(lfno,
     1          '?Missing deviation for DBUNCH.',' ')
          err=.true.
          return
        endif
        if(add)then
          bunch(n)=bunch(n)+r
        else
          bunch(n)=r
        endif
      endif
      err=.false.
      return
      end
