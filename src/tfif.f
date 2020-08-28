      subroutine tfif(word,iflevel,lfno,exist)
      implicit none
      integer*4 ,intent(inout):: iflevel
      integer*4 ,intent(in):: lfno
      character*(*) ,intent(inout):: word
      character ch
      character*13 delim
      integer*4 next,ifl
      real*8 getva
      logical*4 f
      logical*4 ,intent(out):: exist
      data delim/'[*/)}]><=.&^|'/
      exist=.false.
 2    if(word .eq. 'IF')then
        iflevel=iflevel+1  
 1      f=getva(exist) .ne. 0.d0
c        write(*,*)'tfif-1 ',f,exist
        if(.not. exist)then
          call termes(lfno,'?Missing condition for IF',' ')
          call skipline
          iflevel=0
          return
        endif
        if(f)then
          exist=.true.
        else
          ifl=iflevel
          do while(word .ne. 'ELSE' .and. iflevel .ge. ifl)
            call tfgettok(word)
            call peekwd(ch,next)
            if(index(delim,ch) .le. 0)then
              if(word .eq. 'ELSE')then
                if(iflevel .gt. ifl)then
                  word=' '
                endif
              elseif(word .eq. 'IF')then
                iflevel=iflevel+1
              elseif(word .eq. 'ENDIF')then
                iflevel=iflevel-1
              elseif(word .eq. 'ELSEIF')then
                if(iflevel .eq. ifl)then
                  word=' '
                  go to 1
                endif
              endif
            endif
c            write(*,*)'tfif-if  ',word(1:len_trim(word)),' ',ch,' ',iflevel
          enddo
          exist=.true.
        endif
      elseif(word .eq. 'ELSE')then
        if(iflevel .le. 0)then
          call termes(lfno,'?ELSE without IF',' ')
          call skipline
          iflevel=0
          exist=.true.
          return
        endif
        ifl=iflevel
        do while(iflevel .ge. ifl)
          call tfgettok(word)
          call peekwd(ch,next)
          if(index(delim,ch) .le. 0)then
c          write(*,*)'tfif-else ',word(1:len_trim(word)),iflevel
            if(word .eq. 'IF')then
              iflevel=iflevel+1
            elseif(word .eq. 'ENDIF')then
              iflevel=iflevel-1
            endif
          endif
        enddo
        exist=.true.
      elseif(word .eq. 'ELSEIF')then
        if(iflevel .le. 0)then
          call termes(lfno,'?ELSEIF without IF',' ')
          call skipline
          exist=.true.
          iflevel=0
          return
        endif
        word='ELSE'
        go to 2
      elseif(word .eq. 'ENDIF')then
        if(iflevel .le. 0)then
          call termes(lfno,'?ENDIF without IF',' ')
          call skipline
          exist=.true.
          iflevel=0
          return
        endif
        iflevel=iflevel-1
        exist=.true.
      endif
      return
      end
