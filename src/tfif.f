      subroutine tfif(word,iflevel,exist)
      use geto
      implicit none
      integer*4 ,intent(inout):: iflevel
      character*(*) ,intent(inout):: word
      character ch
      character*13 delim
      integer*4 next,ifl
      real*8 getva
      logical*4 f
      logical*4 ,intent(out):: exist
      data delim/'[*/)}]><=.&^|'/
      exist=.false.
 2    if(word == 'IF')then
        iflevel=iflevel+1  
 1      f=getva(exist) .ne. 0.d0
c        write(*,*)'tfif-1 ',f,exist
        if(.not. exist)then
          call termes('?Missing condition for IF',' ')
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
              if(word == 'ELSE')then
                if(iflevel .gt. ifl)then
                  word=' '
                endif
              elseif(word == 'IF')then
                iflevel=iflevel+1
              elseif(word == 'ENDIF')then
                iflevel=iflevel-1
              elseif(word == 'ELSEIF')then
                if(iflevel == ifl)then
                  word=' '
                  go to 1
                endif
              endif
            endif
c            write(*,*)'tfif-if  ',word(1:len_trim(word)),' ',ch,' ',iflevel
          enddo
          exist=.true.
        endif
      elseif(word == 'ELSE')then
        if(iflevel .le. 0)then
          call termes('?ELSE without IF',' ')
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
            if(word == 'IF')then
              iflevel=iflevel+1
            elseif(word == 'ENDIF')then
              iflevel=iflevel-1
            endif
          endif
        enddo
        exist=.true.
      elseif(word == 'ELSEIF')then
        if(iflevel .le. 0)then
          call termes('?ELSEIF without IF',' ')
          call skipline
          exist=.true.
          iflevel=0
          return
        endif
        word='ELSE'
        go to 2
      elseif(word == 'ENDIF')then
        if(iflevel .le. 0)then
          call termes('?ENDIF without IF',' ')
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
