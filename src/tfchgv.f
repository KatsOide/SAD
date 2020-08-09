      subroutine tfchgv(lfno)
      use tfstk
      use ffs
      use ffs_pointer
      use tffitcode
      use tfcsi,only:ipoint
      implicit none
      integer*4 ,intent(in):: lfno
      integer*4 i,k,j,next1
      character*(MAXPNAME) word
      character*8 key,ki,tfkwrd,tfkwrd1
      logical*4 tmatch,apply,exist
      call getwdl(key)
      apply=.false.
 1    call peekwd(word,next1)
      exist=.false.
      LOOP_I: do i=1,nele
        k=idelc(nelvx(i)%klp)
        if(tmatch(pname(k),word))then
          if(.not. exist)then
            ipoint=next1
            exist=.true.
          endif
          j=0
          ki='-'
          do while(ki .ne. ' ')
            j=j+1
            ki=tfkwrd(idtype(k),j)
            if(ki .eq. key)then
              apply=.true.
              nelvx(i)%ival=j
              cycle LOOP_I
            endif
            ki=tfkwrd1(idtype(k),j)
            if(ki .eq. key)then
              apply=.true.
              nelvx(i)%ival=j
              cycle LOOP_I
            endif
          enddo
        endif
      enddo LOOP_I
      apply=apply .or. exist
      if(exist)then
        go to 1
      endif
      if(apply)then
        evarini=.true.
      else
        call termes(lfno,
     $       'No element for VAR_Y with keyword: ',key)
        ipoint=next1
      endif
      return
      end
