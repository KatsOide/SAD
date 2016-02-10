      subroutine tfchgv(latt,ival,klp,lfno)
      use tfstk
      use ffs
      use tffitcode
      implicit none
      integer*4 lfno,i,k,j,next1
      integer*4 latt(2,nlat),ival(nele),klp(nele)
      character*(MAXPNAME) word
      character*8 key,ki,tfkwrd
      logical*4 tmatch,apply,exist
      call getwdl(key)
      apply=.false.
 1    call peekwd(word,next1)
      exist=.false.
      LOOP_I: do i=1,nele
        k=latt(1,klp(i))
        if(tmatch(pname(k),word))then
          if(.not. exist)then
            call cssetp(next1)
            exist=.true.
          endif
          j=0
          ki='-'
          do while(ki .ne. ' ')
            j=j+1
            ki=tfkwrd(idtype(k),j)
            if(ki .eq. key)then
              apply=.true.
              ival(i)=j
              cycle LOOP_I
            endif
          enddo
        endif
      enddo LOOP_I
      apply=apply .or. exist
      if(exist)then
        go to 1
      endif
      if(.not. apply)then
        call termes(lfno,
     $       'No element for VAR_Y with keyword: ',key)
        call cssetp(next1)
      endif
      return
      end
