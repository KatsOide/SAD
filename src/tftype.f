      subroutine tftype(lfno,word)
      use tfstk
      use ffs_pointer
      use tffitcode
      use mackw
      implicit none
      integer*4 ntyp
      parameter (ntyp=18)
      integer*8 kavl
      integer*4 lfno,k1,ltyp(ntyp)
      integer*4 j,kx,lt,kp,notchar,ifany,lpw,itfgetrecl,
     $     nl,kkk,irtc
      character*(*) word
      character*(MAXPNAME) name1,tfkwrd
      character*80 patt
      logical*4 exist,exist1,tmatch,start,mulc
      data ltyp /icDRFT,
     $     icBEND,icQUAD,icSEXT,icOCTU,icDECA,icDODECA,icMULT,icSOL,
     $     icCAVI,icTCAV,
     $     icMAP,icBEAM,icINS,icCOORD,icAPRT,icMONI,icMARK/
      exist1=.false.
      lpw=min(131,itfgetrecl())
      name1=' '
1     call getwdl(word)
      if(word .eq. 'ALL' .or. word .eq. ' ' .and. .not. exist1)then
        patt='*'
      elseif(word .eq. ' ')then
        patt=' '
        exist=.false.
        go to 21
      elseif(notchar(word,'*',1) .eq. 0)then
        patt='*'
        call getwdl(word)
      else
        patt=word
      endif
      exist=.false.
      if(patt .eq. ' ')then
        go to 21
      endif
      call tfgetlineps(patt,len_trim(patt),nl,kavl,1,irtc)
      if(irtc .ne. 0 .or. nl .le. 0)then
        go to 21
      endif
      do j=1,ntyp
        lt=ltyp(j)
        mulc=kytbl(kwMAX,lt) .gt. kytbl(kwMAX, icDRFT)
        start=.true.
        do kkk=1,nl
          k1=int(rlist(kavl+kkk))
          kx=klp(k1)
c     Note: Skip no-head multiple elements
c     *     klp(iele1(kx)) == kx if singlet or head of multipole elements
c          if(klp(iele1(kx)) .ne. kx)cycle
          kp=idelc(kx)
          if(tmatch(pname(kp),patt))then
            exist=.true.
            if(idtype(kp) .eq. lt)then
              if(start)then
                write(lfno,*)';'
                call twbuf(' ',lfno,1,lpw,0,0)
                call twbuf(tfkwrd(lt,0),lfno,1,lpw,7,1)
                start=.false.
              endif
              call tftyp1(iele1(kx),kx,latt(kx),kp,lt,lfno,lpw)
              if(mulc)then
                call twbuf(' ',lfno,10,lpw,0,-1)
              endif
            endif
          endif
        enddo
        if(.not. start .and. .not. mulc)then
          call twbuf(' ',lfno,10,lpw,0,-1)
        endif
      enddo
21    exist1=exist1 .or. exist
c      if(.not. exist1)then
c        patt='*'
c        go to 2
c      endif
      if(exist .and. notchar(patt,'*',1) .ne. 0)then
        go to 1
      else
        exist=ifany(patt,'*%{|',1) .ne. 0
        write(lfno,*)';'
        return
      endif
      end
