      subroutine tfrej(word,nlist,nfc,mfpnta,mfpnta1,
     1    kfit,ndp,ifitp,ifitp1,exist)
      use tfstk
      use ffs
      use tffitcode
      implicit none
      integer*4 i,nfc,mfpnta,mfpnta1,l,lene,j
      character*8 nlist(mfit1),name1
      character*(*) ,intent(out):: word
      integer*4 kfit(*),ndp(*),ifitp(*),ifitp1(*)
      logical*4 exist,tmatch
      exist=.true.
      do while(exist)
        call getwdl(word)
        if(word .eq. 'TOTAL' .or. word .eq. 'TOTALFIT')then
          do i=1,nfc
            ndp(i)=0
          enddo
          if(word .eq. 'TOTAL')then
            call tfinitcalc
          endif
          return
        endif
        exist=.false.
        do i=1,mfit1
          l=lene(nlist(i))
          name1=nlist(i)
          name1(l+1:l+1)='M'
          if(tmatch(nlist(i),word) .or. tmatch(name1,word))then
            exist=.true.
            do j=1,nfc
              if(ifitp(j) .eq. mfpnta .and. ifitp1(j) .eq. mfpnta1)then
                if(kfit(j) .eq. i)then
                  ndp(j)=0
                  exit
                endif
              endif
            enddo
          endif
        enddo
      enddo
      return
      end
