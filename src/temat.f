      logical*4 function temat(i,name,word)
      use tfstk
      use ffs
      use ffs_pointer
      use tffitcode
      implicit none
      integer*4 i
      character*(*) name,word
      integer*4 lenw,ifany1
      logical*4 tmatchl
      character*(MAXPNAME+8) name1
      integer*4 lw,l
c
      lw=lenw(word)
      call elname(i,name)
      l=lenw(name)
c      write(*,*)'temat ',i,latt(i),ilist(2,latt(i)),name(1:l)
      if(ifany1(word,lw,'*{|%',1) .eq. 0)then
        if(name(1:l) .eq. word(1:lw))then
c     Exact match with standard element name(elname)
          temat=.true.
        elseif(i .eq. nlat)then
c     No more test for nlat($$$) 
          temat=.false.
        elseif(mult(i) .eq. 0)then
c     Test implicit match with singlet element
          call elnameK(i,name1)
          temat=name1(1:lenw(name1)) .eq. word(1:lw)
        elseif(nelvx(ilist(i,ifele1))%klp .eq. i)then
c     Test implicit match with head of multiple elements
          temat=pname(idelc(i))(1:lpnamec(i)) .eq. word(1:lw)
        else
          temat=.false.
        endif
      else
        if(tmatchl(name,l,word,lw))then
c     Pattern match with standard element name(elname)
          temat=.true.
        elseif(i .eq. nlat)then
c     No more test for nlat($$$) 
          temat=.false.
        elseif(mult(i) .eq. 0)then
c     Test implicit pattern match with singlet element
          call elnameK(i,name1)
          temat=tmatchl(name1,lenw(name1),word,lw)
        elseif(nelvx(ilist(i,ifele1))%klp .eq. i)then
c     Test implicit pattern match with head of multiple elements
          temat=tmatchl(pnamec(i),lpnamec(i),word,lw)
        else
          temat=.false.
        endif
      endif
      return
      end
