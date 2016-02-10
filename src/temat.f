      logical*4 function temat(latt,i,mult,name,word)
      use tfstk
      use ffs
      use tffitcode
      implicit none
      integer*4 latt(2,nlat),i,mult(nlat)
      character*(*) name,word
      integer*4 lenw,ifany1,lpname
      logical*4 tmatchl
      character*(MAXPNAME+8) name1
      integer*4 lw,l
c
      lw=lenw(word)
      call elname1(latt,i,mult,name,.false.)
      l=lenw(name)
      if(ifany1(word,lw,'*{|%',1) .eq. 0)then
        if(name(1:l) .eq. word(1:lw))then
c     Exact match with standard element name(elname)
          temat=.true.
        elseif(i .eq. nlat)then
c     No more test for nlat($$$) 
          temat=.false.
        elseif(mult(i) .eq. 0)then
c     Test implicit match with singlet element
          call elname1(latt,i,mult,name1,.true.)
          temat=name1(1:lenw(name1)) .eq. word(1:lw)
        elseif(ilist(ilist(i,ifele1),ifklp) .eq. i)then
c     Test implicit match with head of multiple elements
          temat=pname(latt(1,i))(1:lpname(latt(1,i))) .eq. word(1:lw)
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
          call elname1(latt,i,mult,name1,.true.)
          temat=tmatchl(name1,lenw(name1),word,lw)
        elseif(ilist(ilist(i,ifele1),ifklp) .eq. i)then
c     Test implicit pattern match with head of multiple elements
          temat=tmatchl(pname(latt(1,i)),lpname(latt(1,i)),word,lw)
        else
          temat=.false.
        endif
      endif
      return
      end
