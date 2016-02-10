      logical*4 function tmatch(name,patt)
      implicit none
      integer*4 lenw
      character*(*) name,patt
      logical*4 tmatchl
      tmatch=tmatchl(name,lenw(name),patt,lenw(patt))
      return
      end

      recursive logical*4 function tmatchl(name,ln,patt,lp)
     $     result(lx)
      use tfstk
      implicit none
      integer*4 lp,ln,notchar1,ia,regexp,ir,ia1,ia2
      character*(*) name,patt
      character en,ep
      logical*4 tmatchl1,tmatchl3
      integer*8 iafwild
      data iafwild/0/
      if(lp .le. 0)then
        lx=ln .le. 0
        return
      endif
      if(iafwild .eq. 0)then
        iafwild=ktfsymbolz('$WildCardID',11)
      endif
      if(rlist(iafwild-4) .ne. 0.d0)then
        en=name(ln+1:ln+1)
        name(ln+1:ln+1)=char(0)
        ep=patt(lp+1:lp+1)
        patt(lp+1:lp+1)=char(0)
        ir=regexp(name,patt)
        name(ln+1:ln+1)=en
        patt(lp+1:lp+1)=ep
        lx=ir .eq. 1
        return
      endif
      ia=index(patt(1:lp),'|')
      if(ia .gt. 0)then
        lx=tmatchl(name,ln,patt(1:ia-1),ia-1) .or.
     $       tmatchl(name,ln,patt(ia+1:lp),lp-ia)
        return
      else
        ia1=index(patt(1:lp),'{')
        if(ia1 .le. 0)then
          ia2=index(patt(1:lp),'%')
          if(ia2 .le. 0)then
            lx=tmatchl3(name,ln,patt,lp)
            return
          endif
        endif
        if(ln .le. 0)then
          lx=notchar1(patt(1:lp),lp,'*',1) .eq. 0
        else
          lx=tmatchl1(name,ln,patt,lp)
        endif
      endif
      return
      end

      recursive logical*4 function tmatchl1(name,ln,patt,lp)
     $     result(lx)
      implicit none
      integer*4 lp,ln,notchar1,ic,notany1,ifany1,ic1,ic2,id1
      character*(*) name,patt
      logical*4 tmatchl2
      lx=.false.
      id1=index(patt(1:lp),'{')
      if(id1 .le. 0)then
        lx=tmatchl2(name,ln,patt,lp)
        return
      elseif(id1 .eq. 1)then
        ic=index(patt(1:lp),'}')
        if(ic .eq. 0)then
          ic=lp+1
        endif
        if(patt(2:2) .eq. '^')then
          if(index(patt(3:ic-1),name(1:1)) .gt. 0)then
            return
          endif
        elseif(index(patt(2:ic-1),name(1:1)) .eq. 0)then
          return
        endif
        if(lp .le. ic)then
          lx=ln .eq. 1
        elseif(ln .eq. 1)then
          lx=notchar1(patt,lp,'*',ic+1) .eq. 0
        else
          lx=tmatchl1(name(2:ln),ln-1,patt(ic+1:lp),lp-ic)
        endif
        return
      endif
      ic1=notany1(patt,lp,'*%{}',1)
      if(ic1 .eq. 1)then
        ic2=max(ifany1(patt,lp,'*%{',ic1+1)-1,ic1)
        if(ln .lt. ic2)then
          return
        elseif(name(1:ic2) .ne. patt(1:ic2))then
          return
        elseif(ic2 .eq. lp)then
          lx=ln .eq. ic2
        elseif(ln .eq. ic2)then
          lx=notchar1(patt,lp,'*',ic2+1) .eq. 0
        else
          lx=tmatchl1(name(ic2+1:ln),ln-ic2,
     $         patt(ic2+1:lp),lp-ic2)
        endif
        return
      endif
      if(ic1 .ne. 0)then
        if(id1+1 .ne. ic1)then
          ic2=max(ifany1(patt,lp,'*%{',ic1+1)-1,ic1)
          if(index(name(1:ln),patt(ic1:ic2)) .le. 0)then
            return
          endif
        endif
      endif
      if(patt(1:1) .eq. '*')then
        if(lp .eq. 1)then
          lx=.true.
        else
          lx=tmatchl1(name,ln,patt(2:lp),lp-1)
          if(.not. lx .and. ln .ne. 1)then
            lx=tmatchl1(name(2:ln),ln-1,patt(1:lp),lp)
          endif
        endif
        return
      elseif(patt(1:1) .eq. '%')then
        if(lp .eq. 1)then
          lx=ln .eq. 1
        elseif(ln .eq. 1)then
          lx=notchar1(patt,lp,'*',2) .eq. 0
        else
          lx=tmatchl1(name(2:ln),ln-1,patt(2:lp),lp-1)
        endif
      elseif(name(1:1) .eq. '}')then
        if(lp .eq. 1)then
          lx=ln .eq. 1
        elseif(ln .eq. 1)then
          lx=notchar1(patt,lp,'*',2) .eq. 0
        else
          lx=tmatchl1(name(2:ln),ln-1,patt(2:lp),lp-1)
        endif
      endif
      return
      end

      recursive logical*4 function tmatchl2(name,ln,patt,lp)
     $     result(lx)
      implicit none
      integer*4 lp,ln,notchar1,notany1,ifany1,ic1,ic2
      character*(*) name,patt
      lx=.false.
      ic1=notany1(patt,lp,'*%',1)
      if(ic1 .eq. 1)then
        ic2=max(ifany1(patt,lp,'*%',ic1+1)-1,ic1)
        if(ln .lt. ic2)then
          return
        elseif(name(1:ic2) .ne. patt(1:ic2))then
          return
        elseif(ic2 .eq. lp)then
          lx=ln .eq. ic2
        elseif(ln .eq. ic2)then
          lx=notchar1(patt,lp,'*',ic2+1) .eq. 0
        else
          lx=tmatchl2(name(ic2+1:ln),ln-ic2,
     $         patt(ic2+1:lp),lp-ic2)
        endif
        return
      endif
      if(ic1 .ne. 0)then
        ic2=max(ifany1(patt,lp,'*%',ic1+1)-1,ic1)
        if(index(name(1:ln),patt(ic1:ic2)) .le. 0)then
          return
        endif
      endif
      if(patt(1:1) .eq. '*')then
        if(lp .eq. 1)then
          lx=.true.
        else
          lx=tmatchl2(name,ln,patt(2:lp),lp-1)
          if(.not. lx .and. ln .ne. 1)then
            lx=tmatchl2(name(2:ln),ln-1,patt(1:lp),lp)
          endif
        endif
        return
      else
        if(lp .eq. 1)then
          lx=ln .eq. 1
        elseif(ln .eq. 1)then
          lx=notchar1(patt(2:lp),lp-1,'*',1) .eq. 0
        else
          lx=tmatchl2(name(2:ln),ln-1,patt(2:lp),lp-1)
        endif
      endif
      return
      end

      recursive logical*4 function tmatchl3(name,ln,patt,lp)
     $     result(lx)
      implicit none
      integer*4 lp,ln,notchar1,ic1,ic2,id1,id2
      character*(*) name,patt
      lx=.false.
      ic1=notchar1(patt,lp,'*',1)
      if(ic1 .eq. 1)then
        ic2=max(0,index(patt(ic1+1:lp),'*')-1)+ic1
        if(ln .lt. ic2)then
          return
        elseif(name(1:ic2) .ne. patt(1:ic2))then
          return
        elseif(ic2 .eq. lp)then
          lx=ln .eq. ic2
        elseif(ln .eq. ic2)then
          lx=notchar1(patt,lp,'*',ic2+1) .eq. 0
        else
          lx=tmatchl3(name(ic2+1:ln),ln-ic2,
     $         patt(ic2+1:lp),lp-ic2)
        endif
        return
      endif
      if(ic1 .ne. 0)then
        ic2=max(0,index(patt(ic1+1:lp),'*')-1)+ic1
        if(ic2 .eq. lp)then
          if(ln .ge. ic2-ic1+1)then
            lx=name(ln-ic2+ic1:ln) .eq. patt(ic1:ic2)
          endif
          return
        endif
        id1=index(name(1:ln),patt(ic1:ic2))
        if(id1 .le. 0)then
          return
        endif
        id2=id1+ic2-ic1
        if(ln .eq. id2)then
          lx=notchar1(patt,lp,'*',ic2+1) .eq. 0
        else
          lx=tmatchl3(name(id2+1:ln),ln-id2,
     $         patt(ic2+1:lp),lp-ic2) .or.
     $         tmatchl3(name(id2+1:ln),ln-id2,
     $         patt,lp)
        endif
        return
      endif
      if(lp .eq. 1)then
        lx=.true.
      else
        lx=tmatchl3(name,ln,patt(2:lp),lp-1)
        if(.not. lx .and. ln .ne. 1)then
          lx=tmatchl3(name(2:ln),ln-1,patt(1:lp),lp)
        endif
      endif
      return
      end
