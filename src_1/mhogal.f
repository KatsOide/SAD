      logical function mhogal(l1,l2,l)
      if(l1.lt.l2) then
        mhogal=.false.
        if(l.ge.l1) then
          if(l.lt.l2) then
            mhogal=.true.
          endif
        endif
      elseif(l1.gt.l2) then
        mhogal=.false.
        if(l.ge.l1) then
          mhogal=.true.
        elseif(l.lt.l2) then
          mhogal=.true.
        endif
      elseif(l1.eq.l) then
        mhogal=.true.
      else
        mhogal=.false.
      endif
      return
      end
