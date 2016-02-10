      logical function mhogan(l1,l2,l)
      if(l1.lt.l2) then
        mhogan=.false.
        if(l.gt.l1) then
          if(l.lt.l2) then
            mhogan=.true.
          endif
        endif
      else
        mhogan=.false.
        if(l.gt.l1) then
          mhogan=.true.
        elseif(l.lt.l2) then
          mhogan=.true.
        endif
      endif
      return
      end
