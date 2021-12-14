      integer*4 function lrflct(idx)
      integer*4 idx
c.....this function treats (-1)*<element>.
c.....function porduce 1*<-element> ic core.
c.....and returns a value of index for <-element>.
c.....following line is temporaly one.
      lrflct=idx
      return
      end
