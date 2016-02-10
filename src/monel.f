      logical function monel(imon,nmon,n,m)
      implicit real*8 (a-h,o-z)
      dimension imon(nmona,4)
      logical errflg
      common /codcor/eptsol,errflg(4),errval(4),dpshft,optiv(18),
     1               ibckup,nmona,nmonact,nstra,nstract,
     1               itmon,itemon,itstr,itestr
c          mcepst msolvg mcmon mcnrmc mdpmax palgn pcbak pundo prkick
c eptsol     o      o                          o
c errflg                   o
c errval                   o     o             o
c dpshft                         o      o
c ibckup                                             o      o     o
c optiv                                              o
c nmona                    o     o             o
c nmonact
c nstra                                        o                  o
c nstract
c itmon
c itemon
c itstr
c itestr
c          pbump  corinit mcstr mccor mcrcod mcrmda monel pwrite mrecal
c eptsol     o      o             o
c errflg            o
c errval            o
c dpshft            o
c ibckup     o      o
c optiv
c nmona      o                    o     o      o      o     o      o
c nmonact
c nstra                     o     o                         o      o
c nstract
c itmon             o
c itemon            o
c itstr             o       o
c itestr            o       o
c          pkill pstati mstore ptrim petcod mclear twsdrw pbumps monact
c eptsol
c errflg
c errval
c dpshft
c ibckup
c optiv
c nmona      o     o      o      o      o            o
c nmonact                                                          o
c nstra      o     o      o                   o              o
c nstract
c itmon                                                            o
c itemon
c itstr
c itestr
c          mcrmat mbmp msolb pcrmat ptol pcset mstack mweght
c eptsol
c errflg
c errval
c dpshft
c ibckup
c optiv
c nmona                                          o      o
c nmonact
c nstra      o     o     o     o     o     o     o
c nstract
c itmon
c itemon
c itstr
c itestr
c     begin initialize for preventing compiler warning
      mm=0
c     begin initialize for preventing compiler warning
      monel=.false.
      if(nmon.eq.0)return
      m=min(max(1,m),nmon)
      if(imon(imon(m,2),1).lt.n) then
        do 10 i=m+1,nmon
          j=imon(i,2)
          if(imon(j,1).eq.n) then
            monel=.true.
            mm=i
            goto 21
          elseif(imon(j,1).gt.n) then
            return
          endif
   10   continue
      else
        do 20 i=m,1,-1
          j=imon(i,2)
          if(imon(j,1).eq.n) then
            monel=.true.
            mm=i
            goto 21
          elseif(imon(j,1).lt.n) then
            return
          endif
   20   continue
      endif
   21 m=mm
      return
      end
