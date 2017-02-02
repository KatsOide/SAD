      subroutine pwmatq(latt,twiss,gammab,pos,a,b,imon,nmon,nq,io)
      use tfstk
      use ffs
      use tffitcode
      implicit real*8(a-h,o-z)
      parameter (halfpi=pi*0.5d0)
      logical avecmo
      integer*8 latt(nlat)
      dimension twiss(nlat,-ndim:ndim,ntwissfun),gammab(nlat),pos(nlat)
      dimension imon(nmona,4)
      dimension a(nmon,nq,2),b(nq)
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
      iq=0
c     end   initialize for preventing compiler warning
      psix=twiss(nlat,0,3)-twiss(1,0,3)
      psiy=twiss(nlat,0,6)-twiss(1,0,6)
      lq=0
      itemp=italoc(nstra)
      do 20 l=1,nlat-1
        if(idtype(ilist(2,latt(l))).eq.icquad) then
          avecmo=.false.
          t=rlist(idval(ilist(2,latt(l)))+4)
c         .... reject skew quads ....
          if(abs(t-pi/4d0).lt.0.01 .or. abs(t+pi/4d0).lt.0.01) goto 20
          lq=lq+1
          do 21 i=1,nmon
            if(imon(imon(i,2),4).eq.l) then
              avecmo=.true.
              iq=i
            endif
   21     continue
          b(lq)=0.5d0*(pos(l)+pos(l+1))
          ilist(1,itemp)=l
          ilist(mod(nstra,2)+1,itemp+nstra/2)=1
          call mcrmat(latt,twiss,gammab,0,psix,psiy,a(1,lq,1),nmon,
     1               .false.,.false.,ilist(1,itemp),1,imon,nmon,'X')
          do 22 i=1,nmon
            a(i,lq,1)=a(i,lq,1)*rlist(idval(ilist(2,latt(l)))+2)
   22     continue
          dx=rlist(idval(ilist(2,latt(l)))+5)
          rlist(idval(ilist(2,latt(l)))+5)=t+halfpi
          call mcrmat(latt,twiss,gammab,0,psix,psiy,a(1,lq,2),nmon,
     1                .false.,.false.,ilist(1,itemp),1,imon,nmon,'Y')
          rlist(idval(ilist(2,latt(l)))+5)=dx
          do 24 i=1,nmon
            a(i,lq,2)=-a(i,lq,2) * rlist(idval(ilist(2,latt(l)))+2)
   24     continue
          if(avecmo) then
            a(iq,lq,1)=a(iq,lq,1)-1d0
            a(iq,lq,2)=a(iq,lq,2)-1d0
          endif
        endif
   20 continue
      call tfree(int8(itemp))
      write(io,err=99) nmon,nq
      write(*,*) nmon,nq
      write(io,err=99) ((a(i,j,1),i=1,nmon),j=1,nq)
      write(io,err=99) ((a(i,j,2),i=1,nmon),j=1,nq)
      write(io,err=99) (b(i),i=1,nq)
   99 return
      end
