c$Header: /SAD/cvsroot/oldsad/src/tfeval.f,v 1.156.2.3 2012/07/20 03:55:45 oide Exp $
      subroutine tfeval(string,l,ist1,istop,kx,re,irtc)
      implicit none
      include 'inc/TFCBK.inc'
      include 'inc/TFSTK.inc'
      include 'inc/TFCODE.inc'
      logical*4 re,T,F
      parameter (T=.true.,F=.false.)
      character*1023 string
      integer*8 kx,ktfmakelist,kt,ktpfaloc,kax,ka
      integer*4 istart,istop,irtc,isp0,ist10,iop1,
     $     iprior(0:mtfnopc),icsmrk,i,ist1,ishash,
     $     l,ifchar,mopc,itgetfpe,m1,
     $     itfmessage,level1,ist2,icslfni,icslfno,irt
      logical*4 nullfirst(0:mtfnopc),lastfirst(0:mtfnopc)
      real*8 vx
c      include 'DEBUG.inc'
      logical*4 tfreadevalbuf,eol
      data iprior/
     $     9999,
     $     10,  20,  50,  50,  40,  40,  15,  15,  100, 100,
     $     100, 100, 100, 100, 120, 120, 150, 160, 170, 80,
     $     6,   3000,9999,3000,250, 250, 7000,9999,8000,9000,
     $     1000,220, 180, 190, 190, 200, 200, 250, 250, 900,
     $     4,   3,   3,   3,   10,  175, 9,   9,   9,   172,
     $     172, 130, 210, 210, 210, 210, 7,   7,   6,   5,
     $     2,   240, 9999,9999,9999,9999/
c          null
c          m    i    +    -    *    /    v    ^    e    n    
c          >    <    g    l    E    N    ~    &    o    c
c          [    ]    {    }    s    =    C    (    )    ,
c          ;    &    :    r    d    RepA RepR u    U    USet
c          ?    f    #    ##   .    |    M    MA   A    rept 
c          repn ineq AT   SF   TB   DB   Inc  Dec  Part @
c          msgn TagS (*   *)   hold z
      data nullfirst/
     $     T,
     $     T,   T,   F,   F,   F,   F,   F,   F,   F,   F,   
     $     F,   F,   F,   F,   F,   F,   T,   F,   F,   F,   
     $     T,   T,   T,   T,   F,   F,   F,   T,   T,   T,
     $     T,   F,   F,   F,   F,   F,   F,   F,   F,   F,   
     $     T,   T,   T,   T,   F,   F,   F,   F,   F,   F,   
     $     F,   F,   F,   F,   F,   F,   T,   T,   F,   F,
     $     F,   F,   F,   F,   F,   F/
      data lastfirst/
     $     F,
     $     F,   F,   F,   F,   F,   F,   F,   T,   F,   F,   
     $     F,   F,   F,   F,   F,   F,   F,   F,   F,   F,   
     $     F,   F,   F,   F,   T,   T,   F,   F,   F,   F,
     $     F,   F,   F,   T,   T,   F,   F,   F,   F,   F,   
     $     F,   F,   F,   F,   F,   F,   T,   T,   T,   F,   
     $     F,   F,   F,   F,   F,   F,   F,   F,   F,   F,
     $     F,   F,   F,   F,   F,   F/
c     begin initialize for preventing compiler warning
      mopc=0
c     end   initialize for preventing compiler warning
      isp0=isp+1
      isp=isp0
      ierrorf=0
      ierrorprint=0
      jtastk(1,ivstkoffset+isp)=mtfnull
      ktastk(isp)=ktfoper
      istart=ist1
      istop=istart
      ishash=-1
      vx=0.d0
      if(.not. re)then
        jtastk(1,ivstkoffset+isp)=mtfleftparen
        isp=isp+1
        jtastk(1,ivstkoffset+isp)=mtfnull
        ktastk(isp)=ktfoper+mtfnull
      endif
      levele=levele+1
      if(levele .ge. maxlevele-10)then
        go to 8120
      endif
      irtc=0
      eol=.false.
 1    continue
c      write(*,*)'tfeval ',istart,l,string(istart:l)
      call tfetok(string(istart:l),istop,kx,itfcontext,irt)
      istop=istop+istart-1
c        call tfreecheck1('tfeval-0',1,0,0.d0,irtc)
        write(*,*)'tfeval-0 ',irt,kx,istart,istop,
     $       string(istart:istop)
c        if(irtc .ne. 0)then
c          rlist(7)=0.d0
c        endif
      if(irt .ge. 0)then
        go to 2400
      endif
 2    go to (2100,2200,2300,8040),-irt
c            oper eol  cmnt clq
c
      write(*,*)'tfeval ',irt
      stop
 2100 mopc=ktfaddr(kx)
 2101 if(ktastk(isp) .eq. ktfoper+mtfnull)then
        if(isp .gt. isp0)then
          m1=jtastk(1,ivstkoffset+isp-1)
          if(m1 .eq. mtfslot .or. m1 .eq. mtfslotseq)then
          elseif(m1 .eq. mtfrightbra .and. mopc .ne. mtfrightbra)then
            do i=isp-2,isp0,-1
              if(jtastk(1,ivstkoffset+i) .eq. mtfpart)then
                isp=isp-1
                kax=ktfmakelist(i-1)
                klist(kax)=ktfoper+mtfpart
                ktastk(i)=ktflist+kax
                jtastk(1,ivstkoffset+i)=mtfnull
                isp=i
                go to 910
              endif
            enddo
            go to 8050
          elseif(mopc .eq. mtfminus)then
            if(m1 .eq. mtfpower .or. m1 .eq. mtfrevpower)then
              mopc=mtfneg
            else
              mopc=mtftimes
            endif
            rtastk(isp)=-1.d0
          elseif(mopc .eq. mtfplus)then
            go to 1010
          elseif(m1 .eq. mtfdot .and. mopc .eq. mtftimes)then
            kx=ktfoper+mtfnull
            istop=ist1
            irtc=0
            go to 8910
          elseif(m1 .eq. mtffun)then
            isp=isp-1
            ktastk(isp)=ktflist+ktpfaloc(ktastk(isp))
            go to 910
          elseif(mopc .eq. mtfpattest)then
            mopc=mtfflag
          elseif(.not. nullfirst(mopc))then
            if(m1 .eq. mtfpower .and. mopc .eq. mtfpower)then
              go to 8600
            endif
            go to 8700
          endif
        else
          if(mopc .eq. mtfcomp .and.
     $         jtastk(1,ivstkoffset+isp) .eq. mtfnull)then
            kx=ktfoper+mtfnull
            irtc=-1
            if(re)then
              istop=istop-1
            endif
            go to 9000
          elseif(mopc .eq. mtfminus)then
            mopc=mtftimes
            rtastk(isp)=-1.d0
          elseif(mopc .eq. mtfplus)then
            go to 1010
          elseif(mopc .eq. mtfpattest)then
            mopc=mtfflag
          elseif(.not. nullfirst(mopc))then
            go to 8700
          endif
        endif
      else
        if(mopc .eq. mtfleftparen .or. mopc .eq. mtflist
     $       .or. mopc .eq. mtfslot
     $       .or. mopc .eq. mtfslotseq)then
          do i=isp0,isp-1
            if(jtastk(1,ivstkoffset+i) .eq. mtflist
     $           .or. jtastk(1,ivstkoffset+i) .eq. mtfleftbra
     $           .or. jtastk(1,ivstkoffset+i) .eq. mtfleftparen
     $           .or. jtastk(1,ivstkoffset+i) .eq. mtfpart)then
              jtastk(1,ivstkoffset+isp)=mtftimes
              call tfestk(isp0,iprior,lastfirst,irtc)
              if(irtc .ne. 0)then
                go to 8900
              endif
              isp=isp+1
              go to 910
            endif
          enddo
          jtastk(1,ivstkoffset+isp)=mtfnull
          istop=istart
          irt=-2
          go to 2
        endif
      endif
 910  if(mopc .gt. 0)then
        jtastk(1,ivstkoffset+isp)=mopc
        if(isp .gt. isp0)then
          call tfestk(isp0,iprior,lastfirst,irtc)
          if(irtc .ne. 0)then
            go to 8900
          endif
        endif
        if(mopc .eq. mtfcomma)then
          if(isp .eq. isp0)then
            go to 7000
          endif
        elseif(mopc .eq. mtfcomp)then
          if(isp .eq. isp0 .and. re)then
            istop=istop-1
            go to 7000
          endif
        elseif(mopc .eq. mtfminus)then
          jtastk(1,ivstkoffset+isp)=mtfplus
          isp=isp+1
          rtastk(isp)=-1.d0
          jtastk(1,ivstkoffset+isp)=mtftimes
        elseif(mopc .eq. mtfdiv)then
          jtastk(1,ivstkoffset+isp)=mtftimes
          isp=isp+1
          rtastk(isp)=-1.d0
          jtastk(1,ivstkoffset+isp)=mtfrevpower
        elseif(mopc .eq. mtfunset .or.
     $         mopc .eq. mtfrepeated .or.
     $         mopc .eq. mtfrepeatednull .or.
     $         (mopc .eq. mtfincrement .or.
     $         mopc .eq. mtfdecrement) .and.
     $         ktastk(isp) .ne. ktfoper+mtfnull)then
          ka=ktfmakelist(isp-1)
          klist(ka)=ktfoper+mopc
          ktastk(isp)=ktflist+ka
          jtastk(1,ivstkoffset+isp)=mtfnull
        elseif(jtastk(1,ivstkoffset+isp) .eq. mtfrightbrace
     $         .or. jtastk(1,ivstkoffset+isp) .eq. mtfrightparen)then
          go to 8050
        endif
        if(jtastk(1,ivstkoffset+isp) .ne. mtfnull)then
          isp=isp+1
          if(isp .gt. mstk)then
            go to 8110
          endif
          jtastk(1,ivstkoffset+isp)=mtfnull
          ktastk(isp)=ktfoper+mtfnull
        endif
        go to 1010
      else
        go to 8010
      endif
 2200 if(ktastk(isp) .eq. ktfoper+mtfnull)then
        if(isp .eq. isp0 .and.
     $       jtastk(1,ivstkoffset+isp) .eq. mtfnull)then
          irtc=-1
        elseif(isp .gt. isp0)then
          iop1=jtastk(1,ivstkoffset+isp-1)
          if(re)then
            if(string(istop-1:istop-1) .eq. char(10))then
              if(tfreadevalbuf(istart,istop,l,iop1))then
                eol=.true.
                go to 1
              endif
            endif
          elseif(iop1 .eq. mtfcomp)then
            go to 3
          elseif(iop1 .eq. mtffun)then
            isp=isp-1
            ktastk(isp)=ktflist+ktpfaloc(ktastk(isp))
            jtastk(1,ivstkoffset+isp)=mtfnull
            call tfdebugprint(ktastk(isp),'tfeval',3)
            go to 3
          endif
          go to 8700
        endif
        kx=ktastk(isp0)
        go to 9000
      endif
 3    if(isp .gt. isp0)then
        if(.not. re)then
          jtastk(1,ivstkoffset+isp)=mtfrightparen
          call tfestk(isp0,iprior,lastfirst,irtc)
          if(irtc .ne. 0)then
            go to 8900
          endif
          if(isp .le. isp0)then
            go to 7000
          endif
        endif            
        do i=isp0,isp
          if(jtastk(1,ivstkoffset+i) .eq. mtflist
     $         .or. jtastk(1,ivstkoffset+i) .eq. mtfleftbra
     $         .or. jtastk(1,ivstkoffset+i) .eq. mtfleftparen
     $         .or. jtastk(1,ivstkoffset+i) .eq. mtfpart)then
            if(re)then
              if(string(istop-1:istop-1) .eq. char(10))then
                if(tfreadevalbuf(istart,istop,l,
     $               int(jtastk(1,ivstkoffset+i))))then
                  jtastk(1,ivstkoffset+isp)=mtfnull
                  eol=.true.
                  go to 1
                endif
              endif
            endif
            go to 8700
          endif
        enddo
c     Sanity check for SAD stack
c     Stack MIGHT not have mtfright(brace|bra|paren),
c     because stack does not have mtfleft(brace|bra|paren)!!
        do i=isp,isp0,-1
          if(jtastk(1,ivstkoffset+i) .eq. mtfrightbrace
     $         .or. jtastk(1,ivstkoffset+i) .eq. mtfrightbra
     $         .or. jtastk(1,ivstkoffset+i) .eq. mtfrightparen)then
            mopc=1
            go to 8050
          endif
        enddo
        jtastk(1,ivstkoffset+isp)=mtfend
        call tfestk(isp0,iprior,lastfirst,irtc)
        if(irtc .ne. 0)then
          go to 8900
        endif
      endif
      go to 7000
 2300 if(re)then
        if(string(istop-1:istop-1) .eq. char(10))then
 3001     if(tfreadevalbuf(istart,istop,l,
     $         mtfleftcomment))then
            ist2=index(string(istart:l),'*)')
            if(ist2 .le. 0)then
              go to 3001
            endif
            istart=ist2+istart+1
            eol=.false.
            go to 1
          endif
        endif
      endif
      irtc=itfmessage(9999,'General::comment',
     $     '"'//string(ist1:min(istop-1,l))//'"')
      istop=ist1
      go to 8800
 2400 if(ktastk(isp) .ne. ktfoper+mtfnull)then
        do i=isp0,isp-1
          if(jtastk(1,ivstkoffset+i) .eq. mtflist
     $         .or. jtastk(1,ivstkoffset+i) .eq. mtfleftbra
     $         .or. jtastk(1,ivstkoffset+i) .eq. mtfleftparen
     $         .or. jtastk(1,ivstkoffset+i) .eq. mtfpart)then
            if(eol)then
              irtc=itfmessage(9999,'General::missop',
     $             '"'//string(ist1:min(istop-1,l))//'"')
              go to 8800
            endif
            jtastk(1,ivstkoffset+isp)=mtftimes
            call tfestk(isp0,iprior,lastfirst,irtc)
            if(irtc .ne. 0)then
              go to 8900
            endif
            isp=isp+1
            go to 900
          endif
        enddo
c     Sanity check for SAD stack
c     Stack MIGHT not have mtfright(brace|bra|paren),
c     because stack does not have mtfleft(brace|bra|paren)!!
        do i=isp-1,isp0,-1
          if(jtastk(1,ivstkoffset+i) .eq. mtfrightbrace
     $         .or. jtastk(1,ivstkoffset+i) .eq. mtfrightbra
     $         .or. jtastk(1,ivstkoffset+i) .eq. mtfrightparen)then
            mopc=1
            go to 8050
          endif
        enddo
        jtastk(1,ivstkoffset+isp)=mtfnull
        istop=istart
        irt=-2
        go to 2
      endif
 900  ktastk(isp)=kx
      jtastk(1,ivstkoffset+isp)=mtfnull
c
 1010 istart=istop
      eol=.false.
      go to 1
 7000 kt=iand(ktfmask,ktastk(isp))
      if(kt .ne. ktflist .and. kt .ne. ktfpat .and.
     $     kt .ne. ktfsymbol)then
        go to 8001
      endif
      ist10=max(ist1,istop)
      ist1=istop
      call tclrfpe
c      call tfdebugprint(ktastk(isp),'tfeval-8',3)
      call tfeevalref(ktastk(isp),kx,irtc)
c      call tfdebugprint(kx,'tfeval-9',3)
      istop=max(icsmrk(),ist1,ist10)
      if(irtc .eq. -1)then
        kx=ktfoper+mtfnull
        go to 9000
      elseif(irtc .eq. -5)then
        kx=ktfoper+mtfnull
        irtc=itfmessage(999,'General::abort',
     $     '"'//string(ist1:min(istop-1,l))//'"')
        go to 8900
      elseif(irtc .eq. -6)then
        kx=ktfoper+mtfnull
        irtc=itfmessage(9999,'General::abort',
     $     '"'//string(ist1:min(istop-1,l))//'"')
        go to 8901
      elseif(irtc .ne. 0)then
        ist1=ist10
        go to 8900
      endif
      if(itgetfpe() .ne. 0)then
        irtc=itfmessage(9,'General::fpe',
     $     '"'//string(ist1:min(istop-1,l))//'"')
        go to 8900
      endif
      go to 9000
 8001 kx=ktastk(isp)
      go to 9000
 8010 irtc=itfmessage(9999,'General::invop',
     $     '"'//string(ist1:min(istop-1,l))//'"')
      go to 8900
 8040 irtc=itfmessage(9999,'General::clquote','""')
      go to 8900
 8050 irtc=itfmessage(9999,'General::mismatch',
     $     '"'//opcode(mopc)//'"')
      go to 8900
 8110 irtc=itfmessage(9999,'General::stack','""')
      go to 8900
 8120 level1=levele
      levele=0
      irtc=itfmessage(9999,'General::deep','""')
      levele=level1
      go to 8900
 8600 irtc=itfmessage(999,'General::incomplete',
     $     '"'//string(ist1:min(istop-1,l))//'"')
      go to 8710
 8700 irtc=itfmessage(9999,'General::incomplete',
     $     '"'//string(ist1:min(istop-1,l))//'"')
 8710 kx=ktfoper+mtfnull
      istop=ist1
 8800 if(re .and. icslfni() .eq. 5)then
        if(irtc .gt. 0 .and. ierrorprint .ne. 0)then
          call tfreseterror
        endif
        go to 8910
      endif
 8900 if(irtc .lt. 0 .and. irtc .ne. -6)then
        modethrow=-1
        if(irtc .lt. -6)then
          call tflocal(kerror)
        endif
        kx=ktfoper+mtfnull
        irtc=itfmessage(999,'General::unexpbreak',
     $       '"'//string(ist1:min(l,ist1+20))//'"')
      endif
 8901 if(ierrorprint .ne. 0)then
        kx=ktfoper+mtfnull
        ierrorf=0
        if(kerror .ne. 0)then
          call tfaddmessage(kerror,string(1:l),min(istop,l+1),icslfno())
        endif
      endif
 8910 istop=ifchar(string(1:l),char(10),ist1)+1
      if(istop .eq. 1)then
        istop=l+1
      endif
 9000 if(levele .gt. 1)then
        call tfconnectk(kx,irtc)
      endif
      isp=isp0-1
      return
      include 'inc/TFSF.inc'
      end
