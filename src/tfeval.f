      function tfeval(string,ist1,istop,re,irtc) result(kx)
      use tfstk
      use ophash
      use opdata
      use tfcsi
      use eeval
      implicit none
      type (sad_descriptor) kx
      type (sad_dlist), pointer :: kla,klx
      logical*4 , intent(in) ::re
      integer*4 , intent(in) :: ist1
      integer*4 ,intent(out):: istop,irtc
      integer*4 istart,isp0,ist10,iop1,
     $     i,ishash,l,ifchar,mopc,itgetfpe,m1,iste,
     $     itfmessage,itfmessagestr,level1,ist2,irt,l0
      character*(*) , intent(in) :: string
      logical*4 tfreadevalbuf,eol
      type (csiparam) sav
c     begin initialize for preventing compiler warning
      mopc=0
c     end   initialize for preventing compiler warning
      l=len(string)
      l0=l
      isp0=isp+1
      isp=isp0
      ierrorf=0
      ierrorprint=0
      itastk2(1,isp)=mtfnull
      ktastk(isp)=ktfoper
      istart=ist1
      istop=istart
      iste=ist1
      ishash=-1
      if(.not. re)then
        itastk2(1,isp)=mtfleftparen
        isp=isp+1
        itastk2(1,isp)=mtfnull
        ktastk(isp)=ktfoper+mtfnull
      endif
      levele=levele+1
      if(levele >= maxlevele-10)then
        go to 8120
      endif
      irtc=0
      eol=.false.
 1    continue
c      write(*,*)'eval-1 ',string(istart:l)
      call tfetok(string(istart:l),istop,kx,itfcontext,irt)
      istop=min(l+1,istop+istart-1)
      if(irt >= 0)then
        go to 2400
      endif
 2    select case(irt)
c 2    go to (2100,2200,2300,8040),-irt
c            oper eol  cmnt clq
c
      case (-1)
        mopc=int(ktfaddrd(kx))
 2101   if(ktastk(isp) == ktfoper+mtfnull)then
          if(isp > isp0)then
            m1=itastk2(1,isp-1)
            select case(m1)
            case(mtfslot,mtfslotseq)
              go to 910
            case(mtfrightbra)
              if(mopc /= mtfrightbra)then
                do i=isp-2,isp0,-1
                  if(itastk2(1,i) == mtfpart)then
                    isp=isp-1
                    dtastk(i)=kxmakelist(i-1,klx)
                    klx%head%k=ktfoper+mtfpart
                    itastk2(1,i)=mtfnull
                    isp=i
                    go to 910
                  endif
                enddo
                go to 8050
              endif
            case(mtffun)
              dtastk(isp-1)=kxpfaloc(dtastk(isp-1))
              isp=isp-1
              go to 2101
            end select
            select case (mopc)
            case (mtfminus)
              mopc=merge(mtfneg,mtftimes,m1 == mtfpower .or. m1 == mtfrevpower)
              rtastk(isp)=-1.d0
            case (mtfplus)
              go to 1010
            case (mtftimes)
              if(m1 == mtfdot)then
                kx%k=ktfoper+mtfnull
                istop=ist1
                irtc=0
                go to 8910
              endif
              go to 910
            case (mtfpattest)
              mopc=mtfflag
            case default
              if(.not. nullfirst(mopc))then
                if(m1 == mtfpower .and. mopc == mtfpower)then
                  go to 8600
                endif
                go to 8700
              endif
            end select
          else
            if(mopc == mtfcomp .and. itastk2(1,isp) == mtfnull)then
              kx%k=ktfoper+mtfnull
              irtc=-1
              if(re)then
                istop=istop-1
              endif
              go to 9000
            elseif(mopc == mtfminus)then
              mopc=mtftimes
              rtastk(isp)=-1.d0
            elseif(mopc == mtfplus)then
              go to 1010
            elseif(mopc == mtfpattest)then
              mopc=mtfflag
            elseif(.not. nullfirst(mopc))then
              go to 8700
            endif
          endif
        else
          if(mopc == mtfleftparen .or. mopc == mtflist
     $         .or. mopc == mtfslot
     $         .or. mopc == mtfslotseq)then
            do i=isp0,isp-1
              if(itastk2(1,i) == mtflist
     $             .or. itastk2(1,i) == mtfleftbra
     $             .or. itastk2(1,i) == mtfleftparen
     $             .or. itastk2(1,i) == mtfpart)then
                itastk2(1,isp)=mtftimes
                call tfestk(isp0,iprior,lastfirst,irtc)
                if(irtc /= 0)then
                  go to 8900
                endif
                isp=isp+1
                go to 910
              endif
            enddo
            itastk2(1,isp)=mtfnull
            istop=istart
            irt=-2
            go to 2
          endif
        endif
 910    if(mopc > 0)then
          itastk2(1,isp)=mopc
          if(isp > isp0)then
            call tfestk(isp0,iprior,lastfirst,irtc)
            if(irtc /= 0)then
              go to 8900
            endif
          endif
          select case(mopc)
          case (mtfcomma)
            if(isp == isp0)then
              go to 7000
            endif
          case (mtfcomp)
            if(isp == isp0 .and. re)then
              istop=istop-1
              go to 7000
            endif
          case (mtfminus)
            itastk2(1,isp)=mtfplus
            isp=isp+1
            rtastk(isp)=-1.d0
            itastk2(1,isp)=mtftimes
          case (mtfdiv)
            itastk2(1,isp)=mtftimes
            isp=isp+1
            rtastk(isp)=-1.d0
            itastk2(1,isp)=mtfrevpower
          case (mtfunset,mtfrepeated,mtfrepeatednull,mtfincrement,mtfdecrement)
            if(ktastk(isp) /= ktfoper+mtfnull)then
              dtastk(isp)=kxmakelist(isp-1,kla)
              kla%head%k=ktfoper+mopc
              itastk2(1,isp)=mtfnull
            endif
          case default
            if(itastk2(1,isp) == mtfrightbrace .or. itastk2(1,isp) == mtfrightparen)then
              go to 8050
            endif
          end select
          if(itastk2(1,isp) /= mtfnull)then
            isp=isp+1
            if(isp > mstk)then
              go to 8110
            endif
            itastk2(1,isp)=mtfnull
            ktastk(isp)=ktfoper+mtfnull
          endif
          go to 1010
        else
          go to 8010
        endif
      case (-2)
        if(ktastk(isp) == ktfoper+mtfnull)then
          if(isp == isp0 .and.
     $         itastk2(1,isp) == mtfnull)then
            irtc=-1
          elseif(isp > isp0)then
            iop1=itastk2(1,isp-1)
            if(iop1 == mtffun)then
              isp=isp-1
              dtastk(isp)=kxpfaloc(dtastk(isp))
              itastk2(1,isp)=mtfnull
              go to 3
            endif
            if(re)then
                if(tfreadevalbuf(istart,istop,l,iop1))then
                  eol=.true.
                  go to 1
                endif
            elseif(iop1 == mtfcomp .or. iop1 == mtfleftparen)then
              go to 3
            endif
            go to 8700
          endif
          kx=dtastk(isp0)
          go to 9000
        endif
 3      if(isp > isp0)then
          if(.not. re)then
            itastk2(1,isp)=mtfrightparen
            call tfestk(isp0,iprior,lastfirst,irtc)
            if(irtc /= 0)then
              go to 8900
            endif
            if(isp .le. isp0)then
              go to 7000
            endif
          endif
          do i=isp0,isp
            if(itastk2(1,i) == mtflist
     $           .or. itastk2(1,i) == mtfleftbra
     $           .or. itastk2(1,i) == mtfleftparen
     $           .or. itastk2(1,i) == mtfpart)then
              if(re)then
c                if(string(istop-1:istop-1) == char(10))then
                  if(tfreadevalbuf(istart,istop,l,
     $                 int(itastk2(1,i))))then
                    itastk2(1,isp)=mtfnull
                    eol=.true.
                    go to 1
                  endif
c                endif
              endif
              go to 8700
            endif
          enddo
c     Sanity check for SAD stack
c     Stack MIGHT not have mtfright(brace|bra|paren),
c     because stack does not have mtfleft(brace|bra|paren)!!
          do i=isp,isp0,-1
            if(itastk2(1,i) == mtfrightbrace
     $           .or. itastk2(1,i) == mtfrightbra
     $           .or. itastk2(1,i) == mtfrightparen)then
              mopc=1
              go to 8050
            endif
          enddo
          itastk2(1,isp)=mtfend
          call tfestk(isp0,iprior,lastfirst,irtc)
          if(irtc /= 0)then
            go to 8900
          endif
        endif
        go to 7000
      case (-3)
        if(re)then
c          if(string(istop-1:istop-1) == char(10))then
 3001       if(tfreadevalbuf(istart,istop,l,
     $           mtfleftcomment))then
              ist2=index(string(istart:l),'*)')
              if(ist2 .le. 0)then
                go to 3001
              endif
              istart=ist2+istart+1
              eol=.false.
              go to 1
            endif
c          endif
        endif
        irtc=itfmessagestr(9999,'General::comment',
     $       string(ist1:min(istop-1,l)))
        istop=ist1
        go to 8800
      case (-4)
        go to 8040
      case default
        write(*,*)'tfeval undefined return code ',irt
        call abort
      end select

 2400 if(ktastk(isp) /= ktfoper+mtfnull)then
        do i=isp0,isp-1
          if(itastk2(1,i) == mtflist
     $         .or. itastk2(1,i) == mtfleftbra
     $         .or. itastk2(1,i) == mtfleftparen
     $         .or. itastk2(1,i) == mtfpart)then
            if(eol)then
              irtc=itfmessagestr(9999,'General::missop',
     $             string(ist1:min(ist1+80,istop-1,l)))
              go to 8800
            endif
            itastk2(1,isp)=mtftimes
            call tfestk(isp0,iprior,lastfirst,irtc)
            if(irtc /= 0)then
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
          if(itastk2(1,i) == mtfrightbrace
     $         .or. itastk2(1,i) == mtfrightbra
     $         .or. itastk2(1,i) == mtfrightparen)then
            mopc=1
            go to 8050
          endif
        enddo
        itastk2(1,isp)=mtfnull
        istop=istart
        irt=-2
        go to 2
      endif
 900  dtastk(isp)=kx
      itastk2(1,isp)=mtfnull
c
 1010 istart=istop
      eol=.false.
      go to 1
 7000 select case(ktftype(ktastk(isp)))
      case (ktflist,ktfpat,ktfsymbol)
        ist10=max(ist1,istop)
        ist2=istop
        call tclrfpe
        istop=max(ist2,ist10)
        if(re)then
          sav=savep
          ipoint=istop
        endif
        kx=tfeevalref(dtastk(isp),irtc)
        if(re)then
          savep=sav
        endif
        if(irtc == -1)then
          kx%k=ktfoper+mtfnull
          go to 9000
        elseif(irtc == irtcabort)then
          kx%k=ktfoper+mtfnull
          irtc=itfmessagestr(999,'General::abort',
     $         string(ist1:min(istop-1,l)))
          go to 8900
        elseif(irtc /= 0)then
          iste=ist10
          go to 8900
        endif
        if(itgetfpe() /= 0)then
          irtc=itfmessagestr(9,'General::fpe',
     $         string(ist1:min(istop-1,l)))
          go to 8900
        endif
        go to 9000
      case default
        kx=dtastk(isp)
      end select
      go to 9000
 8010 irtc=itfmessagestr(9999,'General::invop',
     $     string(ist1:min(istop-1,l)))
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
 8600 irtc=itfmessagestr(999,'General::incomplete',
     $     string(ist1:min(istop-1,l)))
      go to 8710
 8700 irtc=itfmessagestr(9999,'General::incomplete',
     $     string(ist1:min(istop-1,l)))
 8710 kx%k=ktfoper+mtfnull
      istop=ist1
 8800 if(re .and. icslfni() == 5)then
        go to 8910
      endif
 8900 if(irtc .lt. -1 .and. irtc > irtcabort)then
        modethrow=-1
        if(irtc .le. irtcret)then
          call tfreseterror
        endif
        kx%k=ktfoper+mtfnull
        irtc=itfmessagestr(999,'General::unexpbreak',
     $       string(max(ist1,iste-16):min(l,iste+20)))
      endif
      if(ierrorprint /= 0)then
        kx%k=ktfoper+mtfnull
        ierrorf=0
        if(kerror /= 0)then
          call tfaddmessage(string(1:l),min(istop,l+1),icslfnm())
        endif
      endif
 8910 istop=ifchar(string(1:l),char(10),iste)+1
      if(istop .le. 1)then
        istop=l+1
      endif
 9000 if(levele > 1)then
        call tfconnect(kx,irtc)
      endif
      isp=isp0-1
      return
      end
