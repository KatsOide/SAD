      module strbuf
      use tfstk

      type sad_strbuf
      sequence
      type (sad_string) string(1:0)
      integer*4 indw,llevel,remlines,maxllevel,column,lexp,nch,maxnch
      integer*1 istr(1:0)
      character*(mbody1) str
      end type

      contains
        subroutine strbuf_loc(locp,loc)
        use iso_c_binding
        implicit none
        type (sad_strbuf), pointer, intent(out) :: loc
        integer*8 locp
        call c_f_pointer(c_loc(klist(locp-3)),loc)
        return
        end subroutine

        subroutine descr_strbuf(k,strb)
        use iso_c_binding
        implicit none
        type (sad_descriptor) k
        type (sad_strbuf), pointer, intent(out) :: strb
        call c_f_pointer(c_loc(klist(ktfaddr(k)-3)),strb)
        return
        end subroutine

        subroutine getstringbuf(strb,n,stk)
        use tfmem
        implicit none
        type (sad_strbuf), pointer, intent(out) :: strb
        integer*8 kbuf
        integer*4 n,m,minsize,mexp,l,maxint
        logical*4 stk
        parameter (minsize=1024,maxint=1073741823)
        if(n .lt. 0)then
          m=min(minsize,-n)
          mexp=-n
        elseif(n .eq. 0)then
          m=minsize
          mexp=-1
        else
          m=n
          mexp=n
        endif
        l=m/8+2
        if(stk)then
          kbuf=isp+4+ispbase
          isp=isp+l+3
          if(isp .ge. mstk)then
            isp=isp-l-3
            kbuf=ktzaloc(ktfstring,l)
          endif
        else
          kbuf=ktzaloc(ktfstring,l)
        endif
        if(kbuf .le. 0)then
          write(6,*)'Memory allocation error (getstringbuf), size =',
     $         m
          call abort
        endif
        call strbuf_loc(kbuf,strb)
        strb%indw=maxint
        strb%llevel=0
        strb%remlines=maxint
        strb%maxllevel=0
        strb%column=0
        strb%lexp=mexp
        strb%nch=0
        strb%maxnch=m
c     ilist(1,kbuf)=0         ! nch
c     ilist(2,kbuf)=m         ! max nch
c     ilist(1,kbuf-1)=0       ! column
c     ilist(2,kbuf-1)=mexp    ! lexp; indent = lexp > 0
c     ilist(1,kbuf-2)=maxint  ! remaining lines
c     ilist(2,kbuf-2)=0       ! max level?
c     ilist(1,kbuf-3)=maxint  ! indent width
c     ilist(2,kbuf-3)=0       ! llevel
        return
        end subroutine

        subroutine tfquotestring(strb,string,l,lfno,irtc)
        use ISO_C_BINDING
        implicit none
        type (sad_strbuf), pointer :: strb
        integer*4 l,lfno,irtc,i,jp
        character*(l) string
        character ch
        character*(l*4+2) buff
        character*2 str
        character*4 buf
        jp=1
        buff(1:1)='"'
        LOOP_I: do i=1,l
          ch=string(i:i)
          if(ch .eq. C_NEW_LINE)then
            str='\\n'
            go to 20
          elseif(ch .eq. C_CARRIAGE_RETURN)then
            str='\\r'
            go to 20
          elseif(ch .eq. C_FORM_FEED)then
            str='\\f'
            go to 20
          elseif(ch .eq. C_HORIZONTAL_TAB)then
            str='\\t'
            go to 20
          elseif(ch .eq. '"')then
            str='\\"'
            go to 20
          elseif(ch .eq. '\')then
c'\
            str='\\\\'
            go to 20
          elseif(ichar(string(i:i)) .lt. 32 .or.
     $           ichar(string(i:i)) .gt. 126)then
            write(buf,'(1h\,o3.3)')string(i:i)
            buff(jp+1:jp+1)=buf(1:1)
            buff(jp+2:jp+2)=buf(2:2)
            buff(jp+3:jp+3)=buf(3:3)
            buff(jp+4:jp+4)=buf(4:4)
            jp=jp+4
            cycle LOOP_I
          else          
            jp=jp+1
            buff(jp:jp)=ch
            cycle LOOP_I
          endif
 20       jp=jp+2
          buff(jp-1:jp-1)=str(1:1)
          buff(jp:jp)    =str(2:2)
        enddo LOOP_I
        jp=jp+1
        buff(jp:jp)='"'
        call putstringbufpb(strb,buff,jp,.true.,lfno,irtc)
        return
        end subroutine

        subroutine tfconvbase(strb,v,base)
        implicit none
        type (sad_strbuf), pointer :: strb
        integer*4 lbuf,i,ifrac,lmax
        parameter (lbuf=1040)
        real*8 v,av,av1,af,af1,v1,base
        character*(lbuf) buf
        character ch
        logical*4 fr,full
        if(v .lt. 0.d0)then
          call putstringbufb1(strb,'-')
          v1=-v
        else
          v1=v
        endif
        av=aint(v1)
        fr=av .eq. 0.d0
        af=v1-av
        i=lbuf
        if(av .eq. 0.d0)then
          buf(i:i)='0'
          i=i-1
        else
          do while(av .ne. 0.d0)
            av1=aint(av/base)
            ifrac=int(av-av1*base)
            if(ifrac .lt. 10)then
              ch=char(ichar('0')+ifrac)
            elseif(ifrac .lt. 36)then
              ch=char(ichar('a')+ifrac-10)
            else
              ch=char(ichar('A')+ifrac-36)
            endif
            buf(i:i)=ch
            i=i-1
            av=av1
          enddo
        endif
        call putstringbufb(strb,buf(i+1:lbuf),lbuf-i,full)
        if(af .eq. 0.d0)then
          return
        endif
        lmax=int(56.d0*log(2.d0)/log(base)-(lbuf-i))
        if(lmax .le. 0)then
          return
        endif
        i=1
        do while(af .ne. 0.d0 .and. i .le. lmax)
          af1=aint(af*base)
          ifrac=int(af1)
          if(ifrac .eq. 0)then
            if(fr)then
              lmax=min(lbuf,lmax+1)
            endif
            ch='0'
          else
            fr=.false.
            if(ifrac .lt. 10)then
              ch=char(ichar('0')+ifrac)
            elseif(ifrac .lt. 36)then
              ch=char(ichar('a')+ifrac-10)
            else
              ch=char(ichar('A')+ifrac-36)
            endif
          endif
          buf(i:i)=ch
          i=i+1
          af=af*base-af1
        enddo
        if(i .gt. 1)then
          call putstringbufb(strb,'.'//buf(:i-1),i,full)
        endif
        return
        end subroutine
 
        recursive subroutine tfconvstrb(strb,
     $     k,nc,str,gens,lfno,form,irtc)
        use tfcode
        use iso_c_binding
        implicit none
        type (sad_descriptor) k
        type (sad_pat), pointer :: pat
        type (sad_namtbl), pointer :: loc
        type (sad_strbuf), pointer :: strb
        type (sad_string), pointer :: strx
        type (sad_symbol), pointer :: sym
        integer*8 kv,ktv,kav,ic,i,ip
        integer*4 nc,lenw,lfno,irtc,nc0,nc1,kpat
        real*8 v
        character*27 buff
        character*26 form1,tfgetform,autos1,autofg
        character*(*) form
        character*3 patstr
        data patstr/'___'/
        logical*4 str,gens
        if(ktfstringq(k,strx))then
          nc=strx%nch
          if(str)then
            nc0=strb%nch
            call tfquotestring(strb,strx%str,nc,lfno,irtc)
            if(irtc .ne. 0)then
              return
            endif
            nc=max(strb%nch-nc0,0)
          else
            call putstringbufpb(strb,strx%str,nc,.false.,lfno,irtc)
          endif
        elseif(ktfsymbolq(k,sym))then
          call loc_namtbl(sym%loc,loc)
          ic=loc%cont
          do i=itfcontextpath,
     $         itfcontextpath+ilist(2,itfcontextpath-1)-1
            if(klist(i) .eq. ic)then
              nc=0
              go to 1
            endif
          enddo
          if(ic .eq. itfcontroot)then
            call putstringbufp(strb,'`',lfno,irtc)
            nc=1
          elseif(klist(ic) .ne. 0)then
            call tfconvstrb(strb,dlist(ic),
     $           nc,.false.,gens,lfno,form,irtc)
          endif
          if(irtc .ne. 0)then
            return
          endif
 1        nc1=loc%str%nch
          call putstringbufpb(strb,loc%str%str,nc1,.true.,lfno,irtc)
          if(irtc .ne. 0)then
            return
          endif
          nc=nc+nc1
          if(.not. gens)then
            if(sym%gen .gt. 0 .and.
     $           sym%gen .ne. maxgeneration)then
              buff='$'//autos1(dble(sym%gen))
              nc1=lenw(buff)
              call putstringbufpb(strb,buff,nc1,.true.,lfno,irtc)
              nc=nc+nc1
            endif
          endif
        elseif(ktfoperq(k))then
          call loc_namtbl(klist(klist(ifunbase+ktfaddr(k))),loc)
          nc=loc%str%nch
          call putstringbufpb(strb,loc%str%str,nc,.true.,lfno,irtc)
        elseif(ktfpatq(k,pat))then
          kv=pat%expr%k
          kav=ktfaddr(kv)
          ktv=kv-kav
          kpat=merge(merge(int(kav),0,kav .le. 3),-1,
     $         ktv .eq. ktfref)
          if(pat%default%k .ne. ktfref)then
            if(kpat .lt. 0)then
              call putstringbufp(strb,'((',lfno,irtc)
              nc=2
            else
              call putstringbufp(strb,'(',lfno,irtc)
              nc=1
            endif
          else
            nc=0
          endif
          if(pat%sym%loc .ne. 0)then
            ip=ktfaddr(pat%sym%alloc%k)
            if(ip .gt. 0)then
              call tfconvstrb(strb,transfer(ktfsymbol+ip,k),
     $             nc1,.false.,gens,lfno,form,irtc)
              if(irtc .ne. 0)then
                return
              endif
              nc=nc+nc1
            endif
          endif
          if(kpat .gt. 0)then
            call putstringbufp(strb,patstr(1:kpat),lfno,irtc)
            if(irtc .ne. 0)then
              return
            endif
            nc=nc+kpat
            if(pat%head%k .ne. 0)then
              call tfconvstrb(strb,
     $             pat%head,nc1,.false.,gens,lfno,form,irtc)
              if(irtc .ne. 0)then
                return
              endif
              nc=nc+nc1
            endif
          else
            call putstringbufp(strb,':',lfno,irtc)
            if(irtc .ne. 0)then
              return
            endif
            call tfconvstrb(strb,pat%expr,nc1,.true.,
     $           gens,lfno,form,irtc)
            nc=nc+nc1+1
          endif
          if(pat%default%k .ne. ktfref)then
            if(kpat .lt. 0)then
              call putstringbufp(strb,'):',lfno,irtc)
              nc=nc+2
            else
              call putstringbufp(strb,':',lfno,irtc)
              nc=nc+1
            endif
            call tfconvstrb(strb,pat%default,
     $           nc1,.true.,gens,lfno,form,irtc)
            call putstringbufp(strb,')',lfno,irtc)
            nc=nc+nc1+1
          endif
        elseif(ktflistq(k))then
          nc0=strb%nch
          call tfconvstrl(strb,ktfaddr(k),lfno,form,gens,irtc)
          nc=strb%nch-nc0
        elseif(ktfrefq(k))then
          nc=0
        elseif(ktfenanq(k%x(1)))then
          nc=3
          call putstringbufpb(strb,'NaN',nc,.true.,lfno,irtc)        
        elseif(ktfrealq(k,v))then
          if(form .eq. '*')then
            form1=tfgetform()
          else
            form1=form
          endif
          buff=autofg(v,form1)
          nc=len_trim(buff)
          call putstringbufpb(strb,buff,nc,.true.,lfno,irtc)
        else
          write(*,*)'tfconvstrb ',ktftype(k%k),ktfaddr(k)
          nc=5
          call putstringbufpb(strb,' ??? ',nc,.true.,lfno,irtc)
        endif
        return
        end subroutine

        recursive subroutine tfconvstrl(strb,ka,lfno,form,gens,irtc)
        use ophash
        use opdata
        implicit none
        type (sad_descriptor) k1,ki
        type (sad_dlist), pointer ::list,listi
        type (sad_strbuf), pointer :: strb
        integer*8 ka,kt,kai
        integer*4 lfno,nd,iaaf,ncx,nc,i,irtc,i1,llevel,lenw,le,istep,
     $       iaaf1
        real*8 v1
        character*(*) form
        character*4 opcx
        character*2 opce
        logical*4 gens
        data llevel/0/
c     Initialize to avoid compiler warning
        ncx=-1
c     
        istep=1
        nd=ilist(2,ka-1)
        i1=1
        if(tfexprq(ktflist+ka))then
          iaaf=int(ktfaddr(klist(ka)))
          kt=klist(ka)-iaaf
          if(kt .eq. ktfoper .and. iaaf .le. mtfend
     $         .and. (nd .ge. 2
     $         .or. iaaf .eq. mtfslot
     $         .or. iaaf .eq. mtfslotseq
     $         .or. iaaf .eq. mtfincrement
     $         .or. iaaf .eq. mtfdecrement
     $         .or. iaaf .eq. mtffun
     $         .or. iaaf .eq. mtfrepeated
     $         .or. iaaf .eq. mtfunset
     $         .or. iaaf .eq. mtfrepeatednull
     $         .or. iaaf .eq. mtfflag))then
            select case (iaaf)
            case (mtfcomplex)
              if(ktfrealq(klist(ka+1)))then
                if(rlist(ka+1) .eq. 0.d0)then
                  if(iand(ktrmask,klist(ka+2)) .ne. ktfnr)then
                    if(rlist(ka+2) .eq. 1.d0)then
                      call putstringbufp(strb,'I',lfno,irtc)
                      return
                    elseif(rlist(ka+2) .eq. -1.d0)then
                      call putstringbufp(strb,'(-I)',lfno,irtc)
                      return
                    endif
                  endif
                  call putstringbufp(strb,'(',lfno,irtc)
                  if(irtc .ne. 0)then
                    return
                  endif
                  call tfconvstrb(strb,dlist(ka+2),
     $                 nc,.true.,gens,lfno,form,irtc)
                  if(irtc .ne. 0)then
                    return
                  endif
                  call putstringbufp(strb,' I)',lfno,irtc)
                  return
                endif
              endif
              call putstringbufp(strb,'(',lfno,irtc)
              if(irtc .ne. 0)then
                return
              endif
              call tfconvstrb(strb,dlist(ka+1),
     $             nc,.true.,gens,lfno,form,irtc)
              if(irtc .ne. 0)then
                return
              endif
              if(iand(ktrmask,klist(ka+2)) .ne. ktfnr)then
                if(rlist(ka+2) .ge. 0.d0)then
                  call putstringbufp(strb,'+',lfno,irtc)
                endif
                if(rlist(ka+2) .eq. 1.d0)then
                  call putstringbufp(strb,'I)',lfno,irtc)
                  return
                elseif(rlist(ka+2) .eq. -1.d0)then
                  call putstringbufp(strb,'-I)',lfno,irtc)
                  return
                endif
              else
                call putstringbufp(strb,'+',lfno,irtc)
              endif
              if(irtc .ne. 0)then
                return
              endif
              call tfconvstrb(strb,dlist(ka+2),
     $             nc,.true.,gens,lfno,form,irtc)
              if(irtc .ne. 0)then
                return
              endif
              call putstringbufp(strb,' I)',lfno,irtc)
              return
            case (mtfmult)
              k1=dlist(ka+1)
              if(ktfrealq(k1,v1))then
                if(v1 .eq. -1.d0)then
                  call putstringbufp(strb,'(-',lfno,irtc)
                  if(irtc .ne. 0)then
                    return
                  endif
                  opce=')'
                  opcx=' '
                  ncx=1
                  i1=2
                  go to 101
                endif
              endif
              call putstringbufp(strb,'(',lfno,irtc)
              ncx=1
              opcx=' '
              opce=')'
            case (mtfnull)
              call putstringbufp(strb,'[',lfno,irtc)
              opcx=','
              ncx=1
              opce=']'
            case (mtfpart)
              k1=dlist(ka+1)
              call tfconvstrb(strb,k1,nc,
     $             .true.,gens,lfno,form,irtc)
              if(irtc .ne. 0)then
                return
              endif
              call putstringbufp(strb,'[[',lfno,irtc)
              if(irtc .ne. 0)then
                return
              endif
              opcx=','
              ncx=1
              opce=']]'
              i1=2
              go to 101
            case (mtfslot,mtfslotseq,mtfflag)
              opcx=opcode(iaaf)
              call putstringbufp(strb,opcx(1:lenw(opcx)),lfno,irtc)
              opce=' '
            case (mtfhold)
              call putstringbufp(strb,'Hold[',lfno,irtc)
              opcx=','
              ncx=1
              opce=']'
            case default
              opcx=opcode(iaaf)
              ncx=lenw(opcx)
              if(iaaf .eq. mtfmult)then
                opcx=' '
              endif
              opce=' '
            end select
            if(iaaf .eq. mtfincrement .or. iaaf .eq. mtfdecrement)then
              i1=nd
            elseif(iaaf .eq. mtfinequality)then
              istep=2
            endif
          else
            call tfconvstrb(strb,dlist(ka),nc,
     $           .true.,gens,lfno,form,irtc)
            if(irtc .ne. 0)then
              return
            endif
            iaaf=-1
            call putstringbufp(strb,'[',lfno,irtc)
            if(irtc .ne. 0)then
              return
            endif
            opcx=','
            ncx=1
            opce=']'
          endif
        else
          iaaf=mtflist
          call putstringbufp(strb,'{',lfno,irtc)
          if(irtc .ne. 0)then
            return
          endif
          opcx=','
          ncx=1
          opce='}'
        endif
        if(i1 .gt. 1)then
          call putstringbufp(strb,opcx(1:ncx),lfno,irtc)
          if(irtc .ne. 0)then
            return
          endif
        endif
 101    llevel=llevel+1
        call loc_sad(ka,list)
        do i=i1,nd,istep
          strb%llevel=llevel
          ki=list%dbody(i)
          if(ktflistq(ki,listi))then
            kai=ktfaddrd(ki)
            if(opcx .ne. ',' .and. ktfoperq(listi%head%k))then
              iaaf1=int(ktfaddr(listi%head))
              if(iaaf1 .le. mtfend .and. iaaf1 .ne. mtflist)then
                if(iprior(iaaf1) .ge. iprior(iaaf))then
                  call putstringbufp(strb,'(',lfno,irtc)
                  if(irtc .ne. 0)then
                    go to 8000
                  endif
                  call tfconvstrl(strb,kai,lfno,form,gens,irtc)
                  if(irtc .ne. 0)then
                    go to 8000
                  endif
                  call putstringbufp(strb,')',lfno,irtc)
                  go to 102
                endif
              endif
            endif
            call tfconvstrl(strb,kai,lfno,form,gens,irtc)
          elseif(ki%k .eq. ktfoper+mtfnull)then
          else
            call tfconvstrb(strb,ki,nc,.true.,gens,lfno,form,irtc)
          endif
 102      if(irtc .ne. 0)then
            go to 8000
          endif
          if(iaaf .eq. mtfinequality .and. i .lt. nd)then
            ki=list%dbody(i+1)
            if(ki%k .le. ktfoper+mtfend)then
              opcx=opcode(iand(ktamask,ki%k))
              ncx=lenw(opcx)
            else
              opcx=','
              ncx=1
              call putstringbufp(strb,opcx(1:ncx),lfno,irtc)
              if(irtc .ne. 0)then
                return
              endif
              call tfconvstrb(strb,ki,nc,.true.,gens,lfno,form,irtc)
              if(irtc .ne. 0)then
                go to 8000
              endif
            endif
          elseif(nd .gt. 2 .and.
     $           (iaaf .eq. mtfmap .or. iaaf .eq. mtfapply))then
            opce='])'
            if(i .eq. 1)then
              opcx=opcx(1:2)//'['
              ncx=3
            elseif(i .gt. 1)then
              opcx=','
              ncx=1
            endif
          endif
          if(i .lt. nd .or. iaaf .eq. mtfunset .or.
     $         iaaf .eq. mtfrepeated .or. iaaf .eq. mtfrepeatednull
     $         .or. ((iaaf .eq. mtffun .or. iaaf .eq. mtfincrement
     $         .or. iaaf .eq. mtfdecrement) .and. nd .eq. 1))then
            call putstringbufp(strb,opcx(1:ncx),lfno,irtc)
            if(irtc .ne. 0)then
              go to 8000
            endif
          endif
          if(lfno .ge. 0 .and. strb%lexp .gt. 0)then
            if(strb%indw .ge. llevel)then
              strb%indw=llevel
              strb%column=strb%nch
            endif
          endif
        enddo
        le=lenw(opce)
        if(le .ne. 0)then
          call putstringbufp(strb,opce(1:le),lfno,irtc)
        endif
 8000   llevel=llevel-1
        return
        end subroutine

        subroutine tfconvreal(strb,x)
        implicit none
        type (sad_strbuf), pointer :: strb
        integer*4 l,lenw,ich
        real*8 x
        character*22 str1,autos1
        logical*4 full
        if(strb%nch .ne. 0)then
          ich=strb%istr(strb%nch)
          if(ich .ge. ichar('0') .and. ich .le. ichar('9'))then
            call putstringbufb1(strb,' ')
          endif
        endif
        str1=autos1(x)
        l=lenw(str1)
        call putstringbufb(strb,str1,l,full)
        return
        end subroutine

        subroutine extendstringbuf(strb,lnew)
        use tfmem
        implicit none
        type (sad_strbuf), pointer :: strb
        integer*8 i
        integer*4 lnew,l
        if(strb%lexp .ge. lnew .or. strb%lexp .eq. -1)then
          l=lnew/8+2
          i=ktzaloc(ktfstring,l)
          if(i .le. 0)then
            write(6,*)
     $           'Memory allocation error (extendstringbuf), size =',
     $           lnew
            call abort
          endif
          call tmov(strb%indw,ilist(1,i-3),strb%maxnch/8+5)
          ilist(2,i)=lnew
          call tfreestringbuf(strb)
        else
          i=0
        endif
        call strbuf_loc(i,strb)
        return
        end subroutine

        subroutine putstringbufp(strb,string,lfno,irtc)
        implicit none
        type (sad_strbuf), pointer :: strb
        integer*4 lfno,irtc
        character*(*) string
        call putstringbufpb(strb,string,len(string),.false.,lfno,irtc)
        return
        end subroutine

        subroutine putstringbufpb(strb,string,l,quote,lfno,irtc)
        implicit none
        type (sad_strbuf), pointer :: strb
        integer*4 lfno,irtc,l,lw,i,i1,lexp,lv,iext,itfmessage
        character*(l) string
        logical*4 full,indent,quote
        irtc=0
        if(l .le. 0)then
          return
        endif
        if(lfno .lt. 0)then
          if(l .eq. 1)then
            call putstringbufb1(strb,string)
          else
            call putstringbufb(strb,string,l,full)
            if(full)then
              irtc=itfmessage(9,'General::longstr',' ')
            endif
          endif
          return
        else
          lv=strb%llevel
          iext=merge(3,1,
     $         ichar(string(1:1)) .ge. ichar('0') .and.
     $         ichar(string(1:1)) .le. ichar('9')
     $         .or.
     $         ichar(string(1:1)) .ge. ichar('A') .and.
     $         ichar(string(1:1)) .le. ichar('Z')
     $         .or.
     $         ichar(string(1:1)) .ge. ichar('a') .and.
     $         ichar(string(1:1)) .le. ichar('z')
     $         )
          indent=strb%lexp .ge. 0
          i1=1
          do10: do while(i1 .le. l)
            lexp=strb%lexp
            if(lexp .gt. 0)then
              if(indent .and. strb%nch .gt. 0 .and.
     $             strb%llevel .lt. strb%maxllevel)then
                strb%column=0
                call flushstringbuf(strb,.true.,.true.,lfno,irtc)
                if(irtc .ne. 0)then
                  return
                endif
              endif
              lw=lexp-strb%nch
              if((.not. indent .or.
     $             strb%nch .gt. min(lv,32)) .and.
     $             (lw .le. 0 .or. lw .lt. l-i1+iext .and.
     $             lexp .ge. l-i1+1))then
                call flushstringbuf(strb,indent,.true.,lfno,irtc)
                if(irtc .ne. 0)then
                  return
                endif
                lw=lexp-strb%nch
              endif
            else
              lw=l-i1+1
            endif
            indent=.false.
            do i=i1,min(l,lw+i1-1)
              if(string(i:i) .eq. char(13))then
                call putstringbufb(strb,string(i1:i1),i-i1+1,full)
                call flushstringbuf(strb,indent,.false.,lfno,irtc)
                if(irtc .ne. 0)then
                  return
                endif
                i1=i+1
                cycle do10
              elseif(string(i:i) .eq. char(10))then
                call putstringbufb(strb,string(i1:i1),i-i1,full)
                call flushstringbuf(strb,indent,.true.,lfno,irtc)
                if(irtc .ne. 0)then
                  return
                endif
                i1=i+1
                cycle do10
              endif
            enddo
 20         if(lw .lt. l-i1+1)then
              if(quote)then
                call putstringbufb(strb,string(i1:i1),lw-1,full)
                call putstringbufb(strb,'\',1,full)
                call flushstringbuf(strb,indent,.true.,lfno,irtc)
                if(irtc .ne. 0)then
                  return
                endif
                i1=i1+lw-1
                lw=lexp-strb%nch
                go to 20
              else
                call putstringbufb(strb,string(i1:i1),lw,full)
                i1=i1+lw
              endif
            else
              call putstringbufb(strb,string(i1:i1),l-i1+1,full)
              i1=l+1
            endif
          enddo do10
        endif
        return
        end subroutine

        subroutine putstringbufb(strb,string,l,full)
        implicit none
        type (sad_strbuf), pointer :: strb
        integer*4 l,ip,lexp,ip1,l1
        character*(l) string
        logical*4 full
        ip=strb%nch
        ip1=ip+l
        lexp=strb%lexp
        full=lexp .gt. 0 .and. ip1 .gt. lexp
        if(full)then
          l1=lexp-ip
          if(l1 .le. 0)then
            return
          endif
          ip1=lexp
        else
          l1=l
        endif
        if(ip1 .gt. strb%maxnch)then
          if(lexp .gt. 0)then
            call extendstringbuf(strb,min(lexp,ip1+strb%maxnch))
          else
            call extendstringbuf(strb,ip1+strb%maxnch)
          endif
        endif
        strb%str(ip+1:ip+l1)=string(1:l1)
        strb%nch=ip1
        return
        end subroutine

        subroutine putstringbufb1(strb,string)
        implicit none
        type (sad_strbuf), pointer :: strb
        integer*4 ip,lexp,ip1
        character string
        ip=strb%nch
        ip1=ip+1
        lexp=strb%lexp
        if(lexp .gt. 0 .and. ip1 .gt. lexp)then
          if(lexp .le. ip)then
            return
          endif
          ip1=lexp
        endif
        if(ip1 .gt. strb%maxnch)then
          if(lexp .gt. 0)then
            call extendstringbuf(strb,min(lexp,ip1+strb%maxnch))
          else
            call extendstringbuf(strb,ip1+strb%maxnch)
          endif
        endif
        strb%str(ip1:ip1)=string
        strb%nch=ip1
        return
        end subroutine

        subroutine writestringbufn(strb,nl,lfno)
        implicit none
        type (sad_strbuf) strb
        integer*4 lfno
        logical*4 nl
        if(strb%nch .gt. 0)then
          call writeb(strb%str,strb%nch,nl,lfno)
          strb%nch=0
        endif
        return
        end subroutine

        subroutine writestringbuf(strb,nl,lfno)
        implicit none
        type (sad_strbuf) strb
        integer*4 lfno
        logical*4 nl
        call writeb(strb%str,strb%nch,nl,lfno)
        strb%nch=0
        strb%maxllevel=0
        return
        end subroutine

        type (sad_descriptor) function kxstringbuftostring(strb)
        implicit none
        type (sad_strbuf) strb
        type (sad_string), pointer :: str
        integer*8 kbuf,ktfaloc,ka,kp
        integer*4 n,l,m,i,n1,k
        kbuf=sad_loc(strb%nch)
        n=strb%nch
        if(n .eq. 0)then
          kxstringbuftostring=dxnulls
          if(.not. tfonstackq(kbuf))then
            strb%indw=strb%maxnch/8+5
            call tfree(kbuf-2)
          endif
        elseif(n .eq. 1)then
          k=strb%istr(1)
          if(k .lt. 0)then
            k=k+256
          endif
          kxstringbuftostring%k=ktfstring+iaxschar+k*5+3
          if(.not. tfonstackq(kbuf))then
            strb%indw=strb%maxnch/8+5
            call tfree(kbuf-2)
          endif
        else
          if(tfonstackq(kbuf))then
            m=n/8
c            write(*,*)'strbuftostr ',m,n,kbuf,ispbase,mstk+ispbase
            ka=ktfaloc(-1,ktfstring,m+2)
c            write(*,*)'strbuftostr-ka ',ka
            ilist(2,ka-3)=-1
            ilist(1,ka)=n
c            do i=1,m
              klist(ka+1:ka+m)=klist(kbuf+1:kbuf+m)
c            enddo
            klist(ka+m+1)=0
            do i=m*8+1,n
              jlist(i,ka+1)=jlist(i,kbuf+1)
            enddo
            kxstringbuftostring%k=ktfstring+ka
          else
            call loc_string(kbuf,str)
            str%len=strb%maxnch/8+5
            str%override=-1
            str%ref=1
            str%alloc%k=ktfstring
            str%gen=0
            n1=n/8
            l=n1+5
            n1=(n1+1)*8
c            do i=n+1,n1
            jlist(n+1:n1,kbuf+1)=int(0,1)
c            enddo
            if(ilist(1,kbuf-3)-l .ge. 4)then
              kp=kbuf+l-3
              ilist(1,kp)=ilist(1,kbuf-3)-l
              call tfree(kp+1)
              ilist(1,kbuf-3)=l
            endif
            ilist(2,kbuf)=0
            kxstringbuftostring%k=ktfstring+kbuf
            call tflocal1(kbuf)
          endif
        endif
        return
        end function

      end module
