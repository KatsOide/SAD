      module funs

      contains
      function tfsequence(isp1,isp2) result(kx)
      use tfstk
      implicit none
      type (sad_descriptor) kx
      integer*4 ,intent(in):: isp1,isp2
      kx=merge(dxnull,merge(dtastk(isp2),
     $     kxcrelistm(isp2-isp1,ktastk(isp1+1:isp2),
     $     k_descr(ktfoper+mtfnull)),isp1+1 .eq. isp2),
     $     isp1 .ge. isp2)
      return
      end

      function tfpuref(isp1,kf,irtc) result(kx)
      use tfstk
      use eeval
      implicit none
      type (sad_descriptor) kx,ki,ka
      type (sad_dlist) ,intent(in):: kf
      type (sad_dlist), pointer :: kla
      integer*4 ,intent(in):: isp1
      integer*4 ,intent(out):: irtc
      integer*4 itfmessage,narg,m,j,i,ipf0,nap0,isp0
      logical*4 rep
      if(kf%nl .eq. 1)then
        isp0=isp
        do i=isp1+1,isp
          dtastk(i)=dtfcopy(dtastk(i))
        enddo
        ipf0=ipurefp
        nap0=napuref
        ipurefp=isp1
        napuref=isp-isp1
        isp=isp+1
        itastk(1,isp)=ipf0
        itastk(2,isp)=nap0
        kx=tfeevalref(kf%dbody(1),irtc)
        ipurefp=ipf0
        napuref=nap0
        do i=isp1+1,isp0
          call tflocald(dtastk(i))
        enddo
        isp=isp0
      elseif(kf%nl .eq. 2)then
        narg=isp-isp1
        ka=kf%dbody(1)
        if(ktfsymbolq(ka))then
          if(narg .ne. 1)then
            irtc=itfmessage(9,'General::narg',
     $           '"equal to actual number of args"')
            return
          endif
          dtastk(isp+1)=ka
          isp=isp+2
          dtastk(isp)=dtastk(isp1+1)
        elseif(tflistq(ka,kla))then
          m=kla%nl
          if(m .ne. narg)then
            irtc=itfmessage(9,'General::narg',
     $           '"equal to actual number of args"')
            return
          endif
          if(m .ne. 0)then
            if(ktfreallistq(kla))then
              irtc=itfmessage(9,'General::wrongtype',
     $             '"List of symbols"')
              return
            endif
            do i=1,m
              ki=kla%dbody(i)
              if(.not. ktfsymbolq(ki))then
                irtc=itfmessage(9,'General::wrongtype',
     $               '"List of symbols"')
                return
              endif
              j=isp+i*2
              dtastk(j-1)=ki
              ktastk(j)=ktastk(isp1+i)
            enddo
            isp=isp+2*m
          endif
        else
          irtc=itfmessage(9,'General::wrongtype','"List of symbols"')
          return
        endif
        kx=kf%dbody(2)
        if(narg .ne. 0)then
          call tfreplacesymbolstk(kx,isp1+narg,narg,kx,.true.,rep,irtc)
c          call tfdebugprint(kx,'puref-2',3)
c          write(*,*)irtc
          if(irtc .ne. 0)then
            isp=isp1+narg
            return
          endif
        endif
        kx=tfeevalref(kx,irtc)
        isp=isp1+narg
      else
        irtc=itfmessage(9,'General::narg','"1 or 2"')
      endif
      return
      end

      recursive function tfset(isp1,upvalue,irtc) result(kx)
      use tfstk
      use eeval
      implicit none
      type (sad_descriptor) kx,k1,k2,k110,k11,tfset1
      type (sad_symbol), pointer :: sym
      type (sad_symdef), pointer :: symd
      type (sad_dlist), pointer :: kl11
      integer*4 ,intent(in):: isp1
      integer*4 ,intent(out):: irtc
      integer*4 i,itfmessage,isp11,isp0
      logical*4 ,intent(in):: upvalue
      logical*4 euv
      if(isp1+1 .ge. isp)then
        kx=dxnull
        irtc=itfmessage(9,'General::narg','"2 or more"')
        return
      endif
      k2=dtastk(isp)
      if(isp1+2 .eq. isp)then
        k1=dtastk(isp1+1)
        if(ktfsymbolq(k1,sym))then
          if(sym%override .eq. 0)then
            sym=>tfsydef(sym)
            k1=sad_descr(sym)
          endif
        elseif(ktflistq(k1))then
          call tfeevaldef(k1,k1,irtc)
          if(irtc .ne. 0)then
            return
          endif
        endif
        if(k1%k .ne. ktastk(isp1+1) .and. upvalue)then
          k11=k1
          k110=k1
          if(ktflistq(k11,kl11))then
            k11=kl11%head
            if(ktflistq(k11))then
              k110=k11
              do while(ktflistq(k11,kl11))
                k11=kl11%head
              enddo
            endif
          endif
          if(ktfsymbolq(k11,sym))then
            if(sym%override .eq. 0)then
              sym=>tfsydef(sym)
            endif
            call sym_symdef(sym,symd)
            if(symd%upval .ne. 0)then
              isp=isp+1
              isp11=isp
              dtastk(isp11)=dtastk(isp1)
              isp=isp+1
              dtastk(isp)=k1
              isp=isp+1
              dtastk(isp)=k2
              call tfdeval(isp11,ksad_loc(sym%loc),kx,0,
     $             .false.,euv,irtc)
              isp=isp11-1
              if(euv)then
                return
              endif
            endif
          endif
        endif
        kx=tfset1(k1,k2,ktfaddr(ktastk(isp1)),irtc)
      else
        isp0=isp
        kx=dtastk(isp)
        do i=isp-1,isp1+1,-1
          isp=isp0+3
          dtastk(isp-2)=dtastk(isp1)
          dtastk(isp-1)=dtastk(i)
          dtastk(isp  )=kx
          kx=tfset(isp0+1,upvalue,irtc)
          if(irtc .ne. 0)then
            return
          endif
        enddo
        isp=isp0
      endif
      return
      end

      function tfreplace1(isp1,irtc) result(kx)
      use tfstk
      implicit none
      type (sad_descriptor) kx,k1
      integer*4 ,intent(in):: isp1
      integer*4 ,intent(out):: irtc
      integer*4 i,itfmessage
      kx=dxnullo
      if(isp .le. isp1+1)then
        irtc=itfmessage(9,'General::narg','"2 or more"')
        return
      endif
      k1=dtastk(isp1+1)
      do i=isp1+2,isp
        call tfreplace(k1,dtastk(i),kx,.true.,.true.,.false.,irtc)
        if(irtc .ne. 0)then
          return
        endif
        k1=kx
      enddo
      return
      end

      function tfreplacerepeated1(isp1,irtc) result(kx)
      use tfstk
      implicit none
      type (sad_descriptor) kx
      integer*4 ,intent(in):: isp1
      integer*4 ,intent(out):: irtc
      integer*4 i,itfmessage
      kx=dxnullo
      if(isp .le. isp1+1)then
        irtc=itfmessage(9,'General::narg','"2 or more"')
        return
      endif
      kx%k=ktastk(isp1+1)
      irtc=0
      do i=isp1+2,isp
        call tfreplacerepeated(kx,ktastk(i),kx,.true.,.true.,irtc)
        if(irtc .ne. 0)then
          return
        endif
      enddo
      return
      end

      function tfsameq1(isp1,iopc,irtc) result(kx)
      use tfstk
      implicit none
      type (sad_descriptor) kx,k,k1
      integer*4 ,intent(in):: isp1,iopc
      integer*4 ,intent(out):: irtc
      integer*4 itfmessage
      kx=dxnullo
      if(isp .ne. isp1+2)then
        irtc=itfmessage(9,'General::narg','"2"')
        return
      endif
      k=dtastk(isp)
      k1=dtastk(isp1+1)
      kx%k=0
      if(k%k .eq. k1%k)then
        kx%k=ktftrue
      elseif(tfsameq(k,k1))then
        kx%k=ktftrue
      endif
      if(iopc .eq. mtfunsame)then
        kx%k=ktftrue-kx%k
      endif
      irtc=0
      return
      end

      function tfequal(isp1,iopc,irtc) result(kx)
      use tfstk
      implicit none
      type (sad_descriptor) kx,tfrelation
      integer*4 ,intent(in):: isp1,iopc
      integer*4 ,intent(out):: irtc
      integer*4 itfmessage
      if(isp .lt. isp1+2)then
        irtc=itfmessage(9,'General::narg','"2 or more"')
        return
      endif
      if(isp .eq. isp1+2 .and.
     $     ktfstringq(dtastk(isp)) .and. ktfstringq(dtastk(isp1+1)))then
        kx%k=merge(ktftrue,ktffalse,
     $       tfsamestringq(dtastk(isp),dtastk(isp1+1)))
        if(iopc .eq. mtfunequal)then
          kx%k=ktftrue-kx%k
        endif
        irtc=0
      else
        kx=tfrelation(isp1,iopc,irtc)
      endif
      return
      end

      function tfupset(k1,k2,kas,irtc) result(kx)
      use tfstk
      implicit none
      type (sad_descriptor) ,intent(in):: k1,k2
      type (sad_descriptor) kx,ki,karg
      type (sad_dlist), pointer :: kl,kli
      type (sad_symbol), pointer :: symi
      type (sad_symdef), pointer :: symd
      integer*8 ,intent(in):: kas
      integer*4 ,intent(out):: irtc
      integer*4 i,isp0,isp1,m,itfmessage
      kx=dxnullo
      if(ktfnonlistq(k1,kl))then
        irtc=itfmessage(999,'General::wrongtype','"Expression"')
        return
      endif
      m=kl%nl
      if(m .le. 0)then
        irtc=itfmessage(999,'General::wrongleng',
     $       '"Expression","longer than 0"')
        return
      endif
      isp0=isp
      isp1=isp0+1
      call tfgetllstk(kl,0,-1)
      karg=kxcompose(isp1)
      LOOP_I: do i=isp1+1,isp
        ki=dtastk(i)
        do while(ktflistq(ki,kli))
          ki=kli%head
        enddo
        if(ktfsymbolqdef(ki%k,symd))then
          if(symd%sym%override .ne. 0)then
            if(symd%sym%gen .lt. 0 .and. symd%sym%gen .ne. -3)then
              cycle LOOP_I
            endif
            if(kas .eq. 0 .or. kas .eq. ktfaddr(ki))then
              call tfdset(k2,symd%upval,kx,karg)
              if(kas .ne. 0)then
                cycle LOOP_I
              endif
            endif
          else
            symi=>tfsydef(symd%sym)
            if(symi%gen .lt. 0 .and. symi%gen .ne. -3)then
              cycle LOOP_I
            endif
            if(kas .eq. 0 .or. kas .eq. ksad_loc(symi%loc))then
              call sym_symdef(symi,symd)
              call tfdset(k2,symd%upval,kx,karg)
              if(kas .ne. 0)then
                cycle LOOP_I
              endif
            endif
          endif
        endif
      enddo LOOP_I
      kx=k2
      isp=isp0
      irtc=0
      return
      end

      function tfoverride(isp1,irtc) result(kx)
      use tfstk
      implicit none
      type (sad_descriptor) kx
      type (sad_dlist), pointer :: kli
      integer*8 ki,k1
      integer*4 ,intent(in):: isp1
      integer*4 ,intent(out):: irtc
      integer*4 isp0,isp2,n,itfmessage,isp3,isp4,i,j
      if(isp1 .eq. isp)then
        kx=dxnulll
        irtc=0
        return
      elseif(isp1+1 .eq. isp)then
        if(ktastk(isp) .eq. ktfoper+mtfnull)then
          kx=dxnulll
          irtc=0
          return
        endif
      endif
      kx=dxnullo
      isp=isp+1
      isp0=isp
      do i=isp1+1,isp0-1
        ki=ktastk(i)
        if(ktflistq(ki,kli))then
          k1=kli%head%k
          if(k1 .eq. ktfoper+mtflist)then
            call tfgetllstkall(kli)
          elseif(k1 .eq. ktfoper+mtfrule .or.
     $           k1 .eq. ktfoper+mtfruledelayed)then
            isp=isp+1
            ktastk(isp)=ki
          else
            go to 9000
          endif
        else
          isp=isp+1
          ktastk(isp)=ki
        endif
      enddo
      isp2=isp
      do i=isp0+1,isp2
        isp=isp+1
        dtastk(isp)=merge(kli%dbody(1),dtastk(i),
     $       ktflistq(ktastk(i),kli))
      enddo
      n=isp-isp2
      isp3=isp
      call tfsortl(ktastk(isp2-3),.false.,n,2,ktfref,.true.,irtc)
      if(irtc .ne. 0)then
        isp=isp0-1
        return
      endif
      isp4=isp
      do i=1,n
        j=int(ktastk(isp3+i))
        if(j .ne. 0 .and.
     $       ktastk(isp0+j) .ne. ktfoper+mtfnull)then
          isp=isp+1
          ktastk(isp)=ktastk(isp0+j)
        endif
      enddo
      kx=kxmakelist(isp4)
      isp=isp0-1
      irtc=0
      return
 9000 irtc=itfmessage(9,'General::wrongtype',
     $     '"List, Rule, Symbol, String, Real"')
      isp=isp0-1
      return
      end

      function tfslot(iopc,kls,ref,irtc) result(kx)
      use tfstk
      use eeval
      implicit none
      type (sad_descriptor) kx,ka
      type (sad_dlist) ,intent(in):: kls
      type (sad_symbol), pointer :: sym
      type (sad_namtbl), pointer :: nam
      integer*8 kaa
      integer*4 ,intent(in):: iopc
      integer*4 ,intent(out):: irtc
      integer*4 ind,nc,isp1,isps,
     $     itfmessage,ns,ipf0,naf0,ls,isp2,itfmessagestr
      real*8 ffval,vx
      character*256 name
      character*12 inds
      logical*4 ,intent(in):: ref
      logical*4 exist
      kx=dxnullo
      ns=kls%nl
      if(ns .gt. 1)then
        irtc=itfmessage(9,'General::narg','"0 or 1"')
        return
      endif
      if(ns .eq. 0)then
        ind=1
      else
        ka=kls%dbody(1)
        if(ktfoperq(ka,kaa))then
          if(kaa .eq. mtfnull)then
            ind=1
          else
            irtc=itfmessage(999,'General::invop',' ')
            return
          endif
        elseif(ktfrealq(ka,ind))then
        elseif(ktfsymbolq(ka,sym) .and. iopc .eq. mtfslot)then
          call sym_namtbl(sym,nam)
          nc=nam%str%nch+1
          name(2:nc)=nam%str%str(1:nc-1)
          name(1:1)='#'
          call capita(name(1:nc))
          vx=ffval(name(1:nc),exist)
          if(exist)then
            kx=dfromr(vx)
            irtc=0
          else
            irtc=itfmessage(999,'FFS::undef','"element"')
          endif
          return
        else
          irtc=itfmessage(999,'General::wrongtype','"Number or symbol"')
          return
        endif
        if(ind .lt. 0)then
          ind=napuref+ind+1
        endif
      endif
      isps=ipurefp+ind
      if(iopc .eq. mtfslot)then
        if(ipurefp .eq. 0 .or. ind .le. 0 .or. ind .gt. napuref)then
          call strfromil(ind,inds,ls)
          irtc=itfmessagestr(999,'General::slot',
     $         '#'//inds(:ls))
          return
        endif
        kx=dtastk(isps)
      else
        if(ipurefp .eq. 0 .or. ind .le. 0)then
          call strfromil(ind,inds,ls)
          irtc=itfmessagestr(999,'General::slot',
     $         '##'//inds(:ls))
          return
        endif
        isp1=isp
        isp2=ipurefp+napuref
        kx=tfsequence(isps-1,isp2)
      endif
      ipf0=ipurefp
      naf0=napuref
      ipurefp=itastk(1,ipf0+naf0+1)
      napuref=itastk(2,ipf0+naf0+1)
c      write(*,*)'tfslot ',ipf0,naf0,ipurefp,napuref
c      call tfdebugprint(kx,'puref',1)
      kx=tfeeval(kx,ref,irtc)
      ipurefp=ipf0
      napuref=naf0
      return
      end

      end module
