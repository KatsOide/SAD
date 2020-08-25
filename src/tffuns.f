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

      subroutine tfpuref(isp1,kf,kx,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) kx,ki,ka
      type (sad_dlist) kf
      type (sad_dlist), pointer :: kla
      integer*4 isp1,irtc,itfmessage,narg,m,j,i,ipf0,nap0,isp0
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
        call tfeevalref(kf%dbody(1),kx,irtc)
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
        call tfeevalref(kx,kx,irtc)
        isp=isp1+narg
      else
        irtc=itfmessage(9,'General::narg','"1 or 2"')
      endif
      return
      end

      subroutine tfplus(isp1,kx,iopc,irtc)
      use tfstk
      use eexpr
      implicit none
      type (sad_descriptor) kx,k1,k,ki,tfecmplxl
      type (sad_dlist), pointer :: klx
      integer*4 isp1,irtc,i,iopc,narg
      real*8 v,v1,vx,vi
      narg=isp-isp1
      if(narg .eq. 2)then
        k1=dtastk(isp1+1)
        k =dtastk(isp)
        if(ktfrealq(k1,v1) .and. ktfrealq(k,v))then
          kx=merge(dfromr(v1+v),dfromr(v1*v),
     $         iopc .eq. mtfplus)    
          irtc=0
        else
          if(tfnumberq(k) .and. tfnumberq(k1))then
            call tfcmplx(k1,k,kx,iopc,irtc)
          else
            kx=merge(tfecmplxl(k1,k,iopc),tfeexpr(k1,k,iopc),
     $           tflistq(k1) .or. tflistq(k))
            irtc=0
          endif
        endif
        return
      elseif(narg .eq. 0)then
        kx%k=merge(ktffalse,ktftrue,iopc .eq. mtfplus)
      elseif(narg .eq. 1)then
        if(ktfsymbolq(dtastk(isp)) .or.
     $       ktfpatq(dtastk(isp)))then
          kx=kxmakelist(isp1,klx)
          klx%head%k=ktfoper+iopc
        else
          kx=dtastk(isp1+1)
        endif          
        irtc=0
      else
        kx=dtastk(isp1+1)
        irtc=0
        do i=isp1+2,isp
          ki=dtastk(i)
          if(ktfrealq(ki,vi) .and. ktfrealq(kx,vx))then
            kx=merge(dfromr(vx+vi),dfromr(vx*vi),
     $           iopc .eq. mtfplus)    
          else
            k1=kx
            if(tfnumberq(k1) .and. tfnumberq(ki))then
              call tfcmplx(k1,ki,kx,iopc,irtc)
              if(irtc .ne. 0)then
                return
              endif
            elseif(tflistq(k1) .or. tflistq(ki))then
              kx=tfecmplxl(k1,ki,iopc)
              if(irtc .ne. 0)then
                return
              endif
            else
              kx=tfeexpr(k1,ki,iopc)
            endif
          endif
        enddo
      endif
      return
      end

      subroutine tfpower(isp1,kx,irtc)
      use tfstk
      use eexpr
      implicit none
      type (sad_descriptor) kx,k1,ki,tfecmplxl
      integer*4 isp1,irtc,i,itfmessage
      if(isp1 .eq. isp)then
        irtc=itfmessage(9,'General::narg','"1 or more"')
        return
      endif
      kx=dtastk(isp)
      irtc=0
      if(isp .eq. isp1+1)then
        return
      endif
      do i=isp-1,isp1+1,-1
        ki=dtastk(i)
        k1=kx
        if(tfnumberq(k1) .and. tfnumberq(ki))then
          call tfcmplx(ki,k1,kx,mtfpower,irtc)
          if(irtc .ne. 0)then
            return
          endif
        elseif(tflistq(k1) .or. tflistq(ki))then
          kx=tfecmplxl(ki,k1,mtfpower)
          irtc=0
        else
          kx=tfeexpr(ki,k1,mtfpower)
        endif
      enddo
      return
      end

      subroutine tfrevpower(isp1,kx,irtc)
      use tfstk
      use eexpr
      implicit none
      type (sad_descriptor) kx,k1,ki,tfecmplxl
      integer*4 isp1,irtc,i,itfmessage
      if(isp1 .eq. isp)then
        irtc=itfmessage(9,'General::narg','"1 or more"')
        return
      endif
      kx=dtastk(isp1+1)
      irtc=0
      if(isp .eq. isp1+1)then
        return
      endif
      do i=isp1+2,isp
        ki=dtastk(i)
        k1=kx
        if(tfnumberq(k1) .and. tfnumberq(ki))then
          call tfcmplx(ki,k1,kx,mtfpower,irtc)
          if(irtc .ne. 0)then
            return
          endif
        elseif(tflistq(k1) .or. tflistq(ki))then
          kx=tfecmplxl(ki,k1,mtfpower)
          irtc=0
        else
          kx=tfeexpr(ki,k1,mtfpower)
        endif
      enddo
      return
      end

      recursive function tfset(isp1,upvalue,irtc) result(kx)
      use tfstk
      implicit none
      type (sad_descriptor) kx,k1,k2,k110,k11,tfset1
      type (sad_symbol), pointer :: sym
      type (sad_symdef), pointer :: symd
      type (sad_dlist), pointer :: kl11
      integer*4 isp1,irtc,i,itfmessage,isp11,isp0
      logical*4 euv,upvalue
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

      subroutine tfreplace1(isp1,kx,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) kx,k1
      integer*4 isp1,irtc,i,itfmessage
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

      subroutine tfreplacerepeated1(isp1,kx,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) kx
      integer*4 isp1,irtc,i,itfmessage
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

      recursive subroutine tfrelation(isp1,kx,iopc,irtc)
      use tfstk
      use eexpr
      implicit none
      type (sad_descriptor) kx
      integer*4 isp1,irtc,iopc,itfmessage,isp0,k
      if(isp .lt. isp1+2)then
        irtc=itfmessage(9,'General::narg','"2 or more"')
        return
      endif
      if(isp .eq. isp1+2)then
c        call tfdebugprint(dtastk(isp1+1),'tfrelation',1)
c        call tfdebugprint(dtastk(isp),'and',1)
c        write(*,*)'with: ',iopc
        if(tfnumberq(dtastk(isp1+1)) .and.
     $       tfnumberq(dtastk(isp)))then
          call tfcmplx(dtastk(isp1+1),dtastk(isp),kx,iopc,irtc)
        elseif(tflistq(dtastk(isp1+1))
     $         .or. tflistq(dtastk(isp)))then
          call tfearray(dtastk(isp1+1),dtastk(isp),kx,iopc,irtc)
        else
          kx=tfeexpr(dtastk(isp1+1),dtastk(isp),iopc)
          irtc=0
        endif
      else
        isp0=isp
        do k=1,isp0-isp1-1
          ktastk(isp0+1)=ktastk(isp1+k)
          ktastk(isp0+2)=ktastk(isp1+k+1)
          isp=isp0+2
          call tfrelation(isp0,kx,iopc,irtc)
          if(irtc .ne. 0)then
            isp=isp0
            return
          elseif(ktfnonrealq(kx))then
            irtc=-1
            return
          endif
          if(kx%k .eq. 0)then
            isp=isp0
            return
          endif
        enddo
        isp=isp0
      endif
      return
      end

      subroutine tfsameq1(isp1,kx,iopc,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) kx,k,k1
      integer*4 isp1,irtc,iopc,itfmessage
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

      subroutine tfequal(isp1,kx,iopc,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) kx
      integer*4 isp1,irtc,iopc,itfmessage
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
        call tfrelation(isp1,kx,iopc,irtc)
      endif
      return
      end

      subroutine tfnot(isp1,kx,iopc,irtc)
      use tfstk
      use eexpr
      implicit none
      type (sad_descriptor) kx
      integer*4 isp1,irtc,iopc,itfmessage
      if(isp .ne. isp1+1)then
        irtc=itfmessage(9,'General::narg','"1"')
        return
      endif
      if(tfnumberq(dtastk(isp)))then
        call tfcmplx(dfromr(0.d0),dtastk(isp),kx,iopc,irtc)
      elseif(tflistq(dtastk(isp)))then
        call tfearray(dfromr(0.d0),dtastk(isp),kx,iopc,irtc)
      else
        kx=tfeexpr(dfromr(0.d0),dtastk(isp),iopc)
        irtc=0
      endif
      return
      end

      subroutine tfupset(k1,k2,kas,kx,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) k1,k2,kx,ki,karg
      type (sad_dlist), pointer :: kl,kli
      type (sad_symbol), pointer :: symi
      type (sad_symdef), pointer :: symd
      integer*8 kas
      integer*4 irtc,i,isp0,isp1,m,itfmessage
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
