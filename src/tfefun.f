      module eexpr

      contains
      recursive function tfeexpr(k1,k,iopc1) result(ke)
      use tfstk
      implicit none
      type (sad_descriptor) ,intent(in):: k1,k
      type (sad_descriptor) ke,kx,ky,ky1,k2,kx2,tftake,tfcompose,
     $     tfjoin2
      type (sad_dlist), pointer :: listy,list1,listi,klx,kl1
      type (sad_rlist), pointer :: kle
      type (sad_symbol), pointer :: sym
      type (sad_string), pointer :: str
      type (sad_pat), pointer :: kp1
      type (sad_complex), pointer :: cx
      integer*8 ks
      integer*4 ,intent(in):: iopc1
      integer*4 m,irtc,i,isp0,iopc,itfcanonicalorder
      real*8 vx1,vy,v2,vx,x
      logical*4 eval
      iopc=iopc1
      ky=k
      select case (iopc)
      case(mtfplus,mtftimes)
        if(.not. tfnumberq(k1) .and. tfnumberq(ky))then
          kx=ky
          ky=k1
        else
          kx=k1
        endif
        if(ktflistq(ky,listy))then
          if(listy%head%k .eq. ktfoper+iopc)then
            ke=tfjoine(kx,ky,irtc)
            if(irtc .ne. -1)then
              return
            endif
            irtc=0
            m=listy%nl
            if(tfnumberq(listy%dbody(1)))then
              if(tfnumberq(kx))then
                ky1=tfcmplx(listy%dbody(1),kx,iopc,irtc)
                if(irtc .ne. 0)then
                  ke%k=ktfoper+mtfnull
                  return
                endif
                if(ktfrealq(ky1))then
                  if(iopc .eq. mtfplus)then
                    if(ky1%k .eq. 0)then
                      ke=merge(listy%dbody(2),
     $                     tftake(ky,dfromr(dble(-m+1)),
     $                     .true.,.false.,irtc),m.eq. 2)
                      return
                    endif
                  else
                    if(ky1%k .eq. ktftrue)then
                      ke=merge(listy%dbody(2),
     $                     tftake(ky,dfromr(dble(-m+1)),
     $                     .true.,.false.,irtc),m .eq. 2)
                      return
                    elseif(ky1%k .eq. 0 .and.
     $                     redmath%value%k .ne. 0)then
                      ke%k=0
                      return
                    endif
                  endif
                endif
                listy=>tfclonelist(listy)
                call tfreplist(listy,1,ky1,eval)
                ke=sad_descr(listy)
              else
                ke=tfinsertsort(listy,kx)
              endif
              return
            endif
            if(ktfrealq(kx))then
              if(iopc .eq. mtfplus .and. kx%k .eq. 0)then
                ke=ky
                return
              elseif(iopc .eq. mtftimes)then
                if(kx%k .eq. ktftrue)then
                  ke=ky
                  return
                elseif(kx%k .eq. 0 .and.redmath%value%k .ne. 0)then
                  ke%k=0
                  return
                endif
              endif
            endif
            ke=tfinsertsort(listy,kx)
            return
          endif
        endif
        if(ktflistq(kx,klx))then
          if(klx%head%k .eq. ktfoper+iopc)then
            ke=tfinsertsort(klx,ky)
            return
          endif
        elseif(ktfrealq(kx))then
          if(iopc .eq. mtfplus .and. kx%k .eq. 0)then
            ke=ky
            return
          elseif(iopc .eq. mtftimes)then
            if(kx%k .eq. ktftrue)then
              ke=ky
              return
            elseif(kx%k .eq. 0 .and. redmath%value%k .ne. 0)then
              ke%k=0
              return
            endif
          endif
        endif
        isp=isp+3
        ktastk(isp-2)=ktfoper+iopc
        if(itfcanonicalorder(ky,kx) .ge. 0)then
          dtastk(isp-1)=kx
          dtastk(isp  )=ky
        else
          dtastk(isp-1)=ky
          dtastk(isp  )=kx
        endif
        ke=kxcompose(isp-2)
        isp=isp-3
        return
      case (mtfnot)
        if(ktflistq(ky,listy))then
          if(listy%head%k .eq. ktfoper+mtfnot)then
            ke=listy%dbody(1)
            return
          endif
        endif
        go to 5000
      case (mtfslot,mtfslotseq)
        if(ktfrealq(ky,vy))then
          ks=int8(vy)
c          write(*,*)'tfeexpr-slot ',vy,ks
          if(dble(ks) .ne. vy) then
            go to 5000
          endif
          if(ks .le. 0.or. ks .gt. nslots)then
            go to 5000
          endif
        elseif(ky%k .ne. ktfoper+mtfnull)then
          go to 5000
        else
          ks=1
        endif
        ke=merge(dlist(iaxslotnull+(ks-1)*2),
     $       dlist(iaxslotnull+(ks-1)*2+1),iopc .eq. mtfslot)
        return
      case (mtfflag)
        go to 5000
      case (mtfcomp,mtfconcat,mtfand,mtfor,mtfalt,mtfmessagename)
        if(ktflistq(k1,list1))then
          if(list1%head%k .eq. ktfoper+iopc)then
            if(ktflistq(ky,listy))then
              if(listy%head%k .eq. ktfoper+iopc)then
                ke=tfjoin2(k1,ky,.false.,irtc)
                return
              endif
            endif
            ke=tfappend(k1,ky,.false.,0,irtc)
            return
          endif
        endif
        if(ktflistq(ky,listy))then
          if(listy%head%k .eq. ktfoper+iopc)then
            ke=tfappend(ky,k1,.false.,1,irtc)
            return
          endif
        endif
      case (mtfset,mtfpower)
        if(ktflistq(ky,listy))then
          if(listy%head%k .eq. ktfoper+iopc)then
            ke=tfappend(ky,k1,.false.,1,irtc)
            return
          endif
        endif
        if(iopc .eq. mtfpower)then
          if(ktfrealq(ky))then
            if(ky%k .eq. ktftrue)then
              ke=k1
              return
            elseif(ky%k .eq. 0 .and. redmath%value%k .ne. 0)then
              ke%k=ktftrue
              return
            endif
          endif
          if(ktflistq(k1,list1))then
            if(list1%head%k .eq. ktfoper+mtfpower)then
              m=list1%nl
              if(m .eq. 1)then
                kx=ky
              else
                k2=merge(list1%dbody(2),
     $               tftake(k1,dfromr(dble(-m+1)),.true.,.false.,irtc),
     $               m .eq. 2)
                if(ktfrealq(k2,v2) .and. ktfrealq(ky,vy))then
                  kx=dfromr(v2*vy)
                else
                  kx=tfeexpr(k2,ky,mtftimes)
                endif
              endif
              ky=list1%dbody(1)
              ke=tfeexpr(ky,kx,mtfpower)
              return
            elseif(list1%head%k .eq. ktfoper+mtftimes
     $             .and. ktfrealq(ky))then
              vx1=1.d0
              kx2%k=ktftrue
              m=list1%nl
              isp0=isp
              isp=isp+1
              ktastk(isp)=ktfoper+mtftimes
              do i=1,m
                isp=isp+1
                dtastk(isp)=list1%dbody(i)
                if(ktfrealq(ktastk(isp)))then
                  vx1=vx1*rtastk(isp)
                  isp=isp-1
                elseif(ktflistq(ktastk(isp),listi))then
                  if(listi%head%k .eq. ktfoper+mtfcomplex)then
                    kx=tfeexpr(dtastk(isp),ky,mtfpower)
                    if(ktfrealq(kx,vx))then
                      vx1=vx1*vx
                    else
                      kx2=tfeexpr(kx2,kx,mtftimes)
                    endif
                    isp=isp-1
                  endif
                endif
              enddo
              if(isp .ne. isp0+m+1)then
                if(isp .gt. isp0+2)then
                  kx=tfcompose(isp0+1,ktfoper+mtftimes,irtc)
                elseif(isp .eq. isp0+2)then
                  kx=dtastk(isp)
                else
                  kx%k=ktftrue
                endif
                isp=isp0
                ke=tfeexpr(kx,ky,mtfpower)
                if(ktfrealq(kx2,v2))then
                  vx1=vx1*v2
                else
                  k2=tfeexpr(kx2,ke,mtftimes)
                endif
                if(vx1 .ne. 1.d0)then
                  kx=tfcmplx(sad_descr(vx1),ky,mtfpower,irtc)
                  ke=tfeexpr(kx,ke,mtftimes)
                endif
                return
              else
                isp=isp0
              endif
            endif
          elseif(ktfrealq(k1))then
            if(k1%k .eq. ktftrue .and. redmath%value%k .ne. 0)then
              ke%k=ktftrue
              return
            endif
          endif
        endif
      case (mtfreplace,mtfreplacerepeated)
        if(ktflistq(k1,kl1))then
          if(kl1%head%k .eq. ktfoper+iopc)then
            ke=tfappend(k1,ky,.false.,0,irtc)
            return
          endif
        endif
      case (mtfcolon)
        if(ktfsymbolq(k1,sym))then
          call sym_symstr(sym,str)
          ke=kxpalocb(str%str,str%nch,ky,transfer(ktfref,k))
          return
        elseif(ktfpatq(k1,kp1))then
          kp1%default=dtfcopy(ky)
          ke%k=ktfpat+ktfaddrd(k1)
          return
        endif
      case (mtfrevpower)
        ke=tfeexpr(k,k1,mtfpower)
        return
      case (mtfatt)
c        call tfdebugprint(k1,'eexpr',1)
c        call tfdebugprint(ky,'@',1)
        if(ktfnonrealq(k1))then
          if(ktfrealq(k,x))then
            ke=kxavaloc(-1,1,kle)
            kle%rbody(1)=x
            kle%head=dtfcopy(k1)
            return
          elseif((ktfsymbolq(ky) .or. ktfoperq(ky)) .and.
     $           (ktfsymbolq(k1) .or. ktflistq(k1)) .or.
     $           ktfpatq(ky))then
            go to 4900
          elseif(ktflistq(k1,kl1) .and.
     $           kl1%head%k .eq. ktfoper+mtfatt)then
            isp=isp+1
            isp0=isp
            ktastk(isp0)=ktfoper+mtfatt
            call tfgetllstkall(kl1)
            isp=isp+1
            dtastk(isp)=ky
            ke=kxcompose(isp0)
            isp=isp0-1
            return
          else
            isp=isp+2
            dtastk(isp-1)=k1
            dtastk(isp  )=ky
            ke=kxcompose(isp-1)
            isp=isp-2
            return
          endif
        endif
      case (mtfcomplex)
        if(ktfrealq(k))then
          if(k%k .eq. 0)then
            ke=k1
            return
          endif
        elseif(tfcomplexq(k,cx))then
          ky=tfeexpr(cx%dbody(1),k1,mtfplus)
          ke=tfeexpr(ky,cx%dbody(2),mtfcomplex)
          return
        elseif(tfcomplexq(k1,cx))then
          ky=tfeexpr(cx%dbody(2),k,mtfplus)
          ke=tfeexpr(cx%dbody(1),ky,mtfcomplex)
          return
        endif
      case (mtffun)
        if(k%k .eq. ktfoper+mtfnull)then
          ke=kxpfaloc(k1)
          return
        endif
      end select
 4900 isp=isp+3
      ktastk(isp-2)=ktfoper+iopc
      dtastk(isp-1)=k1
      dtastk(isp  )=ky
c      call tfdebugprint(ktfoper+iopc,'tfeexpr',1)
c      call tfdebugprint(k1,'tfeexpr',1)
c      call tfdebugprint(ky,'tfeexpr',1)
c      write(*,*)isp
      ke=kxcompose(isp-2)
      isp=isp-3
      return
 5000 if(ktfrealq(ky))then
        ke=kxavaloc(-1,1,kle)
        kle%dbody(1)=ky
      else
        ke=kxadaloc(-1,1,klx)
        call descr_rlist(ke,kle)
        klx%dbody(1)=dtfcopy(ky)
      endif
      kle%head%k=ktfoper+iopc
      return
      end function

      function tfinsertsort(kl,ki) result(kx)
      use tfstk
      implicit none
      type (sad_descriptor) kx
      type (sad_descriptor) ,intent(in):: ki
      type (sad_dlist) ,intent(inout):: kl
      integer*4 isp0,isp1,isp2,i,ispm,itfcanonicalorder,isp3
      isp0=isp
      call tfgetllstkall(kl)
      isp1=isp0+1
      isp2=isp+1
      do while(isp2 .gt. isp1)
        ispm=isp1+(isp2-isp1)/2
        i=itfcanonicalorder(ki,dtastk(ispm))
        if(i .gt. 0)then
          isp1=ispm+1
        elseif(i .eq. 0)then
          isp1=ispm
          isp2=ispm
        else
          isp2=ispm
        endif
      enddo
      isp2=isp
      isp3=isp+isp1-isp0-1
      ktastk(isp+1:isp3)=ktastk(isp0+1:isp1-1)
      isp=isp3+1
      dtastk(isp)=ki
      ktastk(isp+1:isp+isp2-isp1+1)=ktastk(isp1:isp2)
      isp=isp+isp2-isp1+1
      kx=kxcrelistm(isp-isp2,ktastk(isp2+1:isp),kl%head)
      isp=isp0
      return
      end function

      function tfappend(kl,k,eval,mode,irtc) result(kx)
      use tfstk
      use eeval
      implicit none
      type (sad_descriptor) kx
      type (sad_descriptor) ,intent(in):: kl,k
      type (sad_dlist), pointer :: list,listx
      integer*4 ,intent(out):: irtc
      integer*4 ,intent(in):: mode
      integer*4 m,itfmessage,i
      logical*4 ,intent(in):: eval
      logical*4 ev
      if(.not. ktflistq(kl,list))then
        kx=dxnull
        irtc=itfmessage(9,'General::wrongtype',
     $       '"List or composition for #1"')
        return
      endif
      ev=eval .and. list%head%k .ne. ktfoper+mtflist .and.
     $     list%head%k .ne. ktfoper+mtfalt .and.
     $     list%head%k .ne. ktfoper+mtfnull
      m=list%nl
      call loc_sad(ktaaloc(-1,m+1),listx)
      listx%attr=list%attr
      if(ktfreallistq(list))then
        if(mode .eq. 0)then
          listx%dbody(1:m)=list%dbody(1:m)
          if(ktfrealq(k))then
            listx%dbody(m+1)=k
          else
            listx%dbody(m+1)=dtfcopy(k)
            listx%attr=ior(listx%attr,lnonreallist)
          endif
        else
          listx%dbody(2:m+1)=list%dbody(1:m)
          if(ktfrealq(k))then
            listx%dbody(1)=k
          else
            listx%dbody(1)=dtfcopy(k)
            listx%attr=ior(listx%attr,lnonreallist)
          endif
        endif
      else
        if(mode .eq. 0)then
          do i=1,m
            listx%dbody(i)=dtfcopy(list%dbody(i))
          enddo
          listx%dbody(m+1)=dtfcopy(k)
        else
          do i=1,m
            listx%dbody(i+1)=dtfcopy(list%dbody(i))
          enddo
          listx%dbody(1)=dtfcopy(k)
        endif
      endif
      listx%head=dtfcopy(list%head)
      if(iand(list%attr,kconstarg) .ne. 0)then
        if(.not. tfconstq(k%k))then
          listx%attr=ior(listx%attr-kconstarg,knoconstarg+lnoconstlist)
        endif
      endif
      if(ev)then
        kx=tfleval(listx,.true.,irtc)
      else
        kx=sad_descr(listx)
        irtc=0
      endif
      return
      end function

      function tfnot(isp1,iopc,irtc) result(kx)
      use tfstk
      implicit none
      type (sad_descriptor) kx
      integer*4 isp1,irtc,iopc,itfmessage
      if(isp .ne. isp1+1)then
        kx=dxnullo
        irtc=itfmessage(9,'General::narg','"1"')
        return
      endif
      if(tfnumberq(dtastk(isp)))then
        kx=tfcmplx(dfromr(0.d0),dtastk(isp),iopc,irtc)
      elseif(tflistq(dtastk(isp)))then
        kx=tfearray(dfromr(0.d0),dtastk(isp),iopc,irtc)
      else
        kx=tfeexpr(dfromr(0.d0),dtastk(isp),iopc)
        irtc=0
      endif
      return
      end

      function tfplus(isp1,iopc,irtc) result(kx)
      use tfstk
      implicit none
      type (sad_descriptor) kx,k1,k,ki,tfecmplxl
      type (sad_dlist), pointer :: klx
      integer*4 ,intent(in):: isp1,iopc
      integer*4 ,intent(out):: irtc
      integer*4 i,narg
      real*8 v,v1,vx,vi
      kx=dxnullo
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
            kx=tfcmplx(k1,k,iopc,irtc)
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
              kx=tfcmplx(k1,ki,iopc,irtc)
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

      function tfpower(isp1,irtc) result(kx)
      use tfstk
      implicit none
      type (sad_descriptor) kx,k1,ki,tfecmplxl
      integer*4 ,intent(in):: isp1
      integer*4 ,intent(out):: irtc
      integer*4 i,itfmessage
      if(isp1 .eq. isp)then
        kx=dxnullo
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
          kx=tfcmplx(ki,k1,mtfpower,irtc)
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

      function tfrevpower(isp1,irtc) result(kx)
      use tfstk
      implicit none
      type (sad_descriptor) kx,k1,ki,tfecmplxl
      integer*4 ,intent(in):: isp1
      integer*4 ,intent(out):: irtc
      integer*4 i,itfmessage
      if(isp1 .eq. isp)then
        kx=dxnullo
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
          kx=tfcmplx(ki,k1,mtfpower,irtc)
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

      function tfequal(isp1,iopc,irtc) result(kx)
      use tfstk
      implicit none
      type (sad_descriptor) kx
      integer*4 ,intent(in):: isp1,iopc
      integer*4 ,intent(out):: irtc
      integer*4 itfmessage
      if(isp .lt. isp1+2)then
        kx=dxnullo
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

      recursive function tfrelation(isp1,iopc,irtc) result(kx)
      use tfstk
      implicit none
      type (sad_descriptor) kx
      integer*4 ,intent(in):: isp1,iopc
      integer*4 ,intent(out):: irtc
      integer*4 itfmessage,isp0,k
      kx=dxnullo
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
          kx=tfcmplx(dtastk(isp1+1),dtastk(isp),iopc,irtc)
        elseif(tflistq(dtastk(isp1+1))
     $         .or. tflistq(dtastk(isp)))then
          kx=tfearray(dtastk(isp1+1),dtastk(isp),iopc,irtc)
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
          kx=tfrelation(isp0,iopc,irtc)
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

      function tfecmplx1(k1,k2i,iopc1,c,d) result(kxi)
      use tfstk
      implicit none
      type (sad_descriptor) kxi,tfecmplxl
      type (sad_descriptor) ,intent(in):: k1,k2i
      integer*4 ,intent(in):: iopc1
      logical*4 ,intent(inout):: c,d
      if(tflistq(k2i))then
        kxi=dtfcopy(tfecmplxl(k1,k2i,iopc1))
        d=.true.
        c=c .and. tfconstq(kxi%k)
        kxi=dtfcopy(kxi)
      else
        kxi=tfeexpr(k1,k2i,iopc1)
        if(ktfnonrealq(kxi))then
          c=c .and. tfconstq(kxi%k)
          d=.true.
          kxi=dtfcopy(kxi)                  
        endif
      endif
      return
      end function

      function tfcmplx(k1,k2,iopc,irtc) result(kx)
      use tfstk
      implicit none
      type (sad_descriptor) kx
      type (sad_descriptor) ,intent(in):: k1,k2
      type (sad_descriptor) tfenum,tfeval1
      type (sad_complex), pointer :: cx1,cx2
      integer*8 ki1,ki2
      integer*4 ,intent(out):: irtc
      integer*4 ,intent(in):: iopc
      real*8 v1,v2
      complex*16 c1
      irtc=0
      if(iopc .eq. mtfnot)then
        if(ktfrealq(k2))then
          if(k2%k .eq. 0)then
            kx%k=ktftrue
          else
            kx%k=0
          endif
        else
          go to 8000
        endif
        return
      endif
      if(ktfrealq(k1,v1))then
        if(ktfrealq(k2,v2))then
          if(iopc .eq. mtfcomplex)then
            if(k2%k .ne. 0)then
              kx=kxcalocv(-1,v1,v2)
            else
              kx=k1
            endif
          elseif(iopc .eq. mtfpower)then
            if(k1%k .lt. 0)then
              ki2=int(v2)
              if(ki2 .ne. v2)then 
                c1=dcmplx(v1,0.d0)**v2
                if(imag(c1) .eq. 0.d0)then
                  kx=sad_descr(dble(c1))
                else
                  kx=kxcalocv(-1,dble(c1),imag(c1))
                endif
              else
                kx=dfromr(v1**ki2)
              endif
              return
            endif
            kx=tfenum(k1%x(1),k2%x(1),mtfpower,irtc)
          elseif(iopc .eq. mtfrevpower)then
            if(k2%k .lt. 0)then
              ki1=int(v1)
              if(ki1 .ne. v1)then 
                c1=dcmplx(v2,0.d0)**v1
                if(imag(c1) .eq. 0.d0)then
                  kx=sad_descr(dble(c1))
                else
                  kx=kxcalocv(-1,dble(c1),imag(c1))
                endif
              else
                kx=sad_descr(v2**ki1)
              endif
              return
            endif
            kx=tfenum(v2,v1,mtfpower,irtc)
          else
            kx=tfenum(k1%x(1),k2%x(1),iopc,irtc)
          endif
          return
         elseif(tfcomplexq(k2,cx2))then
          kx=tfcmplxmath(dcmplx(v1,0.d0),cx2%cx(1),iopc,irtc)
          return
        else
          go to 8000
        endif
      elseif(tfcomplexq(k1,cx1))then
        if(ktfrealq(k2,v2))then
          kx=tfcmplxmath(cx1%cx(1),dcmplx(v2,0.d0),iopc,irtc)
          return
        elseif(tfcomplexq(k2,cx2))then
          kx=tfcmplxmath(cx1%cx(1),cx2%cx(1),iopc,irtc)
          return
        endif
      endif
 8000 kx=tfeval1(k1,k2,iopc,irtc)
      return
      end

      function tfcmplxmath(c1,c2,iopc1,irtc) result(kx)
      use tfstk
      implicit none
      type (sad_descriptor) kx
      integer*4 ,intent(in):: iopc1
      integer*4 ,intent(out):: irtc
      complex*16 ,intent(in):: c1,c2
      complex*16 cx
      if(iopc1 .gt. mtfunequal .and. iopc1 .ne. mtfcomplex)then
        kx=dxnullo
        irtc=-1
      else
        cx=tfcmplxmathv(c1,c2,iopc1)
        kx=merge(dfromr(dble(cx)),
     $       kxcalocv(-1,dble(cx),imag(cx)),
     $       imag(cx) .eq. 0.d0)
        irtc=0
      endif
      return
      end

      complex*16 function tfcmplxmathv(c1,c2,iopc1)
      use tfstk
      implicit none
      integer*4 ,intent(in):: iopc1
      integer*8 i1,i2
      complex*16 ,intent(in):: c1,c2
      select case(iopc1)
      case (mtfneg)
        tfcmplxmathv=-c2
      case (mtfinv)
        tfcmplxmathv=1.d0/c2
      case (mtfplus)
        tfcmplxmathv=c1+c2
      case (mtftimes)
        tfcmplxmathv=c1*c2
      case (mtfrevpower)
        if(imag(c1) .eq. 0.d0)then
          i1=int8(c1)
          tfcmplxmathv=merge(merge(1.d0/c2,c2**i1,i1 .eq. -1),
     $         c2**dble(c1),i1 .eq. dble(c1))
        else
          tfcmplxmathv=c2**c1
        endif
      case(mtfpower)
        if(imag(c2) .eq. 0.d0)then
          i2=int8(c2)
          tfcmplxmathv=merge(merge(1.d0/c1,merge((1.d0,0.d0),c1**i2,
     $         i2 .eq. 0 .and. redmath%value%k .ne. 0),
     $         i2 .eq. -1),c1**dble(c2),i2 .eq. dble(c2))
        else
          tfcmplxmathv=c1**c2
        endif
      case (mtfequal)
        tfcmplxmathv=merge(1.d0,0.d0,c1 .eq. c2)
      case (mtfunequal)
        tfcmplxmathv=merge(1.d0,0.d0,c1 .ne. c2)
      case (mtfcomplex)
        tfcmplxmathv=c1+dcmplx(-imag(c2),dble(c2))
      case default
        tfcmplxmathv=0.d0
      end select
      return
      end

      function tfearray(k1,k,iopc1,irtc) result(kx)
      use tfstk
      implicit none
      type (sad_descriptor) ,intent(in):: k,k1
      type (sad_descriptor) kx
      type (sad_descriptor) ky,tfdot,tfecmplxl
      type (sad_dlist), pointer :: kl,kl1
      integer*4 ,intent(out):: irtc
      integer*4 ne,ne1,i,iopc1,isp0
      logical*4 list1,list
c     begin initialize for preventing compiler warning
c     end   initialize for preventing compiler warning
c$$$      if(tfmatrixqd(k1,kl1))then
c$$$        if(tfmatrixqd(k,kl2))then
c$$$          call tfematrix(kl1,kl2,kx,iopc1,irtc)
c$$$        else
c$$$c     call tfematrix1(kl1,k,kx,iopc1,irtc)
c$$$        endif
c$$$        if(irtc .ne. 0)then
c$$$          go to 101
c$$$        endif
c$$$        return
c$$$      elseif(tfmatrixqd(k,kl2))then
c$$$c     call tfmatrix2(k1,kl2,kx,iopc1,irtc)
c$$$        if(irtc .ne. 0)then
c$$$          exit
c$$$        endif
c$$$        return
c$$$      endif
c      call tfdebugprint(k1,'earray',1)
c      call tfdebugprint(k,'and',1)
c      write(*,*)'with ',iopc1
      irtc=0
      do
        select case(iopc1)
        case (mtfplus:mtfnot,mtfcomplex)
          kx=tfecmplxl(k1,k,iopc1)
          return
        case (mtfsame:mtfunsame)
          exit
        end select
        if(ktflistq(k1,kl1))then
          if(tfcomplexq(k1))then
            ne1=0
            list1=.false.
          else
            if(tfexprq(k1))then
              exit
            endif
            ne1=kl1%nl
            list1=.true.
          endif
        else
          ne1=0
          list1=.false.
        endif
        if(ktflistq(k,kl))then
          if(tfcomplexq(k))then
            list=.false.
            ne=0
          else
            if(tfexprq(k))then
              exit
            endif
            if(iopc1 .eq. mtfdot)then
              kx=tfdot(k1,k,irtc)
              return
            endif
            ne=kl%nl
            list=.true.
            if(list1 .and. ne .ne. ne1)then
              if(iopc1 .eq. mtfequal)then
                kx%k=0
              elseif(iopc1 .eq. mtfunequal)then
                kx%k=ktftrue
              else
                exit
              endif
              return
            endif
          endif
        else
          list=.false.
          ne=0
        endif
        if(list1)then
          if(list)then
            if(iopc1 .eq. mtfequal)then
              kx%k=ktftrue
              do i=1,ne
                ky=tfcmplx(kl1%dbody(i),kl%dbody(i),iopc1,irtc)
                if(irtc .ne. 0)then
                  return
                endif
                if(ky%k .eq. 0)then
                  kx%k=0
                  return
                elseif(.not. ktfrealq(ky))then
                  exit
                endif
              enddo
              return
            elseif(iopc1 .eq. mtfunequal)then
              kx%k=0
              do i=1,ne
                ky=tfcmplx(kl1%dbody(i),kl%dbody(i),iopc1,irtc)
                if(irtc .ne. 0)then
                  return
                endif
                if(ky%k .eq. ktftrue)then
                  kx%k=ktftrue
                  return
                elseif(.not. ktfrealq(ky))then
                  exit
                endif
              enddo
              return
            else
              isp0=isp
              do i=1,ne
                isp=isp+1
                dtastk(isp)=tfcmplx(kl1%dbody(i),kl%dbody(i),iopc1,irtc)
                if(irtc .ne. 0)then
                  kx=dxnullo
                  isp=isp0
                  return
                endif
              enddo
              kx=kxmakelist(isp0)
              isp=isp0
            endif
          else
            if(iopc1 .eq. mtfequal .or. iopc1 .eq. mtfunequal)then
              exit
            endif
            isp0=isp
            do i=1,ne
              isp=isp+1
              dtastk(isp)=tfcmplx(kl1%dbody(i),k,iopc1,irtc)
              if(irtc .ne. 0)then
                kx=dxnullo
                isp=isp0
                return
              endif
            enddo
            kx=kxmakelist(isp0)
            isp=isp0
          endif
        else
          if(iopc1 .eq. mtfequal .or. iopc1 .eq. mtfunequal)then
            exit
          endif
          isp0=isp
          do i=1,ne
            isp=isp+1
            dtastk(isp)=tfcmplx(k1,kl%dbody(i),iopc1,irtc)
            if(irtc .ne. 0)then
              kx=dxnullo
              isp=isp0
              return
            endif
          enddo
          kx=kxmakelist(isp0)
          isp=isp0
        endif
        return
      enddo
      kx=tfeexpr(k1,k,iopc1)
      return
      end

      function tfjoine(k1,k2,irtc) result(kx)
      use tfstk
      implicit none
      type (sad_descriptor) kx,tfjoin2
      type (sad_descriptor) ,intent(in):: k1,k2
      type (sad_descriptor) k10,k20,ky1
      type (sad_dlist), pointer ::kl1,kl2
      integer*4 ,intent(out):: irtc
      integer*4 ma1,ma2,m,iopc,isp1
      if(ktfnonlistq(k1,kl1) .or. ktfnonlistq(k2,kl2)
     $     .or. .not. tfsameheadq(k1,k2))then
        kx=dxnullo
        irtc=-1
        return
      endif
      irtc=0
      iopc=int(ktfaddr(kl1%head%k))
      if(iopc .eq. mtfplus .or. iopc .eq. mtftimes)then
        if(.not. tfnumberq(kl1%dbody(1)))then
          kx=tfjoin2(k2,k1,.false.,irtc)
          go to 1000
        endif
        if(.not. tfnumberq(kl2%dbody(1)))then
          kx=tfjoin2(k1,k2,.false.,irtc)
          go to 1000
        endif
        ma1=kl1%nl
        ma2=kl2%nl
        m=ma1+ma2-1
        k10=kl1%dbody(1)
        k20=kl2%dbody(1)
        ky1=tfcmplx(k10,k20,iopc,irtc)
        if(irtc .ne. 0)then
          kx=dxnullo
          return
        endif
        if(ktfrealq(ky1))then
          if(iopc .eq. mtfplus)then
            if(ky1%k .eq. 0)then
              m=m-1
            endif
          else
            if(ky1%k .eq. 0)then
              kx%k=0
              return
            elseif(ky1%k .eq. ktftrue)then
              m=m-1
            endif
          endif
        endif
        if(m .eq. 1)then
          kx=kl2%dbody(2)
          return
        endif
        isp1=isp
        if(m .eq. ma1+ma2-2)then
          call tfgetllstk(kl1,2,-1)
          call tfgetllstk(kl2,2,-1)
        else
          isp=isp+1
          dtastk(isp)=ky1
          call tfgetllstk(kl1,2,-1)
          call tfgetllstk(kl2,2,-1)
        endif
        kx=kxcrelistm(isp-isp1,ktastk(isp1+1:isp),k_descr(ktfoper+iopc))
        isp=isp1
      else
        kx=tfjoin2(k1,k2,.false.,irtc)
        return
      endif
 1000 isp=isp+1
      dtastk(isp)=kx
      call tfsort(isp-1,kx,0,irtc)
      isp=isp-1
      return
      end

      end module eexpr

      module efun

      contains
      recursive function tfefunref(isp1,upvalue,irtc)
     $     result(kx)
      use tfstk
      use tfmem
      use tfshare
      use tfcsi,only:lfno
      use mathfun
      use eexpr
      use funs
      use eeval
      use readbuf
      use gammaf
      implicit none
      type (sad_descriptor) kx,k1,k,kh,tfmodule,tfjoin2,
     $     tfsolvemember,tftable,tfefun1,tfeval1,tfeintf,
     $     tfeintf2,tfget,tftake,tfeval1to,tfmap,tfminmax,
     $     tfgetcommandline,tfreplacepart,tfpart,tfwrite,
     $     tftemporaryname,tfmapfile,tfunmapfile,tfcmplxf,
     $     tfgaussiancoulomb,tfgaussiancoulombu
      type (sad_dlist), pointer :: kl,kl1,klx,klh
      type (sad_symbol), pointer :: sym1
      type (sad_symdef), pointer :: symd
      type (sad_string), pointer :: str
      integer*8 ka1,ka,kax,km,kop
      integer*4 ,intent(in):: isp1
      integer*4 ,intent(out):: irtc
      integer*4 i,id,narg,nc,isp2,itfdepth,
     $     itfgetrecl,ltr0,iaf,itfopenwrite,itfmessage,
     $     itfopenappend,itfopcode,isp0,nsize
      real*8 vx,v1,v,f
      character*3 opcx
      character*8 char8
      logical*4 ,intent(in):: upvalue
      logical*4 euv,tfnearlysameqf,rep
      type (sad_descriptor), save :: ilog2
      data ilog2%k /0/

c     DOUBLE specific math intrinsic function
c     from Fortran77    77   77   77
      intrinsic dabs,dsqrt,dexp,dlog
c     from Fortran77   77   77    77    77    77     77
      intrinsic dsin,dcos,dtan,dasin,dacos,datan,datan2
c     from Fortran 77    77    77
      intrinsic dsinh,dcosh,dtanh
c     from Fortran 08    08  <-- somehow cannot pass as an argument @ gfortran 5
c      intrinsic derf,derfc

c     DOUBLE COMPLEX specific math intrinsic function
c     from Fortran EX    EX     EX    EX    EX
      intrinsic cdsin,cdcos,cdsqrt,cdexp,cdlog

c     DOUBLE specific `t-prefix' math function proxy
c     for Fortran2008 generic function
      real*8   tasinh,tacosh,tatanh
      external tasinh,tacosh,tatanh

c     DOUBLE COMPLEX specific `t-prefix' math function proxy
c     for vendor extended math intrinsic function
      complex*16 tctan
      external   tctan

c     DOUBLE COMPLEX specific math function implemented by SAD

c      real*8 , external:: aloggamma1,factorial,gammaq,gammap,inverseerf,
c     $     productlog,gamma0,ferf,ferfc
c      complex*16, external:: cloggamma1,cfactorial,cerfc,cerf,
c     $     cproductlog
      if(upvalue)then
        LOOP_I: do i=isp1+1,isp
          k1=dtastk(i)
 12       if(ktflistq(k1,kl1))then
            klh=>kl1
            k1=kl1%head
            do while(ktflistq(k1,kl1))
              k1=kl1%head
            enddo
            if(k1%k .eq. ktfoper+mtfatt .or.
     $           k1%k .eq. ktfoper+mtfslot)then
              k1=tfsolvemember(klh,rep,irtc)
              if(irtc .eq. 0)then
                dtastk(i)=k1
                go to 12
              elseif(irtc .gt. 0)then
                return
              endif
            endif
          endif
          if(ktfsymbolq(k1,sym1))then
            if(sym1%override .eq. 0)then
              sym1=>tfsydef(sym1)
            endif
            call sym_symdef(sym1,symd)
            if(symd%upval .ne. 0)then
              call tfdeval(isp1,sad_loc(symd%sym%loc),kx,
     $             0,.false.,euv,irtc)
              if(euv)then
                return
              endif
            endif
          endif
        enddo LOOP_I
      endif
      k1=dtastk(isp1)
      narg=isp-isp1
      if(ktfoperq(k1,ka1))then
        if(ka1 .le. mtfend)then
          go to 6000
        endif
        if(narg .eq. 0)then
          narg=1
          isp=isp+1
          ktastk(isp)=ktfoper+mtfnull
        endif
        k=dtastk(isp)
        id=iget_fun_id(ka1)
c        write(*,*)'tfefunref-1 ',id
        irtc=-1
        go to (110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 
     $         210, 220, 230, 240, 240, 260, 270, 280, 290, 300,
     $         310, 320, 330, 340, 350, 360, 370, 380, 380, 400,
     $         410, 420, 430, 440, 440, 460, 470, 480, 490, 500,
     $         510, 520, 530, 540, 550, 560, 570, 580, 590, 600,
c              Sin  Cos  Tan  Sinh Cosh Tanh Exp  Log  Atan Det
c              Sqrt Flor Ceil Min  Max  Mod  StrL Arg  Sign Leng
c              Dims RplP ASin ACos ASh  ACh  ATh  Tabl Do   Attr
c              Peek Abs  Revs Modl Blck StrR SwiC Fltn If   Take
c              Selc Whil Join Apnd Ppnd Clr  Prot Unpr Drop MpAt
c
     $         610, 620, 630, 640, 650, 660, 670, 680, 690, 700,
     $         700, 720, 730, 730, 730, 760, 770, 780, 790, 800,
     $         810, 820, 830, 840, 850, 860, 870, 880, 890, 900,
     $         910, 920, 930, 940, 950, 960, 970, 980, 990,1000,
     $        1010,1020,1030,1040,1050,1060,1070,1080,1090,1100,
c              Innr Trps SglV DiaM LinS IdtN Eigs Oper Posi Sum
c              Prod Rang Re   Im   Conj ToSt Dpth Levl Writ Get
c              OpnW OpnA Clos Flsh Prnt WStr Retn Head RLiQ Pttn
c              Thrw Ctch Thrd SetA MpId FrCh ToCh CmpQ Tr   SvSM
c              Stch Sort Uni1 Ordr MChk Scan Iden TimU NumQ VecQ
c
     $        1110,1120,1130,1140,1150,1160,1170,1180,1190,1200,
     $        1210,1220,1230,1240,1250,1260,1270,1280,1290,1300,
     $        1310,1320,1330,1340,1350,1360,1370,1380,1390,1390,
     $        1290,1420,1430,1440,1450,1460,1470,1480,1490,1500,
     $        1510,1520,1530,1540,1540,1560,1570,1580,1590,1600,
c              AtmQ Outr MatQ TrcP Defi READ Ints Cmpl Roun IErf
c              FrDt ToDt ToIS ReaS OpnR ToEx StrM StrP ToUp Brek
c              Cont Goto Four IFou Chek Whic MapF UnmF GetU GetG
c              ToLo Unev Case DelC Vect BDPi Nams GbCl LgGm Fact
c              With WhiC Ovrr AppT PreT FndR GamR GmRP Erf  Erfc
c
     $        1610,1620, 760,1640,1650,1660,1670,1680,1690,1700,
     $        1710,1720,1730,1740,1750,1760,1770,1770,1770,1770,
     $        1810,1820,1830,1840,1850,1860,1870,1870,1870,1900,
     $        1910,1920,1930,1940,1950,1960,1970,1980,1980,1980,
     $        2010,2020,2030,2040,2050,2060,2070,2080,2080,2100,
c              Fit  Symb SyNm Extr Read Skip TmpN Exit StrF Rstr
c              MM   Shrt $SOT Dir  SDir Wait BesJ BesY BesI BesK
c              BasF StrT S2S  Even OddQ DatS Inst Delt FlAt Repl
c              SetE Spl$ FInd SCnt SCnP ToCt Ctxt BAnd BOr  BXor
c              RepM MSca StdF Abrt ChkA RelH NaNQ MapT ScaT Last
c
     $        2100,2100,2100,2140,2150,2160,2170,2180,2190,2200,
     $        2210,2220,2230,2240,2250,2260,2270,2280,2290,2300,
     $        2310,2320,2330,2340,2350,2360,2370,2380
c              Frst Scnd Thrd ObjS PrdL GauC ClMO MAll Dupl GCLn
c              Seek DigQ LetQ ReaQ NSmQ OpSh RdSh WrSh ShSz FBuQ
c              GaCU GaCF Rest RRt1 Diff Gam0 **** XSin
     $       ),id
c
        if(id .gt. 0)then
          if(id .le. 2000)then
            kx=tfefun1(isp1,id,.true.,irtc)
            go to 6900
          elseif(id .le. 3000) then
            call tfefun2(isp1,id,k,kx,irtc)
            go to 6900
          elseif(id .le. 4000) then
            call tfefun3ep(isp1,id,kx,irtc)
            go to 6900
          elseif(id .le. 5000) then
c            call tfdebugprint(ktastk(isp1),'tfefun-ctbl8',1)
c            call tfdebugprint(ktastk(isp),'with',1)
            call tfefunctbl8(isp1,id,kx,irtc)
c            call tfdebugprint(kx,'result',1)
c            write(*,*)'irtc: ',irtc
            go to 6900
          endif
        endif
        go to 100
 100    write(*,*)'Function implementation error: ',id
        call tfdebugprint(k1,'tfefun',1)
        irtc=0
        go to 7000
 110    if(narg .eq. 1)then
          kx=tfeintf(dsin,cdsin,k,
     $         .true.,-dinfinity,dinfinity,irtc)
        else
          go to 6810
        endif
        go to 6900
 120    if(narg .eq. 1)then
          kx=tfeintf(dcos,cdcos,k,
     $         .true.,-dinfinity,dinfinity,irtc)
        else
          go to 6810
        endif
        go to 6900
 130    if(narg .eq. 1)then
          kx=tfeintf(dtan,tctan,k,
     $         .true.,-dinfinity,dinfinity,irtc)
        else
          go to 6810
        endif
        go to 6900
 140    if(narg .eq. 1)then
          kx=tfeintf(dsinh,tcsinh,k,
     $         .true.,-dinfinity,dinfinity,irtc)
        else
          go to 6810
        endif
        go to 6900
 150    if(narg .eq. 1)then
          kx=tfeintf(dcosh,tccosh,k,
     $         .true.,-dinfinity,dinfinity,irtc)
        else
          go to 6810
        endif
        go to 6900
 160    if(narg .eq. 1)then
          kx=tfeintf(dtanh,tctanh,k,
     $         .true.,-dinfinity,dinfinity,irtc)
        else
          go to 6810
        endif
        go to 6900
 170    if(narg .eq. 1)then
          kx=tfeintf(dexp,cdexp,k,
     $         .true.,-dinfinity,dinfinity,irtc)
        else
          go to 6810
        endif
        go to 6900
 180    if(narg .eq. 1)then
          kx=tfeintf(dlog,cdlog,k,
     $         .true.,0.d0,dinfinity,irtc)
        elseif(narg .eq. 2)then
          if(ilog2%k .eq. 0)then
            ilog2=kxsymbolf('Log2$',5,.true.)
          endif
          dtastk(isp1)=ilog2
          call tfefun(isp1,kx,.true.,.false.,irtc)
          isp=isp1+2
          dtastk(isp1)=k1
        else
          irtc=itfmessage(9,'General::narg','"1 or 2"')
        endif
        go to 6900
 190    if(narg .eq. 1)then
          kx=tfeintf(datan,tcatan,k,
     $         .true.,-dinfinity,dinfinity,irtc)
        elseif(narg .eq. 2)then
          dtastk(isp-1)=tfeintf2(datan2,tcatan2,k,k1,.true.,irtc)
        else
          go to 6812
        endif
        go to 6900
 200    call tfdet(isp1,kx,irtc)
        go to 6900
 210    if(narg .eq. 1)then
          kx=tfeintf(dsqrt,cdsqrt,k,
     $         .true.,0.d0,dinfinity,irtc)
        else
          go to 6810
        endif
        go to 6900
 220    if(narg .eq. 1)then
          kx=tfeintf(tfloor,tcfloor,k,
     $         .true.,-dinfinity,dinfinity,irtc)
        else
          go to 6810
        endif
        go to 6900
 230    if(narg .eq. 1)then
          kx=tfeintf(tceiling,tcceiling,k,
     $         .true.,-dinfinity,dinfinity,irtc)
        else
          go to 6810
        endif
        go to 6900
 240    kx=tfminmax(isp1,id-13,irtc)
        go to 6900
 260    call tfmod(isp1,kx,0,irtc)
        go to 6900
 270    if(narg .eq. 1)then
          if(ktfstringq(k,str))then
            kx=dfromr(dble(str%nch))
            go to 8000
          else
            irtc=itfmessage(9,'General::wrongtype','"Character-string"')
          endif
        else
          go to 6810
        endif
        go to 6900
 280    if(narg .eq. 1)then
          kx=tfeintf(tfarg,tfcarg,k,.false.,-dinfinity,dinfinity,irtc)
        else
          go to 6810
        endif
        go to 6900
 290    if(narg .eq. 1)then
          kx=tfeintf(tfsign,tfcsign,k,.true.,-dinfinity,dinfinity,irtc)
        else
          go to 6810
        endif
        go to 6900
 300    if(narg .ne. 1)then
          go to 6810
        endif
        kx%x(1)=merge(dble(kl%nl),0.d0,ktflistq(k,kl))
        go to 8000
 310    call tfdimensions(isp1,kx,irtc)
        go to 6900
 320    kx=tfreplacepart(isp1,0,irtc)
        go to 6900
 330    if(narg .eq. 1)then
          kx=tfeintf(dasin,tcasin,k,.true.,-1.d0,1.d0,irtc)
        else
          go to 6810
        endif
        go to 6900
 340    if(narg .eq. 1)then
          kx=tfeintf(dacos,tcacos,k,.true.,-1.d0,1.d0,irtc)
        else
          go to 6810
        endif
        go to 6900
 350    if(narg .eq. 1)then
          kx=tfeintf(tasinh,tcasinh,k,.true.,-dinfinity,dinfinity,irtc)
        else
          go to 6810
        endif
        go to 6900
 360    if(narg .eq. 1)then
          kx=tfeintf(tacosh,tcacosh,k,.true.,1.d0,dinfinity,
     $         irtc)
        else
          go to 6810
        endif
        go to 6900
 370    if(narg .eq. 1)then
          kx=tfeintf(tatanh,tcatanh,k,.true.,-1.d0,1.d0,
     $         irtc)
        else
          go to 6810
        endif
        go to 6900
 380    isp2=isp
        kx=tftable(isp1,isp1+2,isp2,29-id,irtc)
        go to 6900
 400    call tfattributes(isp1,kx,irtc)
        go to 6900
 410    if(narg .eq. 2)then
          if(ktfnonrealq(dtastk(isp),f))then
            go to 6810
          endif
        elseif(narg .ne. 1)then
          go to 6810
        else
          f=0.d0
        endif
        if(ktfnonrealq(k))then
          irtc=itfmessage(9,'General::wrongtype','"Real number"')
        else
          ka=int8(rtastk(isp1+1))
          if(f .eq. 0.d0 .and. .not. tfchecklastp(ka))then
            irtc=itfmessage(9,'General::wrongnum',
     $           '"within allocated block"')
          else
            kx=kxadaloc(-1,4,klx)
            klx%rbody(1)=rlist(ka)
            klx%rbody(2)=dble(klist(ka))
            klx%rbody(3)=dble(ilist(1,ka))
            klx%rbody(4)=dble(ilist(2,ka))
            klx%dbody(5)=kxsalocb(0,transfer(klist(ka),char8),8)
            irtc=0
          endif
        endif
        go to 6900
 420    if(narg .eq. 1)then
          kx=tfeintf(dabs,ccdabs,k,
     $         .false.,-dinfinity,dinfinity,irtc)
        else
          go to 6810
        endif
        go to 6900
 430    call tfreverse(isp1,kx,irtc)
        go to 6900
 440    kx=tfmodule(isp1,id .eq. nfunmodule,.true.,irtc)
        go to 6900
 460    call tfstringreplace(isp1,kx,irtc)
        go to 6900
 470    call tfswitchcases(isp1,kx,0,irtc)
        go to 6900
 480    call tfflatten(isp1,kx,irtc)
        go to 6900
 490    call tfiff(isp1,kx,irtc)
        go to 6900
 500    if(narg .ne. 2)then
          go to 6812
        else
          kx=tftake(dtastk(isp-1),k,.true.,.true.,irtc)
        endif
        go to 6900
 510    call tfselect(isp1,kx,irtc)
        go to 6900
 520    call tfwhile(isp1,kx,irtc)
        go to 6900
 530    kx=tfjoin(isp1,.true.,irtc)
        go to 6900
 540    if(narg .eq. 2)then
          kx=tfappend(dtastk(isp1+1),k,.true.,0,irtc)
        elseif(narg .ne. 1)then
          go to 6812
        endif
        go to 6900
 550    if(narg .eq. 2)then
          kx=tfappend(dtastk(isp1+1),k,.true.,1,irtc)
        elseif(narg .ne. 1)then
          go to 6812
        endif
        go to 6900
 560    call tfclear(isp1,kx,irtc)
        go to 6900
 570    call tfprotect(isp1,kx,.true.,irtc)
        go to 6900
 580    call tfprotect(isp1,kx,.false.,irtc)
        go to 6900
 590    if(narg .ne. 2)then
          go to 6812
        else
          kx=tftake(dtastk(isp-1),k,.false.,.true.,irtc)
        endif
        go to 6900
 600    kx=tfreplacepart(isp1,1,irtc)
        go to 6900
 610    if(narg .ne. 4)then
          irtc=itfmessage(9,'General::wrongnum','"4"')
        else
          call tfinner(ktastk(isp-2),ktastk(isp-1),
     $         kx,k,ktastk(isp-3),irtc)
        endif
        go to 6900
 620    if(narg .ne. 1)then
          go to 6810
        else
          call tftranspose(k,kx,irtc)
        endif
        go to 6900
 630    call tfsingularvalues(isp1,kx,irtc)
        go to 6900
 640    if(narg .ne. 1)then
          go to 6810
        else
          call tfdiagonalmatrix(k,kx,irtc)
        endif
        go to 6900
 650    call tflinearsolve(isp1,kx,irtc)
        go to 6900
 660    if(narg .ne. 1)then
          go to 6810
        else
          call tfidentitymatrix(k,kx,irtc)
        endif
        go to 6900
 670    if(narg .ne. 1)then
          go to 6810
        else
          call tfeigensystem(k,kx,irtc)
        endif
        go to 6900
 680    if(narg .ne. 1)then
          go to 6810
        elseif(.not. ktfstringq(k))then
          irtc=itfmessage(9,'General::wrongtype','"Character-string"')
        else
          ka=ktfaddr(k)
          nc=ilist(1,ka)
          if(nc .eq. 0 .or. nc .gt. 3)then
            irtc=itfmessage(9,'General::invop',' ')
            go to 6900
          endif
          opcx=' '
          call tmovb(ilist(1,ka+1),opcx,nc)
          kax=itfopcode(opcx)
          if(kax .lt. 0)then
            irtc=itfmessage(9,'General::invop',' ')
          else
            kx%k=ktfoper+kax
            irtc=0
          endif
        endif
        go to 6900
 690    call tfposition(isp1,kx,0,irtc)
        go to 6900
 700    isp2=isp
        kx=tftable(isp1,isp1+2,isp2,id-58,irtc)
        go to 6900
 720    call tfrange(isp1,kx,irtc)
        go to 6900
 730    if(narg .eq. 1)then
          kx=tfcmplxf(k,id-62,int(ka1))
          irtc=0
          return
        else
          go to 6810
        endif
        go to 6900
 760    call tftostring(isp1,kx,id .eq. 153,irtc)
        go to 6900
 770    if(narg .eq. 1)then
          kx=dfromr(dble(itfdepth(k)))
          go to 8000
        else
          go to 6810
        endif
 780    if(narg .eq. 2)then
          call tflevel(ktastk(isp1+1),k,kx,irtc)
        else
          go to 6812
        endif
        go to 6900
 790    kx=tfwrite(isp1,irtc)
        go to 6900
 800    if(narg .ne. 1)then
          go to 6810
        else
          kx=tfget(k,irtc)
        endif
        go to 6900
 810    if(narg .ne. 1)then
          go to 6810
        else
          vx=itfopenwrite(k,irtc)
        endif
        go to 829
 820    if(narg .ne. 1)then
          go to 6810
        else
          vx=itfopenappend(k,irtc)
        endif
 829    if(irtc .ne. 0)then
          if(vx .eq. -2.d0)then
            call tfaddmessage(' ',0,lfno)
            irtc=0
            kx%k=kxfailed
            return
          endif
          go to 6900
        endif
        kx=dfromr(vx)
        go to 8000
 830    if(narg .ne. 1)then
          go to 6810
        else
          call tfclosef(k,irtc)
        endif
        go to 8200
 840    call tfflush(isp1,kx,irtc)
        go to 6900
 850    call tfprintf(isp1,kx,irtc)
        go to 6900
 860    call tfwritestring(isp1,kx,irtc)
        go to 6900
 870    if(narg .eq. 1)then
          call tfthrow(irtcret,k,irtc)
          return
        else
          go to 6810
        endif
        go to 6900
 880    call tfhead(k,kx)
        go to 8100
 890    if(narg .eq. 1)then
          kx%k=merge(ktftrue,i00,tfreallistq(k))
          irtc=0
        else
          irtc=itfmessage(9,'General::narg','"1"')
        endif
        go to 6900
 900    call tfpartition(isp1,kx,irtc)
        go to 6900
 910    if(narg .eq. 1)then
          call tfthrow(irtcthrow,k,irtc)
          return
        else
          go to 6810
        endif
        go to 6900
 920    call tfcatch(isp1,kx,irtc)
        go to 6900
 930    call tfthread(isp1,kx,0,irtc)
        go to 6900
 940    call tfsetattributes(isp1,kx,irtc)
        go to 6900
 950    kx=tfmap(isp1,4,1,irtc)
        go to 6900
 960    call tffromcharactercode(isp1,kx,irtc)
        go to 6900
 970    call tftocharactercode(isp1,kx,irtc)
        go to 6900
 980    call tfcomplexlistqkf(isp1,kx,irtc)
        go to 6900
 990    call tftr(isp1,kx,irtc)
        go to 6900
c 990    irtc=itfmessage(999,'General::unregister',' ')
c        go to 6900
 1000   call tfsavesharedmap()
        irtc=0
        go to 6900
c 1000   irtc=itfmessage(999,'General::unregister',' ')
c        go to 6900
 1010   call tfswitch(isp1,kx,irtc)
        go to 6900
 1020   call tfsort(isp1,kx,0,irtc)
        go to 6900
 1030   call tfsort(isp1,kx,1,irtc)
        go to 6900
 1040   call tforder(isp1,kx,irtc)
        go to 6900
 1050   call tfmemcheck(isp1,kx,irtc)
        go to 6900
 1060   kx=tfmap(isp1,1,1,irtc)
        go to 6900
 1070   if(narg .eq. 1)then
          kx=k
          irtc=0
          return
        else
          go to 6810
        endif
 1080   if(narg .eq. 1)then
          call cputime(vx,irtc)
          kx=dfromr(vx*1.d-6)
          go to 8000
        else
          go to 6810
        endif
        go to 6900
 1090   if(narg .eq. 1)then
          kx%k=merge(ktftrue,ktffalse,tfnumberq(k))
          go to 8000
        else
          go to 6810
        endif
        go to 6900
 1100   call tfvectorqf(isp1,kx,irtc)
        go to 6900
 1110   if(narg .eq. 1)then
          kx%k=merge(merge(ktftrue,ktffalse,
     $         kl%head%k .eq. ktfoper+mtfcomplex),
     $         ktftrue,ktflistq(k,kl))
          go to 8000
        else
          go to 6810
        endif
        go to 6900
 1120   call tfouter(isp1,kx,irtc)
        go to 6900
 1130   call tfmatchqf(isp1,kx,irtc)
        go to 6900
 1140   if(narg .eq. 1)then
          if(ktflistq(k,kl))then
            if(kl%head%k .ne. ktfoper+mtfcomp)then
              call tfprint1(k,
     $             6,-itfgetrecl(),4,.true.,.true.,irtc)
            endif
          else
            call tfprint1(k,6,-itfgetrecl(),4,.true.,.true.,irtc)
          endif
          ltr0=ltrace
          ltrace=6
          kx=tfeevalref(k,irtc)
          ltrace=ltr0
        else
          go to 6810
        endif
        go to 6900
 1150   call tfdefinition(isp1,kx,irtc)
        go to 6900
 1160   if(narg .ne. 1)then
          go to 6810
        else
          call nfread(k,kx,irtc)
        endif
        go to 6900
 1170   call tfintersection(isp1,kx,0,irtc)
        go to 6900
 1180   call tfintersection(isp1,kx,1,irtc)
        go to 6900
 1190   if(narg .eq. 1)then
          kx=tfeintf(tround,tcround,k,.true.,-dinfinity,dinfinity,irtc)
        else
          go to 6810
        endif
        go to 6900
 1200   if(narg .eq. 1)then
          kx=tfeintf(inverseerf,inverseerf,k,.true.,
     $         -dinfinity,dinfinity,irtc)
        else
          go to 6810
        endif
        go to 6900
 1210   irtc=itfmessage(999,'General::unregister',' ')
        go to 6900
 1220   irtc=itfmessage(999,'General::unregister',' ')
        go to 6900
 1230   call tftoinputstring(isp1,kx,irtc)
        go to 6900
 1240   call tfreadstring(isp1,kx,.false.,.false.,irtc)
        go to 6900
 1250   call tfopenread(isp1,kx,irtc)
        go to 6900
 1260   call tftoexpression(isp1,kx,irtc)
        go to 6900
 1270   call tfstringmatchq(isp1,kx,irtc)
        go to 6900
 1280   call tfstringposition(isp1,kx,irtc)
        go to 6900
 1290   call tftouppercase(isp1,kx,id-119,irtc)
        go to 6900
 1300   continue
 1310   call tfbreak(id-123,narg,kx,irtc)
        go to 6900
 1320   if(narg .eq. 1)then
          call tfthrow(irtcgoto,k,irtc)
          return
        else
          go to 6810
        endif
        go to 6900
 1330   continue
 1340   if(narg .eq. 1)then
          call tffourier(id .eq. 124,k,kx,irtc)
        else
          go to 6810
        endif
        go to 6900
 1350   kx=tfcheck(isp1,irtc)
        go to 6900
 1360   call tfwhich(isp1,kx,irtc)
        go to 6900
 1370   kx=tfmapfile(isp1,irtc)
        go to 6900
 1380   kx=tfunmapfile(isp1,irtc)
        go to 6900
 1390   irtc=itfmessage(999,'General::unregister',' ')
        go to 6900
 1420   kx=tfsequence(isp1,isp)
        irtc=0
        return
 1430   call tfcases(isp1,kx,irtc)
        go to 6900
 1440   call tfposition(isp1,kx,2,irtc)
        go to 6900
 1450   call tfvectorize(isp1,kx,irtc)
        go to 6900
 1460   irtc=itfmessage(999,'General::unregister',' ')
        goto 6900
 1470   call tfnames(isp1,kx,irtc)
        go to 6900
 1480   call tfgarbagecollect(isp1,kx,irtc)
        go to 6900
 1490   if(narg .eq. 1)then
          kx=tfeintf(aloggamma1,cloggamma1,k,
     $         .true.,0.d0,dinfinity,irtc)
        else
          go to 6810
        endif
        go to 6900
 1500   if(narg .eq. 1)then
          kx=tfeintf(factorial,cfactorial,k,
     $         .true.,-1.d0+5.56d-17,dinfinity,irtc)
        else
          go to 6810
        endif
        go to 6900
 1510   call tfwith(isp1,kx,.true.,irtc)
        go to 6900
 1520   call tfswitchcases(isp1,kx,1,irtc)
        go to 6900
 1530   kx=tfoverride(isp1,irtc)
        go to 6900
 1540   call tfappendto(isp1,kx,id-143,irtc)
        go to 6900
 1560   call tffindroot(isp1,kx,irtc)
        go to 6900
 1570   if(narg .eq. 2)then
          k=tfeintf2(gammaq,tfdummy,dtastk(isp-1),k1,.false.,irtc)
        else
          go to 6812
        endif
        go to 6900
 1580   if(narg .eq. 2)then
          k=tfeintf2(gammap,tfdummy,dtastk(isp-1),k1,.false.,irtc)
        else
          go to 6812
        endif
        go to 6900
 1590   if(narg .eq. 1)then
          kx=tfeintf(ferf,cerf,k,.true.,-dinfinity,dinfinity,irtc)
        else
          go to 6810
        endif
        go to 6900
 1600   if(narg .eq. 1)then
          kx=tfeintf(ferfc,cerfc,k,.true.,-dinfinity,dinfinity,irtc)
        else
          go to 6810
        endif
        go to 6900
 1610   call tffit(isp1,kx,irtc)
        go to 6900
 1620   call tfsymbol(isp1,kx,irtc)
        go to 6900
 1640   call tfextract(isp1,kx,irtc)
        go to 6900
 1650   call tfread(isp1,kx,irtc)
        go to 6900
 1660   call tfskip(isp1,kx,irtc)
        go to 6900
 1670   kx=tftemporaryname(isp1,irtc)
        go to 6900
 1680   call tfexit(isp1,kx,irtc)
        go to 6900
 1690   call tfstringfill(isp1,kx,irtc)
        go to 6900
 1700   call tfrestrict(isp1,kx,irtc)
        go to 6900
 1710   kx=tfminmax(isp1,0,irtc)
        go to 6900
 1720   call tfshort(isp1,kx,irtc)
        go to 6900
 1730   call setompnumthreads(isp1,kx,irtc)
        go to 6900
 1740   irtc=itfmessage(999,'General::unregister',' ')
        go to 6900
 1750   irtc=itfmessage(999,'General::unregister',' ')
        go to 6900
 1760   irtc=itfmessage(999,'General::unregister',' ')
        go to 6900
 1770   call tfbessel(isp1,kx,id-167,irtc)
        go to 6900
 1810   call tfbaseform(isp1,kx,irtc)
        go to 6900
 1820   call tfstringtrim(isp1,kx,irtc)
        go to 6900
 1830   call tfstringtostream(isp1,kx,irtc)
        go to 6900
 1840   if(narg .eq. 1)then
          kx=tfeintf(tfevenq,tfcevenq,k,
     $         .true.,-dinfinity,dinfinity,irtc)
        else
          go to 6810
        endif
        go to 6900
 1850   if(narg .eq. 1)then
          kx=tfeintf(tfoddq,tfcoddq,k,
     $         .true.,-dinfinity,dinfinity,irtc)
        else
          go to 6810
        endif
        go to 6900
 1860   irtc=itfmessage(999,'General::unregister',' ')
        go to 6900
 1870   kx=tfreplacepart(isp1,id-175,irtc)
        go to 6900
 1900   if(narg .ne. 2)then
          go to 6812
        endif
        call tfreplace(dtastk(isp1+1),dtastk(isp),
     $       kx,.false.,.true.,.false.,irtc)
        go to 6900
 1910   irtc=itfmessage(999,'General::unregister',' ')
        go to 6900
 1920   call tfspline(isp1,kx,irtc)
        go to 6900
 1930   call tffindindex(isp1,kx,irtc)
        go to 6900
 1940   call tfsetcontext(isp1,kx,irtc)
        go to 6900
 1950   call tfsetcontextpath(isp1,kx,irtc)
        go to 6900
 1960   call tftocontext(isp1,kx,irtc)
        go to 6900
 1970   call tfcontext(isp1,kx,irtc)
        go to 6900
 1980   call tfmod(isp1,kx,id-187,irtc)
        go to 6900
 2010   call tfreplacemember(isp1,kx,irtc)
        go to 6900
 2020   call tfmemberscan(isp1,kx,irtc)
        go to 6900
 2030   call tfstandardform(isp1,kx,irtc)
        go to 6900
 2040   irtc=-7
        if(narg .eq. 1)then
          if(ktfrealq(k,v))then
            irtc=int(min(max(-3.d0,v),-1.d0)-3.d0)
          endif
        endif
        go to 6900
 2050   go to 6900
 2060   call tfreleasehold(isp1,kx,irtc)
        go to 6900
 2070   call tfnanqk(isp1,kx,irtc)
        go to 6900
 2080   call tfthread(isp1,kx,id-197,irtc)
        go to 6900
 2100   call tffirst(isp1,kx,id-201,irtc)
        go to 6900
 2140   call tfobjectsymbol(isp1,kx,irtc)
        go to 6900
 2150   if(narg .eq. 1)then
          kx=tfeintf(productlog,cproductlog,k,
     $         .true.,-exp(-1.d0),dinfinity,irtc)
        else
          go to 6810
        endif
        go to 6900
 2160   kx=tfgaussiancoulomb(isp1,irtc)
        go to 6900
 2170   call tfclearmemberobject(isp1,kx,irtc)
        go to 6900
 2180   call tfmalloc(isp1,kx,irtc)
        go to 6900
 2190   if(narg .eq. 1)then
          kx=k
          if(ktflistq(k))then
            kx=kxcopylist(k)
          endif
          irtc=0
          return
        else
          go to 6810
        endif
 2200   kx=tfgetcommandline(isp1,irtc)
        go to 6900
 2210   call tfseek(isp1,kx,irtc)
        go to 6900
 2220   call tfdigitQ(isp1,kx,irtc)
        go to 6900
 2230   call tfletterQ(isp1,kx,irtc)
        go to 6900
 2240   call tfrealqk(isp1,kx,irtc)
        go to 6900
 2250   if(narg .ne. 4)then
          go to 6814
        endif
        if(ktfnonrealq(dtastk(isp-1)) .or.
     $       ktfnonrealq(dtastk(isp)))then
          irtc=itfmessage(9,'General::wrongarg',
     $         '"$NearlySameQ[a,b,relthre,absthre]"')
          go to 7000
        endif
        irtc=0
        kx%k=merge(ktftrue,ktffalse,
     $       tfnearlysameqf(dtastk(isp1+1),dtastk(isp1+2),
     $       rtastk(isp1+3),rtastk(isp)))
        go to 6900
 2260   kx=tfopenshared(isp1,irtc)
        go to 6900
 2270   kx=tfreadshared(isp1,irtc)
        go to 6900
 2280   kx=tfwriteshared(isp1,irtc)
        go to 6900
 2290   if(narg .ne. 1)then
          irtc=itfmessage(9,'General::narg','"1"')
          go to 7000
        endif
        isp0=isp1+1
        call tfsharedsize(isp0,k,nsize,irtc)
        isp=isp1+1
        if(irtc .ne. 0)then
          go to 7000
        endif
        kx=dfromr(dble(max(0,nsize-2)*8))
        go to 6900
 2300   if(narg .ne. 1)then
           irtc=itfmessage(9,'General::narg','"1"')
           go to 7000
        endif
        irtc=0
        kx%k=merge(ktftrue,ktffalse,ktfoperq(k))
        go to 6900
 2310   kx=tfgaussiancoulombu(isp1,irtc)
        go to 6900
 2320   call tfgaussiancoulombfitted(isp1,kx,irtc)
        go to 6900
 2330   call tfrest(isp1,kx,irtc)
        go to 6900
 2340   call tfrotateright1(isp1,kx,irtc)
        go to 6900
 2350   call tfdifference(isp1,kx,irtc)
        go to 6900
 2360   if(narg .eq. 1)then
          kx=tfeintf(gamma0,tfdummy,k,.false.,-dinfinity,dinfinity,irtc)
        else
          go to 6810
        endif
        go to 6900
 2370   go to 6900
 2380   if(narg .eq. 1)then
          kx=tfeintf(xsin,tcxsin,k,.true.,-dinfinity,dinfinity,irtc)
        else
          go to 6810
        endif
        go to 6900
 8000   continue
 8100   irtc=0
        return
 8200   if(irtc .ne. 0)then
          go to 6900
        endif
        kx%k=ktfoper+mtfnull
        return
 6000   iaf=int(ka1)
c        go to (
c     $       6010,
c     $       7000,7000,6600,7000,6600,7000,6620,6630,6690,6690,
c     $       6680,6680,6680,6680,6700,6700,6610,6100,6100,6510,
c     $       7000,7000,6010,7000,6640,6650,6750,7000,7000,7000,
c     $       7000,6002,6010,6010,6010,6660,6670,6560,6560,6710,
c     $       6010,6740,7000,7000,6200,6090,6520,6530,6540,6010,
c     $       6010,6550,6570,6570,6570,6570,6580,6580,6590,6720,
c     $       6730,6010,7000,7000,6010,7000),
c     $       iaf+1
c            null
c            m    i    +    -    *    /    v    ^    e    n    
c            >    <    g    l    E    N    ~    &&   o    c
c            [    ]    {    }    s    =    C    (    )    ,
c            ;    &    :    r    d    RepA RepR u    U    S    
c            ?    f    #    ##   .    |    M    MA   A    rept 
c            repn ineq AT   SF   TB   DB   INc  Dec  Part @
c            msgn TagS (*   *)   Hold z
c        go to 6001
        select case(iaf)
        case (mtfnull,mtflist,mtfcolon,mtfhold,mtftagset,
     $       mtfrepeated,mtfrepeatednull,mtfpattest,
     $       mtfrule,mtfruledelayed)
          kx=kxcompose(isp1)
          irtc=0
          return
        case (mtfset)
          kx=tfset(isp1,.true.,irtc)
          go to 6900
        case (mtfplus,mtftimes)
          kx=tfplus(isp1,iaf,irtc)
          go to 6900
        case (mtfrevpower)
          kx=tfrevpower(isp1,irtc)
          go to 6900
        case (mtfpower)
          kx=tfpower(isp1,irtc)
          go to 6900
        case (mtfgreater,mtfless,mtfgeq,mtfleq)
          kx=tfrelation(isp1,iaf,irtc)
          go to 6900
        case (mtfinequality)
          call tfinequality(isp1,kx,irtc)
          go to 6900
        case (mtfpart)
          if(ktflistq(ktastk(isp1+1)))then
            kx=tfpart(isp1+1,.true.,irtc)
            if(irtc .eq. 0)then
              kx=tfeevalref(kx,irtc)
            endif
          elseif(ktfsymbolq(ktastk(isp1+1)))then
            irtc=-1
          else
            irtc=itfmessage(9,'General::wrongtype',
     $           '"List or composition"')
          endif
          go to 6900
        case (mtfnot)
          kx=tfnot(isp1,iaf,irtc)
          go to 6900
        case (mtfupset,mtfupsetdelayed)
          if(narg .ne. 2)then
            go to 6812
          endif
          kx=tfupset(dtastk(isp1+1),dtastk(isp),i00,irtc)
          go to 6900
        case (mtffun)
          if(isp .eq. isp1+2)then
            if(ktastk(isp) .eq. ktfoper+mtfnull)then
              kx=kxpfaloc(dtastk(isp-1))
              irtc=0
              return
            endif
          endif
        case (mtfalt,mtfand,mtfor)
          if(narg .lt. 2 .and. iaf .eq. mtfalt)then
            irtc=-1
            go to 6900
          endif
          if(narg .eq. 2)then
            kx=tfeval1(dtastk(isp1+1),dtastk(isp),iaf,irtc)
            go to 6900
          endif
          if(narg .eq. 0)then
            kx=dfromr(dble(mtfor-iaf))
            irtc=0
            return
          endif
          kx=dtastk(isp1+1)
          if(narg .eq. 1)then
            irtc=0
            return
          endif
          do i=isp1+2,isp
            k1=kx
            kx=tfeval1(k1,dtastk(i),iaf,irtc)
            if(irtc .ne. 0)then
              go to 6900
            endif
          enddo
          return
        case (mtfdot)
          if(narg .eq. 2)then
            kx=tfeval1(dtastk(isp1+1),dtastk(isp),iaf,irtc)
            go to 6900
          endif
          if(narg .eq. 0)then
            irtc=itfmessage(9,'General::narg','"1 or more"')
            go to 6900
          endif
          kx=dtastk(isp)
          if(narg .eq. 1)then
            irtc=0
            return
          endif
          do i=isp-1,isp1+1,-1
            k1=kx
            kx=tfeval1(dtastk(i),k1,iaf,irtc)
            if(irtc .ne. 0)then
              go to 6900
            endif
          enddo
          return
        case (mtfconcat)
          call tfstringjoin(isp1,kx,irtc)
          go to 6900
        case (mtfmap)
          kx=tfmap(isp1,3,1,irtc)
          go to 6900
        case (mtfmapall)
          call tfmapall(isp1,kx,irtc)
          go to 6900
        case (mtfapply)
          call tfapply(isp1,kx,irtc)
          go to 6900
        case (mtfaddto,mtfsubtractfrom,mtftimesby,mtfdivideby)
          if(narg .ne. 2)then
            go to 6812
          endif
          kx=tfeval1to(dtastk(isp1+1),dtastk(isp),iaf,.false.,irtc)
          go to 6900
        case (mtfincrement,mtfdecrement)
          v1=merge(1.d0,-1.d0,iaf .eq. mtfincrement)
          if(narg .eq. 1)then
            kx=tfeval1to(dtastk(isp1+1),dfromr(v1),mtfaddto,.true.,irtc)
          elseif(narg .eq. 2 .and.
     $           ktastk(isp1+1) .eq. ktfoper+mtfnull)then
            kx=tfeval1to(dtastk(isp),dfromr(v1),mtfaddto,.false.,irtc)
          else
            go to 6810
          endif
          go to 6900
        case (mtfsetdelayed)
          if(narg .ne. 2)then
            go to 6812
          endif
          kx=tfset(isp1,.true.,irtc)
          go to 6900
        case (mtfreplace)
          kx=tfreplace1(isp1,irtc)
          go to 6900
        case (mtfreplacerepeated)
          kx=tfreplacerepeated1(isp1,irtc)
          go to 6900
        case (mtfequal,mtfunequal)
          kx=tfequal(isp1,iaf,irtc)
          go to 6900
        case (mtfsame,mtfunsame)
          kx=tfsameq1(isp1,iaf,irtc)
          go to 6900
        case (mtfunset)
          if(narg .ne. 1)then
            go to 6810
          endif
          isp=isp1+2
          ktastk(isp)=ktfref
          kx=tfset(isp1,.true.,irtc)
          isp=isp1+1
          kx%k=ktfoper+mtfnull
          go to 6900
        case (mtfatt)
          call tfatt(isp1,kx,.true.,irtc)
          go to 6900
        case (mtfmessagename)
          km=klist(ifunbase+mtfmessagename)
          call tfdeval(isp1,km,kx,1,.false.,euv,irtc)
          go to 6900
        case (mtfflag)
c          write(*,*)'tfefun-mtfflag'
          call tfflagordef(isp1,kx,irtc)
          go to 6900
        case (mtfcomplex)
          if(narg .ne. 2)then
            go to 6812
          endif
          kx=tfeval1(dtastk(isp1+1),dtastk(isp),mtfcomplex,irtc)
          go to 6900
        case default
          if(iaf .le. mtfend)then
            go to 7000
          endif
          irtc=itfmessage(999,'General::invop',' ')
          return
        end select
      elseif(ktfsymbolqdef(k1%k,symd))then
        if(symd%sym%override .ne. 0)then
          if(symd%downval .ne. 0)then
            call tfdeval(isp1,ktfaddr(k1),kx,1,.false.,euv,irtc)
            go to 6900
          else
            go to 6800
          endif
        else
          go to 6800
        endif
      elseif(ktflistq(k1,kl1))then
        kx=k1
        if(ktfoperq(kl1%head,kop))then
          id=iget_fun_id(kop)
          select case (id)
          case (-mtffun)
            kx=tfpuref(isp1,kl1,irtc)
            go to 6900
          case (-mtfnull)
            if(kl1%nl .eq. 0)then
              kx=tfsequence(isp1,isp)
              if(ktflistq(kx,klx))then
                kx=tfleval(klx,.true.,irtc)
              elseif(ktfsymbolq(kx) .or. ktfpatq(kx))then
                kx=tfeevalref(kx,irtc)
              endif
              go to 6900
            elseif(kl1%nl .eq. 1 .and.
     $             ktfnonreallistqo(kl1))then
              kx=kxmakelist(isp1,klx)
              kh=klx%dbody(1)
              klx%head=dtfcopy(kh)
              irtc=0
              return
            endif
          case (-mtflist)
            kx=tfpart(isp1,.true.,irtc)
            if(irtc .eq. 0)then
              if(ktflistq(kx,klx))then
                kx=tfleval(klx,.true.,irtc)
              elseif(ktfsymbolq(kx) .or. ktfpatq(kx))then
                kx=tfeevalref(kx,irtc)
              endif
            endif
            go to 6900
          case (-mtfmap,-mtfapply)
            if(kl1%nl .eq. 1)then
              isp2=isp+1
              dtastk(isp2)=kl1%head
              dtastk(isp2+1)=kl1%dbody(1)
              dtastk(isp2+2:isp2+isp-isp1+1)=dtastk(isp1+1:isp)
              isp=isp+isp2-isp1+1
              kx=tfefunref(isp2,upvalue,irtc)
              isp=isp2
              go to 6900
            endif
          case (nfunappend,nfunprepend,nfuncases,nfundelcases,
     $           nfunselcases,nfundelete,nfunposition,nfunselect,
     $           nfunreppart,nfunextract,nfunswicases)
            if(kl1%nl .eq. 1)then
              isp2=isp+1
              dtastk(isp2)=kl1%head
              dtastk(isp2+1)=dtastk(isp1+1)
              dtastk(isp2+2)=kl1%dbody(1)
              if(isp .gt. isp1+1)then
                dtastk(isp2+3:isp2+isp-isp1+1)=dtastk(isp1+2:isp)
              endif
              isp=isp2+isp-isp1+1
              kx=tfefunref(isp2,upvalue,irtc)
              isp=isp2
              go to 6900
            endif
          case (nfuninsert)
            if(kl1%nl .eq. 2)then
              isp2=isp+1
              dtastk(isp2)=kl1%head
              dtastk(isp2+1)=dtastk(isp1+1)
              dtastk(isp2+2)=kl1%dbody(1)
              dtastk(isp2+3)=kl1%dbody(2)
              isp=isp2+3
              kx=tfefunref(isp2,upvalue,irtc)
              isp=isp2
              go to 6900
            endif
          end select
          go to 6800
        endif
        do while(ktflistq(kx,klx))
          kx=klx%head
        enddo
        if(ktfsymbolqdef(kx%k,symd))then
          if(symd%sym%override .ne. 0 .and. symd%downval .ne. 0)then
            call tfdeval(isp1,ktfaddrd(kx),kx,1,.false.,euv,irtc)
            go to 6900
          endif
        endif
        go to 6800
      elseif(ktfstringq(k1))then
        if(narg .le. 0 .or. narg .gt. 2)then
          go to 6800
        elseif(ktfnonrealq(ktastk(isp)) .or.
     $         ktfnonrealq(ktastk(isp1+1)))then
          go to 6800
        endif
        kx=kxsubstring(k1,isp1+1,isp)
        irtc=0
        return
      else
        go to 6800
      endif
      go to 6900
 6800 irtc=0
      go to 7000
 6810 irtc=itfmessage(9,'General::narg','"1"')
      go to 7000
 6812 irtc=itfmessage(9,'General::narg','"2"')
      go to 7000
 6814 irtc=itfmessage(9,'General::narg','"4"')
      go to 7000
 6900 if(irtc .eq. 0 .or. irtc .lt. -1)then
        return
      elseif(irtc .eq. -1)then
        irtc=0
      endif
 7000 isp=isp1+narg
      kx=kxcrelistm(narg,ktastk(isp1+1:isp),dtastk(isp1))
      if(irtc .gt. 0)then
        if(ierrorprint .ne. 0)then
          call tferrorhandle(kx,irtc)
        else
          call tfdebugprint(kx,'... in',3)
        endif
      elseif(irtc .eq. -1)then
        irtc=0
      endif
      return
      end

      function tfcheck(isp1,irtc) result(kx)
      use tfstk
      use eeval
      use tfcsi
      implicit none
      type (sad_descriptor) kx,kf,kxcheckmessage,kxmessagelist
      type (sad_symdef), pointer,save :: symd
      integer*4 ,intent(in):: isp1
      integer*4 ,intent(out):: irtc
      integer*4 isp0,itgetfpe,itfmessage,narg,irtc1
      data kxmessagelist%k,kxcheckmessage%k /0,0/
      narg=isp-isp1
      if(narg .lt. 2)then
        kx=dxnullo
        irtc=itfmessage(9,'General::narg','"2 or more"')
        return
      endif
      isp0=isp
      rlist(ierrorgen)=0.d0
      kx=tfeevalref(dtastk(isp1+1),irtc)
      isp=isp0
      if(irtc .eq. 0)then
        if(itgetfpe() .ne. 0)then
          call tclrfpe
          irtc=itfmessage(9,'General::fpe','""')
        endif
      endif
      if(irtc .gt. 0)then
        if(ierrorprint .ne. 0)then
          call tfaddmessage(' ',0,icslfno())
        endif
      elseif(irtc .eq. irtcabort .and. rlist(ierrorgen) .eq. 0.d0)then
        irtc=itfmessage(999,'General::abort',' ')
      elseif(irtc .le. irtcret)then
        rlist(ierrorgen)=1.d0
      endif
      if(rlist(ierrorgen) .ne. 0.d0)then
        if(irtc .le. irtcret)then
          call tfcatchreturn(modethrow,kx,irtc)
c          write(*,*)'tfcheck-1 ',modethrow,irtc
        endif
        if(narg .gt. 2)then
          if(kxcheckmessage%k .eq. 0)then
            kxcheckmessage=kxsymbolz('Check$Message',13)
          endif
          dtastk(isp1+1)=kxcheckmessage
          kf=tfefunref(isp1+1,.false.,irtc1)
          if(irtc1 .eq. 0 .and. ktfrealq(kf) .and. kf%k .ne. 0)then
            rlist(ierrorgen)=0.d0
            kx=tfeevalref(dtastk(isp1+2),irtc)
          endif
        else
          if(kxmessagelist%k .eq. 0)then
            kxmessagelist=kxsymbolz('$MessageList',12)
            call descr_sad(kxmessagelist,symd)
          endif
          call tflocald(symd%value)
          symd%value=dtfcopy1(dxnulll)
          rlist(ierrorgen)=0.d0
          kx=tfeevalref(dtastk(isp1+2),irtc)
        endif
      endif
      return
      end

      function tfjoin(isp1,eval,irtc) result(kx)
      use tfstk
      implicit none
      type (sad_descriptor) kx
      type (sad_descriptor) kf
      type (sad_dlist), pointer :: kl1,kli
      integer*4 ,intent(in):: isp1
      integer*4 ,intent(out):: irtc
      integer*4 itfmessage,i,narg,isp0
      logical*4 ,intent(in):: eval
      logical*4 ev
      narg=isp-isp1
      if(ktfnonlistq(ktastk(isp1+1),kl1))then
        go to 9010
      endif
      if(narg .le. 1)then
        if(narg .eq. 1)then
          kx=dtastk(isp1+1)
          irtc=0
        else
          kx=dxnullo
          irtc=itfmessage(9,'General::narg','"1 or more"')
        endif
        return
      endif
      kf=kl1%head
      isp=isp+1
      isp0=isp
      call tfgetllstkall(kl1)
      if(isp .ge. mstk)then
        kx=dxnullo
        irtc=itfmessage(9999,'General::stack','"Join"')
        return
      endif
      do i=isp1+2,isp0-1
        if(ktfnonlistq(ktastk(i),kli))then
          isp=isp0-1
          go to 9000
        endif
        if(.not. tfsameq(kli%head,kf))then
          go to 9100
        endif
        call tfgetllstkall(kli)
        if(isp .ge. mstk)then
          kx=dxnullo
          irtc=itfmessage(9999,'General::stack','"Join"')
          return
        endif
      enddo
      dtastk(isp0)=kf
      ev=eval
      if(ev .and. (kf%k .eq. ktfoper+mtflist .or.
     $     kf%k .eq. ktfoper+mtfalt .or.
     $     kf%k .eq. ktfoper+mtfnull))then
          ev=.false.
      endif
      if(ev)then
        kx=tfefunref(isp0,.true.,irtc)
      else
        kx=kxcompose(isp0)
        irtc=0
      endif
      isp=isp0-1
      return
 9000 isp=isp0-1
 9010 irtc=itfmessage(9,'General::wrongtype',
     $     '"List or composition for all args"')
      kx=dxnullo
      return
 9100 irtc=itfmessage(9,'General::samehead',' ')
      kx=dxnullo
      isp=isp0-1
      return
      end

      end module efun

      subroutine tfefun(isp1,kx,ref,upvalue,irtc)
      use tfstk
      use efun
      implicit none
      type (sad_descriptor) kx
      integer*4 isp1,irtc
      logical*4 ref,upvalue
      if(ref)then
        kx=tfefunref(isp1,upvalue,irtc)
      else
        call tfefundef(isp1,kx,irtc)
      endif
      return
      end

      subroutine tfefunrefc(isp1,kx,irtc)
      use tfstk
      use efun
      implicit none
      type (sad_descriptor) kx
      integer*4 isp1,irtc
      levele=levele+1
      kx=tfefunref(isp1,.true.,irtc)
      call tfconnect(kx,irtc)
      return
      end

      subroutine tfefunrefstk(isp1,isp2,irtc)
      use tfstk
      use efun
      implicit none
      type (sad_descriptor) kx
      type (sad_dlist), pointer :: kl
      integer*4 isp1,irtc,isp2
      kx=tfefunref(isp1,.true.,irtc)
      if(irtc .ne. 0)then
        return
      endif
      if(ktflistq(kx,kl))then
        if(kl%head%k .eq. ktfoper+mtfnull)then
          isp=isp2-1
          call tfgetllstkall(kl)
          return
        endif
      endif
      isp=isp2
      dtastk(isp)=kx
      return
      end

      subroutine tfefundef(isp1,kx,irtc)
      use tfstk
      use funs
      use eeval
      implicit none
      type (sad_descriptor) kx,k1,tfsolvemember,tfefun1
      type (sad_dlist), pointer :: kl,kli,klx
      type (sad_symbol), pointer :: sym
      type (sad_symdef), pointer :: symd
      integer*8 ka1
      integer*4 isp1,irtc,narg,id,itfmessageexp
      logical*4 tfrefq,rep
      irtc=0
      narg=isp-isp1
 1    k1=dtastk(isp1)
      if(ktfoperq(k1,ka1))then
        id=iget_fun_id(ka1)
        if(ka1 .gt. mtfend .and.
     $       id .gt. 999 .and. id .le. 2000)then
          kx=tfefun1(isp1,id,.false.,irtc)
          if(irtc .eq. 0 .and. tfrefq(kx))then
          else
            irtc=itfmessageexp(999,'General::invset',k1)
          endif
          go to 6900
        elseif(ka1 .eq. mtfatt)then
          call tfatt(isp1,kx,.false.,irtc)
          go to 6900
        else
          kx=kxcompose(isp1)
          return
        endif
      elseif(ktfsymbolqdef(k1%k,symd))then
        if(symd%sym%override .ne. 0)then
          if(ktfoperq(symd%value))then
            dtastk(isp1)=symd%value
            go to 1
          endif
        endif
        go to 7000
      elseif(ktflistq(k1,kl))then
        kx=k1
        if(kl%head%k .eq. ktfoper+mtffun)then
          kx=tfpuref(isp1,kl,irtc)
          go to 6900
        elseif(kl%head%k .eq. ktfoper+mtfnull)then
          if(kl%nl .eq. 0)then
            kx=tfsequence(isp1,isp)
            if(ktflistq(kx,klx))then
              kx=tfleval(klx,.false.,irtc)
            elseif(ktfsymbolq(kx,sym))then
              if(sym%override .eq. 0)then
                sym=>tfsydef(sym)
                kx=sad_descr(sym)
              endif
            endif
            go to 6900
          elseif(kl%nl .eq. 1 .and.
     $           ktfnonreallistqo(kl))then
            kx=kxmakelist(isp1,klx)
            klx%head=dtfcopy(kl%head)
            return
          endif
        elseif(kl%head%k .eq. ktfoper+mtflist)then
          go to 6900
        elseif(ktflistq(kl%head,kli))then
          k1=tfsolvemember(kli,rep,irtc)
          if(irtc .eq. 0)then
            dtastk(isp1)=k1
          endif
        endif
        go to 7000
      elseif(ktfstringq(k1))then
        irtc=itfmessageexp(999,'General::invset',k1)
      else
        go to 7000
      endif
 6900 if(irtc .eq. 0 .or. irtc .lt. -1)then
        return
      elseif(irtc .eq. -1)then
        irtc=0
      endif
 7000 isp=isp1+narg
      kx=kxcrelistm(narg,ktastk(isp1+1:isp1+narg),dtastk(isp1))
      if(irtc .gt. 0)then
        call tferrorhandle(kx,irtc)
      elseif(irtc .eq. -1)then
        irtc=0
      endif
      return
      end

c     Type fixed function proxy for generic math function/vendor extension
c     DOUBLE instance of ASINH in Fortran2008
      real*8 function tasinh(x)
      implicit none
c      include 'inc/MATH.inc'
      real*8 x
      tasinh=asinh(x)
      return
      end
c     DOUBLE instance of ACOSH in Fortran2008
      real*8 function tacosh(x)
      implicit none
c      include 'inc/MATH.inc'
      real*8 x
      tacosh=acosh(x)
      return
      end
c     DOUBLE instance of ATANH in Fortran2008
      real*8 function tatanh(x)
      implicit none
c      include 'inc/MATH.inc'
      real*8 x
      tatanh=atanh(x)
      return
      end
c     DOUBLE COMPLEX proxy of ZTAN vendor extension
      complex*16 function tctan(z)
      implicit none
      complex*16 z,ztan
      tctan=ztan(z)
      return
      end

      subroutine tfefundummy
      use mackw
      return
      end

      subroutine tfseval(ks,ns,kh,kx,stk,av,ref,irtc)
      use tfstk
      use tfcode
      use iso_c_binding
      use efun
      use funs
      use eeval
      implicit none
      type (sad_descriptor) kx,kf,kh,kl
      type (sad_funtbl), pointer :: fun
      type (sad_symbol), pointer :: sym
      type (sad_dlist), pointer :: list,klf,kls,kls1
      integer*4 maxlevel
      parameter (maxlevel=2**12)
      integer*8 ks,iaat,kaf
      integer*4 irtc,ns,iaf,isp10,isp11,j
      integer*4 isp1,level,itfgetrecl,
     $     isp0,mf,i1,level1,i,itfmessage,lpw,l,itfdownlevel
      real*8 v
      logical*4 ref,ev,evalh,rep,stk,tfconstlistqo,tfgetseqstk,av
      data level/0/
      level=level+1
      if(level .ge. maxlevel-64)then
        level1=level
        level=0
        irtc=itfmessage(999,'General::deep','""')
        level=level1
        go to 9000
      endif
      call loc_sad(ks,kls)
      isp0=isp
      kf=kh
      evalh=.false.
      kaf=ktfaddr(kf%k)
      select case(kf%k-kaf)
      case (ktfoper)
        call c_f_pointer(c_loc(klist(klist(ifunbase+kaf)-9)),fun)
        iaf=-fun%id
        if(fun%narg .lt. 0)then
          go to 100
        endif
      case (ktfsymbol)
        if(ref)then
          kf=tfsyeval(kf,irtc)
          if(irtc .ne. 0)then
            go to 8000
          endif
        else
          call loc_sym(kaf,sym)
          sym=>tfsydef(sym)
          kf=sad_descr(sym)
        endif
      case (ktflist)
        call loc_dlist(kaf,list)
        if(iand(lconstlist,list%attr) .eq. 0)then
          kf=tfleval(list,ref,irtc)
          if(irtc .ne. 0)then
            go to 8000
          endif
        endif
      case (ktfpat)
        if(ref)then
          kf=tfpateval(kf,irtc)
          if(irtc .ne. 0)then
            go to 8000
          endif
        endif
      end select
      ev=.true.
      iaat=0

      do
        if(ktfoperq(kf%k))then
          kaf=ktfaddr(kf)
          call c_f_pointer(c_loc(klist(klist(ifunbase+kaf)-9)),fun)
          iaat=klist(ifunbase+kaf)+1
          mf=fun%narg
          if(mf .lt. 0)then
            iaf=-fun%id
            evalh=.true.
            exit
          endif
          ev=fun%mapeval(2,1) .eq. -2
          if(fun%mapeval(2,1) .eq. -1)then
            iaat=0
          endif
        elseif(ktfsymbolq(kf%k,sym))then
          if(sym%override .ne. 0)then
            iaat=iand(iattrholdall,sym%attr)
            if(iaat .ne. 0)then
              ev=.false.
              iaat=merge(i00,-iaat,iaat .eq. iattrholdall)
            endif
          endif
        elseif(ktflistq(kf,klf))then
          if(klf%head%k .eq. ktfoper+mtfnull)then
            if(klf%nl .eq. 1)then
              kf=klf%dbody(1)
              cycle
            elseif(klf%nl .eq. 0)then
              kf%k=ktfoper+mtfnull
              evalh=.true.
              exit
            endif
          elseif(ktfsymbolq(klf%head,sym) .and. .not. ref)then
            if(sym%gen .eq. -3 .and. ktfreallistq(klf))then
              ev=.false.
              iaat=0
            endif
          endif
        endif
        isp=isp0+1
        isp1=isp
        go to 3000
      enddo

 100  isp=isp0+1
      isp1=isp
      select case(iaf)
      case(mtfnull,mtflist,mtfrule,mtfrepeated,mtfrepeatednull)
        if(stk)then
          ev=.true.
          go to 3000
        endif
        if(evalh .or. .not. ref)then
          if(tfonstackq(ks) .or. kls%ref .gt. 0)then
            kls1=>tfduplist(kls)
            kls=>kls1
          endif
          kls%head%k=ktfoper+iaf
c          call tfloadlstk(lista,lista)
c          lista%head=ktfoper+kaf
c          call tfstk2l(lista,lista)
          if(ref .and. tfconstlistqo(kls))then
            kx=sad_descr(kls)
            irtc=0
            go to 8000
          endif
        endif
        if(ref)then
          kx=tfevallev(kls,irtc)
        else
          call tfevallstkall(kls,.false.,.false.,irtc)
          if(irtc .eq. 0)then
            levele=levele+1
            go to 6000
          endif
        endif
      case(mtfset)
        rep=tfgetseqstk(ks,ns)
        if(isp .gt. isp1)then
          dtastk(isp)=tfeevalref(dtastk(isp),irtc)
          if(irtc .ne. 0)then
            go to 8000
          endif
        endif
        levele=levele+1
        go to 6000
      case (mtfand)
        do i=1,ns
          isp10=isp
          call tfseqevalstk(kls%dbody(1),ns,i,av,irtc)
          if(irtc .ne. 0)then
            go to 8000
          endif
          isp11=isp
          isp=isp10
          do j=isp10+1,isp11
            if(ktfrealq(ktastk(j),v))then
              if(v==0.d0)then
                kx%k=0
                go to 8000
              endif
            else
              isp=isp+1
              ktastk(isp)=ktastk(j)
            endif
          enddo
        enddo
        if(isp .eq. isp1)then
          kx%k=ktftrue
        elseif(isp .eq. isp1+1)then
          kx=dtastk(isp)
        else
          levele=levele+1
          go to 6000
        endif
      case (mtfor)
        do i=1,ns
          isp10=isp
          call tfseqevalstk(kls%dbody(1),ns,i,av,irtc)
          if(irtc .ne. 0)then
            go to 8000
          endif
          isp11=isp
          isp=isp10
          do j=isp10+1,isp11
            if(ktfrealq(ktastk(j),v))then
              if(v/=0.d0)then
                kx%k=ktftrue
                go to 8000
              endif
            else
              isp=isp+1
              ktastk(isp)=ktastk(j)
            endif
          enddo
        enddo
        if(isp .eq. isp1)then
          kx%k=ktffalse
        elseif(isp .eq. isp1+1)then
          kx=dtastk(isp)
        else
          levele=levele+1
          go to 6000
        endif
      case (mtfpart)
        if(ref .or. ns .eq. 0)then
          ev=.true.
          go to 3000
        endif
        isp=isp+1
        dtastk(isp)=tfeevaldef(kls%dbody(1),irtc)
        if(irtc .ne. 0)then
          go to 8000
        endif
        call tfseqevalstkall(kls%dbody(2),ns-1,av,irtc)
        if(irtc .eq. 0)then
          levele=levele+1
          go to 6000
        endif
      case (mtfslot,mtfslotseq)
        kx=tfslot(iaf,kls,ref,irtc)
      case (mtfcomp)
        if(ns .eq. 0)then
          kx%k=ktfoper+mtfnull
          go to 8000
        endif
        i1=1
        do
          if(ltrace .gt. 0)then
            levele=levele+1
            lpw=min(131,itfgetrecl())
            do i=i1,ns
              call tfprint1(kls%dbody(i),6,-lpw,4,.true.,
     $             .true.,irtc)
              kx=tfeevalref(kls%dbody(i),irtc)
              if(irtc .ne. 0)then
                go to 1320
              endif
            enddo
          else
            levele=levele+1
            irtc=0
            do i=i1,ns
              kx=kls%dbody(i)
              if(ktflistq(kx,list))then
                kx=tfleval(list,.true.,irtc)
              elseif(ktfsymbolq(kx))then
                kx=tfsyeval(kx,irtc)
              elseif(ktfpatq(kx))then
                kx=tfpateval(kx,irtc)
              endif
              if(irtc .ne. 0)then
                go to 1320
              endif
            enddo
          endif
          go to 7000
 1320     if(irtc .gt. irtcret)then
            go to 7000
          endif
          call tfcatchreturn(irtcgoto,kl,irtc)
          l=itfdownlevel()
          if(irtc .ne. 0)then
            exit
          endif
          call tffindlabel(kls,ns,i1,kl)
          if(i1 .le. 0)then
            call tfthrow(irtcgoto,kl,irtc)
            exit
          endif
          i1=i1+1
        enddo
      case (mtffun,mtfpattest,mtftagset,mtfhold)
        rep=tfgetseqstk(ks,ns)
        if(rep .or. stk .or. evalh)then
          levele=levele+1
          go to 6000
        endif
        kx%k=ktflist+ks
      case default
        write(*,*)'tfleval-implementation error: ',kaf,iaf
        call abort
      end select
      go to 8000

 3000 if(levele .ge. maxlevele-32)then
        irtc=itfmessage(999,'General::deep','""')
        go to 8000
      endif
      levele=levele+1
      if(ev)then
        call tfseqevalstkall(kls%dbody(1),ns,av,irtc)
        if(irtc .ne. 0)then
          go to 7000
        endif
      elseif(iaat .eq. 0 .or. av)then
        rep=tfgetseqstk(ks,ns)
      else
        call tfargevalstk(isp1,kls,ns,iaat,mf,.false.,irtc)
        if(irtc .ne. 0)then
          go to 7000
        endif
      endif
 6000 dtastk(isp1)=kf
c      call tfmemcheckprint('seval-efun',.false.,irtc)
c      if(irtc .ne. 0)then
c        call tfdebugprint(kf,'tfseval-efunref-in',3)
c        call tfdebugprint(ktastk(isp1+1),' ',3)
c        call tfdebugprint(ktastk(isp),' ',3)
c      endif
      if(ref)then
        kx=tfefunref(isp1,.true.,irtc)
c          call tfdebugprint(kf,'tfseval-efunref-out',3)
c          call tfdebugprint(ktastk(isp1+1),' ',3)
c          call tfdebugprint(ktastk(isp),' ',3)
      else
        call tfefundef(isp1,kx,irtc)
      endif
 7000 continue
c      call tfdebugprint(kx,'tfseval-connect',3)
      call tfconnect(kx,irtc)
 8000 isp=isp0
 9000 level=max(0,level-1)
      return
      end

      function tfjoin2(k1,k2,eval,irtc) result(kx)
      use tfstk
      use efun
      implicit none
      type (sad_descriptor) kx
      type (sad_descriptor) ,intent(in):: k1,k2
      integer*4 ,intent(out):: irtc
      integer*4 isp1
      logical*4 ,intent(in):: eval
      isp1=isp
      isp=isp+1
      dtastk(isp)=k1
      isp=isp+1
      dtastk(isp)=k2
      kx=tfjoin(isp1,eval,irtc)
      isp=isp1
      return
      end function
