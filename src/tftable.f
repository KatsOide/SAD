      recursive function tftable(isp1,isp2,ispa,mode,irtc)
     $     result(kx)
      use tfstk
      use efun
      implicit none
      integer*4 maxint
      parameter (maxint=2**31-1)
      type (sad_descriptor) kx,kj,kxlistcopied,ke,ki,k1,kl
      type (sad_dlist), pointer :: listi,kle,klj
      type (sad_rlist), pointer :: klr
      type (sad_symbol), pointer :: name
      type (sad_symdef), pointer :: symd
      integer*8 ls
      integer*4 isp1,isp2,irtc,mode,ispa,narg,ns,
     $     itfmessage,itfdownlevel,lv,m,j,isp0,ispb,i
      real*8 x0,x1,xstep,xns,xi,ve
      logical*4 var
      kx=dxnull
      narg=ispa-isp1
      if(narg .lt. 2)then
        irtc=itfmessage(9,'General::narg','"2 or more"')
        return
      endif
      if(levele .ge. maxlevele-100)then
        irtc=itfmessage(999,'General::deep',' ')
        return
      endif
      ki=dtastk(isp2)
      if(ktflistq(ki,listi))then
        if(listi%head%k .eq. ktfoper+mtflist)then
          m=listi%nl
          if(m .eq. 1)then
            call tfeevalref(listi%dbody(1),kl,irtc)
            if(.not. ktfrealq(kl,x1))then
              irtc=itfmessage(9,'General::wrongtype','"Real number"')
              return
            endif
            x0=1.d0
            xstep=1.d0
            nullify(name)
          elseif(m .gt. 4)then
            go to 9500
          else
            if(ktfreallistq(listi))then
              go to 9500
            endif
            k1=listi%dbody(1)
            if(.not. ktfsymbolq(k1,name))then
              go to 9500
            endif
            call tfeevalref(listi%dbody(2),k1,irtc)
            if(irtc .ne. 0)then
              return
            endif
            if(ktfnonrealq(k1,x1))then
              go to 9500
            endif
            if(m .eq. 2)then
              x0=1.d0
              xstep=1.d0
            else
              x0=x1
              call tfeevalref(listi%dbody(3),k1,irtc)
              if(irtc .ne. 0)then
                return
              endif
              if(ktfnonrealq(k1,x1))then
                go to 9500
              endif
              if(m .eq. 3)then
                xstep=1.d0
              else
                call tfeevalref(listi%dbody(4),k1,irtc)
                if(irtc .ne. 0)then
                  return
                endif
                if(ktfnonrealq(k1,xstep))then
                  go to 9500
                endif
              endif
            endif
          endif
        else
          go to 9500
        endif
      else
        go to 9500
      endif
      if(xstep .eq. 0.d0)then
        irtc=itfmessage(9,'General::wrongnum','"non-zero"')
        return
      endif
      if(associated(name))then
        call descr_symdef(kxnaloc1(name%gen,name%loc),symd)
        call tflocald(symd%value)
        symd%value=dfromr(x0)
        xi=0.d0
        var=.true.
      else
        nullify(symd)
        var=.false.
      endif
      if(x1 .eq. dinfinity .or. x0 .eq. -dinfinity)then
        if(xstep .gt. 0.d0)then
          ls=maxint
        else
          ls=0
        endif
      elseif(x1 .eq. -dinfinity .or. x0 .eq. dinfinity)then
        if(xstep .gt. 0.d0)then
          ls=0
        else
          ls=maxint
        endif
      else
        xns=(x1-x0)/xstep*(1.d0+1.d-11)
        ls=max(-1,min(int8(xns),maxint-1))+1
      endif
      if(ls .gt. maxint)then
        irtc=itfmessage(9,'General::wrongval',
     $       '"# of elements","less than 2^31"')
        return
      endif
      ns=int(ls)
      ispb=isp+1
      if(mode .ne. 0)then
        isp=ispb
      endif
      if(isp2 .lt. ispa)then
        do j=1,ns
          kj=tftable(isp1,isp2+1,ispa,mode,irtc)
          if(irtc .ne. 0)then
            go to 9000
          endif
          if(mode .ne. 0)then
            isp=isp+1
            dtastk(isp)=dtfcopy(kj)
          endif
          if(var)then
            xi=xi+1.d0
            symd%value=dfromr(xstep*xi+x0)
          endif
        enddo
      else
        ke=dtfcopy(dtastk(isp1+1))
        if(var)then
          if(mode .eq. 0)then
            if(ktflistq(ke,kle))then
              do j=1,ns
                levele=levele+1
                call tfleval(kle,kj,.true.,irtc)
                lv=itfdownlevel()
                if(irtc .ne. 0)then
                  if(irtc .eq. -2)then
                    irtc=0
                  else
                    go to 9200
                  endif
                endif
                xi=xi+1.d0
                symd%value=dfromr(xstep*xi+x0)
              enddo
            else
              do j=1,ns
                levele=levele+1
                call tfeevalref(ke,kj,irtc)
                lv=itfdownlevel()
                if(irtc .ne. 0)then
                  if(irtc .eq. -2)then
                    irtc=0
                  else
                    go to 9200
                  endif
                endif
                xi=xi+1.d0
                symd%value=dfromr(xstep*xi+x0)
              enddo
            endif
          else
            if(ktfrealq(ke,ve))then
              irtc=0
              if(mode .eq. 1)then
                kx=kxavaloc(-1,ns,klr)
                klr%rbody(1:ns)=ve
              elseif(mode .eq. 2)then
                kx=dfromr(ve*ns)
              elseif(mode .eq. 3)then
                kx=dfromr(ve**ns)
              endif
              isp=ispb-1
              call tfdelete(symd,.true.,.false.)
              return
            else
              if(ktflistq(ke,kle))then
                do j=1,ns
                  levele=levele+1
                  call tfleval(kle,kj,.true.,irtc)
                  if(irtc .ne. 0)then
                    if(irtc .eq. -2)then
                      irtc=0
                      go to 9010
                    elseif(irtc .eq. -3)then
                      lv=itfdownlevel()
                      irtc=0
                      exit
                    endif
                    go to 9100
                  endif
                  kj=dtfcopy(kj)
                  lv=itfdownlevel()
                  if(ktflistq(kj,klj))then
                    if(klj%head%k .eq. ktfoper+mtfnull)then
                      isp0=isp
                      call tfgetllstkall(klj)
                      call tflocal1d(kj)
c                      call incr1i(ilist(1,ktaobj(ktastk(isp0+1:isp))-1))
                      call ktfcopym(ktastk(isp0+1:isp))
c                      ilist(1,ktaobj(ktastk(isp0+1:isp))-1)=
c     $                     ilist(1,ktaobj(ktastk(isp0+1:isp))-1)+1
c                      do i=isp0+1,isp
c                        ktastk(i)=ktfcopy(ktastk(i))
c                      enddo
                    else
                      isp=isp+1
                      dtastk(isp)=kj
                    endif
                  else
                    isp=isp+1
                    dtastk(isp)=kj
                  endif
 9010             xi=xi+1.d0
                  symd%value=dfromr(xstep*xi+x0)
                enddo
              else
                do j=1,ns
                  levele=levele+1
                  call tfeevalref(ke,kj,irtc)
                  if(irtc .ne. 0)then
                    if(irtc .eq. -2)then
                      irtc=0
                      go to 9020
                    elseif(irtc .eq. -3)then
                      lv=itfdownlevel()
                      irtc=0
                      exit
                    endif
                    go to 9100
                  endif
                  kj=dtfcopy(kj)
                  lv=itfdownlevel()
                  if(ktflistq(kj,klj))then
                    if(klj%head%k .eq. ktfoper+mtfnull)then
                      isp0=isp
                      call tfgetllstkall(klj)
                      call tflocal1d(kj)
                      call ktfcopym(ktastk(isp0+1:isp))
c                      do i=isp0+1,isp
c                        ktastk(i)=ktfcopy(ktastk(i))
c                      enddo
                    else
                      isp=isp+1
                      dtastk(isp)=kj
                    endif
                  else
                    isp=isp+1
                    dtastk(isp)=kj
                  endif
 9020             xi=xi+1.d0
                  symd%value=dfromr(xstep*xi+x0)
                enddo
              endif
            endif
          endif
        else
          if(mode .eq. 0)then
            if(ktflistq(ke,kle))then
              do j=1,ns
                levele=levele+1
                call tfleval(kle,kj,.true.,irtc)
                lv=itfdownlevel()
                if(irtc .ne. 0)then
                  if(irtc .eq. -2)then
                    irtc=0
                  else
                    go to 9200
                  endif
                endif
              enddo
            else
              do j=1,ns
                levele=levele+1
                call tfeevalref(ke,kj,irtc)
                lv=itfdownlevel()
                if(irtc .ne. 0)then
                  if(irtc .eq. -2)then
                    irtc=0
                  else
                    go to 9200
                  endif
                endif
              enddo
            endif
          else
            if(ktfrealq(ke,ve))then
              irtc=0
              if(mode .eq. 1)then
                kx=kxavaloc(-1,ns,klr)
                klr%rbody(1:ns)=ve
              elseif(mode .eq. 2)then
                kx=dfromr(ve*ns)
              elseif(mode .eq. 3)then
                kx=dfromr(ve**ns)
              endif
              isp=ispb-1
              return
            else
              if(ktflistq(ke,kle))then
                do j=1,ns
                  levele=levele+1
                  call tfleval(kle,kj,.true.,irtc)
                  if(irtc .ne. 0)then
                    if(irtc .eq. -2)then
                      irtc=0
                      cycle
                    elseif(irtc .eq. -3)then
                      irtc=0
                      lv=itfdownlevel()
                      exit
                    endif
                    go to 9100
                  endif
                  kj=dtfcopy(kj)
                  lv=itfdownlevel()
                  if(ktflistq(kj,klj))then
                    if(klj%head%k .eq. ktfoper+mtfnull)then
                      isp0=isp
                      call tfgetllstkall(klj)
                      call tflocal1d(kj)
                      call ktfcopym(ktastk(isp0+1:isp))
c                      do i=isp0+1,isp
c                        ktastk(i)=ktfcopy(ktastk(i))
c                      enddo
                    else
                      isp=isp+1
                      dtastk(isp)=kj
                    endif
                  else
                    isp=isp+1
                    dtastk(isp)=kj
                  endif
                enddo
              else
                do j=1,ns
                  levele=levele+1
                  call tfeevalref(ke,kj,irtc)
                  if(irtc .ne. 0)then
                    if(irtc .eq. -2)then
                      irtc=0
                      cycle
                    elseif(irtc .eq. -3)then
                      irtc=0
                      lv=itfdownlevel()
                      exit
                    endif
                    go to 9100
                  endif
                  kj=dtfcopy(kj)
                  lv=itfdownlevel()
                  if(ktflistq(kj,klj))then
                    if(klj%head%k .eq. ktfoper+mtfnull)then
                      isp0=isp
                      call tfgetllstkall(klj)
                      call tflocal1d(kj)
                      call ktfcopym(ktastk(isp0+1:isp))
c                      do i=isp0+1,isp
c                        ktastk(i)=ktfcopy(ktastk(i))
c                      enddo
                    else
                      isp=isp+1
                      dtastk(isp)=kj
                    endif
                  else
                    isp=isp+1
                    dtastk(isp)=kj
                  endif
                enddo
              endif
            endif
          endif
        endif
        call tflocald(ke)
      endif
      if(mode .eq. 0)then
        kx%k=ktfoper+mtfnull
        irtc=0
      elseif(mode .eq. 1)then
        kx=kxlistcopied(ispb)
        irtc=0
      else
        do i=ispb+1,isp
          call tflocal(ktastk(i))
        enddo
        if(mode .eq. 2)then
          ktastk(ispb)=ktfoper+mtfplus
        else
          ktastk(ispb)=ktfoper+mtfmult
        endif
        kx=tfefunref(ispb,.true.,irtc)
      endif
 9000 isp=ispb-1
      if(var)then
        call tfdelete(symd,.true.,.false.)
      endif
      if(irtc .eq. 0)then
        return
      elseif(irtc .eq. -3)then
        kx%k=ktfoper+mtfnull
        irtc=0
c      else
c        call tfcatchreturn(0,kx,irtc)
      endif
      return
 9100 lv=itfdownlevel()
      do i=ispb+1,isp
        call tflocal(ktastk(i))
      enddo
 9200 call tflocald(ke)
      go to 9000
 9500 irtc=itfmessage(9,'General::wrongtype',
     $     '"{symbol, n}, {symbol, start, stop},'//
     $     ' {symbol, start, stop, step}"')
      return
      end

      type (sad_descriptor) function kxlistcopied(isp1)
      use tfstk
      implicit none
      type (sad_dlist), pointer ::klx
      type (sad_rlist), pointer ::klr
      integer*4 isp1,i,n
      logical*4 nr,re
      if(isp1 .ge. isp)then
        kxlistcopied=dxnulll
      else
        n=isp-isp1
        nr=.false.
        re=.false.
        do i=isp1+1,isp
          if(ktfrealq(dtastk(i)))then
            re=.true.
            if(nr)then
              go to 10
            endif
          else
            nr=.true.
            if(re)then
              go to 10
            endif
          endif
        enddo
 10     if(nr)then
          kxlistcopied=kxadaloc(-1,n,klx)
        else
          kxlistcopied=kxavaloc(-1,n,klr)
          call descr_sad(kxlistcopied,klx)
        endif
        klx%dbody(1:n)=dtastk(isp1+1:isp1+n)
      endif
      return
      end

      subroutine tfrange(isp1,kx,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) kx
      type (sad_rlist), pointer :: kl
      integer*8 n
      integer*4 isp1,irtc,narg,i,itfmessage
      real*8 x0,x1,xs,xi
      narg=isp-isp1
      if(narg .gt. 3)then
        irtc=itfmessage(9,'General::narg','"1, 2, or 3"')
        return
      endif
      do i=isp1+1,isp
        if(ktfnonrealq(ktastk(i)))then
          irtc=itfmessage(9,'General::wrongtype','"Real number"')
          return
        endif
      enddo
      if(narg .eq. 1)then
        x0=1.d0
        x1=rtastk(isp)
        xs=1.d0
      elseif(narg .eq. 2)then
        x0=rtastk(isp-1)
        x1=rtastk(isp)
        xs=1.d0
      else
        x0=rtastk(isp-2)
        x1=rtastk(isp-1)
        xs=rtastk(isp)
      endif
      if(xs .eq. 0.d0)then
        irtc=itfmessage(9,'General::wrongnum','"non-zero"')
        return
      endif
      n=max(int8((x1-x0)/xs*(1.d0+1.d-11)+1.d0),0)
      if(n .gt. 2**31-1)then
        irtc=itfmessage(9,'General::wrongval',
     $       '"# of elements","less than 2^31"')
        return
      endif
      kx=kxavaloc(-1,int(n),kl)
      kl%attr=ior(kl%attr,lconstlist)
      xi=0.d0
      do i=1,int(n)
        kl%rbody(i)=x0+xs*xi
        xi=xi+1.d0
      enddo
      irtc=0
      return
      end
