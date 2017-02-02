      recursive subroutine tfeexpr(k1,k,ke,iopc1)
      use tfstk
      implicit none
      type (sad_descriptor) k1,k,ke,kx,ky,ky1,k2,kx2
      type (sad_list), pointer :: listy,list1,listi,klx,kl1,kle
      type (sad_symbol), pointer :: sym
      type (sad_string), pointer :: str
      type (sad_pat), pointer :: kp1
      type (sad_complex), pointer :: cx
      integer*8 ks,isp0
      integer*4 m,irtc,i,iopc1,iopc,itfcanonicalorder
      real*8 vx1,vy,v2,vx
      logical*4 eval
      iopc=iopc1
      ky=k
      select case (iopc)
      case(mtfplus,mtftimes)
        if(.not. tfnumberqd(k1) .and. tfnumberqd(ky))then
          kx=ky
          ky=k1
        else
          kx=k1
        endif
        if(ktflistqd(ky,listy))then
          if(listy%head .eq. ktfoper+iopc)then
            call tfjoine(kx,ky,ke,irtc)
            if(irtc .ne. -1)then
              return
            endif
            irtc=0
            m=listy%nl
            if(tfnumberqd(listy%dbody(1)))then
              if(tfnumberqd(kx))then
                call tfcmplx(listy%dbody(1),kx,ky1,iopc,irtc)
                if(irtc .ne. 0)then
                  ke%k=ktfoper+mtfnull
                  return
                endif
                if(ktfrealqd(ky1))then
                  if(iopc .eq. mtfplus)then
                    if(ky1%k .eq. 0)then
                      if(m .eq. 2)then
                        ke=listy%dbody(2)
                      else
                        call tftake(ky,dble(-m+1),ke,
     $                       .true.,.false.,irtc)
                      endif
                      return
                    endif
                  else
                    if(ky1%k .eq. ktftrue)then
                      if(m .eq. 2)then
                        ke=listy%dbody(2)
                      else
                        call tftake(ky,dble(-m+1),ke,
     $                       .true.,.false.,irtc)
                      endif
                      return
                    elseif(ky1%k .eq. 0 .and.
     $                     redmath%value%k .ne. 0)then
                      ke%k=0
                      return
                    endif
                  endif
                endif
                call tfclonelist(listy,listy)
                call tfreplist(listy,1,ky1,eval)
                ke%k=ktflist+ksad_loc(listy%head)
              else
                call tfinsertsort(listy,kx,ke)
              endif
              return
            endif
            if(ktfrealqd(kx))then
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
            call tfinsertsort(listy,kx,ke)
            return
          endif
        endif
        if(ktflistqd(kx,klx))then
          if(klx%head .eq. ktfoper+iopc)then
            call tfinsertsort(klx,ky,ke)
            return
          endif
        elseif(ktfrealqd(kx))then
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
        if(ktflistqd(ky,listy))then
          if(listy%head .eq. ktfoper+mtfnot)then
            ke=listy%dbody(1)
            return
          endif
        endif
        go to 5000
      case (mtfslot,mtfslotseq)
        if(ktfrealqd(ky,vy))then
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
        if(iopc .eq. mtfslot)then
          ke=dlist(iaxslotnull+(ks-1)*2)
        else
          ke=dlist(iaxslotnull+(ks-1)*2+1)
        endif
c        call tfdebugprint(ke,'tfeexpr-slot-end',3)
        return
      case (mtfflag)
        go to 5000
      case (mtfcomp,mtfconcat,mtfand,mtfor,mtfalt,mtfmessagename)
        if(ktflistqd(k1,list1))then
          if(list1%head .eq. ktfoper+iopc)then
            if(ktflistqd(ky,listy))then
              if(listy%head .eq. ktfoper+iopc)then
                call tfjoin2(k1,ky,ke,.false.,irtc)
                return
              endif
            endif
            call tfappend(k1,ky,ke,.false.,0,irtc)
            return
          endif
        endif
        if(ktflistqd(ky,listy))then
          if(listy%head .eq. ktfoper+iopc)then
            call tfappend(ky,k1,ke,.false.,1,irtc)
            return
          endif
        endif
      case (mtfset,mtfpower)
        if(ktflistqd(ky,listy))then
          if(listy%head .eq. ktfoper+iopc)then
            call tfappend(ky,k1,ke,.false.,1,irtc)
            return
          endif
        endif
        if(iopc .eq. mtfpower)then
          if(ktfrealqd(ky))then
            if(ky%k .eq. ktftrue)then
              ke=k1
              return
            elseif(ky%k .eq. 0 .and. redmath%value%k .ne. 0)then
              ke%k=ktftrue
              return
            endif
          endif
          if(ktflistqd(k1,list1))then
            if(list1%head .eq. ktfoper+mtfpower)then
              m=list1%nl
              if(m .eq. 1)then
                kx=ky
              else
                if(m .eq. 2)then
                  k2=list1%dbody(2)
                else
                  call tftake(k1,dble(-m+1),k2,.true.,.false.,irtc)
                endif
                if(ktfrealqd(k2,v2) .and. ktfrealqd(ky,vy))then
                  kx=dfromr(v2*vy)
                else
                  call tfeexpr(k2,ky,kx,mtftimes)
                endif
              endif
              ky=list1%dbody(1)
              call tfeexpr(ky,kx,ke,mtfpower)
              return
            elseif(list1%head .eq. ktfoper+mtftimes
     $             .and. ktfrealqd(ky))then
              vx1=1.d0
              kx2%k=ktftrue
              m=list1%nl
              isp0=isp
              isp=isp+1
              ktastk(isp)=ktfoper+mtftimes
              do i=1,m
                isp=isp+1
                ktastk(isp)=list1%body(i)
                if(ktfrealq(ktastk(isp)))then
                  vx1=vx1*rtastk(isp)
                  isp=isp-1
                elseif(ktflistq(ktastk(isp),listi))then
                  if(listi%head .eq. ktfoper+mtfcomplex)then
                    call tfeexpr(dtastk(isp),ky,kx,mtfpower)
                    if(ktfrealqd(kx,vx))then
                      vx1=vx1*vx
                    else
                      call tfeexpr(kx2,kx,kx2,mtftimes)
                    endif
                    isp=isp-1
                  endif
                endif
              enddo
              if(isp .ne. isp0+m+1)then
                if(isp .gt. isp0+2)then
                  call tfcompose(isp0+1,ktfoper+mtftimes,kx,irtc)
                elseif(isp .eq. isp0+2)then
                  kx=dtastk(isp)
                else
                  kx%k=ktftrue
                endif
                isp=isp0
                call tfeexpr(kx,ky,ke,mtfpower)
                if(ktfrealqd(kx2,v2))then
                  vx1=vx1*v2
                else
                  call tfeexpr(kx2,ke,ke,mtftimes)
                endif
                if(vx1 .ne. 1.d0)then
                  call tfcmplx(vx1,ky,kx,mtfpower,irtc)
                  call tfeexpr(kx,ke,ke,mtftimes)
                endif
                return
              else
                isp=isp0
              endif
            endif
          elseif(ktfrealqd(k1))then
            if(k1%k .eq. ktftrue .and. redmath%value%k .ne. 0)then
              ke%k=ktftrue
              return
            endif
          endif
        endif
      case (mtfreplace,mtfreplacerepeated)
        if(ktflistqd(k1,kl1))then
          if(kl1%head .eq. ktfoper+iopc)then
            call tfappend(k1,ky,ke,.false.,0,irtc)
            return
          endif
        endif
      case (mtfcolon)
        if(ktfsymbolqd(k1,sym))then
          call sym_symstr(sym,str)
          ke=kxpalocb(str%str,str%nch,ky,transfer(ktfref,k))
          return
        elseif(ktfpatqd(k1,kp1))then
          kp1%default=dtfcopy(ky)
          ke%k=ktfpat+ktfaddrd(k1)
          return
        endif
      case (mtfrevpower)
        call tfeexpr(k,k1,ke,mtfpower)
        return
      case (mtfatt)
c        call tfdebugprint(k1,'eexpr',1)
c        call tfdebugprint(ky,'@',1)
        if(ktfnonrealqd(k1))then
          if(ktfrealqd(k))then
            ke=kxavaloc(-1,1,kle)
            kle%dbody(1)=k
            kle%dbody(0)=dtfcopy(k1)
            return
          elseif((ktfsymbolqd(ky) .or. ktfoperqd(ky)) .and.
     $           (ktfsymbolqd(k1) .or. ktflistqd(k1)) .or.
     $           ktfpatqd(ky))then
            go to 4900
          elseif(ktflistqd(k1,kl1) .and.
     $           kl1%head .eq. ktfoper+mtfatt)then
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
        if(ktfrealqd(k))then
          if(k%k .eq. 0)then
            ke=k1
            return
          endif
        elseif(tfcomplexqx(k%k,cx))then
          call tfeexpr(cx%dbody(1),k1,ky,mtfplus)
          call tfeexpr(ky,cx%dbody(2),ke,mtfcomplex)
          return
        elseif(tfcomplexqx(k1%k,cx))then
          call tfeexpr(cx%dbody(2),k,ky,mtfplus)
          call tfeexpr(cx%dbody(1),ky,ke,mtfcomplex)
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
      ke=kxcompose(isp-2)
      isp=isp-3
      return
 5000 if(ktfrealqd(ky))then
        ke=kxavaloc(-1,1,kle)
        kle%dbody(1)=ky
      else
        ke=kxadaloc(-1,1,kle)
        kle%dbody(1)=dtfcopy(ky)
      endif
      kle%head=ktfoper+iopc
      return
      end

      subroutine tfinsertsort(kl,ki,kx)
      use tfstk
      implicit none
      type (sad_descriptor) ki,kx
      type (sad_list) kl
      integer*8 isp0,isp1,isp2,ispm,isp3
      integer*4 i,itfcanonicalorder
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
c      do i=isp0+1,isp1-1
c        isp=isp+1
c        ktastk(isp)=ktastk(i)
c      enddo
      isp=isp3+1
      dtastk(isp)=ki
      ktastk(isp+1:isp+isp2-isp1+1)=ktastk(isp1:isp2)
      isp=isp+isp2-isp1+1
c      do i=isp1,isp2
c        isp=isp+1
c        ktastk(isp)=ktastk(i)
c      enddo
      kx=kxcrelistm(int(isp-isp2),ktastk(isp2+1:isp),kl%dbody(0))
      isp=isp0
      return
      end

      subroutine tfcmplx(k1,k2,kx,iopc,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) k1,k2,kx
      type (sad_complex), pointer :: cx1,cx2
      integer*8 ki1,ki2
      integer*4 irtc,iopc
      real*8 v1,v2,tfenum
      complex*16 c1
      irtc=0
      if(iopc .eq. mtfnot)then
        if(ktfrealqd(k2))then
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
      if(ktfrealqd(k1,v1))then
        if(ktfrealqd(k2,v2))then
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
                  kx=dfromr(dble(c1))
                else
                  kx=kxcalocv(-1,dble(c1),imag(c1))
                endif
              else
                kx=dfromr(v1**ki2)
              endif
              return
            endif
            kx=dfromr(tfenum(k1,k2,mtfpower,irtc))
          elseif(iopc .eq. mtfrevpower)then
            if(k2%k .lt. 0)then
              ki1=int(v1)
              if(ki1 .ne. v1)then 
                c1=dcmplx(v2,0.d0)**v1
                if(imag(c1) .eq. 0.d0)then
                  kx=dfromr(dble(c1))
                else
                  kx=kxcalocv(-1,dble(c1),imag(c1))
                endif
              else
                kx=dfromr(v2**ki1)
              endif
              return
            endif
            kx=dfromr(tfenum(v2,v1,mtfpower,irtc))
          else
            kx=dfromr(tfenum(k1,k2,iopc,irtc))
          endif
          return
         elseif(tfcomplexqx(k2%k,cx2))then
          call tfcmplxmath(dcmplx(v1,0.d0),cx2%cx,kx,iopc,irtc)
          return
        else
          go to 8000
        endif
      elseif(tfcomplexqx(k1%k,cx1))then
        if(ktfrealqd(k2,v2))then
          call tfcmplxmath(cx1%cx,dcmplx(v2,0.d0),
     $         kx,iopc,irtc)
          return
        elseif(tfcomplexqx(k2%k,cx2))then
          call tfcmplxmath(cx1%cx,cx2%cx,kx,iopc,irtc)
          return
        endif
      endif
 8000 call tfeval1(k1,k2,kx,iopc,irtc)
      return
      end

      real*8 function tfenum(v1,v,iopc1,irtc)
      use tfstk
      implicit none
      real*8 v1,v
      integer*4 ix,iopc1,irtc,itfmessage
      irtc=0
      go to (
     $     2000,
     $     1010,1020,1030,2000,1050,2000,1070,1080,1090,1100,
     $     1110,1120,1130,1140,1090,1100,1170,1180,1190),iopc1+1
c          m    i    +    -    *    /    v    ^    e    n    
c          >    <    g    l    E    N    ~    &    o    c
      go to 2000
 1010 tfenum=-v
      return
 1020 tfenum=1.d0/v
      return
 1030 tfenum=v1+v
      return
 1050 if(abs(v) .eq. dinfinity .or. abs(v1) .eq. dinfinity)then
        if(v .eq. 0.d0 .or. v1 .eq. 0.d0)then
          tfenum=v*v1
        elseif(v .gt. 0.d0 .and. v1 .gt. 0.d0 .or.
     $         v .lt. 0.d0 .and. v1 .lt. 0.d0)then
          tfenum=dinfinity
        else
          tfenum=-dinfinity
        endif
      else
        tfenum=v1*v
      endif
      return
 1080 if(v .eq. -1.d0)then
        if(abs(v1) .eq. dinfinity)then
          tfenum=0.d0
        else
          tfenum=1.d0/v1
        endif
      elseif(v .eq. 2.d0)then
        tfenum=v1**2
      elseif(v .eq. .5d0)then
        tfenum=sqrt(v1)
      elseif(v .eq. 0.d0 .and. redmath%value%k .ne. 0)then
        tfenum=1.d0
      else
        ix=int(v)
        if(ix .eq. v)then
          tfenum=v1**ix
        else
          tfenum=v1**v
        endif
      endif
      return
 1070 if(v1 .eq. -1.d0)then
        tfenum=1.d0/v
      elseif(v1 .eq. 2.d0)then
        tfenum=v**2
      elseif(v1 .eq. .5d0)then
        tfenum=sqrt(v)
      elseif(v1 .eq. 0.d0 .and. redmath%value%k .ne. 0)then
        tfenum=1.d0
      else
        ix=int(v1)
        if(ix .eq. v1)then
          tfenum=v**ix
        else
          tfenum=v**v1
        endif
      endif
      return
 1090 if(v1 .eq. v)then
        tfenum=1.d0
      else
        tfenum=0.d0
      endif
      return
 1100 if(v1 .ne. v)then
        tfenum=1.d0
      else
        tfenum=0.d0
      endif
      return
 1110 if(v1 .gt. v)then
        tfenum=1.d0
      else
        tfenum=0.d0
      endif
      return
 1120 if(v1 .lt. v)then
        tfenum=1.d0
      else
        tfenum=0.d0
      endif
      return
 1130 if(v1 .ge. v)then
        tfenum=1.d0
      else
        tfenum=0.d0
      endif
      return
 1140 if(v1 .le. v)then
        tfenum=1.d0
      else
        tfenum=0.d0
      endif
      return
 1170 if(v .eq. 0.d0)then
        tfenum=1.d0
      else
        tfenum=0.d0
      endif
      return
 1180 if(v1 .ne. 0.d0 .and. v .ne. 0.d0)then
        tfenum=1.d0
      else
        tfenum=0.d0
      endif
      return
 1190 if(v1 .ne. 0.d0 .or. v .ne. 0.d0)then
        tfenum=1.d0
      else
        tfenum=0.d0
      endif
      return
 2000 irtc=itfmessage(999,'General::invop',' ')
      tfenum=0.d0
      return
      end

      subroutine tfjoine(k1,k2,kx,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) k1,k2,kx,k10,k20,ky1
      type (sad_list), pointer ::kl1,kl2
      integer*8 isp1
      integer*4 irtc,ma1,ma2,m,iopc
      logical*4 tfsameheadqk
      if(ktfnonlistqd(k1,kl1) .or. ktfnonlistqd(k2,kl2))then
        irtc=-1
        return
      endif
      if( .not. tfsameheadqk(k1,k2))then
        irtc=-1
        return
      endif
      irtc=0
      iopc=int(ktfaddr(kl1%head))
      if(iopc .eq. mtfplus .or. iopc .eq. mtftimes)then
        if(.not. tfnumberqd(kl1%dbody(1)))then
          call tfjoin2(k2,k1,kx,.false.,irtc)
          go to 1000
        endif
        if(.not. tfnumberqd(kl2%dbody(1)))then
          call tfjoin2(k1,k2,kx,.false.,irtc)
          go to 1000
        endif
        ma1=kl1%nl
        ma2=kl2%nl
        m=ma1+ma2-1
        k10=kl1%dbody(1)
        k20=kl2%dbody(1)
        call tfcmplx(k10,k20,ky1,iopc,irtc)
        if(irtc .ne. 0)then
          return
        endif
        if(ktfrealqd(ky1))then
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
        kx=kxcrelistm(int(isp-isp1),ktastk(isp1+1:isp),
     $       k_descr(ktfoper+iopc))
        isp=isp1
      else
        call tfjoin2(k1,k2,kx,.false.,irtc)
        return
      endif
 1000 isp=isp+1
      dtastk(isp)=kx
      call tfsort(isp-1,kx,0,irtc)
      isp=isp-1
      return
      end

      subroutine tfjoin2(k1,k2,kx,eval,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) k1,k2,kx
      integer*4 irtc
      integer*8 isp1
      logical*4 eval
      isp1=isp
      isp=isp+1
      dtastk(isp)=k1
      isp=isp+1
      dtastk(isp)=k2
      call tfjoin(isp1,kx,eval,irtc)
      isp=isp1
      return
      end

      subroutine tfjoin(isp1,kx,eval,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) kx,kf
      type (sad_list), pointer :: kl1,kli
      integer*8 isp1,isp0,i
      integer*4 irtc,itfmessage,narg
      logical*4 tfsameqk,eval,ev
      narg=int(isp-isp1)
      if(ktfnonlistq(ktastk(isp1+1),kl1))then
        go to 9010
      endif
      if(narg .le. 1)then
        if(narg .eq. 1)then
          kx=dtastk(isp1+1)
          irtc=0
        else
          irtc=itfmessage(9,'General::narg','"1 or more"')
        endif
        return
      endif
      kf=kl1%dbody(0)
      isp=isp+1
      isp0=isp
      call tfgetllstkall(kl1)
      if(isp .ge. mstk)then
        irtc=itfmessage(9999,'General::stack','"Join"')
        return
      endif
      do i=isp1+2,isp0-1
        if(ktfnonlistq(ktastk(i),kli))then
          isp=isp0-1
          go to 9000
        endif
        if(.not. tfsameqk(kli%head,kf))then
          go to 9100
        endif
        call tfgetllstkall(kli)
        if(isp .ge. mstk)then
          irtc=itfmessage(9999,'General::stack','"Join"')
          return
        endif
      enddo
      dtastk(isp0)=kf
      ev=eval
      if(ev)then
        if(kf%k .eq. ktfoper+mtflist .or.
     $       kf%k .eq. ktfoper+mtfalt .or.
     $       kf%k .eq. ktfoper+mtfnull)then
          ev=.false.
        endif
      endif
      if(ev)then
        call tfefunref(isp0,kx,.true.,irtc)
      else
        kx=kxcompose(isp0)
        irtc=0
      endif
      isp=isp0-1
      return
 9000 isp=isp0-1
 9010 irtc=itfmessage(9,'General::wrongtype',
     $     '"List or composition for all args"')
      return
 9100 irtc=itfmessage(9,'General::samehead',' ')
      isp=isp0-1
      return
      end

      subroutine tfappend(kl,k,kx,eval,mode,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) kl,k,kx
      type (sad_list), pointer :: list,listx
      integer*4 irtc,m,itfmessage,mode,i
      logical*4 eval,ev,tfconstqk
      if(.not. ktflistqd(kl,list))then
        irtc=itfmessage(9,'General::wrongtype',
     $       '"List or composition for #1"')
        return
      endif
      ev=eval .and. list%head .ne. ktfoper+mtflist .and.
     $     list%head .ne. ktfoper+mtfalt .and.
     $     list%head .ne. ktfoper+mtfnull
      m=list%nl
      call loc_list(ktaaloc(-1,m+1),listx)
      listx%attr=list%attr
      if(ktfreallistqo(list))then
        if(mode .eq. 0)then
          listx%body(1:m)=list%body(1:m)
          if(ktfrealqd(k))then
            listx%dbody(m+1)=k
          else
            listx%dbody(m+1)=dtfcopy(k)
            listx%attr=ior(listx%attr,lnonreallist)
          endif
        else
          listx%body(2:m+1)=list%body(1:m)
          if(ktfrealqd(k))then
            listx%dbody(1)=k
          else
            listx%dbody(1)=dtfcopy(k)
            listx%attr=ior(listx%attr,lnonreallist)
          endif
        endif
      else
        if(mode .eq. 0)then
          do i=1,m
            listx%body(i)=ktfcopy(list%body(i))
          enddo
          listx%dbody(m+1)=dtfcopy(k)
        else
          do i=1,m
            listx%body(i+1)=ktfcopy(list%body(i))
          enddo
          listx%dbody(1)=dtfcopy(k)
        endif
      endif
      listx%head=ktfcopy(list%head)
      if(iand(list%attr,kconstarg) .ne. 0)then
        if(.not. tfconstqk(k%k))then
          listx%attr=ior(listx%attr-kconstarg,knoconstarg+lnoconstlist)
        endif
      endif
      if(ev)then
        call tfleval(listx,kx,.true.,irtc)
      else
        kx%k=ktflist+sad_loc(listx%head)
        irtc=0
      endif
      return
      end

      recursive subroutine tfcmplxf(k,kx,mode,iaf)
      use tfstk
      implicit none
      type (sad_descriptor) k,kx
      type (sad_list), pointer :: kl,klx
      integer*8 isp0
      integer*4 m,i,mode,iaf
      if(tfnumberqd(k))then
        if(ktfrealqd(k))then
          if(mode .eq. 1)then
            kx=k
          elseif(mode .eq. 2)then
            kx%k=0
          else
            kx=k
          endif
        else
          call loc_list(ktfaddrd(k),kl)
          if(mode .eq. 1)then
            kx=kl%dbody(1)
          elseif(mode .eq. 2)then
            kx=kl%dbody(2)
          else
            kx=kxcalocv(-1,kl%rbody(1),-kl%rbody(2))
          endif
        endif
        return
      elseif(tfreallistqd(k,kl))then
        if(mode .eq. 1 .or. mode .eq. 3)then
          kx=k
        else
          kx%k=ktflist+ktraaloc(-1,kl%nl)
        endif
        return
      elseif(ktflistqd(k,kl))then
        m=kl%nl
        if(mode .eq. 3 .or. kl%head .eq. ktfoper+mtflist)then
          isp0=isp
          do i=1,m
            isp=isp+1
            call tfcmplxf(kl%dbody(i),dtastk(isp),mode,iaf)
          enddo
          kx=kxmakelist(isp0)
          isp=isp0
          return
        endif
      endif
      kx=kxadaloc(-1,1,klx)
      klx%head=ktfoper+iaf
      klx%dbody(1)=dtfcopy(k)
      return
      end
