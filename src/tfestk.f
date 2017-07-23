
      module ophash
      use tfstk
      implicit none
      integer*4 nhash
      parameter (nhash=4)
      integer*4 iophash(nhash,0:63)
      logical*4 , save :: opini=.true.
      integer*8 ktfcode(0:ntfarg)
      character*4, save :: opcode(0:mtfnopc) =(/
     $     '    ','    ','    ','+   ','-   ',
     $     '*   ','/   ','    ','^   ','==  ',

     $     '<>  ','>   ','<   ','>=  ','<=  ',
     $     '=== ','<=> ','~   ','&&  ','||  ',

     $     '//  ','[   ',']   ','{   ','}   ',
     $     ':=  ','=   ','    ','(   ',')   ',

     $     ',   ',';   ','&   ',':   ','->  ',
     $     ':>  ','/.  ','//. ','^=  ','^:= ',

     $     '=.  ','?   ','?   ','#   ','##  ',
     $     '.   ','|   ','/@  ','//@ ','@@  ',

     $     '..  ','... ','    ','+=  ','-=  ',
     $     '*=  ','/=  ','++  ','--  ','[[  ',

     $     '@   ','::  ','/:  ','(*  ','*)  ',
     $     '    ','    '/)
      logical*4 :: constop(0:mtfnopc) = (/
     $     .false.,.false.,.false.,.false.,.false.,
     $     .false.,.false.,.false.,.false.,.false.,

     $     .false.,.false.,.false.,.false.,.false.,
     $     .false.,.false.,.false.,.false.,.false.,

     $     .false.,.false.,.false.,.true. ,.false.,
     $     .false.,.false.,.true., .false.,.false.,

     $     .false.,.false.,.true. ,.true., .true.,
     $     .true., .false.,.false.,.false.,.false.,

     $     .false.,.true. ,.false.,.false.,.false.,
     $     .false.,.false. ,.false.,.false.,.false.,
c Alternatives is temporarily set .false. due to possible reducution.

     $     .true. ,.true. ,.false.,.false.,.false.,
     $     .false.,.false.,.false.,.false.,.false.,

     $     .false.,.true., .true. ,.false.,.false.,
     $     .true., .false.
     $     /)

      contains
        subroutine tfopcodehash
        integer*4 i,j,k
        character*4 oper1
        do i=0,63
          do k=1,nhash
            iophash(k,i)=-1
          enddo
        enddo
        LOOP_I: do i=0,mtfnopc
          if(opcode(i) .ne. ' ')then
            oper1=opcode(i)
            j=ichar(oper1(1:1))+ichar(oper1(2:2))
     $           +ichar(oper1(3:3))+ichar(oper1(4:4))
            j=iand(j,63)
            do k=1,nhash
              if(iophash(k,j) .lt. 0)then
                iophash(k,j)=i
                cycle LOOP_I
              endif
            enddo
            write(*,*)'tfopcode hash implementation error. ',opcode(i)
            stop
          endif
        enddo LOOP_I
        opini=.false.
        return
        end subroutine

      end module

      module opdata
      use tfstk
      implicit none
      integer*4 :: iprior(0:mtfnopc) = (/
     $     9999,
     $     10,  20,  50,  50,  40,  40,  15,  15,  100, 100,
     $     100, 100, 100, 100, 120, 120, 150, 160, 170, 80,
     $     6,   3000,9999,3000,250, 250, 7000,9999,8000,9000,
     $     1000,220, 180, 190, 190, 200, 200, 250, 250, 900,
     $     4,   3,   3,   3,   10,  175, 9,   9,   9,   172,
     $     172, 130, 210, 210, 210, 210, 7,   7,   6,   5,
     $     2,   240, 9999,9999,1,   9999/)
c          null
c          m    i    +    -    *    /    v    ^    ==   <>
c          >    <    >=   <=   ===  <=>  ~    &&   ||   //
c          [    ]    {    }    :=   =    C    (    )    ,
c          ;    &    :    ->   :>   /.   //.  ^=   ^:=  =.
c          ?    flg  #    ##   .    |    /@   //@  @@   .. 
c          ...  ineq +=   -=   *=   /=   ++   --   [[   @
c          msgn /:   (*   *)   Hold z
      logical*4, parameter :: T=.true.,F=.false.
      logical*4 :: nullfirst(0:mtfnopc) = (/
     $     T,
     $     T,   T,   F,   F,   F,   F,   F,   F,   F,   F,   
     $     F,   F,   F,   F,   F,   F,   T,   F,   F,   F,   
     $     T,   T,   T,   T,   F,   F,   F,   T,   T,   T,
     $     T,   F,   F,   F,   F,   F,   F,   F,   F,   F,   
     $     T,   T,   T,   T,   F,   F,   F,   F,   F,   F,   
     $     F,   F,   F,   F,   F,   F,   T,   T,   F,   F,
     $     F,   F,   F,   F,   F,   F/)
      logical*4 :: lastfirst(0:mtfnopc) = (/
     $     F,
     $     F,   F,   F,   F,   F,   F,   F,   T,   F,   F,   
     $     F,   F,   F,   F,   F,   F,   F,   F,   F,   F,   
     $     F,   F,   F,   F,   T,   T,   F,   F,   F,   F,
     $     F,   F,   F,   T,   T,   F,   F,   F,   F,   F,   
     $     F,   F,   F,   F,   F,   F,   T,   T,   T,   F,   
     $     F,   F,   F,   F,   F,   F,   F,   F,   F,   F,
     $     F,   F,   F,   F,   F,   F/)

      end module

      module tfform
      use tfstk
      integer*8, save :: iaxform=0,iaxpagewidth=0
      type (sad_symbol), pointer, save :: symform,sympw
      contains
        subroutine tfforminit
        implicit none
        iaxform=ktfsymbolz('System`$FORM',12)
        iaxpagewidth=ktfsymbolz('System`PageWidth',16)
        call loc_sym(iaxform,symform)
        call loc_sym(iaxpagewidth,sympw)
        end subroutine
      end module

      subroutine tfestk(isp0,iprior,lastfirst,irtc)
      use tfstk
      use ophash
      implicit none
      type (sad_descriptor) kx,kh
      type (sad_list), pointer :: klx
      integer*4 isp0,iprior(0:mtfnopc),irtc
      integer*4 iop,iop1,isp1,i,itgetfpe,itfmessage
      logical*4 tfconstqk,tfheldqd,lastfirst(0:mtfnopc)
      irtc=0
      do while(isp .gt. isp0)
        iop=itastk2(1,isp)
        if(iop .eq. mtfnull .or.
     $       iop .eq. mtfleftparen .or. iop .eq. mtflist)then
          return
        endif
        isp1=isp-1
        iop1=itastk2(1,isp1)
        if(iop .eq. mtfrightbra)then
          if(iop1 .eq. mtfrightbra)then
            if(ktastk(isp) .eq. ktfoper+mtfnull)then
              isp=isp-1
              do i=isp,isp0,-1
                if(itastk2(1,i) .eq. mtfpart)then
                  isp1=i
                  kx=kxmakelist(isp1-1,klx)
                  klx%head%k=ktfoper+mtfpart
                  iop=mtfnull
                  go to 1010
                endif
              enddo
            endif
            go to 9010
          elseif(iop1 .eq. mtfcomma .or.
     $           iop1 .eq. mtfleftbra)then
            do i=isp1,isp0,-1
              if(itastk2(1,i) .eq. mtfleftparen
     $             .or. itastk2(1,i) .eq. mtflist)then
                go to 9010
              elseif(itastk2(1,i) .eq. mtfpart)then
                return
              elseif(itastk2(1,i) .eq. mtfleftbra)then
                isp1=i
                if(isp .eq. isp1+1 .and.
     $               ktastk(isp) .eq. ktfoper+mtfnull)then
                  isp=isp1
                endif
                kh=dtastk(i)
                call tfcomposefull(i,kh,kx,irtc)
                if(irtc .ne. 0)then
                  return
                endif
                iop=mtfnull
                go to 1010
              endif
            enddo
            go to 9010
          elseif(iop1 .eq. mtfpart)then
            return
          endif
        elseif(iop .eq. mtfrightbrace .and.
     $         (iop1 .eq. mtfcomma .or. iop1 .eq. mtflist))then
          do i=isp1,isp0,-1
            if(itastk2(1,i) .eq. mtfleftparen
     $           .or. itastk2(1,i) .eq. mtfleftbra
     $           .or. itastk2(1,i) .eq. mtfpart)then
              go to 9010
            elseif(itastk2(1,i) .eq. mtflist)then
              isp1=i
              kx=kxmakelist(isp1)
              iop=mtfnull
              go to 1010
            endif
          enddo
          go to 9010
        elseif(iop .eq. mtfrightparen
     $         .and. iop1 .eq. mtfleftparen)then
          ktastk(isp1)=ktastk(isp)
          itastk2(1,isp1)=mtfnull
          isp=isp1
          return
        elseif(iop1 .eq. mtffun .and.
     $         ktastk(isp) .eq. ktfoper+mtfnull)then
          kx=kxpfaloc(dtastk(isp-1))
          go to 1010
        endif
        if(iprior(iop) .lt. iprior(iop1))then
          return
        elseif(iprior(iop) .eq. iprior(iop1) .and.
     $         lastfirst(iop))then
          return
        endif
        if(iop1 .eq. mtfcomma)then
          if(iop .ne. mtfcomma)then
            go to 9020
          endif
          return
        elseif(iop1 .eq. mtfpower)then
          if(ktastk(isp) .eq. ktfoper+mtfnull)then
            return
          endif
        elseif(iop1 .eq. mtfleftbra .or. iop1 .eq. mtfpart)then
          return
        endif
        if(iop1 .ge. mtfequal .and. iop1 .le. mtfunsame)then
          call tfeinequal(dtastk(isp1),dtastk(isp),kx,iop1,
     $         iop .ge. mtfequal .and. iop .le. mtfunsame)
        else
          if(iop1 .eq. mtfneg)then
            iop1=mtftimes
          elseif(iop1 .eq. mtfnot)then
            if(ktastk(isp1) .ne. ktfoper+mtfnull)then
              go to 9030
            endif
          endif
          if(tfconstqk(ktastk(isp1)) .and.
     $         .not. tfheldqd(dtastk(isp1)) .and.
     $         tfconstqk(ktastk(isp)) .and.
     $         .not. tfheldqd(dtastk(isp)) .and.
     $         ktfimmediateq(klist(ifunbase+iop1)))then
            call tfeval1(ktastk(isp1),ktastk(isp),kx%k,iop1,irtc)
            if(irtc .ne. 0)then
              return
            endif
            if(itgetfpe() .gt. 0)then
              irtc=itfmessage(9,'General::fpe','""')
              return
            endif
          else
            call tfeexpr(dtastk(isp1),dtastk(isp),kx,iop1)
          endif
        endif
 1010   isp=isp1
        itastk2(1,isp)=iop
        dtastk(isp)=kx
      enddo
      return
 9010 irtc=itfmessage(9999,'General::mismatch',
     $     '"'//opcode(iop)//'"')
      return
 9020 irtc=itfmessage(9999,'General::incomplete','""')
      return
 9030 irtc=itfmessage(9,'General::narg','"1 for ~"')
      return
      end

      subroutine tfeinequal(k1,k2,kx,iop1,nextrel)
      use tfstk
      implicit none
      type (sad_descriptor) k1,k2,kx
      type (sad_list), pointer :: kl1,kl2
      type (sad_dlist), pointer :: klx
      integer*4 iop1,i,m1,m2,irtc
      logical*4 tfinequalityqk,nextrel,tfconstlistqo
      if(tfinequalityqk(k1%k))then
        call loc_sad(ktfaddrd(k1),kl1)
        m1=kl1%nl
        if(tfinequalityqk(k2))then
          call loc_sad(ktfaddrd(k2),kl2)
          m2=kl2%nl
          kx=kxadaloc(-1,m1+m2+1,klx)
          klx%head%k=ktfoper+mtfinequality
          do i=1,m1
            klx%body(i)=ktfcopy(kl1%body(i))
          enddo
          klx%body(m1+1)=ktfoper+iop1
          do i=1,m2
            klx%body(i+m1+1)=ktfcopy(kl2%body(i))
          enddo
        else
          kx=kxadaloc(-1,m1+2,klx)
          klx%head%k=ktfoper+mtfinequality
          do i=1,m1
            klx%dbody(i)=dtfcopy(kl1%dbody(i))
          enddo
          klx%body(m1+1)=ktfoper+iop1
          klx%dbody(m1+2)=dtfcopy(k2)
        endif
      else
        if(tfinequalityqk(k2%k))then
          call loc_sad(ktfaddrd(k2),kl2)
          m2=kl2%nl
          kx=kxadaloc(-1,m2+2,klx)
          klx%head%k=ktfoper+mtfinequality
          klx%dbody(1)=dtfcopy(k1)
          klx%body(2)=ktfoper+iop1
          do i=1,m2
            klx%dbody(i+2)=dtfcopy(kl2%dbody(i))
          enddo
        else
          if(nextrel)then
            kx=kxadaloc(-1,3,klx)
            klx%head%k=ktfoper+mtfinequality
            klx%dbody(1)=dtfcopy(k1)
            klx%body(2)=ktfoper+iop1
            klx%dbody(3)=dtfcopy(k2)
          else
            kx=kxadaloc(-1,2,klx)
            klx%head%k=ktfoper+iop1
            klx%dbody(1)=dtfcopy(k1)
            klx%dbody(2)=dtfcopy(k2)
          endif
        endif
      endif
      if(.not. nextrel)then
        if(tfconstlistqo(klx))then
          call tfleval(klx,kx%k,.true.,irtc)
        endif
      endif
      return
      end

      subroutine tfinequality(isp1,kx,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) kx
      integer*8 ka
      integer*4 isp1,irtc,i,narg,itfmessage
      narg=isp-isp1
      if(narg .lt. 3 .or. narg .ne. (narg/2)*2+1)then
        irtc=itfmessage(9,'General::narg','"3 or more"')
        return
      endif
      do i=isp1+1,isp,2
        call tfeevalref(dtastk(i),dtastk(i),irtc)
        if(irtc .ne. 0)then
          return
        endif
        if(i .gt. isp1+1)then
          if(ktfnonoperqd(dtastk(i-1),ka) .or.
     $         ka .lt. mtfequal .or. ka .gt. mtfunsame)then
            irtc=itfmessage(9,'General::wrongtype',
     $           '"==, <>, <, >, <=, >=, ===, <=>"')
            return
          endif
          call tfeval1(dtastk(i-2),dtastk(i),kx,int(ka),irtc)
          if(irtc .ne. 0)then
            return
          endif
          if(kx%k .eq. 0)then
            return
          elseif(ktfnonrealq(kx))then
            irtc=-1
            return
          endif
        endif
      enddo
      return
      end

      subroutine tfcompose(isp1,kh,kx,irtc)
      use tfstk
      implicit none
c      include 'DEBUG.inc'
      type (sad_descriptor) kx,kh
      integer*8 kah
      integer*4 isp1,irtc,iah,isp0
      if(rlist(iaximmediate) .ne. 0.d0)then
        if(ktfoperqd(kh,kah))then
          iah=int(kah)
          if(iah .gt. mtfend)then
            call tfcomposefun(isp1,iah,kx,.false.,irtc)
          else
            call tfcomposeoper(isp1,iah,kx,.true.,isp0,irtc)
          endif
          if(irtc .eq. 0)then
            return
          endif
        elseif(ktfstringqd(kh))then
          if(isp .eq. isp1+1 .or. isp .eq. isp1+2)then
            if(ktfrealq(dtastk(isp1+1)) .and.
     $           ktfrealq(dtastk(isp)))then
              kx=kxsubstring(kh,isp1+1,isp)
              irtc=0
              return
            endif
          endif
        endif
      endif
      kx=kxcrelistm(isp-isp1,ktastk(isp1+1:isp),kh)
      irtc=0
      return
      end

      subroutine tfcomposefull(isp1,kh,kx,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) kh,kx
      integer*8 kah
      integer*4 isp1,i,irtc,isp0,iah
      logical*4 tfconstqk,tfheldqd
      if(isp .ne. isp1 .and. rlist(iaximmediate) .ne. 0.d0)then
        if(ktfoperqd(kh,kah))then
          iah=int(kah)
          if(iah .gt. mtfend)then
            call tfcomposefun(isp1,iah,kx,.true.,irtc)
            if(irtc .eq. 0)then
              return
            elseif(irtc .gt. 0)then
              go to 10
            endif
            if(.not. ktfimmediateq(klist(ifunbase+iah)))then
              go to 10
            endif
          else
            call tfcomposeoper(isp1,iah,kx,.true.,isp0,irtc)
            if(irtc .eq. 0)then
              return
            elseif(irtc .gt. 0)then
              if(iah .eq. mtfcomplex)then
                go to 1
              endif
              go to 10
            endif
          endif
        elseif(ktfstringqd(kh))then
          if(isp .eq. isp1+1 .or. isp .eq. isp1+2)then
            if(ktfrealq(ktastk(isp1+1)) .and. ktfrealq(ktastk(isp)))then
              kx=kxsubstring(kh,isp1+1,isp)
              irtc=0
              return
            endif
          endif
          go to 10
        else
          go to 10
        endif
 1      isp0=isp
        isp=isp+1
        dtastk(isp)=kh
        do i=isp1+1,isp0
          if(tfconstqk(ktastk(i)))then
            if(ktastk(i) .eq. ktfoper+mtfnull)then
              isp=isp0
              go to 10
            elseif(tfheldqd(dtastk(i)))then
              isp=isp0
              go to 10
            endif
            isp=isp+1
            ktastk(isp)=ktastk(i)
          else
            isp=isp0
            go to 10
          endif
        enddo
        call tfefunref(isp0+1,kx,.true.,irtc)
        isp=isp0
        return
      endif
 10   kx=kxcrelistm(isp-isp1,ktastk(isp1+1:isp),kh)
      irtc=0
      return
      end

      subroutine tfcomposeoper(isp1,iah,kx,comp,isp0,irtc)
      use tfstk
      use ophash
      implicit none
      type (sad_descriptor) kx
      type (sad_list), pointer :: klx
      integer*8 iee0
      integer*4 isp1,iah,irtc,isp0,i,iy,isp00,iv
      logical*4 tfconstqk,tfsameqk,tfsamestringqk,comp
      real*8 vx,x,y
      irtc=1
      if(constop(iah))then
        return
      elseif(isp1+2 .eq. isp)then
        select case (iah)
        case (mtfplus,mtftimes)
          if(tfconstqk(ktastk(isp1+1)))then
            if(tfconstqk(ktastk(isp)))then
              call tfplus(isp1,kx,iah,irtc)
              return
            endif
          endif
        case (mtfpower)
          if(ktfrealq(ktastk(isp)))then
            y=rtastk(isp)
            if(y .eq. 1.d0)then
              kx=dtastk(isp1+1)
              irtc=0
              return
            elseif(y .eq. 0.d0 .and.
     $             redmath%value%k .ne. 0)then
              kx%k=ktftrue
              irtc=0
              return
            elseif(iand(ktrmask,ktastk(isp1+1)) .ne. ktfnr)then
              x=rtastk(isp1+1)
              if(y .eq. -1.d0)then
                if(x .eq. 0.d0)then
                  return
                else
                  vx=1.d0/x
                endif
              elseif(y .eq. 2.d0)then
                vx=x**2
              elseif(y .eq. .5d0)then
                vx=sqrt(x)
              else
                iy=int(y)
                if(y-iy .eq. 0.d0)then
                  vx=x**iy
                else
                  vx=x**y
                endif
              endif
              kx=dfromr(vx)
              irtc=0
              return
            endif
          endif
        case (mtfequal)
          if(iand(ktrmask,ktastk(isp1+1)) .ne. ktfnr .and.
     $         iand(ktrmask,ktastk(isp)) .ne. ktfnr )then
            if(rtastk(isp1+1) .eq. rtastk(isp))then
              kx%k=ktftrue
            else
              kx%k=0
            endif
            irtc=0
            return
          elseif(ktftype(ktastk(isp1+1)) .eq. ktfstring .and.
     $           ktftype(ktastk(isp)) .eq. ktfstring)then
            if(tfsamestringqk(ktastk(isp1+1),ktastk(isp)))then
              kx%k=ktftrue
            else
              kx%k=0
            endif
            irtc=0
            return
          endif
        case (mtfunequal)
           if(iand(ktrmask,ktastk(isp1+1)) .ne. ktfnr .and.
     $         iand(ktrmask,ktastk(isp)) .ne. ktfnr )then
            if(rtastk(isp1+1) .eq. rtastk(isp))then
              kx%k=0
            else
              kx%k=ktftrue
            endif
            irtc=0
            return
          elseif(ktftype(ktastk(isp1+1)) .eq. ktfstring .and.
     $           ktftype(ktastk(isp)) .eq. ktfstring)then
            if(tfsamestringqk(ktastk(isp1+1),ktastk(isp)))then
              kx%k=0
            else
              kx%k=ktftrue
            endif
            irtc=0
            return
          endif
        case (mtfsame)
          if(tfconstqk(ktastk(isp1+1)) .and. tfconstqk(ktastk(isp)))then
            if(tfsameqk(ktastk(isp1+1),ktastk(isp)))then
              kx%k=ktftrue
            else
              kx%k=0
            endif
            irtc=0
            return
          endif
        case (mtfunsame)
          if(tfconstqk(ktastk(isp1+1)) .and. tfconstqk(ktastk(isp)))then
            if(tfsameqk(ktastk(isp1+1),ktastk(isp)))then
              kx%k=0
            else
              kx%k=ktftrue
            endif
            irtc=0
            return
          endif
        case (mtfgreater)
           if(iand(ktrmask,ktastk(isp1+1)) .ne. ktfnr .and.
     $         iand(ktrmask,ktastk(isp)) .ne. ktfnr )then
            if(rtastk(isp1+1) .gt. rtastk(isp))then
              kx%k=ktftrue
            else
              kx%k=0
            endif
            irtc=0
            return
          endif
        case (mtfgeq)
           if(iand(ktrmask,ktastk(isp1+1)) .ne. ktfnr .and.
     $         iand(ktrmask,ktastk(isp)) .ne. ktfnr )then
            if(rtastk(isp1+1) .ge. rtastk(isp))then
              kx%k=ktftrue
            else
              kx%k=0
            endif
            irtc=0
            return
          endif
        case(mtfless)
           if(iand(ktrmask,ktastk(isp1+1)) .ne. ktfnr .and.
     $         iand(ktrmask,ktastk(isp)) .ne. ktfnr )then
            if(rtastk(isp1+1) .lt. rtastk(isp))then
              kx%k=ktftrue
            else
              kx%k=0
            endif
            irtc=0
            return
          endif
        case (mtfleq)
           if(iand(ktrmask,ktastk(isp1+1)) .ne. ktfnr .and.
     $         iand(ktrmask,ktastk(isp)) .ne. ktfnr )then
            if(rtastk(isp1+1) .le. rtastk(isp))then
              kx%k=ktftrue
            else
              kx%k=0
            endif
            irtc=0
            return
          endif
        case (mtfpart)
          if(ktastk(isp1+1) .eq. klist(iaxslotnull))then
            if(ktfrealq(ktastk(isp)))then
              iv=int(rtastk(isp))
              if(iv .ge. -2 .and. iv .le. nslots)then
                kx=dlist(iaxslotpart+iv+2)
                irtc=0
                return
              endif
            endif
            return
          endif
        case (mtfatt)
          call tfclassmember(dtastk(isp1+1),dtastk(isp),kx,.false.,irtc)
c          call tfdebugprint(dtastk(isp1+1),'estk',1)
c          call tfdebugprint(dtastk(isp),'@',1)
c          call tfdebugprint(kx,'==>',1)
c          write(*,*)'with ',irtc
          return
        case (mtffun)
          if(ktastk(isp) .eq. ktfoper+mtfnull)then
            kx=kxpfaloc(dtastk(isp-1))
            irtc=0
          endif
          return
        case (mtfset,mtfsetdelayed,mtfmap,mtfapply,
     $         mtfmapall,mtfunset,mtfmessagename)
          return
        end select
      endif
      select case (iah)
      case (mtfpart)
        if(tfconstqk(ktastk(isp1+1)))then
          if(tflistq(ktastk(isp1+1)))then
            do i=isp1+2,isp
              if(ktfnonrealq(ktastk(i)))then
                if(ktastk(i) .ne. ktfoper+mtfnull)then
                  return
                endif
              endif
            enddo
            isp0=isp
            iee0=ierrorexp
            ierrorexp=1
            call tfpart(isp1+1,kx,.false.,irtc)
            ierrorexp=iee0
            isp=isp0
            if(irtc .ne. 0)then
              irtc=1
              return
            endif
          endif
        endif
        return
      case (mtfplus)
        vx=0.d0
        isp0=isp
        isp=isp+1
        do i=isp1+1,isp0
          if(ktfnonrealq(ktastk(i)))then
            isp=isp+1
            ktastk(isp)=ktastk(i)
          else
            vx=vx+rtastk(i)
          endif
        enddo
        if(isp .eq. isp0+1)then
          kx=dfromr(vx)
          irtc=0
          return
        elseif(vx .eq. 0.d0)then
          if(isp-isp0-1 .ne. isp0-isp1)then
            isp=isp0
            return
          endif
          isp00=isp0
          isp0=isp0+1
          go to 110
        else
          rtastk(isp0+1)=vx
          isp00=isp0
          go to 110
        endif
      case (mtftimes)
        vx=1.d0
        isp0=isp
        isp=isp+1
        do i=isp1+1,isp0
          if(ktfnonrealq(ktastk(i)))then
            isp=isp+1
            ktastk(isp)=ktastk(i)
          else
            vx=vx*rtastk(i)
          endif
        enddo
        if(isp .eq. isp0+1)then
          kx=dfromr(vx)
          isp=isp0
          irtc=0
          return
        elseif(vx .eq. 1.d0)then
          if(isp-isp0-1 .ne. isp0-isp1)then
            isp=isp0
            return
          endif
          isp00=isp0
          isp0=isp0+1
          go to 110
        else
          rtastk(isp0+1)=vx
          isp00=isp0
          go to 110
        endif
      case (mtfand)
        isp0=isp
        do i=isp1+1,isp
          if(ktfnonrealq(ktastk(i)))then
            isp=isp+1
            ktastk(isp)=ktastk(i)
          elseif(rtastk(i) .eq. 0.d0)then
            kx%k=0
            isp=isp0
            irtc=0
            return
          endif
        enddo
        if(isp .eq. isp0)then
          kx%k=ktftrue
          irtc=0
          return
        else
          isp00=isp0
          go to 110
        endif
      case (mtfor)
        isp0=isp
        do i=isp1+1,isp
          if(ktfnonrealq(ktastk(i)))then
            isp=isp+1
            ktastk(isp)=ktastk(i)
          elseif(rtastk(i) .ne. 0.d0)then
            kx%k=ktftrue
            isp=isp0
            irtc=0
            return
          endif
        enddo
        if(isp .eq. isp0)then
          kx%k=0
          irtc=0
          return
        else
          isp00=isp0
          go to 110
        endif
      case (mtfconcat)
        do i=isp1+1,isp
          if(ktfnonstringq(ktastk(i)))then
            go to 1
          endif
        enddo
        call tfstringjoin(isp1,kx,irtc)
        return
      case (mtfset,mtfsetdelayed,mtfmap,mtfapply,
     $       mtfmapall,mtfatt,mtfunset,mtffun,mtfmessagename)
        return
      end select
 1    irtc=-1
      return
 110  irtc=0
      if(isp .eq. isp0+1)then
        kx=dtastk(isp)
        if(ktfrealq(kx))then
          isp=isp00
          return
        elseif(ktflistq(kx,klx))then
          if(klx%head%k .eq. ktfoper+mtfcomplex)then
            isp=isp00
            return
          endif
        endif
      endif
      if(comp)then
        kx=kxmakelist(isp0,klx)
        klx%head%k=ktfoper+iah
        isp=isp00
      else
        kx%k=ktfref
      endif
      return
      end

      recursive subroutine tfcomposefun(isp1,iah,kx,full,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) kx,dh
      integer*8 ka,kti,kai,i
      integer*4 isp1,iah,irtc,narg,id
      real*8 rimmediate0
      logical*4 tfconstqk,full,re
      irtc=1
      id=iget_fun_id(int8(iah))
      select case (id)
      case (nfunif)
        narg=isp-isp1
        if(narg .ge. 2 .and. narg .le. 4)then
          if(ktastk(isp1+1) .eq. 0)then
            if(narg .ge. 3)then
              kx=dtastk(isp1+3)
            else
              kx%k=ktfoper+mtfnull
            endif
          elseif(ktfrealq(ktastk(isp1+1)))then
            kx=dtastk(isp1+2)
          else
            return
          endif
          irtc=0
          return
        endif
        return
      case (nfunmodule)
        if(full)then
 11       if(isp .eq. isp1+2)then
            call tfmodule(isp1,kx,.true.,.false.,irtc)
            if(irtc .ne. 0)then
              if(irtc .gt. 0 .and. ierrorprint .ne. 0)then
                call tfreseterror
              endif
              irtc=1
            endif
            isp=isp1+2
          elseif(isp .gt. isp1+2)then
            dh=dtastk(isp1+1)
            dtastk(isp1+1)=dtastk(isp1)
            call tfcomposefun(isp1+1,iah,kx,.true.,irtc)
            dtastk(isp1+1)=dh
            if(irtc .ne. 0)then
              if(irtc .gt. 0 .and. ierrorprint .ne. 0)then
                call tfreseterror
              endif
              irtc=1
              return
            endif
            isp=isp1+2
            dtastk(isp)=kx
            go to 11
          endif
        endif
        return

      case (nfunwith)
        if(full)then
 21       if(isp .eq. isp1+2)then
            kx=kxcompose(isp1)
            irtc=0
          elseif(isp .gt. isp1+2)then
            dh=dtastk(isp1+1)
            dtastk(isp1+1)=dtastk(isp1)
            call tfcomposefun(isp1+1,iah,kx,.true.,irtc)
            dtastk(isp1+1)=dh
            if(irtc .ne. 0)then
              if(irtc .gt. 0 .and. ierrorprint .ne. 0)then
                call tfreseterror
              endif
              irtc=1
              return
            endif
            isp=isp1+2
            dtastk(isp)=kx
            go to 21
          endif
          return
        endif

      case (nfunlength)
        if(isp .eq. isp1+1)then
          if(ktftype(ktastk(isp)) .eq. ktflist)then
            ka=iand(ktamask,ktastk(isp))
            if(klist(ka) .eq. ktfoper+mtflist)then
              if(iand(ilist(2,ka-3),lnonreallist) .ne. 0)then
                do i=ka+1,ka+ilist(2,ka-1)
                  kti=ktftype(klist(i))
                  if(kti .eq. ktfsymbol)then
                    return
                  elseif(kti .eq. ktflist)then
                    kai=iand(ktamask,klist(i))
                    if(klist(kai) .eq. ktfoper+mtfnull .or.
     $                   klist(kai) .eq. ktfoper+mtfslotseq)then
                      return
                    endif
                  endif
                enddo
              endif
              kx=dfromr(dble(ilist(2,ka-1)))
              irtc=0
            endif
          endif
        endif
        return
      end select
      irtc=-1
      if(ktfnumericq(klist(ifunbase+iah)))then
        rimmediate0=rlist(iaximmediate)
        rlist(iaximmediate)=-1.d0
        if(isp .eq. isp1+1)then
          if(ktfrealq(ktastk(isp)))then
            call tfefunref(isp1,kx,.false.,irtc)
          elseif(tfconstqk(ktastk(isp)))then
            if(ktflistq(ktastk(isp)))then
              if(ilist(2,ktfaddr(ktastk(isp))-1) .le. 15)then
c                call tfdebugprint(ktastk(isp),'composefun-a',1)
                call tfefunref(isp1,kx,.true.,irtc)
              endif
            else
c              call tfdebugprint(ktastk(isp),'composefun-b',1)
              call tfefunref(isp1,kx,.true.,irtc)
            endif
          endif
        else
          re=.false.
c          call tfdebugprint(ktastk(isp1),'composefun-c',1)
c          write(*,*)'with ',isp1,isp
          do i=isp1+1,isp
            if(ktfrealq(ktastk(i)))then
            elseif(tfconstqk(ktastk(i)))then
c              call tfdebugprint(ktastk(i),'composefun-c-i',1)
              re=.true.
            else
              rlist(iaximmediate)=rimmediate0
              return
            endif
          enddo
          call tfefunref(isp1,kx,re,irtc)
c          call tfdebugprint(kx,'==>',1)
        endif
        rlist(iaximmediate)=rimmediate0
      endif
      return
      end

      subroutine tfevalb(string,nc,kx,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) kx
      integer*4 irtc,istart,istop,nc
      character string(nc)
      istart=1
      levele=levele+1
c      write(*,*)'tfevalb ',string
      call tfeval(string,nc,istart,istop,kx,.false.,irtc)
      call tfconnect(kx,irtc)
      return
      end

      subroutine tfevalc(string,nc)
      use tfstk
      implicit none
      type (sad_descriptor) kx
      integer*4 irtc,istart,istop,l,nc,itfdownlevel
      character string(nc)
      istart=1
      levele=levele+1
c      write(*,*)'tfevalc ',string
      call tfeval(string,nc,istart,istop,kx,.false.,irtc)
      if(irtc .gt. 0 .and. ierrorprint .ne. 0)then
        call tfreseterror
      endif
      l=itfdownlevel()
      return
      end

      subroutine tfevals(string,kx,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) kx
      integer*4 irtc,istart,istop
      character*(*) string
      istart=1
      levele=levele+1
      call tfeval(string,len_trim(string),
     $     istart,istop,kx,.false.,irtc)
c      call tfdebugprint(kx,'==>',1)
      call tfconnect(kx,irtc)
      return
      end

      logical*4 function tfreadevalbuf(istart,istop,l,ipr)
      use tfcsi
      implicit none
      integer*4 ip1,istart,istop,l,ipr
      ip1=icsmrk()
      call cssetp(istop)
      call tprmpt(ipr,-1,0)
      call getbuf
      call tprmpt(0,-1,0)
      if(csrec())then
        call cssetp(ip1)
      else
        istop=icsmrk()
      endif
      tfreadevalbuf=icsstat() .eq. 0
      if(tfreadevalbuf)then
        istart=istop
        l=icslrecl()
      endif
      return
      end

      character*(*) function tfkname(k)
      use tfstk
      implicit none
      integer*8 k
      if(ktfrealq(k))then
        tfkname='Real'
      else
        select case(ktftype(k))
        case(ktfoper)
          tfkname='Function'
        case(ktflist)
          tfkname='List'
        case(ktfstring)
          tfkname='String'
        case(ktfsymbol)
          tfkname='Symbol'
        case(ktfpat)
          tfkname='Pattern'
        case(ktfref)
          tfkname='Reference'
        case default
          tfkname='Unknown'
        end select
      endif
      return
      end

      logical*4 function tfheldqd(k)
      use tfstk
      implicit none
      type (sad_descriptor) k
      type (sad_list), pointer :: kl
      tfheldqd=ktflistq(k,kl) .and.
     $     kl%head%k .eq. ktfoper+mtfhold
      return
      end

