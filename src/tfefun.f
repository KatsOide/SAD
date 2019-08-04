      subroutine tfefun(isp1,kx,ref,upvalue,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) kx
      integer*4 isp1,irtc
      logical*4 ref,upvalue
      if(ref)then
        call tfefunref(isp1,kx,upvalue,irtc)
      else
        call tfefundef(isp1,kx,irtc)
      endif
      return
      end

      subroutine tfefunrefc(isp1,kx,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) kx
      integer*4 isp1,irtc
      levele=levele+1
      call tfefunref(isp1,kx,.true.,irtc)
      call tfconnect(kx,irtc)
      return
      end

      subroutine tfefunrefstk(isp1,isp2,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) kx
      type (sad_dlist), pointer :: kl
      integer*4 isp1,irtc,isp2
      call tfefunref(isp1,kx,.true.,irtc)
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

      recursive subroutine tfefunref(isp1,kx,upvalue,irtc)
      use tfstk
      use tfmem
      use tfshare
      use tfcsi
      use mathfun
      implicit none
      type (sad_descriptor) kx,k1,k,kh
      type (sad_dlist), pointer :: kl,kl1,klx,klh
      type (sad_symbol), pointer :: sym1
      type (sad_symdef), pointer :: symd
      type (sad_string), pointer :: str
      integer*8 ka1,ka,kax,km,kop
      integer*4 isp1,irtc,i,id,narg,nc,isp2,itfdepth,
     $     itfgetrecl,ltr0,iaf,itfopenwrite,itfmessage,
     $     itfopenappend,itfopcode,isp0,nsize
      real*8 vx,v1,v
      character*3 opcx
      character*8 char8
      logical*4 euv,upvalue,tfnearlysameqf,rep
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

      real*8 aloggamma1,factorial,gammaq,gammap,inverseerf,
     $     productlog,gamma0,erf,erfc
      complex*16 cloggamma1,cfactorial,cerfc,cerf,
     $     cproductlog
      external aloggamma1,cloggamma1,factorial,cfactorial,gammaq,
     $     gammap,cerfc,cerf,gamma0,erf,erfc,inverseerf,
     $     cproductlog,productlog
c      call tfreecheck1('tfefun-0',itastk(1,isp),
c     $     itastk(2,isp),vstk(ivstkoffset+isp),irtc)
c      if(irtc .ne. 0)then
c        call tfdebugprint(itastk(1,isp1),itastk(2,isp1),
c     $       vstk(ivstkoffset+isp1),'tfefun-0',1)
c      endif
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
              call tfsolvemember(klh,k1,rep,irtc)
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
              call tfsydef(sym1,sym1)
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
c
     $         610, 620, 630, 640, 650, 660, 670, 680, 690, 700,
     $         700, 720, 730, 730, 730, 760, 770, 780, 790, 800,
     $         810, 820, 830, 840, 850, 860, 870, 880, 890, 900,
     $         910, 920, 930, 940, 950, 960, 970, 980, 990,1000,
     $        1010,1020,1030,1040,1050,1060,1070,1080,1090,1100,
c
     $        1110,1120,1130,1140,1150,1160,1170,1180,1190,1200,
     $        1210,1220,1230,1240,1250,1260,1270,1280,1290,1300,
     $        1310,1320,1330,1340,1350,1360,1370,1380,1380,1380,
     $        1290,1420,1430,1440,1450,1460,1470,1480,1490,1500,
     $        1510,1520,1530,1540,1540,1560,1570,1580,1590,1600,
c
     $        1610,1620, 760,1640,1650,1660,1670,1680,1690,1700,
     $        1710,1720,1730,1740,1750,1760,1770,1770,1770,1770,
     $        1810,1820,1830,1840,1850,1860,1870,1870,1870,1900,
     $        1910,1920,1930,1940,1950,1960,1970,1980,1980,1980,
     $        2010,2020,2030,2040,2050,2060,2070,2080,2080,2100,
c
     $        2100,2100,2100,2140,2150,2160,2170,2180,2190,2200,
     $        2210,2220,2230,2240,2250,2260,2270,2280,2290,2300,
     $        2310,2320,2330,2340,2350,2360,2370,2380
     $       ),id
c              Sin  Cos  Tan  Sinh Cosh Tanh Exp  Log  Atan Det
c              Sqrt Flor Ceil Min  Max  Mod  StrL Arg  Sign Leng
c              Dims RplP ASin ACos ASh  ACh  ATh  Tabl Do   Attr
c              Peek Abs  Revs Modl Blck StrR SwiC Fltn If   Take
c              Selc Whil Join Apnd Ppnd Clr  Prot Unpr Drop MpAt
c
c              Innr Trps SglV DiaM LinS IdtN Eigs Oper Posi Sum
c              Prod Rang Re   Im   Conj ToSt Dpth Levl Writ Get
c              OpnW OpnA Clos Flsh Prnt WStr Retn Head RLiQ Pttn
c              Thrw Ctch Thrd SetA MpId FrCh ToCh CmpQ Tr   SvSM
c              Stch Sort Uni1 Ordr MChk Scan Iden TimU NumQ VecQ
c
c              AtmQ Outr MatQ TrcP Defi READ Ints Cmpl Roun IErf
c              FrDt ToDt ToIS ReaS OpnR ToEx StrM StrP ToUp Brek
c              Cont Goto Four IFou Chek Whic Syst GetP GetU GetG
c              ToLo Unev Case DelC Vect BDPi Nams GbCl LgGm Fact
c              With WhiC Ovrr AppT PreT FndR GamR GmRP Erf  Erfc
c
c              Fit  Symb SyNm Extr Read Skip TmpN Exit StrF Rstr
c              MM   Shrt Fork Dir  SDir Wait BesJ BesY BesI BesK
c              BasF StrT S2S  Even OddQ DatS Inst Delt FlAt Repl
c              SetE Spl$ FInd SCnt SCnP ToCt Ctxt BAnd BOr  BXor
c              RepM MSca StdF Abrt ChkA RelH NaNQ MapT ScaT Last
c
c              Frst Scnd Thrd ObjS PrdL GauC ClMO MAll Dupl GCLn
c              Seek DigQ LetQ ReaQ NSmQ OpSh RdSh WrSh ShSz FBuQ
c              GaCU GaCF Rest RRt1 Diff Gam0 **** XSin
        if(id .gt. 0)then
          if(id .le. 2000)then
            call tfefun1(isp1,id,kx,.true.,irtc)
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
          call tfeintf(dsin,cdsin,k,kx,
     $         .true.,-dinfinity,dinfinity,irtc)
        else
          go to 6810
        endif
        go to 6900
 120    if(narg .eq. 1)then
          call tfeintf(dcos,cdcos,k,kx,
     $         .true.,-dinfinity,dinfinity,irtc)
        else
          go to 6810
        endif
        go to 6900
 130    if(narg .eq. 1)then
          call tfeintf(dtan,tctan,k,kx,
     $         .true.,-dinfinity,dinfinity,irtc)
        else
          go to 6810
        endif
        go to 6900
 140    if(narg .eq. 1)then
          call tfeintf(dsinh,tcsinh,k,kx,
     $         .true.,-dinfinity,dinfinity,irtc)
        else
          go to 6810
        endif
        go to 6900
 150    if(narg .eq. 1)then
          call tfeintf(dcosh,tccosh,k,kx,
     $         .true.,-dinfinity,dinfinity,irtc)
        else
          go to 6810
        endif
        go to 6900
 160    if(narg .eq. 1)then
          call tfeintf(dtanh,tctanh,k,kx,
     $         .true.,-dinfinity,dinfinity,irtc)
        else
          go to 6810
        endif
        go to 6900
 170    if(narg .eq. 1)then
          call tfeintf(dexp,cdexp,k,kx,
     $         .true.,-dinfinity,dinfinity,irtc)
        else
          go to 6810
        endif
        go to 6900
 180    if(narg .eq. 1)then
          call tfeintf(dlog,cdlog,k,kx,
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
          call tfeintf(datan,tcatan,k,kx,
     $         .true.,-dinfinity,dinfinity,irtc)
        elseif(narg .eq. 2)then
          call tfeintf2(datan2,tcatan2,k,
     $         ktastk(isp-1),.true.,kx,irtc)
        else
          go to 6812
        endif
        go to 6900
 200    call tfdet(isp1,kx,irtc)
        go to 6900
 210    if(narg .eq. 1)then
          call tfeintf(dsqrt,cdsqrt,k,kx,
     $         .true.,0.d0,dinfinity,irtc)
        else
          go to 6810
        endif
        go to 6900
 220    if(narg .eq. 1)then
          call tfeintf(tfloor,tcfloor,k,kx,
     $         .true.,-dinfinity,dinfinity,irtc)
        else
          go to 6810
        endif
        go to 6900
 230    if(narg .eq. 1)then
          call tfeintf(tceiling,tcceiling,k,kx,
     $         .true.,-dinfinity,dinfinity,irtc)
        else
          go to 6810
        endif
        go to 6900
 240    call tfminmax(isp1,kx,id-13,irtc)
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
          call tfeintf(tfarg,tfcarg,k,
     $         kx,.false.,-dinfinity,dinfinity,irtc)
        else
          go to 6810
        endif
        go to 6900
 290    if(narg .eq. 1)then
          call tfeintf(tfsign,tfcsign,k,
     $         kx,.true.,-dinfinity,dinfinity,irtc)
        else
          go to 6810
        endif
        go to 6900
 300    if(narg .ne. 1)then
          go to 6810
        endif
        if(ktflistq(k,kl))then
          kx=dfromr(dble(kl%nl))
        else
          kx%k=0
        endif
        go to 8000
 310    call tfdimensions(isp1,kx,irtc)
        go to 6900
 320    call tfreplacepart(isp1,kx,0,irtc)
        go to 6900
 330    if(narg .eq. 1)then
          call tfeintf(dasin,tcasin,k,
     $         kx,.true.,-1.d0,1.d0,
     $         irtc)
        else
          go to 6810
        endif
        go to 6900
 340    if(narg .eq. 1)then
          call tfeintf(dacos,tcacos,k,
     $         kx,.true.,-1.d0,1.d0,
     $         irtc)
        else
          go to 6810
        endif
        go to 6900
 350    if(narg .eq. 1)then
          call tfeintf(tasinh,tcasinh,k,
     $         kx,.true.,-dinfinity,dinfinity,irtc)
        else
          go to 6810
        endif
        go to 6900
 360    if(narg .eq. 1)then
          call tfeintf(tacosh,tcacosh,k,
     $         kx,.true.,1.d0,dinfinity,
     $         irtc)
        else
          go to 6810
        endif
        go to 6900
 370    if(narg .eq. 1)then
          call tfeintf(tatanh,tcatanh,k,
     $         kx,.true.,-1.d0,1.d0,
     $         irtc)
        else
          go to 6810
        endif
        go to 6900
 380    isp2=isp
        call tftable(isp1,isp1+2,isp2,kx,29-id,irtc)
        go to 6900
 400    call tfattributes(isp1,kx,irtc)
        go to 6900
 410    if(narg .ne. 1)then
          go to 6810
        elseif(ktfnonrealq(k))then
          irtc=itfmessage(9,'General::wrongtype','"Real number"')
        else
          ka=int8(rtastk(isp))
          if(.not. tfchecklastp(ka))then
            irtc=itfmessage(9,'General::wrongnum',
     $           '"within allocated block"')
          else
            kx=kxadaloc(-1,4,klx)
            klx%rbody(1)=dble(klist(ka))
            klx%rbody(2)=dble(ilist(1,ka))
            klx%rbody(3)=dble(ilist(2,ka))
            klx%dbody(4)=kxsalocb(0,transfer(klist(ka),char8),8)
            irtc=0
          endif
        endif
        go to 6900
 420    if(narg .eq. 1)then
          call tfeintf(dabs,ccdabs,k,kx,
     $         .false.,-dinfinity,dinfinity,irtc)
        else
          go to 6810
        endif
        go to 6900
 430    call tfreverse(isp1,kx,irtc)
        go to 6900
 440    call tfmodule(isp1,kx,id .eq. nfunmodule,.true.,irtc)
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
          call tftake(ktastk(isp-1),k,kx,.true.,.true.,irtc)
        endif
        go to 6900
 510    call tfselect(isp1,kx,irtc)
        go to 6900
 520    call tfwhile(isp1,kx,irtc)
        go to 6900
 530    call tfjoin(isp1,kx,.true.,irtc)
        go to 6900
 540    if(narg .eq. 2)then
          call tfappend(ktastk(isp1+1),k,kx,.true.,0,irtc)
        elseif(narg .ne. 1)then
          go to 6812
        endif
        go to 6900
 550    if(narg .eq. 2)then
          call tfappend(ktastk(isp1+1),k,kx,.true.,1,irtc)
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
          call tftake(ktastk(isp-1),k,kx,.false.,.true.,irtc)
        endif
        go to 6900
 600    call tfreplacepart(isp1,kx,1,irtc)
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
        call tftable(isp1,isp1+2,isp2,kx,id-58,irtc)
        go to 6900
 720    call tfrange(isp1,kx,irtc)
        go to 6900
 730    if(narg .eq. 1)then
          call tfcmplxf(k,kx,id-62,ka1)
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
 790    call tfwrite(isp1,kx,irtc)
        go to 6900
 800    if(narg .ne. 1)then
          go to 6810
        else
          call tfget(k,kx,irtc)
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
            call tfaddmessage(' ',0,icslfno())
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
          if(tfreallistq(k))then
            kx%k=ktftrue
          else
            kx%k=0
          endif
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
 950    call tfmap(isp1,kx,4,1,irtc)
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
 1060   call tfmap(isp1,kx,1,1,irtc)
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
          if(tfnumberq(k))then
            kx%k=ktftrue
          else
            kx%k=0
          endif
          go to 8000
        else
          go to 6810
        endif
        go to 6900
 1100   call tfvectorqf(isp1,kx,irtc)
        go to 6900
 1110   if(narg .eq. 1)then
          if(ktflistq(k,kl))then
            if(kl%head%k .eq. ktfoper+mtfcomplex)then
              kx%k=ktftrue
            else
              kx%k=0
            endif
          else
            kx%k=ktftrue
          endif
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
          call tfeevalref(k,kx,irtc)
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
          call tfeintf(tround,tcround,k,
     $         kx,.true.,-dinfinity,dinfinity,irtc)
        else
          go to 6810
        endif
        go to 6900
 1200   if(narg .eq. 1)then
          call tfeintf(inverseerf,inverseerf,k,
     $         kx,.true.,-dinfinity,dinfinity,irtc)
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
 1350   call tfcheck(isp1,kx,irtc)
        go to 6900
 1360   call tfwhich(isp1,kx,irtc)
        go to 6900
 1370   irtc=itfmessage(999,'General::unregister',' ')
        go to 6900
 1380   irtc=itfmessage(999,'General::unregister',' ')
        go to 6900
 1420   call tfsequence(isp1,isp,kx)
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
          call tfeintf(aloggamma1,cloggamma1,k,kx,
     $         .true.,0.d0,dinfinity,irtc)
        else
          go to 6810
        endif
        go to 6900
 1500   if(narg .eq. 1)then
          call tfeintf(factorial,cfactorial,k,kx,
     $         .true.,-1.d0+5.56d-17,dinfinity,irtc)
        else
          go to 6810
        endif
        go to 6900
 1510   call tfwith(isp1,kx,.true.,irtc)
        go to 6900
 1520   call tfswitchcases(isp1,kx,1,irtc)
        go to 6900
 1530   call tfoverride(isp1,kx,irtc)
        go to 6900
 1540   call tfappendto(isp1,kx,id-143,irtc)
        go to 6900
 1560   call tffindroot(isp1,kx,irtc)
        go to 6900
 1570   if(narg .eq. 2)then
          call tfeintf2(gammaq,0.d0,ktastk(isp-1),k,.false.,kx,irtc)
        else
          go to 6812
        endif
        go to 6900
 1580   if(narg .eq. 2)then
          call tfeintf2(gammap,0.d0,ktastk(isp-1),k,.false.,kx,irtc)
        else
          go to 6812
        endif
        go to 6900
 1590   if(narg .eq. 1)then
          call tfeintf(erf,cerf,k,
     $         kx,.true.,-dinfinity,dinfinity,irtc)
        else
          go to 6810
        endif
        go to 6900
 1600   if(narg .eq. 1)then
          call tfeintf(erfc,cerfc,k,
     $         kx,.true.,-dinfinity,dinfinity,irtc)
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
 1670   call tftemporaryname(isp1,kx,irtc)
        go to 6900
 1680   call tfexit(isp1,kx,irtc)
        go to 6900
 1690   call tfstringfill(isp1,kx,irtc)
        go to 6900
 1700   call tfrestrict(isp1,kx,irtc)
        go to 6900
 1710   call tfminmax(isp1,kx,0,irtc)
        go to 6900
 1720   call tfshort(isp1,kx,irtc)
        go to 6900
 1730   irtc=itfmessage(999,'General::unregister',' ')
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
          call tfeintf(tfevenq,tfcevenq,k,kx,
     $         .true.,-dinfinity,dinfinity,irtc)
        else
          go to 6810
        endif
        go to 6900
 1850   if(narg .eq. 1)then
          call tfeintf(tfoddq,tfcoddq,k,kx,
     $         .true.,-dinfinity,dinfinity,irtc)
        else
          go to 6810
        endif
        go to 6900
 1860   irtc=itfmessage(999,'General::unregister',' ')
        go to 6900
 1870   call tfreplacepart(isp1,kx,id-175,irtc)
        go to 6900
 1900   if(narg .ne. 2)then
          go to 6812
        endif
        call tfreplace(ktastk(isp1+1),ktastk(isp),
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
          call tfeintf(productlog,cproductlog,k,kx,
     $         .true.,-exp(-1.d0),dinfinity,irtc)
        else
          go to 6810
        endif
        go to 6900
 2160   call tfgaussiancoulomb(isp1,kx,irtc)
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
 2200   call tfgetcommandline(isp1,kx,irtc)
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
        if(tfnearlysameqf(dtastk(isp1+1),dtastk(isp1+2),
     $       rtastk(isp1+3),rtastk(isp)))then
          kx%k=ktftrue
        else
          kx%k=0
        endif
        go to 6900
 2260   call tfopenshared(isp1,kx,irtc)
        go to 6900
 2270   call tfreadshared(isp1,kx,irtc)
        go to 6900
 2280   call tfwriteshared(isp1,kx,irtc)
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
        if(ktfoperq(k))then
           kx%k=ktftrue
        else
           kx%k=0
        endif
        go to 6900
 2310   call tfgaussiancoulombu(isp1,kx,irtc)
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
          call tfeintf(gamma0,tfdummy,k,
     $         kx,.false.,-dinfinity,dinfinity,irtc)
        else
          go to 6810
        endif
        go to 6900
 2370   go to 6900
 2380   if(narg .eq. 1)then
          call tfeintf(xsin,tcxsin,k,
     $         kx,.true.,-dinfinity,dinfinity,irtc)
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
        go to (
     $       6010,
     $       7000,7000,6600,7000,6600,7000,6620,6630,6690,6690,
     $       6680,6680,6680,6680,6700,6700,6610,6100,6100,6510,
     $       7000,7000,6010,7000,6640,6650,6750,7000,7000,7000,
     $       7000,6002,6010,6010,6010,6660,6670,6560,6560,6710,
     $       6010,6740,7000,7000,6200,6090,6520,6530,6540,6010,
     $       6010,6550,6570,6570,6570,6570,6580,6580,6590,6720,
     $       6730,6010,7000,7000,6010,7000),
     $       iaf+1
c            null
c            m    i    +    -    *    /    v    ^    e    n    
c            >    <    g    l    E    N    ~    &&   o    c
c            [    ]    {    }    s    =    C    (    )    ,
c            ;    &    :    r    d    RepA RepR u    U    S    
c            ?    f    #    ##   .    |    M    MA   A    rept 
c            repn ineq AT   SF   TB   DB   INc  Dec  Part @
c            msgn TagS (*   *)   Hold z
        go to 6001
 6001   irtc=itfmessage(999,'General::invop',' ')
        return
 6002   if(isp .eq. isp1+2)then
          if(ktastk(isp) .eq. ktfoper+mtfnull)then
            kx=kxpfaloc(dtastk(isp-1))
            irtc=0
            return
          endif
        endif
 6010   kx=kxcompose(isp1)
        irtc=0
        return
 6090   if(narg .lt. 2)then
          irtc=-1
          go to 6900
        endif
 6100   if(narg .eq. 2)then
          call tfeval1(ktastk(isp1+1),ktastk(isp),kx,iaf,irtc)
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
          call tfeval1(k1,ktastk(i),kx,iaf,irtc)
          if(irtc .ne. 0)then
            go to 6900
          endif
        enddo
        return
 6200   if(narg .eq. 2)then
          call tfeval1(ktastk(isp1+1),ktastk(isp),kx,iaf,irtc)
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
          call tfeval1(ktastk(i),k1,kx,iaf,irtc)
          if(irtc .ne. 0)then
            go to 6900
          endif
        enddo
        return
 6510   call tfstringjoin(isp1,kx,irtc)
        go to 6900
 6520   call tfmap(isp1,kx,3,1,irtc)
        go to 6900
 6530   call tfmapall(isp1,kx,irtc)
        go to 6900
 6540   call tfapply(isp1,kx,irtc)
        go to 6900
 6550   call tfinequality(isp1,kx,irtc)
        go to 6900
 6560   if(narg .ne. 2)then
          go to 6812
        endif
        call tfupset(dtastk(isp1+1),dtastk(isp),int8(0),kx,irtc)
        go to 6900
 6570   if(narg .ne. 2)then
          go to 6812
        endif
        call tfeval1to(ktastk(isp1+1),ktastk(isp),
     $       kx,iaf,.false.,irtc)
        go to 6900
 6580   if(iaf .eq. mtfincrement)then
          v1=1.d0
        else
          v1=-1.d0
        endif
        if(narg .eq. 1)then
          call tfeval1to(ktastk(isp1+1),dfromr(v1),
     $         kx,mtfaddto,.true.,irtc)
        elseif(narg .eq. 2 .and.
     $         ktastk(isp1+1) .eq. ktfoper+mtfnull)then
          call tfeval1to(ktastk(isp),dfromr(v1),
     $         kx,mtfaddto,.false.,irtc)
        else
          go to 6810
        endif
        go to 6900
 6590   if(ktflistq(ktastk(isp1+1)))then
          call tfpart(isp1+1,kx,.true.,irtc)
          if(irtc .eq. 0)then
            call tfeevalref(kx,kx,irtc)
          endif
        elseif(ktfsymbolq(ktastk(isp1+1)))then
          irtc=-1
        else
          irtc=itfmessage(9,'General::wrongtype',
     $         '"List or composition"')
        endif
        go to 6900
 6600   call tfplus(isp1,kx,iaf,irtc)
        go to 6900
 6610   call tfnot(isp1,kx,iaf,irtc)
        go to 6900
 6620   call tfrevpower(isp1,kx,irtc)
        go to 6900
 6630   call tfpower(isp1,kx,irtc)
        go to 6900
 6640   if(narg .ne. 2)then
          go to 6812
        endif
 6650   call tfset(isp1,kx,.true.,irtc)
        go to 6900
 6660   call tfreplace1(isp1,kx,irtc)
        go to 6900
 6670   call tfreplacerepeated1(isp1,kx,irtc)
        go to 6900
 6680   call tfrelation(isp1,kx,iaf,irtc)
        go to 6900
 6690   call tfequal(isp1,kx,iaf,irtc)
        go to 6900
 6700   call tfsameq1(isp1,kx,iaf,irtc)
        go to 6900
 6710   if(narg .ne. 1)then
          go to 6810
        endif
        isp=isp1+2
        ktastk(isp)=ktfref
        call tfset(isp1,kx,.true.,irtc)
        isp=isp1+1
        kx%k=ktfoper+mtfnull
        go to 6900
 6720   call tfatt(isp1,kx,.true.,irtc)
        go to 6900
 6730   km=klist(ifunbase+mtfmessagename)
        call tfdeval(isp1,km,kx,1,.false.,euv,irtc)
        go to 6900
 6740   call tfflagordef(isp1,kx,irtc)
        go to 6900
 6750   if(narg .ne. 2)then
          go to 6812
        endif
        call tfeval1(ktastk(isp1+1),ktastk(isp),kx,mtfcomplex,irtc)
        go to 6900
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
            call tfpuref(isp1,kl1,kx,irtc)
            go to 6900
          case (-mtfnull)
            if(kl1%nl .eq. 0)then
              call tfsequence(isp1,isp,kx)
              if(ktflistq(kx,klx))then
                call tfleval(klx,kx,.true.,irtc)
              elseif(ktfsymbolq(kx) .or. ktfpatq(kx))then
                call tfeevalref(kx,kx,irtc)
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
            call tfpart(isp1,kx,.true.,irtc)
            if(irtc .eq. 0)then
              if(ktflistq(kx,klx))then
                call tfleval(klx,kx,.true.,irtc)
              elseif(ktfsymbolq(kx) .or. ktfpatq(kx))then
                call tfeevalref(kx,kx,irtc)
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
              call tfefunref(isp2,kx,upvalue,irtc)
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
              call tfefunref(isp2,kx,upvalue,irtc)
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
              call tfefunref(isp2,kx,upvalue,irtc)
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
      kx=kxcrelistm(narg,ktastk(isp1+1:isp1+narg),dtastk(isp1))
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

      subroutine tfsequence(isp1,isp2,kx)
      use tfstk
      implicit none
      type (sad_descriptor) kx
      integer*4 isp1,isp2
      if(isp1 .ge. isp2)then
        kx=dxnull
      elseif(isp1+1 .eq. isp2)then
        kx=dtastk(isp2)
      else
        kx=kxcrelistm(isp2-isp1,ktastk(isp1+1:isp2),
     $       k_descr(ktfoper+mtfnull))
      endif
      return
      end

      subroutine tfefundef(isp1,kx,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) kx,k1
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
          call tfefun1(isp1,id,kx,.false.,irtc)
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
          call tfpuref(isp1,kl,kx,irtc)
          go to 6900
        elseif(kl%head%k .eq. ktfoper+mtfnull)then
          if(kl%nl .eq. 0)then
            call tfsequence(isp1,isp,kx)
            if(ktflistq(kx,klx))then
              call tfleval(klx,kx,.false.,irtc)
            elseif(ktfsymbolq(kx,sym))then
              if(sym%override .eq. 0)then
                call tfsydef(sym,sym)
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
          call tfsolvemember(kli,k1,rep,irtc)
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
      implicit none
      type (sad_descriptor) kx,k1,k,ki
      type (sad_dlist), pointer :: klx
      integer*4 isp1,irtc,i,iopc,narg
      real*8 v,v1,vx,vi
      narg=isp-isp1
      if(narg .eq. 2)then
        k1=dtastk(isp1+1)
        k =dtastk(isp)
        if(ktfrealq(k1,v1) .and. ktfrealq(k,v))then
          if(iopc .eq. mtfplus)then
            kx=dfromr(v1+v)
          else
            kx=dfromr(v1*v)
          endif
          irtc=0
        else
          if(tfnumberq(k) .and. tfnumberq(k1))then
            call tfcmplx(k1,k,kx,iopc,irtc)
          elseif(tflistq(k1) .or. tflistq(k))then
            call tfecmplxl(k1,k,kx,iopc)
            irtc=0
          else
            call tfeexpr(k1,k,kx,iopc)
            irtc=0
          endif
        endif
        return
      elseif(narg .eq. 0)then
        if(iopc .eq. mtfplus)then
          kx%k=0
        else
          kx%k=ktftrue
        endif
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
            if(iopc .eq. mtfplus)then
              kx=dfromr(vx+vi)
            else
              kx=dfromr(vx*vi)
            endif
          else
            k1=kx
            if(tfnumberq(k1) .and. tfnumberq(ki))then
              call tfcmplx(k1,ki,kx,iopc,irtc)
              if(irtc .ne. 0)then
                return
              endif
            elseif(tflistq(k1) .or. tflistq(ki))then
              call tfecmplxl(k1,ki,kx,iopc)
              if(irtc .ne. 0)then
                return
              endif
            else
              call tfeexpr(k1,ki,kx,iopc)
            endif
          endif
        enddo
      endif
      return
      end

      subroutine tfpower(isp1,kx,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) kx,k1,ki
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
          call tfecmplxl(ki,k1,kx,mtfpower)
          irtc=0
        else
          call tfeexpr(ki,k1,kx,mtfpower)
        endif
      enddo
      return
      end

      subroutine tfrevpower(isp1,kx,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) kx,k1,ki
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
          call tfecmplxl(ki,k1,kx,mtfpower)
          irtc=0
        else
          call tfeexpr(ki,k1,kx,mtfpower)
        endif
      enddo
      return
      end

      recursive subroutine tfset(isp1,kx,upvalue,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) kx,k1,k2,k110,k11
      type (sad_symbol), pointer :: sym
      type (sad_symdef), pointer :: symd
      type (sad_dlist), pointer :: kl11
      integer*4 isp1,irtc,i,itfmessage,isp11,isp0
      logical*4 euv,upvalue
      if(isp1+1 .ge. isp)then
        irtc=itfmessage(9,'General::narg','"2 or more"')
        return
      endif
      k2=dtastk(isp)
      if(isp1+2 .eq. isp)then
        k1=dtastk(isp1+1)
        if(ktfsymbolq(k1,sym))then
          if(sym%override .eq. 0)then
            call tfsydef(sym,sym)
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
              call tfsydef(sym,sym)
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
        call tfset1(k1,k2,kx,ktfaddr(ktastk(isp1)),irtc)
      else
        isp0=isp
        kx=dtastk(isp)
        do i=isp-1,isp1+1,-1
          isp=isp0+3
          dtastk(isp-2)=dtastk(isp1)
          dtastk(isp-1)=dtastk(i)
          dtastk(isp  )=kx
          call tfset(isp0+1,kx,upvalue,irtc)
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
      implicit none
      type (sad_descriptor) kx
      integer*4 isp1,irtc,iopc,itfmessage,isp0,k
      if(isp .lt. isp1+2)then
        irtc=itfmessage(9,'General::narg','"2 or more"')
        return
      endif
      if(isp .eq. isp1+2)then
        if(tfnumberq(dtastk(isp1+1)) .and.
     $       tfnumberq(dtastk(isp)))then
          call tfcmplx(ktastk(isp1+1),ktastk(isp),kx,iopc,irtc)
        elseif(tflistq(dtastk(isp1+1))
     $         .or. tflistq(dtastk(isp)))then
          call tfearray(ktastk(isp1+1),ktastk(isp),kx,iopc,irtc)
        else
          call tfeexpr(ktastk(isp1+1),ktastk(isp),kx,iopc)
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
     $     ktfstringq(ktastk(isp)) .and. ktfstringq(ktastk(isp1+1)))then
        if(tfsamestringq(ktastk(isp),ktastk(isp1+1)))then
          kx%k=ktftrue
        else
          kx%k=0
        endif
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
      implicit none
      type (sad_descriptor) kx
      integer*4 isp1,irtc,iopc,itfmessage
      if(isp .ne. isp1+1)then
        irtc=itfmessage(9,'General::narg','"1"')
        return
      endif
      if(tfnumberq(dtastk(isp)))then
        call tfcmplx(0.d0,ktastk(isp),kx,iopc,irtc)
      elseif(tflistq(dtastk(isp)))then
        call tfearray(0.d0,ktastk(isp),kx,iopc,irtc)
      else
        call tfeexpr(0.d0,ktastk(isp),kx,iopc)
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
            call tfsydef(symd%sym,symi)
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
