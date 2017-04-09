      module maccode
c Do not forget to update sim/MACCODE.h when you change this module!!!!
      integer*4 , parameter ::
     $     icNULL   =   0, icDRFT   =   1,
     $     icBEND   =   2, icQUAD   =   4, icSEXT   =   6,
     $     icOCTU   =   8, icDECA   =  10, icDODECA =  12,
     $     icUND    =  18, icWIG    =  19, icSOL    =  20,
     $     icST     =  21, icMULT   =  22, icTEST   =  30,
     $     icCAVI   =  31, icTCAV   =  32, icMAP    =  33,
     $     icINS    =  34, icCOORD  =  35, icBEAM   =  36,
     $     icPROT   =  37, icSPCH   =  38,
     $     icMARK   =  41, icMONI   =  42, icAPRT   =  43,
     $     icMXEL   =  99, icLINE   = 100, icCELL   = 100,

     $     icRSVD = 258,
     $     icDEF  = icRSVD + 2, icACT  = icDEF  + 2,
     $     icPROC = icACT  + 2, icVAR  = icPROC + 2,
     $     icKWRD = icVAR  + 2, icUNIT = icKWRD + 2,
     $     icRAND = icUNIT + 2, icENV  = icRAND + 2,
     $     icFLAG = icENV  + 2, icGLI  = icFLAG + 4,
     $     icGLL  = icGLI  + 4, icGLR  = icGLI  + 8,
     $     icGraf = icGLL  + 2, icPART = icGraf + 2
      end module

      module maccbk
      implicit none
      integer*4 FLAGON,FLAGOF
      parameter (FLAGON=-1,FLAGOF=0)
      integer*4 HTMAX,MAXPNAME,LILISTDUMMY
      parameter(MAXPNAME=32,LILISTDUMMY=3)
      character*(MAXPNAME) NULSTR
      parameter(HTMAX=2**16,NULSTR='        ')

      integer*4 pagesz,inipage
      parameter(pagesz=4096/8,inipage=4)
      integer*4 MAXSTR,MAXMEM0,MAXMEM
      parameter (MAXSTR=256)
      parameter (MAXMEM=2*inipage*pagesz)
      parameter (MAXMEM0=6*1024*pagesz)
      character*(MAXPNAME) pname(HTMAX)
      integer*4 idtype(HTMAX)
      integer*8 idval(HTMAX)
      integer*8 ilistroot
      integer*4, parameter :: klistlen=16
      integer*8, pointer, dimension(:) :: klist
      real*8, pointer, dimension(:) :: rlist
      integer*4, pointer, dimension(:,:) :: ilist
      integer*1, pointer, dimension(:,:)  :: jlist

      interface sethtb
        module procedure sethtb4,sethtb8
      end interface

      contains
        integer*4 function sethtb8(token,type,ival)
        use iso_c_binding
        use maccode
        implicit none
        character*(*) token
        integer*8, target, intent(in):: ival
        integer*4 type
        integer*4 idx,hsrch,lenw
        sethtb8=0
     
        idx= hsrch(token(:lenw(token)))
        if(idx .le. 0 .or. idx .gt. HTMAX) then
          call errmsg('sethtb8'
     &         ,'illegal index value for sethashtble'
     &         , 0,16)
        else
          idtype(idx)=type
          if(type .eq. icRSVD) then
            idval(idx)=transfer(c_loc(ival),int8(0))/8
          else
            idval(idx)=ival
          endif
          sethtb8=idx
        endif
        return
        end function

        integer*4 function sethtb4(token,type,ival)
        implicit none
        character*(*) token
        integer*4 type,ival
        sethtb4=sethtb8(token,type,int8(ival))
        return
        end function

        subroutine forcesf
        implicit none
        call abort
        end subroutine

      end module

c     Don't confuse, Emacs. This is -*- fortran -*- mode!
      module tfcbk
      integer*4 maxgeneration,maxlevele,nsymhash,nslots
      parameter (maxgeneration=2**30-1,maxlevele=2**14,nsymhash=2047,
     $     nslots=32)
      integer*4 maxlbuf
      parameter (maxlbuf=2**22)
      real*8 dinfinity,dnotanumber
      integer*8
     $     itfcontroot,itfcontext,itfcontextpath,itflocal,
     $     kxeof,kxfailed,iaxhold,iaximmediate,
     $     iaxline,kxliteral,kxnull,kxnulll,kxnulls,
     $     iaxout,iaxpriority,iaxschar,
     $     iaxslotnull,iaxslotpart,iaxslotseqnull,
     $     kxvect,kxvect1,iavpw,
     $     kerror,ierrorf,ierrorgen,ierrorprint,ierrorth,
     $     ierrorexp,ifunbase,initmessage,levelcompile
      integer*4 
     $     levele,levelp,lgeneration,ltrace,
     $     modethrow,iordless
      end module

      module tfmem
      implicit none
      integer*8, parameter :: mpsize=2**22,kcpklist0=0
      integer*4, parameter :: nindex=64,mhash=32767,
     $     minseg0=9,minseg1=16,minseg2=16
      integer*8 :: icp=0,nmem=0,nnet=0,ich=0,maxic=3,
     $     minic,icsep,nitaloc
      integer*8 kfirstalloc

      type cbkalloc
        integer*8, allocatable :: ca(:)
      end type

      integer*8, parameter :: ncbk = 2**16, kcpoffset = 0*2**28,
     $     kcpthre  = 2**34, nlarge = 2**30, nitaloc0 = 2**24
      integer*4 icbk,jcbk
      type (cbkalloc), target :: sadalloc(ncbk)
      integer*8 kcbk(3,ncbk)

      interface sad_loc
        module procedure ksad_loc,isad_loc,rsad_loc
      end interface

      contains
        integer*8 function ksad_loc(k)
        use iso_c_binding
        implicit none
        integer*8, target :: k
        ksad_loc=(transfer(c_loc(k),int8(0))-kcpklist0)/8
        return
        end function

        integer*8 function isad_loc(i)
        use iso_c_binding
        implicit none
        integer*4, target:: i
        integer*8 k
        isad_loc=(transfer(c_loc(i),k)-kcpklist0)/8
        return
        end function

        integer*8 function rsad_loc(x)
        use iso_c_binding
        implicit none
        real*8, target:: x
        rsad_loc=(transfer(c_loc(x),int8(0))-kcpklist0)/8
        return
        end function

        subroutine tfcbkinit
        use iso_c_binding
        use maccbk
        implicit none
        type (c_ptr) cp
        integer*4, save::lps=0
        integer*4 getpagesize
c        kcpklist0=transfer(c_loc(kdummy),int8(0))
c        write(*,*)'tfcbkinit ',kcpklist0,2**31
c        kcpklist0=transfer(c_loc(kdummy),int8(0))-kcpoffset-8
        if(lps .eq. 0)then
          lps=getpagesize()
        endif
        call c_f_pointer(transfer(kcpklist0+8,cp),klist,[klistlen])
        call lminit(klist(0),lps)
        call c_f_pointer(c_loc(klist(1)),rlist,[klistlen])
        call c_f_pointer(c_loc(klist(1)),ilist,[2,klistlen])
        call c_f_pointer(c_loc(klist(1)),jlist,[8,klistlen])
        return
        end subroutine

        subroutine talocinit
        use maccbk
        use iso_c_binding
        implicit none
        integer*8 ka,ic
        allocate(sadalloc(1)%ca(nindex*2+mhash+16))
        ka=transfer(c_loc(sadalloc(1)%ca(1)),int8(0))
        kfirstalloc=ka
c     kcpklist0=0
        call tfcbkinit
        icp=ksad_loc(sadalloc(1)%ca(1))
        icbk=1
        jcbk=1
        kcbk(1,1)=icp
        kcbk(2,1)=icp+nindex*2+mhash+15
        kcbk(3,1)=kcbk(2,1)
        do ic=icp,icp+nindex*2,2
          klist(ic)=ic
          klist(ic+1)=ic
        enddo
        ich=icp+nindex*2+4
        do ic=ich,ich+mhash
          klist(ic)=ic
        enddo
        icsep=ich+mhash+1
        klist(icsep)=icsep
        nmem=0
        nnet=0
        return
        end subroutine

        integer*8 function ktaloc(n)
        use maccbk
        implicit none
        integer*4 n,m,n1,m1
        integer*8 ic1,i,ic,ic2,ip1,j,i1
        n1=max(n,3)
        m=n1+1
        if(n1 .lt. nindex)then
          ic=icp+n1*2
          i=klist(ic)
          if(i .ne. ic)then
            m=ilist(1,i-1)
            klist(ic)=klist(i)
            klist(klist(i)+1)=ic
            j=ich+iand(i+m+2,mhash)
            do while(klist(j) .ne. i+2)
              j=klist(j)
            enddo
            klist(j)=klist(i+2)
            nnet=nnet+m
            ktaloc=i
            return
          endif
          ic1=ic+min(m,minseg1)*2
          if(ic1 .le. maxic)then
            ic2=min(maxic-2,ic1+10)
            do ic=ic1,ic2,2
              i=klist(ic)
              if(ic .ne. i)then
                m1=ilist(1,i-1)
                klist(ic)=klist(i)
                klist(klist(i)+1)=ic
                j=ich+iand(i+m1+2,mhash)
                do while(klist(j) .ne. i+2)
                  j=klist(j)
                enddo
                klist(j)=klist(i+2)
                call tsetindexhash(i+m,m1-m)
                ilist(1,i-1)=m
                nnet=nnet+m
                ktaloc=i
                return
              endif
            enddo
            do ic=maxic,ic2+2,-2
              i=klist(ic)
              if(ic .ne. i)then
                m1=ilist(1,i-1)
                klist(ic)=klist(i)
                klist(klist(i)+1)=ic
                j=ich+iand(i+m1+2,mhash)
                do while(klist(j) .ne. i+2)
                  j=klist(j)
                enddo
                klist(j)=klist(i+2)
                call tsetindexhash(i+m,m1-m)
                ilist(1,i-1)=m
                nnet=nnet+m
                if(klist(ic) .eq. ic)then
                  maxic=ic-2
                else
                  maxic=ic
                endif
                ktaloc=i
                return
              endif
            enddo
            maxic=ic1-2
          endif
        endif
        ic=icp+nindex*2
 1000   i1=ic
        i=klist(i1)
        do while(i .ne. ic)
          m1=ilist(1,i-1)
          if(m1 .eq. m)then
            klist(i1)=klist(i)
            klist(klist(i)+1)=i1
            j=ich+iand(i+m+2,mhash)
            do while(klist(j) .ne. i+2)
              j=klist(j)
            enddo
            klist(j)=klist(i+2)
            nnet=nnet+m
            ktaloc=i
            return
          elseif(m1 .ge. m+minseg2)then
            klist(i1)=klist(i)
            klist(klist(i)+1)=i1
            j=ich+iand(i+m1+2,mhash)
            do while(klist(j) .ne. i+2)
              j=klist(j)
            enddo
            klist(j)=klist(i+2)
            call tsetindexhash(i+m,m1-m)
            ilist(1,i-1)=m
            nnet=nnet+m
            ktaloc=i
            return
          endif
          i1=i
          i=klist(i)
        enddo
        call talocp(m,ip1)
        if(ip1 .gt. 0)then
          go to 1000
        endif
        ktaloc=-1
        return
        end

        subroutine tsetindexhash(ip,m)
        use maccbk
        implicit none
        integer*4 m
        integer*8 ip,ic,ic1,ia
        if(m .gt. nindex)then
          ic=icp+nindex*2
        else
          ic=icp+(m-1)*2
          maxic=max(ic,maxic)
        endif
        klist(ip  )=klist(ic)
        klist(ip+1)=ic
        klist(klist(ic)+1)=ip
        klist(ic)=ip
        ia=ip+2
        ic1=ich+iand(ia+m,mhash)
        klist(ia  )=klist(ic1)
        klist(ic1 )=ia
        ilist(1,ip-1)=m
c     call tfsetlastp(ip+m-1)
        return
        end

        subroutine tfree(ka)
        use maccbk
        implicit none
        integer*8 ka,ix,ik,ik0,ip,ix1
        integer*4 m,mx
        m=ilist(1,ka-1)
        if(m .lt. 4)then
          if(m .ne. 0)then
            write(*,*)'tfree-too small segment: ',ka,m
            call abort
          endif
          return
        endif
        nnet=nnet-m
        ix=ka+2
        ik=iand(ix,mhash)+ich
        ik0=ik
        ip=klist(ik)
        do while(ip .ne. ik0)
          if(ip .lt. ix)then
            mx=ilist(1,ip-3)
            if(ip+mx .eq. ix)then
              klist(ik  )=klist(ip)
              klist(klist(ip-2)+1)=klist(ip-1)
              klist(klist(ip-1)  )=klist(ip-2)
              m=m+mx
              ix=ip
              exit
            endif
          endif
          ik=ip
          ip=klist(ik)
        enddo
        ix1=ix+m
c     if(tfchecklastp(ix1))then
        ik=iand(ix1+ilist(1,ix1-3),mhash)+ich
        ik0=ik
        ip=klist(ik)
        do while(ip .ne. ik0)
          if(ip .eq. ix1)then
            klist(ik  )=klist(ip)
            klist(klist(ip-2)+1)=klist(ip-1)
            klist(klist(ip-1)  )=klist(ip-2)
            m=m+ilist(1,ix1-3)
            exit
          endif
          ik=ip
          ip=klist(ik)
        enddo
c     endif
        call tsetindexhash(ix-2,m)
        return
        end

      end module

      module tfcode
      real*8 xinfinity
      parameter (xinfinity=1.7976931348623157D308)
      integer*4 ntfoper,ntfreal,ntflist,ntflistr,ntfdef,ntfstkseq,
     $     ntfstring,ntfsymbol,ntfpat,ntffun,ntfstk,ntfarg
      parameter (ntfoper=0,ntfreal=1,ntflist=3,ntflistr=4,ntfstkseq=5,
     $     ntfstk=6,
     $     ntfstring=101,ntfsymbol=201,ntfpat=203,ntfarg=205,
     $     ntffun=ntfoper,ntfdef=ntfsymbol)
      integer*4 nfunif,nfunlength,nfundo,nfunmodule,nfunblock,
     $     nfununeval,nfunwith,nfunthread
      parameter (nfunif=39,nfunlength=20,nfundo=29,nfunmodule=34,
     $     nfunblock=35,nfununeval=132,nfunwith=141,nfunthread=83)
      integer*4 mtfnull,mtfneg,mtfinv,mtfplus,mtfminus,mtfmult,mtfdiv,
     $     mtfdot,mtfpower,mtfequal,mtfunequal,mtfgreater,mtfless,
     $     mtfgeq,mtfleq,mtfsame,mtfunsame,mtfnot,mtfand,mtfor,
     $     mtfconcat,mtfleftbra,mtfrightbra,mtfleftbrace,
     $     mtfrightbrace,mtfsetdelayed,mtfset,mtfcomplex,mtfleftparen,
     $     mtfrightparen,mtfcomma,mtftimes,
     $     mtfcomp,mtffun,mtfcolon,mtfrule,mtfruledelayed,mtfreplace,
     $     mtfreplacerepeated,mtfupset,mtfupsetdelayed,mtfunset,
     $     mtfpattest,mtfflag,mtfslot,mtfslotseq,mtfrevpower,mtfalt,
     $     mtflist,mtfmap,mtfmapall,mtfapply,mtfrepeated,
     $     mtfrepeatednull,mtfinequality,mtfaddto,mtfsubtractfrom,
     $     mtftimesby,mtfdivideby,mtfincrement,mtfdecrement,
     $     mtfpart,mtfatt,mtfmessagename,mtftagset,
     $     mtfleftcomment,mtfrightcomment,mtfhold,mtfend
      parameter (
     $     mtfnull=0,
     $     mtfneg=1,mtfinv=2,mtfplus=3,mtfminus=4,mtfmult=5,mtfdiv=6,
     $     mtfrevpower=7,mtfpower=8,mtfequal=9,mtfunequal=10,
     $     mtfgreater=11,
     $     mtfless=12,mtfgeq=13,mtfleq=14,mtfsame=15,mtfunsame=16,
     $     mtfnot=17,mtfand=18,mtfor=19,mtfconcat=20,mtfleftbra=21,
     $     mtfrightbra=22,mtfleftbrace=23,mtfrightbrace=24,
     $     mtfsetdelayed=25,mtfset=26,mtfcomplex=27,mtfleftparen=28,
     $     mtfrightparen=29,mtfcomma=30,mtfcomp=31,mtffun=32,
     $     mtfcolon=33,mtfrule=34,mtfruledelayed=35,mtfreplace=36,
     $     mtfreplacerepeated=37,mtfupset=38,mtfupsetdelayed=39,
     $     mtfunset=40,mtfpattest=41,mtfflag=42,mtfslot=43,
     $     mtfslotseq=44,mtfdot=45,mtfalt=46,mtfmap=47,mtfmapall=48,
     $     mtfapply=49,mtfrepeated=50,mtfrepeatednull=51,
     $     mtfinequality=52,mtfaddto=53,mtfsubtractfrom=54,
     $     mtftimesby=55,mtfdivideby=56,
     $     mtfincrement=57,mtfdecrement=58,
     $     mtfpart=59,mtfatt=60,mtfmessagename=61,mtftagset=62,
     $     mtfleftcomment=63,mtfrightcomment=64,mtfhold=65,
     $     mtfend=66,
     $     mtflist=mtfleftbrace,mtftimes=mtfmult)
      integer*4 mtfnopc
      parameter (mtfnopc=mtfend)
      integer*4 lsimplepat,lsimplepatlist,lconstlist,lnoconstlist,
     $     larglist,lnopatlist,lmemberlist,lnodefsymbol,lnoseqlist,
     $     lnonreallist
      parameter (lsimplepat=1,lsimplepatlist=2,larglist=4,
     $     lconstlist=8,lnoconstlist=16,lnopatlist=32,lmemberlist=64,
     $     lnodefsymbol=128,lnoseqlist=256,lnonreallist=512)
      integer*4 kconstlist,knoconstlist,kfixedarg,knofixedarg,
     $     kallnofixedarg,knopatarg,kpatarg,kconstarg,knoconstarg,
     $     kseqarg,knoseqarg,ksymbollist,knosymbollist,ktoberebuilt
      parameter (knopatarg=lnonreallist*2,
     $     kpatarg=lnonreallist*4,
     $     kconstlist=lnonreallist*8,
     $     knoconstlist=lnonreallist*16,
     $     kfixedarg=lnonreallist*32,
     $     knofixedarg=lnonreallist*64,
     $     kallnofixedarg=lnonreallist*128,
     $     kconstarg=lnonreallist*256,
     $     knoconstarg=lnonreallist*512,
     $     kseqarg=lnonreallist*1024,
     $     knoseqarg=lnonreallist*2048,
     $     ksymbollist=lnonreallist*4096,
     $     knosymbollist=lnonreallist*8192,
     $     ktoberebuilt=lnonreallist*16384)
      integer*4 iattrholdfirst,iattrholdrest,iattrholdall,
     $     iattrconstant,iattrimmediate,iattrorderless,
     $     iattrdynamic,iattrprotected,iattrnumeric
      parameter (iattrholdfirst=1,iattrholdrest=2,
     $     iattrholdall=iattrholdfirst+iattrholdrest,
     $     iattrimmediate=4,iattrnumeric=8,
     $     iattrorderless=16,iattrdynamic=32,
     $     iattrprotected=64,iattrconstant=128)
      real*8 rtfnull
      parameter (rtfnull=0.d0)
      integer*8 ktfoper,ktflist,ktfstring,ktfsymbol,ktfpat,ktfobj,
     $     ktfmask,ktamask,ktrmask,ktfnull,ktfnr,ktfref,ktfother,
     $     ktomask,ktftrue,ktfnan
      parameter (
     $     ktfnull  =int8(z'fff0000000000000'),
     $     ktfother =int8(z'fff2000000000000'),
     $     ktfnr    =int8(z'7ff2000000000000'),
     $     ktfoper  =int8(z'fff6000000000000'),
     $     ktfref   =int8(z'fffa000000000000'),
     $     ktfobj   =int8(z'7ff2000000000000'),
     $     ktflist  =int8(z'7ff2000000000000'),
     $     ktfpat   =int8(z'7ff6000000000000'),
     $     ktfstring=int8(z'7ffa000000000000'),
     $     ktfsymbol=int8(z'7ffe000000000000'),
     $     ktomask  =int8(z'fff2000000000000'),
     $     ktrmask  =int8(z'7ff2000000000000'),
     $     ktfmask  =int8(z'fffe000000000000'),
     $     ktamask  =int8(z'0001ffffffffffff'),
     $     ktftrue  =int8(z'3ff0000000000000'),
     $     ktfnan   =int8(z'fff8000000000000'))

      type sad_descriptor
      sequence
      integer*8 k
      end type

      type sad_object
      sequence
      integer*4 len,attr
      integer*8 alloc
      integer*4 ref,nl
      integer*8 body(0:2**31-1)
      end type

      type sad_list
      sequence
      integer*2 lenp,lena
      integer*4 attr
      integer*8 alloc
      integer*4 ref,nl
      integer*8 head
      real*8 rbody(1:0)
      complex*16 cbody(1:0)
      type (sad_descriptor) dbody(1:0)
      integer*8 body(1:2**31-1)
      end type

      type sad_complex
      sequence
      integer*2 lenp,lena
      integer*4 attr
      integer*8 alloc
      integer*4 ref,nl
      integer*8 head
      complex*16 cx(1:0)
      integer*8 body(1:0)
      type (sad_descriptor) dbody(1:0)
      real*8 re,im
      end type

      type sad_symbol
      sequence
      integer*4 attr,override
      type (sad_descriptor) alloc
      integer*4 ref,gen
      integer*8 loc,dummy(11:10)
      end type

      type sad_symdef
      sequence
      integer*4 len,attr
      integer*8 next,prev,upval,downval
      type (sad_descriptor) value
      type (sad_symbol) sym
      end type

      type sad_funtbl
      sequence
      type (sad_symdef) def
      integer*4 narg,id,mapeval(2,1)
      end type

      type sad_pat
      sequence
      integer*4 len,mat
      integer*8 alloc
      integer*4 ref,gen
      type (sad_descriptor) expr,head
      type (sad_pat), pointer :: equiv
      type (sad_descriptor) value
      type (sad_symbol) sym
      type (sad_descriptor) default
      end type

      type sad_string
      sequence
      integer*4 len,override
      integer*8 alloc
      integer*4 ref,gen
      integer*4 nch,nc
      integer*1 istr(1:0)
      integer*8 kstr(1:0)
      character*(2**31-1) str
      end type

      type sad_namtbl
      sequence
      integer*4 len,dummy
      integer*8 next,symdef,cont
      type (sad_string) str
      end type

      type sad_deftbl
      sequence
      integer*4 len,pat
      integer*8 next,prev
      integer*4 npat,attr
      type (sad_descriptor) arg,argc,body,bodyc
      real*8 compile
      type (sad_descriptor) pattbl(1:2**10)
      end type

      type sad_defhash
      sequence
      integer*4 len,attr
      integer*8 next,prev
      integer*4 gen,nhash
      integer*8 hash(0:2**10-1)
      end type

      end module

      module tfstk
      use tfcbk
      use tfcode
      use maccbk
      use tfmem, only:sad_loc,ksad_loc,ktaloc,tfree
      integer*8 ispbase
      integer*4 mstk,isp,ivstkoffset,ipurefp,napuref,isporg
      integer*4, pointer, dimension(:,:) :: ivstk,itastk,ivstk2,itastk2
      real*8, pointer, dimension(:) :: vstk,rtastk,vstk2,rtastk2
      integer*8, pointer, dimension(:) :: ktastk,ktastk2
      type (sad_descriptor), pointer, dimension(:) :: dtastk,dtastk2,
     $     dlist
      type (sad_descriptor) kxmatrix,dxliteral,dxeof,dxfailed,
     $     dxvect,dxvect1,dxnull,dxnulll,dxnulls

c      integer*4 ivstk (2,RBASE:RBASE+MAXMEM0-1)
c      real*8    vstk  (  RBASE:RBASE+MAXMEM0-1)
c      integer*4 itastk(2,RBASE:RBASE+MAXMEM0-1)
c      integer*2 jtastk(4,RBASE:RBASE+MAXMEM0-1)
c      integer*8 ktastk(  RBASE:RBASE+MAXMEM0-1)
c      real*8    rtastk(  RBASE:RBASE+MAXMEM0-1)
c      equivalence ( ivstk(1,RBASE),ilist(1,RBASE))
c      equivalence (  vstk(  RBASE),ilist(1,RBASE))
c      equivalence (itastk(1,RBASE),ilist(1,RBASE))
c      equivalence (jtastk(1,RBASE),ilist(1,RBASE))
c      equivalence (rtastk(  RBASE),ilist(1,RBASE))
c      equivalence (ktastk(  RBASE),ilist(1,RBASE))
      logical*4 :: tfstkinit = .false.

      type (sad_symdef), pointer :: redmath

      interface loc_sad
        module procedure loc_list,loc_sym,loc_string,
     $     loc_pat,loc_obj,loc_complex,loc_symdef
      end interface

      interface descr_sad
        module procedure descr_list,descr_sym,descr_string,
     $     descr_pat,descr_obj,descr_complex,descr_symdef
      end interface

      interface ktfaddr
        module procedure ktfaddrk,ktfaddrd
      end interface

      contains
        subroutine tfinitstk
        use iso_c_binding
        implicit none
        integer*4 igetgl1
        if(tfstkinit)then
          return
        endif
        mstk=max(2**18,igetgl1('STACKSIZ'))
        ispbase=ktaloc(mstk*2)-1
        if(ispbase .le. 0)then
          write(*,*)'Stack allocation failed: ',mstk,ispbase
          call abort
        endif
        isp=0
        isporg=isp+1
        ivstkoffset=mstk
        ipurefp=0
        napuref=0
        call c_f_pointer(c_loc(klist(ispbase+1)),vstk,[klistlen])
        call c_f_pointer(c_loc(klist(ispbase+1)),rtastk,[klistlen])
        call c_f_pointer(c_loc(klist(ispbase+1)),ktastk,[klistlen])
        call c_f_pointer(c_loc(klist(ispbase+1)),dtastk,[klistlen])
        call c_f_pointer(c_loc(klist(ispbase+1)),itastk,[2,klistlen])
        call c_f_pointer(c_loc(klist(ispbase+1)),ivstk,[2,klistlen])
        call c_f_pointer(c_loc(klist(ivstkoffset+ispbase+1)),
     $       vstk2,[klistlen])
        call c_f_pointer(c_loc(klist(ivstkoffset+ispbase+1)),
     $       rtastk2,[klistlen])
        call c_f_pointer(c_loc(klist(ivstkoffset+ispbase+1)),
     $       ktastk2,[klistlen])
        call c_f_pointer(c_loc(klist(ivstkoffset+ispbase+1)),
     $       dtastk2,[klistlen])
        call c_f_pointer(c_loc(klist(ivstkoffset+ispbase+1)),
     $       itastk2,[2,klistlen])
        call c_f_pointer(c_loc(klist(ivstkoffset+ispbase+1)),
     $       ivstk2,[2,klistlen])
        call c_f_pointer(c_loc(klist(1)),dlist,[klistlen])
        tfstkinit=.true.
        return
        end subroutine

        integer*8 function ktfsadalloc(n)
        use tfmem
        use iso_c_binding
        implicit none
        integer*8 n
        integer*4 istat,i
        ktfsadalloc=0
        do i=icbk+1,int(ncbk)
          allocate(sadalloc(i)%ca(n),stat=istat)
          if(istat .ne. 0)then
            write(*,*)'ktfsadalloc allocation error in ALLOCATE: ',
     $           istat,n
            call abort
          endif
          icbk=icbk+1
          ktfsadalloc=sad_loc(sadalloc(i)%ca(1))
          if(ktfsadalloc .ge. 0)then
            call tfentercbk(ktfsadalloc,n)
            return
          else
            write(*,*)'Negative allocation - retry: ',i,ktfsadalloc
          endif
        enddo          
        if(icbk .ge. ncbk)then
          write(*,*)'ktfsadalloc too many allocations: ',icbk
          call abort
        endif
        end function

        subroutine tfentercbk(ka,n)
        use tfmem
        implicit none
        integer*8 ka,n
        integer*4 j,j0,k
        j0=jcbk+1
        do j=1,jcbk
          if(kcbk(2,j) .eq. ka-1)then
            kcbk(2,j)=ka+n-1
            kcbk(3,j)=kcbk(2,j)
            do k=1,jcbk
              if(kcbk(1,k) .eq. kcbk(2,j)+1)then
                if(k .lt. j)then
                  kcbk(1,k)=kcbk(1,j)
                  kcbk(1,j)=0
                  kcbk(2,j)=0
                  kcbk(3,j)=0
                  if(j .eq. jcbk)then
                    jcbk=jcbk-1
                  endif
                else
                  kcbk(2,j)=kcbk(2,k)
                  kcbk(3,j)=kcbk(2,k)
                  kcbk(1,k)=0
                  kcbk(2,k)=0
                  kcbk(3,k)=0
                endif
                return
              endif
            enddo
          elseif(kcbk(1,j) .eq. ka+n)then
            kcbk(1,j)=ka
            do k=1,jcbk
              if(kcbk(2,k) .eq. kcbk(1,j)+1)then
                if(k .lt. j)then
                  kcbk(2,k)=kcbk(2,j)
                  kcbk(3,k)=kcbk(3,j)
                  kcbk(1,j)=0
                  kcbk(2,j)=0
                  kcbk(3,j)=0
                  if(j .eq. jcbk)then
                    jcbk=jcbk-1
                  endif
                else
                  kcbk(1,j)=kcbk(1,k)
                  kcbk(1,k)=0
                  kcbk(2,k)=0
                  kcbk(3,k)=0
                endif
                return
              endif
            enddo
          elseif(kcbk(2,j) .eq. 0)then
            j0=min(j0,j)
          endif
        enddo
        kcbk(1,j0)=ka
        kcbk(2,j0)=ka+n-1
        kcbk(3,j0)=kcbk(2,j0)
        jcbk=max(jcbk,j0)
        return
        end subroutine

        integer*4 function itfcbk(k)
        use tfmem
        implicit none
        integer*8 k
        integer*4 i
        do i=1,jcbk
          if(k .le. kcbk(2,i) .and. k .ge. kcbk(1,i))then
            itfcbk=i
            return
          endif
        enddo
        itfcbk=0
        return
        end

        subroutine tfsetlastp(k)
        use tfmem
        implicit none
        integer*8 k
        integer*4 i
        i=itfcbk(k)
        if(i .ne. 0)then
          kcbk(3,i)=max(kcbk(3,i),k)
        endif
        return
        end

        logical*4 function tfchecklastp(k)
        use tfmem
        implicit none
        integer*8 k
        integer*4 i
        i=itfcbk(k)
        if(i .ne. 0)then
          tfchecklastp=k .le. kcbk(3,i)
        else
          tfchecklastp=.false.
        endif
        if(.not. tfchecklastp)then
          if(i .ne. 0)then
            write(*,*)'tfcklastp ',k,i,kcbk(3,i),kcbk(2,i)
          else
            write(*,*)'tfcklastp ',k,0
          endif
        endif
        return
        end

        integer*4 function iget_fun_id(ka)
        use iso_c_binding
        implicit none
        type (sad_funtbl), pointer :: fun
        integer*8 ka
        call c_f_pointer(c_loc(klist(klist(ifunbase+ka)-9)),fun)
        iget_fun_id=fun%id
        return
        end function

        subroutine loc_list(locp,list)
        use iso_c_binding
        implicit none
        type (sad_list), pointer, intent(out) :: list
        integer*8 locp
        call c_f_pointer(c_loc(klist(locp-3)),list)
        return
        end subroutine

        subroutine descr_list(dscr,list)
        implicit none
        type (sad_descriptor) dscr
        type (sad_list), pointer, intent(out) :: list
        call loc_list(ktfaddrd(dscr),list)
        return
        end subroutine

        subroutine loc_complex(locp,cx)
        use iso_c_binding
        implicit none
        type (sad_complex), pointer, intent(out) :: cx
        integer*8 locp
        call c_f_pointer(c_loc(klist(locp-3)),cx)
        return
        end subroutine

        subroutine descr_complex(dscr,complex)
        implicit none
        type (sad_descriptor) dscr
        type (sad_complex), pointer, intent(out) :: complex
        call loc_complex(ktfaddrd(dscr),complex)
        return
        end subroutine

        subroutine loc_obj(locp,obj)
        use iso_c_binding
        implicit none
        type (sad_object), pointer, intent(out) :: obj
        integer*8 locp
        call c_f_pointer(c_loc(klist(locp-3)),obj)
        return
        end subroutine

        subroutine descr_obj(dscr,obj)
        implicit none
        type (sad_descriptor) dscr
        type (sad_object), pointer, intent(out) :: obj
        call loc_obj(ktfaddrd(dscr),obj)
        return
        end subroutine

        subroutine loc_namtbl(locp,loc)
        use iso_c_binding
        implicit none
        type (sad_namtbl), pointer, intent(out) :: loc
        integer*8 locp
        call c_f_pointer(c_loc(klist(locp-1)),loc)
        return
        end subroutine

        subroutine descr_namtbl(dscr,namtbl)
        implicit none
        type (sad_descriptor) dscr
        type (sad_namtbl), pointer, intent(out) :: namtbl
        call loc_namtbl(ktfaddrd(dscr),namtbl)
        return
        end subroutine

        subroutine sym_namtbl(sym,loc)
        implicit none
        type (sad_symbol) sym
        type (sad_namtbl), pointer, intent(out) :: loc
        call loc_namtbl(sym%loc,loc)
        return
        end subroutine

        subroutine sym_symstr(sym,str)
        implicit none
        type (sad_symbol) sym
        type (sad_string), pointer, intent(out) :: str
        call loc_symstr(sym%loc,str)
        return
        end subroutine

        subroutine loc_symstr(loc,str)
        implicit none
        integer*8 loc
        type (sad_string), pointer, intent(out) :: str
        type (sad_namtbl), pointer :: nam
        if(loc .eq. 0)then
          call loc_sad(ktfaddr(kxnulls),str)
        else
          call loc_namtbl(loc,nam)
          str=>nam%str
        endif
        return
        end subroutine

        subroutine loc_pat(locp,pat)
        use iso_c_binding
        implicit none
        type (sad_pat), pointer, intent(out) :: pat
        integer*8 locp
        call c_f_pointer(c_loc(klist(locp-3)),pat)
        return
        end subroutine

        subroutine descr_pat(dscr,pat)
        implicit none
        type (sad_descriptor) dscr
        type (sad_pat), pointer, intent(out) :: pat
        call loc_pat(ktfaddrd(dscr),pat)
        return
        end subroutine

        subroutine loc_string(locp,str)
        use iso_c_binding
        implicit none
        type (sad_string), pointer, intent(out) :: str
        integer*8 locp
        call c_f_pointer(c_loc(klist(locp-3)),str)
        return
        end subroutine

        subroutine descr_string(dscr,string)
        implicit none
        type (sad_descriptor) dscr
        type (sad_string), pointer, intent(out) :: string
        call loc_string(ktfaddrd(dscr),string)
        return
        end subroutine

        subroutine loc_sym(locp,sym)
        use iso_c_binding
        implicit none
        type (sad_symbol), pointer, intent(out) :: sym
        integer*8 locp
        call c_f_pointer(c_loc(klist(locp-3)),sym)
        return
        end subroutine

        subroutine descr_sym(dscr,symbol)
        implicit none
        type (sad_descriptor) dscr
        type (sad_symbol), pointer, intent(out) :: symbol
        call loc_sym(ktfaddrd(dscr),symbol)
        return
        end subroutine

        subroutine loc1_symdef(locp1,symd)
        use iso_c_binding
        implicit none
        type (sad_symdef), pointer, intent(out) :: symd
        integer*8 locp1
        call c_f_pointer(c_loc(klist(locp1-1)),symd)
        return
        end subroutine

        subroutine loc_symdef(locp,symd)
        use iso_c_binding
        implicit none
        type (sad_symdef), pointer, intent(out) :: symd
        integer*8 locp
        call c_f_pointer(c_loc(klist(locp-9)),symd)
        return
        end subroutine

        subroutine descr_symdef(dscr,symdef)
        implicit none
        type (sad_descriptor) dscr
        type (sad_symdef), pointer, intent(out) :: symdef
        call loc_symdef(ktfaddrd(dscr),symdef)
        return
        end subroutine

        subroutine sym_symdef(sym,symd)
        use iso_c_binding
        implicit none
        type (sad_symbol), target :: sym
        type (sad_symdef), pointer, intent(out) :: symd
        call c_f_pointer(c_loc(sym%dummy(1)),symd)
        return
        end subroutine

        subroutine loc_deftbl(locp,loc)
        use iso_c_binding
        implicit none
        type (sad_deftbl), pointer, intent(out) :: loc
        integer*8 locp
        call c_f_pointer(c_loc(klist(locp-1)),loc)
        return
        end subroutine

        subroutine loc_defhash(locp,loc)
        use iso_c_binding
        implicit none
        type (sad_defhash), pointer, intent(out) :: loc
        integer*8 locp
        call c_f_pointer(c_loc(klist(locp-1)),loc)
        return
        end subroutine

        subroutine list_list(kl,kl1)
        use iso_c_binding
        implicit none
        type (sad_list) , target :: kl
        type (sad_list), pointer :: kl1
        kl1=>kl
        return
        end subroutine list_list

        integer*8 function ktfaddrk(k)
        implicit none
        integer*8 k
        ktfaddrk=iand(ktamask,k)
        return
        end function ktfaddrk

        integer*8 function ktfaddrd(k)
        implicit none
        type (sad_descriptor) k
        ktfaddrd=iand(ktamask,k%k)
        return
        end function ktfaddrd

        integer*8 function ktftype(k)
        implicit none
        integer*8 k
        ktftype=iand(ktfmask,k)
        return
        end function ktftype

        logical*4 function ktfobjq(k,obj)
        implicit none
        type (sad_object), pointer, optional, intent(out) :: obj
        integer*8 k
        if(iand(ktomask,k) .eq. ktfobj)then
          ktfobjq=.true.
          if(present(obj))then
            call loc_obj(iand(ktamask,k),obj)
          endif
        else
          ktfobjq=.false.
        endif
        return
        end function ktfobjq

        logical*4 function ktfobjqd(k,obj)
        implicit none
        type (sad_descriptor) k
        type (sad_object), pointer, optional, intent(out) :: obj
        if(iand(ktomask,k%k) .eq. ktfobj)then
          ktfobjqd=.true.
          if(present(obj))then
            call loc_obj(iand(ktamask,k%k),obj)
          endif
        else
          ktfobjqd=.false.
        endif
        return
        end function ktfobjqd

        logical*4 function ktfnonobjq(ka)
        implicit none
        integer*8 ka
        ktfnonobjq=iand(ktomask,ka) .ne. ktfobj
        return
        end function ktfnonobjq

        logical*4 function ktfrealq(k)
        implicit none
        integer*8 k
        ktfrealq=iand(ktrmask,k) .ne. ktfnr
        return
        end function ktfrealq

        logical*4 function ktfnonrealq(k)
        implicit none
        integer*8 k
        ktfnonrealq=iand(ktrmask,k) .eq. ktfnr
        return
        end function ktfnonrealq

        logical*4 function ktfrealqd(k,v)
        implicit none
        type (sad_descriptor) k
        real*8, optional, intent(out) :: v
        real*8 rfromk
        ktfrealqd=iand(ktrmask,k%k) .ne. ktfnr
        if(ktfrealqd .and. present(v))then
          v=rfromk(k%k)
        endif
        return
        end function ktfrealqd

        logical*4 function ktfnonrealqd(k,v)
        implicit none
        type (sad_descriptor) k
        real*8 , optional, intent(out) :: v
        real*8 rfromk
        ktfnonrealqd=iand(ktrmask,k%k) .eq. ktfnr
        if(.not. ktfnonrealqd .and. present(v))then
          v=rfromk(k%k)
        endif
        return
        end function ktfnonrealqd

        logical*4 function ktfrealqdi(k,iv)
        implicit none
        type (sad_descriptor) k
        integer*4, optional, intent(out) :: iv
        real*8 rfromk
        ktfrealqdi=iand(ktrmask,k%k) .ne. ktfnr
        if(ktfrealqdi .and. present(iv))then
          iv=rfromk(k%k)
        endif
        return
        end function ktfrealqdi

        logical*4 function ktfnonrealqdi(k,iv)
        implicit none
        type (sad_descriptor) k
        integer*4 , optional, intent(out) :: iv
        real*8 rfromk
        ktfnonrealqdi=iand(ktrmask,k%k) .eq. ktfnr
        if(.not. ktfnonrealqdi .and. present(iv))then
          iv=rfromk(k%k)
        endif
        return
        end function ktfnonrealqdi

        logical*4 function ktfoperq(k)
        implicit none
        integer*8 k
        ktfoperq=iand(ktfmask,k) .eq. ktfoper
        return
        end function ktfoperq

        logical*4 function ktfnonoperq(k)
        implicit none
        integer*8 k
        ktfnonoperq=iand(ktfmask,k) .ne. ktfoper
        return
        end function ktfnonoperq

        logical*4 function ktfoperqd(k,ka)
        implicit none
        type (sad_descriptor) k
        integer*8, optional, intent(out) :: ka
        ktfoperqd=iand(ktfmask,k%k) .eq. ktfoper
        if(ktfoperqd .and. present(ka))then
          ka=ktfaddr(k%k)
        endif
        return
        end function ktfoperqd

        logical*4 function ktfnonoperqd(k,ka)
        implicit none
        type (sad_descriptor) k
        integer*8, optional, intent(out) :: ka
        ktfnonoperqd=iand(ktfmask,k%k) .ne. ktfoper
        if(.not. ktfnonoperqd .and. present(ka))then
          ka=ktfaddr(k%k)
        endif
        return
        end function ktfnonoperqd

        logical*4 function ktfstringq(k,str)
        implicit none
        type (sad_string), pointer, optional, intent(out) :: str
        integer*8 k
        ktfstringq=iand(ktfmask,k) .eq. ktfstring
        if(present(str) .and. ktfstringq)then
          call loc_string(ktfaddr(k),str)
        endif
        return
        end function ktfstringq

        logical*4 function ktfstringqd(k,str)
        implicit none
        type (sad_descriptor) k
        type (sad_string), pointer, optional, intent(out) :: str
        ktfstringqd=iand(ktfmask,k%k) .eq. ktfstring
        if(present(str) .and. ktfstringqd)then
          call loc_string(ktfaddr(k%k),str)
        endif
        return
        end function ktfstringqd

        logical*4 function ktfnonstringq(k,str)
        implicit none
        type (sad_string), pointer, optional, intent(out) :: str
        integer*8 k
        ktfnonstringq=iand(ktfmask,k) .ne. ktfstring
        if(present(str) .and. .not. ktfnonstringq)then
          call loc_string(ktfaddr(k),str)
        endif
        return
        end function ktfnonstringq

        logical*4 function ktflistq(k,kl)
        implicit none
        type (sad_list), pointer, optional, intent(out) :: kl
        integer*8 k
        if(iand(ktfmask,k) .eq. ktflist)then
          ktflistq=.true.
          if(present(kl))then
            call loc_list(iand(ktamask,k),kl)
          endif
        else
          ktflistq=.false.
        endif          
        return
        end function ktflistq

        logical*4 function ktflistqd(k,kl)
        implicit none
        type (sad_descriptor) k
        type (sad_list), pointer, optional, intent(out) :: kl
        if(iand(ktfmask,k%k) .eq. ktflist)then
          ktflistqd=.true.
          if(present(kl))then
            call loc_list(iand(ktamask,k%k),kl)
          endif
        else
          ktflistqd=.false.
        endif          
        return
        end function ktflistqd

        logical*4 function ktfnonlistq(k,kl)
        implicit none
        type (sad_list), pointer, optional, intent(out) :: kl
        integer*8 k
        if(iand(ktfmask,k) .eq. ktflist)then
          ktfnonlistq=.false.
          if(present(kl))then
            call loc_list(iand(ktamask,k),kl)
          endif
        else
          ktfnonlistq=.true.
        endif          
        return
        end function ktfnonlistq

        logical*4 function ktfnonlistqd(k,kl)
        implicit none
        type (sad_descriptor) k
        type (sad_list), pointer, optional, intent(out) :: kl
        if(iand(ktfmask,k%k) .eq. ktflist)then
          ktfnonlistqd=.false.
          if(present(kl))then
            call loc_list(iand(ktamask,k%k),kl)
          endif
        else
          ktfnonlistqd=.true.
        endif          
        return
        end function ktfnonlistqd

        logical*4 function tfreallistqd(k,kl)
        implicit none
        type (sad_descriptor) k
        type (sad_list), pointer, optional, intent(out) :: kl
        type (sad_list), pointer :: kl1
        if(ktflistqd(k,kl1))then
          tfreallistqd=kl1%head .eq. ktfoper+mtflist
     $         .and. ktfreallistqo(kl1)
          if(tfreallistqd .and. present(kl))then
            kl=>kl1
          endif
        else
          tfreallistqd=.false.
        endif
        return
        end function

        logical*4 function tfreallistq(k,kl)
        implicit none
        type (sad_list), pointer, optional, intent(out) :: kl
        type (sad_list), pointer :: kl1
        integer*8 k
        if(ktflistq(k,kl1))then
          tfreallistq=kl1%head .eq. ktfoper+mtflist
     $         .and. ktfreallistqo(kl1)
          if(tfreallistq .and. present(kl))then
            kl=>kl1
          endif
        else
          tfreallistq=.false.
        endif
        return
        end function

        logical*4 function ktflistqx(k,cx)
        implicit none
        type (sad_complex), pointer, optional, intent(out) :: cx
        integer*8 k
        if(iand(ktfmask,k) .eq. ktflist)then
          ktflistqx=.true.
          if(present(cx))then
            call loc_complex(iand(ktamask,k),cx)
          endif
        else
          ktflistqx=.false.
        endif          
        return
        end function ktflistqx

        logical*4 function tflistqk(k,kl)
        implicit none
        type (sad_list), pointer, optional, intent(out) :: kl
        integer*8 k,ka
        if(iand(ktfmask,k) .eq. ktflist)then
          ka=ktfaddr(k)
          if(klist(ka) .eq. ktfoper+mtflist)then
            tflistqk=.true.
            if(present(kl))then
              call loc_list(ka,kl)
            endif
            return
          endif
        endif
        tflistqk=.false.
        return
        end function tflistqk

        logical*4 function tflistqd(k,kl)
        implicit none
        type (sad_descriptor) k
        type (sad_list), pointer, optional, intent(out) :: kl
        tflistqd=tflistqk(k%k,kl)
        return
        end function tflistqd

        logical*4 function tfnonlistqk(k,kl)
        implicit none
        type (sad_list), pointer, optional, intent(out) :: kl
        integer*8 k,ka
        if(iand(ktfmask,k) .eq. ktflist)then
          ka=ktfaddr(k)
          if(klist(ka) .eq. ktfoper+mtflist)then
            tfnonlistqk=.false.
            if(present(kl))then
              call loc_list(ka,kl)
            endif
            return
          endif
        endif
        tfnonlistqk=.true.
        return
        end function tfnonlistqk

        recursive logical*4 function tfruleqk(k,klx) result(lx)
        implicit none
        type (sad_list), pointer :: kl
        type (sad_list), pointer, optional, intent(out) :: klx
        integer*8 k
        integer*4 i
        lx=.false.
        if(ktflistq(k,kl))then
          select case (kl%head)
          case (ktfoper+mtflist)
            if(.not. ktfnonreallistqo(kl))return
            do i=1,kl%nl
              if(.not. tfruleqk(kl%body(i)))return
            enddo
          case (ktfoper+mtfrule,ktfoper+mtfruledelayed)
            if(kl%nl .ne. 2)return
          case default
            return
          end select
        else
          return
        endif
        lx=.true.
        if(present(klx))then
          klx=>kl
        endif
        return
        end function tfruleqk

        logical*4 function tfruleqd(k,klx)
        implicit none
        type (sad_descriptor) k
        type (sad_list), pointer, optional, intent(out) :: klx
        tfruleqd=tfruleqk(k%k,klx)
        return
        end function tfruleqd

        logical*4 function tfnumberqd(k,c)
        implicit none
        type (sad_descriptor) k
        type (sad_complex), pointer :: cx
        complex*16, optional, intent(out) :: c
        real*8 v
        if(ktfrealqd(k,v))then
          tfnumberqd=.true.
          if(present(c))then
            c=v
          endif
        else
          tfnumberqd=tfcomplexqx(k%k,cx)
          if(tfnumberqd .and. present(c))then
            c=cx%cx(1)
          endif
        endif
        return
        end function

        logical*4 function tfnumlistqnk(k,n,kl1)
        implicit none
        type (sad_descriptor) k
        type (sad_list), pointer, optional, intent(out) :: kl1
        type (sad_list), pointer :: kl
        integer*4 n 
        if(ktflistqd(k,kl))then
          tfnumlistqnk=kl%head .eq. ktfoper+mtflist
     $         .and. ktfreallistqo(kl) .and. kl%nl .eq. n
          if(present(kl1))then
            kl1=>kl
          endif
        else
          tfnumlistqnk=.false.
        endif
        return
        end function

        logical*4 function tfcomplexnumlistqk(k,klx)
        implicit none
        type (sad_list), pointer, optional, intent(out) :: klx
        type (sad_list), pointer :: kl
        integer*8 k
        integer*4 i
        tfcomplexnumlistqk=.false.
        if(tflistqk(k,kl))then
          if(ktfnonreallistqo(kl))then
            do i=1,kl%nl
              if(.not. tfnumberqd(kl%dbody(i)))then
                return
              endif
            enddo
          endif
          tfcomplexnumlistqk=.true.
          if(present(klx))then
            klx=>kl
          endif
        endif
        return
        end function
        
        logical*4 function ktfsymbolq(k,sym)
        implicit none
        type (sad_symbol), pointer, optional, intent(out) :: sym
        integer*8 k
        ktfsymbolq=iand(ktfmask,k) .eq. ktfsymbol
        if(present(sym) .and. ktfsymbolq)then
          call loc_sym(ktfaddr(k),sym)
        endif
        return
        end function ktfsymbolq

        logical*4 function ktfnonsymbolq(k,sym)
        implicit none
        type (sad_symbol), pointer, optional, intent(out) :: sym
        integer*8 k
        ktfnonsymbolq=iand(ktfmask,k) .ne. ktfsymbol
        if(present(sym) .and. .not. ktfnonsymbolq)then
          call loc_sym(ktfaddr(k),sym)
        endif
        return
        end function ktfnonsymbolq

        logical*4 function ktfsymbolqd(k,sym)
        implicit none
        type (sad_descriptor) k
        type (sad_symbol), pointer, optional, intent(out) :: sym
        ktfsymbolqd=iand(ktfmask,k%k) .eq. ktfsymbol
        if(present(sym) .and. ktfsymbolqd)then
          call loc_sym(ktfaddr(k%k),sym)
        endif
        return
        end function ktfsymbolqd

        logical*4 function ktfnonsymbolqd(k,sym)
        implicit none
        type (sad_descriptor) k
        type (sad_symbol), pointer, optional, intent(out) :: sym
        ktfnonsymbolqd=iand(ktfmask,k%k) .ne. ktfsymbol
        if(present(sym) .and. .not. ktfnonsymbolqd)then
          call loc_sym(ktfaddr(k%k),sym)
        endif
        return
        end function ktfnonsymbolqd

        logical*4 function ktfsymbolqdef(k,symd)
        implicit none
        type (sad_symdef), pointer, optional, intent(out) :: symd
        integer*8 k
        ktfsymbolqdef=iand(ktfmask,k) .eq. ktfsymbol
        if(present(symd) .and. ktfsymbolqdef)then
          call loc_symdef(ktfaddr(k),symd)
        endif
        return
        end function ktfsymbolqdef

        logical*4 function ktfpatq(k,pat)
        implicit none
        type (sad_pat), pointer, optional, intent(out) :: pat
        integer*8 k
        if(iand(ktfmask,k) .eq. ktfpat)then
          ktfpatq=.true.
          if(present(pat))then
            call loc_pat(iand(ktamask,k),pat)
          endif
        else
          ktfpatq=.false.
        endif
        return
        end function ktfpatq

        logical*4 function ktfpatqd(k,pat)
        implicit none
        type (sad_descriptor) k
        type (sad_pat), pointer, optional, intent(out) :: pat
        ktfpatqd=ktfpatq(k%k,pat)
        return
        end function ktfpatqd

        logical*4 function ktfnonpatq(k,pat)
        implicit none
        type (sad_pat), pointer, optional, intent(out) :: pat
        integer*8 k
        ktfnonpatq=iand(ktfmask,k) .ne. ktfpat
        if(present(pat) .and. .not. ktfnonpatq)then
          call loc_pat(ktfaddr(k),pat)
        endif
        return
        end function ktfnonpatq

        logical*4 function ktfrefq(k,ka)
        implicit none
        integer*8 k
        integer*8 , optional, intent(out) :: ka
        ktfrefq=iand(ktfmask,k) .eq. ktfref
        if(ktfrefq .and. present(ka))then
          ka=ktfaddr(k)
        endif
        return
        end function ktfrefq

        logical*4 function ktfrefqd(k,ka)
        implicit none
        type (sad_descriptor) k
        integer*8 , optional, intent(out) :: ka
        ktfrefqd=iand(ktfmask,k%k) .eq. ktfref
        if(ktfrefqd .and. present(ka))then
          ka=ktfaddr(k%k)
        endif
        return
        end function ktfrefqd

        logical*4 function ktfnonrefq(k)
        implicit none
        integer*8 k
        ktfnonrefq=iand(ktfmask,k) .ne. ktfref
        return
        end function ktfnonrefq

        logical*4 function ktfreallistq(ka)
        implicit none
        integer*8 ka
        ktfreallistq=iand(ilist(2,ka-3),lnonreallist) .eq. 0
        return
        end function ktfreallistq

        logical*4 function ktfnonreallistq(ka)
        implicit none
        integer*8 ka
        ktfnonreallistq=iand(ilist(2,ka-3),lnonreallist) .ne. 0
        return
        end function ktfnonreallistq

        logical*4 function ktfreallistqo(list)
        implicit none
        type (sad_list) list
        ktfreallistqo=iand(list%attr,lnonreallist) .eq. 0
        return
        end function ktfreallistqo

        logical*4 function ktfnonreallistqo(list)
        implicit none
        type (sad_list) list
        ktfnonreallistqo=iand(list%attr,lnonreallist) .ne. 0
        return
        end function ktfnonreallistqo

        logical*4 function ktftrueq(ka)
        implicit none
        integer*8 ka
        ktftrueq=ka .ne. 0 .and. iand(ktrmask,ka) .ne. ktfnr
        return
        end function ktftrueq

        logical*4 function ktfsequenceq(k,kl)
        implicit none
        type (sad_list), pointer, optional, intent(out) :: kl
        integer*8 k
        ktfsequenceq=iand(ktfmask,k) .eq. ktflist .and.
     $       klist(iand(ktamask,k)) .eq. ktfoper+mtfnull
        if(ktfsequenceq .and. present(kl))then
          call loc_list(iand(ktamask,k),kl)
        endif
        return
        end function ktfsequenceq

        logical*4 function ktfsequenceqd(k,kl)
        implicit none
        type (sad_list), pointer, optional, intent(out) :: kl
        type (sad_descriptor) k
        ktfsequenceqd=iand(ktfmask,k%k) .eq. ktflist .and.
     $       klist(iand(ktamask,k%k)) .eq. ktfoper+mtfnull
        if(ktfsequenceqd .and. present(kl))then
          call loc_list(iand(ktamask,k%k),kl)
        endif
        return
        end function ktfsequenceqd

        logical*4 function ktfprotectedq(ka)
        implicit none
        integer*8 ka
        ktfprotectedq=iand(ilist(1,ka-3),iattrprotected) .ne. 0
        return
        end function ktfprotectedq

        logical*4 function ktfprotectedqo(sym)
        implicit none
        type (sad_symbol) sym
        ktfprotectedqo=iand(sym%attr,iattrprotected) .ne. 0
        return
        end function ktfprotectedqo

        logical*4 function ktfconstantq(ka)
        implicit none
        integer*8 ka
        ktfconstantq=iand(ilist(1,ka-3),iattrconstant) .ne. 0
        return
        end function ktfconstantq

        logical*4 function ktfconstantsymq(sym)
        implicit none
        type (sad_symbol) sym
        ktfconstantsymq=iand(sym%attr,iattrconstant) .ne. 0
        return
        end function ktfconstantsymq

        logical*4 function ktfimmediateq(ka)
        implicit none
        integer*8 ka
        ktfimmediateq=iand(ilist(1,ka-3),iattrimmediate) .ne. 0
        return
        end function ktfimmediateq

        logical*4 function ktfnumericq(ka)
        implicit none
        integer*8 ka
        ktfnumericq=iand(ilist(1,ka-3),iattrnumeric) .ne. 0
        return
        end function ktfnumericq

        logical*4 function ktfovrwrtq(kl)
        implicit none
        type (sad_list) kl
        ktfovrwrtq=kl%ref .le. 0 .or.
     $       kl%ref .eq. 1 .and. ktfaddr(kl%alloc) .eq. 0
        return
        end function ktfovrwrtq

        logical*4 function tfcomplexqd(k,klx)
        implicit none
        type (sad_descriptor) k
        type (sad_list), pointer :: kl
        type (sad_list), pointer, optional, intent(out) :: klx
        tfcomplexqd= ktflistqd(k,kl) .and.
     $       kl%head .eq. ktfoper+mtfcomplex
     $       .and. kl%nl .eq. 2 .and. ktfreallistqo(kl)
        if(tfcomplexqd .and. present(klx))then
          klx=>kl
        endif
        return
        end function

        logical*4 function tfcomplexqk(k,klx)
        implicit none
        type (sad_list), pointer :: kl
        type (sad_list), pointer, optional, intent(out) :: klx
        integer*8 k
        tfcomplexqk= ktflistq(k,kl) .and.
     $       kl%head .eq. ktfoper+mtfcomplex
     $       .and. kl%nl .eq. 2 .and. ktfreallistqo(kl)
        if(tfcomplexqk .and. present(klx))then
          klx=>kl
        endif
        return
        end function

        logical*4 function tfcomplexqx(k,cx)
        implicit none
        type (sad_complex), pointer :: c
        type (sad_complex), pointer, optional, intent(out) :: cx
        integer*8 k
        tfcomplexqx=ktflistqx(k,c) .and.
     $       c%head .eq. ktfoper+mtfcomplex
     $       .and. c%nl .eq. 2 .and.
     $       iand(lnonreallist,c%attr) .eq. 0
        if(tfcomplexqx .and. present(cx))then
          cx=>c
        endif
        return
        end function

        logical*4 function tfmatrixqd(k,kl,klind,klbody)
        implicit none
        type (sad_descriptor) k
        type (sad_list), pointer, optional,
     $       intent(out) :: kl,klind,klbody
        type (sad_list), pointer :: kl1,kli1,klb1
        logical*4 tfsameqd
        if(ktflistqd(k,kl1))then
          if(tfsameqd(kl1%dbody(0),kxmatrix) .and. kl1%nl .eq. 2)then
            if(ktflistqd(kl1%dbody(1),kli1) .and.
     $           ktflistqd(kl1%dbody(2),klb1) .and.
     $           ktfreallistqo(kli1))then
              tfmatrixqd=.true.
              if(present(kl))then
                kl=>kl1
                if(present(klind))then
                  klind=>kli1
                  if(present(klbody))then
                    klbody=>klb1
                  endif
                endif
              endif
              return
            endif
          endif
        endif
        tfmatrixqd=.false.
        return
        end function

        complex*16 function cfromr(r)
        implicit none
        real*8 r(2)
        cfromr=dcmplx(r(1),r(2))
        return
        end function cfromr

        integer*8 function ktfcopy1(k)
        implicit none
        integer*8 k,ka
        ka=iand(ktamask,k)
        ilist(1,ka-1)=ilist(1,ka-1)+1
        ktfcopy1=k
        return
        end function ktfcopy1

        integer*8 function ktfcopyd(k,d)
        implicit none
        integer*8 k,ka
        logical*4 d
        d=ktfobjq(k)
        if(d)then
          ka=iand(ktamask,k)
          ilist(1,ka-1)=ilist(1,ka-1)+1
        else
          d=ktfnonrealq(k)
        endif
        ktfcopyd=k
        return
        end function ktfcopyd

        integer*8 function ktfcopy(k)
        implicit none
        integer*8 k,ka
        if(ktfobjq(k))then
          ka=ktfaddr(k)
          ilist(1,ka-1)=ilist(1,ka-1)+1
        endif
        ktfcopy=k
        return
        end function

        type (sad_descriptor) function dtfcopy(d)
        implicit none
        type (sad_descriptor) d
        dtfcopy%k=ktfcopy(d%k)
        return
        end function

        type (sad_descriptor) function dtfcopy1(d)
        implicit none
        type (sad_descriptor) d
        dtfcopy1%k=ktfcopy1(d%k)
        return
        end function

        subroutine tflocald(k)
        implicit none
        type (sad_descriptor) k
        type (sad_object), pointer :: obj
        integer*8 itfroot
        if(ktfobjqd(k,obj))then
          obj%ref=obj%ref-1
          if(obj%ref .le. 0)then
            obj%ref=0
            if(ktfaddr(obj%alloc) .eq. 0)then
              itfroot=itflocal+levele
              obj%alloc=obj%alloc+ktfaddr(klist(itfroot))
              klist(itfroot)=ktfaddr(k%k)-2
            endif
          endif
        endif
        return
        end subroutine 

        subroutine tflocal1d(k)
        implicit none
        type (sad_descriptor) k
        type (sad_object), pointer :: obj
        integer*8 ka,itfroot
        ka=ktfaddr(k)
        call loc_obj(ka,obj)
        obj%ref=obj%ref-1
        if(obj%ref .le. 0)then
          obj%ref=0
          if(ktfaddr(obj%alloc) .eq. 0)then
            itfroot=itflocal+levele
            obj%alloc=obj%alloc+ktfaddr(klist(itfroot))
            klist(itfroot)=ka-2
          endif
        endif
        return
        end subroutine

        subroutine tfconnect(k,irtc)
        implicit none
        type (sad_descriptor) k
        integer*4 irtc
        call tfconnectk(k%k,irtc)
        return
        end subroutine

        subroutine tfconnectk(k,irtc)
        implicit none
        type (sad_object), pointer :: obj
        integer*4 l,itfdownlevel,irtc
        integer*8 k,ka,j
        if(levele .gt. 0)then
          if(irtc .ne. 0)then
            l=itfdownlevel()
          elseif(ktfobjq(k))then
            ka=ktfaddr(k)
            call loc_obj(ktfaddr(ka),obj)
            obj%ref=obj%ref+1
            l=itfdownlevel()
            call tflocal1(ka)
            if(ktfaddr(obj%alloc) .eq. 0)then
              j=itflocal+levele
              obj%alloc=ktftype(obj%alloc)+klist(j)
              klist(j)=sad_loc(obj%alloc)
            endif
c     call tfdebugprint(ktftype(klist(ka-2))+ka,'tfconnectk',1)
c     write(*,*)'with ',ilist(1,ka-1),ktfaddr(klist(ka-2))
          else
            l=itfdownlevel()
          endif
        endif
        return
        end subroutine

        subroutine tfconnectk1(k,irtc)
        implicit none
        type (sad_object), pointer :: obj
        integer*4 l,itfdownlevel,irtc
        integer*8 k,ka
        if(levele .gt. 0)then
          if(irtc .ne. 0)then
            l=itfdownlevel()
          elseif(ktfobjq(k))then
            ka=ktfaddr(k)
            call loc_obj(ktfaddr(ka),obj)
            obj%ref=obj%ref+1
            l=itfdownlevel()
            call tflocal1(ka)
          else
            l=itfdownlevel()
          endif
        endif
        return
        end subroutine

        type (sad_descriptor) function kxaaloc(mode,nd,kl)
        implicit none
        type (sad_list), pointer, optional, intent(out) :: kl
        integer*4 mode,nd
        kxaaloc%k=ktflist+ktaaloc(mode,nd,kl)
        return
        end function

        type (sad_descriptor) function kxadaloc(mode,nd,kl)
        implicit none
        type (sad_list), pointer, optional, intent(out) :: kl
        integer*4 mode,nd
        kxadaloc%k=ktflist+ktadaloc(mode,nd,kl)
        return
        end function

        type (sad_descriptor) function kxadalocnull(mode,nd,kl)
        implicit none
        type (sad_list), pointer, optional, intent(out) :: kl
        integer*4 mode,nd
        kxadalocnull%k=ktflist+ktadalocnull(mode,nd,kl)
        return
        end function

        type (sad_descriptor) function kxavaloc(mode,nd,kl)
        implicit none
        type (sad_list), pointer, optional, intent(out) :: kl
        integer*4 mode,nd
        kxavaloc%k=ktflist+ktavaloc(mode,nd,kl)
        return
        end function

        type (sad_descriptor) function kxmakelist(isp1,kl)
        implicit none
        type (sad_list), pointer, optional, intent(out) :: kl
        integer*4 isp1
        kxmakelist%k=ktflist+ktfmakelist(isp1,kl)
        return
        end function

        type (sad_descriptor) function kxmakelist0(isp1,kl)
        implicit none
        type (sad_list), pointer, optional, intent(out) :: kl
        integer*4 isp1
        kxmakelist0%k=ktflist+ktfmakelist0(isp1,kl)
        return
        end function

        type (sad_descriptor) function kxcompose(isp1,kl)
        implicit none
        type (sad_list), pointer, optional, intent(out) :: kl
        integer*4 isp1
        kxcompose%k=ktflist+ktfcompose(isp1,kl)
        return
        end function

        type (sad_descriptor) function kxcomposer(isp1)
        implicit none
        integer*4 isp1
        integer*8 ktfcrelistr
        kxcomposer%k=ktflist+
     $       ktfcrelistr(isp-isp1,ktastk(isp1+1),ktastk(isp1))
        return
        end function

        type (sad_descriptor) function kxcomposev(isp0,kl)
        implicit none
        type (sad_list), pointer, optional, intent(out) :: kl
        integer*4 isp0
        kxcomposev%k=ktflist+ktfcomposev(isp0,kl)
        return
        end function

        type (sad_descriptor) function kxcalocv(mode,x,y)
        implicit none
        integer*4 mode
        real*8 x,y
        kxcalocv%k=ktflist+ktcalocv(mode,x,y)
        return
        end function

        type (sad_descriptor) function kxsalocb(mode,string,leng,str)
        implicit none
        type (sad_string), pointer, optional, intent(out) :: str
        integer*4 mode,leng
        character string(leng)
        kxsalocb%k=ktfstring+ktsalocb(mode,string,leng)
        if(present(str))then
          call descr_sad(kxsalocb,str)
        endif
        return
        end function

        type (sad_descriptor) function kxsalocbb(mode,leng,str)
        implicit none
        type (sad_string), pointer, optional, intent(out) :: str
        integer*4 mode,leng
        kxsalocbb%k=ktfstring+ktsalocbb(mode,leng)
        if(present(str))then
          call descr_sad(kxsalocbb,str)
        endif
        return
        end function

        type (sad_descriptor) function kxsymbolz(name,l,symd)
        implicit none
        type (sad_symdef), pointer, optional, intent(out) :: symd
        integer*4 l
        character name(l)
        kxsymbolz%k=ktfsymbol+ktfsymbolz(name,l,symd)
        return
        end function

        type (sad_descriptor) function kxsymbolv(name,l,symd0)
        implicit none
        type (sad_symdef), pointer, optional, intent(out) :: symd0
        type (sad_symdef), pointer :: symd
        integer*8 k
        integer*4 l
        character name(l)
        k=ktfsymbolz(name,l,symd)
        kxsymbolv=symd%value
        if(present(symd0))then
          symd0=>symd
        endif
        return
        end function

        type (sad_descriptor) function kxsymbolf(name,l,const,symd0)
        implicit none
        type (sad_symdef), pointer, optional, intent(out) :: symd0
        type (sad_symdef), pointer :: symd
        integer*4 l
        character name(l)
        logical*4 const
        kxsymbolf=kxsymbolz(name,l,symd)
        if(present(symd0))then
          symd0=>symd
        endif        
        if(const .and. symd%sym%gen .le. 0)then
          symd%sym%attr=ior(symd%sym%attr,iattrconstant+iattrprotected)
        endif
        return
        end function

        subroutine tfsetsymbolr(name,l,val)
        implicit none
        integer*4 l
        character name(l)
        real*8 val
        type (sad_symdef), pointer :: symd
        type (sad_descriptor) kx
        kx=kxsymbolz(name,l,symd)
        symd%value=dfromr(val)
        return
        end subroutine

        type (sad_descriptor) function k_descr(k)
        implicit none
        integer*8 k
        k_descr%k=k
        return
        end function

        type (sad_descriptor) function kxscopy(ka,leng,str)
        implicit none
        type (sad_string), pointer, optional, intent(out) :: str
        integer*8 ka,ktfaloc,l
        integer*4 leng,nw,n
        if(leng .eq. 0)then
          kxscopy=dxnulls
          go to 10
        endif
        nw=leng/8+2
        l=ktfaloc(-1,ktfstring,nw)
        ilist(2,l-3)=-1
        ilist(1,l)=leng
        klist(l+nw-1)=0
        n=(min(leng,ilist(1,ka))+7)/8
        klist(l+1:l+n)=klist(ka+1:ka+n)
        kxscopy%k=ktfstring+l
 10     if(present(str))then
          call descr_sad(kxscopy,str)
        endif
        return
        end function

        type (sad_descriptor) function kxnaloc1(lg,locp)
        implicit none
        integer*8 locp
        integer*4 lg
        kxnaloc1=kxnaloc(lg,locp,0)
        return
        end function

        type (sad_descriptor) function kxnaloc(lg,locp,n)
        implicit none
        integer*8 kp,locp,kp1, kp0,ktalocr
        integer*4 lg,ipg,n
        type (sad_symdef), pointer :: def,def0
        type (sad_namtbl), pointer :: loc
        call loc_namtbl(locp,loc)
        kp0=sad_loc(loc%symdef)
        kp=loc%symdef
        if(kp .le. 0)then
          kp=ktalocr(9+n)
          call loc1_symdef(kp,def)
          def%next=0
          def%prev=kp0
          def%upval=0
          def%downval=0
          def%value%k=ktfsymbol+kp+8
          def%sym%attr=0
          def%sym%override=-2
          def%sym%alloc%k=ktfsymbol
          def%sym%ref=2
          def%sym%gen=lg
          def%sym%loc=locp
          loc%symdef=kp
          loc%str%ref=loc%str%ref+1
        else
          call loc1_symdef(kp,def0)
          ipg=def0%sym%gen
          do while(lg .lt. ipg)
            kp0=kp
            kp=def0%next
            if(kp .eq. 0)then
              exit
            endif
            call loc1_symdef(kp,def0)
            ipg=def0%sym%gen
          enddo
          kp1=ktaloc(9+n)
          call loc1_symdef(kp1,def)
          def%next=kp
          def%prev=kp0
          if(kp .ne. 0)then
            def0%prev=kp1
            if(lg .eq. ipg)then
              def0%sym%override=0
            endif
          endif
          klist(kp0)=kp1
          def%upval=0
          def%downval=0
          def%value%k=ktfsymbol+kp1+8
          def%sym%attr=0
          def%sym%override=-2
          def%sym%alloc%k=ktfsymbol
          def%sym%ref=2
          def%sym%gen=lg
          def%sym%loc=locp
        endif
        kxnaloc%k=ktfsymbol+sad_loc(def%sym%loc)
        return
        end function

        real*8 function rfromd(d)
        implicit none
        type (sad_descriptor) d
        real*8 rfromk
        rfromd=rfromk(d%k)
        return
        end function

        integer*4 function ifromd(d)
        implicit none
        type (sad_descriptor) d
        real*8 rfromk
        ifromd=int(rfromk(d%k))
        return
        end function

        type (sad_descriptor) function dfromr(x)
        implicit none
        integer*8 kfromr
        real*8 x
        dfromr%k=kfromr(x)
        return
        end function

        type (sad_descriptor) function dfromk(k)
        implicit none
        integer*8 k
        dfromk%k=k
        return
        end function

        integer*8 function ktfmakelist0(isp1,kl)
        implicit none
        type (sad_list), pointer, optional, intent(out) :: kl
        integer*4 isp1,narg
        narg=isp-isp1
        if(narg .le. 0)then
          ktfmakelist0=ktfaddr(kxnulll)
          go to 10
        elseif(narg .eq. 1)then
          if(ktastk(isp) .eq. ktfoper+mtfnull)then
            ktfmakelist0=ktfaddr(kxnulll)
            go to 10
          endif
        endif
        ktfmakelist0=ktfaddr(
     $       kxcrelistm(narg,ktastk(isp1+1:isp1+narg),
     $       k_descr(ktfoper+mtflist)))
 10     if(present(kl))then
          call loc_sad(ktfmakelist0,kl)
        endif
        return
        end function

        type (sad_descriptor) function kxcrelistm(m,ks,kh)
        implicit none
        type (sad_descriptor) kh
        type (sad_list), pointer ::kl
        integer*4 m
        integer*8 ks(m)
        kxcrelistm=kxaaloc(-1,m,kl)
        call tfcrelista(m,ks,kh,kl)
        return
        end function

        integer*8 function ktsalocb(mode,string,leng)
        integer*8 ktfaloc,l
        integer*4 leng,nw,mode,ik
        character string(leng)
        if(leng .eq. 1)then
          ik=ichar(string(1))
          ktsalocb=iaxschar+ik*5+3
          if(mode .eq. 0)then
            ilist(1,ktsalocb-1)=ilist(1,ktsalocb-1)+1
          endif
        elseif(leng .eq. 0)then
          ktsalocb=ktfaddr(kxnulls)
          if(mode .eq. 0)then
            ilist(1,ktsalocb-1)=ilist(1,ktsalocb-1)+1
          endif
        else
          nw=leng/8+2
          l=ktfaloc(mode,ktfstring,nw)
          call tfpadstr(string,l+1,leng)
          ilist(1,l)=leng
          ilist(2,l)=0
          ilist(2,l-3)=-1
          ktsalocb=l
        endif
        return
        end function
      
        integer*8 function ktsalocbb(mode,leng)
        implicit none
        integer*8 ktfaloc,l
        integer*4 leng,nw,mode
        if(leng .eq. 0)then
          ktsalocbb=kxnulls
          if(mode .eq. 0)then
            ilist(1,ktfaddr(kxnulls)-1)=
     $           ilist(1,ktfaddr(kxnulls)-1)+1
          endif
          return
        endif
        nw=leng/8+2
        l=ktfaloc(mode,ktfstring,nw)
        ilist(2,l-3)=-1
        ilist(1,l)=leng
        klist(l+nw-1)=0
        ktsalocbb=l
        return
        end function

        integer*8 function ktcalocv(mode,x,y)
        implicit none
        type (sad_list), pointer :: kl
        integer*4 mode
        real*8 x,y
        ktcalocv=ktavaloc(mode,2,kl)
        kl%attr=lconstlist
        kl%head=ktfoper+mtfcomplex
        kl%rbody(1)=x
        kl%rbody(2)=y
        return
        end function

        integer*8 function ktadalocnull(mode,nd,kl)
        implicit none
        integer*4 mode,nd
        type (sad_list), pointer :: kl1
        type (sad_list), pointer, optional, intent(out) :: kl
        ktadalocnull=ktadaloc(mode,nd,kl1)
        kl1%body(1:nd)=ktfoper+mtfnull
        if(present(kl))then
          kl=>kl1
        endif
        return
        end function

        integer*8 function ktaaloc(mode,nd,kl)
        implicit none
        type (sad_list), pointer, optional, intent(out) :: kl
        integer*8 ka,ktfaloc
        integer*4 nd,mode
        ka=ktfaloc(mode,ktflist,nd+1)
        ilist(1,ka-3)=0
        ilist(2,ka-1)=nd
        klist(ka)=ktfoper+mtflist
        ktaaloc=ka
        if(present(kl))then
          call loc_list(ka,kl)
        endif
        return
        end function

        integer*8 function ktavaloc(mode,nd,kl)
        implicit none
        type (sad_list), pointer, optional, intent(out) :: kl
        integer*4 mode,nd
        integer*8 k1,itfroot
        k1=ktaloc(nd+3)
        ilist(1,k1-1)=0
        ilist(2,k1-1)=kconstarg
        if(mode .eq. 0)then
          klist(k1)=ktflist
        else
          itfroot=itflocal+levele
          klist(k1)=ktflist+klist(itfroot)
          klist(itfroot)=k1
        endif
        ilist(1,k1+1)=mode+1
        ilist(2,k1+1)=nd
        klist(k1+2)=ktfoper+mtflist
        ktavaloc=k1+2
        if(present(kl))then
          call loc_list(k1+2,kl)
        endif
        return
        end function

        integer*8 function ktadaloc(mode,nd,kl)
        implicit none
        type (sad_list), pointer, optional, intent(out) :: kl
        integer*4 mode,nd
        integer*8 k1,itfroot
        k1=ktaloc(nd+3)
        klist(k1-1)=0
        ilist(2,k1-1)=lnonreallist
        if(mode .eq. 0)then
          klist(k1)=ktflist
        else
          itfroot=itflocal+levele
          klist(k1)=ktflist+klist(itfroot)
          klist(itfroot)=k1
        endif
        ilist(1,k1+1)=mode+1
        ilist(2,k1+1)=nd
        klist(k1+2)=ktfoper+mtflist
        ktadaloc=k1+2
        if(present(kl))then
          call loc_list(k1+2,kl)
        endif
        return
        end function

        integer*8 function ktaalocr(mode,nd,kl)
        implicit none
        type (sad_list), pointer, optional, intent(out) :: kl
        integer*8 ka,ktfalocr
        integer*4 nd,mode
        ka=ktfalocr(mode,ktflist,nd+1)
        ilist(1,ka-3)=0
        ilist(2,ka-1)=nd
        klist(ka)=ktfoper+mtflist
        ktaalocr=ka
        if(present(kl))then
          call loc_list(ka,kl)
        endif
        return
        end function

        integer*8 function ktraaloc(mode,nd,kl)
        implicit none
        type (sad_list), pointer, optional, intent(out) :: kl
        integer*4 mode,nd
        integer*8 ka
        ka=ktavaloc(mode,nd,kl)
        klist(ka+1:ka+nd)=0
        ktraaloc=ka
        return
        end function

        integer*8 function ktaalocsp(nd,lp,la,kl1)
        implicit none
        type (sad_list), pointer, optional, intent(out) :: kl1
        type (sad_list), pointer :: kl
        integer*8 ka
        integer*4 nd
        integer*2 lp,la
        ka=ktaloc(nd+lp+la+3)+lp+2
        call loc_list(ka,kl)
        kl%lenp=lp
        kl%lena=la
        kl%attr=0
        kl%alloc=ktflist
        kl%ref=1
        kl%nl=nd
        call tflocal1(ka)
        ktaalocsp=ka
        if(present(kl1))then
          kl1=>kl
        endif
        return
        end function

        integer*8 function ktfcomposev(isp0,kl)
        implicit none
        type (sad_list), pointer, optional, intent(out) :: kl
        integer*8 kax,kh
        integer*4 isp0
        if(isp .gt. isp0)then
          ktfcomposev=ktfcompose(isp0)
        else
          kh=ktastk(isp0)
          if(kh .eq. ktfoper+mtflist)then
            ktfcomposev=ktfaddr(kxnulll)
          elseif(kh .eq. ktfoper+mtfnull)then
            ktfcomposev=ktfaddr(kxnull)
          elseif(kh.eq. ktfoper+mtfslot)then
            ktfcomposev=ktfaddr(klist(iaxslotnull))
          elseif(ktastk(isp0) .eq. ktfoper+mtfslotseq)then
            ktfcomposev=ktfaddr(klist(iaxslotnull+1))
          else
            kax=ktaaloc(-1,0,kl)
            klist(kax)=ktfcopy(kh)
            ktfcomposev=kax
            return
          endif
        endif
        if(present(kl))then
          call loc_list(ktfcomposev,kl)
        endif
        return
        end function

        integer*8 function ktfcompose(isp1,kl)
        implicit none
        type (sad_list), pointer, optional, intent(out) :: kl
        integer*4 isp1
        ktfcompose=ktfaddr(
     $       kxcrelistm(isp-isp1,ktastk(isp1+1:isp),dtastk(isp1)))
        if(present(kl))then
          call loc_list(ktfcompose,kl)
        endif
        return
        end function

        integer*8 function ktfmakelist(isp1,kl)
        implicit none
        type (sad_list), pointer, optional, intent(out) :: kl
        integer*4 isp1,narg
        narg=isp-isp1
        if(narg .eq. 1)then
          if(ktastk(isp) .eq. ktfoper+mtfnull)then
            ktfmakelist=ktaaloc(-1,0,kl)
            return
          endif
        endif
        ktfmakelist=ktfaddr(kxcrelistm(narg,ktastk(isp1+1:isp1+narg),
     $       k_descr(ktfoper+mtfleftbrace)))
        if(present(kl))then
          call loc_list(ktfmakelist,kl)
        endif
        return
        end function

        integer*8 function ktfsymbolz(name,l,symd)
        implicit none
        type (sad_symdef), pointer, optional, intent(out) :: symd
        integer*4 l
        integer*8 ktfsymbolc
        character name(l)
        ktfsymbolz=ktfsymbolc(name,l,int8(0))
        if(present(symd))then
          call loc_symdef(ktfsymbolz,symd)
        endif
        return
        end function

        subroutine tfsydef(sym,symx)
        implicit none
        type (sad_descriptor) kx
        type (sad_symbol) sym
        type (sad_symbol), pointer, intent(out) :: symx
        call tfsydefg(sym%loc,kx,sym%gen)
        call loc_sad(ktfaddrd(kx),symx)
        return
        end subroutine

        type (sad_descriptor) function dxsycopy(sym)
        implicit none
        type (sad_symbol) sym
        integer*8 kax,ktfaloc
        kax=ktfaloc(-1,ktfsymbol,1)
        ilist(2,kax-1)=maxgeneration
        klist(kax)=sym%loc
        dxsycopy%k=ktfsymbol+kax
        return
        end function

        type (sad_descriptor) function kxpfaloc(k)
        implicit none
        type (sad_descriptor) k
        type (sad_list), pointer :: klx
        if(ktfrealqd(k))then
          kxpfaloc=kxavaloc(-1,1,klx)
          klx%dbody(1)=k
        else
          kxpfaloc=kxadaloc(-1,1,klx)
          klx%dbody(1)=dtfcopy(k)
        endif
        klx%head=ktfoper+mtffun
        return
        end function

        type (sad_descriptor) function kxpaloc(string)
        implicit none
        character*(*) string
        type (sad_descriptor) kh
        integer*8 ipk
        integer*4 lenw,l,ip,k
        l=lenw(string)
        ip=index(string(1:l),'_')
        if(ip .gt. 0)then
          if(string(ip:min(l,ip+2)) .eq. '___')then
            k=3
          elseif(string(ip:min(l,ip+1)) .eq. '__')then
            k=2
          else
            k=1
          endif
        else
          k=0
          ip=l+1
        endif
        ipk=ip+k
        if(ipk .gt. l)then
          kh%k=ktfref
        else
          kh%k=ktfsymbol+ktfsymbolz(string(ipk:l),int(l-ipk+1))
        endif
        kxpaloc=kxpalocb(string(1:ip-1),ip-1,transfer(ktfref+k,kh),kh)
        return
        end function

        type (sad_descriptor) function kxpalocb(symb,ls,kp,kh)
        implicit none
        type (sad_descriptor) kp,kh,ks
        integer*4 ls
        character symb(ls)
        if(ls .gt. 0)then
          ks=kxsymbolz(symb,ls)
        else
          ks%k=0
        endif
        kxpalocb=kxpcopyss(kp,kh,ks,transfer(ktfref,kp))
        return
        end function

        type (sad_descriptor) function kxpcopyss(kp,kh,ks,kd)
        use tfcode
        use iso_c_binding
        implicit none
        type (sad_descriptor) kp,kh,kd,ks
        type (sad_pat), pointer :: pat
        type (sad_symbol), pointer :: sym
        integer*8 kax,ktfaloc
        kax=ktfaloc(-1,ktfpat,9)
        call loc_pat(kax,pat)      
        pat%expr=dtfcopy(kp)
        pat%head=dtfcopy(kh)
        nullify(pat%equiv)
        pat%default=dtfcopy(kd)
        pat%value%k=ktfref
        pat%sym%attr=pat%len-7
        pat%sym%override=1
        pat%sym%ref=1
        if(ktfsymbolqd(ks,sym))then
          pat%sym%alloc=dtfcopy1(ks)
          pat%sym%gen=sym%gen
          pat%sym%loc=sym%loc
        else
          pat%sym%alloc%k=ktfsymbol
          pat%sym%gen=0
          pat%sym%loc=0
        endif
        kxpcopyss%k=ktfpat+kax
        return
        end function

        type (sad_descriptor) function kxsubstring(kh,isp1,isp2)
        implicit none
        type (sad_descriptor) kh,kx
        type (sad_string), pointer :: str
        integer*4 n,ic1,ic2,isp1,isp2
        call descr_sad(kh,str)
        n=str%nch
        ic1=int(rtastk(isp1))
        if(ic1 .lt. 0)then
          ic1=n+1+ic1
        endif
        ic2=int(rtastk(isp2))
        if(ic2 .lt. 0)then
          ic2=n+1+ic2
        endif
        if(ic2 .lt. ic1 .or. ic2 .le. 0 .or. ic1 .gt. n)then
          kx%k=kxnulls
        else
          ic1=max(1,ic1)
          ic2=min(n,ic2)
          kx=kxsalocb(-1,str%str(ic1:ic2),ic2-ic1+1)
        endif
        kxsubstring=kx
        return
        end function

        subroutine tfclonelist(list,listc)
        implicit none
        type (sad_list), target :: list
        type (sad_list), pointer, intent(out) :: listc
        integer*4 i
        if(ktfovrwrtq(list))then
          listc=>list
        else
          call loc_list(ktaaloc(-1,list%nl),listc)
          listc%attr=list%attr
          if(ktfreallistqo(list))then
            listc%head=ktfcopy(list%head)
            listc%body(1:list%nl)=list%body(1:list%nl)
          else
            do i=0,list%nl
              listc%body(i)=ktfcopy(list%body(i))
            enddo
          endif
        endif
        listc%attr=ior(listc%attr,ktoberebuilt)
        return
        end subroutine

        subroutine tfduplist(list,listc)
        implicit none
        type (sad_list), target :: list
        type (sad_list), pointer, intent(out) :: listc
        integer*4 i
        call loc_list(ktaaloc(-1,list%nl),listc)
        listc%attr=list%attr
        if(ktfreallistqo(list))then
          listc%head=ktfcopy(list%head)
          listc%body(1:list%nl)=list%body(1:list%nl)
        else
          do i=0,list%nl
            listc%body(i)=ktfcopy(list%body(i))
          enddo
        endif
        return
        end subroutine

        logical*4 function tfonstackq(ka)
        implicit none
        integer*8 ka
        tfonstackq=ka .ge. isporg+ispbase
     $       .and. ka .le. isporg+ivstkoffset*2+ispbase
        return
        end function

        subroutine tfgetdefargp(kl,kas,kp,ev,irtc)
        implicit none
        type (sad_list) kl
        type (sad_descriptor) kv
        integer*8 kp,kas
        integer*4 isp0,irtc
        logical*4 ev
        isp=isp+1
        isp0=isp
        ktastk(isp)=kl%head
        call tfgetllstkall(kl)
        call tfdeval(isp0,kas,kv,1,.true.,ev,irtc)
        kp=ktfaddrd(kv)
        isp=isp0-1
        return
        end subroutine

        recursive subroutine tfrebuildl(kl,klx,rep)
        implicit none
        type (sad_list), target :: kl
        type (sad_list), pointer :: kli,klxi
        type (sad_list), pointer, intent(out) :: klx
        integer*8 kax
        integer*4 i,isp0
        logical*4 rep,rep1
        rep=.false.
        if(iand(kl%attr,ktoberebuilt) .eq. 0)then
          klx=>kl
          return
        endif
        kl%attr=kl%attr-ktoberebuilt
        if(ktfreallistqo(kl))then
          klx=>kl
          return
        endif
        isp=isp+1
        isp0=isp
        ktastk(isp)=kl%head
        do i=1,kl%nl
          if(ktfsequenceq(kl%body(i),kli))then
            call tfgetllstkall(kli)
            rep=.true.
          elseif(ktflistq(kl%body(i),kli))then
            call tfrebuildl(kli,klxi,rep1)
            isp=isp+1
            if(rep1)then
              ktastk(isp)=ktflist+sad_loc(klxi%head)
              rep=.true.
            else
              ktastk(isp)=kl%body(i)
            endif
          else
            isp=isp+1
            ktastk(isp)=kl%body(i)
          endif
        enddo
        if(rep)then
          kax=ktfcompose(isp0,klx)
        else
          klx=>kl
        endif
        isp=isp0-1
        return
        end subroutine 

        subroutine tfmatrixmaybeq(k,cmplm,realm,vec,n,m,kl)
        implicit none
        type (sad_descriptor) k
        type (sad_list), pointer, intent(out) :: kl
        type (sad_list), pointer :: kli
        integer*4 i,n,m
        logical*4 cmplm,realm,vec
        n=0
        m=0
        realm=.false.
        cmplm=.false.
        vec=.false.
        if(ktflistqd(k,kl))then
          if(kl%head .eq. ktfoper+mtflist)then
            n=kl%nl
            if(ktfnonreallistqo(kl))then
              do i=1,n
                if(ktfnonlistq(kl%body(i)))then
                  return
                endif
                call loc_list(ktfaddr(kl%body(i)),kli)
                if(kli%head .ne. ktfoper+mtflist)then
                  return
                endif
                if(i .eq. 1)then
                  m=kli%nl
                elseif(m .ne. kli%nl)then
                  m=0
                  return
                endif
                if(ktfnonreallistqo(kli))then
                  cmplm=.true.
                  return
                endif
              enddo
              realm=.true.
            else
              vec=.true.
            endif
          endif
        endif
        return
        end subroutine

        type (sad_descriptor) function kxm2l(a,n,m,nd,trans)
        implicit none
        type (sad_descriptor) kx,ki
        type (sad_list), pointer :: kli,klx
        integer*4 n,m,nd,i
        logical*4 trans
        real*8 a(nd,m)
        if(n .eq. 0)then
          kx=kxavaloc(-1,m,klx)
          klx%rbody(1:m)=a(1:m,1)
        else
          if(trans)then
            kx=kxadaloc(-1,m,klx)
            do i=1,m
              ki=kxavaloc(0,n,kli)
              kli%rbody(1:n)=a(1:n,i)
              kli%attr=ior(lconstlist,kli%attr)
              klx%dbody(i)=ki
            enddo
          else
            kx=kxadaloc(-1,n,klx)
            do i=1,n
              ki=kxavaloc(0,m,kli)
              kli%rbody(1:m)=a(i,:)
              kli%attr=ior(lconstlist,kli%attr)
              klx%dbody(i)=ki
            enddo
          endif
        endif
        klx%attr=ior(klx%attr,lconstlist)
        kxm2l=kx
        return
        end function

        type (sad_descriptor) function kxcopylist(k)
        implicit none
        type (sad_descriptor) k
        type (sad_list), pointer :: kl,klx
        integer*4 m,i
        call descr_list(k,kl)
        m=kl%nl
        kxcopylist=kxaaloc(-1,m,klx)
        if(ktfreallistqo(kl))then
          klx%head=ktfcopy(kl%head)
          klx%body(1:m)=kl%body(1:m)
        else
          do i=0,m
            klx%dbody(i)=dtfcopy(kl%dbody(i))
          enddo
        endif
        klx%attr=kl%attr
        return
        end function

        type (sad_descriptor) function kxargsym(n0)
        implicit none
        integer*4, parameter :: nsym=1024 
        type (sad_descriptor), save :: ksym(nsym)
        integer*4 n0,n,l,n1,ls,ifrac
        character ch
        character*32 name,buf
        data ksym%k /nsym*0/
        data name /'`Class`s                        '/
        if(ksym(n0)%k .ne. 0)then
          kxargsym=ksym(n0)
          return
        endif
        n=n0
        l=32
        do while(n .ne. 0)
          n1=n/62
          ifrac=n-n1*62
          if(ifrac .lt. 10)then
            ch=char(ichar('0')+ifrac)
          elseif(ifrac .lt. 36)then
            ch=char(ichar('a')+ifrac-10)
          else
            ch=char(ichar('A')+ifrac-36)
          endif
          buf(l:l)=ch
          l=l-1
          n=n1
        enddo
        name(9:9+31-l)=buf(l+1:32)
        ls=9+32-l
        name(ls:ls)='$'
        kxargsym=kxsymbolz(name,ls)
        ksym(n0)=kxargsym
        return
        end function

        real*8 pure function p2h(p)
        implicit none
        real*8, intent(in) :: p
        real*8 p2
        real*8, parameter:: pth=1.d3;
        if(p .gt. pth)then
          p2=1.d0/p**2
          p2h=p*(1.d0+p2*(0.5d0-p2*.125d0))
        else
          p2h=sqrt(1.d0+p**2)
        endif
        return
        end function

        real*8 pure function h2p(h)
        implicit none
        real*8, intent(in) :: h
        real*8 h2
        real*8, parameter:: hth=1.d3;
        if(h .gt. hth)then
          h2=-1.d0/h**2
          h2p=h*(1.d0+h2*(0.5d0-h2*.125d0))
        else
          h2p=sqrt(h**2-1.d0)
        endif
        return
        end function

        real*8 pure function pxy2dpz(px,py)
        implicit none
        real*8, intent(in) :: px,py
        real*8 x
        real*8, parameter:: xth=1.d-6,xmin=1.d-100
        x=px**2+py**2
        if(x .lt. xth)then
          pxy2dpz=-x*(0.5d0+x*(0.125d0+x*0.0625d0))
        else
          pxy2dpz=-x/(1.d0+sqrt(max(xmin,1.d0-x)))
        endif
        return
        end function

        real*8 pure function sqrt1(x)
         implicit none
        real*8, intent(in) :: x
        real*8, parameter:: xth=1.d-6,xmin=1.d-100
        if(abs(x) .lt. xth)then
          sqrt1=x*(0.5d0-x*(0.125d0-x*0.0625d0))
        else
          sqrt1=x/(1.d0+sqrt(max(xmin,1.d0+x)))
        endif
        return
        end function

        real*8 function sqrt1n(x)
        implicit none
        real*8, intent(in) :: x
        sqrt1n=x*(0.5d0-x*(0.125d0-x*0.0625d0))
        sqrt1n=(sqrt1n**2+x)/(2.d0+2.d0*sqrt1n)
        sqrt1n=(sqrt1n**2+x)/(2.d0+2.d0*sqrt1n)
        return
        end function

        subroutine resetnan(a)
        implicit none
        real*8 a(:)
        integer*4 i
        do i=1,size(a)
          if(isnan(a(i)))then
            a(i)=0.d0
          endif
        enddo
        return
        end subroutine

        real*8 function sqrtl(x)
        implicit none
        real*8 x
        real*8 ,parameter :: am=1.d-20
        sqrtl=sqrt(max(x,am))
        return
        end function

      end module

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
