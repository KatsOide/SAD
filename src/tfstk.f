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
      public
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
c$$$      character*(MAXPNAME) pname(HTMAX)
c$$$      integer*4 lpname(HTMAX)
c$$$      integer*4 idtype(HTMAX)
c$$$      integer*8 idval(HTMAX)
      character*(MAXPNAME), dimension(:), 
     $     allocatable, target:: pname
      integer*4 , dimension(:), 
     $     allocatable, target :: lpname
      integer*4 , dimension(:), 
     $     allocatable, target :: idtype
      integer*8 , dimension(:), 
     $     allocatable, target :: idval
c$$$      character*(MAXPNAME), pointer, dimension(:) :: ppname
c$$$      integer*4 , pointer, dimension(:) :: plpname(:)
c$$$      integer*4 , pointer, dimension(:) :: pidtype(:)
c$$$      integer*8 , pointer, dimension(:) :: pidval(:)
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
        character*(*) , intent(in) :: token
        integer*8, target, intent(in):: ival
        integer*4 , intent(in) :: type
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
        character*(*) , intent(in) ::token
        integer*4 , intent(in) :: type,ival
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
      integer*4, parameter:: maxgeneration=2**30-1,maxlevele=2**14,
     $     nsymhash=2047,nslots=32,maxlbuf=2**22
      real*8 dinfinity,dnotanumber
      integer*8
     $     itfcontroot,itfcontext,itfcontextpath,itflocal,
     $     kxeof,kxfailed,iaxhold,iaximmediate,
     $     iaxline,kxliteral,kxnull,kxnulll,kxnulls,
     $     iaxout,iaxpriority,iaxschar,
     $     iaxslotnull,iaxslotpart,iaxslotseqnull,
     $     kxvect,kxvect1,iavpw,
     $     kerror,ierrorf,ierrorgen,ierrorprint,ierrorth,
     $     ierrorexp,ifunbase,initmessage,levelcompile,
     $     kinfinity,kminfinity,knotanumber
      integer*4 
     $     levele,levelp,lgeneration,ltrace,
     $     modethrow,iordless
      end module

      module tfmem
      implicit none
      integer*8, parameter :: mpsize=2**22,kcpklist0=0,maxstack=2**25,
     $     minstack=2**18
      integer*4, parameter :: nindex=64,mhash=32767,
     $     minseg0=9,minseg1=16,minseg2=16
      integer*4, parameter :: ncbk = 2**16
      integer*8, parameter :: kcpoffset = 0
      type cbkalloc
        integer*8, allocatable :: ca(:)
      end type

      type (cbkalloc), target, allocatable :: sadalloc(:)
      integer*8 ,allocatable :: kcbk(:,:)
      integer*8 :: icp=0,nmem=0,nnet=0,ich=0,maxic=3,
     $     minic,icsep,nitaloc
      integer*8 kfirstalloc
      integer*4 :: icbk=1,jcbk=1

      type sad_descriptor
        sequence
        real*8 x(1:0)
        integer*8 k
      end type

      interface sad_loc
        module procedure ksad_loc,dsad_loc,isad_loc,rsad_loc
      end interface

      contains
        integer*8 function dsad_loc(k)
        use iso_c_binding
        implicit none
        type (sad_descriptor), target, intent(in) :: k
        dsad_loc=(transfer(c_loc(k),int8(0))-kcpklist0)/8
        return
        end function

        integer*8 function ksad_loc(k)
        use iso_c_binding
        implicit none
        integer*8, target, intent(in) :: k
        ksad_loc=(transfer(c_loc(k),int8(0))-kcpklist0)/8
        return
        end function

        integer*8 function isad_loc(i)
        use iso_c_binding
        implicit none
        integer*4, target, intent(in):: i
        isad_loc=(transfer(c_loc(i),int8(0))-kcpklist0)/8
        return
        end function

        integer*8 function rsad_loc(x)
        use iso_c_binding
        implicit none
        real*8, target, intent(in):: x
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
        allocate(pname(0:HTMAX))
        allocate(lpname(0:HTMAX))
        allocate(idtype(0:HTMAX))
        allocate(idval(0:HTMAX))
c$$$        ppname=>pname
c$$$        plpname=>lpname
c$$$        pidtype=>idtype
c$$$        pidval=>idval
        call c_f_pointer(transfer(kcpklist0+8,cp),klist,[klistlen])
        call lminit(klist(0),lps,pname(0),lpname(0),idtype(0),idval(0))
        call c_f_pointer(c_loc(klist(1)),rlist,[klistlen])
        call c_f_pointer(c_loc(klist(1)),ilist,[2,klistlen])
        call c_f_pointer(c_loc(klist(1)),jlist,[8,klistlen])
        return
        end subroutine

        subroutine talocinit
        use maccbk
        use iso_c_binding, only:c_loc
        implicit none
        integer*8 ka,ic
        allocate(sadalloc(ncbk))
        allocate(sadalloc(1)%ca(nindex*2+mhash+16))
        ka=transfer(c_loc(sadalloc(1)%ca(1)),int8(0))
        kfirstalloc=ka
c     kcpklist0=0
        call tfcbkinit
        icp=ksad_loc(sadalloc(1)%ca(1))
        allocate (kcbk(3,ncbk))
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
        integer*4 , intent(in):: n
        integer*4 m,n1,m1
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
          elseif(m1-minseg2 .ge. m)then
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
        end function

        subroutine tsetindexhash(ip,m)
        use maccbk
        implicit none
        integer*4 , intent(in) :: m
        integer*8 , intent(in) ::ip
        integer*8 ic,ic1,ia
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
        end subroutine

        subroutine tfree(ka)
        use maccbk
c        use ffs0, only:ifinext
        implicit none
        integer*8, intent(in):: ka
        integer*8 ix,ik,ik0,ip,ix1
        integer*4 m,mx
        m=ilist(1,ka-1)
c        if(associated(ifinext) .and. ka .eq. ifinext)then
c          write(*,*)'tfree ',ka,m
c        endif
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
        end subroutine
      
        integer*8 function ktzaloc(ktype,nw)
        use maccbk
        implicit none
        integer*4, intent(in):: nw
        integer*8, intent(in):: ktype
        integer*8 k1
        k1=ktaloc(nw+2)
        ilist(2,k1-1)=0
        klist(k1)=ktype
        ilist(1,k1+1)=0
        ktzaloc=k1+2
        return
        end function

      end module

      module tfcode
      use tfmem, only:sad_descriptor
      real*8, parameter :: xinfinity=1.7976931348623157D308
      integer*4, parameter :: ntfoper=0,ntfreal=1,ntflist=3,ntflistr=4,
     $     ntfstkseq=5,
     $     ntfstk=6,ntfstring=101,ntfsymbol=201,ntfpat=203,ntfarg=205,
     $     ntffun=ntfoper,ntfdef=ntfsymbol
      integer*4, parameter :: nfunif=39,nfunlength=20,nfunreppart=22,
     $     nfundo=29,
     $     nfunmodule=34,nfunblock=35,nfunswicases=37,
     $     nfunselect=41,nfunappend=44,
     $     nfunprepend=45,nfunposition=59,nfunthread=83,nfunscan=96,
     $     nfununeval=132,nfuncases=133,nfundelcases=134,
     $     nfunwith=141,nfunselcases=142,nfunextract=154,
     $     nfuninsert=177,nfundelete=178
      integer*4, parameter ::
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
     $     mtflist=mtfleftbrace,mtftimes=mtfmult
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
      integer*4, parameter :: irtcret=-4,irtcthrow=-5,
     $     irtcgoto=-6,irtcabort=-7
      integer*8 ktfoper,ktflist,ktfstring,ktfsymbol,ktfpat,ktfobj,
     $     ktfmask,ktamask,ktrmask,ktfnull,ktfnr,ktfref,ktfother,
     $     ktomask,ktftrue,ktfnan,ktfenan
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
     $     ktfnan   =int8(z'fff8000000000000'),
     $     ktfenan  =int8(z'7ff0000000000000'),
     $     ktfenanb =int8(z'000fffffffffffff')
     $     )
      integer*4 , parameter :: mbody = 2**12

      type sad_object
      sequence
      integer*4 len,attr
      type (sad_descriptor) alloc
      integer*4 ref,nl
      integer*8 body(1:0)
      type (sad_descriptor) dbody(0:mbody)
      end type

      type sad_list
      sequence
      integer*2 lenp,lena
      integer*4 attr
      integer*8 alloc
      integer*4 ref,nl
      type (sad_descriptor) head
      real*8 rbody(1:0)
      complex*16 cbody(1:0)
      type (sad_descriptor) dbody(1:0)
      integer*8 body(1:mbody)
      end type

      type sad_dlist
      sequence
      type (sad_list) list(1:0)
      integer*2 lenp,lena
      integer*4 attr
      type (sad_descriptor) alloc
      integer*4 ref,nl
      type (sad_descriptor) head
      real*8 rbody(1:0)
      complex*16 cbody(1:0)
      integer*8 body(1:0)
      type (sad_descriptor) dbody(1:mbody)
      end type

      type sad_rlist
      sequence
      type (sad_dlist) dlist(1:0)
      integer*2 lenp,lena
      integer*4 attr
      type (sad_descriptor) alloc
      integer*4 ref,nl
      type (sad_descriptor) head
      complex*16 cbody(1:0)
      integer*8 body(1:0)
      type (sad_descriptor) dbody(1:0)
      real*8 rbody(1:mbody)
      end type

      type sad_complex
      sequence
      integer*2 lenp,lena
      integer*4 attr
      integer*8 alloc
      integer*4 ref,nl
      type (sad_descriptor) head
      integer*8 body(1:0)
      type (sad_descriptor) dbody(1:0)
      complex*16 cx(1:0)
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
      type (sad_descriptor) alloc
      integer*4 ref,gen
      integer*4 nch,nc
      integer*1 istr(1:0)
      integer*8 kstr(1:0)
c size limitation due to gfortran 7 on macOS ???
      character*(mbody) str
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
      integer*8 hash(0:-1)
      type (sad_descriptor) dhash(0:2**10-1)
      end type

      end module

      module tfstk
      use tfcbk
      use tfcode
      use maccbk
      use tfmem, only:sad_loc,ksad_loc,ktaloc,tfree
      public
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
        module procedure loc_sym,loc_string,
     $     loc_pat,loc_obj,loc_complex,loc_symdef,
     $     loc_dlist,loc_rlist
      end interface

      interface descr_sad
        module procedure descr_sym,descr_string,
     $     descr_pat,descr_obj,descr_complex,descr_symdef,
     $     descr_dlist,descr_rlist
      end interface

      interface ktfaddr
        module procedure ktfaddrk,ktfaddrd
      end interface

      interface kxaaloc
        module procedure kxaaloc_dlist,kxaaloc_rlist
      end interface

      interface ktaaloc
        module procedure ktaaloc_dlist,ktaaloc_rlist
      end interface

      interface ktfreallistqo
        module procedure ktfreallistqo_rlist,
     $     ktfreallistqo_dlist
      end interface

      interface ktfreallistq
        module procedure ktfreallistqo_rlist,
     $     ktfreallistqo_dlist,ktfreallistqk,ktfreallistqd,
     $     ktfreallistqk_rlist
      end interface

      interface ktfnonreallistqo
        module procedure ktfnonreallistqo_rlist,
     $     ktfnonreallistqo_dlist
      end interface

      interface kxadaloc
        module procedure kxadaloc_dlist
      end interface

      interface ktadaloc
        module procedure ktadaloc_dlist
      end interface

      interface ktflistq
        module procedure ktflistq_dlist,ktflistqd_rlist,
     $     ktflistqd_dlist,ktflistq_rlist
      end interface

      interface ktfnonlistq
        module procedure ktfnonlistq_dlist,ktfnonlistq_rlist,
     $     ktfnonlistqd_dlist,ktfnonlistqd_rlist
      end interface

      interface tflistq
        module procedure tflistq_dlist,
     $     tflistqd_dlist,tflistq_rlist
      end interface

      interface tfnonlistq
        module procedure tfnonlistq_dlist,
     $     tfnonlistqd_dlist,tfnonlistq_rlist,tfnonlistqd_rlist
      end interface

      interface ktfrealq
        module procedure ktfrealq_k,ktfrealq_d,ktfrealq_ki,ktfrealq_di
      end interface

      interface ktfnonrealq
        module procedure ktfnonrealq_k,ktfnonrealq_d,
     $     ktfnonrealq_ki,ktfnonrealq_di
      end interface

      interface tfruleq
        module procedure tfruleqk_dlist,tfruleqd_dlist
      end interface

      interface tfnumlistqn
        module procedure tfnumlistqnk,tfnumlistqnd
      end interface

      interface tfcomplexq
        module procedure tfcomplexqc,tfcomplexqx,
     $     tfcomplexqdx,tfcomplexqdc
      end interface

      interface tfreallistq
        module procedure tfreallistqk,tfreallistqd
      end interface

      interface tfnonreallistq
        module procedure tfnonreallistqk,tfnonreallistqd
      end interface

      interface tfgetllstkall
        module procedure tfgetllstkall_dlist,
     $     tfgetllstkall_rlist
      end interface

      interface tfgetdefargp
        module procedure tfgetdefargp_dlist
      end interface

      interface ktfstringq
        module procedure ktfstringqk,ktfstringqd
      end interface

      interface kxmakelist
        module procedure kxmakelist_dlist,kxmakelist_rlist
      end interface

      interface sad_descr
        module procedure list_descr,rlist_descr,dlist_descr,
     $     symbol_descr,dfromr,dfromk,string_descr,pat_descr
      end interface

      interface ktfpatq
        module procedure ktfpatqk,ktfpatqd
      end interface

      interface ktfsymbolq
        module procedure ktfsymbolqk,ktfsymbolqd
      end interface

      interface ktfsequenceq
        module procedure ktfsequenceq_dlist,ktfsequenceq_rlist,
     $     ktfsequenceqd_dlist,ktfsequenceqd_rlist
      end interface

      interface tfcomplexnumlistqk
        module procedure 
     $     tfcomplexnumlistqk_dlist
      end interface

      interface ktfmakelist
        moduleprocedure ktfmakelist_dlist,ktfmakelist_rlist
      end interface

      interface ktfobjq
        module procedure ktfobjqk,ktfobjqd
      end interface

      interface tfmakerulestk
        module procedure tfmakerulestk_dd,tfmakerulestk_dr
      end interface

      interface ktfnonoperq
        module procedure ktfnonoperqk,ktfnonoperqd
      end interface

      interface ktfoperq
        module procedure ktfoperqk,ktfoperqd
      end interface

      interface ktfrefq
        module procedure ktfrefqk,ktfrefqd
      end interface

      interface tfnumberq
        module procedure tfnumberqk,tfnumberqd
      end interface

      interface ktfnonsymbolq
        module procedure ktfnonsymbolqk,ktfnonsymbolqd
      end interface

      interface tfsameq
        module procedure tfsameqk,tfsameqd
      end interface

      interface tfsamesymbolq
        module procedure tfsamesymbolqk,tfsamesymbolqd,tfsamesymbolqo
      end interface

      interface tfsamestringq
        module procedure tfsamestringqk,tfsamestringqd,tfsamestringqo
      end interface

      interface tfconstq
        module procedure tfconstqk,tfconstqd
      end interface

      interface tfconstpatternq
        module procedure tfconstpatternqk,tfconstpatternqd
      end interface

      interface tfsameheadq
        module procedure tfsameheadqk,tfsameheadqd
      end interface

      interface tfexprq
        module procedure tfexprqk,tfexprqd
      end interface

      interface tfinequalityq
        module procedure tfinequalityqk,tfinequalityqd
      end interface

      contains
        subroutine tfinitstk
        use iso_c_binding
        use tfmem, only:maxstack,minstack
        implicit none
        integer*4 idummy
        if(tfstkinit)then
          return
        endif
c        mstk=max(2**18,int(rgetgl1('STACKSIZ')))
        mstk=int(maxstack*2)
        ispbase=0
        do while (ispbase .le. 0)
          mstk=mstk/2
          if(mstk .lt. minstack)then
            write(*,*)'Stack allocation failed: ',mstk,ispbase
            call abort
          endif
          ispbase=ktaloc(mstk*2)-1
        enddo
        call rsetGL('STACKSIZ',dble(mstk/2),idummy)
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
        integer*8 , intent(in) :: n
        integer*4 istat,i
        ktfsadalloc=0
        do i=icbk+1,ncbk
          allocate(sadalloc(i)%ca(n),stat=istat)
          if(istat .ne. 0)then
c            write(*,*)'ktfsadalloc allocation error in ALLOCATE: ',
c     $           n,i,istat
            ktfsadalloc=-1
            return
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
        integer*8 , intent(in) :: ka,n
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
                  kcbk(:,j)=0
                  if(j .eq. jcbk)then
                    jcbk=jcbk-1
                  endif
                else
                  kcbk(2,j)=kcbk(2,k)
                  kcbk(3,j)=kcbk(2,k)
                  kcbk(:,k)=0
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
                  kcbk(:,j)=0
                  if(j .eq. jcbk)then
                    jcbk=jcbk-1
                  endif
                else
                  kcbk(1,j)=kcbk(1,k)
                  kcbk(:,k)=0
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
        integer*8 , intent(in) ::k
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
        integer*8 , intent(in) ::k
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
        integer*8, intent(in):: k
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
        integer*8 , intent(in):: ka
        call c_f_pointer(c_loc(klist(klist(ifunbase+ka)-9)),fun)
        iget_fun_id=fun%id
        return
        end function

        subroutine dlist_list(dlist,list)
        use iso_c_binding
        implicit none
        type (sad_dlist), target, intent(in) :: dlist
        type (sad_list), pointer, intent(out) :: list
        call c_f_pointer(c_loc(dlist),list)
        return
        end subroutine

        subroutine dlist_rlist(dlist,rlist)
        use iso_c_binding
        implicit none
        type (sad_dlist), target, intent(in) :: dlist
        type (sad_rlist), pointer, intent(out) :: rlist
        call c_f_pointer(c_loc(dlist),rlist)
        return
        end subroutine

        subroutine list_dlist(list,dlist)
        use iso_c_binding
        implicit none
        type (sad_list), target, intent(in) :: list
        type (sad_dlist), pointer, intent(out) :: dlist
        call c_f_pointer(c_loc(list),dlist)
        return
        end subroutine

        subroutine loc_list(locp,list)
        use iso_c_binding
        implicit none
        type (sad_list), pointer, intent(out) :: list
        integer*8 , intent(in) :: locp
        call c_f_pointer(c_loc(klist(locp-3)),list)
        return
        end subroutine

        subroutine loc_dlist(locp,dlist)
        use iso_c_binding
        implicit none
        type (sad_dlist), pointer, intent(out) :: dlist
        integer*8 , intent(in) :: locp
        call c_f_pointer(c_loc(klist(locp-3)),dlist)
        return
        end subroutine

        subroutine loc_rlist(locp,rlist)
        use iso_c_binding
        implicit none
        type (sad_rlist), pointer, intent(out) :: rlist
        integer*8 , intent(in) :: locp
        call c_f_pointer(c_loc(klist(locp-3)),rlist)
        return
        end subroutine

        subroutine descr_list(dscr,list)
        implicit none
        type (sad_descriptor) , intent(in)::dscr
        type (sad_list), pointer, intent(out) :: list
        call loc_list(ktfaddrd(dscr),list)
        return
        end subroutine

        subroutine descr_rlist(dscr,list)
        implicit none
        type (sad_descriptor) , intent(in)::dscr
        type (sad_rlist), pointer, intent(out) :: list
        call loc_rlist(ktfaddrd(dscr),list)
        return
        end subroutine

        subroutine descr_dlist(dscr,dlist)
        implicit none
        type (sad_descriptor) , intent(in)::dscr
        type (sad_dlist), pointer, intent(out) :: dlist
        call loc_dlist(ktfaddrd(dscr),dlist)
        return
        end subroutine

        subroutine loc_complex(locp,cx)
        use iso_c_binding
        implicit none
        type (sad_complex), pointer, intent(out) :: cx
        integer*8 , intent(in)::locp
        call c_f_pointer(c_loc(klist(locp-3)),cx)
        return
        end subroutine

        subroutine descr_complex(dscr,complex)
        implicit none
        type (sad_descriptor) , intent(in)::dscr
        type (sad_complex), pointer, intent(out) :: complex
        call loc_complex(ktfaddrd(dscr),complex)
        return
        end subroutine

        subroutine loc_obj(locp,obj)
        use iso_c_binding
        implicit none
        type (sad_object), pointer, intent(out) :: obj
        integer*8 , intent(in)::locp
        call c_f_pointer(c_loc(klist(locp-3)),obj)
        return
        end subroutine

        subroutine descr_obj(dscr,obj)
        implicit none
        type (sad_descriptor) , intent(in)::dscr
        type (sad_object), pointer, intent(out) :: obj
        call loc_obj(ktfaddrd(dscr),obj)
        return
        end subroutine

        subroutine loc_namtbl(locp,loc)
        use iso_c_binding
        implicit none
        type (sad_namtbl), pointer, intent(out) :: loc
        integer*8 , intent(in)::locp
        call c_f_pointer(c_loc(klist(locp-1)),loc)
        return
        end subroutine

        subroutine descr_namtbl(dscr,namtbl)
        implicit none
        type (sad_descriptor) , intent(in):: dscr
        type (sad_namtbl), pointer, intent(out) :: namtbl
        call loc_namtbl(ktfaddrd(dscr),namtbl)
        return
        end subroutine

        subroutine sym_namtbl(sym,loc)
        implicit none
        type (sad_symbol) , intent(in)::sym
        type (sad_namtbl), pointer, intent(out) :: loc
        call loc_namtbl(sym%loc,loc)
        return
        end subroutine

        subroutine sym_symstr(sym,str)
        implicit none
        type (sad_symbol) , intent(in)::sym
        type (sad_string), pointer, intent(out) :: str
        call loc_symstr(sym%loc,str)
        return
        end subroutine

        subroutine loc_symstr(loc,str)
        implicit none
        integer*8 , intent(in)::loc
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
        integer*8 , intent(in)::locp
        call c_f_pointer(c_loc(klist(locp-3)),pat)
        return
        end subroutine

        subroutine descr_pat(dscr,pat)
        implicit none
        type (sad_descriptor) , intent(in)::dscr
        type (sad_pat), pointer, intent(out) :: pat
        call loc_pat(ktfaddrd(dscr),pat)
        return
        end subroutine

        subroutine loc_string(locp,str)
        use iso_c_binding
        implicit none
        type (sad_string), pointer, intent(out) :: str
        integer*8 , intent(in)::locp
        call c_f_pointer(c_loc(klist(locp-3)),str)
        return
        end subroutine

        subroutine descr_string(dscr,string)
        implicit none
        type (sad_descriptor) , intent(in)::dscr
        type (sad_string), pointer, intent(out) :: string
        call loc_string(ktfaddrd(dscr),string)
        return
        end subroutine

        subroutine loc_sym(locp,sym)
        use iso_c_binding
        implicit none
        type (sad_symbol), pointer, intent(out) :: sym
        integer*8 , intent(in)::locp
        call c_f_pointer(c_loc(klist(locp-3)),sym)
        return
        end subroutine

        subroutine descr_sym(dscr,symbol)
        implicit none
        type (sad_descriptor) , intent(in)::dscr
        type (sad_symbol), pointer, intent(out) :: symbol
        call loc_sym(ktfaddrd(dscr),symbol)
        return
        end subroutine

        subroutine loc1_symdef(locp1,symd)
        use iso_c_binding
        implicit none
        type (sad_symdef), pointer, intent(out) :: symd
        integer*8 , intent(in)::locp1
        call c_f_pointer(c_loc(klist(locp1-1)),symd)
        return
        end subroutine

        subroutine loc_symdef(locp,symd)
        use iso_c_binding
        implicit none
        type (sad_symdef), pointer, intent(out) :: symd
        integer*8 , intent(in)::locp
        call c_f_pointer(c_loc(klist(locp-9)),symd)
        return
        end subroutine

        subroutine descr_symdef(dscr,symdef)
        implicit none
        type (sad_descriptor) , intent(in)::dscr
        type (sad_symdef), pointer, intent(out) :: symdef
        call loc_symdef(ktfaddrd(dscr),symdef)
        return
        end subroutine

        subroutine sym_symdef(sym,symd)
        use iso_c_binding
        implicit none
        type (sad_symbol), target, intent(in) :: sym
        type (sad_symdef), pointer, intent(out) :: symd
        call c_f_pointer(c_loc(sym%dummy(1)),symd)
        return
        end subroutine

        subroutine loc_deftbl(locp,loc)
        use iso_c_binding
        implicit none
        type (sad_deftbl), pointer, intent(out) :: loc
        integer*8 , intent(in)::locp
        call c_f_pointer(c_loc(klist(locp-1)),loc)
        return
        end subroutine

        subroutine loc_defhash(locp,loc)
        use iso_c_binding
        implicit none
        type (sad_defhash), pointer, intent(out) :: loc
        integer*8 , intent(in)::locp
        call c_f_pointer(c_loc(klist(locp-1)),loc)
        return
        end subroutine

        subroutine dlist_dlist(kl,kl1)
        implicit none
        type (sad_dlist) , target, intent(in) :: kl
        type (sad_dlist), pointer, intent(out) :: kl1
        kl1=>kl
        return
        end subroutine dlist_dlist

        type (sad_descriptor) function dlist_descr(kl)
        implicit none
        type (sad_dlist) , intent(in)::kl
        dlist_descr%k=ktflist+sad_loc(kl%head)
        return
        end

        type (sad_descriptor) function list_descr(kl)
        implicit none
        type (sad_list) , intent(in)::kl
        list_descr%k=ktflist+sad_loc(kl%head)
        return
        end

        type (sad_descriptor) function rlist_descr(kl)
        implicit none
        type (sad_rlist) , intent(in)::kl
        rlist_descr%k=ktflist+sad_loc(kl%head)
        return
        end

        type (sad_descriptor) function symbol_descr(s)
        implicit none
        type (sad_symbol) , intent(in)::s
        symbol_descr%k=ktfsymbol+sad_loc(s%loc)
        return
        end

        type (sad_descriptor) function string_descr(s)
        implicit none
        type (sad_string) , intent(in)::s
        string_descr%k=ktfstring+sad_loc(s%nch)
        return
        end

        type (sad_descriptor) function pat_descr(p)
        implicit none
        type (sad_pat) , intent(in)::p
        pat_descr%k=ktfpat+sad_loc(p%len)+3
        return
        end

        integer*8 function ktfaddrk(k)
        implicit none
        integer*8 , intent(in)::k
        ktfaddrk=iand(ktamask,k)
        return
        end function ktfaddrk

        integer*8 function ktfaddrd(k)
        implicit none
        type (sad_descriptor) , intent(in)::k
        ktfaddrd=iand(ktamask,k%k)
        return
        end function ktfaddrd

        integer*8 function ktftype(k)
        implicit none
        integer*8 , intent(in)::k
        ktftype=iand(ktfmask,k)
        return
        end function ktftype

        logical*4 function ktfobjqk(k,obj)
        implicit none
        type (sad_object), pointer, optional, intent(out) :: obj
        integer*8 , intent(in)::k
        if(iand(ktomask,k) .eq. ktfobj)then
          ktfobjqk=.true.
          if(present(obj))then
            call loc_obj(iand(ktamask,k),obj)
          endif
        else
          ktfobjqk=.false.
        endif
        return
        end function ktfobjqk

        logical*4 function ktfobjqd(k,obj)
        implicit none
        type (sad_descriptor) , intent(in)::k
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
        integer*8 , intent(in)::ka
        ktfnonobjq=iand(ktomask,ka) .ne. ktfobj
        return
        end function ktfnonobjq

        logical*4 function ktfnanq(x)
        implicit none
        real*8 , intent(in)::x
        integer*8 kfromr
        ktfnanq=kfromr(x) .eq. ktfnan
        return
        end

        logical*4 function ktfenanq(x)
        implicit none
        real*8 , intent(in)::x
        integer*8 kfromr,k
        k=kfromr(x)
        ktfenanq=k .eq. knotanumber .or.
     $       iand(k,ktfenan) .eq. ktfenan .and.
     $       iand(k,ktfenanb) .ne. 0 .and. k .ne. kinfinity .and.
     $       k .ne. kminfinity
        return
        end

        logical*4 function ktfrealq_k(k,v)
        implicit none
        integer*8 , intent(in)::k
        real*8, optional, intent(out) :: v
        real*8 rfromk
        ktfrealq_k=iand(ktrmask,k) .ne. ktfnr
        if(ktfrealq_k .and. present(v))then
          v=rfromk(k)
        endif
        return
        end function ktfrealq_k

        logical*4 function ktfrealq_d(k,v)
        implicit none
        type (sad_descriptor) , intent(in)::k
        real*8, optional, intent(out) :: v
        real*8 rfromk
        ktfrealq_d=iand(ktrmask,k%k) .ne. ktfnr
        if(ktfrealq_d .and. present(v))then
          v=rfromk(k%k)
        endif
        return
        end function ktfrealq_d

        logical*4 function ktfrealq_ki(k,iv)
        implicit none
        integer*8 , intent(in)::k
        integer*4 , intent(out) :: iv
        real*8 rfromk
        ktfrealq_ki=iand(ktrmask,k) .ne. ktfnr
        if(ktfrealq_ki)then
          iv=int(rfromk(k))
        endif
        return
        end function ktfrealq_ki

        logical*4 function ktfrealq_di(k,iv)
        implicit none
        type (sad_descriptor) , intent(in)::k
        integer*4 , intent(out) :: iv
        real*8 rfromk
        ktfrealq_di=iand(ktrmask,k%k) .ne. ktfnr
        if(ktfrealq_di)then
          iv=int(rfromk(k%k))
        endif
        return
        end function ktfrealq_di

        logical*4 function ktfnonrealq_k(k,v)
        implicit none
        integer*8 , intent(in)::k
        real*8, optional, intent(out) :: v
        real*8 rfromk
        ktfnonrealq_k=iand(ktrmask,k) .eq. ktfnr
        if(.not. ktfnonrealq_k .and. present(v))then
          v=rfromk(k)
        endif
        return
        end function ktfnonrealq_k

        logical*4 function ktfnonrealq_d(k,v)
        implicit none
        type (sad_descriptor) , intent(in)::k
        real*8, optional, intent(out) :: v
        real*8 rfromk
        ktfnonrealq_d=iand(ktrmask,k%k) .eq. ktfnr
        if(.not. ktfnonrealq_d .and. present(v))then
          v=rfromk(k%k)
        endif
        return
        end function ktfnonrealq_d

        logical*4 function ktfnonrealq_ki(k,iv)
        implicit none
        integer*8 , intent(in)::k
        integer*4 , intent(out) :: iv
        real*8 rfromk
        ktfnonrealq_ki=iand(ktrmask,k) .eq. ktfnr
        if(.not. ktfnonrealq_ki)then
          iv=int(rfromk(k))
        endif
        return
        end function ktfnonrealq_ki

        logical*4 function ktfnonrealq_di(k,iv)
        implicit none
        type (sad_descriptor) , intent(in)::k
        integer*4 , intent(out) :: iv
        real*8 rfromk
        ktfnonrealq_di=iand(ktrmask,k%k) .eq. ktfnr
        if(.not. ktfnonrealq_di)then
          iv=int(rfromk(k%k))
        endif
        return
        end function ktfnonrealq_di

        logical*4 function ktfoperqk(k,ka)
        implicit none
        integer*8 , intent(in)::k
        integer*8 , optional, intent(out) :: ka
        ktfoperqk=iand(ktfmask,k) .eq. ktfoper
        if(ktfoperqk .and. present(ka))then
          ka=ktfaddr(k)
        endif
        return
        end function ktfoperqk

        logical*4 function ktfnonoperqk(k,ka)
        implicit none
        integer*8 , intent(in)::k
        integer*8 , optional, intent(out) :: ka
        ktfnonoperqk=iand(ktfmask,k) .ne. ktfoper
        if(.not. ktfnonoperqk .and. present(ka))then
          ka=ktfaddr(k)
        endif
        return
        end function ktfnonoperqk

        logical*4 function ktfoperqd(k,ka)
        implicit none
        type (sad_descriptor) , intent(in)::k
        integer*8, optional, intent(out) :: ka
        ktfoperqd=iand(ktfmask,k%k) .eq. ktfoper
        if(ktfoperqd .and. present(ka))then
          ka=ktfaddr(k)
        endif
        return
        end function ktfoperqd

        logical*4 function ktfnonoperqd(k,ka)
        implicit none
        type (sad_descriptor) , intent(in)::k
        integer*8, optional, intent(out) :: ka
        ktfnonoperqd=iand(ktfmask,k%k) .ne. ktfoper
        if(.not. ktfnonoperqd .and. present(ka))then
          ka=ktfaddr(k)
        endif
        return
        end function ktfnonoperqd

        logical*4 function ktfstringqk(k,str)
        implicit none
        type (sad_string), pointer, optional, intent(out) :: str
        integer*8 , intent(in)::k
        ktfstringqk=iand(ktfmask,k) .eq. ktfstring
        if(present(str) .and. ktfstringqk)then
          call loc_string(ktfaddr(k),str)
        endif
        return
        end function ktfstringqk

        logical*4 function ktfstringqd(k,str)
        implicit none
        type (sad_descriptor) , intent(in)::k
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
        integer*8 , intent(in)::k
        ktfnonstringq=iand(ktfmask,k) .ne. ktfstring
        if(present(str) .and. .not. ktfnonstringq)then
          call loc_string(ktfaddr(k),str)
        endif
        return
        end function ktfnonstringq

        logical*4 function ktflistq_rlist(k,kl)
        implicit none
        type (sad_rlist), pointer, intent(out) :: kl
        integer*8 , intent(in)::k
        if(iand(ktfmask,k) .eq. ktflist)then
          ktflistq_rlist=.true.
          call loc_sad(iand(ktamask,k),kl)
        else
          ktflistq_rlist=.false.
        endif          
        return
        end function ktflistq_rlist

        logical*4 function ktflistq_dlist(k,kl)
        implicit none
        type (sad_dlist), pointer, optional, intent(out) :: kl
        integer*8 , intent(in)::k
        if(iand(ktfmask,k) .eq. ktflist)then
          ktflistq_dlist=.true.
          if(present(kl))then
            call loc_dlist(iand(ktamask,k),kl)
          endif
        else
          ktflistq_dlist=.false.
        endif          
        return
        end function ktflistq_dlist

        logical*4 function ktflistqd_rlist(k,kl)
        implicit none
        type (sad_rlist), pointer, intent(out) :: kl
        type (sad_descriptor), intent(in):: k
        if(iand(ktfmask,k%k) .eq. ktflist)then
          ktflistqd_rlist=.true.
          call loc_sad(iand(ktamask,k%k),kl)
        else
          ktflistqd_rlist=.false.
        endif          
        return
        end function ktflistqd_rlist

        logical*4 function ktflistqd_dlist(k,kl)
        implicit none
        type (sad_dlist), pointer, optional, intent(out) :: kl
        type (sad_descriptor), intent(in):: k
        if(iand(ktfmask,k%k) .eq. ktflist)then
          ktflistqd_dlist=.true.
          if(present(kl))then
            call loc_dlist(iand(ktamask,k%k),kl)
          endif
        else
          ktflistqd_dlist=.false.
        endif          
        return
        end function ktflistqd_dlist

        logical*4 function ktfnonlistq_rlist(k,kl)
        implicit none
        type (sad_rlist), pointer,  intent(out) :: kl
        integer*8 k
        if(iand(ktfmask,k) .eq. ktflist)then
          ktfnonlistq_rlist=.false.
          call loc_sad(iand(ktamask,k),kl)
        else
          ktfnonlistq_rlist=.true.
        endif          
        return
        end function ktfnonlistq_rlist

        logical*4 function ktfnonlistq_dlist(k,kl)
        implicit none
        type (sad_dlist), pointer, optional, intent(out) :: kl
        integer*8 , intent(in)::k
        if(iand(ktfmask,k) .eq. ktflist)then
          ktfnonlistq_dlist=.false.
          if(present(kl))then
            call loc_dlist(iand(ktamask,k),kl)
          endif
        else
          ktfnonlistq_dlist=.true.
        endif          
        return
        end function ktfnonlistq_dlist

        logical*4 function ktfnonlistqd_rlist(k,kl)
        implicit none
        type (sad_rlist), pointer, intent(out) :: kl
        type (sad_descriptor), intent(in):: k
        if(iand(ktfmask,k%k) .eq. ktflist)then
          ktfnonlistqd_rlist=.false.
          call loc_sad(iand(ktamask,k%k),kl)
        else
          ktfnonlistqd_rlist=.true.
        endif          
        return
        end function ktfnonlistqd_rlist

        logical*4 function ktfnonlistqd_dlist(k,kl)
        implicit none
        type (sad_dlist), pointer, optional, intent(out) :: kl
        type (sad_descriptor), intent(in):: k
        if(iand(ktfmask,k%k) .eq. ktflist)then
          ktfnonlistqd_dlist=.false.
          if(present(kl))then
            call loc_dlist(iand(ktamask,k%k),kl)
          endif
        else
          ktfnonlistqd_dlist=.true.
        endif          
        return
        end function ktfnonlistqd_dlist

        logical*4 function tflistq_rlist(k,kl)
        implicit none
        type (sad_rlist), pointer, optional, intent(out) :: kl
        integer*8 , intent(in)::k
        if(iand(ktfmask,k) .eq. ktflist .and.
     $       klist(ktfaddr(k)) .eq. ktfoper+mtflist)then
          tflistq_rlist=.true.
          if(present(kl))then
            call loc_sad(iand(ktamask,k),kl)
          endif
        else
          tflistq_rlist=.false.
        endif          
        return
        end function tflistq_rlist

        logical*4 function tflistq_dlist(k,kl)
        implicit none
        type (sad_dlist), pointer, intent(out) :: kl
        integer*8 , intent(in)::k
        if(iand(ktfmask,k) .eq. ktflist .and.
     $       klist(ktfaddr(k)) .eq. ktfoper+mtflist)then
          tflistq_dlist=.true.
          call loc_dlist(iand(ktamask,k),kl)
        else
          tflistq_dlist=.false.
        endif          
        return
        end function tflistq_dlist

        logical*4 function tflistqd_rlist(k,kl)
        implicit none
        type (sad_rlist), pointer, intent(out) :: kl
        type (sad_descriptor), intent(in):: k
        if(iand(ktfmask,k%k) .eq. ktflist .and.
     $       klist(ktfaddr(k%k)) .eq. ktfoper+mtflist)then
          tflistqd_rlist=.true.
          call loc_sad(iand(ktamask,k%k),kl)
        else
          tflistqd_rlist=.false.
        endif          
        return
        end function tflistqd_rlist

        logical*4 function tflistqd_dlist(k,kl)
        implicit none
        type (sad_dlist), pointer, optional, intent(out) :: kl
        type (sad_descriptor), intent(in):: k
        if(iand(ktfmask,k%k) .eq. ktflist .and.
     $       klist(ktfaddr(k%k)) .eq. ktfoper+mtflist)then
          tflistqd_dlist=.true.
          if(present(kl))then
            call loc_dlist(iand(ktamask,k%k),kl)
          endif
        else
          tflistqd_dlist=.false.
        endif          
        return
        end function tflistqd_dlist

        logical*4 function tfnonlistq_rlist(k,kl)
        implicit none
        type (sad_rlist), pointer, optional, intent(out) :: kl
        integer*8 , intent(in)::k
        if(iand(ktfmask,k) .eq. ktflist .and.
     $       klist(ktfaddr(k)) .eq. ktfoper+mtflist)then
          tfnonlistq_rlist=.false.
          if(present(kl))then
            call loc_sad(iand(ktamask,k),kl)
          endif
        else
          tfnonlistq_rlist=.true.
        endif          
        return
        end function tfnonlistq_rlist

        logical*4 function tfnonlistq_dlist(k,kl)
        implicit none
        type (sad_dlist), pointer, intent(out) :: kl
        integer*8 , intent(in)::k
        if(iand(ktfmask,k) .eq. ktflist .and.
     $       klist(ktfaddr(k)) .eq. ktfoper+mtflist)then
          tfnonlistq_dlist=.false.
          call loc_dlist(iand(ktamask,k),kl)
        else
          tfnonlistq_dlist=.true.
        endif          
        return
        end function tfnonlistq_dlist

        logical*4 function tfnonlistqd_rlist(k,kl)
        implicit none
        type (sad_rlist), pointer, optional, intent(out) :: kl
        type (sad_descriptor), intent(in):: k
        if(iand(ktfmask,k%k) .eq. ktflist .and.
     $       klist(ktfaddr(k%k)) .eq. ktfoper+mtflist)then
          tfnonlistqd_rlist=.false.
          if(present(kl))then
            call loc_sad(iand(ktamask,k%k),kl)
          endif
        else
          tfnonlistqd_rlist=.true.
        endif          
        return
        end function tfnonlistqd_rlist

        logical*4 function tfnonlistqd_dlist(k,kl)
        implicit none
        type (sad_dlist), pointer, intent(out) :: kl
        type (sad_descriptor), intent(in):: k
        if(iand(ktfmask,k%k) .eq. ktflist .and.
     $       klist(ktfaddr(k%k)) .eq. ktfoper+mtflist)then
          tfnonlistqd_dlist=.false.
          call loc_dlist(iand(ktamask,k%k),kl)
        else
          tfnonlistqd_dlist=.true.
        endif          
        return
        end function tfnonlistqd_dlist

        logical*4 function tfreallistqd(k,kl)
        implicit none
        type (sad_descriptor), intent(in):: k
        type (sad_rlist), pointer, optional, intent(out) :: kl
        type (sad_dlist), pointer :: kl1
        if(ktflistq(k,kl1))then
          tfreallistqd=kl1%head%k .eq. ktfoper+mtflist
     $         .and. ktfreallistq(kl1)
          if(tfreallistqd .and. present(kl))then
            call descr_sad(k,kl)
          endif
        else
          tfreallistqd=.false.
        endif
        return
        end function

        logical*4 function tfreallistqk(k,kl)
        use iso_c_binding
        implicit none
        type (sad_rlist), pointer, optional, intent(out) :: kl
        type (sad_rlist), pointer :: kl1
        integer*8 k
        if(ktflistq(k,kl1))then
          tfreallistqk=kl1%head%k .eq. ktfoper+mtflist
     $         .and. ktfreallistq(kl1)
          if(tfreallistqk .and. present(kl))then
            kl=>kl1
          endif
        else
          tfreallistqk=.false.
        endif
        return
        end function

        logical*4 function tfnonreallistqd(k,kl)
        implicit none
        type (sad_descriptor) , intent(in)::k
        type (sad_rlist), pointer, optional, intent(out) :: kl
        tfnonreallistqd=.not. tfreallistqd(k,kl)
        return
        end function

        logical*4 function tfnonreallistqk(k,kl)
        use iso_c_binding
        implicit none
        type (sad_rlist), pointer, optional, intent(out) :: kl
        integer*8 , intent(in)::k
        tfnonreallistqk=.not. tfreallistqk(k,kl)
        return
        end function

        logical*4 function ktflistqx(k,cx)
        implicit none
        type (sad_complex), pointer, optional, intent(out) :: cx
        integer*8 , intent(in)::k
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

        logical*4 function tfnonlistqk_rlist(k,kl)
        implicit none
        type (sad_rlist), pointer, intent(out) :: kl
        integer*8 , intent(in)::k
        integer*8 ka
        if(iand(ktfmask,k) .eq. ktflist)then
          ka=ktfaddr(k)
          if(klist(ka) .eq. ktfoper+mtflist)then
            tfnonlistqk_rlist=.false.
            call loc_rlist(ka,kl)
            return
          endif
        endif
        tfnonlistqk_rlist=.true.
        return
        end function tfnonlistqk_rlist

        logical*4 function tfnonlistqk_dlist(k,kl)
        implicit none
        type (sad_dlist), pointer, intent(out) :: kl
        integer*8 , intent(in)::k
        integer*8 ka
        if(iand(ktfmask,k) .eq. ktflist)then
          ka=ktfaddr(k)
          if(klist(ka) .eq. ktfoper+mtflist)then
            tfnonlistqk_dlist=.false.
            call loc_dlist(ka,kl)
            return
          endif
        endif
        tfnonlistqk_dlist=.true.
        return
        end function tfnonlistqk_dlist

        recursive logical*4 function tfruleqk_dlist(k,klx) result(lx)
        implicit none
        type (sad_dlist), pointer :: kl
        type (sad_dlist), pointer, optional, intent(out) :: klx
        integer*8 , intent(in)::k
        integer*4 i
        lx=.false.
        if(ktflistq(k,kl))then
          select case (kl%head%k)
          case (ktfoper+mtflist)
            if(.not. ktfnonreallistqo(kl))return
            do i=1,kl%nl
              if(.not. tfruleqk_dlist(kl%dbody(i)%k))then
                return
              endif
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
        end function tfruleqk_dlist

        logical*4 function tfruleqd_dlist(k,klx)
        implicit none
        type (sad_descriptor) , intent(in)::k
        type (sad_dlist), pointer, optional, intent(out) :: klx
        tfruleqd_dlist=tfruleqk_dlist(k%k,klx)
        return
        end function tfruleqd_dlist

        logical*4 function tfnumberqk(k,c)
        implicit none
        integer*8 , intent(in)::k
        type (sad_complex), pointer :: cx
        complex*16, optional, intent(out) :: c
        real*8 v
        if(ktfrealq(k,v))then
          tfnumberqk=.true.
          if(present(c))then
            c=v
          endif
        else
          tfnumberqk=tfcomplexqx(k,cx)
          if(tfnumberqk .and. present(c))then
            c=cx%cx(1)
          endif
        endif
        return
        end function

        logical*4 function tfnumberqd(k,c)
        implicit none
        type (sad_descriptor) , intent(in)::k
        type (sad_complex), pointer :: cx
        complex*16, optional, intent(out) :: c
        real*8 v
        if(ktfrealq(k,v))then
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

        logical*4 function tfnumlistqnd(k,n,kl1)
        implicit none
        type (sad_descriptor) , intent(in)::k
        type (sad_rlist), pointer, optional, intent(out) :: kl1
        type (sad_dlist), pointer :: kl
        integer*4 , intent(in)::n 
        if(ktflistq(k,kl))then
          tfnumlistqnd=kl%head%k .eq. ktfoper+mtflist
     $         .and. ktfreallistq(kl) .and. kl%nl .eq. n
          if(present(kl1))then
            call descr_rlist(k,kl1)
          endif
        else
          tfnumlistqnd=.false.
        endif
        return
        end function

        logical*4 function tfnumlistqnk(k,n,kl1)
        implicit none
        integer*8 , intent(in)::k
        type (sad_rlist), pointer, optional, intent(out) :: kl1
        type (sad_dlist), pointer :: kl
        integer*4 , intent(in)::n 
        if(ktflistq(k,kl))then
          tfnumlistqnk=kl%head%k .eq. ktfoper+mtflist
     $         .and. ktfreallistq(kl) .and. kl%nl .eq. n
          if(present(kl1))then
            call loc_rlist(ktfaddr(k),kl1)
          endif
        else
          tfnumlistqnk=.false.
        endif
        return
        end function
        
        logical*4 function tfcomplexnumlistqk_dlist(k,klx)
        implicit none
        type (sad_dlist), pointer, optional, intent(out) :: klx
        type (sad_dlist), pointer :: kl
        integer*8 , intent(in)::k
        integer*4 i
        tfcomplexnumlistqk_dlist=.false.
        if(tflistq(k,kl))then
          if(ktfnonreallistqo(kl))then
            do i=1,kl%nl
              if(.not. tfnumberqd(kl%dbody(i)))then
                return
              endif
            enddo
          endif
          tfcomplexnumlistqk_dlist=.true.
          if(present(klx))then
            klx=>kl
          endif
        endif
        return
        end function
        
        logical*4 function ktfsymbolqk(k,sym)
        implicit none
        type (sad_symbol), pointer, optional, intent(out) :: sym
        integer*8 , intent(in)::k
        ktfsymbolqk=iand(ktfmask,k) .eq. ktfsymbol
        if(present(sym) .and. ktfsymbolqk)then
          call loc_sym(ktfaddr(k),sym)
        endif
        return
        end function ktfsymbolqk

        logical*4 function ktfnonsymbolqk(k,sym)
        implicit none
        type (sad_symbol), pointer, optional, intent(out) :: sym
        integer*8 , intent(in)::k
        ktfnonsymbolqk=iand(ktfmask,k) .ne. ktfsymbol
        if(present(sym) .and. .not. ktfnonsymbolqk)then
          call loc_sym(ktfaddr(k),sym)
        endif
        return
        end function ktfnonsymbolqk

        logical*4 function ktfsymbolqd(k,sym)
        implicit none
        type (sad_descriptor) , intent(in)::k
        type (sad_symbol), pointer, optional, intent(out) :: sym
        ktfsymbolqd=iand(ktfmask,k%k) .eq. ktfsymbol
        if(present(sym) .and. ktfsymbolqd)then
          call loc_sym(ktfaddr(k%k),sym)
        endif
        return
        end function ktfsymbolqd

        logical*4 function ktfnonsymbolqd(k,sym)
        implicit none
        type (sad_descriptor) , intent(in)::k
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
        integer*8 , intent(in)::k
        ktfsymbolqdef=iand(ktfmask,k) .eq. ktfsymbol
        if(present(symd) .and. ktfsymbolqdef)then
          call loc_symdef(ktfaddr(k),symd)
        endif
        return
        end function ktfsymbolqdef

        logical*4 function ktfpatqk(k,pat)
        implicit none
        type (sad_pat), pointer, optional, intent(out) :: pat
        integer*8 , intent(in)::k
        if(iand(ktfmask,k) .eq. ktfpat)then
          ktfpatqk=.true.
          if(present(pat))then
            call loc_pat(iand(ktamask,k),pat)
          endif
        else
          ktfpatqk=.false.
        endif
        return
        end function ktfpatqk

        logical*4 function ktfpatqd(k,pat)
        implicit none
        type (sad_descriptor) , intent(in)::k
        type (sad_pat), pointer, optional, intent(out) :: pat
        ktfpatqd=ktfpatqk(k%k,pat)
        return
        end function ktfpatqd

        logical*4 function ktfnonpatq(k,pat)
        implicit none
        type (sad_pat), pointer, optional, intent(out) :: pat
        integer*8 , intent(in)::k
        ktfnonpatq=iand(ktfmask,k) .ne. ktfpat
        if(present(pat) .and. .not. ktfnonpatq)then
          call loc_pat(ktfaddr(k),pat)
        endif
        return
        end function ktfnonpatq

        logical*4 function ktfrefqk(k,ka)
        implicit none
        integer*8 , intent(in)::k
        integer*8 , optional, intent(out) :: ka
        ktfrefqk=iand(ktfmask,k) .eq. ktfref
        if(ktfrefqk .and. present(ka))then
          ka=ktfaddr(k)
        endif
        return
        end function ktfrefqk

        logical*4 function ktfrefqd(k,ka)
        implicit none
        type (sad_descriptor) , intent(in)::k
        integer*8 , optional, intent(out) :: ka
        ktfrefqd=iand(ktfmask,k%k) .eq. ktfref
        if(ktfrefqd .and. present(ka))then
          ka=ktfaddr(k%k)
        endif
        return
        end function ktfrefqd

        logical*4 function ktfnonrefq(k)
        implicit none
        integer*8 , intent(in)::k
        ktfnonrefq=iand(ktfmask,k) .ne. ktfref
        return
        end function ktfnonrefq

        logical*4 function ktfreallistqk(ka)
        implicit none
        integer*8 , intent(in)::ka
        ktfreallistqk=iand(ilist(2,ka-3),lnonreallist) .eq. 0
        return
        end function ktfreallistqk

        logical*4 function ktfreallistqk_rlist(ka,kl)
        implicit none
        type (sad_rlist), pointer, intent(out) :: kl
        integer*8 , intent(in)::ka
        ktfreallistqk_rlist=iand(ilist(2,ka-3),lnonreallist) .eq. 0
        if(ktfreallistqk_rlist)then
          call loc_sad(ka,kl)
        endif
        return
        end function ktfreallistqk_rlist

        logical*4 function ktfreallistqd(ka)
        implicit none
        type (sad_descriptor) , intent(in)::ka
        ktfreallistqd=iand(ilist(2,ka%k-3),lnonreallist) .eq. 0
        return
        end function ktfreallistqd

        logical*4 function ktfnonreallistq(ka)
        implicit none
        integer*8 , intent(in)::ka
        ktfnonreallistq=iand(ilist(2,ka-3),lnonreallist) .ne. 0
        return
        end function ktfnonreallistq

        logical*4 function ktfreallistqo_rlist(list)
        implicit none
        type (sad_rlist) , intent(in)::list
        ktfreallistqo_rlist=iand(list%attr,lnonreallist) .eq. 0
        return
        end function ktfreallistqo_rlist

        logical*4 function ktfreallistqo_dlist(list)
        implicit none
        type (sad_dlist) , intent(in)::list
        ktfreallistqo_dlist=iand(list%attr,lnonreallist) .eq. 0
        return
        end function ktfreallistqo_dlist

        logical*4 function ktfnonreallistqo_rlist(list)
        implicit none
        type (sad_rlist) , intent(in)::list
        ktfnonreallistqo_rlist=iand(list%attr,lnonreallist) .ne. 0
        return
        end function ktfnonreallistqo_rlist

        logical*4 function ktfnonreallistqo_dlist(list)
        implicit none
        type (sad_dlist) , intent(in)::list
        ktfnonreallistqo_dlist=iand(list%attr,lnonreallist) .ne. 0
        return
        end function ktfnonreallistqo_dlist

        logical*4 function ktftrueq(ka)
        implicit none
        integer*8 , intent(in)::ka
        ktftrueq=ka .ne. 0 .and. iand(ktrmask,ka) .ne. ktfnr
        return
        end function ktftrueq

        logical*4 function ktfsequenceq_rlist(k,kl) result(v)
        implicit none
        type (sad_rlist), pointer, intent(out) :: kl
        integer*8 , intent(in)::k
        v=iand(ktfmask,k) .eq. ktflist .and.
     $       klist(iand(ktamask,k)) .eq. ktfoper+mtfnull
        if(v)then
          call loc_sad(iand(ktamask,k),kl)
        endif
        return
        end function ktfsequenceq_rlist

        logical*4 function ktfsequenceq_dlist(k,kl) result(v)
        implicit none
        type (sad_dlist), pointer, optional, intent(out) :: kl
        integer*8 , intent(in)::k
        v=iand(ktfmask,k) .eq. ktflist .and.
     $       klist(iand(ktamask,k)) .eq. ktfoper+mtfnull
        if(v .and. present(kl))then
          call loc_dlist(iand(ktamask,k),kl)
        endif
        return
        end function ktfsequenceq_dlist

        logical*4 function ktfsequenceqd_rlist(k,kl) result(v)
        implicit none
        type (sad_rlist), pointer, intent(out) :: kl
        type (sad_descriptor) , intent(in)::k
        v=iand(ktfmask,k%k) .eq. ktflist .and.
     $       klist(iand(ktamask,k%k)) .eq. ktfoper+mtfnull
        if(v)then
          call loc_sad(iand(ktamask,k%k),kl)
        endif
        return
        end function ktfsequenceqd_rlist

        logical*4 function ktfsequenceqd_dlist(k,kl) result(v)
        implicit none
        type (sad_dlist), pointer, optional, intent(out) :: kl
        type (sad_descriptor) , intent(in)::k
        v=iand(ktfmask,k%k) .eq. ktflist .and.
     $       klist(iand(ktamask,k%k)) .eq. ktfoper+mtfnull
        if(v .and. present(kl))then
          call loc_dlist(iand(ktamask,k%k),kl)
        endif
        return
        end function ktfsequenceqd_dlist

        logical*4 function ktfprotectedq(ka)
        implicit none
        integer*8 , intent(in)::ka
        ktfprotectedq=iand(ilist(1,ka-3),iattrprotected) .ne. 0
        return
        end function ktfprotectedq

        logical*4 function ktfprotectedqo(sym)
        implicit none
        type (sad_symbol) , intent(in)::sym
        ktfprotectedqo=iand(sym%attr,iattrprotected) .ne. 0
        return
        end function ktfprotectedqo

        logical*4 function ktfconstantq(ka)
        implicit none
        integer*8 , intent(in)::ka
        ktfconstantq=iand(ilist(1,ka-3),iattrconstant) .ne. 0
        return
        end function ktfconstantq

        logical*4 function ktfconstantsymq(sym)
        implicit none
        type (sad_symbol) , intent(in)::sym
        ktfconstantsymq=iand(sym%attr,iattrconstant) .ne. 0
        return
        end function ktfconstantsymq

        logical*4 function ktfimmediateq(ka)
        implicit none
        integer*8 , intent(in)::ka
        ktfimmediateq=iand(ilist(1,ka-3),iattrimmediate) .ne. 0
        return
        end function ktfimmediateq

        logical*4 function ktfnumericq(ka)
        implicit none
        integer*8 , intent(in)::ka
        ktfnumericq=iand(ilist(1,ka-3),iattrnumeric) .ne. 0
        return
        end function ktfnumericq

        logical*4 function ktfovrwrtq(kl)
        implicit none
        type (sad_dlist) , intent(in)::kl
        ktfovrwrtq=kl%ref .le. 0 .or.
     $       kl%ref .eq. 1 .and. ktfaddr(kl%alloc) .eq. 0
        return
        end function ktfovrwrtq

        logical*4 function tfcomplexqdx(k,cx)
        implicit none
        type (sad_complex), pointer :: c
        type (sad_complex), pointer, optional, intent(out) :: cx
        type (sad_descriptor) , intent(in)::k
        tfcomplexqdx=ktflistqx(k%k,c) .and.
     $       c%head%k .eq. ktfoper+mtfcomplex
     $       .and. c%nl .eq. 2 .and.
     $       iand(lnonreallist,c%attr) .eq. 0
        if(tfcomplexqdx .and. present(cx))then
          cx=>c
        endif
        return
        end function

        logical*4 function tfcomplexqdc(k,cv)
        implicit none
        type (sad_complex), pointer :: c
        complex*16  cv
        type (sad_descriptor) , intent(in)::k
        tfcomplexqdc=ktflistqx(k%k,c) .and.
     $       c%head%k .eq. ktfoper+mtfcomplex
     $       .and. c%nl .eq. 2 .and.
     $       iand(lnonreallist,c%attr) .eq. 0
        if(tfcomplexqdc)then
          cv=c%cx(1)
        endif
        return
        end function

        logical*4 function tfcomplexqx(k,cx)
        implicit none
        type (sad_complex), pointer :: c
        type (sad_complex), pointer, optional, intent(out) :: cx
        integer*8 , intent(in)::k
        tfcomplexqx=ktflistqx(k,c) .and.
     $       c%head%k .eq. ktfoper+mtfcomplex
     $       .and. c%nl .eq. 2 .and.
     $       iand(lnonreallist,c%attr) .eq. 0
        if(tfcomplexqx .and. present(cx))then
          cx=>c
        endif
        return
        end function

        logical*4 function tfcomplexqc(k,cv)
        implicit none
        type (sad_complex), pointer :: c
        complex*16 cv
        integer*8 , intent(in)::k
        tfcomplexqc=ktflistqx(k,c) .and.
     $       c%head%k .eq. ktfoper+mtfcomplex
     $       .and. c%nl .eq. 2 .and.
     $       iand(lnonreallist,c%attr) .eq. 0
        if(tfcomplexqc)then
          cv=c%cx(1)
        endif
        return
        end function

        recursive logical*4 function tfsameqk(ka,kp) result(lx)
        use tfcode
        use iso_c_binding
        implicit none
        type (sad_dlist), pointer :: kla,klp
        type (sad_symbol), pointer :: syma,symp
        type (sad_string), pointer :: stra,strp
        type (sad_pat), pointer :: pata,patp
        integer*8 , intent(in)::ka,kp
        integer*8 nc
        if(ka .eq. kp)then
          lx=.true.
          return
        endif
        lx=.false.
        if(ktfrealq(ka) .or. ktfoperq(ka) .or.
     $       ktftype(ka) .ne. ktftype(kp))then
          return
        endif
        if(ktfsymbolq(ka,syma))then
          call loc_sad(ktfaddr(kp),symp)
          lx=tfsamesymbolqo(syma,symp)
        elseif(ktfstringq(ka,stra))then
          call loc_sad(ktfaddr(kp),strp)
          nc=stra%nch
          if(nc .eq. strp%nch)then
            lx=stra%str(1:nc) .eq. strp%str(1:nc)
          endif
        elseif(ktflistq(ka,kla))then
          call loc_sad(ktfaddr(kp),klp)
          if(kla%nl .eq. klp%nl)then
            lx=tfsamelistqo(kla,klp)
          endif
        elseif(ktfpatq(ka,pata))then
          call loc_sad(ktfaddr(kp),patp)
          if(.not. tfsameqk(pata%expr%k,patp%expr%k))then
            return
          endif
          if(.not. tfsameqk(pata%head%k,patp%head%k))then
            return
          endif
          if(.not. tfsameqk(pata%default%k,patp%default%k))then
            return
          endif
          lx=tfsamesymbolqk(ktfaddr(pata%sym%alloc),
     $         ktfaddr(patp%sym%alloc))
        endif
        return
        end function

        logical*4 function tfsameqd(ka,kp)
        implicit none
        type (sad_descriptor) , intent(in)::ka,kp
        tfsameqd=tfsameqk(ka%k,kp%k)
        return
        end function

        logical*4 function tfsamelistqo(lista,listp)
        implicit none
        type (sad_dlist) , intent(inout)::lista,listp
        type (sad_descriptor) kai,kpi
        integer*8 kaai,kapi
        integer*4 i,m
        tfsamelistqo=.false.
        m=lista%nl
        if(m .ne. listp%nl)then
          return
        endif
        do i=0,m
          kai=lista%dbody(i)
          kpi=listp%dbody(i)
          if(kai%k .ne. kpi%k)then
            if(ktfobjq(kai))then
              if(tfsameq(kai,kpi))then
                kaai=ktfaddr(kai%k)
                kapi=ktfaddr(kpi%k)
                if(ilist(1,kapi-1) .ge. ilist(1,kaai-1))then
                  call tflocal1(kai)
                  lista%dbody(i)=dtfcopy1(kpi)
                else
                  call tflocal1(kpi)
                  listp%dbody(i)=dtfcopy1(kai)
                endif
                cycle
              endif
            endif
            return
          endif
        enddo
        tfsamelistqo=.true.
        return
        end

        logical*4 function tfsamesymbolqo(sa,sp)
        use tfcode
        implicit none
        type(sad_symbol) , intent(in)::sa,sp
        tfsamesymbolqo=sa%loc .eq. sp%loc .and.
     $       max(0,sa%gen) .eq. max(0,sp%gen)
        return
        end function

        logical*4 function tfsamesymbolqk(ka1,kp1)
        implicit none
        type (sad_symbol) ,pointer :: syma,symp
        integer*8 , intent(in)::ka1,kp1
        integer*8 ka,kp
        ka=ktfaddr(ka1)
        kp=ktfaddr(kp1)
        if(ka .eq. kp)then
          tfsamesymbolqk=.true.
        elseif(ka .eq. 0 .or. kp .eq. 0)then
          tfsamesymbolqk=.false.
        else
          call loc_sad(ka,syma)
          call loc_sad(kp,symp)
          tfsamesymbolqk=tfsamesymbolqo(syma,symp)
        endif
        return
        end function

        logical*4 function tfsamesymbolqd(k1,k2)
        implicit none
        type (sad_descriptor) k1,k2
        tfsamesymbolqd=tfsamesymbolqk(k1%k,k2%k)
        return
        end function

        logical*4 function tfsamestringqk(ka1,kp1)
        implicit none
        type (sad_string), pointer :: stra,strp
        integer*8 ka1,kp1
        if(ktfaddr(ka1) .eq. ktfaddr(kp1))then
          tfsamestringqk=.true.
        else
          call loc_sad(ktfaddr(ka1),stra)
          call loc_sad(ktfaddr(kp1),strp)
          tfsamestringqk=tfsamestringqo(stra,strp)
        endif
        return
        end function

        logical*4 function tfsamestringqd(k1,k2)
        implicit none
        type (sad_descriptor) k1,k2
        tfsamestringqd=tfsamestringqk(k1%k,k2%k)
        return
        end function

        logical*4 function tfsamestringqo(sa,sp)
        use tfcode
        implicit none
        type(sad_string) , intent(in)::sa,sp
        tfsamestringqo=sa%nch .eq. sp%nch .and.
     $       sa%str(:sa%nch) .eq. sp%str(:sp%nch)
        return
        end function

        logical*4 function tfexprqd(k)
        implicit none
        type (sad_descriptor) , intent(in)::k
        type (sad_dlist), pointer :: kl
        tfexprqd=ktflistq(k,kl) .and. kl%head%k .ne. ktfoper+mtflist
        return
        end function

        logical*4 function tfexprqk(k)
        implicit none
        integer*8 , intent(in)::k
        type (sad_dlist), pointer :: kl
        tfexprqk=ktflistq(k,kl) .and. kl%head%k .ne. ktfoper+mtflist
        return
        end function

        logical*4 function tfsameheadqk(k1,k2)
        implicit none
        integer*8 , intent(in)::k1,k2
        tfsameheadqk=tfsameq(klist(ktfaddr(k1)),klist(ktfaddr(k2)))
        return
        end function

        logical*4 function tfsameheadqd(k1,k2)
        implicit none
        type (sad_descriptor) , intent(in)::k1,k2
        tfsameheadqd=tfsameq(klist(ktfaddr(k1)),klist(ktfaddr(k2)))
        return
        end function

        logical function tfinequalityqk(k)
        implicit none
        integer*8 , intent(in)::k
        type (sad_dlist), pointer :: kl
        tfinequalityqk=ktflistq(k,kl) .and.
     $       kl%head%k .eq. ktfoper+mtfinequality
        return
        end function

        logical function tfinequalityqd(k)
        implicit none
        type (sad_descriptor) , intent(in)::k
        type (sad_dlist), pointer :: kl
        tfinequalityqd=ktflistq(k,kl) .and.
     $       kl%head%k .eq. ktfoper+mtfinequality
        return
        end function

        recursive logical*4 function tfconstqk(k) result(lx)
        implicit none
        integer*8 , intent(in)::k
        type (sad_dlist), pointer :: kl
        type (sad_symdef), pointer ::symd
        type (sad_pat), pointer :: pat
        logical*4 tfconstlistqo
        lx=.true.
        if(ktfsymbolqdef(k,symd))then
          lx=ktfconstantsymq(symd%sym) .and.
     $         symd%value%k .eq. ktfsymbol+ktfaddr(k)
     $         .and. symd%upval .eq. 0
        elseif(ktflistq(k,kl))then
          lx=tfconstlistqo(kl)
        elseif(ktfpatq(k,pat))then
          lx=tfconstqk(pat%expr%k)
        endif
        return
        end function

        logical*4 function tfconstqd(k)
        implicit none
        type (sad_descriptor) , intent(in)::k
        tfconstqd=tfconstqk(k%k)
        return
        end function

        logical*4 function tfconstpatternqk(k)
        implicit none
        integer*8 , intent(in)::k
        type (sad_dlist), pointer :: kl
        logical*4 tfconstpatternheadqk,tfconstpatternlistbodyqo
        tfconstpatternqk=.true.
        if(ktfpatq(k))then
          tfconstpatternqk=.false.
        elseif(ktflistq(k,kl))then
          tfconstpatternqk=tfconstpatternheadqk(kl%head) .and.
     $         tfconstpatternlistbodyqo(kl)
        endif
        return
        end function

        logical*4 function tfconstpatternqd(k)
        implicit none
        type (sad_descriptor) , intent(in)::k
        type (sad_dlist), pointer :: kl
        logical*4 tfconstpatternheadqk,tfconstpatternlistbodyqo
        tfconstpatternqd=.true.
        if(ktfpatq(k))then
          tfconstpatternqd=.false.
        elseif(ktflistq(k,kl))then
          tfconstpatternqd=tfconstpatternheadqk(kl%head) .and.
     $         tfconstpatternlistbodyqo(kl)
        endif
        return
        end function

        logical*4 function tfheldqd(k)
        implicit none
        type (sad_descriptor) , intent(in)::k
        type (sad_dlist), pointer :: kl
        tfheldqd=ktflistq(k,kl) .and.
     $       kl%head%k .eq. ktfoper+mtfhold
        return
        end function

        logical*4 function tfmatrixqd(k,kl,klind,klbody)
        implicit none
        type (sad_descriptor) , intent(in)::k
        type (sad_dlist), pointer, optional,
     $       intent(out) :: kl,klbody
        type (sad_rlist), pointer, optional,
     $       intent(out) :: klind
        type (sad_rlist), pointer :: kli1
        type (sad_dlist), pointer :: kl1,klb1
        if(ktflistq(k,kl1))then
          if(tfsameq(kl1%head,kxmatrix) .and. kl1%nl .eq. 2)then
            if(tfreallistq(kl1%dbody(1),kli1) .and.
     $           ktflistq(kl1%dbody(2),klb1))then
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
        real*8 , intent(in)::r(2)
        cfromr=dcmplx(r(1),r(2))
        return
        end function cfromr

        integer*8 function ktfcopy1(k)
        implicit none
        integer*8 , intent(in)::k
        integer*8 ka
        ka=iand(ktamask,k)
        ilist(1,ka-1)=ilist(1,ka-1)+1
        ktfcopy1=k
        return
        end function ktfcopy1

        integer*8 function ktfcopyd(k,d)
        implicit none
        integer*8 , intent(in)::k
        integer*8 ka
        logical*4 , intent(out)::d
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
        integer*8 , intent(in)::k
        integer*8 ka
        if(ktfobjq(k))then
          ka=ktfaddr(k)
          ilist(1,ka-1)=ilist(1,ka-1)+1
        endif
        ktfcopy=k
        return
        end function

        type (sad_descriptor) function dtfcopy(d)
        implicit none
        type (sad_descriptor) , intent(in)::d
        dtfcopy%k=ktfcopy(d%k)
        return
        end function

        type (sad_descriptor) function dtfcopy1(d)
        implicit none
        type (sad_descriptor) , intent(in)::d
        dtfcopy1%k=ktfcopy1(d%k)
        return
        end function

        subroutine tflocald(k)
        implicit none
        type (sad_descriptor) , intent(in)::k
        type (sad_object), pointer :: obj
        integer*8 itfroot
        if(ktfobjqd(k,obj))then
          obj%ref=obj%ref-1
          if(obj%ref .le. 0)then
            obj%ref=0
            if(ktfaddr(obj%alloc) .eq. 0)then
              itfroot=itflocal+levele
              obj%alloc%k=obj%alloc%k+ktfaddr(klist(itfroot))
              klist(itfroot)=ktfaddr(k%k)-2
            endif
          endif
        endif
        return
        end subroutine 

        subroutine tflocal1d(k)
        implicit none
        type (sad_descriptor) , intent(in)::k
        type (sad_object), pointer :: obj
        integer*8 ka,itfroot
        ka=ktfaddr(k)
        call loc_obj(ka,obj)
        obj%ref=obj%ref-1
        if(obj%ref .le. 0)then
          obj%ref=0
          if(ktfaddr(obj%alloc) .eq. 0)then
            itfroot=itflocal+levele
            obj%alloc%k=obj%alloc%k+ktfaddr(klist(itfroot))
            klist(itfroot)=ka-2
          endif
        endif
        return
        end subroutine

        subroutine tfconnect(k,irtc)
        implicit none
        type (sad_descriptor) , intent(in)::k
        integer*4 irtc
        call tfconnectk(k%k,irtc)
        return
        end subroutine

        subroutine tfconnectk(k,irtc)
        implicit none
        type (sad_object), pointer :: obj
        integer*4 l,itfdownlevel
        integer*4 , intent(out)::irtc
        integer*8 , intent(in)::k
        integer*8 ka,j
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
              obj%alloc%k=ktftype(obj%alloc%k)+klist(j)
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
        integer*4 l,itfdownlevel
        integer*4 , intent(out)::irtc
        integer*8 , intent(in)::k
        integer*8 ka
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

        type (sad_descriptor) function kxaaloc_dlist(mode,nd,kl)
        implicit none
        type (sad_dlist), pointer, optional, intent(out) :: kl
        integer*4 , intent(in)::mode,nd
        kxaaloc_dlist%k=ktflist+ktaaloc_dlist(mode,nd,kl)
        return
        end function

        type (sad_descriptor) function kxaaloc_rlist(mode,nd,kl)
        implicit none
        type (sad_rlist), pointer, intent(out) :: kl
        integer*4 , intent(in)::mode,nd
        kxaaloc_rlist%k=ktflist+ktaaloc_rlist(mode,nd,kl)
        return
        end function

        type (sad_descriptor) function kxadaloc_dlist(mode,nd,kl)
        implicit none
        type (sad_dlist), pointer, optional, intent(out) :: kl
        integer*4 , intent(in)::mode,nd
        kxadaloc_dlist%k=ktflist+ktadaloc(mode,nd,kl)
        return
        end function

        type (sad_descriptor) function kxadalocnull(mode,nd,kl)
        implicit none
        type (sad_dlist), pointer, optional, intent(out) :: kl
        integer*4 , intent(in)::mode,nd
        kxadalocnull%k=ktflist+ktadalocnull(mode,nd,kl)
        return
        end function

        type (sad_descriptor) function kxavaloc(mode,nd,kl)
        implicit none
        type (sad_rlist), pointer, optional, intent(out) :: kl
        integer*4 , intent(in)::mode,nd
        kxavaloc%k=ktflist+ktavaloc(mode,nd,kl)
        return
        end function

        type (sad_descriptor) function kxmakelist_dlist(isp1,kl)
        implicit none
        type (sad_dlist), pointer, optional, intent(out) :: kl
        integer*4 , intent(in)::isp1
        kxmakelist_dlist%k=ktflist+ktfmakelist(isp1,kl)
        return
        end function

        type (sad_descriptor) function kxmakelist_rlist(isp1,kl)
        implicit none
        type (sad_rlist), pointer, intent(out) :: kl
        integer*4 , intent(in)::isp1
        kxmakelist_rlist%k=ktflist+ktfmakelist(isp1,kl)
        return
        end function

        type (sad_descriptor) function kxmakelist0(isp1,kl)
        implicit none
        type (sad_dlist), pointer, optional, intent(out) :: kl
        integer*4 , intent(in)::isp1
        kxmakelist0%k=ktflist+ktfmakelist0(isp1,kl)
        return
        end function

        type (sad_descriptor) function kxcompose(isp1,kl)
        implicit none
        type (sad_dlist), pointer, optional, intent(out) :: kl
        integer*4 , intent(in)::isp1
        kxcompose%k=ktflist+ktfcompose(isp1,kl)
        return
        end function

        type (sad_descriptor) function kxcomposer(isp1)
        implicit none
        integer*4 , intent(in)::isp1
        integer*8 ktfcrelistr
        kxcomposer%k=ktflist+
     $       ktfcrelistr(isp-isp1,ktastk(isp1+1),ktastk(isp1))
        return
        end function

        type (sad_descriptor) function kxcomposev(isp0,kl)
        implicit none
        type (sad_dlist), pointer, optional, intent(out) :: kl
        integer*4 , intent(in)::isp0
        kxcomposev%k=ktflist+ktfcomposev(isp0,kl)
        return
        end function

        type (sad_descriptor) function kxcalocv(mode,x,y)
        implicit none
        integer*4 , intent(in)::mode
        real*8 , intent(in)::x,y
        kxcalocv%k=ktflist+ktcalocv(mode,x,y)
        return
        end function

        type (sad_descriptor) function kxsalocb(mode,string,leng,str)
        implicit none
        type (sad_string), pointer, optional, intent(out) :: str
        integer*4 , intent(in)::mode,leng
        character , intent(in)::string(leng)
        kxsalocb%k=ktfstring+ktsalocb(mode,string,leng)
        if(present(str))then
          call descr_sad(kxsalocb,str)
        endif
        return
        end function

        type (sad_descriptor) function kxsalocbb(mode,leng,str)
        implicit none
        type (sad_string), pointer, optional, intent(out) :: str
        integer*4 , intent(in)::mode,leng
        kxsalocbb%k=ktfstring+ktsalocbb(mode,leng)
        if(present(str))then
          call descr_sad(kxsalocbb,str)
        endif
        return
        end function

        type (sad_descriptor) function kxsymbolz(name,l,symd)
        implicit none
        type (sad_symdef), pointer, optional, intent(out) :: symd
        integer*4 , intent(in)::l
        character , intent(in)::name(l)
        kxsymbolz%k=ktfsymbol+ktfsymbolz(name,l,symd)
        return
        end function

        type (sad_descriptor) function kxsymbolv(name,l,symd0)
        implicit none
        type (sad_symdef), pointer, optional, intent(out) :: symd0
        type (sad_symdef), pointer :: symd
        integer*8 k
        integer*4 , intent(in)::l
        character , intent(in)::name(l)
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
        integer*4 , intent(in)::l
        character , intent(in)::name(l)
        logical*4 , intent(in)::const
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
        integer*4 , intent(in)::l
        character , intent(in)::name(l)
        real*8 , intent(in)::val
        type (sad_symdef), pointer :: symd
        type (sad_descriptor) kx
        kx=kxsymbolz(name,l,symd)
        symd%value=dfromr(val)
        return
        end subroutine

        type (sad_descriptor) function k_descr(k)
        implicit none
        integer*8 , intent(in)::k
        k_descr%k=k
        return
        end function

        type (sad_descriptor) function kxscopy(ka,leng,str)
        implicit none
        type (sad_string), pointer, optional, intent(out) :: str
        integer*8 , intent(in)::ka
        integer*8 ktfaloc,l
        integer*4 , intent(in)::leng
        integer*4 nw,n
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
        integer*8 , intent(in)::locp
        integer*4 , intent(in)::lg
        kxnaloc1=kxnaloc(lg,locp,0)
        return
        end function

        type (sad_descriptor) function kxnaloc(lg,locp,n)
        implicit none
        integer*8 kp,kp1, kp0,ktalocr
        integer*8 , intent(in)::locp
        integer*4 , intent(in)::lg
        integer*4 ipg,n
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
        kxnaloc=sad_descr(def%sym)
        return
        end function

        real*8 function rfromd(d)
        implicit none
        type (sad_descriptor) , intent(in)::d
        real*8 rfromk
        rfromd=rfromk(d%k)
        return
        end function

        integer*4 function ifromd(d)
        implicit none
        type (sad_descriptor) , intent(in)::d
        real*8 rfromk
        ifromd=int(rfromk(d%k))
        return
        end function

        type (sad_descriptor) function dfromr(x)
        implicit none
        integer*8 kfromr
        real*8 , intent(in)::x
        dfromr%k=kfromr(x)
        return
        end function

        type (sad_descriptor) function dfromk(k)
        implicit none
        integer*8 , intent(in)::k
        dfromk%k=k
        return
        end function

        integer*8 function ktfmakelist0(isp1,kl)
        implicit none
        type (sad_dlist), pointer, optional, intent(out) :: kl
        integer*4 , intent(in)::isp1
        integer*4 narg
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
        type (sad_descriptor) , intent(in)::kh
        type (sad_dlist), pointer ::kl
        integer*4 , intent(in)::m
        integer*8 , intent(in)::ks(m)
        kxcrelistm=kxaaloc(-1,m,kl)
        call tfcrelista(m,ks,kh,kl)
        return
        end function

        integer*8 function ktsalocb(mode,string,leng)
        integer*8 ktfaloc,l
        integer*4 , intent(in)::leng,mode
        integer*4 nw,ik
        character , intent(in)::string(leng)
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
        type (sad_rlist), pointer :: kl
        integer*4 mode
        real*8 x,y
        ktcalocv=ktavaloc(mode,2,kl)
        kl%attr=lconstlist
        kl%head%k=ktfoper+mtfcomplex
        kl%rbody(1)=x
        kl%rbody(2)=y
        return
        end function

        integer*8 function ktadalocnull(mode,nd,kl)
        implicit none
        integer*4 , intent(in)::mode,nd
        type (sad_dlist), pointer :: kl1
        type (sad_dlist), pointer, optional, intent(out) :: kl
        ktadalocnull=ktadaloc(mode,nd,kl1)
        kl1%body(1:nd)=ktfoper+mtfnull
        if(present(kl))then
          kl=>kl1
        endif
        return
        end function

        integer*8 function ktaaloc_dlist(mode,nd,kl)
        implicit none
        type (sad_dlist), pointer, optional, intent(out) :: kl
        integer*8 ka,ktfaloc
        integer*4 , intent(in)::nd,mode
        ka=ktfaloc(mode,ktflist,nd+1)
        ilist(1,ka-3)=0
        ilist(2,ka-1)=nd
        klist(ka)=ktfoper+mtflist
        ktaaloc_dlist=ka
        if(present(kl))then
          call loc_dlist(ka,kl)
        endif
        return
        end function

        integer*8 function ktaaloc_rlist(mode,nd,kl)
        implicit none
        type (sad_rlist), pointer, intent(out) :: kl
        integer*8 ka,ktfaloc
        integer*4 , intent(in)::nd,mode
        ka=ktfaloc(mode,ktflist,nd+1)
        ilist(1,ka-3)=0
        ilist(2,ka-1)=nd
        klist(ka)=ktfoper+mtflist
        ktaaloc_rlist=ka
        call loc_rlist(ka,kl)
        return
        end function

        integer*8 function ktavaloc(mode,nd,kl)
        implicit none
        type (sad_rlist), pointer, optional, intent(out) :: kl
        integer*4 , intent(in)::mode,nd
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
          call loc_rlist(k1+2,kl)
        endif
        return
        end function

        integer*8 function ktadaloc_dlist(mode,nd,kl)
        implicit none
        type (sad_dlist), pointer, optional, intent(out) :: kl
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
        ktadaloc_dlist=k1+2
        if(present(kl))then
          call loc_dlist(k1+2,kl)
        endif
        return
        end function

        integer*8 function ktaalocr(mode,nd,kl)
        implicit none
        type (sad_dlist), pointer, optional, intent(out) :: kl
        integer*8 ka,ktfalocr
        integer*4 , intent(in)::nd,mode
        ka=ktfalocr(mode,ktflist,nd+1)
        ilist(1,ka-3)=0
        ilist(2,ka-1)=nd
        klist(ka)=ktfoper+mtflist
        ktaalocr=ka
        if(present(kl))then
          call loc_sad(ka,kl)
        endif
        return
        end function

        integer*8 function ktraaloc(mode,nd,kl)
        implicit none
        type (sad_rlist), pointer, optional, intent(out) :: kl
        type (sad_rlist), pointer :: kl1
        integer*4 , intent(in)::mode,nd
        integer*8 ka
        ka=ktavaloc(mode,nd,kl1)
        kl1%rbody(1:nd)=0.d0
        if(present(kl))then
          kl=>kl1
        endif
        ktraaloc=ka
        return
        end function

        type (sad_descriptor) function kxraaloc(mode,nd,kl)
        implicit none
        type (sad_rlist), pointer, optional, intent(out) :: kl
        type (sad_rlist), pointer :: kl1
        integer*4 , intent(in)::mode,nd
        integer*8 ka
        ka=ktavaloc(mode,nd,kl1)
        kl1%rbody(1:nd)=0.d0
        if(present(kl))then
          kl=>kl1
        endif
        kxraaloc%k=ktflist+ka
        return
        end function

        integer*8 function ktaalocsp(nd,lp,la,kl1)
        implicit none
        type (sad_dlist), pointer, optional, intent(out) :: kl1
        type (sad_dlist), pointer :: kl
        integer*8 ka
        integer*4 , intent(in)::nd
        integer*2 , intent(in)::lp,la
        ka=ktaloc(nd+lp+la+3)+lp+2
        call loc_sad(ka,kl)
        kl%lenp=lp
        kl%lena=la
        kl%attr=0
        kl%alloc%k=ktflist
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
        type (sad_dlist), pointer, optional, intent(out) :: kl
        integer*8 kax,kh
        integer*4 , intent(in)::isp0
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
          call loc_sad(ktfcomposev,kl)
        endif
        return
        end function

        integer*8 function ktfcompose(isp1,kl)
        implicit none
        type (sad_dlist), pointer, optional, intent(out) :: kl
        integer*4 , intent(in)::isp1
        ktfcompose=ktfaddr(
     $       kxcrelistm(isp-isp1,ktastk(isp1+1:isp),dtastk(isp1)))
        if(present(kl))then
          call loc_sad(ktfcompose,kl)
        endif
        return
        end function

        integer*8 function ktfmakelist_dlist(isp1,kl)
        implicit none
        type (sad_dlist), pointer, optional, intent(out) :: kl
        integer*4 , intent(in)::isp1
        integer*4 narg
        narg=isp-isp1
        if(narg .eq. 1)then
          if(ktastk(isp) .eq. ktfoper+mtfnull)then
            ktfmakelist_dlist=ktaaloc(-1,0,kl)
            return
          endif
        endif
        ktfmakelist_dlist=
     $       ktfaddr(kxcrelistm(narg,ktastk(isp1+1:isp1+narg),
     $       k_descr(ktfoper+mtfleftbrace)))
        if(present(kl))then
          call loc_dlist(ktfmakelist_dlist,kl)
        endif
        return
        end function

        integer*8 function ktfmakelist_rlist(isp1,kl)
        implicit none
        type (sad_rlist), pointer, intent(out) :: kl
        integer*4 , intent(in)::isp1
        integer*4 narg
        narg=isp-isp1
        if(narg .eq. 1)then
          if(ktastk(isp) .eq. ktfoper+mtfnull)then
            ktfmakelist_rlist=ktavaloc(-1,0,kl)
            return
          endif
        endif
        ktfmakelist_rlist=
     $       ktfaddr(kxcrelistm(narg,ktastk(isp1+1:isp1+narg),
     $       k_descr(ktfoper+mtfleftbrace)))
        call loc_rlist(ktfmakelist_rlist,kl)
        return
        end function

        integer*8 function ktfsymbolz(name,l,symd)
        implicit none
        type (sad_symdef), pointer, optional, intent(out) :: symd
        integer*4 , intent(in)::l
        integer*8 ktfsymbolc
        character , intent(in)::name(l)
        ktfsymbolz=ktfsymbolc(name,l,int8(0))
        if(present(symd))then
          call loc_symdef(ktfsymbolz,symd)
        endif
        return
        end function

        subroutine tfsydef(sym,symx)
        implicit none
        type (sad_descriptor) kx
        type (sad_symbol) , intent(in)::sym
        type (sad_symbol), pointer, intent(out) :: symx
        call tfsydefg(sym%loc,kx,sym%gen)
        call loc_sad(ktfaddrd(kx),symx)
        return
        end subroutine

        type (sad_descriptor) function dxsycopy(sym)
        implicit none
        type (sad_symbol) , intent(in)::sym
        integer*8 kax,ktfaloc
        kax=ktfaloc(-1,ktfsymbol,1)
        ilist(2,kax-1)=maxgeneration
        klist(kax)=sym%loc
        dxsycopy%k=ktfsymbol+kax
        return
        end function

        type (sad_descriptor) function kxpfaloc(k)
        implicit none
        type (sad_descriptor) , intent(in)::k
        type (sad_dlist), pointer :: klx
        type (sad_rlist), pointer :: klr
        real*8 x
        if(ktfrealq(k,x))then
          kxpfaloc=kxavaloc(-1,1,klr)
          klr%rbody(1)=x
          klr%head%k=ktfoper+mtffun
        else
          kxpfaloc=kxadaloc(-1,1,klx)
          klx%dbody(1)=dtfcopy(k)
          klx%head%k=ktfoper+mtffun
        endif
        return
        end function

        type (sad_descriptor) function kxpaloc(string)
        implicit none
        character*(*) , intent(in)::string
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
        type (sad_descriptor) , intent(in)::kp,kh
        type (sad_descriptor) ks
        integer*4 , intent(in)::ls
        character , intent(in)::symb(ls)
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
        type (sad_descriptor) , intent(in)::kp,kh,kd,ks
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
        if(ktfsymbolq(ks,sym))then
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
        type (sad_descriptor) , intent(in)::kh
        type (sad_descriptor) kx
        type (sad_string), pointer :: str
        integer*4 n,ic1,ic2
        integer*4 , intent(in)::isp1,isp2
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
        type (sad_dlist), target, intent(in) :: list
        type (sad_dlist), pointer, intent(out) :: listc
        integer*4 i
        if(ktfovrwrtq(list))then
          listc=>list
        else
          call loc_dlist(ktaaloc(-1,list%nl),listc)
          listc%attr=list%attr
          if(ktfreallistq(list))then
            listc%head=dtfcopy(list%head)
            listc%dbody(1:list%nl)=list%dbody(1:list%nl)
          else
            do i=0,list%nl
              listc%dbody(i)=dtfcopy(list%dbody(i))
            enddo
          endif
        endif
        listc%attr=ior(listc%attr,ktoberebuilt)
        return
        end subroutine

        subroutine tfduplist(list,listc)
        implicit none
        type (sad_dlist), target, intent(in) :: list
        type (sad_dlist), pointer, intent(out) :: listc
        integer*4 i
        call loc_dlist(ktaaloc(-1,list%nl),listc)
        listc%attr=list%attr
        if(ktfreallistq(list))then
          listc%head=dtfcopy(list%head)
          listc%dbody(1:list%nl)=list%dbody(1:list%nl)
        else
          do i=0,list%nl
            listc%dbody(i)=dtfcopy(list%dbody(i))
          enddo
        endif
        return
        end subroutine

        recursive subroutine tfgetllstkall_dlist(list)
        implicit none
        type (sad_dlist) , intent(inout)::list
        type (sad_dlist),pointer :: listi
        integer*4 i,m
        logical*4 noseq
        m=list%nl
        if(iand(list%attr,lnoseqlist) .ne. 0)then
          ktastk(isp+1:isp+m)=list%body(1:m)
c     call tmov(klist(ka+1),ktastk(isp+1),m)
          isp=isp+m
          return
        endif
        noseq=.true.
        if(ktfreallistq(list))then
          dtastk(isp+1:isp+m)=list%dbody(1:m)
          isp=isp+m
        else
          do i=1,m
            isp=isp+1
            dtastk(isp)=list%dbody(i)
            if(ktfsequenceq(ktastk(isp)))then
              noseq=.false.
              isp=isp-1
              call loc_dlist(ktfaddr(list%dbody(i)),listi)
              call tfgetllstkall_dlist(listi)
            endif
          enddo
        endif
        if(noseq)then
          list%attr=ior(list%attr,lnoseqlist)
        endif
        return
        end subroutine

        recursive subroutine tfgetllstkall_rlist(list)
        implicit none
        type (sad_rlist) , intent(inout)::list
        type (sad_rlist),pointer :: listi
        integer*4 i,m
        logical*4 noseq
        m=list%nl
        if(iand(list%attr,lnoseqlist) .ne. 0)then
          ktastk(isp+1:isp+m)=list%body(1:m)
c     call tmov(klist(ka+1),ktastk(isp+1),m)
          isp=isp+m
          return
        endif
        noseq=.true.
        if(ktfreallistq(list))then
          dtastk(isp+1:isp+m)=list%dbody(1:m)
          isp=isp+m
        else
          do i=1,m
            isp=isp+1
            dtastk(isp)=list%dbody(i)
            if(ktfsequenceq(dtastk(isp)))then
              noseq=.false.
              isp=isp-1
              call loc_rlist(ktfaddr(list%dbody(i)),listi)
              call tfgetllstkall_rlist(listi)
            endif
          enddo
        endif
        if(noseq)then
          list%attr=ior(list%attr,lnoseqlist)
        endif
        return
        end subroutine

        logical*4 function tfonstackq(ka)
        implicit none
        integer*8 , intent(in)::ka
        tfonstackq=ka .ge. isporg+ispbase
     $       .and. ka .le. isporg+ivstkoffset*2+ispbase
        return
        end function

        subroutine tfgetdefargp_dlist(kl,kas,kp,ev,irtc)
        implicit none
        type (sad_dlist) , intent(inout)::kl
        type (sad_descriptor) kv
        integer*8 , intent(in)::kas
        integer*8 kp
        integer*4 isp0
        integer*4 , intent(out)::irtc
        logical*4 , intent(in)::ev
        isp=isp+1
        isp0=isp
        dtastk(isp)=kl%head
        call tfgetllstkall(kl)
        call tfdeval(isp0,kas,kv,1,.true.,ev,irtc)
        kp=ktfaddrd(kv)
        isp=isp0-1
        return
        end subroutine

        subroutine tfmakerulestk_dd(ks,kx)
        implicit none
        type (sad_descriptor) , intent(in)::kx,ks
        type (sad_dlist), pointer :: kl1
        isp=isp+1
        ktastk(isp)=ktflist+ktadaloc(-1,2,kl1)
        kl1%head%k=ktfoper+mtfrule
        kl1%dbody(1)=dtfcopy1(ks)
        kl1%dbody(2)=dtfcopy(kx)
        return
        end

        subroutine tfmakerulestk_dr(ks,x)
        implicit none
        type (sad_descriptor) , intent(in)::ks
        type (sad_dlist), pointer :: kl1
        real*8 , intent(in)::x
        isp=isp+1
        ktastk(isp)=ktflist+ktadaloc(-1,2,kl1)
        kl1%head%k=ktfoper+mtfrule
        kl1%dbody(1)=dtfcopy1(ks)
        kl1%rbody(2)=x
        return
        end

        recursive subroutine tfrebuildl(kl,klx,rep)
        implicit none
        type (sad_dlist), target, intent(inout) :: kl
        type (sad_dlist), pointer :: kli,klxi
        type (sad_dlist), pointer, intent(out) :: klx
        integer*8 kax
        integer*4 i,isp0
        logical*4, intent(out):: rep
        logical*4 rep1
        rep=.false.
        if(iand(kl%attr,ktoberebuilt) .eq. 0)then
          klx=>kl
          return
        endif
        kl%attr=kl%attr-ktoberebuilt
        if(ktfreallistq(kl))then
          klx=>kl
          return
        endif
        isp=isp+1
        isp0=isp
        dtastk(isp)=kl%head
        do i=1,kl%nl
          if(ktfsequenceq(kl%body(i),kli))then
            call tfgetllstkall(kli)
            rep=.true.
          elseif(ktflistq(kl%body(i),kli))then
            call tfrebuildl(kli,klxi,rep1)
            isp=isp+1
            if(rep1)then
              dtastk(isp)=sad_descr(klxi)
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
        type (sad_descriptor), intent(in):: k
        type (sad_dlist), pointer, intent(out) :: kl
        type (sad_dlist), pointer :: kli
        integer*4 i
        integer*4 , intent(out)::n,m
        logical*4 , intent(out)::cmplm,realm,vec
        n=0
        m=0
        realm=.false.
        cmplm=.false.
        vec=.false.
        if(ktflistq(k,kl))then
          if(kl%head%k .eq. ktfoper+mtflist)then
            n=kl%nl
            if(ktfnonreallistqo(kl))then
              do i=1,n
                if(ktfnonlistq(kl%dbody(i)))then
                  return
                endif
                call descr_sad(kl%dbody(i),kli)
                if(kli%head%k .ne. ktfoper+mtflist)then
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
        type (sad_dlist), pointer :: klx
        type (sad_rlist), pointer :: klr,klri
        integer*4 , intent(in)::n,m,nd
        integer*4 i
        logical*4 , intent(in)::trans
        real*8 , intent(in)::a(nd,m)
        if(n .eq. 0)then
          kx=kxavaloc(-1,m,klr)
          klr%rbody(1:m)=a(1:m,1)
          klr%attr=ior(klr%attr,lconstlist)
        else
          if(trans)then
            kx=kxadaloc(-1,m,klx)
            do i=1,m
              ki=kxavaloc(0,n,klri)
              klri%rbody(1:n)=a(1:n,i)
              klri%attr=ior(lconstlist,klri%attr)
              klx%dbody(i)=ki
            enddo
          else
            kx=kxadaloc(-1,n,klx)
            do i=1,n
              ki=kxavaloc(0,m,klri)
              klri%rbody(1:m)=a(i,:)
              klri%attr=ior(lconstlist,klri%attr)
              klx%dbody(i)=ki
            enddo
          endif
          klx%attr=ior(klx%attr,lconstlist)
        endif
        kxm2l=kx
        return
        end function

        type (sad_descriptor) function kxcopylist(k)
        implicit none
        type (sad_descriptor) , intent(in)::k
        type (sad_dlist), pointer :: kl,klx
        integer*4 m,i
        call descr_sad(k,kl)
        m=kl%nl
        kxcopylist=kxaaloc(-1,m,klx)
        if(ktfreallistq(kl))then
          klx%head=dtfcopy(kl%head)
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
        integer*4 , intent(in)::n0
        integer*4 n,l,n1,ls,ifrac
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

        subroutine resetnan(a,xl)
        implicit none
        real*8, intent(in), optional:: xl
        real*8, intent(inout) ::a(:)
        real*8 x
        integer*4 i
        if(present(xl))then
          x=xl
        else
          x=0.d0
        endif
        do i=1,size(a)
          if(ktfenanq(a(i)))then
            a(i)=x
          endif
        enddo
        return
        end subroutine

        subroutine limitnan(a,xl,xh,xn)
        implicit none
        real*8, intent(in), optional :: xn
        real*8, intent(in) :: xl,xh
        real*8 , intent(inout)::a(:)
        real*8 x
        integer*4 i
        if(present(xn))then
          x=xn
        else
          x=xh
        endif
        do i=1,size(a)
          if(ktfenanq(a(i)))then
            a(i)=x
          else
            a(i)=max(xl,min(xh,a(i)))
          endif
        enddo
        return
        end subroutine

      end module
