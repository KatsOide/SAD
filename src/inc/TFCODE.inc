c$$$c     Don't confuse, Emacs. This is -*- fortran -*- mode!
c$$$      real*8 xinfinity
c$$$      parameter (xinfinity=1.7976931348623157D308)
c$$$      integer*4 ntfoper,ntfreal,ntflist,ntflistr,ntfdef,ntfstkseq,
c$$$     $     ntfstring,ntfsymbol,ntfpat,ntffun,ntfstk,ntfarg
c$$$      parameter (ntfoper=0,ntfreal=1,ntflist=3,ntflistr=4,ntfstkseq=5,
c$$$     $     ntfstk=6,
c$$$     $     ntfstring=101,ntfsymbol=201,ntfpat=203,ntfarg=205,
c$$$     $     ntffun=ntfoper,ntfdef=ntfsymbol)
c$$$      integer*4 nfunif,nfunlength,nfundo,nfunmodule,nfunblock,
c$$$     $     nfununeval,nfunwith,nfunthread
c$$$      parameter (nfunif=39,nfunlength=20,nfundo=29,nfunmodule=34,
c$$$     $     nfunblock=35,nfununeval=132,nfunwith=141,nfunthread=83)
c$$$      integer*4 mtfnull,mtfneg,mtfinv,mtfplus,mtfminus,mtfmult,mtfdiv,
c$$$     $     mtfdot,mtfpower,mtfequal,mtfunequal,mtfgreater,mtfless,
c$$$     $     mtfgeq,mtfleq,mtfsame,mtfunsame,mtfnot,mtfand,mtfor,
c$$$     $     mtfconcat,mtfleftbra,mtfrightbra,mtfleftbrace,
c$$$     $     mtfrightbrace,mtfsetdelayed,mtfset,mtfcomplex,mtfleftparen,
c$$$     $     mtfrightparen,mtfcomma,mtftimes,
c$$$     $     mtfcomp,mtffun,mtfcolon,mtfrule,mtfruledelayed,mtfreplace,
c$$$     $     mtfreplacerepeated,mtfupset,mtfupsetdelayed,mtfunset,
c$$$     $     mtfpattest,mtfflag,mtfslot,mtfslotseq,mtfrevpower,mtfalt,
c$$$     $     mtflist,mtfmap,mtfmapall,mtfapply,mtfrepeated,
c$$$     $     mtfrepeatednull,mtfinequality,mtfaddto,mtfsubtractfrom,
c$$$     $     mtftimesby,mtfdivideby,mtfincrement,mtfdecrement,
c$$$     $     mtfpart,mtfatt,mtfmessagename,mtftagset,
c$$$     $     mtfleftcomment,mtfrightcomment,mtfhold,mtfend
c$$$      parameter (
c$$$     $     mtfnull=0,
c$$$     $     mtfneg=1,mtfinv=2,mtfplus=3,mtfminus=4,mtfmult=5,mtfdiv=6,
c$$$     $     mtfrevpower=7,mtfpower=8,mtfequal=9,mtfunequal=10,
c$$$     $     mtfgreater=11,
c$$$     $     mtfless=12,mtfgeq=13,mtfleq=14,mtfsame=15,mtfunsame=16,
c$$$     $     mtfnot=17,mtfand=18,mtfor=19,mtfconcat=20,mtfleftbra=21,
c$$$     $     mtfrightbra=22,mtfleftbrace=23,mtfrightbrace=24,
c$$$     $     mtfsetdelayed=25,mtfset=26,mtfcomplex=27,mtfleftparen=28,
c$$$     $     mtfrightparen=29,mtfcomma=30,mtfcomp=31,mtffun=32,
c$$$     $     mtfcolon=33,mtfrule=34,mtfruledelayed=35,mtfreplace=36,
c$$$     $     mtfreplacerepeated=37,mtfupset=38,mtfupsetdelayed=39,
c$$$     $     mtfunset=40,mtfpattest=41,mtfflag=42,mtfslot=43,
c$$$     $     mtfslotseq=44,mtfdot=45,mtfalt=46,mtfmap=47,mtfmapall=48,
c$$$     $     mtfapply=49,mtfrepeated=50,mtfrepeatednull=51,
c$$$     $     mtfinequality=52,mtfaddto=53,mtfsubtractfrom=54,
c$$$     $     mtftimesby=55,mtfdivideby=56,
c$$$     $     mtfincrement=57,mtfdecrement=58,
c$$$     $     mtfpart=59,mtfatt=60,mtfmessagename=61,mtftagset=62,
c$$$     $     mtfleftcomment=63,mtfrightcomment=64,mtfhold=65,
c$$$     $     mtfend=66,
c$$$     $     mtflist=mtfleftbrace,mtftimes=mtfmult)
c$$$      integer*4 mtfnopc
c$$$      parameter (mtfnopc=mtfend)
c$$$      integer*4 lsimplepat,lsimplepatlist,lconstlist,lnoconstlist,
c$$$     $     larglist,lnopatarg,lmemberlist,lnodefsymbol,lnoseqlist,
c$$$     $     lnonreallist
c$$$      parameter (lsimplepat=1,lsimplepatlist=2,larglist=4,
c$$$     $     lconstlist=8,lnoconstlist=16,lnopatarg=32,lmemberlist=64,
c$$$     $     lnodefsymbol=128,lnoseqlist=256,lnonreallist=512)
c$$$      integer*4 kconstlist,knoconstlist,kfixedarg,knofixedarg,
c$$$     $     kallnofixedarg,knopatarg,kpatarg,
c$$$     $     kconstarg,knoconstarg,kseqarg,knoseqarg
c$$$      parameter (knopatarg=lnonreallist*2,
c$$$     $     kpatarg=lnonreallist*4,
c$$$     $     kconstlist=lnonreallist*8,
c$$$     $     knoconstlist=lnonreallist*16,
c$$$     $     kfixedarg=lnonreallist*32,
c$$$     $     knofixedarg=lnonreallist*64,
c$$$     $     kallnofixedarg=lnonreallist*128,
c$$$     $     kconstarg=lnonreallist*256,
c$$$     $     knoconstarg=lnonreallist*512,
c$$$     $     kseqarg=lnonreallist*1024,
c$$$     $     knoseqarg=lnonreallist*2048)
c$$$      integer*4 iattrholdfirst,iattrholdrest,iattrholdall,
c$$$     $     iattrconstant,iattrimmediate,iattrorderless,
c$$$     $     iattrdynamic,iattrprotected,iattrnumeric
c$$$      parameter (iattrholdfirst=1,iattrholdrest=2,
c$$$     $     iattrholdall=iattrholdfirst+iattrholdrest,
c$$$     $     iattrimmediate=4,iattrnumeric=8,
c$$$     $     iattrorderless=16,iattrdynamic=32,
c$$$     $     iattrprotected=64,iattrconstant=128)
c$$$      real*8 rtfnull
c$$$      parameter (rtfnull=0.d0)
c$$$      integer*8 ktfoper,ktflist,ktfstring,ktfsymbol,ktfpat,ktfobj,
c$$$     $     ktfmask,ktamask,ktrmask,ktfnull,ktfnr,ktfref,ktfother,
c$$$     $     ktomask,ktftrue,ktfnan
c$$$      parameter (
c$$$     $     ktfnull  =int8(z'fff0000000000000'),
c$$$     $     ktfother =int8(z'fff2000000000000'),
c$$$     $     ktfnr    =int8(z'7ff2000000000000'),
c$$$     $     ktfoper  =int8(z'fff6000000000000'),
c$$$     $     ktfref   =int8(z'fffa000000000000'),
c$$$     $     ktfobj   =int8(z'7ff2000000000000'),
c$$$     $     ktflist  =int8(z'7ff2000000000000'),
c$$$     $     ktfpat   =int8(z'7ff6000000000000'),
c$$$     $     ktfstring=int8(z'7ffa000000000000'),
c$$$     $     ktfsymbol=int8(z'7ffe000000000000'),
c$$$     $     ktomask  =int8(z'fff2000000000000'),
c$$$     $     ktrmask  =int8(z'7ff2000000000000'),
c$$$     $     ktfmask  =int8(z'fffe000000000000'),
c$$$     $     ktamask  =int8(z'0001ffffffffffff'),
c$$$     $     ktftrue  =int8(z'3ff0000000000000'),
c$$$     $     ktfnan   =int8(z'fff8000000000000'))
c$$$
c$$$      type sad_symbol
c$$$      sequence
c$$$      integer*4 attr,override
c$$$      integer*8 alloc
c$$$      integer*4 ref,gen
c$$$      integer*8 loc
c$$$      end type
c$$$
c$$$      type sad_string
c$$$      sequence
c$$$      integer*4 len,override
c$$$      integer*8 alloc
c$$$      integer*4 ref,gen
c$$$      integer*4 nch,nc
c$$$      character*(2**31-1) str
c$$$      end type
