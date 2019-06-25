c     Don't confuse, Emacs. This is -*- fortran -*- mode!
c   Obsolete, moved to module tfcbk, 9/1/2012 by K.O.
c$$$      include 'inc/MACCBK.inc'
c$$$      integer*4 maxgeneration,maxlevele,nsymhash,nslots
c$$$      parameter (maxgeneration=2**30,maxlevele=4096,nsymhash=2047,
c$$$     $     nslots=32)
c$$$      integer*4 maxlbuf
c$$$      parameter (maxlbuf=32768)
c$$$      real*8 dinfinity,dnotanumber
c$$$      integer*8
c$$$     $     itfcontroot,itfcontext,itfcontextpath,itflocal,
c$$$     $     kxeof,kxfailed,iaxhold,iaximmediate,
c$$$     $     iaxline,kxliteral,kxnull,kxnulll,kxnulls,
c$$$     $     iaxout,iaxpriority,iaxschar,
c$$$     $     iaxslotnull,iaxslotpart,iaxslotseqnull,
c$$$     $     kxvect,kxvect1,iavpw,
c$$$     $     kerror,ierrorf,ierrorgen,ierrorprint,ierrorth,
c$$$     $     ierrorexp,ifunbase,initmessage
c$$$      integer*4 
c$$$     $     levelcompile,levele,levelp,lgeneration,ltrace,
c$$$     $     modethrow,iordless
c$$$      common /tffvp/
c$$$     $     itfcontroot,itfcontext,itfcontextpath,itflocal,
c$$$     $     kxeof,kxfailed,iaxhold,iaximmediate,
c$$$     $     iaxline,kxliteral,kxnull,kxnulll,kxnulls,
c$$$     $     iaxout,iaxpriority,iaxschar,
c$$$     $     iaxslotnull,iaxslotpart,iaxslotseqnull,
c$$$     $     kxvect,kxvect1,iavpw,
c$$$     $     kerror,ierrorf,ierrorgen,ierrorprint,ierrorth,
c$$$     $     ierrorexp,ifunbase,initmessage,
c$$$     $     dinfinity,dnotanumber,
c$$$     $     levelcompile,levele,levelp,lgeneration,ltrace,
c$$$     $     modethrow,iordless
