      subroutine  corinit(newcor,nstr,nmon,ipstr,ipestr,ipmon,
     1                    ipemon)
      use ffs
      use tffitcode
c ---- Initialize correction programs
      implicit real*8(a-h,o-z)
      parameter (kfiles=3)
      parameter (kstack=3)
      logical*4 start
      common /mcfiles/icomf,nof(kfiles),iifnam(kfiles),iidata(kfiles),
     1                memax(kfiles)
      common /mcstack/icoma,ipnt(kstack),iistck(kstack),mstkmx(kstack)
      common /pstatis/na,nb,nc,isa,isb,isc,start
      include 'inc/common.inc'
      include 'inc/coroper.inc'
      common /cordefv/ ipvbmp,nvbmp
      common /corbtune/ipyv,ippv,ipobs,ipiv,ipbckup1,ipemibuf,nemitd,
     $     nememax,iter
      common /corcoup/ itcoupste
c
      newcor=0
      ipmon=0
      ipemon=0
      ipstr=0
      ipestr=0
      itmon=0
      itemon=0
      itstr=0
      itestr=0
      nmon=0
      nmona=0
      nmonact=0
      nstr=0
      nstra=0
      eptsol=1d-8
      do 10 i=1,4
        errflg(i)=.false.
        errval(i)=0d0
   10 continue
      dpshft=1d-3
      ibckup=0
      icomf=0
      do 20 i=1,kfiles
        nof(i)=0
        memax(i)=0
        iifnam(i)=0
        iidata(i)=0
   20 continue
      icoma=0
      do 21 i=1,kstack
        ipnt(i)=0
        iistck(i)=0
        mstkmx(i)=0
   21 continue
      na=0
      nb=0
      nc=0
      start=.false.
      istope=0
      nstope=0
      ipvbmp=0
      nvbmp=0
      ipyv=0
      ippv=0
      ipobs=0
      ipiv=0
      ipbckup1=0
      ipemibuf=0
      nemitd=0
      nememax=0
      iter=0
      itcoupste=0
      return
      end
