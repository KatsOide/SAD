      subroutine corfree(newcor,nstr,nmon,ipstr,ipestr,ipmon,ipemon)
      use tfstk
      use ffs
      use tffitcode
c      Release memories
      parameter (kfiles=3)
      parameter (kstack=3)
      integer*8 isa,isb,isc,ipyv,ippv,ipobs,ipiv,ipbckup1,ipemibuf,
     $     ipmon,ipstr,ipestr,ipemon,k,ip,ipdefv,itcoupste
      integer*8 iifnam,iidata
      common /mcfiles/icomf,nof(kfiles),iifnam(kfiles),iidata(kfiles),
     1                memax(kfiles)
      common /mcstack/icoma,ipnt(kstack),iistck(kstack),mstkmx(kstack)
      common /pstatis/na,nb,nc,isa,isb,isc
      common /cordefv/ ipdefv,idefvar
      include 'inc/common.inc'
      include 'inc/coroper.inc'
      common /corbtune/ipyv,ippv,ipobs,ipiv,ipbckup1,ipemibuf,nemitd,
     $     nememax,iter
      common /corcoup/ itcoupste
      newcor=0
      if(ibckup.ne.0) call tfree(ibckup)
      if(itmon.ne.0) then
        call tfree(ipmon)
        call tfree(ipemon)
        nmon=0
        nmona=0
        nmonact=0
        ipmon=0
        ipemon=0
        itmon=0
        itemon=0
      endif
      if(itstr.ne.0) then
        call tfree(ipstr)
        call tfree(ipestr)
        nstr=0
        nstra=0
        ipstr=0
        ipestr=0
        itstr=0
        itestr=0
      endif
      icomf=0
      do 11 i=1,kfiles
        do 10 j=1,nof(i)
          k=ilist(mod(j-1,2)+1,iifnam(i)+(j-1)/2)
          if(k.ne.0) call tfree(k)
          k=ilist(mod(j-1,2)+1,iidata(i)+(j-1)/2)
          if(k.ne.0) call tfree(k)
   10   continue
        if(iifnam(i).ne.0) call tfree(iifnam(i))
        if(iidata(i).ne.0) call tfree(iidata(i))
   11 continue
      icoma=0
      do 13 i=1,kstack
        do 12 j=1,ipnt(i)
          k=ilist(mod(j-1,2)+1,iistck(i)+(j-1)/2)
          if(k.ne.0) call tfree(k)
   12   continue
        if(iistck(i).ne.0) call tfree(int8(iistck(i)))
   13 continue
      na=0
      nb=0
      nc=0
      if(istope.ne.0) call tfree(istope)
      istope=0
      nstope=0
      if(ipdefv.ne.0) then
        do 14 i=1,idefvar
          ip=ilist(mod(i-1,2)+1,ipdefv+(i-1)/2) 
          if(ip.ne.0) call tfree(ip)
 14     continue
        call tfree(ipdefv)
      endif
      if(ipemibuf.ne.0) then
        call tfree(ipemibuf)
        nemitd=0
      endif
      if(ipyv.ne.0) call tfree(ipyv)
      if(ippv.ne.0) call tfree(ippv)
      if(ipobs.ne.0) call tfree(ipobs)
      if(ipiv.ne.0) call tfree(ipiv)
      if(ipbckup1.ne.0) call tfree(ipbckup1)
      if(itcoupste.ne.0) call tfree(itcoupste)
      return
      end
