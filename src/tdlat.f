      subroutine tdlat(latt,ls,le,pos,posa,
     1                                icomp,patt,exist,lfnd)
      use tfstk
      use ffs
      use tffitcode
      use ffs_pointer, only:pnamec,lpnamec,idtypec,idvalc
      implicit real*8 (a-h,o-z)
      parameter (dh=.08,qh=.8,bh=.4,sh=.6,ch=.3,oh=.2,width=.2d0)
      real*4 xc,yc,dx,dy,xs
      integer*8 latt(nlat)
      dimension pos(nlat),posa(nlat),icomp(nlat)
      dimension height(7)
      logical tmatch,exist,fin
      character*8 name
      character*(*) patt
      data height /dh,bh,dh,qh,dh,sh,ch/
      exist=.false.
      if(ls.gt.le) then
        call tmov(pos(ls),posa(ls),nlat-ls+1)
        do 20 i=1,min(le+1,nlat)
20        posa(i)=pos(i)+pos(nlat)
      else
        call tmov(pos(ls),posa(ls),le-ls+1)
      endif
      xlen=posa(le)-posa(ls)
      write(lfnd,*)
     1'SET WINDOW X 2.5 11.0 Y 1.30 2.2'
      write(lfnd,*)
     1'SET LIMITS X ',sngl(posa(ls)),sngl(posa(le)),' Y -1 1'
      write(lfnd,*)
     1'SET OUTLINE ALL OFF;SET TICK ALL OFF;SET LABEL ALL OFF'
      write(lfnd,*)
     1'SET TITLE SIZE -1.6'
      xscale=8.5d0/xlen
      yscale=0.90d0/2.d0
      pa=posa(ls)
      ps=-width/xscale
      istart=ls
      if(ls.gt.le) then
        iend=nlat-1
      else
        iend=le-1
      endif
      do i=istart,iend
        if(idtypec(i) .eq. icMARK)then
          posa(i)=posa(i-1)
        endif
      enddo
      fin=.false.
11    do 10 i=istart,iend
        id=idtypec(i)
        if(id .ge. 41)then
          go to 12
        endif
        if(id .eq. 19)then
          yc=0.d0
          dy=height(2)*2.d0*yscale
        elseif(id .eq. 31)then
          yc=0.d0
          dy=height(7)*2.d0*yscale
        elseif(id .eq. icmult)then
          yc=sign(height(4)-oh,rlist(latt(i)+13))*.5d0
          dy=(height(4)+oh)*yscale
        elseif(icomp(i) .eq. 0 .or. id .eq. 1 .or.
     1          id .gt. icdodeca)then
          yc=0.d0
          dy=dh*2.d0*yscale
        else
          yc=sign(height(min(id,6))-oh,rlist(latt(i)+2))*.5d0
          dy=(height(min(id,6))+oh)*yscale
        endif
 12     if(kytbl(kwL,id) .gt. 0)then
          al=rlist(latt(i)+kytbl(kwL,id))
        else
          al=0.d0
        endif
        xc=pa+al*.5d0
        dx=al*xscale
        if( tmatch(pnamec(i),patt) )then
          exist=.true.
          ln=min(8,lene(pnamec(i)))
          name=' '
          name(9-ln:8)=pnamec(i)
c          xs=max(xc,REAL(ps))
          xs=REAL(ps)
          xs=max(xc,xs)
          if(xs .lt. posa(i+1)+width/xscale)then
            write(lfnd,9002)
     1      'TITLE ',xs,' 0.4  XDATA ANGLE 90 ''',name,''''
9002        format(1x,a,1p,g14.6,a,a,a)
            ps=xs+width/xscale
          endif
        endif
        if(id .ne. icMARK)then
          pa=max(pa+al,posa(i+1))
c     write(*,*)'tdlat ',pnamec(i),xc,yc,dx,dy
          write(lfnd,9001)
     1         'BOX ',xc,yc,' DATA SIZE ',dx,dy
 9001     format(1x,a,1p,2g14.6,a,2g14.6)
        endif
10    continue
      if(fin) return
      if(le.lt.ls) then
        istart=1
        iend=le-1
        fin=.true.
        goto 11
      endif
      return
      end

