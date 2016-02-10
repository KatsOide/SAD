      subroutine tspini(iparam,kptbl,spect)
      use tfstk
      implicit none
      include 'inc/TMACRO1.inc'
      integer, parameter :: nkptbl = 6
      integer*4 iparam,kptbl(np0,nkptbl),nexp2(16),i,j,
     $     lp,lpa,itype,ip,kp,nt,nd1,m,italoc,ispp,nd
      real*8 sp1,c1,phi
      logical*4 spect
      data nexp2  /1,2,4,8,16,32,64,128,256,512,1024,2048,4096,8192,
     1             16382,32768/

      kptbl(1:np0,1:nkptbl)=0

      do i=1,np0
        kptbl(i,1)=i
        kptbl(i,2)=i
      enddo

      if(.not. spect)then
        return
      endif

      lspect=ilist(2,iparam)
      nspect=ilist(1,lspect)-1
      do j=1,nspect
        lp=italoc(7)
        lpa=ilist(2,lspect+j)
        itype=ilist(1,lpa)
        ip=ilist(1,lpa+1)
        if(kptbl(ip,3) .gt. 0)then
          if(itype .lt. 0)then
            kptbl(ip,3)=-kptbl(ip,3)
          endif
        elseif(kptbl(ip,3) .eq. 0)then
          kp=italoc(12)
          if(itype .gt. 0)then
            kptbl(ip,3)=kp
          else
            kptbl(ip,3)=-kp
          endif
        endif
        ilist(2,lpa)=lp
        nt=nturn-ilist(1,lpa+2)+2
        nd1=min(nint(rlist(lpa+3)),nt)
        do i=1,16
          if(nexp2(i) .le. nd1)then
            nd=nexp2(i)
          endif
        enddo
        ilist(1,lp)=nd
        sp1=1.d0/nd
        c1=mod(rlist(lpa+4)+sp1/2.d0,1.d0)
        phi=pi*(2.d0*c1-1.d0)
        rlist(lp+3)=cos(phi)
        rlist(lp+4)=sin(phi)
        rlist(lpa+4)=0.5d0
        rlist(lpa+3)=1.d0
        ilist(1,lp+1)=0
        ilist(1,lp+2)=0
        ilist(1,lp+5)=italoc(nd*2)
        call tclr(rlist(ilist(1,lp+5)),nd*2)
        m=nt/nd
        rlist(lp+6)=1.d0/sqrt(dble(nd*m))
        ispp=italoc(nd)-1
        ilist(1,lspect+j)=ispp
        call tclr(rlist(ispp+1),nd)
      enddo
      return
      end

      subroutine tplini(iparam,kptbl)
      use tfstk
      implicit none
      include 'inc/TMACRO1.inc'
      integer*4 iparam,kptbl(np0,6),j,itype,kp,
     $     lpa,ip,nt,ispp,italoc
      lplot=ilist(2,iparam)
      nplot=ilist(1,lplot)-1
      do j=1,nplot
        lpa=ilist(2,lplot+j)
        itype=ilist(1,lpa)
        ip=ilist(1,lpa+1)
        if(kptbl(ip,3) .gt. 0)then
          if(itype .lt. 0)then
            kptbl(ip,3)=-kptbl(ip,3)
          endif
        elseif(kptbl(ip,3) .eq. 0)then
          kp=italoc(12)
          if(itype .gt. 0)then
            kptbl(ip,3)=kp
          else
            kptbl(ip,3)=-kp
          endif
        endif
        nt=nturn-ilist(1,lpa+2)+2
        ispp=italoc(nt*2)-1
        call tclr(rlist(ispp+1),nt*2)
        ilist(1,lplot+j)=ispp
      enddo
      return
      end
