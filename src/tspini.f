      subroutine tspini(iparam,kptbl,spect)
      use tfstk
      use tmacro
      implicit none
      integer, parameter :: nkptbl = 6
      integer*8 iparam,lp
      integer*4 kptbl(np0,nkptbl),nexp2(16),i,j,
     $     lpa,itype,ip,kp,nt,nd1,m,ispp,nd
      real*8 sp1,c1,phi
      logical*4 spect
      data nexp2  /1,2,4,8,16,32,64,128,256,512,1024,2048,4096,8192,
     1             16382,32768/

      kptbl(1:np0,1:nkptbl)=0

      do i=1,np0
        kptbl(i,1)=i
        kptbl(i,2)=i
      enddo

c      if(spect)then
c        write(*,*)'Use functions within FFS.'
c      endif
      return
c$$$
c$$$      lspect=klist(iparam)
c$$$      nspect=ilist(1,lspect)-1
c$$$      do j=1,nspect
c$$$        lp=ktaloc(7)
c$$$        lpa=ilist(2,lspect+j)
c$$$        itype=ilist(1,lpa)
c$$$        ip=ilist(1,lpa+1)
c$$$        if(kptbl(ip,3) .gt. 0)then
c$$$          if(itype .lt. 0)then
c$$$            kptbl(ip,3)=-kptbl(ip,3)
c$$$          endif
c$$$        elseif(kptbl(ip,3) .eq. 0)then
c$$$          kp=ktaloc(12)
c$$$          if(itype .gt. 0)then
c$$$            kptbl(ip,3)=kp
c$$$          else
c$$$            kptbl(ip,3)=-kp
c$$$          endif
c$$$        endif
c$$$        ilist(2,lpa)=lp
c$$$        nt=nturn-ilist(1,lpa+2)+2
c$$$        nd1=min(nint(rlist(lpa+3)),nt)
c$$$        do i=1,16
c$$$          if(nexp2(i) .le. nd1)then
c$$$            nd=nexp2(i)
c$$$          endif
c$$$        enddo
c$$$        ilist(1,lp)=nd
c$$$        sp1=1.d0/nd
c$$$        c1=mod(rlist(lpa+4)+sp1/2.d0,1.d0)
c$$$        phi=pi*(2.d0*c1-1.d0)
c$$$        rlist(lp+3)=cos(phi)
c$$$        rlist(lp+4)=sin(phi)
c$$$        rlist(lpa+4)=0.5d0
c$$$        rlist(lpa+3)=1.d0
c$$$        ilist(1,lp+1)=0
c$$$        ilist(1,lp+2)=0
c$$$        ilist(1,lp+5)=ktaloc(nd*2)
c$$$        call tclr(rlist(ilist(1,lp+5)),nd*2)
c$$$        m=nt/nd
c$$$        rlist(lp+6)=1.d0/sqrt(dble(nd*m))
c$$$        ispp=ktaloc(nd)-1
c$$$        ilist(1,lspect+j)=ispp
c$$$        call tclr(rlist(ispp+1),nd)
c$$$      enddo
c$$$      return
      end

      subroutine tplini(iparam,kptbl)
      use tfstk
      use tmacro
      implicit none
      integer*8 iparam,ispp
      integer*4 kptbl(np0,6),j,itype,kp,lpa,ip,nt
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
          kp=ktaloc(12)
          if(itype .gt. 0)then
            kptbl(ip,3)=kp
          else
            kptbl(ip,3)=-kp
          endif
        endif
        nt=nturn-ilist(1,lpa+2)+2
        ispp=ktaloc(nt*2)-1
        call tclr(rlist(ispp+1),nt*2)
        ilist(1,lplot+j)=ispp
      enddo
      return
      end
