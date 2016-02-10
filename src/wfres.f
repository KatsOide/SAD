      subroutine wfres(zbuf,zwork,lsp,
     1                 minord,lcoeff,nres,
     1                 iord,aord,ndim,f1,f2,f3)
      implicit none
      integer*4 lsp,minord,maxord,nres,iord(3,*),ndim,i,j,k,jm,i1
      integer*4 lcoeff, ip
      parameter (maxord=100)
      real*8 wtune,f1,f2,f3,aord(3,*),fmax,rgetgl1
      complex*16 zbuf(lsp,3),zwork(lsp)
      f1=wtune(zbuf(1,1),zwork,lsp)
      f2=wtune(zbuf(1,2),zwork,lsp)
      if(ndim .ge. 3)then
        f3=wtune(zbuf(1,3),zwork,lsp)
      else
        f3=0.d0
      endif
      if(nres .le. 0)then
        return
      endif
      fmax=nres/rgetgl1('ADDDENSE')/lsp
      ip=lcoeff-nres+1
      do 10 i=minord,maxord,2
        do 11 k=i,1-i,-1
          call wfbin(f1,f2,f3,k,i-abs(k),0,fmax,ip,iord)
          if(ip .gt. lcoeff)then
            go to 99
          endif
11      continue
        i1=i+1
        if(ndim .ge. 3)then
          do 20 j=1,i
            jm=i-j
            do 21 k=-jm,jm
              call wfbin(f1,f2,f3,k,jm-abs(k),j,fmax,ip,iord)
              if(ip .gt. lcoeff)then
                go to 99
              endif
21          continue
            do 31 k=jm-1,1-jm,-1
              call wfbin(f1,f2,f3,k,abs(k)-jm,j,fmax,ip,iord)
              if(ip .gt. lcoeff)then
                go to 99
              endif
31          continue
20        continue
          do 40 j=i1,1,-1
            jm=i1-j
            do 51 k=1-jm,jm-1
              call wfbin(f1,f2,f3,k,abs(k)-jm,j,fmax,ip,iord)
              if(ip .gt. lcoeff)then
                go to 99
              endif
51          continue
            do 61 k=jm,-jm,-1
              call wfbin(f1,f2,f3,k,jm-abs(k),j,fmax,ip,iord)
              if(ip .gt. lcoeff)then
                go to 99
              endif
61          continue
40        continue
        endif
        do 71 k=1-i1,i1
          call wfbin(f1,f2,f3,k,i-abs(k),0,fmax,ip,iord)
          if(ip .gt. lcoeff)then
            go to 99
          endif
71      continue
10    continue
99    do 100 i=lcoeff-nres+1,lcoeff
        aord(1,i)=iord(1,i)
        aord(2,i)=iord(2,i)
        aord(3,i)=iord(3,i)
100   continue
c     write(*,*)(iord(k,lcoeff),k=1,3)
c     write(*,'(10I5)')(k,nbin(k),k=1,nfreq)
      return
      end
