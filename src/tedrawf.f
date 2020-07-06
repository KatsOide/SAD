      subroutine tedrawf(lfno,bc,bs,hc,hs,
     $     sige,anuz,mphi2,ndp)
      implicit none
      integer*4 mphi2,ndp,lfno,lene
      real*8 bc(10,mphi2,ndp),bs(10,mphi2,ndp),
     $     hc(4,mphi2,ndp),hs(4,mphi2,ndp),
     $     sige,anuz
      character*32 buf,autofg
      write(lfno,'(a,i5,a)')'mphi2=',mphi2,';',
     $     'ndp=',ndp,';'
      buf=autofg(sige,'12.9')
      call mathexp(buf)
      write(lfno,'(a)')'sige='//buf(1:lene(buf))//';'
      buf=autofg(anuz,'12.9')
      call mathexp(buf)
      write(lfno,'(a)')'nuz='//buf(1:lene(buf))//';'
      call matharray(bc,10,mphi2,ndp,mphi2,'bc=',lfno)
      call matharray(bs,10,mphi2,ndp,mphi2,'bs=',lfno)
      call matharray(hc,4,mphi2,ndp,mphi2,'hc=',lfno)
      call matharray(hs,4,mphi2,ndp,mphi2,'hs=',lfno)
      return
      end

      subroutine matharray(bc,n,mphi,ndp,mphi2,s,lfno)
      implicit none
      integer*4 n,mphi2,mphi,ndp,lfno,lene,l
      real*8 bc(n,mphi,ndp)
      integer*4 i,m,j
      character*(*) s
      character*64 s1
      l=lene(s)
      s1=s(:l)
      write(lfno,'(a)')s1(:l)//'{{'
      do j=1,ndp
        do m=1,mphi2
          call mathinit('{')
          do i=1,n
            call mathput(lfno,bc(i,m,j))
          enddo
          if(m .eq. mphi2)then
            call mathterm(lfno,'}')
          else
            call mathterm(lfno,'},')
          endif
        enddo
        if(j .eq. ndp)then
          write(lfno,'(a)')'}};'
        else
          write(lfno,'(a)')'},{'
        endif
      enddo
      return
      end

      subroutine mathinit(s)
      implicit none
      integer*4 lfno,ip,lene
      character*80 mb
      character*(*) s
      common /math/ip,mb
      ip=lene(s)+1
      mb=s(1:ip-1)
      return
      end

      subroutine mathterm(lfno,s)
      implicit none
      character*(*) s
      integer*4 lfno,ip,lene,l
      character*80 mb
      common /math/ip,mb
      l=lene(s)
      if(ip+l .gt. 81)then
        write(lfno,'(a)')mb(1:ip-1)
        ip=1
      endif
      mb(ip:80)=s(1:l)
      write(lfno,'(a)')mb(1:ip+l-1)
      return
      end

      subroutine mathput(lfno,a)
      implicit none
      character*32 buf,autofg
      real*8 a
      integer*4 lfno,ip,l,lene
      character*80 mb
      common /math/ip,mb
      buf=autofg(a,'S12.9')
      call mathexp(buf)
      l=lene(buf)
      if(ip .eq. 1)then
        mb=buf(1:l)
        ip=l+1
      elseif(mb(ip-1:ip-1) .eq. '{')then
        mb(ip:ip+l-1)=buf(1:l)
        ip=ip+l
      elseif(ip+l .lt. 80)then
        mb(ip:ip+l)=','//buf(1:l)
        ip=ip+l+1
      else
        mb(ip:ip)=','
        write(lfno,'(a)')mb(1:ip)
        mb=buf(1:l)
        ip=l+1
      endif
      return
      end

      subroutine mathexp(buf)
      implicit none
      character*(*) buf
      character*256 buf1
      integer*4 lene,i,j
      buf1=buf
      buf=' '
      j=1
      do i=1,lene(buf1)
        if(buf1(i:i) .eq. 'E')then
          buf(j:j+3)=' 10^'
          j=j+4
        else
          buf(j:j)=buf1(i:i)
          j=j+1
        endif
      enddo
      return
      end
