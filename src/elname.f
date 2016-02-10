      subroutine elname(latt,i,mult,name)
      implicit none
      integer*4 latt(2,*),i,mult(*)
      character*(*) name
      call elname1(latt,i,mult,name,.false.)
      return
      end

      subroutine elnameK(latt,i,mult,name)
      implicit none
      integer*4 latt(2,*),i,mult(*)
      character*(*) name
      call elname1(latt,i,mult,name,.true.)
      return
      end

      subroutine elname1(latt,i,mult,name,comp)
      use tfstk
      use ffs
      use tffitcode
      implicit none 
      integer*4 latt(2,nlat),i,mult(nlat)
      character*(*) name
      logical*4 comp
c
      integer*4 lenw
c
      integer ltyp,idx,l
      character*(MAXPNAME) buff
c
      if(i .eq. nlat)then
        name='$$$'
      else
        name=pname(latt(1,i))
        idx=max(0,abs(mult(i)))
        if(comp .and. idx .eq. 0)then
c     Case: Force append suffix number by comp flag
          ltyp=idtype(latt(1,i))
          if(ltyp .gt. icNULL .and. ltyp .lt. icMXEL)then
            idx=kytbl(kwINDX,ltyp)
            if(idx .ne. 0)then
              idx=max(1,int(rlist(idval(latt(1,i))+idx)))
            endif
          else
            idx=0
          endif
c     idx  0: no defaults suffix number
c     *   >0: suffix number starts `idx'
          l=lenw(name)
          if(idx .gt. 1)then
            write(buff,*)idx
            name(l+1:)='.'//adjustl(buff)
          else
            name(l+1:)='.1'
          endif
        elseif(idx .ne. 0)then
c     Case: Append suffix number for multiple elements
          l=lenw(name)
          if(idx .gt. 1)then
            write(buff,*)idx
            name(l+1:)='.'//adjustl(buff)
          else
            name(l+1:)='.1'
          endif
        endif
      endif
      return
      end
