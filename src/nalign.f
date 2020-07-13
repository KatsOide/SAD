      subroutine nalign(latt,mult,master,word,exist)
      use tfstk
      use ffs
      use tffitcode
      use sad_main
      use ffs_pointer, only:elatt
      integer*8 latt(nlat),k,jd
      dimension mult(nlat),master(nlat)
      integer*4 master
      character*(*) word
      character*20  name
      logical exist,temat
      integer*4 ie,is,i,j,id
      real*8 dx,dy
c
      exist=.false.
        do 1010 i=1,nlat-1
          if(temat(i,name,word))then
            is=i
            go to 1020
          end if
1010    continue
        return
1020    continue
        call getwdl(word)
        do 1100 i=1,nlat-1
          if(temat(i,name,word))then
            exist=.true.
            ie=i
            go to 1120
          end if
1100    continue
        ie=nlat
1120    continue
        if( is .gt. ie) then
          i=is
          is=ie
          ie=i
        end if
c
        if (master(is) .lt. 0) then is=-master(is)
        if (master(ie) .gt. 0) ie=master(ie)
c
c       call elname(is,name)
c       write(lfno,*)' nalign start', name
c       call elname(ie,name)
c       write(lfno,*)' nalign end  ', name
c
        j=idcomp(elatt,is)
          k=latt(is)
          id=idtype(j)
c
          if(id .eq. 2) then
            dx=rlist(k+9)
            dy=rlist(k+10)
          else if(id .eq. 4 .or. id .eq. 6 .or. id .eq.8) then
            dx=rlist(k+5)
            dy=rlist(k+6)
          else if(id .eq. 20) then
            jd=idval(j)
            dx=rlist(jd+3)-rlist(k+3)
            dy=rlist(jd+4)-rlist(k+4)
          else
            dx=0.d0
            dy=0.d0
          end if
c
c       write (lfno,*)' align :',dx,dy,is,ie
c
        do 3010 i=is+1,ie
          j=idcomp(elatt,i)
          k=latt(i)
          id=idtype(j)
          go to (3010,3120,3010,3140,3010,3140,3010,3140,3010,3010,
     &           3010,3010,3010,3010,3010,3010,3010,3010,3010,3160),id
            go to 3010
3120        continue
            rlist(k+9)  =dx
            rlist(k+10) =dy
            go to 3010
3140        continue
            rlist(k+5)=dx
            rlist(k+6)=dy
            go to 3010
3160        continue
            jd=idval(j)
            rlist(k+3)=rlist(jd+3)-dx
            rlist(k+4)=rlist(jd+4)-dy
3010    continue
c
c      write(lfno,*) ' align ends'
c
      return
c
      end
