      subroutine tshow(flv,latt,scale,nlist,mult,kx,irtc,ret,lfno)
      use tfstk
      use ffs
      use ffslocal, only: ffslocalv
      use tffitcode
      implicit none
      type (ffslocalv) flv
      type (sad_descriptor) kx
      type (sad_list), pointer :: klx,klxi
      integer*4 lfno,nc,i,j,kf
      integer*4 latt(2,nlat),mult(nlat),irtc
      real*8 scale(mfit1),x
      logical*4 ret
      character*8 nlist(mfit1),fun
      character*12 name,name1
      character*32 namer,name1r
      character*12 autofg,vout,sout
      irtc=0
      if(lfno .gt. 0)then
        write(lfno,'(2a)')'!    component1   component2',
     $       '   fun       goal-value  np     scale'
      endif
      nc=0
      do i=1,flv%nfc
        if(flv%mfitp(i) .ne. 0)then
          if(lfno .gt. 0)then
            call elname(latt,flv%ifitp(i),mult,name)
            kf=flv%kfit(i)
            x=flv%fitval(i)/scale(kf)
            vout=autofg(x,'12.9')
            sout=autofg(scale(kf),'12.9')
            if(flv%ifitp(i) .ne. flv%ifitp1(i))then
              call elname(latt,flv%ifitp1(i),mult,name1)
            else
              name1=' '
            endif
            fun=nlist(kf)
            if(flv%mfitp(i) .lt. 0)then
              fun(len_trim(fun)+1:)='M'
            endif
            write(lfno,9001)name,name1,fun,vout,
     $           abs(flv%mfitp(i))-1,sout
 9001       format(1x,'FIT ',a,' ',a,' ',a,' ',a,i3,' ! *',a)
          endif
          nc=nc+1
        endif
      enddo
      if(ret)then
        kx=kxadaloc(-1,nc,klx)
        j=0
        do i=1,flv%nfc
          if(flv%mfitp(i) .ne. 0)then
            call elname(latt,flv%ifitp(i),mult,namer)
            kf=flv%kfit(i)
            x=flv%fitval(i)/scale(kf)
            if(flv%ifitp(i) .ne. flv%ifitp1(i))then
              call elname(latt,flv%ifitp1(i),mult,name1r)
            else
              name1r=' '
            endif
            fun=nlist(kf)
            if(flv%mfitp(i) .lt. 0)then
              fun(len_trim(fun)+1:)='M'
            endif
            j=j+1
            klx%dbody(j)=kxadaloc(0,6,klxi)
            klxi%dbody(1)=kxsalocb(0,namer,len_trim(namer))
            klxi%dbody(2)=kxsalocb(0,name1r,len_trim(name1r))
            klxi%dbody(3)=kxsalocb(0,fun,len_trim(fun))
            klxi%rbody(4)=x
            klxi%rbody(5)=dble(abs(flv%mfitp(i))-1)
            klxi%rbody(6)=scale(kf)
          endif
        enddo
      endif
      return
      end
