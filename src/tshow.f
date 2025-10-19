      subroutine tshow(kx,irtc,ret,lfno)
      use tfstk
      use ffs
      use ffs_pointer
      use ffs_fit
      use tffitcode
      implicit none
      type (sad_descriptor) ,intent(out):: kx
      type (sad_dlist), pointer :: klx,klxi
      type (sad_rlist), pointer :: krl,krl1
      type (sad_descriptor) ,save :: k1
      data k1%k /0/
      type (sad_rlist) , pointer, save :: kl1
      integer*4 ,intent(out):: irtc
      integer*4 ,intent(in):: lfno
      integer*4 nc,i,j,kf,nc1
      real*8 x,f
      logical*4 ,intent(in):: ret
      character*8 fun
      character*12 name,name1
      character*32 namer,name1r
      character*17 autofg,vout,tfconvstr
      character*12 sout
      irtc=0
      if(lfno > 0)then
        write(lfno,'(2a)')'!    component1   component2   fun         goal-value     np     scale'
      endif
      nc=0
      if(k1%k == 0)then
        k1=kxavaloc(0,2,kl1)
      endif
      do i=1,flv%nfc
        if(flv%mfitp(i) /= 0)then
          if(lfno > 0)then
            call elname(flv%ifitp(i),name)
            kf=flv%kfit(i)
            if(ktfrealq(flv%fitval(i),f))then
              x=f/scale(kf)
              vout=autofg(x,'12.9')
            elseif(ktfreallistq(flv%fitval(i),krl))then
              if(krl%nl /= 2)then
                cycle
              endif
              kl1%rbody(1)=krl%rbody(1)/scale(kf)
              kl1%rbody(2)=krl%rbody(2)/scale(kf)              
c              call tfdebugprint(k1,'tshow',1)
              vout=tfconvstr(k1,nc1,'S7.5')
            else
              cycle
            endif
            sout=autofg(scale(kf),'12.9')
            if(flv%ifitp(i) /= flv%ifitp1(i))then
              call elname(flv%ifitp1(i),name1)
            else
              name1=' '
            endif
            fun=nlist(kf)
            if(flv%mfitp(i) .lt. 0)then
              fun(len_trim(fun)+1:)='M'
            endif
            write(lfno,9001)name,name1,fun,vout,abs(flv%mfitp(i))-1,sout
 9001       format(1x,'FIT ',a,' ',a,' ',a,' ',a,i3,' ! *',a)
          endif
          nc=nc+1
        endif
      enddo
      if(ret)then
        kx=kxadaloc(-1,nc,klx)
        j=0
        do i=1,flv%nfc
          if(flv%mfitp(i) /= 0)then
            call elname(flv%ifitp(i),namer)
            kf=flv%kfit(i)
            if(flv%ifitp(i) /= flv%ifitp1(i))then
              call elname(flv%ifitp1(i),name1r)
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
            if(ktfrealq(flv%fitval(i),f))then
              klxi%rbody(4)=flv%fitval(i)%x(1)/scale(kf)
            elseif(ktfreallistq(flv%fitval(i),krl))then
              klxi%dbody(4)=kxavaloc(0,2,krl1)
              krl1%rbody(1)=krl%rbody(1)/scale(kf)
              krl1%rbody(2)=krl%rbody(2)/scale(kf)
            else
              klxi%rbody(i)=0.d0
            endif
            klxi%rbody(5)=dble(abs(flv%mfitp(i))-1)
            klxi%rbody(6)=scale(kf)
          endif
        enddo
      endif
      return
      end
