      subroutine tfattr(word,lfno,exist,kx,irtc,ret)
      use tfstk
      use ffs
      use ffs_pointer
      use tffitcode
      implicit none
      type (sad_descriptor) kx
      type (sad_dlist), pointer :: klxi
      integer*4 lfno,i,j
      integer*4 irtc,isp1,lenw
      real*8 v
      character*(*) word
      character*(MAXPNAME+16) namc,name
      character*10 autofg
      character*8 tfkwrd,key
      logical*4 exist,exist1,all,temat,ret
      write(lfno,'(a)')
     $     'Element     Keyword    Value     Mimimum  Maximum'
     1           //'  Couple         Coefficient'
      all=.false.
      exist=.false.
      isp1=isp
 1    exist1=.false.
      call getwdl(word)
      if(word .eq. ' ')then
        if(exist)then
          go to 9000
        else
          all=.true.
        endif
      endif
 2    LOOP_J: do j=1,nlat-1
        if(all)then
          call elname(j,name)
        elseif(temat(j,name,word))then
        else
          cycle LOOP_J
        endif
        exist1=.true.
        i=iele1(icomp(j))
        if(icomp(j) .eq. j .or. nelvx(iele1(j))%klp .eq. j .or.
     $       icomp(j) .ne. nelvx(iele1(j))%klp .and.
     $       icomp(j) .ne. icomp(nelvx(iele1(j))%klp))then
          call elname(icomp(j),namc)
          if(nelvx(i)%ival .eq. 0)then
            v=0.d0
          else
            v=rlist(latt(j)+nelvx(i)%ival)/errk(1,j)
          endif
          key=tfkwrd(idtypec(j),nelvx(i)%ival)
          if(ret)then
            dtastk(isp)=kxadaloc(-1,7,klxi)
            klxi%dbody(1)=kxsalocb(0,name,lenw(name))
            klxi%dbody(2)=kxsalocb(0,key,lenw(key))
            klxi%rbody(3)=v
            klxi%rbody(4)=nelvx(i)%vlim(1)
            klxi%rbody(5)=nelvx(i)%vlim(2)
            klxi%dbody(6)=kxsalocb(0,namc,lenw(namc))
            klxi%rbody(7)=couple(j)
            isp=isp+1
          endif
          if(nelvx(i)%klp .eq. j)then
            namc='<--'
          endif
          write(lfno,'(5a,1x,2a)')name(1:12),
     $         key,autofg(v,'10.6'),
     $         autofg(nelvx(i)%vlim(1),'10.6'),
     $         autofg(nelvx(i)%vlim(2),'10.6'),
     1         namc(1:12),autofg(couple(j),'10.6')
        endif
      enddo LOOP_J
      if(.not. all)then
        if(exist)then
          if(exist1)then
            go to 1
          else
            exist=.false.
          endif
        else
          if(exist1)then
            exist=.true.
          else
            all=.true.
            go to 2
          endif
          go to 1
        endif
      endif
 9000 if(ret)then
        kx=kxmakelist(isp1)
        isp=isp1
        irtc=0
      endif
      return
      end
