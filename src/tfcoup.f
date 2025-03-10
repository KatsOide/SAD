      subroutine tfcoup(exist)
      use tfstk
      use ffs
      use ffs_pointer
      use tffitcode
      use ffs_seg
      use ffsa, only:tffssetupcouple
      implicit none
      type (sad_comp), pointer :: cmp
      logical*4 ,intent(out):: exist
      integer*4 kk1,k1,i,kk2,k2,ielm,lenw
      real*8 co,v,getva
      character*(MAXPNAME+16) ele1,ele2,name
      logical*4 comp
      call getwdlp(ele1)
      if(ele1 == '*')then
        call tffssetupcouple()
        exist=.true.
        return
      endif
      call getwdlp(ele2)
      co=getva(exist)
      if( .not. exist)then
        call termes('Missing coefficient for COUP_LE.',' ')
        go to 9000
      endif
      comp=index(ele1,'.') .gt. 0
      kk1=ielm(ele1,exist)
      k1=iele1(kk1)
      if(.not. exist)then
        call termes('Undefined slave element for COUP_LE: ',
     $       ele1)
        go to 9000
      endif
      kk2=ielm(ele2,exist)
      if(.not. exist)then
        call termes('Undefined master component for COUP_LE: ',
     $       ele2)
        go to 9000
      endif
      k2=iele1(kk2)
      if(idtypec(kk1) /=
     $     idtypec(kk2))then
        call termes('Different types of element to COUP_LE.',' ')
        go to 9000
      endif
      if(nelvx(k1)%ival /= nelvx(k2)%ival)then
        call termes('Different keywords of element to COUP_LE.',' ')
        go to 9000
      endif
      evarini=.true.
      if(icomp(kk2) /= kk2)then
        if(kk2 == nelvx(k2)%klp)then
          call elnameK(kk2,name)
          call termes('Info-COUPLEs of Components '//
     $         pname(idelc(kk2))(1:lpnamec(kk2))//'.*'//
     $         ' have been reset to ',name(1:lenw(name))//' .')
          do i=1,nlat-1
            if(iele1(i) == k2)then
              icomp(i)=kk2
              couple(i)=1.d0
            endif
          enddo
        else
          call termes('Info-Component '//ele2(1:lenw(ele2))//
     $         ' has been made uncoupled.',' ')
          icomp(kk2)=kk2
          couple(kk2)=1.d0
        endif
      endif
      v=tfvalvar(kk2,nelvx(k2)%ival)/errk(1,kk2)
c      v=rlist(latt(kk2)+ival(k2))/errk(1,kk2)
      if(comp)then
        if(kk1 == kk2 .and. nelvx(iele1(kk1))%klp == kk1)then
          co=1.d0
          do i=1,nlat-1
            if(icomp(i) == kk1 .and. i /= kk1)then
              call elname(i,name)
              call termes('Info-Component '//name(1:lenw(name))//
     $             ' has been made uncoupled.',' ')
              icomp(i)=i
              couple(i)=1.d0
            endif
          enddo
        else
          do i=1,nlat-1
            if(icomp(i) == kk1 .and. i /= kk1)then
              call tfdecoupcomp(i)
            endif
          enddo
        endif
        icomp(kk1)=kk2
        couple(kk1)=co
        call compelc(kk1,cmp)
c        call tfsetcmp(v*errk(1,kk1)*co,cmp,ival(k1))
        cmp%value(nelvx(k1)%ival)=v*errk(1,kk1)*co
      else
        if(k1 /= k2)then
          do i=1,nlat-1
            if(iele1(i) /= k1)then
              if(iele1(icomp(i)) == k1)then
                call tfdecoupcomp(i)
              endif
            else
              icomp(i)=kk2
              couple(i)=co
              rlist(latt(i)+nelvx(k1)%ival)=v*errk(1,i)*co
            endif
          enddo
        else
          do i=1,nlat-1
            if(iele1(i) == k1)then
              icomp(i)=kk2
              couple(i)=co
              rlist(latt(i)+nelvx(k1)%ival)=v*errk(1,i)*co
            endif
          enddo
        endif
      endif
      if( .not. exist)then
        call termes('Undefined element for COUP_LE: ',ele1)
      endif
      return
 9000 call termes('Usage: COUP_LE slave master coefficient',' ')
      return
      end

      subroutine tfdecoup()
      use tfstk
      use ffs
      use ffs_pointer
      use tffitcode
      use tfcsi,only:ipoint
      implicit none
      integer*4 next,i,j
      character*(MAXPNAME+16) ele,name
      logical*4 mat,temat
 1    call peekwdp(ele,next)
      if(ele == ' ')then
        return
      endif
      mat=.false.
      do i=1,nlat-1
        if(.not. mat)then
          if(.not. temat(i,name,ele))cycle
          mat=.true.
          ipoint=next
        endif
        j=icomp(i)
        if(j /= nelvx(iele1(i))%klp)then
          if(temat(i,name,ele))then
            icomp(i)=nelvx(iele1(i))%klp
            couple(i)=1.d0
          elseif(nelvx(iele1(j))%klp /= j)then
            if(temat(j,name,ele))then
              call tfdecoupcomp(i)
            endif
          endif
        endif
      enddo
      if(mat)then
        evarini=.true.
        go to 1
      endif
      return
      end

      subroutine tfindep
      use tfstk
      use ffs
      use ffs_pointer
      use tffitcode
      use tfcsi,only:ipoint
      implicit none
      integer*4 next,i
      character*(MAXPNAME+16) ele,name
      logical*4 mat,temat
 1    call peekwdp(ele,next)
      if(ele == ' ')then
        return
      endif
      mat=.false.
      do i=1,nlat-1
        if(temat(i,name,ele))then
          if(.not. mat)then
            mat=.true.
            ipoint=next
          endif
          icomp(i)=i
        endif
      enddo
      if(mat)then
        evarini=.true.
        go to 1
      endif
      return
      end

      subroutine tfdecoupcomp(i)
      use tfstk
      use ffs
      use ffs_pointer
      use tffitcode
      implicit none
      integer*4 ,intent(in):: i
      integer*4 lenw
      character*(MAXPNAME+16)name,name1
      call elname(i,name)
      call elnameK(nelvx(iele1(i))%klp,name1)
      call termes('Info-COUPLE of component '//name(1:lenw(name))//
     $     ' has been reset to ',name1(1:lenw(name1))//' .')
      icomp(i)=nelvx(iele1(i))%klp
      couple(i)=1.d0
      return
      end
