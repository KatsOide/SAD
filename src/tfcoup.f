      subroutine tfcoup(latt,
     $     couple,iele,iele1,klp,ival,errk,mult,
     $     lfno,exist)
      use tfstk
      use ffs
      use tffitcode
      implicit none
      integer*4 lfno,kk1,k1,i,kk2,k2,ielm,lenw,lpname
      integer*4 latt(2,nlat),iele(nlat),iele1(nlat),
     $     ival(nele),klp(nele),mult(nlat)
      real*8 errk(2,nlat),couple(nlat),co,v,getva
      character*(MAXPNAME+16) ele1,ele2,name
      logical*4 exist,comp
      call getwdlp(ele1)
      call getwdlp(ele2)
      co=getva(exist)
      if( .not. exist)then
        call termes(lfno,'Missing coefficient for COUP_LE.',' ')
        go to 9000
      endif
      comp=index(ele1,'.') .gt. 0
      kk1=ielm(latt,ele1,1,mult,exist)
      k1=iele1(kk1)
      if(.not. exist)then
        call termes(lfno,'Undefined slave element for COUP_LE: ',
     $       ele1)
        go to 9000
      endif
      kk2=ielm(latt,ele2,1,mult,exist)
      if(.not. exist)then
        call termes(lfno,'Undefined master component for COUP_LE: ',
     $       ele2)
        go to 9000
      endif
      k2=iele1(kk2)
      if(idtype(latt(1,kk1)) .ne. idtype(latt(1,kk2)))then
        call termes(lfno,
     $       'Different types of element to COUP_LE.',' ')
        go to 9000
      endif
      if(ival(k1) .ne. ival(k2))then
        call termes(lfno,
     $       'Different keywords of element to COUP_LE.',' ')
        go to 9000
      endif
      if(iele(kk2) .ne. kk2)then
        if(kk2 .eq. klp(k2))then
          call elname1(latt,kk2,mult,name,.true.)
          call termes(lfno,'Info-COUPLEs of Components '//
     $         pname(latt(1,kk2))(1:lpname(latt(1,kk2)))//'.*'//
     $         ' have been reset to ',name(1:lenw(name))//' .')
          do i=1,nlat-1
            if(iele1(i) .eq. k2)then
              iele(i)=kk2
              couple(i)=1.d0
            endif
          enddo
        else
          call termes(lfno,'Info-Component '//ele2(1:lenw(ele2))//
     $         ' has been made uncoupled.',' ')
          iele(kk2)=kk2
          couple(kk2)=1.d0
        endif
      endif
      v=rlist(latt(2,kk2)+ival(k2))/errk(1,kk2)
      if(comp)then
        if(kk1 .eq. kk2 .and. klp(iele1(kk1)) .eq. kk1)then
          co=1.d0
          do i=1,nlat-1
            if(iele(i) .eq. kk1 .and. i .ne. kk1)then
              call elname(latt,i,mult,name)
              call termes(lfno,'Info-Component '//name(1:lenw(name))//
     $             ' has been made uncoupled.',' ')
              iele(i)=i
              couple(i)=1.d0
            endif
          enddo
        else
          do i=1,nlat-1
            if(iele(i) .eq. kk1 .and. i .ne. kk1)then
              call tfdecoupcomp(latt,i,lfno,iele,iele1,
     $             mult,klp,couple)
            endif
          enddo
        endif
        iele(kk1)=kk2
        couple(kk1)=co
        rlist(latt(2,kk1)+ival(k1))=v*errk(1,kk1)*co
      else
        if(k1 .ne. k2)then
          do i=1,nlat-1
            if(iele1(i) .ne. k1)then
              if(iele1(iele(i)) .eq. k1)then
                call tfdecoupcomp(
     $               latt,i,lfno,iele,iele1,mult,klp,couple)
              endif
            else
              iele(i)=kk2
              couple(i)=co
              rlist(latt(2,i)+ival(k1))=v*errk(1,i)*co
            endif
          enddo
        else
          do i=1,nlat-1
            if(iele1(i) .eq. k1)then
              iele(i)=kk2
              couple(i)=co
              rlist(latt(2,i)+ival(k1))=v*errk(1,i)*co
            endif
          enddo
        endif
      endif
      if( .not. exist)then
        call termes(lfno,'Undefined element for COUP_LE: ',ele1)
      endif
      return
 9000 call termes(lfno,'Usage: COUP_LE slave master coefficient',' ')
      return
      end

      subroutine tfdecoup(latt,
     $     couple,iele,iele1,klp,mult,
     $     lfno)
      use tfstk
      use ffs
      use tffitcode
      implicit none
      integer*4 latt(2,nlat),iele(nlat),iele1(nlat),
     $     klp(nele),mult(nlat),lfno,next,i,j
      real*8 couple(nlat)
      character*(MAXPNAME+16) ele,name
      logical*4 mat,temat
 1    call peekwdp(ele,next)
      if(ele .eq. ' ')then
        return
      endif
      mat=.false.
      do i=1,nlat-1
        if(.not. mat)then
          if(.not. temat(latt,i,mult,name,ele))cycle
          mat=.true.
          call cssetp(next)
        endif
        j=iele(i)
        if(j .ne. klp(iele1(i)))then
          if(temat(latt,i,mult,name,ele))then
            iele(i)=klp(iele1(i))
            couple(i)=1.d0
          elseif(klp(iele1(j)) .ne. j)then
            if(temat(latt,j,mult,name,ele))then
              call tfdecoupcomp(
     $             latt,i,lfno,iele,iele1,mult,klp,couple)
            endif
          endif
        endif
      enddo
      if(mat)then
        go to 1
      endif
      return
      end

      subroutine tfindep(latt,iele,mult)
      use tfstk
      use ffs
      use tffitcode
      implicit none
      integer*4 latt(2,nlat),iele(nlat),mult(nlat),next,i
      character*(MAXPNAME+16) ele,name
      logical*4 mat,temat
 1    call peekwdp(ele,next)
      if(ele .eq. ' ')then
        return
      endif
      mat=.false.
      do i=1,nlat-1
        if(temat(latt,i,mult,name,ele))then
          if(.not. mat)then
            mat=.true.
            call cssetp(next)
          endif
          iele(i)=i
        endif
      enddo
      if(mat)then
        go to 1
      endif
      return
      end

      subroutine tfdecoupcomp(latt,i,lfno,iele,iele1,mult,klp,couple)
      use tfstk
      use ffs
      use tffitcode
      implicit none
      integer*4 latt(2,nlat),iele(nlat),iele1(nlat),klp(nlat),
     $     mult(nlat),i,lfno,lenw
      real*8 couple(nlat)
      character*(MAXPNAME+16)name,name1
      call elname(latt,i,mult,name)
      call elname1(latt,klp(iele1(i)),mult,name1,.true.)
      call termes(lfno,'Info-COUPLE of component '//name(1:lenw(name))//
     $     ' has been reset to ',name1(1:lenw(name1))//' .')
      iele(i)=klp(iele1(i))
      couple(i)=1.d0
      return
      end
