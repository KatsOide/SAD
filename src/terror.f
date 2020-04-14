      subroutine terror(word,new,lfno,exist,errflg)
      use tfstk
      use ffs
      use ffs_pointer
      use tffitcode
      implicit none
      integer*4 lfno,ierr,i,id,kw,i0,istep,
     $     ifany,j,ie
      integer*8 k,jj,le,ldv,jd,ld
      real*8 errg,dl,dk,ddk,dtheta,dx,dy,err,er,v,
     $     getva,tgauss
      character*(*) word
      character*20  name,word1
      logical exist,new,uni,abbrev,temat,all,errflg,abb,
     1        cohe,add
      errflg=.false.
      new=.false.
1     uni=.false.
      cohe=.false.
      add=.false.
      errg=0.d0
      if(word .eq. 'DELX')then
        ierr=1
      elseif(word .eq. 'DELY')then
        ierr=2
      elseif(word .eq. 'DTHETA')then
        ierr=4
      elseif(word .eq. 'DL')then
        ierr=5
      elseif(word .eq. 'DK')then
        ierr=3
      elseif(word .eq. 'DDK')then
        ierr=6
      elseif(word .eq. 'DBEAM')then
        ierr=10
      elseif(word .eq. 'DUMP')then
        write(lfno,'(1x,a)')
     1  ' element     DL(mm)    DK          DDK      DTHETA(mrad)'//
     1  ' DELX(mm) DELY(mm)'
        call getwdl(word1)
        all=.false.
1012    exist=.false.
        if(all)then
          word=word1
          word1='*'
        endif
        do 1010 i=1,nlat-1
c     temat(ilist(2,latt(1))... MUST need because of latt defined by (2,0:nlatt)
          if(temat(i,name,word1))then
            exist=.true.
            j=idelc(i)
            k=latt(i)
            id=idtype(j)
            jd=idval(j)
            kw=kytbl(kwL,id)
            if(kw .ne. 0)then
              dl=rlist(k+kw)-rlist(jd+kw)
            else
              dl=0.d0
            endif
            dk=errk(1,i)-1.d0
            kw=kytbl(kwK0,id)
            if(kw .ne. 0)then
              ddk=rlist(k+kw)-rlist(jd+kw)
            else
              ddk=0.d0
            endif
            kw=kytbl(kwDROT,id)
            if(kw .ne. 0)then
              dtheta=rlist(k+kw)-rlist(jd+kw)
            else
              kw=kytbl(kwROT,id)
              if(kw .ne. 0)then
                dtheta=rlist(k+kw)-rlist(jd+kw)
              else
                dtheta=0.d0
              endif
            endif
            kw=kytbl(kwDX,id)
            if(kw .ne. 0)then
              dx=rlist(k+kw)-rlist(jd+kw)
            else
              dx=0.d0
            endif
            kw=kytbl(kwDY,id)
            if(kw .ne. 0)then
              dy=rlist(k+kw)-rlist(jd+kw)
            else
              dy=0.d0
            endif
            if(dl .ne. 0.d0 .or. dk .ne. 0.d0 .or.
     1         ddk .ne. 0.d0 .or. dtheta .ne. 0.d0 .or.
     1         dx .ne. 0.d0 .or. dy .ne. 0.d0)then
              write(lfno,9001)name(1:10),dl,dk,ddk,dtheta,dx,dy
9001          format(1x,a,3p,f10.5,1p,2g12.5,3p,3f10.5)
            endif
          endif
1010    continue
        if(exist)then
          new=all
          exist=.not. all
          return
        endif
        all=.true.
        go to 1012
      else
        exist=.false.
        return
      endif
      new=.true.
      call tfsetparam
3     err=getva(exist)
      if(.not. exist)then
        call getwdl(word1)
        if(abbrev(word1,'U_NIFORM','_'))then
          uni=.true.
          go to 3
        elseif(abbrev(word1,'R_ANDOM','_'))then
          uni=.false.
          go to 3
        elseif(abbrev(word1,'C_OHERENT','_'))then
          cohe=.true.
          go to 3
        elseif(abbrev(word1,'INC_OHERENT','_'))then
          cohe=.false.
          go to 3
        elseif(abbrev(word1,'A_DD','_'))then
          add=.true.
          go to 3
        elseif(abbrev(word1,'P_UT','_'))then
          add=.false.
          go to 3
        else
          call termes(lfno,'?Missing amount of error ',word)
          errflg=.true.
          return
        endif
      endif
      if(ierr .eq. 10)then
        jj=latt(0)
        k=idval(ilist(2,jj))
        rlist(j+1)=rlist(k+1)+sqrt(1.d0+rlist(k+1)**2)*err*tgauss()
        rlist(j+4)=rlist(k+4)+sqrt(1.d0+rlist(k+4)**2)*err*tgauss()
        rlist(j+2)=rlist(k+2)*abs(1.d0+err*tgauss())
        rlist(j+5)=rlist(k+5)*abs(1.d0+err*tgauss())
        exist=.true.
        return
      endif
2     call getwdl(word)
      if(abbrev(word,'U_NIFORM','_'))then
        uni=.true.
        go to 2
      elseif(abbrev(word,'R_ANDOM','_'))then
        uni=.false.
        go to 2
      elseif(abbrev(word,'C_OHERENT','_'))then
        cohe=.true.
        errg=err*tgauss()
        go to 2
      elseif(abbrev(word,'INC_OHERENT','_'))then
        cohe=.false.
        go to 2
      elseif(abbrev(word,'A_DD','_'))then
        add=.true.
        go to 2
      elseif(abbrev(word,'P_UT','_'))then
        add=.false.
        go to 2
      endif
      if(cohe)then
        if(errg .eq. 0.d0)then
          errg=err*tgauss()
        endif
      endif
      if(word .eq. ' ')then
        return
      endif
      abb=ifany(word,'*%{}',1) .gt. 0
      exist=.false.
      do 10 i=1,nlat-1
        id=idtypec(i)
c     temat(ilist(2,latt(1))... MUST need because of latt defined by (2,0:nlatt)
        if(temat(i,name,word))then
          if(id .eq. icDRFT)then
            ie=i
          else
            ie=master(i)
            if(ie .le. 0)then
              go to 10
            endif
          endif
          if(uni)then
            errg=err
          elseif(.not. cohe)then
            errg=err*tgauss()
          endif
          exist=.true.
          i0=i
          istep=1
          do 110 j=i0,ie,istep
            if((id .ne. icSOL .and. master(j) .ne. 0)
     $           .or. id .eq. icDRFT .or.id .eq. icSOL .and.
     $           idtypec(j) .eq. icSOL)then
              le=latt(j)
              ld=idval(idelc(j))
              if(add)then
                ldv=le
              else
                ldv=ld
              endif
              ie=j
              k=iele1(icomp(j))
              if(ierr .eq. 5)then
                kw=kytbl(kwL,id)
                if(kw .ne. 0)then
                  rlist(le+1)=max(0.d0,rlist(ldv+kw)+errg)
                endif
              elseif(ierr .eq. 3)then
                if(nelvx(k)%ival .gt. 0 .and. id .ne. icSOL
     $               .and. id .ne. icMARK)then
                  if(j .eq. nelvx(k)%klp)then
                    rlist(latt(j)+nelvx(k)%ival)
     $                   =rlist(latt(j)+nelvx(k)%ival)/errk(1,j)
                    er=1.d0
                  else
                    er=errk(1,nelvx(k)%klp)
                  endif
                  if(add)then
                    v=errk(1,j)
                  else
                    v=1.d0
                  endif
                  errk(1,j)=v+errg
                  if(id .ne. 20)then
                    rlist(le+nelvx(k)%ival)
     1                   =rlist(latt(nelvx(k)%klp)+nelvx(k)%ival)
     $                   /er*errk(1,j)*couple(j)
                  endif
                endif
              elseif(ierr .eq. 6)then
                kw=kytbl(kwK0,id)
                if(kw .ne. 0)then
                  if(add)then
                    rlist(le+kw)=rlist(le+kw)+errg
                  else
                    rlist(le+kw)=rlist(ldv+kw)+errg
                  endif                    
                endif
              elseif(ierr .eq. 1)then
                kw=kytbl(kwDX,id)
                if(kw .ne. 0)then
                  if(add)then
                    rlist(le+kw)=rlist(le+kw)+errg
                  else
                    rlist(le+kw)=rlist(ldv+kw)+errg
                  endif                    
                endif
              elseif(ierr .eq. 2)then
                kw=kytbl(kwDY,id)
                if(kw .ne. 0)then
                  if(add)then
                    rlist(le+kw)=rlist(le+kw)+errg
                  else
                    rlist(le+kw)=rlist(ldv+kw)+errg
                  endif                    
                endif
              elseif(ierr .eq. 4)then
                kw=kytbl(kwDROT,id)
                if(kw .ne. 0)then
                  if(add)then
                    rlist(le+kw)=rlist(le+kw)+errg
                  else
                    rlist(le+kw)=rlist(ldv+kw)+errg
                  endif                    
                else
                  kw=kytbl(kwROT,id)
                  if(kw .ne. 0)then
                    if(add)then
                      rlist(le+kw)=rlist(le+kw)+errg
                    else
                      rlist(le+kw)=rlist(ldv+kw)+errg
                    endif
                  endif                 
                endif
              endif
            endif
110       continue
          if(.not. abb)then
            go to 11
          endif
        endif
10    continue
11    if(.not. uni)then
c     tfadjust(ilist(2,latt(1))... MUST need because of latt defined by (2,0:nlatt)
        call tfadjst
      endif
      if(exist)then
        go to 2
      else
        go to 1
      endif
      end
