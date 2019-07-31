      subroutine tffsfreefix(frefix,nvar,lfno)
      use tfstk
      use ffs, only:nve,nele,nlat
      use ffs_pointer
      use ffs_fit
      use tffitcode
      use mackw
      use tfcsi,only:cssetp
      use tflinepcom, only:tftouch
      implicit none
      integer*8 kal
      integer*4 nvar,lfno,l,ifany,
     $     i,j,k,it,iv,ivi,next,itk,ivk,lenw,kk,jj,ivck,
     $     irtc,nl,kkk
      logical*4 frefix,tmatch,wild,found,comp,temat
      character*256 word1,keyword,tfkwrd,tfkwrd1,word,nlist1
      character*(MAXPNAME+16) name,name1
 1    call peekwdp(word,next)
      if(word .eq. ' ')then
        return
      endif
      wild=ifany(word,'%*<{|',1) .gt. 0
      comp=index(word,'.') .gt. 0
      found=.false.
      it=0
      iv=0
      if(frefix)then
        LOOP_K_1: do k=1,ntwissfun
          l=lenw(nlist(k))
          nlist1=nlist(k)(:l)//'I'
          if(word .eq. nlist1(:l+1))then
            call cssetp(next)
            found=.true.
            it=41
            word=pnamec(1)
            iv=k
            exit LOOP_K_1
          endif
        enddo LOOP_K_1
        call tfgetlineps(word,len_trim(word),nl,kal,0,irtc)
        if(irtc .ne. 0 .or. nl .le. 0)then
          go to 3000
        endif
        LOOP_K_2: do kkk=1,nl
          k=int(rlist(kal+kkk))
          if(comp)then
            if(.not. temat(k,name,word))then
              cycle LOOP_K_2
            endif
            ivck=k
          else
c     Note: Skip no-head multiple elements
c     *     klp(iele1(k)) == k if singlet or head of multipole elements
            if(klp(iele1(k)) .ne. k .or.
     $           .not. tmatch(pnamec(k),word))then
              cycle LOOP_K_2
            endif
            ivck=0
          endif
          if(.not. found)then
            call cssetp(next)
            call peekwd(word1,next)
          endif
          itk=idtypec(k)
          if(itk .ne. it)then
            if(word1 .ne. ' ')then
              ivk=1
 1010         continue
              keyword=tfkwrd(itk,ivk)
              if(keyword .eq. ' ')then
                if(it .ne. 0)then
                  cycle LOOP_K_2
                endif
                word1=' '
                iv=0
                go to 1020
              else
                if(keyword .eq. word1)then
                  call cssetp(next)
                  iv=ivk
                  go to 1020
                endif
                keyword=tfkwrd1(itk,ivk)
                if(keyword .eq. word1)then
                  call cssetp(next)
                  iv=ivk
                  go to 1020
                endif
              endif
              ivk=ivk+1
              go to 1010
 1020         continue
            endif
            it=itk
          endif
          found=.true.
          kk=iele1(k)
          if(iv .eq. 0)then
            if(ival(kk) .eq. 0)then
              call termes(lfno,'Can''t use as variable ',
     $             pnamec(k))
              return
            endif
          endif
          LOOP_I_1: do i=1,nvar
            if(ivarele(i) .eq. kk)then
              if(iv .eq. 0)then
                ivi=ival(kk)
              else
                ivi=iv
              endif
              if(ivvar(i) .eq. ivi)then
                if(comp)then
                  if(ivcomp(i) .eq. 0)then
                    call termes(lfno,
     $                   'Element already used as variable: ',
     $                   pnamec(k))
                    return
                  elseif(ivcomp(i) .ne. k)then
                    cycle LOOP_I_1
                  endif
                else
                  if(ivcomp(i) .ne. 0)then
                    call elname(ivcomp(i),name)
                    call termes(lfno,
     $        'A component has been already used as variable: ',
     $                   name)
                    return
                  endif
                endif
                go to 10
              endif
            endif
          enddo LOOP_I_1
          if(nvar .ge. nve)then
            call termes(lfno,'Too many variables',' ')
            return
          endif
          LOOP_I_2: do i=1,nvar
            if(ivarele(i) .ge. kk)then
              do j=nvar,i,-1
                ivarele(j+1)=ivarele(j)
                ivvar(j+1)=ivvar(j)
                ivcomp(j+1)=ivcomp(j)
                valvar2(j+1,1)=valvar2(j,1)
                valvar2(j+1,2)=valvar2(j,2)
              enddo
              go to 11
            endif
          enddo LOOP_I_2
          i=nvar+1
 11       nvar=nvar+1
          ivarele(i)=kk
          if(iv .eq. 0)then
            ivi=ival(kk)
          else
            ivi=iv
          endif
          ivvar(i)=ivi
          valvar2(i,1)=tfvalvar(k,ivi)
          ivcomp(i)=ivck
          if(ivi .eq. ival(kk))then
            if(comp)then
              call elnameK(iele(k),name1)
              if(iele(k) .eq. k)then
                do jj=1,nlat-1
                  if(jj .ne. k .and. iele(jj) .eq. k)then
                    call elnameK(jj,name)
                    call termes(lfno,'Info-Component '//
     $                   name(1:lenw(name))//' is coupled to ',
     $                   name1(1:lenw(name1))//' .')
                  endif
                enddo
              else
                call elnameK(k,name)
                call termes(lfno,'Info-Component '//
     $               name(1:lenw(name))//' is coupled to ',
     $               name1(1:lenw(name1))//' .')
              endif
            else
              do jj=1,nlat-1
                if(klp(iele1(jj)) .eq. k
     $               .and. iele(jj) .ne. k)then
c     `jj' is same family but different master with `k',
c     where klp(iele1(k)) == k
                  call elnameK(jj,name)
                  call elnameK(iele(jj),name1)
                  call termes(lfno,'Info-Component '//
     $                 name(1:lenw(name))//' is coupled to ',
     $                 name1(1:lenw(name1)))
                endif
              enddo
            endif
            valvar2(i,1)=valvar2(i,1)/errk(1,k)
            valvar2(i,2)=valvar2(i,1)
          else
            valvar2(i,2)=valvar2(i,1)
            if(.not. comp)then
              call tftouch(kk,ivi)
            endif
          endif
c          write(*,*)'tffsfreefix ',i,k,ivi,valvar2(i,1),valvar2(i,2)
 10       if(.not. wild)then
            go to 1
          endif
        enddo LOOP_K_2
      else
        iv=0
        it=0
        LOOP_K_3: do k=1,ntwissfun
          l=lenw(nlist(k))
          nlist1=nlist(k)(:l)//'I'
          if(word .eq. nlist1(:l+1))then
            call cssetp(next)
            found=.true.
            it=41
            word=pnamec(1)
            iv=k
            exit LOOP_K_3
          endif
        enddo LOOP_K_3
        LOOP_I_3: do i=1,nvar
 1210     continue
          kk=ivarele(i)
          if(ivcomp(i) .ne. 0)then
            if(.not. temat(ivcomp(i),name,word))then
              cycle LOOP_I_3
            endif
          elseif(comp)then
            cycle LOOP_I_3
          elseif(.not. tmatch(pnamec(klp(kk)),word))then
            cycle LOOP_I_3
          endif
          if(.not. found)then
            call cssetp(next)
            call peekwd(word1,next)
          endif
          found=.true.
          itk=idtypec(klp(kk))
          if(itk .ne. it)then
            if(word1 .ne. ' ')then
              ivk=1
 1110         continue
              keyword=tfkwrd(itk,ivk)
              if(keyword .eq. ' ')then
                if(it .ne. 0)then
                  cycle LOOP_I_3
                endif
                word1=' '
                iv=0
                go to 1120
              else
                if(keyword .eq. word1)then
                  call cssetp(next)
                  iv=ivk
                  go to 1120
                endif
                keyword=tfkwrd1(itk,ivk)
                if(keyword .eq. word1)then
                  call cssetp(next)
                  iv=ivk
                  go to 1120
                endif
              endif                
              ivk=ivk+1
              go to 1110
 1120         continue
            endif
            it=itk
          endif
          if(iv .eq. 0 .or. iv .eq. ivvar(i))then
            do j=i,nvar-1
              ivarele(j)=ivarele(j+1)
              ivvar(j)=ivvar(j+1)
              ivcomp(j)=ivcomp(j+1)
              valvar2(j,1)=valvar2(j+1,1)
              valvar2(j,2)=valvar2(j+1,2)
            enddo
            nvar=nvar-1
            if(i .gt. nvar .or.
     $           .not. wild .and. iv .ne. 0)then
              go to 1
            endif
          else
            go to 1220
          endif
          go to 1210
 1220     continue
        enddo LOOP_I_3
        if(.not. found .and. .not. wild)then
          do k=1,nele
            if(tmatch(pnamec(k),word))then
              call cssetp(next)
              found=.true.
              go to 1
            endif
          enddo
        endif
      endif
 3000 if(wild .and. .not. found)then
        call cssetp(next)
      endif
      if(found .or. wild)then
        go to 1
      endif
      return
      end
