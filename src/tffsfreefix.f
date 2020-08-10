      subroutine tffsfreefix(frefix,nvar,lfno)
      use tfstk
      use ffs, only:nele,nlat,nvevx,nelvx,ffv
      use ffs_pointer
      use ffs_fit
      use tffitcode
      use mackw
      use tfcsi,only:ipoint
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
            ipoint=next
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
            if(nelvx(iele1(k))%klp .ne. k .or.
     $           .not. tmatch(pnamec(k),word))then
              cycle LOOP_K_2
            endif
            ivck=0
          endif
          if(.not. found)then
            ipoint=next
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
                  ipoint=next
                  iv=ivk
                  go to 1020
                endif
                keyword=tfkwrd1(itk,ivk)
                if(keyword .eq. word1)then
                  ipoint=next
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
            if(nelvx(kk)%ival .eq. 0)then
              call termes(lfno,'Can''t use as variable ',
     $             pnamec(k))
              return
            endif
          endif
          LOOP_I_1: do i=1,nvar
            if(nvevx(i)%ivarele .eq. kk)then
              if(iv .eq. 0)then
                ivi=nelvx(kk)%ival
              else
                ivi=iv
              endif
              if(nvevx(i)%ivvar .eq. ivi)then
                if(comp)then
                  if(nvevx(i)%ivcomp .eq. 0)then
                    call termes(lfno,
     $                   'Element already used as variable: ',
     $                   pnamec(k))
                    return
                  elseif(nvevx(i)%ivcomp .ne. k)then
                    cycle LOOP_I_1
                  endif
                else
                  if(nvevx(i)%ivcomp .ne. 0)then
                    call elname(nvevx(i)%ivcomp,name)
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
          call tffsnvealloc(nvar+1)
          LOOP_I_2: do i=1,nvar
            if(nvevx(i)%ivarele .ge. kk)then
              do j=nvar,i,-1
                nvevx(j+1)%ivarele=nvevx(j)%ivarele
                nvevx(j+1)%ivvar=nvevx(j)%ivvar
                nvevx(j+1)%ivcomp=nvevx(j)%ivcomp
                nvevx(j+1)%valvar=nvevx(j)%valvar
                nvevx(j+1)%valvar2=nvevx(j)%valvar2
              enddo
              go to 11
            endif
          enddo LOOP_I_2
          i=nvar+1
 11       nvar=nvar+1
          ffv%evarini=.true.
          nvevx(i)%ivarele=kk
          if(iv .eq. 0)then
            ivi=nelvx(kk)%ival
          else
            ivi=iv
          endif
          nvevx(i)%ivvar=ivi
          nvevx(i)%valvar=tfvalvar(k,ivi)
          nvevx(i)%ivcomp=ivck
          if(ivi .eq. nelvx(kk)%ival)then
            if(comp)then
              call elnameK(icomp(k),name1)
              if(icomp(k) .eq. k)then
                do jj=1,nlat-1
                  if(jj .ne. k .and. icomp(jj) .eq. k)then
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
                if(nelvx(iele1(jj))%klp .eq. k
     $               .and. icomp(jj) .ne. k)then
c     `jj' is same family but different master with `k',
c     where klp(iele1(k)) == k
                  call elnameK(jj,name)
                  call elnameK(icomp(jj),name1)
                  call termes(lfno,'Info-Component '//
     $                 name(1:lenw(name))//' is coupled to ',
     $                 name1(1:lenw(name1)))
                endif
              enddo
            endif
            nvevx(i)%valvar=nvevx(i)%valvar/errk(1,k)
            nvevx(i)%valvar2=nvevx(i)%valvar
          else
            nvevx(i)%valvar2=nvevx(i)%valvar
            if(.not. comp)then
              call tftouch(kk,ivi)
            endif
          endif
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
            ipoint=next
            found=.true.
            it=41
            word=pnamec(1)
            iv=k
            exit LOOP_K_3
          endif
        enddo LOOP_K_3
        LOOP_I_3: do i=1,nvar
 1210     continue
          kk=nvevx(i)%ivarele
          if(nvevx(i)%ivcomp .ne. 0)then
            if(.not. temat(nvevx(i)%ivcomp,name,word))then
              cycle LOOP_I_3
            endif
          elseif(comp)then
            cycle LOOP_I_3
          elseif(.not. tmatch(pnamec(nelvx(kk)%klp),word))then
            cycle LOOP_I_3
          endif
          if(.not. found)then
            ipoint=next
            call peekwd(word1,next)
          endif
          found=.true.
          itk=idtypec(nelvx(kk)%klp)
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
                  ipoint=next
                  iv=ivk
                  go to 1120
                endif
                keyword=tfkwrd1(itk,ivk)
                if(keyword .eq. word1)then
                  ipoint=next
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
          if(iv .eq. 0 .or. iv .eq. nvevx(i)%ivvar)then
            do j=i,nvar-1
              nvevx(j)%ivarele=nvevx(j+1)%ivarele
              nvevx(j)%ivvar=nvevx(j+1)%ivvar
              nvevx(j)%ivcomp=nvevx(j+1)%ivcomp
              nvevx(j)%valvar=nvevx(j+1)%valvar
              nvevx(j)%valvar2=nvevx(j+1)%valvar2
            enddo
            nvar=nvar-1
            ffv%evarini=.true.
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
              ipoint=next
              found=.true.
              go to 1
            endif
          enddo
        endif
      endif
 3000 if(wild .and. .not. found)then
        ipoint=next
      endif
      if(found .or. wild)then
        go to 1
      endif
      return
      end
