      subroutine tffsfreefix(latt,frefix,
     $     ivarele,ivvar,ivcomp,valvar,nvar,nve,
     $     itouchele,itouchv,ntouch,
     $     ival,klp,nele,errk,mult,iele,iele1,
     $     nlat,nlist,lfno)
      use tfstk
      use tffitcode
      implicit none
      include 'inc/MACCODE.inc'
      integer*8 kal
      integer*4 nvar,nele,lfno,nlat,nve,
     $     ivarele(nve),ivvar(nve),l,
     $     itouchele(nve),itouchv(nve),ntouch,
     $     ivcomp(nve),mult(nlat),iele(nlat),iele1(nlat),
     $     klp(nele),ival(nele),ifany,latt(2,nlat),
     $     i,j,k,it,iv,ivi,next,itk,ivk,lenw,kk,jj,ivck,
     $     irtc,nl,kkk
      real*8 valvar(nve,2),errk(2,nlat)
      logical*4 frefix,tmatch,wild,found,comp,temat
      character*(*) nlist(ntwissfun)
      character*80 word1,keyword,tfkwrd
      character*80 word,nlist1
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
            word=pname(latt(1,1))
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
            if(.not. temat(latt,k,mult,name,word))then
              cycle LOOP_K_2
            endif
            ivck=k
          else
c     Note: Skip no-head multiple elements
c     *     klp(iele1(k)) == k if singlet or head of multipole elements
            if(klp(iele1(k)) .ne. k .or.
     $           .not. tmatch(pname(latt(1,k)),word))then
              cycle LOOP_K_2
            endif
            ivck=0
          endif
          if(.not. found)then
            call cssetp(next)
            call peekwd(word1,next)
          endif
          itk=idtype(latt(1,k))
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
              elseif(keyword .eq. word1)then
                call cssetp(next)
                iv=ivk
                go to 1020
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
     $             pname(latt(1,k)))
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
     $                   pname(latt(1,k)))
                    return
                  elseif(ivcomp(i) .ne. k)then
                    cycle LOOP_I_1
                  endif
                else
                  if(ivcomp(i) .ne. 0)then
                    call elname(latt,ivcomp(i),mult,name)
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
                valvar(j+1,1)=valvar(j,1)
                valvar(j+1,2)=valvar(j,2)
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
          valvar(i,1)=rlist(latt(2,k)+ivi)
          ivcomp(i)=ivck
          if(ivi .eq. ival(kk))then
            if(comp)then
              call elname1(latt,iele(k),mult,name1,.true.)
              if(iele(k) .eq. k)then
                do jj=1,nlat-1
                  if(jj .ne. k .and. iele(jj) .eq. k)then
                    call elname1(latt,jj,mult,name,.true.)
                    call termes(lfno,'Info-Component '//
     $                   name(1:lenw(name))//' is coupled to ',
     $                   name1(1:lenw(name1))//' .')
                  endif
                enddo
              else
                call elname1(latt,k,mult,name,.true.)
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
                  call elname1(latt,jj,mult,name,.true.)
                  call elname1(latt,iele(jj),mult,name1,.true.)
                  call termes(lfno,'Info-Component '//
     $                 name(1:lenw(name))//' is coupled to ',
     $                 name1(1:lenw(name1)))
                endif
              enddo
            endif
            valvar(i,1)=valvar(i,1)/errk(1,k)
            valvar(i,2)=valvar(i,1)
          else
            valvar(i,2)=valvar(i,1)
            if(.not. comp)then
              do j=1,ntouch
                if(itouchele(j) .eq. kk .and.
     $               itouchv(j) .eq. ivi)then
                  go to 10
                endif
              enddo
              ntouch=ntouch+1
              itouchele(ntouch)=kk
              itouchv(ntouch)=ivi
            endif
          endif
c          write(*,*)'tffsfreefix ',i,k,ivi,valvar(i,1),valvar(i,2)
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
            word=pname(latt(1,1))
            iv=k
            exit LOOP_K_3
          endif
        enddo LOOP_K_3
        LOOP_I_3: do i=1,nvar
 1210     continue
          kk=ivarele(i)
          if(ivcomp(i) .ne. 0)then
            if(.not. temat(latt,ivcomp(i),mult,name,word))then
              cycle LOOP_I_3
            endif
          elseif(comp)then
            cycle LOOP_I_3
          elseif(.not. tmatch(pname(latt(1,klp(kk))),word))then
            cycle LOOP_I_3
          endif
          if(.not. found)then
            call cssetp(next)
            call peekwd(word1,next)
          endif
          found=.true.
          itk=idtype(latt(1,klp(kk)))
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
              elseif(keyword .eq. word1)then
                call cssetp(next)
                iv=ivk
                go to 1120
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
              valvar(j,1)=valvar(j+1,1)
              valvar(j,2)=valvar(j+1,2)
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
            if(tmatch(pname(latt(1,klp(k))),word))then
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
