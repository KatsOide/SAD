      subroutine tfvectorize(isp1,kx,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) kx,kx1,kf,ki,k1
      type (sad_list), pointer :: kli,kli1,klx1,kl,klx
      integer*8 kaf,kaf1,kai,isp1,isp0,isp2,i,idsp
      integer*4 irtc,j,nv,itfmessage
      kf=dtastk(isp1+1)
      if(ktfoperqd(kf,kaf))then
        if(kaf .eq. mtfplus .or. kaf .eq. mtftimes
     $       .or. kaf .eq. mtfpower .or. kaf .eq. mtfrevpower)then
          nv=0
          do i=isp1+2,isp
            if(ktflistqd(dtastk(i),kli) .and. (kli%head .eq. kxvect
     $           .or. kli%head .eq. kxvect1))then
              ktastk(i)=kli%body(1)
              if(tflistqd(dtastk(i),kli1))then
                if(nv .eq. 0)then
                  nv=kli1%nl
                elseif(nv .ne. kli1%nl)then
                  irtc=itfmessage(999,'General::equalleng',
     $                 '"two Vectors"')
                  return
                endif
              endif
            endif
          enddo
          call tfefunref(isp1+1,kx1,.true.,irtc)
          if(irtc .ne. 0)then
            return
          endif
c          call tfdebugprint(kx1,'vectorize',1)
          if(tflistqd(kx1,klx1) .and. klx1%nl .eq. nv)then
            go to 9000
          endif
          kx=kx1
          return
        elseif(kaf .eq. mtfset .or. kaf .eq. mtfsetdelayed)then
          do i=isp1+2,isp
            if(ktflistq(ktastk(i),kli) .and. kli%head .eq. kxvect1)then
              k1=kli%dbody(1)
              dtastk(i)=kxadaloc(-1,1,kli1)
              kli1%head=ktfcopy1(kxvect)
              kli1%dbody(1)=dtfcopy(k1)
            endif
          enddo
          call tfefunref(isp1+1,kx,.false.,irtc)
          return
        elseif(kaf .gt. mtfend)then
          kaf1=klist(ifunbase+kaf)+1
          if(ilist(1,kaf1) .eq. 1 .and. ilist(1,kaf1+1) .ne. 0 .and.
     $         isp .eq. isp1+2 .and. ktflistq(ktastk(isp),kl) .and.
     $         (kl%head .eq. kxvect .or. kl%head .eq. kxvect1))then
            ktastk(isp)=kl%body(1)
            if(tflistqk(ktastk(isp)))then
              call tfefunref(isp1+1,kx1,.false.,irtc)
              go to 9000
            else
              irtc=itfmessage(999,'General::wrongtype',
     $             "Vector[{...}]")
              return
            endif
          endif
        endif
      endif
      nv=0
      do i=isp1+2,isp
        ki=dtastk(i)
        if(ktflistqd(ki,kli) .and.
     $       (kli%head .eq. kxvect .or. kli%head .eq. kxvect1))then
          k1=kli%dbody(1)
          if(tflistqd(k1,kli1))then
            if(nv .gt. 0)then
              if(nv .ne. kli1%nl)then
                irtc=itfmessage(999,'General::equalleng',
     $               '"two Vectors"')
                return
              endif
            else
              nv=kli1%nl
            endif
            ktastk(i)=ktfref+ktfaddr(k1)
          else
            irtc=itfmessage(999,'General::wrongtype',
     $           "Vector[{...}]")
            return
          endif
        endif
      enddo
      if(nv .eq. 0)then
        call tfefunref(isp1+1,kx,.false.,irtc)
        return
      endif
      isp2=isp
      isp0=isp+nv+1
      idsp=nv+isp2-isp1
      do j=1,nv
        dtastk(isp0)=kf
        do i=isp1+2,isp2
          ki=dtastk(i)
          if(ktfrefqd(ki,kai))then
            ktastk(idsp+i)=klist(kai+j)
          else
            ktastk(idsp+i)=ktastk(i)
          endif
        enddo
        isp=isp2+idsp
        call tfefunref(isp0,ktastk(isp2+j),.true.,irtc)
        if(irtc .ne. 0)then
          return
        endif
      enddo
      isp=isp0-1
      kx1=kxmakelist(isp2)
      isp=isp2
      irtc=0
 9000 if(irtc .ne. 0)then
        return
      endif
      kx=kxadaloc(-1,1,klx)
      klx%head=ktfcopy1(kxvect)
      klx%dbody(1)=dtfcopy(kx1)
      irtc=0
      return
      end
