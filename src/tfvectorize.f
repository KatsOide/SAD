      subroutine tfvectorize(isp1,kx,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) ,intent(out):: kx
      type (sad_descriptor) kx1,kf,ki,k1,tfefunrefu,tfefunrefd
      type (sad_dlist), pointer :: kli,kli1,klx1,kl,klx
      integer*8 kaf,kaf1,kai
      integer*4 ,intent(in):: isp1
      integer*4 ,intent(out):: irtc
      integer*4 i,j,isp0,isp2,nv,idsp,itfmessage
      kf=dtastk(isp1+1)
      if(ktfoperq(kf,kaf))then
        select case (kaf)
        case (mtfplus,mtftimes,mtfpower,mtfrevpower)
          nv=0
          do i=isp1+2,isp
            if(ktflistq(dtastk(i),kli) .and. (kli%head%k == kxvect
     $           .or. kli%head%k == kxvect1))then
              dtastk(i)=kli%dbody(1)
              if(tflistq(dtastk(i),kli1))then
                if(nv == 0)then
                  nv=kli1%nl
                elseif(nv /= kli1%nl)then
                  irtc=itfmessage(999,'General::equalleng',
     $                 '"two Vectors"')
                  return
                endif
              endif
            endif
          enddo
          kx1=tfefunrefu(isp1+1,irtc)
          if(irtc /= 0)then
            return
          endif
c          call tfdebugprint(kx1,'vectorize',1)
          if(tflistq(kx1,klx1) .and. klx1%nl == nv)then
            go to 9000
          endif
          kx=kx1
          return
        case (mtfset,mtfsetdelayed)
          do i=isp1+2,isp
            if(ktflistq(ktastk(i),kli) .and.
     $           kli%head%k == kxvect1)then
              k1=kli%dbody(1)
              dtastk(i)=kxadaloc(-1,1,kli1)
              kli1%head%k=ktfcopy1(kxvect)
              kli1%dbody(1)=dtfcopy(k1)
            endif
          enddo
          kx=tfefunrefd(isp1+1,irtc)
          return
        case (mtfend)
          kaf1=klist(ifunbase+kaf)+1
          if(ilist(1,kaf1) == 1 .and. ilist(1,kaf1+1) /= 0 .and.
     $         isp == isp1+2 .and. ktflistq(ktastk(isp),kl) .and.
     $         (kl%head%k == kxvect .or. kl%head%k == kxvect1))then
            dtastk(isp)=kl%dbody(1)
            if(tflistq(ktastk(isp)))then
              kx1=tfefunrefd(isp1+1,irtc)
              go to 9000
            else
              irtc=itfmessage(999,'General::wrongtype',
     $             "Vector[{...}]")
              return
            endif
          endif
        end select
      endif
      nv=0
      do i=isp1+2,isp
        ki=dtastk(i)
        if(ktflistq(ki,kli) .and.
     $       (kli%head%k == kxvect .or. kli%head%k == kxvect1))then
          k1=kli%dbody(1)
          if(tflistq(k1,kli1))then
            if(nv .gt. 0)then
              if(nv /= kli1%nl)then
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
      if(nv == 0)then
        kx=tfefunrefd(isp1+1,irtc)
        return
      endif
      isp2=isp
      isp0=isp+nv+1
      idsp=nv+isp2-isp1
      do j=1,nv
        dtastk(isp0)=kf
        do i=isp1+2,isp2
          ki=dtastk(i)
          if(ktfrefq(ki,kai))then
            ktastk(idsp+i)=klist(kai+j)
          else
            ktastk(idsp+i)=ktastk(i)
          endif
        enddo
        isp=isp2+idsp
        dtastk(isp2+j)=tfefunrefu(isp0,irtc)
        if(irtc /= 0)then
          return
        endif
      enddo
      isp=isp0-1
      kx1=kxmakelist(isp2)
      isp=isp2
      irtc=0
 9000 if(irtc /= 0)then
        return
      endif
      kx=kxadaloc(-1,1,klx)
      klx%head%k=ktfcopy1(kxvect)
      klx%dbody(1)=dtfcopy(kx1)
      irtc=0
      return
      end
