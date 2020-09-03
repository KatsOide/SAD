      module funs

      contains
      function tfsequence(isp1,isp2) result(kx)
      use tfstk
      implicit none
      type (sad_descriptor) kx
      integer*4 ,intent(in):: isp1,isp2
      kx=merge(dxnull,merge(dtastk(isp2),
     $     kxcrelistm(isp2-isp1,ktastk(isp1+1:isp2),
     $     k_descr(ktfoper+mtfnull)),isp1+1 .eq. isp2),
     $     isp1 .ge. isp2)
      return
      end
      function tfreplace1(isp1,irtc) result(kx)
      use tfstk
      implicit none
      type (sad_descriptor) kx,k1
      integer*4 ,intent(in):: isp1
      integer*4 ,intent(out):: irtc
      integer*4 i,itfmessage
      kx=dxnullo
      if(isp .le. isp1+1)then
        irtc=itfmessage(9,'General::narg','"2 or more"')
        return
      endif
      k1=dtastk(isp1+1)
      do i=isp1+2,isp
        call tfreplace(k1,dtastk(i),kx,.true.,.true.,.false.,irtc)
        if(irtc .ne. 0)then
          return
        endif
        k1=kx
      enddo
      return
      end

      function tfreplacerepeated1(isp1,irtc) result(kx)
      use tfstk
      implicit none
      type (sad_descriptor) kx
      integer*4 ,intent(in):: isp1
      integer*4 ,intent(out):: irtc
      integer*4 i,itfmessage
      kx=dxnullo
      if(isp .le. isp1+1)then
        irtc=itfmessage(9,'General::narg','"2 or more"')
        return
      endif
      kx%k=ktastk(isp1+1)
      irtc=0
      do i=isp1+2,isp
        call tfreplacerepeated(kx,ktastk(i),kx,.true.,.true.,irtc)
        if(irtc .ne. 0)then
          return
        endif
      enddo
      return
      end

      function tfsameq1(isp1,iopc,irtc) result(kx)
      use tfstk
      implicit none
      type (sad_descriptor) kx,k,k1
      integer*4 ,intent(in):: isp1,iopc
      integer*4 ,intent(out):: irtc
      integer*4 itfmessage
      kx=dxnullo
      if(isp .ne. isp1+2)then
        irtc=itfmessage(9,'General::narg','"2"')
        return
      endif
      k=dtastk(isp)
      k1=dtastk(isp1+1)
      kx%k=0
      if(k%k .eq. k1%k)then
        kx%k=ktftrue
      elseif(tfsameq(k,k1))then
        kx%k=ktftrue
      endif
      if(iopc .eq. mtfunsame)then
        kx%k=ktftrue-kx%k
      endif
      irtc=0
      return
      end

      function tfupset(k1,k2,kas,irtc) result(kx)
      use tfstk
      implicit none
      type (sad_descriptor) ,intent(in):: k1,k2
      type (sad_descriptor) kx,ki,karg
      type (sad_dlist), pointer :: kl,kli
      type (sad_symbol), pointer :: symi
      type (sad_symdef), pointer :: symd
      integer*8 ,intent(in):: kas
      integer*4 ,intent(out):: irtc
      integer*4 i,isp0,isp1,m,itfmessage
      kx=dxnullo
      if(ktfnonlistq(k1,kl))then
        irtc=itfmessage(999,'General::wrongtype','"Expression"')
        return
      endif
      m=kl%nl
      if(m .le. 0)then
        irtc=itfmessage(999,'General::wrongleng',
     $       '"Expression","longer than 0"')
        return
      endif
      isp0=isp
      isp1=isp0+1
      call tfgetllstk(kl,0,-1)
      karg=kxcompose(isp1)
      LOOP_I: do i=isp1+1,isp
        ki=dtastk(i)
        do while(ktflistq(ki,kli))
          ki=kli%head
        enddo
        if(ktfsymbolqdef(ki%k,symd))then
          if(symd%sym%override .ne. 0)then
            if(symd%sym%gen .lt. 0 .and. symd%sym%gen .ne. -3)then
              cycle LOOP_I
            endif
            if(kas .eq. 0 .or. kas .eq. ktfaddr(ki))then
              call tfdset(k2,symd%upval,kx,karg)
              if(kas .ne. 0)then
                cycle LOOP_I
              endif
            endif
          else
            symi=>tfsydef(symd%sym)
            if(symi%gen .lt. 0 .and. symi%gen .ne. -3)then
              cycle LOOP_I
            endif
            if(kas .eq. 0 .or. kas .eq. ksad_loc(symi%loc))then
              call sym_symdef(symi,symd)
              call tfdset(k2,symd%upval,kx,karg)
              if(kas .ne. 0)then
                cycle LOOP_I
              endif
            endif
          endif
        endif
      enddo LOOP_I
      kx=k2
      isp=isp0
      irtc=0
      return
      end

      function tfoverride(isp1,irtc) result(kx)
      use tfstk
      implicit none
      type (sad_descriptor) kx
      type (sad_dlist), pointer :: kli
      integer*8 ki,k1
      integer*4 ,intent(in):: isp1
      integer*4 ,intent(out):: irtc
      integer*4 isp0,isp2,n,itfmessage,isp3,isp4,i,j
      if(isp1 .eq. isp)then
        kx=dxnulll
        irtc=0
        return
      elseif(isp1+1 .eq. isp)then
        if(ktastk(isp) .eq. ktfoper+mtfnull)then
          kx=dxnulll
          irtc=0
          return
        endif
      endif
      kx=dxnullo
      isp=isp+1
      isp0=isp
      do i=isp1+1,isp0-1
        ki=ktastk(i)
        if(ktflistq(ki,kli))then
          k1=kli%head%k
          if(k1 .eq. ktfoper+mtflist)then
            call tfgetllstkall(kli)
          elseif(k1 .eq. ktfoper+mtfrule .or.
     $           k1 .eq. ktfoper+mtfruledelayed)then
            isp=isp+1
            ktastk(isp)=ki
          else
            go to 9000
          endif
        else
          isp=isp+1
          ktastk(isp)=ki
        endif
      enddo
      isp2=isp
      do i=isp0+1,isp2
        isp=isp+1
        dtastk(isp)=merge(kli%dbody(1),dtastk(i),
     $       ktflistq(ktastk(i),kli))
      enddo
      n=isp-isp2
      isp3=isp
      call tfsortl(ktastk(isp2-3),.false.,n,2,ktfref,.true.,irtc)
      if(irtc .ne. 0)then
        isp=isp0-1
        return
      endif
      isp4=isp
      do i=1,n
        j=int(ktastk(isp3+i))
        if(j .ne. 0 .and.
     $       ktastk(isp0+j) .ne. ktfoper+mtfnull)then
          isp=isp+1
          ktastk(isp)=ktastk(isp0+j)
        endif
      enddo
      kx=kxmakelist(isp4)
      isp=isp0-1
      irtc=0
      return
 9000 irtc=itfmessage(9,'General::wrongtype',
     $     '"List, Rule, Symbol, String, Real"')
      isp=isp0-1
      return
      end
      end module
