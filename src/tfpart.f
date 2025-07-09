      module part
      use tfstk

      contains
      recursive function tfpart1(list,isp1,isp2,err,irtc)
     $     result(kx)
      implicit none
      type (sad_descriptor) kx
      integer*4 ,intent(in):: isp1,isp2
      integer*4 ,intent(out):: irtc
      type (sad_descriptor) kl,ki
      type (sad_dlist) ,intent(in):: list
      type (sad_dlist), pointer :: listi,listl
      integer*4 narg,ivi,iv,ma,m,i,isp0,itfmessage,itfmessageexp
      logical*4 ,intent(in):: err
      narg=isp2-isp1
      ma=list%nl
      kl=dtastk(isp1+1)
      isp0=isp
      irtc=0
      kx%k=ktfoper+mtfnull
      if(ktfrealq(kl,iv))then
        if(iv .lt. 0)then
          iv=ma+iv+1
          if(iv .lt. 0)then
            ivi=iv
            go to 9030
          endif
        elseif(iv .gt. ma)then
          ivi=iv
          go to 9030
        endif
        kx=list%dbody(iv)
        if(narg /= 1)then
          if(ktflistq(kx,listi))then
            kx=tfpart1(listi,isp1+1,isp2,err,irtc)
          else
            go to 9040
          endif
        endif
        return
      endif
      if(tflistq(kl,listl))then
        if(ktfnonreallistqo(listl))then
          if(err)then
            irtc=itfmessage(9,'General::wrongtype','"List of reals as index"')
          else
            irtc=-1
          endif
          return
        endif
        m=listl%nl
        if(narg .gt. 1)then
          do i=1,m
            ivi=int(listl%rbody(i))
            if(ivi .lt. 0)then
              ivi=ma+ivi+1
              if(ivi .lt. 0)then
                go to 9030
              endif
            elseif(ivi .gt. ma)then
              go to 9030
            endif
            ki=list%dbody(ivi)
            if(ktflistq(ki,listi))then
              isp=isp+1
              dtastk(isp)=tfpart1(listi,isp1+1,isp2,err,irtc)
              if(irtc /= 0)then
                isp=isp0
                return
              endif
            else
              go to 9040
            endif
          enddo
        else
          do i=1,m
            ivi=int(listl%rbody(i))
            if(ivi .lt. 0)then
              ivi=ma+ivi+1
              if(ivi .lt. 0)then
                go to 9030
              endif
            elseif(ivi .gt. ma)then
              go to 9030
            endif
            dtastk(isp0+i)=list%dbody(ivi)
          enddo
          isp=isp0+m
        endif
        kx=kxmakelist(isp0)
        isp=isp0
      elseif(kl%k == ktfoper+mtfnull)then
        if(narg == 1)then
          kx=sad_descr(list)
        else
          if(ktfnonreallistqo(list))then
            do i=1,ma
              if(ktflistq(list%dbody(i),listi))then
                isp=isp+1
                dtastk(isp)=tfpart1(listi,isp1+1,isp2,err,irtc)
                if(irtc /= 0)then
                  isp=isp0
                  return
                endif
              else
                go to 9040
              endif
            enddo
          endif
          kx=kxmakelist(isp0)
          isp=isp0
        endif
      elseif(ierrorexp .gt. 0)then
        irtc=-1
      elseif(err)then
        irtc=itfmessage(9,'General::wrongtype','"Real, List of Reals, or Null as index"')
      else
        irtc=-1
      endif
      return
 9030 if(err)then
        irtc=itfmessageexp(9,'General::index',sad_descr(dble(ivi)))
      else
        irtc=-1
      endif
      isp=isp0
      return
 9040 if(err)then
        irtc=itfmessage(9,'General::toomany','"indices"')
      else
        irtc=-1
      endif
      isp=isp0
      return
      end

      recursive subroutine tfpartrstk(lar,isp1,isp2,
     $     list,last,write,eval,err,irtc)
      implicit none
      integer*4 ,intent(in):: isp1,isp2
      integer*4 ,intent(out):: irtc
      type (sad_descriptor) ki,kl,k1
      type (sad_dlist) ,intent(inout):: lar
      type (sad_dlist),pointer :: lari,listl
      integer*8 ka
      integer*4 narg,ivi,iv,ma,m,i,isp0,itfmessage,
     $     itfmessageexp
      logical*4 ,intent(out):: eval,list
      logical*4 ,intent(in):: err,last,write
      logical*4 list1,eval1
      irtc=0
      eval=.false.
      narg=isp2-isp1
      ma=lar%nl
      if(narg .le. 0)then
        return
      endif
      kl=dtastk(isp1+1)
      if(ktfrealq(kl,iv))then
        if(iv .lt. 0)then
          iv=ma+iv+1
          if(iv .lt. 0)then
            ivi=iv
            go to 9030
          endif
        elseif(last .and. iv == ma+1)then
        elseif(iv .gt. ma)then
          ivi=iv
          go to 9030
        endif
        if(narg .gt. 1)then
          k1=lar%dbody(iv)
          if(ktflistq(k1,lari))then
            if(write)then
              lari=>tfclonelist(lari)
              call tfreplist(lar,iv,sad_descr(lari),eval1)
              eval=eval .or. eval1
            endif
            call tfpartrstk(lari,isp1+1,isp2,list,last,write,eval1,err,irtc)
            eval=eval .or. eval1
          else
            go to 9040
          endif
        elseif(write)then
          isp=isp+1
          ktastk(isp)=sad_loc(lar%head)
          itastk2(1,isp)=iv
          list=.false.
        else
          isp=isp+1
          dtastk(isp)=lar%dbody(iv)
        endif
      else
        list=.true.
        isp0=isp
        if(tflistq(kl,listl))then
          if(ktfnonreallistqo(listl))then
            if(err)then
              irtc=itfmessage(9,'General::wrongtype',
     $             '"List of reals as index"')
            else
              irtc=-1
            endif
            return
          endif
          m=listl%nl
          if(narg .gt. 1)then
            do i=1,m
              ivi=int(listl%rbody(i))
              if(ivi .lt. 0)then
                ivi=ma+ivi+1
                if(ivi .lt. 0)then
                  isp=isp0
                  go to 9030
                endif
              elseif(ivi .gt. ma)then
                isp=isp0
                go to 9030
              endif
              ki=lar%dbody(ivi)
              if(ktflistq(ki,lari))then
                if(write)then
                  lari=>tfclonelist(lari)
                  call tfreplist(lar,ivi,sad_descr(lari),eval1)
                  eval=eval .or. eval1
                endif
                call tfpartrstk(lari,isp1+1,isp2,
     $               list1,last,write,eval1,err,irtc)
                if(irtc /= 0)then
                  return
                endif
                eval=eval .or. eval1
              else
                isp=isp0
                go to 9040
              endif
            enddo
          elseif(write)then
            do i=1,m
              ivi=int(listl%rbody(i))
              if(ivi .lt. 0)then
                ivi=ma+ivi+1
                if(ivi .lt. 0)then
                  isp=isp0
                  go to 9030
                endif
              elseif(last .and. ivi == ma+1)then
              elseif(ivi .gt. ma)then
                isp=isp0
                go to 9030
              endif
              isp=isp+1
              ktastk(isp)=sad_loc(lar%head)
              itastk2(1,isp)=ivi
            enddo
          else
            do i=1,m
              ivi=int(listl%rbody(i))
              if(ivi .lt. 0)then
                ivi=ma+ivi+1
                if(ivi .lt. 0)then
                  isp=isp0
                  go to 9030
                endif
              elseif(last .and. ivi == ma+1)then
              elseif(ivi .gt. ma)then
                isp=isp0
                go to 9030
              endif
              isp=isp+1
              dtastk(isp)=lar%dbody(ivi)
            enddo
          endif
        elseif(kl%k == ktfoper+mtfnull)then
          if(narg == 1)then
            if(write)then
              ka=sad_loc(lar%head)
              ktastk(isp+1:isp+ma)=ka
              itastk2(1,isp+1:isp+ma)=[(i,i=1,ma)]
              isp=isp+ma
            else
              call tfgetllstkall(lar)
            endif
          else
            do i=1,ma
              ki=lar%dbody(i)
              if(ktflistq(ki,lari))then
                if(write)then
                  lari=>tfclonelist(lari)
                  call tfreplist(lar,i,sad_descr(lari),eval1)
                  eval=eval .or. eval1
                endif
                call tfpartrstk(lari,isp1+1,isp2,list1,
     $               last,write,eval1,err,irtc)
                if(irtc /= 0)then
                  isp=isp0
                  return
                endif
                eval=eval .or. eval1
              else
                isp=isp0
                go to 9040
              endif
            enddo
          endif
        elseif(err)then
          irtc=itfmessage(9,'General::wrongtype','"Real, List of Reals, or Null as index"')
        else
          irtc=-1
        endif
      endif
      return
 9030 if(err)then
        irtc=itfmessageexp(9,'General::index',sad_descr(dble(ivi)))
      else
        irtc=-1
      endif
      return
 9040 if(err)then
        irtc=itfmessage(9,'General::toomany','"indices"')
      else
        irtc=-1
      endif
      return
      end

      recursive subroutine tfreprulestk(kn,irtc)
      implicit none
      type (sad_descriptor) ,intent(in):: kn
      type (sad_dlist) , pointer :: knl
      integer*4 ,intent(out):: irtc
      integer*4 itfmessage,i
      if(ktfnonlistq(kn,knl))then
        go to 9000
      endif
      irtc=0
      select case (knl%head%k)
        case (ktfoper+mtflist)
          do i=1,knl%nl
            call tfreprulestk(knl%dbody(i),irtc)
            if(irtc /= 0)then
              return
            endif
          enddo
        case (ktfoper+mtfrule,ktfoper+mtfruledelayed)
          if(knl%nl /= 2)then
            go to 9000
          endif
          isp=isp+1
          dtastk(isp)=knl%dbody(1)
          dtastk2(isp)=knl%dbody(2)
        case default
          go to 9000
      end select
      return
 9000 irtc=itfmessage(9,'General::wrongtype',
     $     '"Rule or List of rules for #2"')
      return
      end

      subroutine tfflattenstk(kl,level0,kh,irtc)
      implicit none
      type (sad_descriptor) ,intent(in):: kh
      type (sad_dlist), intent(in):: kl
      type (sad_dlist), pointer :: list,kli
      integer*4 ,intent(in):: level0
      integer*4 ,intent(out):: irtc
      integer*4 level,i,m,i0,itfmessage,mstk0
      call dlist_dlist(kl,list)
      mstk0=mstk
      level=level0
      i0=1

      irtc=0
 1    m=list%nl
      if(isp+m-i0+1 .gt. mstk)then
        mstk=mstk0
        irtc=itfmessage(999,'General::stack',' ')
        return
      endif
      if(ktfreallistq(list))then
        call tfcopyarray(list%dbody(i0:m),dtastk(isp+i0:isp+m),m-i0+1)
c        dtastk(isp+i0:isp+m)=list%dbody(i0:m)
        isp=isp+m-i0+1
      else
        LOOP_I: do i=i0,m
          isp=isp+1
          dtastk(isp)=list%dbody(i)
          if(ktflistq(ktastk(isp),kli))then
            if(tfsameq(kli%head,kh))then
              if(level /= 0)then
                go to 100
              endif
              cycle LOOP_I
            endif
            if(kli%head%k == ktfoper+mtfnull)then
              go to 100
            endif
          endif
        enddo LOOP_I
      endif
      if(mstk /= mstk0)then
        mstk=mstk+1
        call loc_dlist(ktastk(mstk),list)
        i0=itastk2(1,mstk)+1
        level=level+1
        go to 1
      endif
      return
 100  isp=isp-1
      ktastk(mstk)=ksad_loc(list%head%k)
      itastk2(1,mstk)=i
      mstk=mstk-1
      level=level-1
      list=>kli
      i0=1
      go to 1
      end

      subroutine tfreplist(list,iv,k,eval)
      implicit none
      type (sad_descriptor) ,intent(in):: k
      type (sad_dlist) ,intent(inout):: list
      integer*4 ,intent(in):: iv
      integer*4 i,lattr
      logical*4 ,intent(out):: eval
      eval=.false.
      if(iv == 0)then
        call tflocald(list%head)
        list%head=dtfcopy(k)
        list%attr=iand(list%attr,lnonreallist+ktoberebuilt)
      elseif(ktfreallistq(list))then
        if(ktfrealq(k))then
          list%dbody(iv)=k
          list%attr=ior(list%attr,ktoberebuilt)-ktoberebuilt
        else
          list%dbody(iv)=dtfcopy(k)
          eval=ktfsequenceq(k%k)
          list%attr=merge(ktoberebuilt+lnonreallist,lnonreallist,eval)
        endif
      else
        call tflocald(list%dbody(iv))
        if(ktfnonrealq(k))then
          list%dbody(iv)=dtfcopy(k)
          eval=ktfsequenceq(k)
          lattr=iand(list%attr,ktoberebuilt)+lnonreallist
          if(tfconstlistqo(list))then
            if(.not. tfconstq(k%k) .or. eval)then
              list%attr=lattr
            endif
          else
            list%attr=lattr
          endif
        else
          list%dbody(iv)=k
          do i=1,list%nl
            if(ktfnonrealq(list%dbody(i)))then
              list%attr=iand(list%attr,ktoberebuilt)+lnonreallist
              return
            endif
          enddo
          list%attr=iand(list%attr,lconstlist+lnoconstlist)
        endif
      endif
      return
      end

      recursive function tfpartitionstk(isp10,isp20,kl,irtc)
     $     result(kx)
      implicit none
      integer*4 ,intent(in):: isp10,isp20
      integer*4 ,intent(out):: irtc
      type (sad_descriptor) kx,kk,tfefunrefu
      type (sad_descriptor) ,save::kxthread
      type (sad_dlist) ,intent(out):: kl
      type (sad_dlist), pointer :: klj
      integer*8 kai
      integer*4 isp1,ma,ne,no,i,mp,isp2,isp3,isp0,j,itfmessage
      data kxthread%k /0/
      if(kxthread%k == 0)then
        kxthread=kxsymbolz('`System`Thread',14)
      endif
      kx%k=ktfoper+mtfnull
      isp1=isp10
      isp2=isp20
      ma=kl%nl
      ne=itastk(1,isp1)
      if(ne .le. 0)then
        irtc=itfmessage(9,'General::wrongval',
     $       '"#2","positive number or list of them"')
        return
      endif
      no=itastk(2,isp1)
      if(no .le. 0)then
        irtc=itfmessage(9,'General::wrongval',
     $       '"#3","positive number or list of them"')
        return
      endif
      mp=(ma-ne+no)/no
      isp0=isp
      isp=isp+1
      dtastk(isp)=kl%head
      if(mp .gt. 0)then
        if(isp1 == isp2)then
          do i=1,mp
            kai=(i-1)*no
            isp=isp+1
            isp3=isp
            dtastk(isp)=kl%head
            call tfcopyarray(kl%dbody(kai+1:kai+ne),dtastk(isp+1:isp+ne),ne)
c            dtastk(isp+1:isp+ne)=kl%dbody(kai+1:kai+ne)
            isp=isp+ne
            call tfefunrefstk(isp3,isp3,irtc)
            if(irtc /= 0)then
              isp=isp0
              return
            endif
          enddo
        else
          do i=1,mp
            kai=(i-1)*no
            isp=isp+1
            isp3=isp
            dtastk(isp)=kxthread
            isp=isp+1
            dtastk(isp)=kl%head
            do j=1,ne
              isp=isp+1
              dtastk(isp)=kl%dbody(kai+j)
              if(ktflistq(ktastk(isp),klj))then
                kk=tfpartitionstk(isp1+1,isp2,klj,irtc)
                if(irtc /= 0)then
                  isp=isp0
                  return
                endif
                dtastk(isp)=kk
              endif
            enddo
            dtastk(isp3+1)=kxcompose(isp3+1)
            isp=isp3+1
            call tfefunrefstk(isp3,isp3,irtc)
            if(irtc /= 0)then
              isp=isp0
              return
            endif
          enddo
        endif
      endif
      kx=tfefunrefu(isp0+1,irtc)
      isp=isp0
      return
      end

      subroutine tfpartition(isp1,kx,irtc)
      implicit none
      type (sad_descriptor) ,intent(out):: kx
      integer*4 ,intent(in):: isp1
      integer*4 ,intent(out):: irtc
      type (sad_descriptor) kp,k,ks
      type (sad_dlist), pointer :: klp,kl,kls
      integer*4 narg,m,id,itfmessage
      real*8 vp,vs
      narg=isp-isp1
      k=dtastk(isp1+1)
      if(ktfnonlistq(k,kl))then
        irtc=itfmessage(9,'General::wrongtype','"List or composition"')
        return
      endif
      kp=dtastk(isp1+2)
      if(narg == 2)then
        if(ktfrealq(kp,itastk(1,isp+1)))then
          isp=isp+1
          itastk(2,isp)=itastk(1,isp)
          kx=tfpartitionstk(isp,isp,kl,irtc)
          isp=isp1+2
        elseif(tflistq(kp,klp) .and. ktfreallistq(klp))then
          itastk(1,isp+1:isp+klp%nl)=int(klp%rbody(1:klp%nl))
          itastk(2,isp+1:isp+klp%nl)=itastk(1,isp+1:isp+klp%nl)
          isp=isp+klp%nl
c          do i=1,klp%nl
c            isp=isp+1
c            itastk(1,isp)=int(klp%rbody(i))
c            itastk(2,isp)=itastk(1,isp)
c          enddo
          kx=tfpartitionstk(isp1+3,isp,kl,irtc)
          isp=isp1+2
        else
          irtc=itfmessage(9,'General::wrongtype',
     $         '"Real or List of Reals for index"')
        endif
        return
      elseif(narg == 3)then
        ks=dtastk(isp)
        if(ktfrealq(kp,vp) .and. ktfrealq(ks,vs))then
          isp=isp+1
          itastk(1,isp)=int(vp)
          itastk(2,isp)=int(vs)
          kx=tfpartitionstk(isp,isp,kl,irtc)
          isp=isp1+3
          return
        elseif(tflistq(kp,klp) .and. ktfreallistq(klp))then
          m=klp%nl
          if(tflistq(ks,kls) .and. ktfreallistq(kls) .and.
     $         kls%nl == m)then
            itastk(1,isp+1:isp+m)=int(klp%rbody(1:m))
            itastk(2,isp+1:isp+m)=int(kls%rbody(1:m))
            isp=isp+m
c            do i=1,m
c              isp=isp+1
c              itastk(1,isp)=int(klp%rbody(i))
c              itastk(2,isp)=int(kls%rbody(i))
c            enddo
          elseif(ktfrealq(ks,id))then
            itastk(1,isp+1:isp+m)=int(klp%rbody(1:m))
            itastk(2,isp+1:isp+m)=id
            isp=isp+m
c            do i=1,m
c              isp=isp+1
c              itastk(1,isp)=int(klp%rbody(i))
c              itastk(2,isp)=id
c            enddo
          else
            irtc=itfmessage(9,'General::wrongtype',
     $           '"Real or List of Reals for index"')
          endif
          kx=tfpartitionstk(isp1+4,isp,kl,irtc)
          isp=isp1+3
          return
        endif
        irtc=itfmessage(9,'General::wrongtype',
     $       '"Real or List of Reals for index"')
      else
        irtc=itfmessage(9,'General::narg','"2 or 3"')
      endif
      return
      end

      end module

      subroutine tfflatten(isp1,kx,irtc)
      use tfstk
      use part,only:tfflattenstk
      implicit none
      type (sad_descriptor) ,intent(out):: kx
      integer*4 ,intent(in):: isp1
      integer*4 ,intent(out):: irtc
      type (sad_descriptor) kh,kh0
      type (sad_dlist), pointer :: klx,kl
      integer*4 narg,level,itfmessage
      real*8 ,parameter ::amaxl=1.d9
      narg=isp-isp1
      irtc=-1
      if(narg .gt. 3)then
        irtc=itfmessage(9,'General::narg','"Less than 4"')
        return
      endif
      if(ktfnonlistq(ktastk(isp1+1),kl))then
        irtc=itfmessage(9,'General::wrongtype','"List or composition"')
        return
      endif
      kh=kl%head
      kh0=kh
      if(narg == 1)then
        level=-1
      else
        if(ktfnonrealq(dtastk(isp1+2)))then
          irtc=itfmessage(9,'General::wrongtype','"Real number"')
          return
        endif
        level=int(min(amaxl,rtastk(isp1+2)))
        if(level .le. 0)then
          irtc=itfmessage(9,'General::wrongnum','"Positive"')
          return
        endif
        if(narg == 3)then
          kh=dtastk(isp)
        endif
      endif
      call tfflattenstk(kl,level,kh,irtc)
      kx=kxmakelist(isp1+narg,klx)
      klx%head=dtfcopy(kh0)
      isp=isp1+narg
      return
      end

      function  tfpart(isp1,err,irtc) result(kx)
      use tfstk
      use part,only:tfpart1
      implicit none
      type (sad_descriptor) kx
      integer*4 ,intent(in):: isp1
      integer*4 ,intent(out):: irtc
      type (sad_dlist), pointer ::kl
      integer*4 isp2
      logical*4 ,intent(in):: err
      if(isp .le. isp1)then
        kx=dtastk(isp1)
        irtc=0
        return
      endif
      isp2=isp
      call loc_sad(ktfaddr(ktastk(isp1)),kl)
      kx=tfpart1(kl,isp1,isp2,err,irtc)
      return
      end

      function tfreplacepart(isp1,mode,irtc) result(kx)
      use tfstk
      use part,only:tfreprulestk
      use eeval
      implicit none
      type (sad_descriptor) kx
      integer*4 ,intent(in):: isp1,mode
      integer*4 ,intent(out):: irtc
      type (sad_descriptor) k,kn,kf
      type (sad_dlist), pointer :: klx,list
      integer*4 narg,i,itfmessage,isp0
      logical*4 seq,rep,rule
      kx=dxnullo
      narg=isp-isp1
      if(narg == 1 .and. (mode == 0 .or. mode == 3)
     $     .or. narg == 2 .and. mode == 2)then
        irtc=-1
        return
      endif
      rule=.false.
      if(mode == 0)then
        if(narg == 2)then
          rule=.true.
        elseif(narg /= 3)then
          irtc=itfmessage(9,'General::narg','"3"')
          return
        endif
      elseif(mode .lt. 3)then
        if(narg /= 3)then
          irtc=itfmessage(9,'General::narg','"3"')
          return
        endif
      else
        if(narg /= 2)then
          irtc=itfmessage(9,'General::narg','"2"')
          return
        endif
      endif
      kn=dtastk(isp)
      isp0=isp
      if(rule)then
        call tfreprulestk(kn,irtc)
        if(irtc /= 0)then
          isp=isp0
          return
        elseif(isp == isp0)then
          kx=dtastk(isp0-1)
          return
        endif
      endif
      if(mode == 1)then
        k=dtastk(isp0-1)
        if(ktfnonlistq(k,list))then
          irtc=itfmessage(9,'General::wrongtype',
     $         '"List or composition"')
          return
        endif
        kf=dtfcopy(dtastk(isp0-2))
      else
        k=dtastk(isp1+1)
        if(ktfnonlistq(k,list))then
          irtc=itfmessage(9,'General::wrongtype',
     $         '"List or composition"')
          return
        endif
      endif
      list=>tfclonelist(list)
      if(rule)then
        do i=isp0+1,isp
c          call tfdebugprint(dtastk(i),'reppart',1)
c          call tfdebugprint(dtastk2(i),' -> ',1)
          call tfreplacepart0(mode,list,dtastk(i),dtastk2(i),seq,irtc)
          if(irtc /= 0)then
            isp=isp0
            go to 9000
          endif
        enddo
        isp=isp0
      else
        if(mode == 0)then
          kf=dtastk(isp0-1)
        endif
        call tfreplacepart0(mode,list,kn,kf,seq,irtc)
        if(irtc /= 0)then
          go to 9000
        endif
      endif
      if(seq)then
        call tfrebuildl(list,klx,rep)
      else
        list%attr=ior(list%attr,ktoberebuilt)-ktoberebuilt
        klx=>list
      endif
      kx=tfleval(klx,.true.,irtc)
 9000 if(mode == 1)then
        call tflocald(kf)
      endif
      return
      end

      subroutine tfreplacepart0(mode,list,kn,kf,seq,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) ,intent(in):: kn,kf
      type (sad_descriptor) k1
      type (sad_dlist) ,intent(inout):: list
      type (sad_dlist), pointer :: kln
      integer*4 ,intent(in):: mode
      integer*4 ,intent(out):: irtc
      integer*4 itfmessage,i
      logical*4 ,intent(out):: seq
      logical*4 single
      if(ktfrealq(kn))then
        single=.true.
      elseif(tflistq(kn,kln))then
        k1=kln%dbody(1)
        single=ktfnonlistq(k1)
      else
        irtc=itfmessage(9,'General::wrongtype',
     $       '"Real or List of Reals for index"')
        return
      endif
      if(single)then
        call tfreplacepart1(mode,list,kn,kf,seq,irtc)
      else
        do i=1,kln%nl
          call tfreplacepart1(mode,list,kln%dbody(i),kf,seq,irtc)
          if(irtc /= 0)then
            return
          endif
        enddo
      endif
      return
      end

      subroutine tfreplacepart1(mode,kln,kn,kf,seq,irtc)
      use tfstk
      use part,only:tfpartrstk,tfreplist
      implicit none
      type (sad_descriptor) ,intent(in):: kn,kf
      type (sad_descriptor) ki,kxi,tfefunrefu
      type (sad_dlist) ,intent(inout):: kln
      type (sad_dlist), pointer :: kl,kli,klxi
      integer*4 ,intent(in):: mode
      integer*4 ,intent(out):: irtc
      integer*4 isp0,ivi,itfmessage,isp2,isp3,i
      logical*4 ,intent(out):: seq
      logical*4 list,seq1
      isp0=isp
      if(ktfrealq(kn))then
        isp=isp+1
        dtastk(isp)=kn
      elseif(tflistq(kn,kl))then
        call tfgetllstkall(kl)
        if(isp == isp0)then
          irtc=0
          return
        endif
      else
        isp=isp0
        irtc=itfmessage(9,'General::wrongtype',
     $       '"Real or List of Reals for index"')
        return
      endif
      isp2=isp
      call tfpartrstk(kln,isp0,isp2,list,
     $     mode == 2,.true.,seq,.true.,irtc)
      if(irtc /= 0)then
        go to 9000
      endif
      if(mode == 0)then
        do i=isp2+1,isp
          call loc_sad(ktastk(i),kli)
          call tfreplist(kli,itastk2(1,i),kf,seq1)
          seq=seq .or. seq1
        enddo
      elseif(mode == 1)then
        isp3=isp
        do i=isp2+1,isp3
          call loc_dlist(ktastk(i),kli)
          ivi=itastk2(1,i)
          isp=isp3+1
          dtastk(isp)=kf
          isp=isp+1
          dtastk(isp)=kli%dbody(ivi)
          ki=tfefunrefu(isp3+1,irtc)
          if(irtc /= 0)then
            go to 9000
          endif
          call tfreplist(kli,ivi,ki,seq1)
          seq=seq .or. seq1
          isp=isp3
        enddo
      elseif(mode == 2)then
        isp3=isp
        do i=isp2+1,isp3
          call loc_dlist(ktastk(i),kli)
          ivi=itastk2(1,i)
          if(ivi .gt. kli%nl)then
            if(ivi == 1)then
              kli%nl=1
              call tfreplist(kli,1,dtastk(isp0-1),seq1)
            else
              isp=isp3
              isp=isp+1
              dtastk(isp)=kli%dbody(kli%nl)
              isp=isp+1
              ktastk(isp)=ktastk(isp0-1)
              kxi=kxmakelist(isp3,klxi)
              klxi%head%k=ktfoper+mtfnull
              call tfreplist(kli,kli%nl,kxi,seq1)
c              call tfdebugprint(kxi,'reppart1-kxi',1)
            endif
          else
            isp=isp3+1
            ktastk(isp)=ktastk(isp0-1)
            isp=isp+1
            dtastk(isp)=kli%dbody(ivi)
            kxi=kxmakelist(isp3,klxi)
            klxi%head%k=ktfoper+mtfnull
            call tfreplist(kli,ivi,kxi,seq1)
          endif
          seq=seq .or. seq1
        enddo
      elseif(mode == 3)then
        do i=isp2+1,isp
          call loc_dlist(ktastk(i),kli)
          call tfreplist(kli,itastk2(1,i),dxnull,seq1)
          seq=seq .or. seq1
        enddo
      elseif(mode == 4)then
        isp3=isp
        do i=isp2+1,isp3
          call loc_dlist(ktastk(i),kli)
          ivi=itastk2(1,i)
          ki=kli%dbody(ivi)
          if(ktflistq(ki,kl))then
            call tfgetllstkall(kl)
            kxi=kxmakelist(isp3,klxi)
            klxi%head%k=ktfoper+mtfnull
            isp=isp3
            call tfreplist(kli,ivi,kxi,seq1)
            seq=seq .or. seq1
          endif
        enddo
      endif
 9000 isp=isp0
      return
      end

      subroutine tffirst(isp1,kx,mode,irtc)
      use tfstk
      use eeval
      implicit none
      type (sad_descriptor) ,intent(out):: kx
      integer*4 ,intent(in):: isp1,mode
      integer*4 ,intent(out):: irtc
      type (sad_descriptor) k
      type (sad_dlist), pointer :: kl
      integer*4 itfmessage,m
      if(isp /= isp1+1)then
        if(rlist(iaximmediate) .lt. 0.d0)then
          irtc=-1
        else
          irtc=itfmessage(9,'General::narg','"1"')
        endif
        return
      endif
      k=dtastk(isp)
      if(ktfnonlistq(k,kl))then
        if((iaximmediate) .lt. 0.d0)then
          irtc=-1
        else
          irtc=itfmessage(9,'General::wrongtype','"List"')
        endif
        return
      endif
      m=kl%nl
      if(m .le. max(0,mode))then
        if(rlist(iaximmediate) .lt. 0.d0)then
          irtc=-1
        else
          irtc=itfmessage(9,'General::index','""')
        endif
        return
      endif
      if(mode .lt. 0)then
        kx=kl%dbody(m)
      else
        kx=kl%dbody(mode+1)
      endif
      irtc=0
      kx=tfeevalref(kx,irtc)
      return
      end
