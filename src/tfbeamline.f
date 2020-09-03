      module beamline
      use maccbk, only:MAXPNAME
      implicit none
      character*(MAXPNAME) :: lname=' '
      end module

      subroutine tfbeamline(k,idx,ename,irtc)
      use tfstk
      use mackw
      use maccbk, only:MAXPNAME
      implicit none
      type (sad_descriptor) ,intent(in):: k
      type (sad_descriptor) ki,k1
      type (sad_dlist), pointer :: kl,kli
      integer*8 kdx1
      integer*4 ,intent(out):: idx,irtc
      integer*4 hsrchz,n,lenw,idxi,
     $     i,idir,idti,nc,itfmessage,itfmessagestr
      character*(MAXPNAME) tfgetstrs
      character*(MAXPNAME) ,intent(out):: ename
      type (sad_descriptor) , save:: kxbl=sad_descriptor(1,i00)
      integer*4 ,save :: lid=0
      if(kxbl%k .eq. 0)then
        kxbl=kxsymbolz('BeamLine',8)
      endif
      idx=0
      ename=' '
      if(ktfnonlistq(k,kl))then
        go to 9900
      endif
      if(.not. tfsameq(kl%head,kxbl))then
        go to 9900
      endif
      n=kl%nl
      if(n .lt. 1)then
        go to 9900
      endif
      kdx1=ktaloc(n+1)
      ilist(1,kdx1)=n
      ilist(2,kdx1)=0
      do i=1,n
        ki=kl%dbody(i)
        idir=1
 1      if(ktflistq(ki,kli))then
          if(kli%head%k .eq. ktfoper+mtfmult)then
            k1=kli%dbody(1)
            if(ktfrealq(k1))then
              if(k1%x(1) .eq. -1.d0)then
                idir=-idir
                ki=kli%dbody(2)
                go to 1
              endif
            endif
            irtc=itfmessage(9,'General::wrongval',
     $           '"element or -element is ",""')
            go to 9000
          else
            irtc=itfmessage(9,'General::wrongval',
     $           '"element or -element is ",""')
            go to 9000
          endif
        else
          ename=tfgetstrs(ki%k,nc)
          if(nc .gt. 0)then
            idxi=hsrchz(ename)
            idti=idtype(idxi)
            if(idti .eq. icNULL .or. idti .gt. icMXEL)then
              irtc=itfmessagestr(9,'MAIN::wrongtype',ename(1:nc))
              go to 9000
            elseif(i .eq. 1 .and. idti .ne. icMARK)then
              irtc=itfmessage(9,'FFS::firstmark',' ')
              go to 9000
            endif
            ilist(1,kdx1+i)=idir
            ilist(2,kdx1+i)=idxi
          else
            irtc=itfmessage(9,'General::wrongleng',
     $           '"name of element","nonzero"')
            go to 9000
          endif
        endif
      enddo
      lid=lid+1
      write(ename,'(''L'',i6.6,''$'')')lid
      idx=hsrchz(ename)
      if(idtype(idx) .ne. icNULL)then
        lid=lid-1
        irtc=itfmessagestr(9,'MAIN::exist',ename(1:lenw(ename)))
        go to 9000
      endif
      idtype(idx)=icLINE
      idval(idx)=kdx1
      pname(idx)=ename
      irtc=0
      return
 9000 ilist(1,kdx1)=n+1
      call tfree(kdx1)
      return
 9900 irtc=itfmessage(9,'General::wrongtype','"BeamLine[ ... ]"')
      return
      end

      integer function itfdummyline()
      use kyparam
      use tfstk
      use mackw
      implicit none
      integer*8 kdx1,idxm1,idxd1
      integer*4 itfdummyptr,idx,hsrchz,idxm,n,idxd
      data itfdummyptr /0/
      if(itfdummyptr .eq. 0)then
        idx=hsrchz('$DUMMYLINE')
        idtype(idx)=icLINE
        pname(idx)='$DUMMYLINE'
        kdx1=ktaloc(3)
        ilist(1,kdx1)=2
        ilist(2,kdx1)=0
        idval(idx)=kdx1
        idxm=hsrchz('$DUMMYMARK')
        idtype(idxm)=icMARK
        pname(idxm)='$DUMMYMARK'
        n=ky_MAX_MARK
        idxm1=ktaloc(n+1)
        idval(idxm)=idxm1
        ilist(1,idxm1)=n
        ilist(2,idxm1)=0
        rlist(idxm1+1:idxm1+n-1)=0.d0
        idxd=hsrchz('$DUMMYDRIFT')
        idtype(idxd)=icDRFT
        pname(idxd)='$DUMMYDRIFT'
        n=ky_MAX_DRFT
        idxd1=ktaloc(n+1)
        idval(idxd)=idxd1
        ilist(1,idxd1)=n
        ilist(2,idxd1)=0
        rlist(idxd1+1)=1.d0
        rlist(idxd1+2:idxd1+n-1)=0.d0
        ilist(1,kdx1+1)=1
        ilist(2,kdx1+1)=idxm
        ilist(1,kdx1+2)=1
        ilist(2,kdx1+2)=idxd
        itfdummyptr=idx
      endif
      itfdummyline=itfdummyptr
      return
      end

      function tfsetelement(isp1,irtc) result(kx)
      use tfstk
      use mackw
      use funs
      implicit none
      type (sad_descriptor) kx,kr
      type (sad_dlist), pointer :: klxi,klx
      integer*8 ka1,kas,kdx1,ktcaloc
      integer*4 ,intent(in):: isp1
      integer*4 ,intent(out):: irtc
      integer*4 nc,lenw,narg,idx,itype,
     $     idt,n,i,nce,m,hsrchz,isp0, itfmessage,
     $     itfdownlevel,l,itfmessagestr
      character*(MAXPNAME) ename,type,tfgetstrs,key,tfkwrd
      kx=dxnullo
      ename=tfgetstrs(ktastk(isp1+1),nce)
      if(nce .lt. 0)then
        irtc=itfmessage(9,'General::wrongleng',
     $       '"name of element","nonzero"')
        return
      endif
      narg=isp-isp1
      idx=hsrchz(ename)
      itype=idtype(idx)
      if(narg .le. 1)then
        if(itype .eq. icNULL)then
          type=' '
        else
          type=pname(kytbl(0,itype))(2:)
        endif
      else
        type=tfgetstrs(ktastk(isp1+2),nc)
        call capita(type(1:nc))
 1      if(nc .le. 0)then
          if(ktastk(isp1+2) .eq. ktfoper+mtfnull .or.
     $         ktastk(isp1+2) .eq. kxnulls)then
            type=pname(kytbl(0,itype))(2:)
          elseif(ktfrealq(ktastk(isp1+2)))then
            idt=int(rtastk(isp1+2))
            type=pname(kytbl(0,idt))(2:)
            nc=lenw(type)
            go to 1
          else
            irtc=itfmessage(9,'General::wrongval',
     $           '"#2 is","type-string, number, or Null"')
            return
          endif
        else
          idt=hsrchz('$'//type)
          if(idtype(idt) .ne. icDEF)then
            irtc=itfmessagestr(9,'MAIN::wrongtype',type(1:lenw(type)))
            return
          endif
          if(idval(idt) .eq. icNULL)then
          elseif(itype .eq. icNULL)then
            itype=int(idval(idt))
            idtype(idx)=itype
            if(itype .ne. icNULL)then
              n=kytbl(kwMAX,itype)
              kdx1=ktcaloc(n+1)
              idval(idx)=kdx1
              ilist(2,kdx1-1)=idx
              ilist(1,kdx1)=n
              ilist(2,kdx1)=0
            endif
          elseif(idval(idt) .ne. itype)then
            irtc=itfmessagestr(9,'FFS::equaltype',type(1:lenw(type)))
            return
          endif
        endif
      endif
      ka1=ktfaddr(ktastk(isp1+1))
      kas=merge(ka1,merge(klist(klist(ifunbase+ka1)),klist(ka1),
     $     ktfoperq(ktastk(isp1+1))),ktfstringq(ktastk(isp1+1)))
      if(itype .eq. icNULL)then
        if(narg .gt. 2)then
          irtc=itfmessage(9,'General::narg','"1 or 2"')
          return
        endif
        kx=kxadaloc(-1,2,klx)
        klx%dbody(1)%k=ktfstring+ktfcopy1(kas)
        klx%dbody(2)=dtfcopy1(dxnulls)
      else
        if(isp .gt. isp1+2)then
          kr=tfoverride(isp1+2,irtc)
          if(irtc .ne. 0)then
            return
          endif
          levele=levele+1
          call tfsetelementkey(idx,kr,irtc)
          l=itfdownlevel()
          if(irtc .ne. 0)then
            return
          endif
        endif
        m=kytbl(kwMAX,itype)-1
        kx=kxadaloc(-1,max(2,min(3,2+m)),klx)
        klx%dbody(1)%k=ktfstring+ktfcopy1(kas)
        klx%dbody(2)=kxsalocb(0,type,lenw(type))
        if(m .gt. 0)then
          isp0=isp
          do i=1,m
            key=tfkwrd(itype,i)
            if(key .ne. ' ' .and. key .ne. '-')then
              isp=isp+1
              dtastk(isp)=kxadaloc(-1,2,klxi)
              klxi%head%k=ktfoper+mtfrule
              klxi%dbody(1)=kxsalocb(0,key,lenw(key))
              klxi%dbody(2)=dtfcopy(dlist(idval(idx)+i))
            endif
          enddo
          klx%dbody(3)=dtfcopy1(kxmakelist(isp0))
          isp=isp0
        endif
      endif    
      irtc=0
      return
      end

      recursive subroutine tfsetelementkey(idx,k,irtc)
      use tfstk
      use mackw
      use eeval
      implicit none
      type (sad_dlist), pointer :: kr
      type (sad_dlist), pointer :: kl
      type (sad_descriptor) k,ki,kk,kv
      integer*4 irtc,idx,i,idt,ioff,nc,itfmessage,itfmessagestr
      character*(MAXPNAME) tfgetstrs,key
      if(tflistq(k,kl))then
        do i=1,kl%nl
          ki=kl%dbody(i)
          call tfsetelementkey(idx,ki,irtc)
          if(irtc .ne. 0)then
            return
          endif
        enddo
      elseif(tfruleq(k,kr))then
        kk=kr%dbody(1)
        key=tfgetstrs(kk%k,nc)
        if(nc .le. 0)then
          irtc=itfmessage(9,'General::wrongtype','"Character-string"')
          return
        endif
        idt=idtype(idx)
        do ioff=1,kytbl(kwMAX,idt)-1
          i=kyindex(ioff,idt)
          if(i .ne. 0)then
            if(pname(kytbl(i,0))(2:) .eq. key(1:nc))then
              go to 10
            endif
            i=kyindex1(ioff,idt)
            if(i .ne. 0 .and.
     $           pname(kytbl(i,0))(2:) .eq. key(1:nc))then
              go to 10
            endif
          endif
        enddo
        irtc=itfmessagestr(9,'FFS::undefkey',key(1:nc))
        return
 10     kv=kr%dbody(2)
        if(kr%head%k .eq. ktfoper+mtfruledelayed)then
          kv=tfeevalref(kv,irtc)
          if(irtc .ne. 0)then
            return
          endif
        endif
        if(ktfrealq(kv))then
          call tflocald(dlist(idval(idx)+ioff))
          dlist(idval(idx)+ioff)=kv
        elseif(tfnonlistq(kv))then
          irtc=itfmessage(9,'General::wrongtype',
     $         '"Keyword -> value"')
          return
        else
          call tflocald(dlist(idval(idx)+ioff))
          dlist(idval(idx)+ioff)=dtfcopy(kv)
        endif
      else
        irtc=itfmessage(9,'General::wrongtype',
     $       '"Keyword -> value or List of them"')
        return
      endif
      irtc=0
      return
      end

      subroutine tfbeamlinename(isp1,kx,irtc)
      use beamline, only: lname
      use tfstk
      implicit none
      type (sad_descriptor) kx
      integer*4 isp1,irtc,itfmessage
      if(isp .gt. isp1+1 .or.
     $     ktastk(isp) .ne. ktfoper+mtfnull)then
        irtc=itfmessage(9,'General::wrongtype','"[]"')
        return
      endif
      irtc=0
      kx=kxsalocb(-1,lname,len_trim(lname))
      return
      end

      subroutine tfsetbeamlinename(name)
      use beamline
      implicit none
      character*(*) name
      lname=name
      return
      end

      character*(*) function tfgetbeamlinename()
      use beamline
      implicit none
      tfgetbeamlinename=lname
      return
      end

      function tfextractbeamline(isp1,irtc) result(kx)
      use tfstk
      use sad_main
      use mackw
      use eeval
      implicit none
      type (sad_descriptor) kx
      type (sad_dlist), pointer :: klx,kli
      type (sad_el), pointer ::el
      integer*8 itfilattp,idx
      integer*4 ,intent(in):: isp1
      integer*4 ,intent(out):: irtc
      integer*4 lenw,i,n,hsrchz,idl,nc,itfmessage,itfmessagestr
      character*(MAXPNAME) ename,tfgetstrs
      logical*4 eval
      data eval /.true./
      kx=dxnullo
      if(eval)then
        kx=tfsyeval(kxsymbolz('BeamLine',8),irtc)
        if(irtc .ne. 0)then
          return
        endif
        eval=.false.
      endif
      if(isp .gt. isp1+1)then
        irtc=itfmessage(9,'General::narg','"0 or 1"')
        return
      endif
      if(ktastk(isp) .eq. ktfoper+mtfnull)then
        idx=itfilattp()
      else
        ename=tfgetstrs(ktastk(isp),nc)
        if(nc .lt. 0)then
          irtc=itfmessage(9,'General::wrongtype','"Character-string"')
          return
        elseif(nc .eq. 0 .or. ename .eq. '*')then
          idx=itfilattp()
        else
          idl=hsrchz(ename)
          if(idtype(idl) .ne. icLINE)then
            irtc=itfmessagestr(9,'MAIN::wrongtype',ename(1:lenw(ename)))
            return
          endif
          call expln(idl)
          idx=idval(ilist(2,idval(idl)))
          idx=idval(ilist(2,idval(idl)))
        endif
      endif
      call loc_el(idx,el)
c      write(*,*)'extractbeamline-0 ',idx,idl,el%comp(1),el%comp(2),ename
      n=el%nlat0
      kx=kxadaloc(-1,n,klx)
      klx%head=dtfcopy1(kxsymbolz('BeamLine',8))
      do i=1,n
c        write(*,*)'extractbeamline ',i,n,el%comp(i)
        ename=pname(idcomp(el,i))
        if(dircomp(el,i) .ge. 0.d0)then
          klx%dbody(i)=dtfcopy1(kxsymbolf(ename,lenw(ename),.true.))
        else
          klx%dbody(i)=kxadaloc(0,2,kli)
          kli%head%k=ktfoper+mtftimes
          kli%rbody(1)=-1.d0
          kli%dbody(2)=dtfcopy1(kxsymbolf(ename,lenw(ename),.true.))
        endif
      enddo
      irtc=0
      return
      end

      recursive subroutine tfexpandbeamline(isp1,kx,irtc)
      use tfstk
      implicit none
      type (sad_descriptor) kx,kx1
      type (sad_dlist), pointer :: kl,kll,klx,klx1
      integer*8 kal
      integer*4 ,intent(in):: isp1
      integer*4 ,intent(out):: irtc
      integer*4 i,j,isp0,m,n,isp2
      integer*8 ifbeamline
      data ifbeamline/0/
      if(ifbeamline .eq. 0)then
        call tfevals('BeamLine',kx,irtc)
        if(irtc .ne. 0)then
          return
        endif
        ifbeamline=ktfaddrd(kx)
      endif
      isp0=isp
      LOOP_I: do i=isp1+1,isp
        if(ktflistq(ktastk(i),kl))then
          m=kl%nl
          if(m .eq. 2 .and. kl%head%k .eq. ktfoper+mtftimes)then
            if(ktfnonreallistqo(kl))then
              do j=1,2
                if(ktfrealq(kl%dbody(j)))then
                  if(kl%rbody(j) .ne. -1.d0)then
                    n=int(kl%rbody(j))
c                    write(*,*)'expandbeamline ',j,kl%rbody(j),n
                    if(n .gt. 0)then
                      dtastk(isp+1:isp+n)=kl%dbody(3-j)
                      isp=isp+n
c                      do k=1,n
c                        isp=isp+1
c                        dtastk(isp)=kl%dbody(3-j)
c                      enddo
                    else
                      kal=ktadaloc(-1,2,kll)
                      kll%head%k=ktfoper+mtftimes
                      kll%rbody(1)=-1.d0
                      kll%dbody(2)=dtfcopy(kl%dbody(3-j))
                      ktastk(isp+1:isp-n)=ktflist+kal
                      isp=isp-n
c                      do k=1,-n
c                        isp=isp+1
c                        ktastk(isp)=ktflist+kal
c                      enddo
                    endif
                    cycle LOOP_I
                  endif
                endif
              enddo
            endif
          elseif(kl%head%k .eq. ktfsymbol+ifbeamline)then
            isp2=isp
            call tfgetllstkall(kl)
            call tfexpandbeamline(isp2,kx1,irtc)
            if(irtc .ne. 0)then
              isp=isp0
              return
            endif
            isp=isp2
            call loc_sad(ktfaddrd(kx1),klx1)
            call tfgetllstkall(klx1)
            cycle LOOP_I
          endif
        endif
        isp=isp+1
        ktastk(isp)=ktastk(i)
      enddo LOOP_I
      kx=kxmakelist(isp0,klx)
      isp=isp0
      klx%head%k=ktfsymbol+ktfcopy1(ifbeamline)
      irtc=0
      return
      end
