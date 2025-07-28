      subroutine tgetfv(word,nfc,icalc,ncalc,
     1    kfit,fitval,mfitp,ifitp,ifitp1,exist,err)
      use geto
      use ffs
      use ffs_pointer
      use ffs_fit
      use tffitcode
      use tfcsi,only:ipoint
      implicit none
      integer*4 ,intent(inout):: kfit(*),mfitp(*),ifitp(*),ifitp1(*),icalc(3,maxcond)
      type (sad_descriptor),intent(inout)::  fitval(*)
      character*(*) ,intent(inout):: word
      integer*4 ,intent(inout):: nfc,ncalc
      logical*4 ,intent(out):: exist,err
      type (sad_descriptor) kx1
      type (sad_rlist) ,pointer::kl1
      type (sad_dlist) ,pointer::kl1d
      integer*4 i,l,lenw,ix2,kp,j,next,lw,lnf,lne,ir
      real*8 sc,x1,x,getva,sig
      character*8 name1
      character*(MAXPNAME) tname
      character*(MAXPNAME+8) name4,name3,namee
      character peekch,ch
      logical*4 maxfit,err1,exist1,realv,obj
      exist=.true.
      err=.false.
      lw=lenw(word)
      do i=1,mfit1
        l=lenw(nlist(i))
c        write(*,*)'tgetfv ',i,word(:lw),nlist(i)(:l)
        if(word(:l) /= nlist(i)(:l))then
          cycle
        endif
        name1=nlist(i)
        name3=name1
        name4=name1
        name1(l+1:l+1)='M'
        name3(l+1:l+1)='I'
        name4(l+1:l+5)='SCALE'
        if(word(:lw) == nlist(i) .or. word(:l+1) == name1)then
          maxfit=word == name1
          if(mfpnt /= mfpnt1)then
            if(i > mfit .or. i == mfittrx .or. i == mfittry)then
              call termes('Invalid function for range fitting: ',word)
              err=.true.
              return
            endif
          endif
          if((i == mfittrx .or. i == mfittry) .and. mfpnt /= nlat)then
            call termes('Trace fit is only valid at the end of line.',' ')
            err=.true.
            return
          endif
          sc=scale(i)
 11       continue
          realv=.false.
          obj=.false.
          ch=peekch(next)
          if(peekch(next) == ':')then
            ipoint=next
            sc=1.d0
            go to 11
          elseif(ch == '*')then
            ipoint=next
            word='*'
            x1=1.d0
            kx1%x(1)=1.d0
            realv=.true.
          elseif(ch == '@')then
            ipoint=next
            ch=peekch(next)
            if(ch == '-')then
              ipoint=next
              sig=-1.d0
              word='@-'
            elseif(ch == '*')then
              ipoint=next
              select case (i)
              case (mfitax,mfitay,mfitaz,mfitepx,mfitepy,
     $             mfitzpx,mfitzpy,mfitr2,mfitr3,
     $             mfitdpx,mfitdpy,
     $             mfitpepx,mfitpepy,mfitpzpx,mfitpzpy)
                sig=direlc(mfpnt)
              case default
                sig=1.d0
              end select
              word='@*'
            else
              sig=1.d0
              word='@'
            endif
            if(idtypecx(mfpnt) == icMARK)then
              if(i <= ntwissfun)then
                realv=.true.
                x1=rlist(idvalc(mfpnt)+i)*sig
                kx1%x(1)=x1
              else
                call termes('No Marked value for '//nlist(i),
     1               ' at '//tname(mfpnt))
                err=.true.
                return
              endif
            else
              call termes('No MARK element at ',tname(mfpnt))
              err=.true.
              return
            endif
          else
            ir=itfgeto(kx1)
            if(ir /= 0)then
              exist=.false.
              return
            endif
            obj=.true.
            realv=ktfrealq(kx1,x1)
            if(realv)then
              x1=x1*sc
            elseif(tfreallistqd(kx1,kl1))then
              if(kl1%nl == 1)then
                realv=.true.
                x1=kl1%rbody(1)*sc
              elseif(kl1%nl == 2)then
                if(kl1%rbody(1) == kl1%rbody(2))then
                  realv=.true.
                  x1=kl1%rbody(1)*sc
                else
                  call descr_dlist(kx1,kl1d)
                  kl1d=>tfclonelist(kl1d)
                  kl1d%rbody(1)=kl1d%rbody(1)*sc
                  kl1d%rbody(2)=kl1d%rbody(2)*sc
                  kx1=dlist_descr(kl1d)
                endif
              elseif(kx1%k == dxnull%k)then
                obj=.false.
                exist=.false.
                return
              else
                obj=.false.
                exist=.false.
                return
              endif
            else
              obj=.false.
              exist=.false.
              return
            endif
          endif
          exist=.true.
          ix2=int(max(-1.d0,getva(exist1)))
          do j=1,nfc
c            write(*,*)'tgetfv ',j,ifitp(j),ifitp1(j),mfpnt,mfpnt1,nfc
            if(ifitp(j) == mfpnt .and. ifitp1(j) == mfpnt1)then
              if(kfit(j) == i)then
                kp=j
                go to 111
              endif
            endif
          enddo
          if(word == '*')then
            call termes('No default value of ',nlist(i))
            err=.true.
            return
          endif
          if(nfc == maxcond)then
            if(ix2 .lt. 0)then
              call termes(' The value of '//nlist(i),' is not stored.')
            else
              do j=1,nfc
                if(mfitp(j) == 0)then
                  kp=j
                  ifitp(j)=min(mfpnt,mfpnt1)
                  ifitp1(j)=max(mfpnt,mfpnt1)
                  kfit(j)=i
                  go to 111
                endif
              enddo
              call termes('Too many conditions ',nlist(i))
            endif
            err=.true.
            return
          endif
          nfc=nfc+1
          kp=nfc
          ifitp(kp)=min(mfpnt,mfpnt1)
          ifitp1(kp)=max(mfpnt,mfpnt1)
          kfit(kp)=i
          mfitp(kp)=0
111       if(realv .or. obj)then
            if((i == mfitbx .or. i == mfitby) .and. x1 <= 0.d0)then
              lnf=len_trim(nlist(i))
              call elname(ifitp(kp),namee)
              lne=len_trim(namee)
              call termes('Non-positive value for',
     $             nlist(i)(1:lnf)//" @ "//namee(1:lne))
              err=.true.
              return
            endif
            if(realv)then
              call tflocald(fitval(kp))
              fitval(kp)%x(1)=x1
            elseif(obj)then
              call tfdebugprint(kx1,'tgetfv',1)
              call tflocald(fitval(kp))
              fitval(kp)=dtfcopy(kx1)
            endif
          else
            exist=.false.
            return
          endif
          if(ix2 /= 0)then
            mfitp(kp)=ix2+1
          elseif(mfitp(kp) == 0)then
            mfitp(kp)=2
          endif
          mfitp(kp)=merge(-abs(mfitp(kp)),abs(mfitp(kp)),maxfit)
          if(mfitp(kp) == 0 .and. mfpnt /= mfpnt1)then
            call txcalc(icalc,ncalc,mfpnt,mfpnt1,i,.false.,err1)
          elseif(mfitp(kp) /= 0 .and. i <= mfit .and.
     1       (mfpnt == mfpnt1 .or. .not. maxfit
     1        .or. i == mfitnx .or. i == mfitny .or.
     1        (i .ge. mfitleng .and. i <= mfitchi3)))then
            call txcalc(icalc,ncalc,mfpnt,mfpnt1,i,.true.,err1)
          endif
          return
        elseif(word(:lw) == name3)then
          x=getva(exist)
          if(exist)then
            if((i == mfitbx .or. i == mfitby) .and. x <= 0.d0)then
              call termes('Zero or negative value for ',word)
            else
              rlist(latt(1)+i)=x
            endif
          else
            call termes('?Missing number for ',word)
          endif
          return
        elseif(word(:lw) == name4)then
          x=getva(exist)
          if(exist)then
            if(x <= 0.d0)then
              call termes('Zero or negative scale ',word)
            else
              scale(i)=x
            endif
          else
            call termes('Missing scale ',word)
          endif
          return
        endif
      enddo
      exist=.false.
      return
      end

      character*(*) function tname(i)
      use tmacro, only:nlat
      use ffs_pointer, only:pnamec
      implicit none
      integer*4 i
      if(i /= nlat)then
        tname=pnamec(i)
      else
        tname='$$$'
      endif
      return
      end
