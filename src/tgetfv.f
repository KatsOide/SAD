      subroutine tgetfv(word,nfc,lfno,icalc,ncalc,
     1    kfit,fitval,mfitp,ifitp,ifitp1,
     $    exist,err)
      use tfstk
      use ffs
      use ffs_pointer
      use ffs_fit
      use tffitcode
      use tfcsi,only:ipoint
      implicit none
      integer*4 nfc,lfno,i,l,lenw,ix2,kp,j,next,ncalc
      real*8 sc,x1,x,getva,sig
      character*8 name1
      character*(MAXPNAME) tname
      character*(MAXPNAME+8) name4,name3,namee
      character*(*) word
      character peekch,ch
      integer*4 kfit(*),mfitp(*),ifitp(*),ifitp1(*),
     $     icalc(3,maxcond),lw,lnf,lne
      real*8 fitval(*)
      logical*4 exist,maxfit,err,err1,exist1
      exist=.true.
      err=.false.
      lw=lenw(word)
      do 10 i=1,mfit1
        l=lenw(nlist(i))
c        write(*,*)'tgetfv ',i,word(:lw),nlist(i)(:l)
        if(word(:l) .ne. nlist(i)(:l))then
          cycle
        endif
        name1=nlist(i)
        name3=name1
        name4=name1
        name1(l+1:l+1)='M'
        name3(l+1:l+1)='I'
        name4(l+1:l+5)='SCALE'
        if(word(:lw) .eq. nlist(i) .or. word(:l+1) .eq. name1)then
          maxfit=word .eq. name1
          if(mfpnt .ne. mfpnt1)then
            if(i .gt. mfit .or.
     $           i .eq. mfittrx .or. i .eq. mfittry)then
              call termes(lfno,
     $             'Invalid function for range fitting: ',word)
              err=.true.
              return
            endif
          endif
          if((i .eq. mfittrx .or. i .eq. mfittry)
     $         .and. mfpnt .ne. nlat)then
            call termes(lfno,
     $           'Trace fit is only valid at the end of line.',' ')
            err=.true.
            return
          endif
          sc=scale(i)
 11       continue
          ch=peekch(next)
          if(peekch(next) .eq. ':')then
            ipoint=next
            sc=1.d0
            go to 11
          elseif(ch .eq. '*')then
            ipoint=next
            word='*'
            x1=1.d0
          elseif(ch .eq. '@')then
            ipoint=next
            ch=peekch(next)
            if(ch .eq. '-')then
              ipoint=next
              sig=-1.d0
              word='@-'
            elseif(ch .eq. '*')then
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
            if(idtypecx(mfpnt) .eq. icMARK)then
              if(i .le. ntwissfun)then
                x1=rlist(idvalc(mfpnt)+i)*sig
              else
                call termes(lfno,'No Marked value for '//nlist(i),
     1               ' at '//tname(mfpnt))
                err=.true.
                return
              endif
            else
              call termes(lfno,'No MARK element at ',
     $             tname(mfpnt))
              err=.true.
              return
            endif
          else
            x1=getva(exist)*sc
            if(.not. exist)then
              return
            endif
          endif
          exist=.true.
          ix2=int(max(-1.d0,getva(exist1)))
          do 110 j=1,nfc
            if(ifitp(j) .eq. mfpnt .and. ifitp1(j) .eq. mfpnt1)then
              if(kfit(j) .eq. i)then
                kp=j
                go to 111
              endif
            endif
110       continue
          if(word .eq. '*')then
            call termes(lfno,'No default value of ',nlist(i))
            err=.true.
            return
          endif
          if(nfc .eq. maxcond)then
            if(ix2 .lt. 0)then
              call termes(lfno,
     1             ' The value of '//nlist(i),' is not stored.')
            else
              do 130 j=1,nfc
                if(mfitp(j) .eq. 0)then
                  kp=j
                  ifitp(j)=min(mfpnt,mfpnt1)
                  ifitp1(j)=max(mfpnt,mfpnt1)
                  kfit(j)=i
                  go to 111
                endif
130           continue
              call termes(lfno,
     $             'Too many conditions ',nlist(i))
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
111       if(word .ne. '*')then
            if((i .eq. mfitbx .or. i .eq. mfitby)
     $           .and. x1 .le. 0.d0)then
              lnf=len_trim(nlist(i))
              call elname(ifitp(kp),namee)
              lne=len_trim(namee)
              call termes(lfno,'Non-positive value for',
     $             nlist(i)(1:lnf)//" @ "//namee(1:lne))
              err=.true.
              return
            endif
            fitval(kp)=x1
          endif
          if(ix2 .ne. 0)then
            mfitp(kp)=ix2+1
          elseif(mfitp(kp) .eq. 0)then
            mfitp(kp)=2
          endif
          mfitp(kp)=merge(-abs(mfitp(kp)),abs(mfitp(kp)),maxfit)
          if(mfitp(kp) .eq. 0 .and. mfpnt .ne. mfpnt1)then
            call txcalc(icalc,ncalc,mfpnt,mfpnt1,i,.false.,err1)
          elseif(mfitp(kp) .ne. 0 .and. i .le. mfit .and.
     1       (mfpnt .eq. mfpnt1 .or. .not. maxfit
     1        .or. i .eq. mfitnx .or. i .eq. mfitny .or.
     1        (i .ge. mfitleng .and. i .le. mfitchi3)))then
            call txcalc(icalc,ncalc,mfpnt,mfpnt1,i,.true.,err1)
          endif
          return
        elseif(word(:lw) .eq. name3)then
          x=getva(exist)
          if(exist)then
            if((i .eq. mfitbx .or. i .eq. mfitby) .and. x .le. 0.d0)then
              call termes(lfno,'Zero or negative value for ',word)
            else
              rlist(latt(1)+i)=x
            endif
          else
            call termes(lfno,'?Missing number for ',word)
          endif
          return
        elseif(word(:lw) .eq. name4)then
          x=getva(exist)
          if(exist)then
            if(x .le. 0.d0)then
              call termes(lfno,'Zero or negative scale ',word)
            else
              scale(i)=x
            endif
          else
            call termes(lfno,'Missing scale ',word)
          endif
          return
        endif
10    continue
      exist=.false.
      return
      end

      character*(*) function tname(i)
      use tmacro, only:nlat
      use ffs_pointer, only:pnamec
      implicit none
      integer*4 i
      if(i .ne. nlat)then
        tname=pnamec(i)
      else
        tname='$$$'
      endif
      return
      end
