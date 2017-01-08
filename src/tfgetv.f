      subroutine tfgetv(word,ntouch,lfno,next,exist)
      use tfstk
      use ffs
      use ffs_pointer
      use tffitcode
      use tfcsi,only:cssetp
      use tflinepcom
      implicit none
      integer*8 kav
      integer*4 ii,i,id,iv,next,lfno,j,ivi,kv,next1,lw1
      real*8 v,getva,va,vx
      integer*4 ntouch,nl,irtc,lw,iii,lenw
      character*(*) word
      character*128 word1
      logical*4 exist,get,rel,maxf,minf,
     $     abbrev,var,exist1,diff,vcomp, cont

c     Initialize to avoid compiler warning
      if(iflinep .eq. 0)then
        call tfinitlinep(irtc)
      endif
      v=0
      kv=-1
      iv=-1
      id=0
c
      vcomp=index(word,'.') .gt. 0
      cont=.false.
 1    exist=.false.
      get=.true.
      rel=.false.
      maxf=.false.
      minf=.false.
      diff=.true.
      id=0
      lw=len_trim(word)
      call tfgetlineps(word,lw,nl,kav,1,irtc)
      if(irtc .ne. 0 .or. nl .le. 0)then
        go to 9000
      endif
      exist=.true.
      LOOP_II: do iii=1,nl
        i=int(rlist(kav+iii))
        ii=klp(i)
        if(get)then
          id=idtypec(ii)
          get=.false.
          iv=ival(i)
          kv=0
          var=.true.
          call peekwd(word1,next1)
          lw1=lenw(word1)
          if(lw1 .gt. MAXPNAME)then
            go to 912
          endif
          if(abbrev(word1,'R_ELATIVE','_'))then
            rel=.true.
          elseif(word1 .eq. 'MAX')then
            maxf=.true.
          elseif(word1 .eq. 'MIN')then
            minf=.true.
          elseif(word1 .eq. 'MAXABS' .or. word1 .eq. 'ABSMAX' .or.
     1           word1 .eq. 'MAXMIN' .or. word1 .eq. 'MINMAX')then
            maxf=.true.
            minf=.true.
          elseif(word1 .eq. ' ' .or. word1 .eq. '-')then
            next=next1
            call cssetp(next)
            exist=.true.
            if(cont)then
              go to 9000
            else
              call termes(lfno,'?Missing value for ',word)
              return
            endif
          else
            kv=itftypekey(id,word1,lw1)
            if(kv .eq. 0)then
              go to 912
            else
              iv=kv
            endif
          endif
          next=next1
          call cssetp(next)
          exist=.true.
 912      if(iv .eq. 0)then
            call termes(lfno,'?No default keyword for ',word)
            return
          endif
          v=getva(exist1)
          if(.not. exist1)then
            if(cont)then
              if(.not. vcomp)then
                call tffsadjust(ntouch)
              endif
            else
              call termes(lfno,'?Missing value for ',word)
            endif
            return
          endif
        endif
        if(idtypec(ii) .ne. id)then
          kv=itftypekey(idtypec(ii),word1,lw1)
          if(kv .ne. 0)then
            ivi=kv
          else
            ivi=ival(i)
          endif
          if(ivi .eq. 0)then
            cycle LOOP_II
          endif
          if(diff)then
            call termes(lfno,
     1           'Info-Different types of elements match ',
     $           word(1:lw)//" "//word1(1:lw1))
            diff=.false.
          endif
        else
          ivi=iv
        endif
        if(minf .or. maxf)then
          ivi=ival(i)
        endif
        if(ivi .eq. 0)then
          call termes(lfno,'?No default keyword for ',word)
          return
        endif
        var=ivi .eq. ival(i)
        if(rel)then
          va=rlist(idvalc(ii)+ivi)*(1.d0+v)
        else
          va=v
        endif
        if(var)then
          vx=rlist(latt(ii)+ivi)/errk(1,ii)
          if(minf)then
            if(maxf)then
              vlim(i,1)=-abs(va)
              vlim(i,2)=abs(va)
              vx=min(abs(va),max(-abs(va),vx))
            else
              vlim(i,1)=va
              vx=max(va,vx)
            endif
          elseif(maxf)then
            vlim(i,2)=va
            vx=min(va,vx)
          else
            vx=va
          endif
          rlist(latt(ii)+ivi)=vx*errk(1,ii)
        else
          vx=va
          rlist(latt(ii)+ivi)=vx
          if(.not. vcomp)then
            do j=1,ntouch
              if(itouchele(j) .eq. i .and. itouchv(j) .eq. ivi)then
                cycle LOOP_II
              endif
            enddo
            ntouch=ntouch+1
            itouchele(ntouch)=i
            itouchv(ntouch)=ivi
          endif
        endif
      enddo LOOP_II
      if(exist)then
        cont=.true.
        go to 1
      endif
 9000 if(exist .and. .not. vcomp)then
        call tffsadjust(ntouch)
      endif
      return
      end
