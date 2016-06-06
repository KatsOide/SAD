      subroutine tfgetv(word,ntouch,lfno,exist)
      use tfstk
      use ffs
      use ffs_pointer
      use tffitcode
      implicit none
      integer*8 kav
      integer*4 ii,i,id,iv,next,k,lfno,j,ivi,kv,kk
      real*8 v,getva,va,vx
      integer*4 ntouch,nl,irtc,lw, iii
      character*(*) word
      character*(MAXPNAME) keywrd,tfkwrd1
      character*64 word1
      logical*4 exist,get,rel,maxf,minf,unitf,
     $abbrev,var,exist1,diff,vcomp, cont

c     Initialize to avoid compiler warning
      v=0
      kv=-1
      iv=-1
c
 2    vcomp=index(word,'.') .gt. 0
      cont=.false.
 1    exist=.false.
      get=.true.
      rel=.false.
      maxf=.false.
      minf=.false.
      unitf=.false.
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
          id=idtype(latt(1,ii))
          get=.false.
          iv=ival(i)
          kv=0
          var=.true.
 101      call peekwd(word1,next)
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
          else
            do k=1,1000
              keywrd=tfkwrd1(id,k,kk)
              if(keywrd .eq. ' ')then
                go to 909
              elseif(keywrd .eq. word1)then
                iv=k
                kv=kk
                go to 912
              endif
            enddo
 909        if(cont)then
              if(word1 .eq. ' ')then
                exist=.true.
                go to 9000
              else
                if(.not. vcomp)then
                  call tffsadjust(ntouch)
                endif
                call cssetp(next)
                word=word1
                go to 2
              endif
            endif
            v=getva(exist1)
            if(.not. exist1)then
              call termes(lfno,'?Missing value for ',word)
              return
            endif
            go to 913
          endif
 912      call cssetp(next)
          cont=.false.
          go to 101
        endif
 913    if(idtype(latt(1,ii)) .ne. id)then
          if(diff)then
            call termes(lfno,
     1           'Info-Different types of variables match ',word)
            diff=.false.
          endif
          if(kv .ne. 0)then
            ivi=kytbl(kv,idtype(latt(1,ii)))
          else
            ivi=ival(i)
          endif
          if(ivi .eq. 0)then
            cycle LOOP_II
          endif
        else
          ivi=iv
        endif
        if(minf .or. maxf)then
          ivi=ival(i)
        endif
        var=ivi .eq. ival(i)
        if(rel)then
          va=rlist(idval(latt(1,ii))+ivi)*(1.d0+v)
        else
          va=v
        endif
        if(var)then
          vx=rlist(latt(2,ii)+ivi)/errk(1,ii)
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
          rlist(latt(2,ii)+ivi)=vx*errk(1,ii)
        else
          vx=va
          rlist(latt(2,ii)+ivi)=vx
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
