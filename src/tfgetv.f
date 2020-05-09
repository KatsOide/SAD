      subroutine tfgetv(word,lfno,next,exist)
      use tfstk
      use ffs
      use ffs_pointer
      use tffitcode
      use tfcsi,only:ipoint
      use tflinepcom
      use ffs_seg
      implicit none
      type (sad_comp), pointer :: cmpd
      type (sad_descriptor) kx
      type (sad_rlist),pointer:: kl
      integer*8 kav
      integer*4 ii,i,id,iv,next,lfno,ivi,kv,next1,lw1,j,lv,
     $     nl,irtc,lw,iii,lenw,isp0,itfuplevel,itfdownlevel,
     $     ivj
      real*8 v,getva,va,vx
      character*(*) word
      character*128 word1
      logical*4 exist,get,rel,maxf,minf,
     $     abbrev,var,exist1,diff,vcomp, cont,vs,nvs

c     Initialize to avoid compiler warning
      lv=itfuplevel()
      if(iflinep .eq. 0)then
        call tfinitlinep(irtc)
      endif
      vs=.false.
      nvs=.false.
      v=0
      kv=-1
      iv=-1
      id=0
      isp0=isp
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
        ii=nelvx(i)%klp
        if(get)then
          id=idtypec(ii)
          get=.false.
          iv=nelvx(i)%ival
          kv=0
          var=.true.
          call peekwd(word1,next1)
          lw1=lenw(word1)
          if(lw1 .gt. MAXPNAME)then
            go to 912
          endif
          call capita1(word1(1:lw1))
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
            ipoint=next
            exist=.true.
            if(cont)then
              go to 9000
            else
              call termes(lfno,'?Missing value for ',word)
              exist=.false.
              exit
            endif
          else
            kv=itftypekey(id,word1,lw1,kx)
            if(kv .eq. 0)then
              go to 912
            elseif(tfreallistq(kx,kl))then
              iv=-1
            else
              iv=kv
            endif
          endif
          next=next1
          ipoint=next
          exist=.true.
 912      if(iv .eq. 0)then
            call termes(lfno,'?No default keyword for ',word)
            exist=.false.
            exit
          endif
          v=getva(exist1)
          if(.not. exist1)then
            if(cont)then
            else
              call termes(lfno,'?Missing value for ',word)
            endif
            exist=.false.
            exit
          endif
        endif
        if(idtypec(ii) .ne. id)then
          id=idtypec(ii)
          kv=itftypekey(id,word1,lw1,kx)
          if(kv .gt. 0)then
            ivi=kv
          elseif(tfreallistq(kx,kl))then
            ivi=-1
          else
            ivi=nelvx(i)%ival
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
          ivi=nelvx(i)%ival
        endif
        if(ivi .eq. 0)then
          call termes(lfno,'?No default keyword for ',word)
          exist=.false.
          exit
        endif
        var=ivi .eq. nelvx(i)%ival
        if(rel)then
          call loc_comp(idvalc(ii),cmpd)
          va=cmpd%value(ivi)*(1.d0+v)
        else
          va=v
        endif
        call compelc(ii,cmpd)
        cmpd%update=cmpd%nparam .le. 0
        if(var)then
          vx=cmpd%value(ivi)/errk(1,ii)
          if(minf)then
            if(maxf)then
              nelvx(i)%vlim(1)=-abs(va)
              nelvx(i)%vlim(2)=abs(va)
              vx=min(abs(va),max(-abs(va),vx))
            else
              nelvx(i)%vlim(1)=va
              vx=max(va,vx)
            endif
          elseif(maxf)then
            nelvx(i)%vlim(2)=va
            vx=min(va,vx)
          else
            vx=va
          endif
c          call tfsetcmp(vx*errk(1,ii),cmpd,ivi)
          cmpd%value(ivi)=vx*errk(1,ii)
c          rlist(latt(ii)+ivi)=vx*errk(1,ii)
        else
          vx=va
          if(ivi .gt. 0)then
            if(ktfnonrealq(cmpd%dvalue(ivi)))then
              call tfree(ktfaddr(cmpd%dvalue(ivi)))
            endif
            cmpd%value(ivi)=vx
            if(.not. vcomp)then
              call tftouch(i,ivi)
              isp=isp+1
              itastk(1,isp)=i
              itastk(2,isp)=ivi
            endif
          else
            do j=1,kl%nl
              ivj=int(kl%rbody(j))
              if(ktfnonrealq(cmpd%dvalue(ivj)))then
                call tfree(ktfaddr(cmpd%dvalue(ivj)))
              endif
              cmpd%value(ivj)=vx
              if(.not. vcomp)then
c                call tftouch(i,ivi)
                isp=isp+1
                itastk(1,isp)=i
                itastk(2,isp)=ivj
              endif
            enddo
          endif
        endif
        if(.not. vcomp)then
          isp=isp+1
          itastk(1,isp)=i
          itastk(2,isp)=ivi
          vs=vs .or. var
          nvs=nvs .or. .not. var
        endif
      enddo LOOP_II
      if(exist)then
        cont=.true.
        go to 1
      endif
 9000 if(isp .gt. isp0)then
        call tffsadjust1(isp0,vs,nvs)
      endif
      isp=isp0
      lv=itfdownlevel()
      return
      end
