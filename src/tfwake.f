      subroutine tfwake(word,lfno,err)
      use kyparam
      use tfstk
      use ffs
      use ffs_pointer
      use tffitcode
      integer*8 ktaloc,lwp
      integer*4 minsize,memmax
      parameter (minsize=2000,memmax=32768)
      character*(*) word
      character*16 name,word1
      logical*4 err,abbrev,temat,exist
      err=.false.
      call getwdl(word)
      if(abbrev(word,'L_ONGITUDINAL','_'))then
        lt=0
        ioff=ky_LWAK_CAVI
      elseif(abbrev(word,'T_RANSVERSE','_'))then
        lt=1
        ioff=ky_TWAK_CAVI
      else
        err=.true.
        return
      endif
      lwp=0
      call getwdl(word)
      do 10 i=1,nlat-1
        if(idtypec(i) .eq. icCAVI)then
          if(temat(i,name,word))then
            ip=ilist(1,latt(i)+ioff)
            if(ip .gt. 0)then
              do 20 j=1,nlat
                if(idtypec(i) .eq. 31)then
                  if(ilist(1,latt(j)+ioff) .eq. -ip)then
                    ilist(1,latt(j)+ioff)=ip
                    go to 21
                  endif
                endif
20            continue
              call tfree(int8(ilist(2,ip)))
              call tfree(int8(ip))
            endif
21          if(lwp .eq. 0)then
              na=max((memmax-2)/4,minsize)
              lwp=ktaloc(na*2+1)
              if(lwp .le. 0)then
                call termes(lfno,
     1          '?Insufficient memory space for WAKE.',' ')
                err=.true.
                return
              endif
              do 30 j=1,na
32              call peekwd(word1,next)
                if(word1 .eq. 'ENDWAKE')then
                  ipoint=next
                  nw=j-1
                  go to 31
                elseif(word1 .eq. ' ')then
                  call skiplnget
                  go to 32
                endif
                s=getva(exist)
                if(.not. exist)then
                  err=.true.
                  call termes(lfno,'?Missing ENDWAKE.',' ')
                  call tfree(lwp)
                  lwp=0
                  return
                endif
                w=getva(exist)
                if(.not. exist)then
                  err=.true.
                  call termes(lfno,'?Missing value of wake.',' ')
                  call tfree(lwp)
                  return
                endif
                rlist(lwp+j*2-1)=s
                rlist(lwp+j*2  )=w
30            continue
              err=.true.
              call termes(lfno,'?Too many wake data.',' ')
              call tfree(lwp)
              lwp=0
              return
31            continue
              if(na .gt. nw)then
                call tfreem(int8(lwp+nw*2+1),na*2-nw*2)
c                call freeme(lwp+nw*2+1,na*2-nw*2)
                ilist(1,lwp-1)=nw*2+2
              endif
              ilist(1,latt(i)+ioff)=lwp
              ilist(2,lwp)=0
            else
              ilist(1,latt(i)+ioff)=-lwp
            endif
          endif
        endif
10    continue
      return
      end
