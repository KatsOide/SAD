      subroutine tcorr(word,latt,pos,master,lfno)
      use tfstk
      use ffs
      use tffitcode
      use ffs_pointer, only:idelc,idvalc,idtypec
      use tfcsi,only:cssetp
      implicit real*8 (a-h,o-z)
      integer*8 latt(nlat),le,le1,ix,ktaloc,i1,i2
      dimension pos(nlat)
      integer*4 master(nlat)
      character*(*) word
      character*20 name
      logical exist,exist1,abb,temat,abbrev,add
      add=.false.
      call peekwd(word,next)
      if(abbrev(word,'A_DD','_'))then
        call cssetp(next)
        add=.true.
      endif
      call tfgetr(cl,1.d0,'DELCOR cor. length',lfno,exist)
      if(.not. exist)then
        go to 9990
      endif
      call tfgetr(dx,1.d0,'DELCOR delx',lfno,exist)
      if(.not. exist)then
        go to 9990
      endif
      call tfgetr(dy,1.d0,'DELCOR dely',lfno,exist)
      if(.not. exist)then
        go to 9990
      endif
      al=pos(nlat)
      if(.not. cell)then
        al=al*2.d0
      endif
      nc=128
      do while(nc .le. 32768 .and. nc*.5d0 .lt. al)
        nc=nc*2
      enddo
      ix=ktaloc(nc*2)
      call tvcorr(rlist(ix),cl,al,nc)
      exist1=.false.
 1    call peekwd(word,next)
      if(word .eq. ' ')then
        go to 9980
      endif
      abb=ifany(word,'*%{}',1) .gt. 0
      exist=.false.
      do 10 i=1,nlat-1
        id=idtypec(i)
        if(id .eq. 1 .or.
     1     (id .gt. 8 .and. id .ne. 20 .and.
     $       id .ne. 31 .and. id .ne. 41))then
          go to 10
        endif
        if(temat(i,name,word))then
          exist=.true.
          call cssetp(next)
          if(id .eq. 41)then
            ie=i
          else
            ie=master(i)
            if(ie .le. 0)then
              go to 10
            endif
          endif
          i0=i
          s=(pos(i0)+pos(ie+1))*.5d0
          ns=mod(int(nc*s/al)+nc/2,nc)+1
          ns1=mod(ns,nc)+1
          w=mod(nc*s/al+nc/2-(ns-1),1.d0)
          i1=ix+(ns-1)*2
          i2=ix+(ns1-1)*2
          delx=(rlist(i1)*(1.d0-w)+rlist(i2)*w)*dx
          dely=(rlist(i1+1)*(1.d0-w)+rlist(i2+1)*w)*dy
          if(id .eq. 20)then
            istep=ie-i0
          else
            istep=1
          endif
          do 110 j=i0,ie,istep
            if(master(j) .ne. 0 .or. id .eq. 41)then
              le=latt(j)+5
              le1=le+1
              dx1=delx
              dy1=dely
              if(id .eq. 2)then
                le=le+4
                le1=le+1
              elseif(id .eq. 20)then
                le=le-2
                le1=le+1
                if(.not. add)then
                  dx1=rlist(idvalc(j)+3)-delx
                  dy1=rlist(idvalc(j)+4)-dely
                else
                  dx1=-delx
                  dy1=-dely
                endif
              elseif(id .eq. 31)then
                le=le-2
                le1=le+1
                if(.not. add)then
                  dx1=rlist(idvalc(j)+13)-delx
                  dy1=rlist(idvalc(j)+14)-dely
                else
                  dx1=-delx
                  dy1=-dely
                endif
              elseif(id .eq. 41)then
                le=latt(j)+15
                le1=latt(j)+17
              endif
              if(add)then
                rlist(le)=rlist(le)+dx1
                rlist(le1)=rlist(le1)+dy1
              else
                rlist(le)=dx1
                rlist(le1)=dy1
              endif
              if(id .eq. 41 .and. i .eq. 1)then
                dxi=rlist(le)
                dyi=rlist(le1)
              endif
            endif
110       continue
          if(.not. abb)then
            go to 11
          endif
        endif
10    continue
11    call tfadjst(latt,pos)
      exist1=exist1 .or. exist
      if(exist)then
        go to 1
      else
        if(abb)then
          call cssetp(next)
          go to 1
        else
          if(.not. exist1)then
            call termes(lfno,'?Undefined element ',word)
          endif
          go to 9980
        endif
      endif
9980  call tfree(ix)
      return
9990  call termes(lfno,
     1'Syntax: DELCOR corlength delx dely element [element1...]',' ')
      return
      end
