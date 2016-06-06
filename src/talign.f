      subroutine talign(latt,word,wordp,pos,mult,lfno,exist)
      use tfstk
      use ffs
      use tffitcode
      implicit real*8 (a-h,o-z)
      dimension latt(2,nlat),mult(nlat),pos(nlat)
      character*(*) word,wordp
      logical exist
      i1=igelm(word,exist)
      if(.not. exist)then
        write(lfno,*)' Missing element for ALIGN.'
        return
      endif
      do 210 i=i1,nlat-1
        id=idtype(latt(1,i))
        if(id .eq. 4 .or. id .eq. 6 .or. id .eq. 8)then
          i1=i
          go to 211
        endif
210   continue
      write(lfno,*)' No quad or sext after ',word
      return
211   dx1=getva(exist)
      if(.not. exist)then
        call getwdl2(word,wordp)
        if(word .eq. '*')then
          dx1=rlist(latt(2,i1)+5)
        else
          dy1=0.d0
          go to 212
        endif
      endif
      dy1=getva(exist)
      call getwdl2(word,wordp)
      if(.not. exist)then
        if(word .eq. '*')then
          dy1=rlist(latt(2,i1)+6)
          call getwdl2(word,wordp)
        else
          dy1=0.d0
        endif
      endif
212   i2=ielm(wordp,exist)
      if(.not. exist)then
        i2=i1
        dx2=dx1
        dy2=dy1
      else
        word=' '
        dx2=getva(exist)
        if(.not. exist)then
          dx2=0.d0
          dy2=0.d0
        else
          dy2=getva(exist)
          if(.not. exist)then
            dy2=0.d0
          endif
        endif
      endif
      ds=pos(i2)-pos(i1)
      if(ds .eq. 0.d0)then
        r=0.d0
      else
        r=1.d0/ds
      endif
      ie=0
      do 10 i=i1,min(nlat-1,i2)
        if(i .le. ie)then
          go to 10
        endif
        k=idtype(latt(1,i))
        if(k .eq. 4 .or. k .eq. 6)then
          dx=((pos(i)+pos(i+1))*.5d0-pos(i1))*r*(dx2-dx1)+dx1
          dy=((pos(i)+pos(i+1))*.5d0-pos(i1))*r*(dy2-dy1)+dy1
          do 20 j=i1,nlat-1
            if(pname(latt(1,j)) .eq. pname(latt(1,i)))then
              ie=j
              rlist(latt(2,j)+5)=dx
              rlist(latt(2,j)+6)=dy
            elseif(pos(j) .ne. pos(j+1))then
              go to 10
            endif
20        continue
        endif
10    continue
      return
      end
