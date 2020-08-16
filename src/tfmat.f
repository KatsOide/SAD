      subroutine tfmat(lfno,word,wordp,exist)
      use tfstk
      use ffs
      use ffs_pointer
      use tffitcode
      implicit none
      integer*4 lfno,lenw,idp,l1,ielme,l2,i,j,lene
      real*8 r
      real*8 trans(4,5),trans6(6,6)
      character*(*) word,wordp
      character*12 autofg,tt(6),title
      character*(MAXPNAME+16) name1,name2
      logical exist,symp,abbrev
      symp=.true.
      title='Symplectic'
      idp=0
      call getwdl2(word,wordp)
      if(abbrev(word,'S_YMPLECTIC','_'))then
        symp=.true.
        title='Symplectic'
        call getwdl2(word,wordp)
      elseif(abbrev(word,'P_HYSICAL','_'))then
        symp=.false.
        title='Physical'
        call getwdl2(word,wordp)
      endif
      l1=ielme(wordp,exist,lfno)
      if(.not. exist)then
        l1=1
        l2=nlat
      else
        call getwdl2(word,wordp)
        l2=ielme(wordp,exist,lfno)
        if(.not. exist)then
          l2=nlat
        endif
      endif
      r=merge(1.d0,sqrt(gammab(l1)/gammab(l2)),symp)
      call elname(l1,name1)
      call elname(l2,name2)
      write(lfno,*)title(1:lene(title))//' transfer matrix from '//
     1             name1(1:lenw(name1)),' to ',name2(1:lenw(name2))
      if(calc6d)then
        call tftmat6(trans6,l1,l2)
        do i=1,6
          do j=1,6
            tt(j)=autofg(r*trans6(i,j),'11.6')
          enddo
          write(lfno,9001)(tt(j),j=1,6)
        enddo
      else
        call tftmat(trans,l1,l2,idp,.true.)
        do i=1,4
          do j=1,5
            tt(j)=autofg(r*trans(i,j),'11.6')
          enddo
          write(lfno,9001)(tt(j),j=1,5)
        enddo
      endif
      return
 9001 format(1x,6a)
      end

      subroutine tftmat6(trans,l1,l2)
      use ffs
      use ffs_pointer, only:twiss
      use temw,only:etwiss2ri,tinv6
      implicit none
      integer*4 l1,l2
      real*8 trans(6,6),tw1(ntwissfun),tw2(ntwissfun),
     $     ri1(6,6),v(6),
     $     dnx,dny,dnz,cx,sx,cy,sy,cz,sz
      logical*4 normal1,normal2
      tw1=twiss(l1,0,1:ntwissfun)
      tw2=twiss(l2,0,1:ntwissfun)
      ri1=etwiss2ri(tw1,normal1)
      dnx=tw2(mfitnx)-tw1(mfitnx)
      dny=tw2(mfitny)-tw1(mfitny)
      dnz=tw2(mfitnz)-tw1(mfitnz)
      cx=cos(dnx)
      sx=sin(dnx)
      cy=cos(dny)
      sy=sin(dny)
      cz=cos(dnz)
      sz=sin(dnz)
      trans(1,:)= cx*ri1(1,:)+sx*ri1(2,:)
      trans(2,:)=-sx*ri1(1,:)+cx*ri1(2,:)
      trans(3,:)= cy*ri1(3,:)+sy*ri1(4,:)
      trans(4,:)=-sy*ri1(3,:)+cy*ri1(4,:)
      trans(5,:)= cz*ri1(5,:)+sz*ri1(6,:)
      trans(6,:)=-sz*ri1(5,:)+cz*ri1(6,:)
      trans=matmul(tinv6(etwiss2ri(tw2,normal2)),trans)
c      call tinv6(ri2,ri2i)
c      call tmultr(trans,ri2i,6)
      if(.not. normal2)then
        normal1=.not. normal1
      endif
      if(.not. normal1)then
        v=trans(3,:)
        trans(3,:)=trans(1,:)
        trans(1,:)=v
        v=trans(4,:)
        trans(4,:)=trans(2,:)
        trans(2,:)=v
      endif
      return
      end

