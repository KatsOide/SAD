      subroutine tfmat(latt,mult,twiss,gammab,
     1                 lfno,word,wordp,exist)
      use tfstk
      use ffs
      use tffitcode
      implicit none
      integer*4 lfno,lenw,idp,l1,ielm,l2,i,j,lene
      integer*4 latt(2,nlat),mult(nlat),
     $     twiss(nlat,-ndim:ndim,ntwissfun)
      real*8 gammab(nlat),r
      real*8 trans(4,5)
      character*(*) word,wordp
      character*12 autofg,tt(5),title
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
      l1=ielm(latt,wordp,1,mult,exist)
      if(.not. exist)then
        l1=1
        l2=nlat
        go to 10
      else
        call getwdl2(word,wordp)
        l2=ielm(latt,wordp,1,mult,exist)
        if(.not. exist)then
          l2=nlat
        endif
      endif
10    call tftmat(twiss,gammab,trans,l1,l2,idp,.true.)
      if(symp)then
        r=1.d0
      else
        r=sqrt(gammab(l1)/gammab(l2))
      endif
      call elname(latt,l1,mult,name1)
      call elname(latt,l2,mult,name2)
      write(lfno,*)title(1:lene(title))//' transfer matrix from '//
     1             name1(1:lenw(name1)),' to ',name2(1:lenw(name2))
      do 20 i=1,4
        do 30 j=1,5
          tt(j)=autofg(r*trans(i,j),'11.6')
30      continue
        write(lfno,9001)(tt(j),j=1,5)
9001    format(1x,5a)
20    continue
      return
      end
