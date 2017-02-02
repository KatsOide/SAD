      subroutine twelm(lfno,latt,mult,idisp1,idisp2,title,irmgn,itab)
      use ffs
      use tffitcode
      integer*8 latt(nlat)
      dimension mult(nlat)
      character*(*) title
      character*16 name,name1,title1
      call elname(idisp1,name)
      if(idisp2 .gt. 0)then
        call elname(idisp2,name1)
      else
        name1=' '
      endif
      title1=title
      call twbuf(
     1    title1(1:lene(title1))//' '//name(1:lene(name))//' '//name1,
     1    lfno,1,irmgn,itab,1)
      return
      end
