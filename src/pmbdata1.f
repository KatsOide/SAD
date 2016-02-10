      subroutine pmbdata1(latt,mult,datas,iobs,nobs,observ,ndata,yhi,
     $     print,lfno)
      use ffs
      use tffitcode
c----- Write table of beam size. -------      
      implicit real*8 (a-h,o-z)
      parameter (item=20,itemn=item/2)
      logical print,single
      character*(*) observ
      character line*131,lepton(2)*4,autofg*8,name*8
      dimension latt(*),mult(*)
      dimension iobs(nobs+1),datas(itemn,2,nobs,*)
      data lepton/'(e+)','(e-)'/
c
      if(print) then
c....... register ylo in the data buffer.
        datas(itemn,1,1,ndata)=sign(datas(itemn,1,1,ndata),-1d0)
      else
        datas(itemn,2,1,ndata)=yhi
        return
      endif
      single=observ.eq.'EMIY'.or. observ.eq.'SIGY'
      write(lfno,'(12x,9a8)')'   sigx ','  ex*dp ','   sigy ',
     $     '  ey*dp ','phi/phi0','    dx ','    dy ','  sminx ',
     $     '  sminy '
      do 12 i=1,nobs
        j=iobs(i+1)
        call elname(latt,j,mult,name)
        line(1:8)=name
        do 11 l=1,2
          if(single) then
            if(charge*(dble(l)-1.5d0).gt.0) goto 11
          endif
          line(9:12)=lepton(l)
          do 10 k=1,itemn-1
            if((k.eq.6.or.k.eq.7).and. .not.single) then
              if(l.eq.1) then
                line(13+(k-1)*8:)=autofg(datas(k,l,i,ndata)-
     $               datas(k,2,i,ndata),'8.5')
              else
                line(13+(k-1)*8:)='     :=0'
              endif
            else
              line(13+(k-1)*8:)=autofg(datas(k,l,i,ndata),'8.5')
            endif
 10       continue
          write(lfno,'(a)') line(1:lene(line))
          line=' '
 11     continue
 12   continue
      return
      end
