      subroutine prkick1(latt,mult,istr,nstr,temp,push,
     1                   print,lfno)
      use tfstk
      use ffs
      use tffitcode
      implicit real*8 (a-h,o-z)
      logical push,print
cHP   character*8 vout(2)*80,autofg,name(2)
      character*(8) vout(2)*79,autofg*8,name(2)
      dimension latt(2,nlat),mult(*),
     1          istr(nstra,4),temp(nstr),rc(6)
      data rc /6*0d0/
      data vout(1)/' kick new (mrd)'/,
     1     vout(2)/'      totl(mrd)'/
      data name/2*'        '/
      include 'inc/common.inc'
      save vout,rc,name
c
      if(push) then
        vout(1)(48:55)=autofg(rc(1)*1d3,'S8.5')
        vout(1)(57:64)=name(1)
        vout(1)(64:71)=autofg(rc(2)*1d3,'S8.5')
        vout(1)(72:79)=autofg(rc(3)*1d3,'S8.5')
        vout(2)(48:55)=autofg(rc(4)*1d3,'S8.5')
        vout(2)(57:64)=name(2)
        vout(2)(64:71)=autofg(rc(5)*1d3,'S8.5')
        vout(2)(72:79)=autofg(rc(6)*1d3,'S8.5')
      endif
      if(nstr.eq.0) return
      do 10 i=1,nstr
        j=istr(i,2)
   10   temp(i)=rlist(latt(2,istr(j,1))+11)-rlist(ibckup-1+istr(j,1))
      call mstatp(temp,nstr,ceil,floor,rc(2),rc(3),imax)
      rc(1)=max(ceil,-floor)*sign(1d0,ceil+floor)
      call elname(latt,istr(istr(imax,2),1),mult,name(1))
      do 12 i=1,nstr
   12   temp(i)=rlist(latt(2,istr(istr(i,2),1))+11)
      call mstatp(temp,nstr,ceil,floor,rc(5),rc(6),imax)
      rc(4)=max(ceil,-floor)*sign(1d0,ceil+floor)
      call elname(latt,istr(istr(imax,2),1),mult,name(2))
      if(print) then
        vout(1)(16:23)=autofg(rc(1)*1d3,'S8.5')
        vout(1)(25:32)=name(1)
        vout(1)(32:39)=autofg(rc(2)*1d3,'S8.5')
        vout(1)(40:47)=autofg(rc(3)*1d3,'S8.5')
        vout(2)(16:23)=autofg(rc(4)*1d3,'S8.5')
        vout(2)(25:32)=name(2)
        vout(2)(32:39)=autofg(rc(5)*1d3,'S8.5')
        vout(2)(40:47)=autofg(rc(6)*1d3,'S8.5')
        write(lfno,'(t16,a,t24,a,t32,a,t40,a,t48,a,t56,a,t64,a,t72,a)')
     1  ' max    ',' @      ',' rms    ',' ave    ',
     1  ' max_old',' @_old  ',' rms_old',' ave_old'
        write(lfno,'(2(a:/))') (vout(i),i=1,2)
      endif
      return
      end
