      Subroutine prnflg(fname)
      use maccbk
      implicit real*8 (a-h,o-z)
      character*(*) fname
      include 'inc/MACCODE.inc'
      include 'inc/MACFILE.inc'
      integer*4 fval,idummy,slen,lenw
      character*(MAXPNAME) cwork
      slen=lenw(fname)
      if(fname(1:1) .ne. '$') then
        cwork='$'//fname(:slen)
      else
        cwork=fname(:slen)
      endif
      if(fname(slen:slen) .ne. '$') then
        cwork=cwork(1:lenw(cwork))//'$'
      endif
      fval=IgetGL(cwork,idummy)
      if (fval .eq. FLAGOF)then
        write(outfl,*)' OFF ', cwork(2:lenw(cwork)-1), ';'
      else
        write(outfl,*)' ON  ', cwork(2:lenw(cwork)-1), ';'
      endif
      return
      end
