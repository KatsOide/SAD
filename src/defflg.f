      subroutine defflg(fname,fval)
      use maccbk
      use maccode
      implicit real*8 (a-h,o-z)
      integer*4 fval
      character*(*) fname
      integer idummy,lenw
      character*(MAXPNAME) cwork
      cwork='$'//fname(:lenw(fname))//'$'
      call defglb(cwork,icFLAG,idummy)
      call IsetGL(cwork,fval,idummy)
      return
      end
