      subroutine NewGRF(id)
      use maccbk
      use maccode
      use macvar
      use tfmem, only:ktaloc
      implicit none
      integer*4 id,LPARM,LBUF,ktcaloc
       parameter(LPARM=10,LBUF=1+2*1024)
c
       idtype(id)=icGRAF
       idval(id)=ktcaloc(LPARM)
c       idval(id)=mcfallo(LPARM)
       klist(idval(id))=ktaloc(LBUF)
c       ilist(2,idval(id))=mfalloc(LBUF)
       ilist(2,idval(id)-1)=LBUF
       ilist(1,klist(idval(id)))=LBUF
       ilist(1,idval(id)+1)=50+id
c      ilist(1,idval(id)+1)=50+gid
c      ilist(1,idval(id)+1)=fopen(gid)
       ilist(2,idval(id)+1)=ilist(2,idval(id))+1
       return
       end
