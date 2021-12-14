      subroutine ActGRA(argp)
      use maccbk
      use maccode
      use macvar
      implicit real*8 (a-h,o-z)
      integer*4 argp
      integer id,gtype
      integer slen
      character*(MAXSTR) option
c
      save
c
c.....for debug
c     call ptrace('ActGRA',1)
c.....end debug
      id=ilist(2,argp+1)
c     " pointer to graphic name
      gtype=ilist(2,argp+2)
c
      if (idtype(id) .ne. icGraf) then
        call NewGRF(id)
      endif
      slen=min(80,ilist(1,gtype))
      do 1100 i=1,(slen+7)/8
         write(option(i*8-7:i*8),'(A8)')rlist(ilist(2,gtype)+i-1)
 1100 continue
c.....for debug
c     print * ,'gtype = ',option(:slen)
c     print * ,'id = ',pname(id),id,idval(id)
c.....end debug
c.....for debug
c     call ptrace('ActGRA',-1)
c.....end debug
      return
      end
