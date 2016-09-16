      Subroutine prSad(idx)
      use maccbk
      use maccode
      use macfile
      implicit real*8 (a-h,o-z)
      integer idx
c
      integer*4 len,oldfl
      integer*8 idxl,plist,mkplst,i
c for debug
c     call ptrace('prSad',1)
c end debug
      oldfl=outfl
      outfl=21
      call sprlin(idx)
      idxl=idval(idx)
      plist=mkplst(idxl)
      len=ilist(1,plist)-1
c     print *,'prsad',plist,len
      do 1100 i=plist+1,plist+len
        ival=idtype(ilist(2,i))
        if(ival .eq. icLINE) then
           call sprlin(ilist(2,i))
        else if (ival .lt. icMXEL) then
           call prelem(ilist(2,i),' ')
        else
          call errmsg('prSad','Invalid element'//pname(ilist(2,i))
     &               ,0,0)
        endif
 1100 continue
      call tfree(plist)
c      call tfreem(int8(plist),ilist(1,plist))
c      call freeme(plist,ilist(1,plist))
      do 2110 i=1,HTMAX
        if((idtype(i) .ge. icGLI) .and.
     &       (idtype(i) .le. icGLR)) then
           call iprglb(i)
        else if(idtype(i) .eq. icFLAG) then
           call prnflg(pname(i))
        endif
 2110 continue
c     write(outfl,*)'! End of print '//pname(idx)//' sad'
      outfl=oldfl
c for debug
c     call ptrace('prSad',1)
c end debug
      return
      end
