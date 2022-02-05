      Subroutine push(arg1,arg2,isp,istack,isize)
      integer arg1,arg2,isp,istack(2,isize),isize
c
      call ptrace('push',1)
c
      if(isp .ge. isize) then
        call errmsg('push','stack over flow',0,0)
        stop
      endif
        isp=isp+1
        istack(1,isp)=arg1
        istack(2,isp)=arg2
      call ptrace('push',-1)
        return
c
      Entry pop(arg1,arg2,isp,istack,isize)
      call ptrace('pop',1)
      if(isp .le. 0) then
        call errmsg('ipop','too many pop operation',0,0)
        stop
      endif
      arg1=istack(1,isp)
      arg2=istack(2,isp)
      isp=isp-1
      call ptrace('pop',-1)
        return
c
      end
