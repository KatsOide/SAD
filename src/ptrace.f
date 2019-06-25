      subroutine ptrace(name,ind)
      use maccbk
      implicit real*8 (a-h,o-z)
      character*(*) name
      integer ind
      integer indent
      common /ptrac/indent
      character*80 msg,prmpt*1
c
      if (IgetGL('$DEBUG$',idummy) .EQ. FLAGOF) return
      if(indent .lt. 1) indent=1
      if (ind .gt. 0) then
        indent=indent+2*ind
        msg(indent:)='enter into '//name
        prmpt='>'
      else
        msg(indent:)='leave from '//name
        prmpt='<'
      endif
      do 100 i=1,indent-1
        msg(i:i)=prmpt
  100 continue
      if (ind .lt. 0 ) indent=indent+2*ind
      print *,msg
      return
      end
