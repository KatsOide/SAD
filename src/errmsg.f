      subroutine errmsg(rtn,msg,ercode,erlvl)
      use macttyp
      use macfile
      use macmisc
      implicit none
      integer*4 ercode,erlvl
      character*(*) rtn,msg
c
      integer*4 lene,lst,ln
      character*(MAXLLEN) cwork
      character*(5) ,parameter ::astr=' *** '
c     
c      print *,"errmsg",rtn(:8)
c
      if(msglvl .le. erlvl) then
         ln=lene(rtn)
         lst=lene(msg)
         lst=min(lst,80-2*len(astr)-ln)
         cwork=astr//rtn(:ln)//astr//msg(:lst)
         write(errfl,*)cwork(:ln+lst+2*len(astr)),' ercode = ',ercode
      endif
c
      if (erlvl .lt. 16) then
         return
      else if (erlvl .lt. 32) then
         write(errfl,*) '16<= error level < 32'
         call exit(8000)
c         call myfflush
c         call toplvl
      else 
         write(errfl,*) 'error level >=32.'
         call exit(9000)
      endif
      return
      end
