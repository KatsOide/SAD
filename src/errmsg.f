      subroutine errmsg(rtn,msg,ercode,erlvl)
      implicit none
      integer*4 ercode,erlvl
      character*(*) rtn,msg
      include 'inc/MACFILE.inc'
      include 'inc/MACTTYP.inc'
      include 'inc/MACMISC.inc'
c
      integer*4 lene,lst,ln
      character*(MAXLLEN) cwork
      character*(*) astr
      parameter (astr=' ??? ')
c     
c      print *,"errmsg",rtn(:8)
c
      if(msglvl .le. erlvl) then
         ln=lene(rtn)
         lst=lene(msg)
         lst=min(lst,80-2*len(astr)-ln)
         cwork=astr//rtn(:ln)//astr//msg(:lst)
         write(errfl,*)cwork(:ln+lst+2*len(astr))
      endif
c
      if (erlvl .lt. 16) then
         return
      else if (erlvl .lt. 32) then
         write(errfl,*) '16<= error level < 32'
         call myfflush
         call toplvl
      else 
         write(errfl,*) 'error level >=32.'
         stop 9000
      endif
      return
      end
