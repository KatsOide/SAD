      integer*4 function lstchk(x,list,mpos,lsiz)
      integer*4 x,list(2,lsiz),lsiz,mpos
      do 1000 i=1,mpos
         if(x .eq. list(2,i)) then
           lstchk=i
           return
         endif
 1000 continue
      lstchk=0
      return
      end
