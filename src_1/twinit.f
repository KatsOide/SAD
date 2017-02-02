      subroutine twinit(waked,wakes,lw,nwx,
     1                  nb,ns,sb,ds)
      implicit none
      integer*4 lw,nwx,nb,ns,i,j,iu,il,im,n,nsi,nsb
      real*8 waked(2,lw),wakes(nwx+ns*nb),sb,ds,s,si,dws
      il=1
      do 10 i=1,nb
        if(i .eq. 1)then
          si=0.d0
          nsi=ns
          nsb=0
        else
          si=(i-1)*sb-(ns-1)*ds
          nsi=ns*2-1
          nsb=(i-2)*nsi+ns
        endif
        do 20 j=1,nsi
          s=si+(j-1)*ds
          n=nsb+j
          iu=lw
          if(s .gt. waked(1,lw))then
            il=iu
            wakes(n)=0.d0
          elseif(s .lt. waked(1,1))then
            wakes(n)=0.d0
          else
211         if(iu-il .gt. 1)then
              im=il+(iu-il)/2
              if(s .ge. waked(1,im))then
                il=im
              else
                iu=im
              endif
              go to 211
            endif
            dws=waked(1,iu)-waked(1,il)
            wakes(n)=((s-waked(1,il))*waked(2,iu)+
     1                (waked(1,iu)-s)*waked(2,il))/dws
          endif
20      continue
10    continue
      wakes(1)=wakes(1)*.5d0
      return
      end
