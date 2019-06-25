      subroutine pvbump1(kfit,ifitp,mfitp,fitval,nfc,
     &                   kfitd,ifitd,fitvd)
      implicit real*8 (a-h,o-z)
      parameter (mfitc1=32,mfitc2=28)
      dimension kfit(nfc),ifitp(nfc),mfitp(nfc),fitval(nfc)
      dimension kfitd(*),ifitd(*),fitvd(*)
c
      nc=0
      do 10 i=1,nfc
        if(mfitp(i).ne.0) then
          if(kfit(i) .ge. mfitc1 .and. kfit(i) .le. mfitc1+3 .or.
     &       kfit(i) .ge. mfitc2 .and. kfit(i) .le. mfitc2+3  ) then
            nc=nc+1
            kfitd(nc)=kfit(i)
            ifitd(nc)=ifitp(i)
            fitvd(nc)=fitval(i)
          endif
        endif
 10   continue
      return
      end
