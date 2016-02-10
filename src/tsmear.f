      subroutine tsmear(n,nsmear,np,kptbl,
     1     sjx,sjxjx,sjy,sjyjy,sjz,sjzjz,lfno)
      include 'inc/TMACRO.inc'
      dimension sjx(np0),sjxjx(np0),sjy(np0),sjyjy(np0),
     1          sjz(np0),sjzjz(np0)
      dimension kptbl(np0,6)
      if(np .gt. 0)then
        write(lfno,9002)n
9002    format(' Smear at ',i6,'th turn:')
        write(lfno,'(a)')
     1' p. #   2Jx(m)      2Jy(m)      2Jz(m)    '//
     1      '  sigma2Jx  sigma2Jy  sigma2Jz'
      endif
      ajx=0.d0
      ajy=0.d0
      ajz=0.d0
      do i=1,np0
        ip=kptbl(i,1)
        if(ip .le. np)then
          sx=sjx(i)/nsmear
          if(sx .gt. 0.d0)then
            sxx=sqrt(abs(sjxjx(i)/nsmear-sx**2))/sx
          else
            sxx=0.d0
          endif
          sy=sjy(i)/nsmear
          if(sy .gt. 0.d0)then
            syy=sqrt(abs(sjyjy(i)/nsmear-sy**2))/sy
          else
            syy=0.d0
          endif
          sz=sjz(i)/nsmear
          if(sz .gt. 0.d0)then
            szz=sqrt(abs(sjzjz(i)/nsmear-sz**2))/sz
          else
            szz=0.d0
          endif
          ajx=ajx+sx
          ajy=ajy+sy
          ajz=ajz+sz
          write(lfno,9001)i,sx,sy,sz,sxx,syy,szz
9001      format(i6,1p,3g12.4,0p,3f10.6)
        endif
      enddo
      write(lfno,9003)ajx/np,ajy/np,ajz/np
 9003 format(' Avrg ',:,1p,3g12.4)
      return
      end

