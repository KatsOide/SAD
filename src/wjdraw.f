      subroutine wjdraw(ajx,lsp,ndim,np,maxord)
      implicit none
      integer*4 lsp,ndim,i,j,lfno,np,maxord
      real*8 ajx(lsp,ndim),h,ajmax,ajmin,rgetgl1
      character*4 label(3)
      character*30 dat
      data label/'D2JX','D2JY','D2JZ'/
      lfno=rgetgl1('DJPLOT')
      if(lfno .le. 0)then
        return
      endif
      call fdate1(dat)
      write(lfno,9002)
9002  format(' NEWFRAME;SET FONT DUPLEX')
      write(lfno,*)'TITLE TOP ''DEVIATION OF J'''
      write(lfno,*)'CASE      '' LLLLLLLL LL  '''
      write(lfno,9003)'TITLE 10 9.7 SIZE 2 ''PARTICLE #',np,''''
9003  format(1x,a,i3,a)
      write(lfno,*)   'CASE                '' LLLLLLL  '''
      write(lfno,9004)'TITLE 10 9.3 SIZE 2 ''   ORDER  ',maxord,''''
9004  format(1x,a,i3,a)
      write(lfno,*)   'CASE                ''    LLLL  '''
      write(lfno,*)'TITLE 9 1 SIZE 1.75 ''',dat,''''
      h=7.d0/ndim
      do 10 i=1,ndim
        ajmax=-1.d30
        ajmin= 1.d30
        do 110 j=1,lsp
          ajmax=max(ajx(j,i),ajmax)
          ajmin=min(ajx(j,i),ajmin)
110     continue
        write(lfno,9001)9.d0-h*i,9.d0-h*(i-1),lsp,
     1                  sngl(ajmin),sngl(ajmax)
9001    format(' SET WINDOW X 3 12 Y ',2f10.6,/,
     1         ' SET LIMIT X 0 ',i5,' Y ',1p,2g15.7)
        if(i .eq. ndim)then
          write(lfno,*)'SET LABELS BOTTOM ON'
          write(lfno,*)'TITLE BOTTOM ''TURNS'''
        else
          write(lfno,*)'SET LABELS BOTTOM OFF'
        endif
        call tdinit(lfno,'JOIN 1','10')
        do 20 j=1,lsp
          call tdput(dble(j),ajx(j,i))
20      continue
        call tdterm
        write(lfno,*)'SET WINDOW X 2 12'
        write(lfno,*)'TITLE LEFT ''',label(i),''''
        write(lfno,*)'CASE       ''F  L'''
10    continue
      return
      end
