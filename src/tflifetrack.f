      subroutine tflifetrack(isp1,kx,irtc)
      use tfstk
      implicit none
      integer*8 kx
      integer*4 isp1,irtc,itfmessage,nip,i,lstr
      parameter (lstr=512)
      real*8 gw
      character*(lstr) file,comm
      if(isp .ne. isp1+13)then
        irtc=itfmessage(9,'General::narg','"13"')
        return
      endif
      if(ktfnonrealq(ktastk(isp1+10)))then
        go to 9000
      endif
      nip=int(rtastk(isp1+10))
      if(nip .le. 0)then
        irtc=itfmessage(9,'General::wrongval',
     $       '" # of IP must be positive"')
        return
      endif
      do i=1,9
        if(.not. tfreallistq(ktastk(isp1+i)))then
          go to 9000
        endif
      enddo
      if(ilist(2,ktfaddr(ktastk(isp1+1))-1) .ne. 3)then
        irtc=itfmessage(9,'General::wrongval',
     $       '"EmittanceStrong={ex,ey,ez}"')
        return
      endif
      if(ilist(2,ktfaddr(ktastk(isp1+2))-1) .ne. 36*nip)then
        irtc=itfmessage(9,'General::wrongval',
     $       '"Dim[BeamMatrixStrong]={#IP,6,6}"')
        return
      endif
      if(ilist(2,ktfaddr(ktastk(isp1+3))-1) .ne. 3)then
        irtc=itfmessage(9,'General::wrongval',
     $       '"EmittanceWeak={ex,ey,ez}"')
        return
      endif
      if(ilist(2,ktfaddr(ktastk(isp1+4))-1) .ne. 36*nip)then
        irtc=itfmessage(9,'General::wrongval',
     $       '"Dim[BeamMatrixWeak]={#IP,6,6}"')
        return
      endif
      if(ilist(2,ktfaddr(ktastk(isp1+5))-1) .ne. 6*nip)then
        irtc=itfmessage(9,'General::wrongval',
     $       '"Dim[OrbitWeak]={#IP,6}"')
        return
      endif
      if(ilist(2,ktfaddr(ktastk(isp1+6))-1) .ne. 6*nip)then
        irtc=itfmessage(9,'General::wrongval',
     $       '"Dim[Separations]={#IP,6}"')
        return
      endif
      if(ilist(2,ktfaddr(ktastk(isp1+7))-1) .ne. 3)then
        irtc=itfmessage(9,'General::wrongval',
     $       '"DampingRate={dx,dy,dz}"')
        return
      endif
      if(ilist(2,ktfaddr(ktastk(isp1+8))-1) .ne. 3)then
        irtc=itfmessage(9,'General::wrongval',
     $       '"Apertures={ax,ay,az}"')
        return
      endif
      if(ilist(2,ktfaddr(ktastk(isp1+9))-1) .ne. nip)then
        irtc=itfmessage(9,'General::wrongval',
     $       '"Dim[Particles]={#IP}"')
        return
      endif
      if(ktfnonrealq(ktastk(isp1+11)))then
        irtc=itfmessage(9,'General::wrongval',
     $       '"GammaWeak=Real"')
        return
      endif
      if(ktfnonstringq(ktastk(isp1+12)))then
        irtc=itfmessage(9,'General::wrongval',
     $       '"OutputFile=String"')
        return
      endif
      if(ktfnonstringq(ktastk(isp1+13)))then
        irtc=itfmessage(9,'General::wrongval',
     $       '"Comment=String"')
        return
      endif
      gw=rtastk(isp1+11)
      file=""
      call tmovb(rlist(ktfaddr(ktastk(isp1+12))+1),file,
     $     min(ilist(1,ktfaddr(ktastk(isp1+12))),lstr))
      comm=""
      call tmovb(rlist(ktfaddr(ktastk(isp1+12))+1),comm,
     $     min(ilist(1,ktfaddr(ktastk(isp1+12))),lstr))
      call sad_lifetrack()
c$$$     $     rlist(ktfaddr(ktastk(isp1+1))),
c$$$     $     rlist(ktfaddr(ktastk(isp1+2))),
c$$$     $     rlist(ktfaddr(ktastk(isp1+3))),
c$$$     $     rlist(ktfaddr(ktastk(isp1+4))),
c$$$     $     rlist(ktfaddr(ktastk(isp1+5))),
c$$$     $     rlist(ktfaddr(ktastk(isp1+6))),
c$$$     $     rlist(ktfaddr(ktastk(isp1+7))),
c$$$     $     rlist(ktfaddr(ktastk(isp1+8))),
c$$$     $     rlist(ktfaddr(ktastk(isp1+9))),
c$$$     $     nip,gw,file,comm)
      irtc=0
      kx=ktfoper+mtfnull
      return
 9000 irtc=itfmessage(9,'General::wrongtype',
     $       '"9*NumList,nip,GammaWeak,Out_String,Comment_String"')
      return
      end

      subroutine sad_lifetrack
      return
      end
