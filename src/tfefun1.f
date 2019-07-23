      recursive subroutine tfefun1(isp1,id0,kx,ref,irtc)
      use tfstk
      use tflinepcom
      use temw, only:tfnormalcoord,tfinitemip
      implicit none
      type (sad_descriptor) kx,kxj,kispi
      type (sad_dlist), pointer :: kli
      type (sad_rlist), pointer :: klx
      integer*8 ka
      integer*4 isp1,id,irtc,na,i,ispi,m,j,narg,id0,
     $     itfmessage,l
      real*8 rgetgl1
      logical*4 ref
      irtc=-1
      id=id0-1000
      narg=isp-isp1
      ka=klist(ifunbase+ktfaddr(ktastk(isp1)))+1
      na=ilist(1,ka)
      do i=1,min(na,isp-isp1)
        ispi=isp1+i
        if(ilist(1,ka+i) .ne. 0)then
          if(tflistq(ktastk(ispi),kli))then
            kispi=dtastk(ispi)
            m=kli%nl
            kx=kxaaloc(-1,m,klx)
            irtc=0
            do j=1,m
              dtastk(ispi)=kli%dbody(j)
              call tfefun1(isp1,id0,kxj,ref,irtc)
              if(irtc .ne. 0)then
                do l=j,m
                  klx%rbody(l)=0.d0
                enddo
                return
              endif
              call tfsetlist(dtfcopy(kxj),kx,j)
            enddo
            dtastk(ispi)=kispi
            return
          endif
        endif
      enddo
c-------next 4 lines were modified by Kikuchi---------------
      go to (5010,5020,5030,5040,5050,5060,5070,5080,5090,5100,
     $       5110,5120,5130,5140,5150,5160,5170,5180,5190,5200,
     $       5210,5220,5230,5240,5250,5260,5270,5280,5290,5300,
     $       5310,5320,5330,5340,5350,5360,5370,5380,5390,5400,
     $       5410,5420,5430,5440,5450,5460,5470
     $     ),id
c            ELEM TWIS LINE CalE TraP CalO DyAp RspM Mast FLAG
c            Mcad exDA InDA MAP  FFS  RadF RadS Flag ExBL BLNm
c            SetE TKey NCD6 Tcl  FFSH ExpB GetM ClTr LiTr RGBC
c            CPro TcA1 CaSy TkOA CaSD TcSR LifT SyBE CSRI CSRM
c            CSRC CSRT CSRH CSOS CSOM AliP BBBR
      write(*,*)
     $'Wrong implementation of function (tfefun1).  ID = ',id0
      return
 5010 call tfelement(isp1,kx,ref,irtc)
      return
 5020 call tftwiss(isp1,kx,ref,irtc)
      return
 5030 call tfline(isp1,kx,ref,irtc)
      return
 5040 call tfemit(isp1,kx,irtc)
      return
 5050 call tftrack(isp1,kx,irtc)
      return
 5060 call tfoptics(isp1,kx,irtc)
      return
 5070 call tfdapert(isp1,kx,irtc)
      return
c-------Kikuchi added-----     
 5080 call pgrmat(isp1,kx,irtc)
      return
 5090 call pgmast(isp1,kx,irtc)
      return
 5100 call pgflag(isp1,kx,irtc)
      return
 5110 call pgsolvcond(isp1,kx,irtc)
      return
c-------Kikuchi addition end-----     
 5120 call execdapackage(isp1,kx,irtc)
      return
 5130 call initdainterface(isp1,kx,irtc)
      return
 5140 call getdamap(isp1,kx,irtc)
      return
 5150 call tfffs(isp1,kx,irtc)
      return
 5160 call tlfield(isp1,kx,irtc)
      return
 5170 call tlspect(isp1,kx,int(rgetgl1('NPARA')),irtc)
      return
 5180 call tfflags(isp1,kx,irtc)
      return
 5190 call tfextractbeamline(isp1,kx,irtc)
      return
 5200 call tfbeamlinename(isp1,kx,irtc)
      return
 5210 call tfsetelement(isp1,kx,irtc)
      return
 5220 call tftypekey(isp1,kx,irtc)
      return
 5230 call tfnormalcoord(isp1,kx,irtc)
      return
 5240 call tfinitemip
      irtc=0
      return
 5250 if(isp .gt. isp1+1)then
        irtc=itfmessage(99,'General::narg','"0"')
      else
        call tshow(kx,irtc,.true.,.false.)
      endif
      return
 5260 call tfexpandbeamline(isp1,kx,irtc)
      return
 5270 call tfmain(isp1,kx,irtc)
      return
 5280 irtc=itfmessage(999,'General::unregister',' ')
      return
 5290 irtc=itfmessage(999,'General::unregister',' ')
      return
 5300 call tfrgbcolor(isp1,kx,irtc)
      return
 5310 irtc=itfmessage(999,'General::unregister',' ')
      return
 5320 irtc=itfmessage(999,'General::unregister',' ')
      return
 5330 irtc=itfmessage(999,'General::unregister',' ')
      return
 5340 irtc=itfmessage(999,'General::unregister',' ')
      return
 5350 irtc=itfmessage(999,'General::unregister',' ')
      return
 5360 irtc=itfmessage(999,'General::unregister',' ')
      return
 5370 call tflifetrack(isp1,kx,irtc)
      return
 5380 call tfemits(isp1,kx,irtc)
      return
 5390 call tfcsrinit(isp1,kx,irtc)
      return
 5400 call tfcsrmat(isp1,kx,irtc)
      return
 5410 call tfcsrvec(isp1,kx,irtc)
      return
 5420 call tfcsrtrack(isp1,kx,irtc)
      return
 5430 call tfcsrhaissin(isp1,kx,irtc)
      return
 5440 call tfcsroysetup(isp1,kx,irtc)
      return
 5450 call tfcsroymat(isp1,kx,irtc)
      return
 5460 call tfsurvivedparticles(isp1,kx,irtc)
      return
 5470 call tbcube1(isp1,kx,irtc)
      return
      end
