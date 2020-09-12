      module mackw
        use maccode
        integer*4 , parameter ::
     $     kwL=1,kwANGL=kwL+1,kwROT=kwANGL+1,kwK0=kwROT+1,
     $     kwK1=kwK0+1,kwDK1=kwK1+1,kwK2=kwDK1+1,kwDK2=kwK2+1,
     $     kwK3=kwDK2+1,kwK4=kwK3+1,kwK5=kwK4+1,
     $     kwK6=kwK5+1,kwK7=kwK6+1,kwK8=kwK7+1,
     $     kwK9=kwK8+1,kwK10=kwK9+1,kwDK3=kwK10+1,kwA3=kwDK3+1,
     $     kwA5=kwA3+1,kwA7=kwA5+1,kwA9=kwA7+1,
     $     kwA11=kwA9+1,kwA13=kwA11+1,kwA15=kwA13+1,
     $     kwA17=kwA15+1,kwE1=kwA17+1,kwE2=kwE1+1,
     $     kwTILT=kwE2+1,kwKICK=kwTILT+1,kwDX=kwKICK+1,kwDY=kwDX+1,
     $     kwVOLT=kwDY+1,kwPHI=kwVOLT+1,kwHARM=kwPHI+1,
     $     kwFREQ=kwHARM+1,kwAX=kwFREQ+1,kwAY=kwAX+1,
     $     kwBX=kwAY+1,kwBY=kwBX+1,kwEX=kwBY+1,kwEY=kwEX+1,
     $     kwDP=kwEY+1,kwRAD=kwDP+1,kwCHRO=kwRAD+1,kwDZ=kwCHRO+1,
     $     kwEPX=kwDZ+1,kwEPY=kwEPX+1,kwPX=kwEPY+1,kwPY=kwPX+1,
     $     kwSIGZ=kwPY+1,kwSIGE=kwSIGZ+1 ,kwAZ=kwSIGE+1,
     $     kwDPX=kwAZ+1,kwDPY=kwDPX+1,
     $     kwR1=kwDPY+1,kwR2=kwR1+1,kwR3=kwR2+1,kwR4=kwR3+1,
     $     kwDETR=kwR4+1,kwEMIX=kwDETR+1,kwEMIY=kwEMIX+1,
     $     kwINDX=kwEMIY+1,kwBMAX=kwINDX+1,kwPRD=kwBMAX+1,
     $     kwRHO0=kwBMAX+1,kwRHO=kwRHO0+1,kwCHI1=kwRHO+1,
     $     kwCHI2=kwCHI1+1,kwCHI3=kwCHI2+1,kwGEO=kwCHI3+1,
     $     kwBZ=kwGEO+1,kwBND=kwBZ+1,kwFRIN=kwBND+1,kwK0FR=kwFRIN+1,
     $     kwEPS=kwK0FR+1,kwF1=kwEPS+1,kwF2=kwF1+1,kwFRMD=kwF2+1,
     $     kwRANK=kwFRMD+1,kwRANP=kwRANK+1,kwRANV=kwRANP+1,
     $     kwKIN=kwRANV+1,kwLWAK=kwKIN+1,kwTWAK=kwLWAK+1,
     $     kwSK1=kwTWAK+1,kwSK2=kwSK1+1,
     $     kwSK3=kwSK2+1,kwSK4=kwSK3+1,kwSK5=kwSK4+1,kwSK6=kwSK5+1,
     $     kwSK7=kwSK6+1,kwSK8=kwSK7+1,kwSK9=kwSK8+1,kwSK10=kwSK9+1,
     $     kwSLI=kwSK10+1,kwNP=kwSLI+1,kwDIR=kwNP+1,
     $     kwLAG=kwDIR+1,kwSK0=kwLAG+1,
     $     kwJDX=kwSK0+1,kwJDPX=kwJDX+1,kwJDY=kwJDPX+1,kwJDPY=kwJDY+1,
     $     kwK11=kwJDPY+1,kwK12=kwK11+1,kwK13=kwK12+1,
     $     kwK14=kwK13+1,kwSK11=kwK14+1,kwSK12=kwSK11+1,kwSK13=kwSK12+1,
     $     kwSK14=kwSK13+1,kwV1=kwSK14+1,
     $     kwZX=kwV1+1,kwZPX=kwZX+1,kwZY=kwZPX+1,kwZPY=kwZY+1,
     $     kwXANGLE=kwZPY+1,kwPZ=kwXANGLE+1,kwEMIZ=kwPZ+1,
     $     kwSTURN=kwEMIZ+1,kwBSTRL=kwSTURN+1,
     $     kwDX1=kwBSTRL+1,kwDX2=kwDX1+1,kwDY1=kwDX2+1,kwDY2=kwDY1+1,
     $     kwD11=kwDY2+1,kwD12=kwD11+1,kwD13=kwD11+2,
     $     kwD14=kwD11+3,kwD15=kwD11+4,kwD16=kwD11+5,
     $     kwD21=kwD11+6,kwD22=kwD21+1,kwD23=kwD21+2,
     $     kwD24=kwD21+3,kwD25=kwD21+4,kwD26=kwD21+5,
     $     kwD31=kwD21+6,kwD32=kwD31+1,kwD33=kwD31+2,
     $     kwD34=kwD31+3,kwD35=kwD31+4,kwD36=kwD31+5,
     $     kwD41=kwD31+6,kwD42=kwD41+1,kwD43=kwD41+2,
     $     kwD44=kwD41+3,kwD45=kwD41+4,kwD46=kwD41+5,
     $     kwD51=kwD41+6,kwD52=kwD51+1,kwD53=kwD51+2,
     $     kwD54=kwD51+3,kwD55=kwD51+4,kwD56=kwD51+5,
     $     kwD61=kwD51+6,kwD62=kwD61+1,kwD63=kwD61+2,
     $     kwD64=kwD61+3,kwD65=kwD61+4,kwD66=kwD61+5,
     $     kwB11=kwD61+6,kwB12=kwB11+1,kwB13=kwB11+2,
     $     kwB14=kwB11+3,kwB15=kwB11+4,kwB16=kwB11+5,
     $     kwB22=kwB11+6,kwB23=kwB11+7,kwB24=kwB11+8,
     $     kwB25=kwB11+9,kwB26=kwB11+10,kwB33=kwB11+11,kwB34=kwB11+12,
     $     kwB35=kwB11+13,kwB36=kwB11+14,kwB44=kwB11+15,kwB45=kwB11+16,
     $     kwB46=kwB11+17,kwB55=kwB11+18,kwB56=kwB11+19,kwB66=kwB11+20,
     $     kwR11=kwB11+21,kwR12=kwR11+1,kwR13=kwR11+2,
     $     kwR14=kwR11+3,kwR15=kwR11+4,kwR16=kwR11+5,
     $     kwR22=kwR11+6,kwR23=kwR11+7,kwR24=kwR11+8,
     $     kwR25=kwR11+9,kwR26=kwR11+10,kwR33=kwR11+11,kwR34=kwR11+12,
     $     kwR35=kwR11+13,kwR36=kwR11+14,kwR44=kwR11+15,kwR45=kwR11+16,
     $     kwR46=kwR11+17,kwR55=kwR11+18,kwR56=kwR11+19,kwR66=kwR11+20,
     $     kwKx=kwR11+21,kwQy=kwkx+1,kwpole=kwQy+1,
     $     kwFBx=kwpole+1,kwFBy=kwFBx+1,kwDBZ=kwFBy+1,
     $     kwK15=kwDBZ+1,kwK16=kwK15+1,kwK17=kwK16+1,
     $     kwK18=kwK17+1,kwK19=kwK18+1,kwK20=kwK19+1,
     $     kwK21=kwK20+1,kwSK15=kwK21+1,kwSK16=kwSK15+1,kwSK17=kwSK16+1,
     $     kwSK18=kwSK17+1,kwSK19=kwSK18+1,kwSK20=kwSK19+1,
     $     kwSK21=kwSK20+1,kwV20=kwSK21+1,kwV11=kwV20+1,kwV02=kwV11+1,
     $     kwJDZ=kwV02+1,kwJDPZ=kwJDZ+1,kwOFFSET=kwJDPZ+1,
     $     kwCOUPLE=kwOFFSET+1,kwDDP=kwCOUPLE+1,kwRADI=kwDDP+1,
     $     kwDPHI=kwRADI+1,kwW1=kwDPHI+1,kwDROT=kwW1+1,
     $     kwAE1=kwDROT+1,kwAE2=kwAE1+1,kwFB1=kwAE2+1,kwFB2=kwFB1+1,
     $     kwLDEV=kwFB2+1,kwLRAD=kwLDEV+1,kwFL=kwLRAD+1,kwAPHI=kwFL+1,
     $     kwF1K1F=kwAPHI+1, kwF2K1F=kwF1K1F+1,
     $     kwF1K1B=kwF2K1F+1,kwF2K1B=kwF1K1B+1,
     $     kwDVOLT=kwF2K1B+1,kwPROF=kwDVOLT+1,
     $       kwNPARAM=kwPROF+1,kwMAX=kwNPARAM+1

      integer*4 , allocatable :: kytbl(:,:)
      integer*4, pointer :: kyindex(:,:),kyindex1(:,:)
      integer*4 INDMAX

      contains
        subroutine initkyindex
        use tfcode
        use maccbk,only:pname
        implicit none
        integer*4 i,k,id
        INDMAX=0
        do i=icDRFT,icMXEL
          INDMAX=max(INDMAX,kytbl(kwMAX,i))
        enddo
        allocate(kyindex(0:INDMAX,0:icMXEL))
        allocate(kyindex1(0:INDMAX,0:icMXEL))
        kyindex=0
        kyindex1=0
        do i=icDRFT,icMXEL
          do k=1,kwMAX-2
            id=kytbl(k,i)
            if(id .ne. 0)then
              if(kyindex(id,i) .eq. 0)then
                kyindex(id,i)=k
              elseif(kyindex1(id,i) .eq. 0)then
                kyindex1(id,i)=k
              else
                write(*,*)'Too many aliases ',pname(kytbl(0,i)),
     $               ' ',pname(kytbl(k,0))
                call abort
              endif
            endif
          enddo
        enddo
        return
        end subroutine

      end module

      module macttyp
        integer*4, parameter :: ttypNM=4, ttypID=2, ttypEF=-1,
     $       ttypDM=1,ttypST=8,ttypLS=16
      end module

      module macvar
        integer*4, parameter ::
     $     VarPt=1, VarLog=2, VarInt=4, VarRl=8, VarStr=128, VarLst=256,
     $     maxv=20
      end module

      module macfile
        integer*4, parameter :: STDERR=6,STDIN=5,STDOUT=6,STDPLT=8,
     $     STDLST=21, MAXLLEN=256
c
        integer*8 inflpt
        integer*4 errfl,infl,outfl,pltfl,msglvl,lstfl
        integer*4 pbuf
        character buf(MAXLLEN)
      end module

      module macmisc
        real*8, parameter :: EPS = 1.0d-7
      character*1, parameter ::
     $       LCURL='{',RCURL='}',LPAR='(',RPAR=')',COMMA=',',
     $       SEMIC=';',EQCHR='=',MINUS='-',PLUS='+',STAR='*'
      end module

      module cbkmac
        use macfile
        use macmisc
        use macttyp
        use maccode
        use mackw
      end module

      module tfcsi
        use tfcbk, only:maxlbuf
        implicit none
        integer*4, parameter :: nbmax=maxlbuf,nsav=5,ipbase=1
        type csiparam
          sequence
          integer*4 isav(1:0)
          integer*4 lfni,lfn1,lfno,ios
          logical*4 rep
        end type
        type (csiparam) , target :: savep
        character*16 delim,cmnt
        integer*8 ibcloc
        integer*4, pointer:: ipoint,lrecl,lfni=>savep%lfni,
     $       lfn1=>savep%lfn1,lfno=>savep%lfno,ios=>savep%ios
        logical*4 , pointer :: rep=>savep%rep
        integer*4 iconv,ldel,lcmnt,lastln,ibegt,lastt
        character*(nbmax) , target  :: buffer0
        character*(nbmax) , pointer :: buffer

        contains
        subroutine cssetp(ip)
        implicit none
        integer*4 ip
        ipoint=ip
        return
        end subroutine

        subroutine cssets(ip)
        implicit none
        integer*4 ip
        ios=ip
        return
        end subroutine

        subroutine cssetl(ip)
        implicit none
        integer*4 ip
        lrecl=ip
        return
        end subroutine

        subroutine cssetlfni(ip)
        implicit none
        integer*4 ip
        lfni=ip
        return
        end subroutine

        subroutine cssetlfno(ip)
        implicit none
        integer*4 ip
        lfno=ip
        return
        end subroutine

        subroutine cssetlfn1(ip)
        implicit none
        integer*4 ip
        lfn1=ip
        return
        end subroutine

        integer*4 function icsmrk()
        implicit none
        icsmrk=ipoint
        return
        end function

        integer*4 function icsstat()
        implicit none
        icsstat=ios
        return
        end function

        integer*4 function icslrecl()
        implicit none
        icslrecl=lrecl
        return
        end function

        integer*4 function icslfni()
        implicit none
        icslfni=lfni
        return
        end function

        integer*4 function icslfno()
        implicit none
        icslfno=lfno
        return
        end function

        integer*4 function icslfn1()
        implicit none
        icslfn1=lfn1
        return
        end function

        integer*4 function igetrecl() result(nc)
        implicit none
        integer*4 i
        nc=max(lrecl-ipoint,0)
        if(nc .gt. 0)then
          do i=ipoint,ipoint+nc-1
            if(buffer(i:i) .eq. char(10))then
              nc=i-ipoint
              return
            endif
          enddo
        endif
        return
        end function

      end module

      module tfrbuf
      use tfstk
      use maccbk, only:i00
      implicit none
      ;integer*4 ,parameter :: irbnofile=-999,irbunread=-1,
     $     irbeof=-99
      integer*4 , parameter ::
     $     modeclose=0,moderead=1,modewrite=2,modestring=3,
     $     modeshared=4,modemapped=5
      integer*4 , parameter :: nbuf=1024
      integer*4 ncprolog
      character*128 prolog
      integer*4 :: ifd(0:nbuf)=0
      integer*4 , target :: lbuf(0:nbuf)=0,mbuf(0:nbuf)=0,
     $     itbuf(0:nbuf)=0,lenbuf(0:nbuf)=0
      integer*8 :: ibuf(0:nbuf)=0
      type (sad_descriptor) :: ntable(nbuf)
      type cbkshared
        integer*8, allocatable :: ca(:)
      end type
      type (cbkshared) rbshared(0:nbuf)
      data ntable(:)%k /nbuf*0/

      contains
      integer*4 function nextfn(mode)
      use macfile, only:MAXLLEN
      implicit none
      integer*4 mode
      integer*4 ios,f,is
      logical*4 od
      character*(MAXLLEN) msg
c
      f=0
      is=11
      do f=is,nbuf
        if(itbuf(f) .eq. modeclose)then
          inquire(f,IOSTAT=ios,err=9000,OPENED=od)
          if( .not. od) then
            itbuf(f)=mode
            ibuf(f)=0
            nextfn=f
            return
          endif
        endif
      end do
c
 1000 continue
      nextfn=0
      return
c
 9000 continue
      call perror(msg)
      go to 1000
c
      end function

      subroutine trbassign(lfn)
      use tfstk
      use tfcsi
      use iso_c_binding
      implicit none
      integer*4 lfn
      if(ibuf(lfn) .ne. 0)then
        call c_f_pointer(c_loc(jlist(1,ibuf(lfn))),buffer)
      else
        buffer=>buffer0
        ibuf(lfn)=transfer(c_loc(buffer0),int8(0))/8
      endif
      ipoint=>mbuf(lfn)
      lrecl=>lbuf(lfn)
      lfni=lfn
      return
      end subroutine

      subroutine trbreset(lfn)
      implicit none
      integer*4 lfn
      if(itbuf(lfn) .le. modewrite)then
        lbuf(lfn)=0
        mbuf(lfn)=1
      else
        mbuf(lfn)=lbuf(lfn)+1
      endif
      return
      end subroutine

      subroutine trbclose(lfn)
      use tfstk
      use tfshare
      implicit none
      integer*4 lfn,irtc
      select case (itbuf(lfn))
      case (modewrite)
        close(lfn)
        if(ibuf(lfn) .gt. 0)then
          if(ilist(2,ibuf(lfn)-1) .ne. 0)then
            call unixclose(ilist(2,ibuf(lfn)-1))
            ilist(2,ibuf(lfn)-1)=0
          endif
          call tfree(ibuf(lfn))
        endif
      case(modeclose,moderead)
        close(lfn)
      case (modestring)
        call tflocal1(ibuf(lfn)-1)
        ibuf(lfn)=0
      case (modeshared)
        call tfreeshared(ibuf(lfn),irtc)
        ibuf(lfn)=0
      case default
        call unmapfile(klist(ibuf(lfn)),int8(lenbuf(lfn)))
        call unixclose(ifd(lfn))
      end select
      lbuf(lfn)=0
      mbuf(lfn)=1
      itbuf(lfn)=modeclose
      return
      end subroutine

      subroutine trbnextl(lfn)
      implicit none
      integer*4 lfn
      if(lfn .gt. 0)then
        mbuf(lfn)=lbuf(lfn)+1
      endif
      return
      end subroutine

      subroutine trbeor2bor(lfn)
      implicit none
      integer*4 lfn
      if(lfn .gt. 0)then
        if(mbuf(lfn) .eq. lbuf(lfn))then
          mbuf(lfn)=lbuf(lfn)+1
        endif
      endif
      return
      end subroutine

      integer*8 function itrbibuf(lfn,mode) result(ia)
      implicit none
      integer*4 , intent(in) :: lfn,mode
      if(lfn .gt. 0 .and. itbuf(lfn) .eq. mode)then
        ia=ibuf(lfn)
      else
        ia=0
      endif
      return
      end function

      subroutine trbmovepoint(lfn,nc)
      implicit none
      integer*4 , intent(in) :: lfn,nc
      if(lfn .gt. 0)then
        mbuf(lfn)=min(lbuf(lfn)+1,max(1,mbuf(lfn)+nc))
      endif
      return
      end subroutine

      subroutine trbinit(lfn,ib)
      use tfstk
      implicit none
      integer*4 , intent(in) ::lfn,ib
      itbuf(lfn)=ib
      if(ib .le. modewrite)then
        if(ibuf(lfn) .eq. 0)then
          ibuf(lfn)=ktaloc(maxlbuf/8)
          lenbuf(lfn)=maxlbuf
        endif
        ilist(2,ibuf(lfn)-1)=0
      endif
      lbuf(lfn)=0
      mbuf(lfn)=1
      return
      end

      end module

      module ffsfile
        integer*4 , parameter :: maxlfn=128
        integer*4 :: lfnp=0,lfnbase=0
        integer*4 lfnstk(maxlfn)
      end module

      module trackbypass
        implicit none
        logical*4 :: bypasstrack=.false.
        integer*8 :: lattuse=0, lattredef=0
      end module

      module readbuf

      contains
      subroutine tfreadbuf(lfn,ib,nc)
      use tfstk
      use tfshare
      use tfcsi
      use tfrbuf
c      use iso_c_binding
      implicit none
      integer*8 ls1,mapresizefile,lenfile,ib1
      integer*4 ,intent(in):: lfn,ib
      integer*4 nc
      integer*4 itfgetbuf,irtc,ls,ie,i
      if(lfn .le. 0 .or. ibuf(lfn) .eq. 0)then
        go to 9000
      endif
      if(itbuf(lfn) .le. modewrite)then
c        write(*,*)'tfreadbuf ',lfn,itbuf(lfn),modewrite
        nc=itfgetbuf(lfn,jlist(ib,ibuf(lfn)),
     $       maxlbuf-ib-256,irtc)
        if(irtc .ne. 0)then
          return
        endif
        lbuf(lfn)=ib-1
        lenbuf(lfn)=lbuf(lfn)+nc
c        write(*,*)'tfreadbuf ',lfn,nc,lbuf(lfn),mbuf(lfn),ib
        mbuf(lfn)=ib
      else
 11     ls=lenbuf(lfn)
c        if(lbuf(lfn) .lt. ls .and.
c     $       jlist(lbuf(lfn)+1,ibuf(lfn)) .eq. 10)then
c          lbuf(lfn)=lbuf(lfn)+1
c        endif
        if(lbuf(lfn) .lt. ls)then
          ie=ls
          do i=lbuf(lfn)+1,ls
            if(jlist(i,ibuf(lfn)) .eq. 10)then
              if(i .eq. 1 .or. jlist(i-1,ibuf(lfn)) .ne.
     $             ichar('\\'))then
                ie=i
                exit
              endif
            endif
          enddo
          nc=ie-lbuf(lfn)
          mbuf(lfn)=lbuf(lfn)+1
          lbuf(lfn)=ie
        else
          if(itbuf(lfn) .ge. modemapped)then
            ls1=lenfile(ifd(lfn))
            if(ls .lt. ls1)then
              ib1=mapresizefile(klist(ibuf(lfn)),ifd(lfn),
     $             ls,ls1)/8
              if(ib1 .lt. 0)then
                go to 9000
              endif
              mbuf(lfn)=ls
              ibuf(lfn)=ib1
              lenbuf(lfn)=int(ls1)
              go to 11
            endif
          elseif(mbuf(lfn) .lt. lbuf(lfn))then
            nc=1
            mbuf(lfn)=lbuf(lfn)
            return
          endif
          mbuf(lfn)=ls+1
          go to 9000
        endif
      endif
      return
 9000 nc=irbeof
      return
      end subroutine

      subroutine irbopen1(j,ib,is,jfd)
      use tfrbuf
      use tfstk
      implicit none
      integer*4 ,intent(in):: j,jfd
      integer*8 ,intent(in):: ib,is
      if(itbuf(j) .eq. 0)then
        itbuf(j)=int(is)
        lenbuf(j)=0
        ifd(J)=0
        select case (int(is))
        case (moderead,modewrite)
        case (modestring)
          ibuf(j)=ktfcopy1(ib)+1
          lenbuf(j)=ilist(1,ib)
        case (modeshared)
          ibuf(j)=ib
          lenbuf(j)=ilist(1,ib-1)
        case default
          ibuf(j)=ib
          lenbuf(j)=int(is-modemapped)
          ifd(j)=jfd
        end select
        lbuf(j)=0
        mbuf(j)=1
c        write(*,*)'irbopen1 ',j,ibuf(j)
        return
      endif
      return
      end

      subroutine readstr(in,str,irtc)
      use tfrbuf
      use tfcsi, only:buffer
      implicit none
      integer*4 ,intent(in):: in
      integer*4 ,intent(out):: irtc
      integer*4 nc
      character*(*) str
c      write(*,*)'reststr ',in
      call tfreadbuf(in,lbuf(in)+1,nc)
c      write(*,*)': ',nc,'''',buffer(mbuf(in):mbuf(in)+nc-1),''''
      irtc=0
      if(nc .gt. 0)then
        str=buffer(mbuf(in):mbuf(in)+nc-1)
        mbuf(in)=mbuf(in)+nc
        if(str(nc:nc) .eq. char(10))then
          str(nc:)=' '
        else
          str(nc+1:)=' '
        endif
      else
        str=' '
        if(nc .lt. 0)then
          irtc=-1
        endif
      endif
      return
      end

      function tfopenshared(isp1,irtc) result(kx)
      use tfstk
      use tfrbuf
      use tfshare
      implicit none
      type (sad_descriptor) kx
      integer*4 ,intent(in):: isp1
      integer*4 ,intent(out):: irtc
      integer*8 ia
      integer*4 itfmessage,n,m,iu,nc
      if(isp .ne. isp1+1)then
        kx=dxnullo
        irtc=itfmessage(9,'General::narg','"1"')
        return
      elseif(ktfnonrealq(ktastk(isp)))then
        kx=dxnullo
        irtc=itfmessage(9,'General::wrongtype','"Real"')
        return
      endif
      n=int(rtastk(isp))
      if(n .lt. 0)then
        irtc=itfmessage(9,'General::wrongval','"Non-negative"')
        kx%k=kxfailed
        irtc=0
        return
      endif
      m=(n+7)/8+2
      irtc=0
      ia=ktfallocshared(m+1)
c      ia=mapallocfixed8(rlist(0), m+1, 8, irtc)
      if(irtc .ne. 0)then
        irtc=itfmessage(9,'General::mmap','""')
        kx%k=kxfailed
        irtc=0
        return
      endif
      call trbopen(iu,ia,int8(modeshared),nc)
      if(iu .le. 0)then
        kx%k=kxfailed
        call tfreeshared(ia)
        irtc=itfmessage(9,'General::fileopen','"(Shared)"')
        return
      endif
      ilist(1,ia)=m
      ilist(2,ia)=0
      klist(ia+3)=ktfoper+mtfnull
      kx=dfromr(dble(iu))
      return
      end

      function tfreadshared(isp1,irtc) result(kx)
      use tfstk
      use tfrbuf
      implicit none
      type (sad_descriptor) kx
      integer*4 ,intent(in):: isp1
      integer*4 ,intent(out):: irtc
      type (sad_string), pointer :: str
      integer*8 ia
      integer*4 itfmessage,isp0,iu,ist
      logical*4 tfcheckelement
      if(isp .ne. isp1+1)then
        kx=dxnullo
        irtc=itfmessage(9,'General::narg','"1"')
        return
      elseif(ktfnonrealq(ktastk(isp),iu))then
        kx=dxnullo
        irtc=itfmessage(9,'General::wrongtype','"Real"')
        return
      endif
      irtc=0
      ia=itrbibuf(iu,modeshared)
      if(ia .eq. 0)then
        kx%k=kxeof
        return
      endif
      ist=ilist(2,ia)
      do while(ist .ne. 0)
        if(ist .ne. 1 .and. ist .ne. -1)then
          write(*,*)'tfreadshared shared memory destructed ',ist,ia
          call abort
        endif
        call tpause(10000)
        ist=ilist(2,ia)
      enddo
      ilist(2,ia)=1
      kx=dlist(ia+1)
c      write(*,*)'readshared ',ia,kx%k
      if(ktfobjq(kx))then
c        write(*,*)'readshared-obj '
        if(ktfsymbolq(kx))then
c          write(*,*)'readshared-symbol '
          if( .not. tfcheckelement(kx,.false.))then
            irtc=itfmessage(99,'General::wrongtype',
     $           '"undefined symbol returned in Shared"')
            kx%k=ktfoper+mtfnull
            ilist(2,ia)=0
            return
          endif
        elseif(ktfstringq(kx,str))then
c          call tfdebugprint(kx,'readshared-string',1)
c          write(*,*)'at ',ia
          kx=kxsalocb(-1,str%str(1:str%nch),str%nch)
        elseif(ktflistq(kx))then
          isp0=isp
          call tfrecallshared(isp0,ktflist+ia+3,kx,irtc)
          isp=isp0
          if(irtc .ne. 0)then
            kx%k=ktfoper+mtfnull
          endif
        else
c          write(*,*)'readshared-other '
          kx%k=ktfoper+mtfnull
        endif
      endif
      ilist(2,ia)=0
      return
      end

      function tfwriteshared(isp1,irtc) result(kx)
      use tfstk
      use tfrbuf
      implicit none
      type (sad_descriptor) kx
      integer*4 ,intent(in):: isp1
      integer*4 ,intent(out):: irtc
      integer*8 kas,ka,kt,kap,k
      integer*4 itfmessage,itfmessageexp,isp0,n,i,iu,ist
      kx=dxnullo
      if(isp .ne. isp1+2)then
        irtc=itfmessage(9,'General::narg','"2"')
        return
      elseif(ktfnonrealq(ktastk(isp1+1),iu))then
        irtc=itfmessage(9,'General::wrongtype','"Real"')
        return
      endif
      kas=itrbibuf(iu,modeshared)
      if(kas .eq. 0)then
        irtc=itfmessage(99,'Shared::notopen','""')
        return
      endif
      ist=ilist(2,kas)
      do while(ist .ne. 0)
        if(ist .ne. 1 .and. ist .ne. -1)then
          write(*,*)'writeshared shared memory destructed: ',ist,kas
          call abort
        endif
        call tpause(10000)
        ist=ilist(2,kas)
      enddo
      ilist(2,kas)=-1
      k=ktastk(isp)
      if(ktfnonobjq(k) .or. ktfsymbolq(k))then
        klist(kas+1)=k
      else
        ka=ktfaddr(k)
        kt=k-ka
        if(kt .eq. ktfstring)then
          if(ilist(1,ka) .gt. ilist(1,kas)*8)then
            ilist(2,kas)=0
            irtc=itfmessageexp(9,'Shared::toolarge',
     $         dble(ilist(1,ka)-ilist(1,kas)*8))
            return
          endif
          call tmov(ilist(1,ka+1),ilist(1,kas+3),
     $         (ilist(1,ka)+7)/8)
          ilist(1,kas+2)=ilist(1,ka)
          klist(kas+1)=ktfstring+kas+2
c          write(*,*)'writeshared-string ',kas,klist(kas+1)
        elseif(kt .eq. ktflist)then
          isp0=isp
          call tfsharedsize(isp0,k,n,irtc)
          if(irtc .ne. 0)then
            ilist(2,kas)=0
            isp=isp1+2
            return
          endif
          if(n .gt. ilist(1,kas))then
            ilist(2,kas)=0
            irtc=itfmessageexp(9,'Shared::toolarge',
     $           dble((n-ilist(1,kas))*8))
            return
          endif
          kap=kas+2
          do i=isp1+3,isp
            n=itastk2(1,i)
            ktastk2(i)=kap
            kap=kap+n
          enddo
          do i=isp1+3,isp
            call tfstoreshared(isp1+3,ktastk(i),ktastk2(i))
          enddo
          isp=isp1+2
          klist(kas+1)=ktflist+kas+3
        else
          ilist(2,kas)=0
          irtc=itfmessage(9,'General::wrongtype',
     $         '"No Patterns"')
          return
        endif
      endif
      ilist(2,kas)=0
      irtc=0
      return
      end

      subroutine trbopen(lfn,ib,is,ifd)
      use tfrbuf,only:itbuf,modeclose,nbuf
      implicit none
      integer*8, intent(in) :: ib,is
      integer*4, intent(in) :: ifd
      integer*4, intent(out):: lfn
      integer*4 j
      do j=nbuf,11,-1
        if(itbuf(j) .eq. modeclose)then
          lfn=j
          call irbopen1(lfn,ib,is,ifd)
          return
        endif
      enddo
      lfn=0
      return
      end subroutine

      subroutine trbopenmap(str,kx,irtc)
      use tfstk
      use tfrbuf,only:modemapped
      implicit none
      character*(*) , intent(in)::str
      integer*4 , intent(out)::irtc
      type (sad_descriptor) , intent(out)::kx
      integer*8 kfile,ksize,mapallocfile
      integer*4 lfn,ifd
      kfile=mapallocfile(str,ifd,ksize,irtc)
      if(irtc .eq. 0)then
        call trbopen(lfn,kfile/8,ksize+modemapped,ifd)
        kx%x(1)=dble(lfn)
      else
        kx%x(1)=-1.d0
      endif
      return
      end subroutine

      end module readbuf

      subroutine toplvl
      use trackbypass, only: bypasstrack
      use tfrbuf
      use maccbk
      use mackw
      use macttyp
      use macfile
      use macmisc
      use gammaf, only:aginit
      use tfstk, only:tfinitstk
      use ffsfile, only:lfnbase
      implicit none
      integer*8 argp
      integer*4 jslen,jttype,jival
      integer*4 slen,ttype,ival,hsrch,hsrchz,status
      integer*4 index,IgetGL
      real*8  jrval,rval
      character*(MAXSTR) jtoken,token
c
c      call omp_set_dynamic(.true.)
c      call omp_set_num_threads(1)
      call aginit
 1000 continue
      if (IgetGL('$CTIME$') .eq. FLAGON) call cputix
      call gettok(token,slen,ttype,rval,ival)
c     for debug
c       print *,'toplvl-0 ',token(:slen),slen,infl,ttype,ttypEF
c     end debug
c
 1100 continue
      if (ttype .eq. ttypEF) go to 9000
      if(ttype .eq. ttypID) then
c..........System defined name.
         index=hsrch(token(:slen))
c         write(*,*)'toplvl-1 ',index,idtype(index),idval(index),icRSVD
         if (idtype(index) .eq. icDEF) then
           call tfinitstk
            if (idval(index) .lt. icMXEL) then
              call doelem(idval(index))
            else if((idval(index) .ge. icLINE)
     &              .and. (idval(index) .lt. icRSVD))  then
              call doline(idval(index))
            else
              call errmsg('toplvl','system bug: check initbl',0,0)
            endif
            go to 1000
         else if (idtype(index) .eq. icACT) then
            call tfinitstk
            if(bypasstrack)then
              return
            endif
            call doACT(index)
            go to 1000
         else if (idtype(index) .eq. icRSVD) then
           call tfinitstk
           call funcall1(idval(index),argp)
           go to 1000
         endif
c..........User defined variable.
c     print *,token(:slen)
         index=hsrchz(token(:slen))
         if ((idtype(index) .lt. icMXEL) .and.
     &        (idtype(index) .ne. icNULL)) then
           call tfinitstk
           call doelm2(idtype(index),pname(index),slen,ttype)
         else if ((idtype(index) .eq. icLINE)) then
           call tfinitstk
           call dolin2(index,slen,ttype)
         else
            call tfinitstk
           call dAssgn(token,slen,status)
           if (status .ne. 0) call errmsg('toplvl',
     &          'Unsupported function '//token(:slen)//'!',
     &          0,16)
c     for debug
c      print *,token(:slen),
c     &     index,pname(index),idval(index),idtype(index)
c     end debug
         endif
      else if(token(:slen) .eq. ';') then
         go to 1000
      else if(token(:slen) .eq. RCURL) then
         go to 1000
       elseif(token(:slen) .eq. char(13))then
         go to 1000
      else
         call errmsg('main'
     &        ,'syntax error: invalid input for toplvl '//
     &        '"'//token(:slen)//'"'
     &        ,0,0)
      endif
      go to 1000
c
 9000 continue
      if(itbuf(infl) .ne. 0 .or. lfnbase .gt. 1)then
        return
      endif
c      print *," SAD1 reads EOF."
      call errmsg("main","Stop execution.(READ EOF)" ,0,0)
c
      stop
c.....Entry point for error handling.
      entry bigjmp(jtoken,jslen,jttype,jrval,jival)
c........for debug
c     print *,'Big Jump !!!'
c
      token=jtoken
      slen=jslen
      ttype=jttype
      rval=jrval
      ival=jival
      go to 1100
      end

