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
            if(id /= 0)then
              if(kyindex(id,i) == 0)then
                kyindex(id,i)=k
              elseif(kyindex1(id,i) == 0)then
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
          integer*4 lfni,lfnm,lfno,lfne,ios
          logical*4 rep
        end type
        type (csiparam) , target :: savep
        character*16 delim,cmnt
        integer*8 ibcloc
        integer*4, pointer:: ipoint,lrecl,lfni=>savep%lfni,
     $       lfnm=>savep%lfnm,lfno=>savep%lfno,lfne=>savep%lfne,ios=>savep%ios
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

        subroutine cssetlfnm(ip)
        implicit none
        integer*4 ip
        lfnm=ip
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

        integer*4 function icslfnm()
        implicit none
        icslfnm=lfnm
        return
        end function

        integer*4 function igetrecl() result(nc)
        implicit none
        integer*4 i
        nc=max(lrecl-ipoint,0)
        if(nc > 0)then
          do i=ipoint,ipoint+nc-1
            if(buffer(i:i) == char(10))then
              nc=i-ipoint
              return
            endif
          enddo
        endif
        return
        end function

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

      module tmacro
        use mackw
        use macphys
        use macfile
        use tfstk
        real*8, parameter :: c=cveloc,hp=plankr,e=elemch,epsrad=1.d-6,
     $       emminv=1.d-15,eps00m=0.005d0,ampmaxm=0.05d0
        integer*4 ,parameter :: ndivmaxm=1000
        real*8 amass,charge,h0,p0,omega0,trf0,crad,erad,
     $       codin(6),dleng,anrad,urad,u0,vc0,wrfeff,dp0,brho,
     $       ccintr,cintrb,pbunch,coumin,re0,pgev,emidiv,
     $       emidib,emidiq,emidis,ctouck,dvemit,h1emit,
     $       anbunch,tdummy(6),zlost,alost,
     $       taurdx,taurdy,taurdz,fridiv,beamin(21),
     $       vcalpha,vceff,vcacc,dvcacc,ddvcacc,alphap,
     $       pspac_dx,pspac_dy,pspac_dz,dvfs,rcratio,rclassic,brhoz,
     $       bradprev,amom0,circ,hvc0,cuc,dptaper
        integer*8 ilattp,lspect,ipoltr,ipolb,ipoll,ipolid,ipolo,
     $       intffsm
        integer*4 nflag0,nlat,np0,nturn,isynch,nspect,
     $       lplot,nplot,nuse,nclas,irad,novfl,npelm,ipelm,
     $       nparallel,pspac_nx,pspac_ny,pspac_nz,
     $       pspac_nturn,pspac_nturncalc,l_track
        logical*4 oldflagsdummy,tparaed

        contains
        subroutine tinitintm
        implicit none
        intffsm=ktfsymbolz('System`FFS$InterruptMask',24)-4
        call tsetintm(0.d0)
        return
        end

        subroutine tsetintm(x)
        implicit none
        real*8 ,intent(in):: x
        rlist(intffsm)=x
        return
        end

      end module tmacro

      module tfshare
      use tfstk
      integer*4, parameter:: nshmax=1024,maxwait=1200
      integer*8, save:: kashare(nshmax)=0,lshare(nshmax)=0
      integer*4, save :: ishared(nshmax)=0,kstshare(0:nshmax)=0
      integer*4 ,parameter :: sem_init=0,sem_destroy=1,
     $     sem_post=2, sem_wait=3, sem_trywait=4
      
      contains
      integer*8 function ktfallocshared(n)
      use iso_c_binding
      implicit none
      integer*4, save :: lps=0
      integer*4 ,intent(in):: n
      integer*4 irtc,getpagesize,nsh1,i,na
      integer*8 k,kpb,kcp
      if(lps == 0)then
        lps=getpagesize()/8
      endif
      na=((n+3)/lps+2)*lps
      nsh1=0
      do i=1,kstshare(0)
        if(kstshare(i) == 0 .and. lshare(i) .ge. na)then
          kstshare(i)=1
          ktfallocshared=kashare(i)
          return
        elseif(kstshare(i) .ge. 0)then
          nsh1=i
        endif
      enddo
      kstshare(0)=nsh1
      nsh1=0
      if(kstshare(0) == nshmax)then
        do i=1,kstshare(0)
          if(kstshare(i) .le. 0)then
            if(lshare(i) > 0)then
              kstshare(i)=1
              if(nsh1 /= 0)then
                kstshare(nsh1)=1-kstshare(i)
              endif
              call tfreleaseshared(kashare(i))
              nsh1=i
              exit
            elseif(nsh1 == 0)then
              kstshare(i)=1-kstshare(i)
              nsh1=i
            endif
          endif
        enddo
      else
        nsh1=kstshare(0)+1
        kstshare(nsh1)=1
        kstshare(0)=nsh1
      endif
      k=ktaloc(na)
      kcp=transfer(c_loc(klist(k)),k)/8
      kpb=k+((kcp+1+lps)/lps)*lps-kcp
c      write(*,*)'ktfallocshared-mmap ',nsh1,kpb,na-lps
      call mapallocshared8(klist(kpb),int8(na-lps),8,irtc)
      if(irtc /= 0)then
        write(*,*)'ktfallocshared failed: ',kpb,na-lps
        call abort
      endif
      ktfallocshared=kpb+2
      klist(kpb)=k
      klist(kpb+1)=na-lps
      if(nsh1 /= 0)then
        ishared(nsh1)=1
        kstshare(nsh1)=1
        kashare(nsh1)=kpb+2
        lshare(nsh1)=na
      endif
c      write(*,*)'allocshared ',kpb,na,lps,
c     $     transfer(c_loc(klist(kpb)),k)/8
      return
      end function

      subroutine tfreeshared(kpb,ist)
      implicit none
      integer*4, optional ,intent(in):: ist
      integer*4 i,is
      integer*8 ,intent(in):: kpb
      is=0
      if(present(ist))then
        is=ist
      endif
      do i=1,kstshare(0)
        if(kashare(i) == kpb)then
          kstshare(i)=is
          if(is .lt. 0)then
            lshare(i)=0
          endif
          return
        endif
      enddo
      call tfreleaseshared(kpb)
      return
      end subroutine

      subroutine tfreleaseshared(kpb)
      implicit none
      integer*8 ,intent(in):: kpb
      integer*8 k
      integer*4 irtc
      k=klist(kpb-2)
      call mapallocfixed8(klist(kpb-2),klist(kpb-1),8,irtc)
      if(irtc /= 0)then
        write(*,*)'tffreecshared failed: ',kpb,klist(kpb-1)
        call abort
      endif
c      write(*,*)'tfreeshared ',kpb,klist(kpb-1),irtc
      if(itfcbk(k) == 0)then
        call tfentercbk(k,klist(kpb-1)+1)
      endif
      call tfree(k)
      return
      end subroutine

      subroutine tfsavesharedmap()
      use tmacro, only:tsetintm
      implicit none
      ishared=0
      kstshare(0)=0
c      write(*,*)'savesharedmap'
      return
      end subroutine

      subroutine tfresetsharedmap()
      implicit none
      integer*4 i
      do i=1,kstshare(0)
        if(ishared(i) /= 0)then
          kstshare(i)=-1
        endif
      enddo
      return
      end subroutine

      subroutine ktfinitshare
      use iso_c_binding
      implicit none
      integer*4, save :: lps=0
      integer*4 getpagesize
      if(lps == 0)then
        lps=getpagesize()/8
      endif
      return
      end

      recursive subroutine tfrecallshared(isp0,k,kx,irtc)
      implicit none
      type (sad_descriptor) ,intent(out):: kx
      type (sad_descriptor) ,intent(in):: k
      type (sad_descriptor) k0,ki
      type (sad_string), pointer :: str
      type (sad_dlist), pointer :: kl
      integer*8 kax
      integer*4 ,intent(in):: isp0
      integer*4 ,intent(out):: irtc
      integer*4 m,itfmessage,i
      logical*4 tfcheckelement
      do i=isp0+1,isp
        if(ktastk(i) == k%k)then
          kx=dtastk2(i)
          return
        endif
      enddo
      if(ktfstringq(k,str))then
        kx=kxsalocb(-1,str%str,str%nch)
      elseif(ktfsymbolq(k))then
        if( .not. tfcheckelement(k,.false.))then
          irtc=itfmessage(9,'General::wrongtype',
     $         '"defined symbol returned in Shared"')
          return
        endif
        kx=k
      elseif(ktflistq(k,kl))then
c        call tfdebugprint(k,'recallshared',3)
        k0=kl%head
        if(ktfobjq(k0))then
          call tfrecallshared(isp0,k0,k0,irtc)
          if(irtc /= 0)then
            return
          endif
        endif
        m=kl%nl
        if(kl%ref == 0)then
          kax=ktavaloc(-1,m)
          dlist(kax+1:kax+m)=kl%dbody(1:m)
        else
          kax=ktadaloc(-1,m)
          do i=1,m
            ki=kl%dbody(i)
            if(ktfobjq(ki))then
              call tfrecallshared(isp0,ki,ki,irtc)
              if(irtc /= 0)then
                klist(kax+1:kax+m)=ktfoper+mtfnull
                exit
              endif
              dlist(kax+i)=dtfcopy1(ki)
            else
              dlist(kax+i)=ki
            endif
          enddo
        endif
        dlist(kax)=dtfcopy(k0)
        kx%k=ktflist+kax
      else
        kx%k=ktfoper+mtfnull
      endif
      isp=isp+1
      dtastk(isp)=k
      dtastk2(isp)=kx
      irtc=0
      return
      end

      subroutine tfstoreshared(isp0,k,kap)
      implicit none
      integer*8 ,intent(in):: k,kap
      integer*8 ka,kh,ki,kt
      integer*4 ,intent(in):: isp0
      integer*4 i,j,m
      ka=ktfaddr(k)
      kt=k-ka
      if(kt == ktfstring)then
        klist(kap)=klist(ka)
        do i=1,(ilist(1,ka)+7)/8
          rlist(kap+i)=rlist(ka+i)
        enddo
      elseif(kt == ktflist)then
        kh=klist(ka)
        klist(kap+1)=kh
        if(ktfobjq(kh) .and. ktfnonsymbolq(kh))then
          do j=isp0,isp
            if(ktastk(j) == kh)then
              klist(kap+1)=merge(ktfstring+ktastk2(j),
     $             ktflist+ktastk2(j)+1,ktfstringq(kh))
              exit
            endif
          enddo
        endif
        m=ilist(2,ka-1)
        ilist(2,kap)=m
        if(ktfreallistq(ka))then
          ilist(1,kap)=0
          rlist(kap+2:kap+m+1)=rlist(ka+1:ka+m)
c          do i=1,m
c            rlist(kap+i+1)=rlist(ka+i)
c          enddo
        else
          ilist(1,kap)=1
          do i=1,m
            ki=klist(ka+i)
            if(ktfnonobjq(ki) .or. ktfsymbolq(ki))then
              klist(kap+i+1)=ki
            else
              do j=isp0,isp
                if(ktastk(j) == ki)then
                  if(ktfstringq(ki))then
                    klist(kap+i+1)=ktfstring+ktastk2(j)
                  else
                    klist(kap+i+1)=ktflist+ktastk2(j)+1
                  endif
                  exit
                endif
              enddo
            endif
          enddo
        endif
      endif
      return
      end

      recursive subroutine tfsharedsize(isp0,k,n,irtc)
      implicit none
      integer*8 ,intent(in):: k
      integer*8 ka,kt,ki,kh
      integer*4 ,intent(in):: isp0
      integer*4 ,intent(out):: irtc,n
      integer*4 i,itfmessage,ni
      irtc=0
      if(ktfnonobjq(k) .or. ktfsymbolq(k))then
        n=1
      else
        ka=ktfaddr(k)
        kt=k-ka
        if(kt == ktfstring)then
          n=3+(ilist(1,ka)+7)/8
          do i=isp0+1,isp
            if(ktastk(i) == k)then
              n=1
              exit
            endif
          enddo
          isp=isp+1
          ktastk(isp)=k
          itastk2(1,isp)=n
        elseif(kt == ktflist)then
          do i=isp0+1,isp
            if(ktastk(i) == k)then
              n=1
              return
            endif
          enddo
          n=3+ilist(2,ka-1)
          isp=isp+1
          ktastk(isp)=k
          itastk2(1,isp)=n
          kh=klist(ka)
          if(ktfobjq(kh) .and. ktfnonsymbolq(kh))then
            call tfsharedsize(isp0,klist(ka),ni,irtc)
            if(irtc /= 0)then
              return
            endif
            n=n+ni-1
          endif
          if(ktfnonreallistq(ka))then
            do i=1,ilist(2,ka-1)
              ki=klist(ka+i)
              if(ktfobjq(ki) .and. ktfnonsymbolq(ki))then
                call tfsharedsize(isp0,ki,ni,irtc)
                if(irtc /= 0)then
                  return
                endif
                n=n+ni-1
              endif
            enddo
          endif
        else
          go to 9000
        endif
      endif
      return
 9000 irtc=itfmessage(9,'General::wrongtype',
     $     '"No Symbols nor Patterns"')
      return
      end

      integer*4 function itffork() result(ip)
      use tmacro, only:tsetintm
      implicit none
      integer*4 :: fork_worker,i,npause=5000
      integer*4 ,parameter :: nrpt=10
      do i=1,nrpt
        ip=fork_worker()
        if(ip /= -1)then
          exit
        endif
        call tpause(npause)
        npause=npause*2
      enddo
      if(ip == 0)then
        call tsetintm(-1.d0)
        call tfsavesharedmap()
      endif
      return
      end function

      subroutine tffswait(ipr,npa,npr,kash,nwait,tag,irtc)
      implicit none
      integer*4 ,intent(out):: irtc
      integer*4 ,intent(inout):: ipr
      integer*4 ,intent(in):: npa,nwait
      integer*4 ,intent(inout):: npr(npa)
      integer*4 ist,i,j,waitpid,waitpid_nohang,iwait,lw
      integer*8 ,intent(in):: kash
      character*(*) ,intent(in):: tag
      if(ipr == 0)then
        if(kash /= 0)then
          call tfreeshared(kash,-1)
        endif
c        call tfresetsharedmap()
        call exit_without_hooks(0)
      elseif(ipr > 0)then
        iwait=-1
        if(nwait /= 0)then
          lw=maxwait*(1000000/nwait)
        endif
        irtc=0
        do i=1,npa-1
          dowait: do
            if(nwait == 0)then
              ipr=waitpid(-1,ist)
            else
              ipr=waitpid_nohang(-1,ist)
              if(ipr == 0)then
                call tpause(nwait)
                if(iwait >= 0)then
                  iwait=iwait+1
                  if(iwait > lw)then
                    do j=1,npa-1
                      if(npr(j) /= 0)then
                        call tkill(npr(j))
                        write(*,*)'???-'//tag//' timeout :',j
                      endif
                    enddo
                    irtc=20004
                    return
                  endif
                endif
                cycle dowait
              endif
            endif
            do j=1,npa-1
              if(npr(j) == ipr)then
                iwait=0
                npr(j)=0
                exit dowait
              endif
            enddo
            if(ipr /= -1)then
              write(*,*)'???-'//tag//' Unexpected process:',ipr
            endif
          enddo dowait
          if(ist /= 0)then
            irtc=20003
          endif
        enddo 
      endif
      return
      end subroutine

      integer*8 function itmmapp(n)
      use tfmem
      use tmacro
      implicit none
      integer*4 ,intent(in):: n
      integer*4 irtc
      if(nparallel > 1)then
        irtc=1
        itmmapp=ktfallocshared(n)
      else
        itmmapp=ktaloc(n)
      endif
      return
      end function

      subroutine tmunmapp(i)
      use tmacro
      implicit none
      integer*8 ,intent(in):: i
      if(nparallel > 1)then
        call tfreeshared(i)
      else
        call tfree(i)
      endif
      return
      end subroutine

      end module tfshare

      function ktfallocsharedf(n) result(k)
      use tfshare
      implicit none
      integer*8 k
      integer*4 ,intent(in):: n
      k=ktfallocshared(n)
      return
      end

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
c        write(*,*)'nextfn ',f,mode,itbuf(f)
        if(itbuf(f) == modeclose)then
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
      if(ibuf(lfn) /= 0)then
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
      type (sad_string) ,pointer :: str
      integer*4 lfn,irtc,unixclose
      character*256 cm
      select case (itbuf(lfn))
      case (modewrite)
        close(lfn)
        if(ibuf(lfn) > 0)then
          if(ilist(2,ibuf(lfn)-1) /= 0)then
            irtc=unixclose(ilist(2,ibuf(lfn)-1))
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
        call tfreeshared(ibuf(lfn),0)
        ibuf(lfn)=0
      case default
        call unmapfile(klist(ibuf(lfn)),int8(lenbuf(lfn)))
        irtc=unixclose(ifd(lfn))
c        write(*,'(a,5i10)')'trbclose ',lfn,itbuf(lfn),
c     $       lenbuf(lfn),ifd(lfn),irtc
      end select
      ifd(lfn)=0
      lbuf(lfn)=0
      mbuf(lfn)=1
      itbuf(lfn)=modeclose
      if(ktfstringq(ntable(lfn),str))then
c        call tfdebugprint(ntable(lfn),'trbclose',1)
        cm(:str%nch+3)='rm '//str%str(:str%nch)
c        write(*,*)'trbclose ',cm(1:str%nch+3)
        call system(cm(:str%nch+3))
        call tflocald(ntable(lfn))
      endif
      ntable(lfn)%k=0
      return
      end subroutine

      subroutine trbnextl(lfn)
      implicit none
      integer*4 lfn
      if(lfn > 0)then
        mbuf(lfn)=lbuf(lfn)+1
      endif
      return
      end subroutine

      subroutine trbeor2bor(lfn)
      implicit none
      integer*4 lfn
      if(lfn > 0)then
        if(mbuf(lfn) == lbuf(lfn))then
          mbuf(lfn)=lbuf(lfn)+1
        endif
      endif
      return
      end subroutine

      integer*8 function itrbibuf(lfn,mode) result(ia)
      implicit none
      integer*4 , intent(in) :: lfn,mode
      if(lfn > 0 .and. itbuf(lfn) == mode)then
        ia=ibuf(lfn)
      else
        ia=0
      endif
      return
      end function

      subroutine trbmovepoint(lfn,nc)
      implicit none
      integer*4 , intent(in) :: lfn,nc
      if(lfn > 0)then
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
        if(ibuf(lfn) == 0)then
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
      if(lfn .le. 0 .or. ibuf(lfn) == 0)then
        go to 9000
      endif
      if(itbuf(lfn) .le. modewrite)then
c        write(*,*)'tfreadbuf ',lfn,itbuf(lfn),modewrite
        nc=itfgetbuf(lfn,jlist(ib,ibuf(lfn)),
     $       maxlbuf-ib-256,irtc)
        if(irtc /= 0)then
          return
        endif
        lbuf(lfn)=ib-1
        lenbuf(lfn)=lbuf(lfn)+nc
c        write(*,*)'tfreadbuf ',lfn,nc,lbuf(lfn),mbuf(lfn),ib
        mbuf(lfn)=ib
      else
 11     ls=lenbuf(lfn)
c        if(lbuf(lfn) .lt. ls .and.
c     $       jlist(lbuf(lfn)+1,ibuf(lfn)) == 10)then
c          lbuf(lfn)=lbuf(lfn)+1
c        endif
        if(lbuf(lfn) .lt. ls)then
          ie=ls
          do i=lbuf(lfn)+1,ls
            if(jlist(i,ibuf(lfn)) == 10)then
              if(i == 1 .or. jlist(i-1,ibuf(lfn)) /=
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
      if(itbuf(j) == 0)then
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
      if(nc > 0)then
        str=buffer(mbuf(in):mbuf(in)+nc-1)
        mbuf(in)=mbuf(in)+nc
        if(str(nc:nc) == char(10))then
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
      if(isp /= isp1+1)then
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
      if(irtc /= 0)then
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
      use tfshare
      implicit none
      type (sad_descriptor) kx
      integer*4 ,intent(in):: isp1
      integer*4 ,intent(out):: irtc
      type (sad_string), pointer :: str
      integer*8 ia
      integer*4 itfmessage,isp0,iu,ist
      logical*4 tfcheckelement
      if(isp /= isp1+1)then
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
      if(ia == 0)then
        kx%k=kxeof
        return
      endif
      ist=ilist(2,ia)
      do while(ist /= 0)
        if(ist /= 1 .and. ist /= -1)then
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
          call tfrecallshared(isp0,dfromk(ktflist+ia+3),kx,irtc)
          isp=isp0
          if(irtc /= 0)then
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
      use tfshare, only:tfsharedsize,tfstoreshared
      implicit none
      type (sad_descriptor) kx
      integer*4 ,intent(in):: isp1
      integer*4 ,intent(out):: irtc
      integer*8 kas,ka,kt,kap,k
      integer*4 itfmessage,itfmessageexp,isp0,n,i,iu,ist
      kx=dxnullo
      if(isp /= isp1+2)then
        irtc=itfmessage(9,'General::narg','"2"')
        return
      elseif(ktfnonrealq(ktastk(isp1+1),iu))then
        irtc=itfmessage(9,'General::wrongtype','"Real"')
        return
      endif
      kas=itrbibuf(iu,modeshared)
      if(kas == 0)then
        irtc=itfmessage(99,'Shared::notopen','""')
        return
      endif
      ist=ilist(2,kas)
      do while(ist /= 0)
        if(ist /= 1 .and. ist /= -1)then
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
        if(kt == ktfstring)then
          if(ilist(1,ka) > ilist(1,kas)*8)then
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
        elseif(kt == ktflist)then
          isp0=isp
          call tfsharedsize(isp0,k,n,irtc)
          if(irtc /= 0)then
            ilist(2,kas)=0
            isp=isp1+2
            return
          endif
          if(n > ilist(1,kas))then
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
        if(itbuf(j) == modeclose)then
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
      if(irtc == 0)then
        call trbopen(lfn,kfile/8,ksize+modemapped,ifd)
        kx%x(1)=dble(lfn)
      else
c        write(*,*)'trbopenmap ',irtc,ifd,ksize,str
        kx%x(1)=-1.d0
      endif
      return
      end subroutine

      end module readbuf

      module sad_basics
      use mathfun, only: sqrtl

      contains
      subroutine tdrife1(trans1,cod,dv,dvdp,al)
      implicit none
      real*8 ,intent(in):: al,dv,dvdp
      real*8 ,intent(out):: trans1(6,6)
      real*8 ,intent(out):: cod(6)
      real*8 pxi,pyi,a,pzi,ale,alz,pr
      pr=1.d0+cod(6)
      pxi=cod(2)
      pyi=cod(4)
      a=pxi**2+pyi**2
      pzi=pr*sqrtl(1.d0-a/pr**2)
      ale=al/pzi
      alz=ale/pzi**2
      call tinitr(trans1)
      trans1(1,2)=ale+pxi**2*alz
      trans1(1,4)=pxi*pyi*alz
      trans1(1,6)=-pxi*pr*alz
      trans1(3,2)=trans1(1,4)
      trans1(3,4)=ale+pyi**2*alz
      trans1(3,6)=-pyi*pr*alz
      trans1(5,2)=trans1(1,6)
      trans1(5,4)=trans1(3,6)
      trans1(5,6)=dvdp*al+a*alz
c     trans(1,1:irad)=trans(1,1:irad)+trans1(1,2)*trans(2,1:irad)
c     $       +trans1(1,4)*trans(4,1:irad)+trans1(1,6)*trans(6,1:irad)
c     trans(3,1:irad)=trans(3,1:irad)+trans1(3,2)*trans(2,1:irad)
c     $       +trans1(3,4)*trans(4,1:irad)+trans1(3,6)*trans(6,1:irad)
c     trans(5,1:irad)=trans(5,1:irad)+trans1(5,2)*trans(2,1:irad)
c     $       +trans1(5,4)*trans(4,1:irad)+trans1(5,6)*trans(6,1:irad)
      cod(1)=cod(1)+pxi/pzi*al
      cod(3)=cod(3)+pyi/pzi*al
      cod(5)=cod(5)-(a/(pr+pzi)/pzi+dv)*al
      return
      end

      subroutine tinitr(trans)
      implicit none
      real*8 ,intent(out):: trans(36)
      trans=(/
     $     1.d0,0.d0,0.d0,0.d0,0.d0,0.d0,
     $     0.d0,1.d0,0.d0,0.d0,0.d0,0.d0,
     $     0.d0,0.d0,1.d0,0.d0,0.d0,0.d0,
     $     0.d0,0.d0,0.d0,1.d0,0.d0,0.d0,
     $     0.d0,0.d0,0.d0,0.d0,1.d0,0.d0,
     $     0.d0,0.d0,0.d0,0.d0,0.d0,1.d0/)
      return
      end subroutine

      subroutine tinitr12(trans)
      implicit none
      real*8 ,intent(out):: trans(72)
      trans=0.d0
      trans(1 )=1.d0
      trans(8 )=1.d0
      trans(15)=1.d0
      trans(22)=1.d0
      trans(29)=1.d0
      trans(36)=1.d0
      return
      end subroutine

      subroutine tmultr(trans,trans1,n)
      implicit none
      integer*4 n
      real*8,intent(inout):: trans(6,n)
      real*8 ,intent(in)::trans1(6,6)
c the warning generated by gfortran up tp 10 is said to be a bug in gcc....
      trans=matmul(trans1,trans)
      return
      end subroutine

      subroutine tmultr5(trans,trans1,n)
      implicit none
      integer*4 n
      real*8 , intent(inout) ::trans(6,n)
      real*8 , intent(in)::trans1(6,6)
      real*8 v1(n),v2(n),v3(n),v4(n),v6(n)
      v1=trans(1,:)
      v2=trans(2,:)
      v3=trans(3,:)
      v4=trans(4,:)
      v6=trans(6,:)
      trans(1,:)=trans1(1,1)*v1+trans1(1,2)*v2+trans1(1,3)*v3
     1     +trans1(1,4)*v4+trans1(1,6)*v6
      trans(2,:)=trans1(2,1)*v1+trans1(2,2)*v2+trans1(2,3)*v3
     1     +trans1(2,4)*v4+trans1(2,6)*v6
      trans(3,:)=trans1(3,1)*v1+trans1(3,2)*v2+trans1(3,3)*v3
     1     +trans1(3,4)*v4+trans1(3,6)*v6
      trans(4,:)=trans1(4,1)*v1+trans1(4,2)*v2+trans1(4,3)*v3
     1     +trans1(4,4)*v4+trans1(4,6)*v6
      trans(5,:)=trans1(5,1)*v1+trans1(5,2)*v2+trans1(5,3)*v3
     1     +trans1(5,4)*v4+trans(5,:)+trans1(5,6)*v6
      return
      end subroutine

      real*8 function tinv6(ra) result(rinv)
      implicit none
      real*8, intent(in):: ra(6,6)
      dimension rinv(6,6)
      rinv(1,1)= ra(2,2)
      rinv(1,2)=-ra(1,2)
      rinv(1,3)= ra(4,2)
      rinv(1,4)=-ra(3,2)
      rinv(1,5)= ra(6,2)
      rinv(1,6)=-ra(5,2)
      rinv(2,1)=-ra(2,1)
      rinv(2,2)= ra(1,1)
      rinv(2,3)=-ra(4,1)
      rinv(2,4)= ra(3,1)
      rinv(2,5)=-ra(6,1)
      rinv(2,6)= ra(5,1)
      rinv(3,1)= ra(2,4)
      rinv(3,2)=-ra(1,4)
      rinv(3,3)= ra(4,4)
      rinv(3,4)=-ra(3,4)
      rinv(3,5)= ra(6,4)
      rinv(3,6)=-ra(5,4)
      rinv(4,1)=-ra(2,3)
      rinv(4,2)= ra(1,3)
      rinv(4,3)=-ra(4,3)
      rinv(4,4)= ra(3,3)
      rinv(4,5)=-ra(6,3)
      rinv(4,6)= ra(5,3)
      rinv(5,1)= ra(2,6)
      rinv(5,2)=-ra(1,6)
      rinv(5,3)= ra(4,6)
      rinv(5,4)=-ra(3,6)
      rinv(5,5)= ra(6,6)
      rinv(5,6)=-ra(5,6)
      rinv(6,1)=-ra(2,5)
      rinv(6,2)= ra(1,5)
      rinv(6,3)=-ra(4,5)
      rinv(6,4)= ra(3,5)
      rinv(6,5)=-ra(6,5)
      rinv(6,6)= ra(5,5)
      return
      end function

      real*8 function tsymp(trans) result(rx)
      implicit none
      real*8 , intent(in)::trans(6,6)
      dimension rx(6,6)
c  gfortran up to 10 generates a warning on uninitialized..., which is said to be a bug in gcc...
      rx=matmul(trans,tinv6(trans))
      rx=-0.5d0*matmul(rx,trans)+1.5d0*trans
      return
      end function

      real*8 function tmultr45(a,b) result(c)
      implicit none
      real*8 ,intent(in)::a(4,5),b(4,5)
      dimension c(4,5)
c      real*8 v1,v2,v3,v4
      c=matmul(b(:,1:4),a)
      c(:,5)=c(:,5)+b(:,5)
      return 
      end function

      subroutine qcopymat(trans,transe,dir)
      implicit none
      real*8 ,intent(inout):: trans(4,5),transe(6,6)
      logical*4 ,intent(in):: dir
      if(dir)then
        transe(1,1)=trans(1,1)
        transe(1,2)=trans(1,2)
        transe(1,3)=trans(1,3)
        transe(1,4)=trans(1,4)
        transe(1,5)=0.d0
        transe(1,6)=trans(1,5)
        transe(2,1)=trans(2,1)
        transe(2,2)=trans(2,2)
        transe(2,3)=trans(2,3)
        transe(2,4)=trans(2,4)
        transe(2,5)=0.d0
        transe(2,6)=trans(2,5)
        transe(3,1)=trans(3,1)
        transe(3,2)=trans(3,2)
        transe(3,3)=trans(3,3)
        transe(3,4)=trans(3,4)
        transe(3,5)=0.d0
        transe(3,6)=trans(3,5)
        transe(4,1)=trans(4,1)
        transe(4,2)=trans(4,2)
        transe(4,3)=trans(4,3)
        transe(4,4)=trans(4,4)
        transe(4,5)=0.d0
        transe(4,6)=trans(4,5)
        transe(5,1)=0.d0
        transe(5,2)=0.d0
        transe(5,3)=0.d0
        transe(5,4)=0.d0
        transe(5,5)=1.d0
        transe(5,6)=0.d0
        transe(6,1)=0.d0
        transe(6,2)=0.d0
        transe(6,3)=0.d0
        transe(6,4)=0.d0
        transe(6,5)=0.d0
        transe(6,6)=1.d0
      else
        trans(1,1)=transe(1,1)
        trans(1,2)=transe(1,2)
        trans(1,3)=transe(1,3)
        trans(1,4)=transe(1,4)
        trans(1,5)=transe(1,6)
        trans(2,1)=transe(2,1)
        trans(2,2)=transe(2,2)
        trans(2,3)=transe(2,3)
        trans(2,4)=transe(2,4)
        trans(2,5)=transe(2,6)
        trans(3,1)=transe(3,1)
        trans(3,2)=transe(3,2)
        trans(3,3)=transe(3,3)
        trans(3,4)=transe(3,4)
        trans(3,5)=transe(3,6)
        trans(4,1)=transe(4,1)
        trans(4,2)=transe(4,2)
        trans(4,3)=transe(4,3)
        trans(4,4)=transe(4,4)
        trans(4,5)=transe(4,6)
      endif
      return
      end subroutine

      subroutine tnorm(r,ceig,lfno)
      implicit none
      integer*4 i,j
      integer*4 ,intent(in):: lfno
      real*8 ,intent(inout):: r(6,6)
      real*8 sa,sb,a,s,cost,sint,s1,s2,s3,b,v(6,2)
      complex*16 ,intent(out):: ceig(6)
      complex*16 cc
      do i=1,5,2
        if(imag(ceig(i)) == 0.d0)then
          sa=sum(r(:,i)**2)
          sb=sum(r(:,i+1)**2)
          a=sqrt(sqrt(sa/sb))
          r(:,i  )=r(:,i  )/a
          r(:,i+1)=r(:,i+1)*a
        endif
        s=0.d0
        do j=1,5,2
          s=s+r(j,i)*r(j+1,i+1)-r(j,i+1)*r(j+1,i)
        enddo
        sa=sqrt(abs(s))
        if(s > 0.d0)then
          sb=sa
        elseif(s == 0.d0)then
          if(lfno /= 0)then
            write(lfno,*)'???-tnorm-Unstable transfer matrix.'
          endif
          sa=1.d0
          sb=1.d0
        else
          sb=-sa
          ceig(i  )=conjg(ceig(i  ))
          ceig(i+1)=conjg(ceig(i+1))
        endif
        r(:,i  )=r(:,i  )/sa
        r(:,i+1)=r(:,i+1)/sb
      enddo
      s1=(r(5,1)*r(6,2)-r(5,2)*r(6,1))**2
      s2=(r(5,3)*r(6,4)-r(5,4)*r(6,3))**2
      s3=(r(5,5)*r(6,6)-r(5,6)*r(6,5))**2
      j=maxloc((/s1,s2,s3/),1)*2-1
      if(j /= 5)then
        v=r(:,5:6)
        r(:,5:6)=r(:,j:j+1)
        r(:,j:j+1)=v
        cc=ceig(5)
        ceig(5)=ceig(j)
        ceig(j)=cc
        cc=ceig(6)
        ceig(6)=ceig(j+1)
        ceig(j+1)=cc
      endif
      s1=(r(1,1)*r(2,2)-r(1,2)*r(2,1))**2
      s2=(r(1,3)*r(2,4)-r(1,4)*r(2,3))**2
      if(s2 > s1)then
        v=r(:,1:2)
        r(:,1:2)=r(:,3:4)
        r(:,3:4)=v
        cc=ceig(1)
        ceig(1)=ceig(3)
        ceig(3)=cc
        cc=ceig(2)
        ceig(2)=ceig(4)
        ceig(4)=cc
      endif
      do i=1,5,2
        if(imag(ceig(i)) /= 0.d0)then
          a=hypot(r(i,i),r(i,i+1))
          cost=r(i,i  )/a
          sint=r(i,i+1)/a
          v(:,1)=r(:,i)
          r(:,i  )= v(:,1)*cost+r(:,i+1)*sint
          r(:,i+1)=-v(:,1)*sint+r(:,i+1)*cost
        else
          a=r(i,i)
          if(a /= 0.d0)then
            b=r(i,i+1)
            r(:,i+1)=r(:,i+1)*a-b*r(:,i)
            r(:,i  )=r(:,i  )/a
          endif
        endif
      enddo
c     call tsymp(r)
      return
      end subroutine

      end module sad_basics

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
      if (IgetGL('$CTIME$') == FLAGON) call cputix
      call gettok(token,slen,ttype,rval,ival)
c     for debug
c       print *,'toplvl-0 ',token(:slen),slen,infl,ttype,ttypEF
c     end debug
c
 1100 continue
      if (ttype == ttypEF) go to 9000
      if(ttype == ttypID) then
c..........System defined name.
         index=hsrch(token(:slen))
c         write(*,*)'toplvl-1 ',index,idtype(index),idval(index),icRSVD
         if (idtype(index) == icDEF) then
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
         else if (idtype(index) == icACT) then
            call tfinitstk
            if(bypasstrack)then
              return
            endif
            call doACT(index)
            go to 1000
         else if (idtype(index) == icRSVD) then
           call tfinitstk
           call funcall1(idval(index),argp)
           go to 1000
         endif
c..........User defined variable.
c     print *,token(:slen)
         index=hsrchz(token(:slen))
         if ((idtype(index) .lt. icMXEL) .and.
     &        (idtype(index) /= icNULL)) then
           call tfinitstk
           call doelm2(idtype(index),pname(index),slen,ttype)
         else if ((idtype(index) == icLINE)) then
           call tfinitstk
           call dolin2(index,slen,ttype)
         else
            call tfinitstk
           call dAssgn(token,slen,status)
           if (status /= 0) call errmsg('toplvl',
     &          'Unsupported function '//token(:slen)//'!',
     &          0,16)
c     for debug
c      print *,token(:slen),
c     &     index,pname(index),idval(index),idtype(index)
c     end debug
         endif
      else if(token(:slen) == ';') then
         go to 1000
      else if(token(:slen) == RCURL) then
         go to 1000
       elseif(token(:slen) == char(13))then
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
      if(itbuf(infl) /= 0 .or. lfnbase > 1)then
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
