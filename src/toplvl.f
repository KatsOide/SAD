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
     $     kwSIGZ=kwPY+1,kwAZ=kwSIGZ+1,kwDPX=kwAZ+1,kwDPY=kwDPX+1,
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

      module macmath
c     Math constants

c     Napier's constant(Exp[1])
c     \lim_{n\rightarrow\infty}(1 + \frac{1}{n})^n
      real*8, parameter :: m_e        = 2.7182818284590452354d0

c     Euler-Mascheroni constant(Euler's gamma)
c     \lim_{n\rightarrow\infty}(\sum_{k=1}^{n}\frac{1}{k} - \ln{n})
      real*8, parameter :: m_euler    = 0.57721566490153286061d0

c     Log[2,  Exp[1]]
      real*8, parameter :: m_log2e    = 1.4426950408889634074d0

c     Log[10, Exp[1]]
      real*8, parameter :: m_log10e   = 0.43429448190325182765d0

c     Log[Exp[1],  2]
      real*8, parameter :: m_ln2      = 0.69314718055994530942d0

c     Log[Exp[1], 10]
      real*8, parameter :: m_ln10     = 2.30258509299404568402d0

c         Pi
      real*8, parameter :: m_pi       = 3.14159265358979323846d0

c     2 * Pi
      real*8, parameter :: m_2pi      = 6.28318530717958647693d0

c     4 * Pi
      real*8, parameter :: m_4pi      = 12.5663706143591729539d0

c         Pi / 2
      real*8, parameter :: m_pi_2     = 1.57079632679489661923d0

c         Pi / 4
      real*8, parameter :: m_pi_4     = 0.78539816339744830962d0

c     1 / Pi
      real*8, parameter :: m_1_pi     = 0.31830988618379067154d0

c     2 / Pi
      real*8, parameter :: m_2_pi     = 0.63661977236758134308d0

c     4 / Pi
      real*8, parameter :: m_4_pi     = 1.27323954473516268615d0

c         Sqrt[Pi]
      real*8, parameter :: m_sqrtpi   = 1.77245385090551602730d0

c     2 / Sqrt[Pi]
      real*8, parameter :: m_2_sqrtpi = 1.12837916709551257390d0

c     Sqrt[2]
      real*8, parameter :: m_sqrt2    = 1.41421356237309504880d0

c     Sqrt[1/2]
      real*8, parameter :: m_1_sqrt2  = 0.70710678118654752440d0

c     Sqrt[3]
      real*8, parameter :: m_sqrt3    = 1.73205080756887729353d0

c     Sqrt[1/3]
      real*8, parameter :: m_1_sqrt3  = 0.57735026918962576451d0

c     Standard alias for SAD sources
      real*8, parameter :: napier = m_e, euler = m_euler
      real*8, parameter :: pi  = m_pi
      real*8, parameter :: pi2 = m_2pi,  pi4 = m_4pi
      real*8, parameter :: hpi = m_pi_2, qpi = m_pi_4

      end module

      module macphys
      use macmath
c DO not forget to update sim/MACPHYS.h!!!
c
c     Update History:
c     2016/01/12:	from PDG 2014 https://www.google.co.jp/url?sa=t&rct=j&q=&esrc=s&source=web&cd=2&ved=0ahUKEwjHw4unsKLKAhUEoJQKHaiPCLgQFggjMAE&url=http%3A%2F%2Fpdg.lbl.gov%2F2014%2Freviews%2Frpp2014-rev-phys-constants.pdf&usg=AFQjCNGXw6APvgkpoTdIRLeC-XV651ZbJg&sig2=-wjwSn5pfl6Juvok0hFqow
c     2008/03/27:	from CODATA 2006
c     			http://physics.nist.gov/cuu/Constants/Table/allascii.txt
c
c     2003/07/15:	from Particle Data Book
c

c     Including math constant for permeability of vacuum

c     Speed of light in vacuum:		c
c     Definition:	299 792 458. m/sec
      real*8, parameter :: cveloc = 299792458d0

c     Planck's constant:		h
c     CODATA 2006:	6.626 068 96(33)   x 10^-34 Js
c     PDG 2014::	6.626 069 57(29)  x 10^-34 Js
c     CODATA 2018       6.626 070 15       x 10^-34 Js exact
      real*8, parameter :: plank  = 6.62607015d-34

c     Dirac's constant:			hbar = h / (2 Pi)
c     PDB 2003:		1.054 571 596(82)  x 10^-34 Js
c     CODATA 2006:	1.054 571 628(53)  x 10^-34 Js
c     PDG 2014: 	1.054 571 726(47)  x 10^-34 Js
c     CODATA 2018:      h (exact) / 2pi
      real*8, parameter :: plankr = plank/m_2pi

c     Elementary charge:		e
c     PDB 2003:		1.602 176 462(63)  x 10^-19 C
c     CODATA 2006:	1.602 176 487(40)  x 10^-19 C
c     PDG 2014: 	1.602 176 565(35)  x 10^-19 C
c     CODATA 2018:      1.602 176 634      x 10^-19 C exact
      real*8, parameter :: elemch = 1.602176634d-19

c     Electron charge in elementary charge unit
      real*8, parameter :: echarg = 1.0d0

c     Fine-structure constant		\alpha = \mu_0 e^2 c / (2 h)
c     PDB 2003:		1 / 137.035 999 76(50)
c     CODATA 2006:	1 / 137.035 999 679(94)
c     PDG 2014: 	1 / 137.035 999 074(44)
c     CODATA 2018: 	1 / 137.035 999 084(21)
      real*8, parameter :: finest = 1.d0 / 137.035999084d0

c     Permeability of vacuum:		\mu_0
c     Definition:	4 Pi x 10^-7 N A^-2
c     CODATA 2018:      1.000 000 000 55 * above
c     CODATA 2018:    4 pi alpha hbar/(e^2 c)
c      real*8, parameter :: mu0 = 1.00000000055d0 * pi * 4.d-7
      real*8, parameter :: mu0=m_4pi*finest*plankr/elemch**2/cveloc

c     Permittivity of vacuum:		\epsilon_0 = 1 / (\mu_0 c^2)
c     Delivered from defintion
      real*8, parameter :: ep0 = 1.d0 / (mu0 * cveloc**2)

c     Electron mass energy equivalent in eV:	m_e c^2 / e
c     PDB 2003:		.510 998 902(21) MeV
c     CODATA 2006:	.510 998 910(13) MeV
c     PDG 2014: 	.510 998 928(11) MeV
c     CODATA 2018: 	.510 998 950 00(15) MeV
      real*8, parameter :: elmass =   0.51099895000d6

c     Proton mass energy equivalent in eV:	m_p c^2 / e
c     PDB 2003:		938. 271 998(38) MeV
c     CODATA 2006:	938. 272 013(23) MeV
c     PDG 2014: 	938. 272 046(21) MeV
c     CODATA 2018: 	938. 272 088 16(29) MeV
      real*8, parameter :: prmass = 938.27208816d6

c     Classical electron radius:	r_e = e^2 / (4Pi \epsilon_0 m_e c^2)
c     					    = (\alpha hbar c) / (m_e c^2)
c     					    = (e^2 c^2 \mu_0 / 4Pi) / (m_e c^2)
c     					    = e c^2 * 10^-7 * (e / (m_e c^2))
c     PDB 2003:		2.817 940 285(??)  x 10^-15 m
c     CODATA 2006:	2.817 940 2894(58) x 10^-15 m
c     PDG 2014: 	2.817 940 3267(27) x 10^-15 m
c      parameter (elradi = finest * plankr * cveloc / (elemch * elmass))
c      parameter (elradi = elemch * cveloc**2 * 1.d-7 / elmass)
      real*8, parameter :: elradi = elemch * cveloc**2 * 1.d-7 / elmass

c     Classical proton radius:		r_p = e^2 / (4Pi \epsilon_0 m_p c^2)
c     					    = (\alpha hbar c) / (m_p c^2)
c     					    = (e^2 c^2 \mu_0 / 4Pi) / (m_p c^2)
c     					    = e c^2 * 10^-7 * (e / (m_p c^2))
c      parameter (prradi = finest * plankr * cveloc / (elemch * prmass))
      real*8, parameter :: prradi = elemch * cveloc**2 * 1.d-7 / prmass

c     Boltzmann Constant:
c     PDG2014:          1.380 6488(13) x 10^-23 J/K
c     COD2018:		1.380 649 × 10−23 J/K exact
      real*8, parameter :: kboltzman = 1.380649d-23

c     Spin precession coefficient (ge-2)/2
c     NIST 2014 0.00115965218091
c     CODATA 2018: (2.002 319 304 362 56(35))/2-1 = 0.001159652181280002
c
      real*8 , parameter :: gspin = 0.001159652181280002

      end module

      module tfcsi
        use tfcbk, only:maxlbuf
        implicit none
        integer*4, parameter :: nbmax=maxlbuf,nsav=6,ipbase=1
        type csiparam
          sequence
          integer*4 isav(1:0)
          integer*4 lfni,linep,lfn1,lfno,ios
          logical*4 rec
        end type
        type (csiparam) , target :: savep
        character*16 delim,cmnt
        integer*8 ibcloc
        integer*4, pointer:: ipoint,lrecl,lfni=>savep%lfni,
     $       linep=>savep%linep,lfn1=>savep%lfn1,lfno=>savep%lfno,
     $       ios=>savep%ios
        logical*4 , pointer :: rec=>savep%rec
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

        subroutine cssetrec(f)
        implicit none
        logical*4 f
        rec=f
        return
        end subroutine

        subroutine cssetlinep(ip)
        implicit none
        integer*4 ip
c        write(*,*)'setlinep ',ip,linep,lrecl
        linep=ip
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

        integer*4 function icslinep()
        implicit none
        icslinep=linep
        return
        end function

        logical*4 function csrec()
        implicit none
        csrec=rec
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

      subroutine trbopen(lfn,ib,is,nc)
      implicit none
      integer*8, intent(in) :: ib,is
      integer*4, intent(out):: lfn,nc
      integer*4 j
      do j=nbuf,11,-1
        if(itbuf(j) .eq. modeclose)then
          lfn=j
          call irbopen1(lfn,ib,is,nc)
          return
        endif
      enddo
      lfn=0
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

      subroutine trbopenmap(str,kx,irtc)
      use tfstk
      implicit none
      character*(*) , intent(in)::str
      integer*4 , intent(out)::irtc
      type (sad_descriptor) , intent(out)::kx
      integer*8 kfile,ksize,mapallocfile,kfromr
      integer*4 itfmessage,lfn,nc
      kfile=mapallocfile(str,nc,ksize,irtc)
      if(irtc .ne. 0)then
        irtc=itfmessage(999,'General::fileopen',str)
        kx=dxfailed
      else
        call trbopen(lfn,kfile/8,ksize+modemapped,nc)
        kx%k=kfromr(dble(lfn))
      endif
      return
      end subroutine

      end module

      module ffsfile
        integer*4 , parameter :: maxlfn=128
        integer*4 :: lfnp=0,lfnbase=0
        integer*4 lfnstk(maxlfn),lfret(maxlfn),
     $     lfrecl(maxlfn),lflinep(maxlfn)
        logical*4 lfopen(maxlfn)
      end module

      module trackbypass
        implicit none
        logical*4 :: bypasstrack=.false.
        integer*8 :: lattuse=0, lattredef=0
      end module

      subroutine toplvl
      use trackbypass, only: bypasstrack
      use tfrbuf
      use maccbk
      use mackw
      use macttyp
      use macfile
      use macmisc
      use tfstk, only:tfinitstk
      use ffsfile, only:lfnbase
      implicit none
      integer*8 argp
      integer*4 jslen,jttype,jival
      integer*4 slen,ttype,ival,hsrch,hsrchz,status
      integer*4 index,IgetGL,idummy
      real*8  jrval,rval
      character*(MAXSTR) jtoken,token
c
 1000 continue
      if (IgetGL('$CTIME$',idummy) .eq. FLAGON) call cputix
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
