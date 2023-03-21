c     Don't confuse, Emacs. This is -*- fortran -*- mode!
c$$$      integer kwL,kwANGL,kwROT,kwK1,kwDK1,kwK2,kwDK2
c$$$      integer kwK3,kwDK3,kwE1,kwE2,kwTILT,kwDX,kwDY,kwKICK
c$$$      integer kwk0,kwk4,kwk5,kwk6,kwk7,kwk8,kwk9,kwk10,
c$$$     $     kwK11,kwK12,kwK13,kwK14,kwK15,kwK16,kwK17,kwK18,
c$$$     $     kwK19,kwK20,kwK21,kwW1
c$$$      integer kwsk0,kwsk4,kwsk5,kwsk6,kwsk7,kwsk8,kwsk9,kwsk10
c$$$      integer kwsk1,kwsk2,kwsk3,kwSK11,kwSK12,kwSK13,kwSK14
c$$$      integer kwSK15,kwSK16,kwSK17,kwSK18,kwSK19,kwSK20,kwSK21
c$$$      integer kwA3,kwA5,kwA7,kwA9,kwA11,kwA13,kwA15,kwA17
c$$$      integer kwVOLT,kwLAG,kwFREQ
c$$$      integer kwAX,kwAY,kwBX,kwBY,kwEX,kwEY,kwDP,kwDDP
c$$$      integer kwBMAX,kwPRD,kwRHO0,kwFRIN,kwEPS,kwK0FR
c$$$      integer kwRAD,kwCHRO,kwDZ
c$$$      integer kwEPX,kwEPY,kwAZ,kwSIGZ,kwPX,kwPY,kwEMIX,kwEMIY
c$$$      integer kwDPX,kwDPY,kwR1,kwR2,kwR3,kwR4,kwDETR
c$$$      integer kwRHO,kwCHI1,kwCHI2,kwCHI3,kwGEO,kwBZ,kwBND
c$$$      integer kwF1,kwF2,kwFRMD,kwRANK,kwRANP,kwRANV,kwKIN
c$$$      integer kwJDX, kwJDPX, kwJDY, kwJDPY, kwJDZ, kwJDPZ
c$$$      integer kwINDX,kwSLI,kwNP,kwDIR,kwOFFSET,kwMAX
c$$$      integer kwphi, kwharm, kwlwak, kwtwak, kwDPHI
c$$$      integer kwV1,kwV20,kwV11,kwV02
c$$$      integer kwZX,  kwZPX,  kwZY,  kwZPY,  kwXANGLE, kwPZ
c$$$      integer kwEMIZ,kwSTURN,kwBSTRL
c$$$      integer*4 kwDX1,kwDX2,kwDY1,kwDY2
c$$$      integer kwD11,  kwD12,  kwD13,  kwD14,  kwD15, kwD16
c$$$      integer kwD21,  kwD22,  kwD23,  kwD24,  kwD25, kwD26
c$$$      integer kwD31,  kwD32,  kwD33,  kwD34,  kwD35, kwD36
c$$$      integer kwD41,  kwD42,  kwD43,  kwD44,  kwD45, kwD46
c$$$      integer kwD51,  kwD52,  kwD53,  kwD54,  kwD55, kwD56
c$$$      integer kwD61,  kwD62,  kwD63,  kwD64,  kwD65, kwD66
c$$$      integer kwB11,  kwB12,  kwB13,  kwB14,  kwB15, kwB16
c$$$      integer         kwB22,  kwB23,  kwB24,  kwB25, kwB26
c$$$      integer                 kwB33,  kwB34,  kwB35, kwB36
c$$$      integer                         kwB44,  kwB45, kwB46
c$$$      integer                                 kwB55, kwB56
c$$$      integer                                        kwB66
c$$$      integer kwR11,  kwR12,  kwR13,  kwR14,  kwR15, kwR16
c$$$      integer         kwR22,  kwR23,  kwR24,  kwR25, kwR26
c$$$      integer                 kwR33,  kwR34,  kwR35, kwR36
c$$$      integer                         kwR44,  kwR45, kwR46
c$$$      integer                                 kwR55, kwR56
c$$$      integer                                        kwR66
c$$$      integer kwKx,kwQy,kwpole,kwFBx,kwFBy,kwDBZ,kwCOUPLE,
c$$$     $     kwRADI,kwDROT,kwAE1,kwAE2,kwFB1,kwFB2,kwLDEV,kwLRAD,
c$$$     $     kwFL,kwAPHI
c$$$      parameter (kwL=1,kwANGL=kwL+1,kwROT=kwANGL+1)
c$$$      parameter (kwK0=kwROT+1)
c$$$      parameter (kwK1=kwK0+1,kwDK1=kwK1+1)
c$$$      parameter (kwK2=kwDK1+1,kwDK2=kwK2+1)
c$$$      parameter (kwK3=kwDK2+1)
c$$$      parameter (kwK4=kwK3+1)
c$$$      parameter (kwK5=kwK4+1)
c$$$      parameter (kwK6=kwK5+1)
c$$$      parameter (kwK7=kwK6+1)
c$$$      parameter (kwK8=kwK7+1)
c$$$      parameter (kwK9=kwK8+1)
c$$$      parameter (kwK10=kwK9+1,kwDK3=kwK10+1)
c$$$      parameter (kwA3=kwDK3+1)
c$$$      parameter (kwA5=kwA3+1)
c$$$      parameter (kwA7=kwA5+1)
c$$$      parameter (kwA9=kwA7+1)
c$$$      parameter (kwA11=kwA9+1)
c$$$      parameter (kwA13=kwA11+1)
c$$$      parameter (kwA15=kwA13+1)
c$$$      parameter (kwA17=kwA15+1)
c$$$      parameter (kwE1=kwA17+1,kwE2=kwE1+1)
c$$$      parameter (kwTILT=kwE2+1,kwKICK=kwTILT+1)
c$$$      parameter (kwDX=kwKICK+1,kwDY=kwDX+1)
c$$$      parameter (kwVOLT=kwDY+1)
c$$$      parameter (kwPHI=kwVOLT+1,kwHARM=kwPHI+1)
c$$$      parameter (kwFREQ=kwHARM+1)
c$$$      parameter (kwAX=kwFREQ+1,kwAY=kwAX+1)
c$$$      parameter (kwBX=kwAY+1,kwBY=kwBX+1)
c$$$      parameter (kwEX=kwBY+1,kwEY=kwEX+1)
c$$$      parameter (kwDP=kwEY+1)
c$$$      parameter (kwRAD=kwDP+1,kwCHRO=kwRAD+1,kwDZ=kwCHRO+1)
c$$$      parameter (kwEPX=kwDZ+1,kwEPY=kwEPX+1)
c$$$      parameter (kwPX=kwEPY+1,kwPY=kwPX+1)
c$$$      parameter (kwSIGZ=kwPY+1,kwAZ=kwSIGZ+1)
c$$$      parameter (kwDPX=kwAZ+1,kwDPY=kwDPX+1,
c$$$     $    kwR1=kwDPY+1,kwR2=kwR1+1,kwR3=kwR2+1,kwR4=kwR3+1,
c$$$     $     kwDETR=kwR4+1,
c$$$     $    kwEMIX=kwDETR+1,kwEMIY=kwEMIX+1,
c$$$     $    kwINDX=kwEMIY+1,kwBMAX=kwINDX+1,kwPRD=kwBMAX+1,
c$$$     $     kwRHO0=kwBMAX+1)
c$$$      parameter (kwRHO=kwRHO0+1,
c$$$     $    kwCHI1=kwRHO+1,kwCHI2=kwCHI1+1,kwCHI3=kwCHI2+1,kwGEO=kwCHI3+1,
c$$$     $    kwBZ=kwGEO+1,kwBND=kwBZ+1,kwFRIN=kwBND+1,kwK0FR=kwFRIN+1,
c$$$     $    kwEPS=kwK0FR+1,
c$$$     $    kwF1=kwEPS+1,kwF2=kwF1+1,kwFRMD=kwF2+1,kwRANK=kwFRMD+1,
c$$$     $    kwRANP=kwRANK+1,kwRANV=kwRANP+1,kwKIN=kwRANV+1,
c$$$     $    kwLWAK=kwKIN+1,kwTWAK=kwLWAK+1)
c$$$      parameter (kwSK1=kwTWAK+1,kwSK2=kwSK1+1,
c$$$     $     kwSK3=kwSK2+1,kwSK4=kwSK3+1,kwSK5=kwSK4+1,kwSK6=kwSK5+1,
c$$$     $     kwSK7=kwSK6+1,kwSK8=kwSK7+1,kwSK9=kwSK8+1,kwSK10=kwSK9+1)
c$$$      parameter (kwSLI=kwSK10+1,kwNP=kwSLI+1)
c$$$      parameter (kwDIR=kwNP+1,kwLAG=kwDIR+1, kwSK0=kwLAG+1,
c$$$     $     kwJDX=kwSK0+1,kwJDPX=kwJDX+1,kwJDY=kwJDPX+1,kwJDPY=kwJDY+1)
c$$$      parameter (kwK11=kwJDPY+1,kwK12=kwK11+1,kwK13=kwK12+1,
c$$$     $     kwK14=kwK13+1)
c$$$      parameter (kwSK11=kwK14+1,kwSK12=kwSK11+1,kwSK13=kwSK12+1,
c$$$     $     kwSK14=kwSK13+1)
c$$$      parameter (kwV1=kwSK14+1)
c$$$      parameter (kwZX=kwV1+1,kwZPX=kwZX+1,kwZY=kwZPX+1,kwZPY=kwZY+1,
c$$$     $     kwXANGLE=kwZPY+1,kwPZ=kwXANGLE+1,kwEMIZ=kwPZ+1,
c$$$     $     kwSTURN=kwEMIZ+1,kwBSTRL=kwSTURN+1)
c$$$      parameter (kwDX1=kwBSTRL+1,kwDX2=kwDX1+1,
c$$$     $     kwDY1=kwDX2+1,kwDY2=kwDY1+1)
c$$$      parameter (kwD11=kwDY2+1,kwD12=kwD11+1,kwD13=kwD11+2,
c$$$     $           kwD14=kwD11+3,kwD15=kwD11+4,kwD16=kwD11+5)
c$$$      parameter (kwD21=kwD11+6,kwD22=kwD21+1,kwD23=kwD21+2,
c$$$     $           kwD24=kwD21+3,kwD25=kwD21+4,kwD26=kwD21+5)
c$$$      parameter (kwD31=kwD21+6,kwD32=kwD31+1,kwD33=kwD31+2,
c$$$     $           kwD34=kwD31+3,kwD35=kwD31+4,kwD36=kwD31+5)
c$$$      parameter (kwD41=kwD31+6,kwD42=kwD41+1,kwD43=kwD41+2,
c$$$     $           kwD44=kwD41+3,kwD45=kwD41+4,kwD46=kwD41+5)
c$$$      parameter (kwD51=kwD41+6,kwD52=kwD51+1,kwD53=kwD51+2,
c$$$     $           kwD54=kwD51+3,kwD55=kwD51+4,kwD56=kwD51+5)
c$$$      parameter (kwD61=kwD51+6,kwD62=kwD61+1,kwD63=kwD61+2,
c$$$     $           kwD64=kwD61+3,kwD65=kwD61+4,kwD66=kwD61+5)
c$$$c
c$$$      parameter (kwB11=kwD61+6,kwB12=kwB11+1,kwB13=kwB11+2,
c$$$     $           kwB14=kwB11+3,kwB15=kwB11+4,kwB16=kwB11+5)
c$$$      parameter (kwB22=kwB11+6,kwB23=kwB11+7,kwB24=kwB11+8,
c$$$     $           kwB25=kwB11+9,kwB26=kwB11+10)
c$$$      parameter (kwB33=kwB11+11,kwB34=kwB11+12,
c$$$     $           kwB35=kwB11+13,kwB36=kwB11+14)
c$$$      parameter (kwB44=kwB11+15,kwB45=kwB11+16,kwB46=kwB11+17)
c$$$      parameter (kwB55=kwB11+18,kwB56=kwB11+19,kwB66=kwB11+20)
c$$$      parameter (kwR11=kwB11+21,kwR12=kwR11+1,kwR13=kwR11+2,
c$$$     $           kwR14=kwR11+3,kwR15=kwR11+4,kwR16=kwR11+5)
c$$$      parameter (kwR22=kwR11+6,kwR23=kwR11+7,kwR24=kwR11+8,
c$$$     $           kwR25=kwR11+9,kwR26=kwR11+10)
c$$$      parameter (kwR33=kwR11+11,kwR34=kwR11+12,
c$$$     $           kwR35=kwR11+13,kwR36=kwR11+14)
c$$$      parameter (kwR44=kwR11+15,kwR45=kwR11+16,kwR46=kwR11+17)
c$$$      parameter (kwR55=kwR11+18,kwR56=kwR11+19,kwR66=kwR11+20)
c$$$      parameter (kwKx=kwR11+21,kwQy=kwkx+1,kwpole=kwQy+1)
c$$$      parameter (kwFBx=kwpole+1,kwFBy=kwFBx+1)
c$$$      parameter (kwDBZ=kwFBy+1)
c$$$      parameter (kwK15=kwDBZ+1,kwK16=kwK15+1,kwK17=kwK16+1,
c$$$     $           kwK18=kwK17+1,kwK19=kwK18+1,kwK20=kwK19+1,
c$$$     $           kwK21=kwK20+1)
c$$$      parameter (kwSK15=kwK21+1,kwSK16=kwSK15+1,kwSK17=kwSK16+1,
c$$$     $           kwSK18=kwSK17+1,kwSK19=kwSK18+1,kwSK20=kwSK19+1,
c$$$     $           kwSK21=kwSK20+1)
c$$$      parameter (kwV20=kwSK21+1,kwV11=kwV20+1,kwV02=kwV11+1,
c$$$     $     kwJDZ=kwV02+1,kwJDPZ=kwJDZ+1,kwOFFSET=kwJDPZ+1,
c$$$     $     kwCOUPLE=kwOFFSET+1,kwDDP=kwCOUPLE+1,kwRADI=kwDDP+1,
c$$$     $     kwDPHI=kwRADI+1,kwW1=kwDPHI+1,kwDROT=kwW1+1,
c$$$     $     kwAE1=kwDROT+1,kwAE2=kwAE1+1,
c$$$     $     kwFB1=kwAE2+1,kwFB2=kwFB1+1,
c$$$     $     kwLDEV=kwFB2+1,kwLRAD=kwLDEV+1,
c$$$     $     kwFL=kwLRAD+1,kwAPHI=kwFL+1)
c$$$      parameter (kwMAX=kwAPHI+1)
c$$$c
c$$$      integer kytbl
c$$$      common /kywdtb/ kytbl(0:kwMAX,0:icMXEL)
