      subroutine initb1
      use maccbk
      implicit none
      include 'inc/MACCODE.inc'
      include 'inc/MACKW.inc'
      include 'inc/MACVAR.inc'
      include 'inc/MACPHYS.inc'
      include 'inc/MACFILE.inc'
c
      integer*4 idummy,sethtb,hsrch,mtaloc
      integer*8 ktcaloc
c     external doline
      external doprin, doexpn, doread, dolist, docod, dostop, dotwis
      external dooffl, doonfl,dorvrs
      external ActLie,ActTra,ActPlt,ActGRA

       call defglb('$PLOT$',icGLI,idummy)
       call IsetGL('$PLOT$',0,idummy)
       idummy=sethtb('PLOT   ',icACT,mtaloc(8))
       ilist(1,idval(idummy))=7
       call setfnp(ilist(1,idval(idummy)+1),ActPlt)
       ilist(1,idval(idummy)+2)=hsrch('PTYPE')
       ilist(1,idval(idummy)+3)=hsrch('NPART')
       ilist(1,idval(idummy)+4)=hsrch('TURNS')
       ilist(1,idval(idummy)+5)=hsrch('SPAN    ')
       ilist(1,idval(idummy)+6)=hsrch('CENTER  ')
c
       idummy=sethtb('ID      ',icVAR,VarPt+VarLst)
       idummy=sethtb('GTYPE   ',icVAR,VarStr)
c
c      idummy=sethtb('DRAW    ',icACT,mtaloc(8))
c      idummy=sethtb('graw    ',icACT,idval(idummy))
c      ilist(1,idval(idummy))=7
c      call setfnp(ilist(1,idval(idummy)+1),Act)
c
       idummy=sethtb('GRAPH   ',icACT,mtaloc(8))
       idummy=sethtb('graph   ',icACT,idval(idummy))
       ilist(1,idval(idummy))=7
       call setfnp(ilist(1,idval(idummy)+1),ActGra)
       ilist(1,idval(idummy)+2)=hsrch('ID')
       ilist(1,idval(idummy)+3)=hsrch('GTYPE')
       ilist(1,idval(idummy)+4)=hsrch('TURNS')
       ilist(1,idval(idummy)+5)=hsrch('SPAN    ')
       ilist(1,idval(idummy)+6)=hsrch('CENTER  ')
c
       kytbl(kwL,0)   =sethtb('L       ',icKWRD,kwL   )
       kytbl(kwKIN,0) =sethtb('DISKIN  ',icKWRD,kwKIN )
       kytbl(kwANGL,0)=sethtb('ANGLE   ',icKWRD,kwANGL)
       kytbl(kwROT ,0)=sethtb('ROTATE  ',icKWRD,kwROT)
             idummy   =sethtb('ROT     ',icKWRD,kwROT)
c      kytbl(kwTILT,0)=sethtb('TILT    ',icKWRD,kwROT )
             idummy   =sethtb('TILT    ',icKWRD,kwROT )
       kytbl(kwK0  ,0)=sethtb('K0      ',icKWRD,kwK0  )
       kytbl(kwK1  ,0)=sethtb('K1      ',icKWRD,kwK1  )
c      kytbl(kwDK1 ,0)=sethtb('DK1     ',icKWRD,kwDK1 )
       kytbl(kwK2  ,0)=sethtb('K2      ',icKWRD,kwK2  )
c      kytbl(kwDK2 ,0)=sethtb('DK2     ',icKWRD,kwDK2 )
       kytbl(kwK3  ,0)=sethtb('K3      ',icKWRD,kwK3  )
c      kytbl(kwDK3 ,0)=sethtb('DK3     ',icKWRD,kwDK3 )
       kytbl(kwK4  ,0)=sethtb('K4      ',icKWRD,kwK4  )
       kytbl(kwK5  ,0)=sethtb('K5      ',icKWRD,kwK5  )
       kytbl(kwK6  ,0)=sethtb('K6      ',icKWRD,kwK6  )
       kytbl(kwK7  ,0)=sethtb('K7      ',icKWRD,kwK7  )
       kytbl(kwK8  ,0)=sethtb('K8      ',icKWRD,kwK8  )
       kytbl(kwK9  ,0)=sethtb('K9      ',icKWRD,kwK9  )
       kytbl(kwK10 ,0)=sethtb('K10     ',icKWRD,kwK10  )
       kytbl(kwK10 ,0)=sethtb('K10    ',icKWRD,kwK10 )
       kytbl(kwK11 ,0)=sethtb('K11    ',icKWRD,kwK11 )
       kytbl(kwK12 ,0)=sethtb('K12    ',icKWRD,kwK12 )
       kytbl(kwK13 ,0)=sethtb('K13    ',icKWRD,kwK13 )
       kytbl(kwK14 ,0)=sethtb('K14    ',icKWRD,kwK14 )
       kytbl(kwK15 ,0)=sethtb('K15    ',icKWRD,kwK15 )
       kytbl(kwK16 ,0)=sethtb('K16    ',icKWRD,kwK16 )
       kytbl(kwK17 ,0)=sethtb('K17    ',icKWRD,kwK17 )
       kytbl(kwK18 ,0)=sethtb('K18    ',icKWRD,kwK18 )
       kytbl(kwK19 ,0)=sethtb('K19    ',icKWRD,kwK19 )
       kytbl(kwK20 ,0)=sethtb('K20    ',icKWRD,kwK20 )
       kytbl(kwK21 ,0)=sethtb('K21    ',icKWRD,kwK21 )
c
       kytbl(kwSK0 ,0)=sethtb('SK0     ',icKWRD,kwSK0 )
       kytbl(kwSK1 ,0)=sethtb('SK1     ',icKWRD,kwSK1 )
       kytbl(kwSK2 ,0)=sethtb('SK2     ',icKWRD,kwSK2 )
       kytbl(kwSK3 ,0)=sethtb('SK3     ',icKWRD,kwSK3 )
       kytbl(kwSK4 ,0)=sethtb('SK4     ',icKWRD,kwSK4 )
       kytbl(kwSK5 ,0)=sethtb('SK5     ',icKWRD,kwSK5 )
       kytbl(kwSK6 ,0)=sethtb('SK6     ',icKWRD,kwSK6 )
       kytbl(kwSK7 ,0)=sethtb('SK7     ',icKWRD,kwSK7  )
       kytbl(kwSK8 ,0)=sethtb('SK8     ',icKWRD,kwSK8  )
       kytbl(kwSK9 ,0)=sethtb('SK9     ',icKWRD,kwSK9  )
       kytbl(kwSK10,0)=sethtb('SK10    ',icKWRD,kwSK10 )
       kytbl(kwSK11,0)=sethtb('SK11    ',icKWRD,kwSK11 )
       kytbl(kwSK12,0)=sethtb('SK12    ',icKWRD,kwSK12 )
       kytbl(kwSK13,0)=sethtb('SK13    ',icKWRD,kwSK13 )
       kytbl(kwSK14,0)=sethtb('SK14    ',icKWRD,kwSK14 )
       kytbl(kwSK15 ,0)=sethtb('SK15    ',icKWRD,kwSK15 )
       kytbl(kwSK16 ,0)=sethtb('SK16    ',icKWRD,kwSK16 )
       kytbl(kwSK17 ,0)=sethtb('SK17    ',icKWRD,kwSK17 )
       kytbl(kwSK18 ,0)=sethtb('SK18    ',icKWRD,kwSK18 )
       kytbl(kwSK19 ,0)=sethtb('SK19    ',icKWRD,kwSK19 )
       kytbl(kwSK20 ,0)=sethtb('SK20    ',icKWRD,kwSK20 )
       kytbl(kwSK21 ,0)=sethtb('SK21    ',icKWRD,kwSK21 )
c
       kytbl(kwA3  ,0)=sethtb('A3      ',icKWRD,kwA3  )
       kytbl(kwA5  ,0)=sethtb('A5      ',icKWRD,kwA5  )
       kytbl(kwA7  ,0)=sethtb('A7      ',icKWRD,kwA7  )
       kytbl(kwA9  ,0)=sethtb('A9      ',icKWRD,kwA9  )
       kytbl(kwA11 ,0)=sethtb('A11     ',icKWRD,kwA11 )
       kytbl(kwA13 ,0)=sethtb('A13     ',icKWRD,kwA13 )
       kytbl(kwA15 ,0)=sethtb('A15     ',icKWRD,kwA15 )
       kytbl(kwA17 ,0)=sethtb('A17     ',icKWRD,kwA17 )
       kytbl(kwE1  ,0)=sethtb('E1      ',icKWRD,kwE1  )
       kytbl(kwE2  ,0)=sethtb('E2      ',icKWRD,kwE2  )
       kytbl(kwDX  ,0)=sethtb('DX      ',icKWRD,kwDX  )
       kytbl(kwDY  ,0)=sethtb('DY      ',icKWRD,kwDY  )
       kytbl(kwKICK,0)=sethtb('KICK    ',icKWRD,kwKICK)
       kytbl(kwVOLT,0)=sethtb('VOLT    ',icKWRD,kwVOLT)
       kytbl(kwPHI ,0)=sethtb('PHI     ',icKWRD,kwPHI )
       kytbl(kwDPHI,0)=sethtb('DPHI    ',icKWRD,kwDPHI)
       kytbl(kwFREQ,0)=sethtb('FREQ    ',icKWRD,kwFREQ )
       kytbl(kwHARM,0)=sethtb('HARM    ',icKWRD,kwHARM)
       kytbl(kwLWAK,0)=sethtb('LWAKE   ',icKWRD,kwLWAK)
       kytbl(kwTWAK,0)=sethtb('TWAKE   ',icKWRD,kwTWAK)
       kytbl(kwAX  ,0)=sethtb('AX      ',icKWRD,kwAX  )
                idummy=sethtb('ALPHAX  ',icKWRD,kwAX  )
       kytbl(kwAY  ,0)=sethtb('AY      ',icKWRD,kwAY  )
                idummy=sethtb('ALPHAY  ',icKWRD,kwAY  )
       kytbl(kwBX  ,0)=sethtb('BX      ',icKWRD,kwBX  )
                idummy=sethtb('BETAX   ',icKWRD,kwBX  )
       kytbl(kwBY  ,0)=sethtb('BY      ',icKWRD,kwBY  )
                idummy=sethtb('BETAY   ',icKWRD,kwBY  )
       kytbl(kwEMIX,0)=sethtb('EMITX   ',icKWRD,kwEMIX)
                idummy=sethtb('EMIX    ',icKWRD,kwEMIX)
       kytbl(kwEMIY,0)=sethtb('EMITY   ',icKWRD,kwEMIY)
                idummy=sethtb('EMIY    ',icKWRD,kwEMIY)
       kytbl(kwEMIZ,0)=sethtb('EMITZ   ',icKWRD,kwEMIZ)
                idummy=sethtb('EMIZ    ',icKWRD,kwEMIZ)
       kytbl(kwPX  ,0)=sethtb('PSIX    ',icKWRD,kwPX  )
       kytbl(kwPY  ,0)=sethtb('PSIY    ',icKWRD,kwPY  )
       kytbl(kwPZ  ,0)=sethtb('PSIZ    ',icKWRD,kwPZ  )
       kytbl(kwDP  ,0)=sethtb('DP      ',icKWRD,kwDP  )
       kytbl(kwSIGZ,0)=sethtb('SIGMAZ  ',icKWRD,kwSIGZ)
                idummy=sethtb('SIGZ    ',icKWRD,kwSIGZ)
       kytbl(kwGEO ,0)=sethtb('GEO     ',icKWRD,kwGEO )
       kytbl(kwR1  ,0)=sethtb('R1      ',icKWRD,kwR1  )
       kytbl(kwR2  ,0)=sethtb('R2      ',icKWRD,kwR2  )
       kytbl(kwR3  ,0)=sethtb('R3      ',icKWRD,kwR3  )
       kytbl(kwR4  ,0)=sethtb('R4      ',icKWRD,kwR4  )
       kytbl(kwDETR,0)=sethtb('DETR    ',icKWRD,kwR4  )
       kytbl(kwEX  ,0)=sethtb('EX      ',icKWRD,kwEX  )
       kytbl(kwEPX ,0)=sethtb('EPX     ',icKWRD,kwEPX )
       kytbl(kwEY  ,0)=sethtb('EY      ',icKWRD,kwEY  )
       kytbl(kwEPY ,0)=sethtb('EPY     ',icKWRD,kwEPY )
       kytbl(kwZX  ,0)=sethtb('ZX      ',icKWRD,kwZX  )
       kytbl(kwZPX ,0)=sethtb('ZPX     ',icKWRD,kwZPX )
       kytbl(kwZY  ,0)=sethtb('ZY      ',icKWRD,kwZY  )
       kytbl(kwZPY ,0)=sethtb('ZPY     ',icKWRD,kwZPY )
       kytbl(kwDPX ,0)=sethtb('DPX     ',icKWRD,kwDPX )
       kytbl(kwDPY ,0)=sethtb('DPY     ',icKWRD,kwDPY )
       kytbl(kwAZ  ,0)=sethtb('AZ      ',icKWRD,kwAZ  )
       kytbl(kwD11  ,0)=sethtb('D11      ',icKWRD,kwD11 )
       kytbl(kwD12  ,0)=sethtb('D12      ',icKWRD,kwD12 )
       kytbl(kwD13  ,0)=sethtb('D13      ',icKWRD,kwD13 )
       kytbl(kwD14  ,0)=sethtb('D14      ',icKWRD,kwD14 )
       kytbl(kwD15  ,0)=sethtb('D15      ',icKWRD,kwD15 )
       kytbl(kwD16  ,0)=sethtb('D16      ',icKWRD,kwD16 )
       kytbl(kwD21  ,0)=sethtb('D21      ',icKWRD,kwD21 )
       kytbl(kwD22  ,0)=sethtb('D22      ',icKWRD,kwD22 )
       kytbl(kwD23  ,0)=sethtb('D23      ',icKWRD,kwD23 )
       kytbl(kwD24  ,0)=sethtb('D24      ',icKWRD,kwD24 )
       kytbl(kwD25  ,0)=sethtb('D25      ',icKWRD,kwD25 )
       kytbl(kwD26  ,0)=sethtb('D26      ',icKWRD,kwD26 )
       kytbl(kwD31  ,0)=sethtb('D31      ',icKWRD,kwD31 )
       kytbl(kwD32  ,0)=sethtb('D32      ',icKWRD,kwD32 )
       kytbl(kwD33  ,0)=sethtb('D33      ',icKWRD,kwD33 )
       kytbl(kwD34  ,0)=sethtb('D34      ',icKWRD,kwD34 )
       kytbl(kwD35  ,0)=sethtb('D35      ',icKWRD,kwD35 )
       kytbl(kwD36  ,0)=sethtb('D36      ',icKWRD,kwD36 )
       kytbl(kwD41  ,0)=sethtb('D41      ',icKWRD,kwD41 )
       kytbl(kwD42  ,0)=sethtb('D42      ',icKWRD,kwD42 )
       kytbl(kwD43  ,0)=sethtb('D43      ',icKWRD,kwD43 )
       kytbl(kwD44  ,0)=sethtb('D44      ',icKWRD,kwD44 )
       kytbl(kwD45  ,0)=sethtb('D45      ',icKWRD,kwD45 )
       kytbl(kwD46  ,0)=sethtb('D46      ',icKWRD,kwD46 )
       kytbl(kwD51  ,0)=sethtb('D51      ',icKWRD,kwD51 )
       kytbl(kwD52  ,0)=sethtb('D52      ',icKWRD,kwD52 )
       kytbl(kwD53  ,0)=sethtb('D53      ',icKWRD,kwD53 )
       kytbl(kwD54  ,0)=sethtb('D54      ',icKWRD,kwD54 )
       kytbl(kwD55  ,0)=sethtb('D55      ',icKWRD,kwD55 )
       kytbl(kwD56  ,0)=sethtb('D56      ',icKWRD,kwD56 )
       kytbl(kwD61  ,0)=sethtb('D61      ',icKWRD,kwD61 )
       kytbl(kwD62  ,0)=sethtb('D62      ',icKWRD,kwD62 )
       kytbl(kwD63  ,0)=sethtb('D63      ',icKWRD,kwD63 )
       kytbl(kwD64  ,0)=sethtb('D64      ',icKWRD,kwD64 )
       kytbl(kwD65  ,0)=sethtb('D65      ',icKWRD,kwD65 )
       kytbl(kwD66  ,0)=sethtb('D66      ',icKWRD,kwD66 )
       kytbl(kwB11  ,0)=sethtb('B11      ',icKWRD,kwB11 )
       kytbl(kwB12  ,0)=sethtb('B12      ',icKWRD,kwB12 )
       kytbl(kwB13  ,0)=sethtb('B13      ',icKWRD,kwB13 )
       kytbl(kwB14  ,0)=sethtb('B14      ',icKWRD,kwB14 )
       kytbl(kwB15  ,0)=sethtb('B15      ',icKWRD,kwB15 )
       kytbl(kwB16  ,0)=sethtb('B16      ',icKWRD,kwB16 )
       kytbl(kwB22  ,0)=sethtb('B22      ',icKWRD,kwB22 )
       kytbl(kwB23  ,0)=sethtb('B23      ',icKWRD,kwB23 )
       kytbl(kwB24  ,0)=sethtb('B24      ',icKWRD,kwB24 )
       kytbl(kwB25  ,0)=sethtb('B25      ',icKWRD,kwB25 )
       kytbl(kwB26  ,0)=sethtb('B26      ',icKWRD,kwB26 )
       kytbl(kwB33  ,0)=sethtb('B33      ',icKWRD,kwB33 )
       kytbl(kwB34  ,0)=sethtb('B34      ',icKWRD,kwB34 )
       kytbl(kwB35  ,0)=sethtb('B35      ',icKWRD,kwB35 )
       kytbl(kwB36  ,0)=sethtb('B36      ',icKWRD,kwB36 )
       kytbl(kwB44  ,0)=sethtb('B44      ',icKWRD,kwB44 )
       kytbl(kwB45  ,0)=sethtb('B45      ',icKWRD,kwB45 )
       kytbl(kwB46  ,0)=sethtb('B46      ',icKWRD,kwB46 )
       kytbl(kwB55  ,0)=sethtb('B55      ',icKWRD,kwB55 )
       kytbl(kwB56  ,0)=sethtb('B56      ',icKWRD,kwB56 )
       kytbl(kwB66  ,0)=sethtb('B66      ',icKWRD,kwB66 )
       kytbl(kwR11  ,0)=sethtb('R11      ',icKWRD,kwR11 )
       kytbl(kwR12  ,0)=sethtb('R12      ',icKWRD,kwR12 )
       kytbl(kwR13  ,0)=sethtb('R13      ',icKWRD,kwR13 )
       kytbl(kwR14  ,0)=sethtb('R14      ',icKWRD,kwR14 )
       kytbl(kwR15  ,0)=sethtb('R15      ',icKWRD,kwR15 )
       kytbl(kwR16  ,0)=sethtb('R16      ',icKWRD,kwR16 )
       kytbl(kwR22  ,0)=sethtb('R22      ',icKWRD,kwR22 )
       kytbl(kwR23  ,0)=sethtb('R23      ',icKWRD,kwR23 )
       kytbl(kwR24  ,0)=sethtb('R24      ',icKWRD,kwR24 )
       kytbl(kwR25  ,0)=sethtb('R25      ',icKWRD,kwR25 )
       kytbl(kwR26  ,0)=sethtb('R26      ',icKWRD,kwR26 )
       kytbl(kwR33  ,0)=sethtb('R33      ',icKWRD,kwR33 )
       kytbl(kwR34  ,0)=sethtb('R34      ',icKWRD,kwR34 )
       kytbl(kwR35  ,0)=sethtb('R35      ',icKWRD,kwR35 )
       kytbl(kwR36  ,0)=sethtb('R36      ',icKWRD,kwR36 )
       kytbl(kwR44  ,0)=sethtb('R44      ',icKWRD,kwR44 )
       kytbl(kwR45  ,0)=sethtb('R45      ',icKWRD,kwR45 )
       kytbl(kwR46  ,0)=sethtb('R46      ',icKWRD,kwR46 )
       kytbl(kwR55  ,0)=sethtb('R55      ',icKWRD,kwR55 )
       kytbl(kwR56  ,0)=sethtb('R56      ',icKWRD,kwR56 )
       kytbl(kwR66  ,0)=sethtb('R66      ',icKWRD,kwR66 )
       kytbl(kwRAD ,0)=sethtb('DISRAD  ',icKWRD,kwRAD )
       kytbl(kwCHRO,0)=sethtb('ACHROMA ',icKWRD,kwCHRO)
       kytbl(kwFRIN,0)=sethtb('DISFRIN ',icKWRD,kwFRIN)
       kytbl(kwF1  ,0)=sethtb('F1      ',icKWRD,kwF1  )
       kytbl(kwF2  ,0)=sethtb('F2      ',icKWRD,kwF2  )
       kytbl(kwFRMD,0)=sethtb('FRINGE  ',icKWRD,kwFRMD)
       kytbl(kwK0FR,0)=sethtb('DISK0FR ',icKWRD,kwK0FR)
       kytbl(kwEPS ,0)=sethtb('EPS     ',icKWRD,kwEPS )
       kytbl(kwRANK,0)=sethtb('RANKICK ',icKWRD,kwRANK)
       kytbl(kwRANV,0)=sethtb('RANVOLT ',icKWRD,kwRANV)
       kytbl(kwRANP,0)=sethtb('RANPHASE',icKWRD,kwRANP)
       kytbl(kwDZ  ,0)=sethtb('DZ      ',icKWRD,kwDZ  )
       kytbl(kwDDP ,0)=sethtb('DDP     ',icKWRD,kwDZ  )
       kytbl(kwINDX,0)=sethtb('INDEX   ',icKWRD,kwINDX)
       kytbl(kwBMAX,0)=sethtb('BMAX    ',icKWRD,kwBMAX)
       kytbl(kwBND ,0)=sethtb('BOUND   ',icKWRD,kwBND )
       kytbl(kwPRD ,0)=sethtb('PERIOD  ',icKWRD,kwPRD )
       kytbl(kwBZ  ,0)=sethtb('BZ      ',icKWRD,kwBZ  )
       kytbl(kwDBZ ,0)=sethtb('DBZ     ',icKWRD,kwDBZ )
       kytbl(kwCHI1,0)=sethtb('CHI1    ',icKWRD,kwCHI1)
       kytbl(kwCHI2,0)=sethtb('CHI2    ',icKWRD,kwCHI2)
       kytbl(kwCHI3,0)=sethtb('CHI3    ',icKWRD,kwCHI3)
       kytbl(kwDIR ,0)=sethtb('DIR     ',icKWRD,kwDIR )
       kytbl(kwSLI ,0)=sethtb('SLICE   ',icKWRD,kwSLI )
       kytbl(kwSTURN ,0)=sethtb('STURN   ',icKWRD,kwSTURN )
       kytbl(kwXANGLE ,0)=sethtb('XANGLE  ',icKWRD,kwXANGLE )
       kytbl(kwNP  ,0)=sethtb('NP      ',icKWRD,kwNP  )
       kytbl(kwKx  ,0)=sethtb('KX      ',icKWRD,kwKx  )
       kytbl(kwQy  ,0)=sethtb('QY      ',icKWRD,kwQy  )
       kytbl(kwFBx ,0)=sethtb('B0X     ',icKWRD,kwFBx  )
       kytbl(kwFBy ,0)=sethtb('B0Y     ',icKWRD,kwFBy  )
       kytbl(kwPole,0)=sethtb('POLE    ',icKWRD,kwPole  )
       kytbl(kwJDX ,0)=sethtb('JDX     ',icKWRD,kwJDX  )
       kytbl(kwJDY  ,0)=sethtb('JDY     ',icKWRD,kwJDY  )
       kytbl(kwJDZ  ,0)=sethtb('JDZ     ',icKWRD,kwJDZ  )
       kytbl(kwJDPX ,0)=sethtb('JDPX    ',icKWRD,kwJDPX  )
       kytbl(kwJDPY ,0)=sethtb('JDPY    ',icKWRD,kwJDPY  )
       kytbl(kwJDPZ ,0)=sethtb('JDPZ    ',icKWRD,kwJDPZ  )
       kytbl(kwOFFSET ,0)=sethtb('OFFSET    ',icKWRD,kwOFFSET  )
       kytbl(kwCOUPLE ,0)=sethtb('COUPLE    ',icKWRD,kwCOUPLE  )
       kytbl(kwV1   ,0)=sethtb('V1      ',icKWRD,kwV1  )
       kytbl(kwV20  ,0)=sethtb('V20     ',icKWRD,kwV20  )
       kytbl(kwV11  ,0)=sethtb('V11     ',icKWRD,kwV11  )
       kytbl(kwV02  ,0)=sethtb('V02     ',icKWRD,kwV02  )
       kytbl(kwDX1  ,0)=sethtb('DX1     ',icKWRD,kwDX1  )
       kytbl(kwDX2  ,0)=sethtb('DX2     ',icKWRD,kwDX2  )
       kytbl(kwDY1  ,0)=sethtb('DY1     ',icKWRD,kwDY1  )
       kytbl(kwDY2  ,0)=sethtb('DY2     ',icKWRD,kwDY2  )
       kytbl(kwRADI ,0)=sethtb('RADIUS  ',icKWRD,kwRADI )
       kytbl(kwW1   ,0)=sethtb('W1      ',icKWRD,kwW1 )
       kytbl(kwDROT ,0)=sethtb('DROTATE ',icKWRD,kwDROT)
             idummy    =sethtb('DROT    ',icKWRD,kwDROT)
       kytbl(kwAE1 ,0) =sethtb('AE1     ',icKWRD,kwAE1)
       kytbl(kwAE2 ,0) =sethtb('AE2     ',icKWRD,kwAE2)
       kytbl(kwFB1 ,0) =sethtb('FB1     ',icKWRD,kwFB1)
       kytbl(kwFB2 ,0) =sethtb('FB2     ',icKWRD,kwFB2)
       kytbl(kwLDEV,0) =sethtb('LDEV    ',icKWRD,kwLDEV)
       kytbl(kwLRAD,0) =sethtb('LRAD    ',icKWRD,kwLRAD)
       kytbl(kwFL,0)   =sethtb('FLAT    ',icKWRD,kwLRAD)
       kytbl(kwAPHI,0) =sethtb('AUTOPHI ',icKWRD,kwAPHI)
c  for drift 
      idummy=sethtb('drift   ',icDEF,icDRFT)
      kytbl(0,icDRFT)=sethtb('DRIFT   ',icDEF,icDRFT)
      kytbl(kwL     ,icDRFT)=1
      kytbl(kwKIN   ,icDRFT)=2
      kytbl(kwCOUPLE,icDRFT)=3
      kytbl(kwRADI  ,icDRFT)=4
c
      kytbl(kwAX    ,icDRFT)=5
      kytbl(kwAY    ,icDRFT)=6
      kytbl(kwRAD,   icDRFT)=7
c
      kytbl(kwMAX,   icDRFT)=8
cc for STeering
      idummy=sethtb('st      ',icDEF,icBEND)
      idummy=sethtb('ST      ',icDEF,icBEND)
c                         number of parameters
c  for bend
      idummy=sethtb('bend    ',icDEF,icBEND)
      kytbl(0,icBEND)=sethtb('BEND    ',icDEF,icBEND)
      kytbl(kwL   ,icBEND)=1
      kytbl(kwANGL,icBEND)=2
      kytbl(kwE1  ,icBEND)=3
      kytbl(kwE2  ,icBEND)=4
      kytbl(kwROT ,icBEND)=5
c     kytbl(kwTILT,icBEND)=6
      kytbl(kwRAD ,icBEND)=7
      kytbl(kwK1  ,icBEND)=8
      kytbl(kwDX  ,icBEND)=9
      kytbl(kwDY  ,icBEND)=10
      kytbl(kwK0  ,icBEND)=11
      kytbl(kwFRIN,icBEND)=12
      kytbl(kwEPS ,icBEND)=13
      kytbl(kwRANK,icBEND)=14
      kytbl(kwF1  ,icBEND)=15
      kytbl(kwFRMD,icBEND)=16
      kytbl(kwCOUPLE ,icBEND)=17
      kytbl(kwDROT,icBEND)=18
      kytbl(kwAE1, icBEND)=19
      kytbl(kwAE2, icBEND)=20
      kytbl(kwFB1, icBEND)=21
      kytbl(kwFB2, icBEND)=22
      kytbl(kwINDX,icBEND)=23
c
      kytbl(kwAX  ,icBEND)=24
      kytbl(kwAY  ,icBEND)=25
      kytbl(kwLDEV,icBEND)=26
      kytbl(kwLRAD,icBEND)=27
c
      kytbl(kwMAX ,icBEND)=28
c  for quad
      idummy=sethtb('quad    ',icDEF,icQUAD)
      kytbl(0,icQUAD)=sethtb('QUAD    ',icDEF,icQUAD)
      kytbl(kwL   ,icQUAD)=1
      kytbl(kwK1  ,icQUAD)=2
c     kytbl(kwDK1 ,icQUAD)=3
      kytbl(kwROT ,icQUAD)=4
      kytbl(kwDX  ,icQUAD)=5
      kytbl(kwDY  ,icQUAD)=6
      kytbl(kwRAD     ,icQUAD)=7
      kytbl(kwCHRO    ,icQUAD)=8
      kytbl(kwFRIN    ,icQUAD)=9
      kytbl(kwF1      ,icQUAD)=10
      kytbl(kwF2      ,icQUAD)=11
      kytbl(kwFRMD    ,icQUAD)=12
      kytbl(kwEPS ,icQUAD)=13
      kytbl(kwKIN ,icQUAD)=14
      kytbl(kwCOUPLE     ,icQUAD)=15
      kytbl(kwINDX    ,icQUAD)=16
c
      kytbl(kwAX      ,icQUAD)=17
      kytbl(kwAY      ,icQUAD)=18
      kytbl(kwLDEV    ,icQUAD)=19
      kytbl(kwLRAD    ,icQUAD)=20
c
      kytbl(kwMAX     ,icQUAD)=21
c  for sextu
      idummy=sethtb('sext    ',icDEF,icsext)
      kytbl(0,icSEXT)=sethtb('SEXT    ',icDEF,icsext)
      kytbl(kwL   ,icSEXT)=1
      kytbl(kwK2  ,icSEXT)=2
c     kytbl(kwDK2 ,icSEXT)=3
      kytbl(kwROT ,icSEXT)=4
      kytbl(kwDX  ,icSEXT)=5
      kytbl(kwDY  ,icSEXT)=6
      kytbl(kwRAD     ,icSEXT)=7
      kytbl(kwFRIN    ,icSEXT)=8
      kytbl(kwCOUPLE     ,icSEXT)=9
      kytbl(kwINDX    ,icSEXT)=10
c
      kytbl(kwAX      ,icSEXT)=11
      kytbl(kwAY      ,icSEXT)=12
      kytbl(kwLDEV    ,icSEXT)=13
      kytbl(kwLRAD    ,icSEXT)=14
c
      kytbl(kwMAX     ,icSEXT)=15
c  for octu
      idummy=sethtb('oct     ',icDEF,icoctu)
      kytbl(0,icOCTU)=sethtb('OCT     ',icDEF,icoctu)
      kytbl(kwL   ,icOCTU)=1
      kytbl(kwK3  ,icOCTU)=2
c     kytbl(kwDK3 ,icOCTU)=3
      kytbl(kwROT ,icOCTU)=4
      kytbl(kwDX  ,icOCTU)=5
      kytbl(kwDY  ,icOCTU)=6
      kytbl(kwRAD     ,icOCTU)=7
      kytbl(kwFRIN    ,icOCTU)=8
      kytbl(kwCOUPLE     ,icOCTU)=9
      kytbl(kwINDX    ,icOCTU)=10
c
      kytbl(kwAX      ,icOCTU)=11
      kytbl(kwAY      ,icOCTU)=12
      kytbl(kwLDEV    ,icOCTU)=13
      kytbl(kwLRAD    ,icOCTU)=14
c
      kytbl(kwMAX     ,icOCTU)=15
c  for deca
      idummy=sethtb('deca    ',icDEF,icdeca)
      kytbl(0,icDECA)=sethtb('DECA    ',icDEF,icdeca)
      kytbl(kwL   ,icDECA)=1
      kytbl(kwK4  ,icDECA)=2
c     kytbl(kwDK4 ,icDECA)=3
      kytbl(kwROT ,icDECA)=4
      kytbl(kwDX  ,icDECA)=5
      kytbl(kwDY  ,icDECA)=6
      kytbl(kwRAD     ,icDECA)=7
      kytbl(kwFRIN    ,icDECA)=8
      kytbl(kwCOUPLE     ,icDECA)=9
      kytbl(kwMAX     ,icDECA)=10
c  for dodeca
      idummy=sethtb('dodeca  ',icDEF,icdodeca)
      kytbl(0,icDODECA)=sethtb('DODECA  ',icDEF,icdodeca)
      kytbl(kwL   ,icDODECA)=1
      kytbl(kwK5  ,icDODECA)=2
c     kytbl(kwDK4 ,icDODECA)=3
      kytbl(kwROT ,icDODECA)=4
      kytbl(kwDX  ,icDODECA)=5
      kytbl(kwDY  ,icDODECA)=6
      kytbl(kwRAD     ,icDODECA)=7
      kytbl(kwFRIN    ,icDODECA)=8
      kytbl(kwCOUPLE     ,icDODECA)=9
      kytbl(kwMAX     ,icDODECA)=10
c  for MULT
      idummy=sethtb('mult    ',icDEF,icMULT)
      kytbl(0,icMULT)=sethtb('MULT    ',icDEF,icMULT)
      kytbl(kwL   ,icMULT)=1
      kytbl(kwDX  ,icMULT)=3
      kytbl(kwDY  ,icMULT)=4
      kytbl(kwDZ  ,icMULT)=5
      kytbl(kwCHI1,icMULT)=6
      kytbl(kwCHI2,icMULT)=7
      kytbl(kwCHI3,icMULT)=8
      kytbl(kwROT ,icMULT)=8
      kytbl(kwEPS ,icMULT)=9
      kytbl(kwRAD ,icMULT)=10
      kytbl(kwFRIN,icMULT)=11
      kytbl(kwF1  ,icMULT)=12
      kytbl(kwF2  ,icMULT)=13
      kytbl(kwFRMD,icMULT)=14
      kytbl(kwVOLT,icMULT)=15
      kytbl(kwHARM,icMULT)=16
      kytbl(kwPHI ,icMULT)=17
      kytbl(kwFREQ,icMULT)=18
      kytbl(kwCOUPLE,icMULT)=19
      kytbl(kwRADI,icMULT)=20
      kytbl(kwDPHI,icMULT)=21
      kytbl(kwW1  ,icMULT)=22
c
      kytbl(kwANGL,icMULT)=23
      kytbl(kwE1,  icMULT)=24
      kytbl(kwE2,  icMULT)=25
      kytbl(kwDROT,icMULT)=26
      kytbl(kwK0FR,icMULT)=27
      kytbl(kwFB1,icMULT)=28
      kytbl(kwFB2,icMULT)=29
c
      kytbl(kwK0  ,icMULT)=30
      kytbl(kwSK0 ,icMULT)=31
      kytbl(kwK1  ,icMULT)=32
      kytbl(kwSK1 ,icMULT)=33
      kytbl(kwK2  ,icMULT)=34
      kytbl(kwSK2 ,icMULT)=35
      kytbl(kwK3  ,icMULT)=36
      kytbl(kwSK3 ,icMULT)=37
      kytbl(kwK4  ,icMULT)=38
      kytbl(kwSK4 ,icMULT)=39
      kytbl(kwK5  ,icMULT)=40
      kytbl(kwSK5 ,icMULT)=41
      kytbl(kwK6  ,icMULT)=42
      kytbl(kwSK6 ,icMULT)=43
      kytbl(kwK7  ,icMULT)=44
      kytbl(kwSK7 ,icMULT)=45
      kytbl(kwK8  ,icMULT)=46
      kytbl(kwSK8 ,icMULT)=47
      kytbl(kwK9  ,icMULT)=48
      kytbl(kwSK9 ,icMULT)=49
      kytbl(kwK10 ,icMULT)=50
      kytbl(kwSK10,icMULT)=51
      kytbl(kwK11 ,icMULT)=52
      kytbl(kwSK11,icMULT)=53
      kytbl(kwK12 ,icMULT)=54
      kytbl(kwSK12,icMULT)=55
      kytbl(kwK13 ,icMULT)=56
      kytbl(kwSK13,icMULT)=57
      kytbl(kwK14 ,icMULT)=58
      kytbl(kwSK14,icMULT)=59
      Kytbl(kwK15 ,icMULT)=60
      kytbl(kwSK15,icMULT)=61
      kytbl(kwK16 ,icMULT)=62
      kytbl(kwSK16,icMULT)=63
      kytbl(kwK17 ,icMULT)=64
      kytbl(kwSK17,icMULT)=65
      kytbl(kwK18 ,icMULT)=66
      kytbl(kwSK18,icMULT)=67
      kytbl(kwK19 ,icMULT)=68
      kytbl(kwSK19,icMULT)=69
      kytbl(kwK20 ,icMULT)=70
      kytbl(kwSK20,icMULT)=71
      kytbl(kwK21 ,icMULT)=72
      kytbl(kwSK21,icMULT)=73
      kytbl(kwAE1, icMULT)=74
      kytbl(kwAE2, icMULT)=75
c
      kytbl(kwAX  ,icMULT)=76
      kytbl(kwAY  ,icMULT)=77
      kytbl(kwLDEV,icMULT)=78
      kytbl(kwLRAD,icMULT)=79
      kytbl(kwAPHI,icMULT)=80
c
      kytbl(kwMAX ,icMULT)=81
c  for UNDULATOR
      idummy=sethtb('und     ',icDEF,icUND)
      kytbl(0,icUND)=sethtb('UND    ',icDEF,icUND)
      kytbl(kwL   ,icUND)=1
      kytbl(kwFBX ,icUND)=2
      kytbl(kwFBY ,icUND)=3
      kytbl(kwKx  ,icUND)=4
      kytbl(kwQy  ,icUND)=5
      kytbl(kwPHI ,icUND)=6
      kytbl(kwSLI ,icUND)=7
      kytbl(kwPole,icUND)=8
      kytbl(kwCOUPLE ,icUND)=9
      kytbl(kwMAX ,icUND)=10
c  for WIG
      idummy=sethtb('wig     ',icDEF,icWIG)
      kytbl(0,icWIG)=sethtb('WIG    ',icDEF,icWIG)
      kytbl(kwL   ,icWIG)=1
      kytbl(kwBMAX,icWIG)=2
      kytbl(kwPRD ,icWIG)=3
      kytbl(kwROT ,icWIG)=4
      kytbl(kwDX  ,icWIG)=5
      kytbl(kwDY  ,icWIG)=6
      kytbl(kwA3  ,icWIG)=7
      kytbl(kwA5  ,icWIG)=8
      kytbl(kwA7  ,icWIG)=9
      kytbl(kwA9  ,icWIG)=10
      kytbl(kwA11 ,icWIG)=11
      kytbl(kwA13 ,icWIG)=12
      kytbl(kwA15 ,icWIG)=13
      kytbl(kwA17 ,icWIG)=14
      kytbl(kwCOUPLE ,icWIG)=15
      kytbl(kwMAX ,icWIG)=16
c  for solenoid
      idummy=sethtb('sol     ',icDEF,icSOL)
      kytbl(0,icSOL)=sethtb('SOL     ',icDEF,icSOL)
      kytbl(kwL   ,icSOL)=1
      kytbl(kwBZ  ,icSOL)=2
      kytbl(kwDX  ,icSOL)=3
      kytbl(kwDY  ,icSOL)=4
      kytbl(kwDZ  ,icSOL)=5
      kytbl(kwDPX ,icSOL)=6
      kytbl(kwDPY ,icSOL)=7
      kytbl(kwBND ,icSOL)=8
      kytbl(kwCHI1,icSOL)=9
      kytbl(kwCHI2,icSOL)=10
      kytbl(kwCHI3,icSOL)=11
      kytbl(kwGEO ,icSOL)=12
      kytbl(kwF1  ,icSOL)=13
      kytbl(kwDBZ ,icSOL)=14
      kytbl(kwCOUPLE ,icSOL)=18
      kytbl(kwFRIN,icSOL)=19
      kytbl(kwRAD, icSOL)=20
      kytbl(kwFL, icSOL)=21
      kytbl(kwMAX  ,icSOL)=22
cc for TEST
      idummy=sethtb('test    ',icDEF,icTEST)
      kytbl(0,icTEST)=sethtb('TEST    ',icDEF,icTEST)
      kytbl(kwL   ,icTEST)=1
      kytbl(kwANGL,icTEST)=2
      kytbl(kwMAX ,icTEST)=3
cc for Cavity
      idummy=sethtb('cavi    ',icDEF,icCavi)
      kytbl(0,icCAVI)=sethtb('CAVI    ',icDEF,icCavi)
      kytbl(kwL   ,icCAVI)=1
      kytbl(kwVOLT,icCAVI)=2
      kytbl(kwHARM,icCAVI)=3
      kytbl(kwPHI,icCAVI)=4
      kytbl(kwFREQ,icCAVI)=5
      kytbl(kwRANV,icCAVI)=9
      kytbl(kwRANP,icCAVI)=10
      kytbl(kwLWAK,icCAVI)=11
      kytbl(kwTWAK,icCAVI)=12
      kytbl(kwDX  ,icCAVI)=13
      kytbl(kwDY  ,icCAVI)=14
      kytbl(kwROT ,icCAVI)=15
      kytbl(kwV1  ,icCAVI)=16
      kytbl(kwV20  ,icCAVI)=17
      kytbl(kwV11  ,icCAVI)=18
      kytbl(kwV02  ,icCAVI)=19
      kytbl(kwCOUPLE ,icCAVI)=20
      kytbl(kwDPHI,icCAVI)=21
      kytbl(kwFRIN,icCAVI)=22
      kytbl(kwFRMD,icCAVI)=23
c
      kytbl(kwAX  ,icCAVI)=24
      kytbl(kwAY  ,icCAVI)=25
      kytbl(kwLDEV,icCAVI)=26
      kytbl(kwAPHI,icCAVI)=27
      kytbl(kwMAX ,icCAVI)=28
cc for t-Cavity
      idummy=sethtb('tcavi   ',icDEF,icTCAV)
      kytbl(0,icTCAV)=sethtb('TCAVI   ',icDEF,icTCAV)
      kytbl(kwL   ,icTCAV)=1
      kytbl(kwK0  ,icTCAV)=2
      kytbl(kwHARM,icTCAV)=3
      kytbl(kwPHI ,icTCAV)=4
      kytbl(kwFREQ,icTCAV)=5
      kytbl(kwDX  ,icTCAV)=6
      kytbl(kwDY  ,icTCAV)=7
      kytbl(kwROT ,icTCAV)=8
      kytbl(kwRANK,icTCAV)=9
      kytbl(kwRANP,icTCAV)=10
      kytbl(kwLWAK,icTCAV)=11
      kytbl(kwTWAK,icTCAV)=12
      kytbl(kwCOUPLE     ,icTCAV)=13
c
      kytbl(kwAX  ,icTCAV)=14
      kytbl(kwAY  ,icTCAV)=15
      kytbl(kwLDEV,icTCAV)=16
c     kytbl(kwLRAD,icTCAV)=17
c
      kytbl(kwMAX ,icTCAV)=18
cc for MAP
      idummy=sethtb('map    ',icDEF,icMAP)
      kytbl(0,icMAP)=sethtb('MAP    ',icDEF,icMAP)
      kytbl(kwL       ,icMAP)=1
      kytbl(kwMAX     ,icMAP)=2
cc for INS
      idummy=sethtb('ins    ',icDEF,icINS)
      kytbl(0,icINS)=sethtb('INS    ',icDEF,icINS)
      kytbl(kwAX      ,icINS)=1
      kytbl(kwBX      ,icINS)=2
      kytbl(kwPX      ,icINS)=3
      kytbl(kwAY      ,icINS)=4
      kytbl(kwBY      ,icINS)=5
      kytbl(kwPY      ,icINS)=6
      kytbl(kwEX      ,icINS)=7
      kytbl(kwEPX     ,icINS)=8
      kytbl(kwEY      ,icINS)=9
      kytbl(kwEPY     ,icINS)=10
      kytbl(kwR1      ,icINS)=11
      kytbl(kwR2      ,icINS)=12
      kytbl(kwR3      ,icINS)=13
      kytbl(kwR4      ,icINS)=14
      kytbl(kwDX      ,icINS)=15
      kytbl(kwDPX     ,icINS)=16
      kytbl(kwDY      ,icINS)=17
      kytbl(kwDPY     ,icINS)=18
      kytbl(kwDIR     ,icINS)=19
      kytbl(kwCOUPLE     ,icINS)=62
      kytbl(kwMAX     ,icINS)=63
cc for Coordinate transformation
      idummy=sethtb('COORD   ',icDEF,icCOORD)
      kytbl(0,     icCOORD  )=sethtb('COORD   ',icDEF,icCOORD)
      kytbl(kwDX  ,icCOORD  )=1
      kytbl(kwDY  ,icCOORD  )=2
      kytbl(kwDZ  ,icCOORD  )=3
      kytbl(kwCHI1,icCOORD  )=4
      kytbl(kwCHI2,icCOORD  )=5
      kytbl(kwCHI3,icCOORD  )=6
      kytbl(kwDIR ,icCOORD  )=7
      kytbl(kwCOUPLE ,icCOORD  )=8
      kytbl(kwMAX ,icCOORD  )=9
cc for BEAMBEAM
      idummy=sethtb('beambeam',icDEF,icBEAM)
      kytbl(0,icBEAM)=sethtb('BEAMBEAM',icDEF,icBEAM)
      kytbl(kwAX      ,icBEAM)=1
      kytbl(kwBX      ,icBEAM)=2
      kytbl(kwAY      ,icBEAM)=3
      kytbl(kwBY      ,icBEAM)=4
      kytbl(kwR1      ,icBEAM)=5
      kytbl(kwR2      ,icBEAM)=6
      kytbl(kwR3      ,icBEAM)=7
      kytbl(kwR4      ,icBEAM)=8
      kytbl(kwEX      ,icBEAM)=9
      kytbl(kwEPX     ,icBEAM)=10
      kytbl(kwEY      ,icBEAM)=11
      kytbl(kwEPY     ,icBEAM)=12
      kytbl(kwZX      ,icBEAM)=13
      kytbl(kwZPX     ,icBEAM)=14
      kytbl(kwZY      ,icBEAM)=15
      kytbl(kwZPY     ,icBEAM)=16
      kytbl(kwDX      ,icBEAM)=17
      kytbl(kwDPX     ,icBEAM)=18
      kytbl(kwDY      ,icBEAM)=19
      kytbl(kwDPY     ,icBEAM)=20
      kytbl(kwXANGLE  ,icBEAM)=21
      kytbl(kwEMIX    ,icBEAM)=22
      kytbl(kwEMIY    ,icBEAM)=23
      kytbl(kwDP      ,icBEAM)=24
      kytbl(kwAZ      ,icBEAM)=25
      kytbl(kwDZ      ,icBEAM)=26
      kytbl(kwSIGZ    ,icBEAM)=27
      kytbl(kwSLI     ,icBEAM)=28
      kytbl(kwNP      ,icBEAM)=29
      kytbl(kwSTURN   ,icBEAM)=30
      kytbl(kwR11     ,icBEAM)=31
      kytbl(kwR12     ,icBEAM)=32
      kytbl(kwR13     ,icBEAM)=33
      kytbl(kwR14     ,icBEAM)=34
      kytbl(kwR15     ,icBEAM)=35
      kytbl(kwR16     ,icBEAM)=36
      kytbl(kwR22     ,icBEAM)=37
      kytbl(kwR23     ,icBEAM)=38
      kytbl(kwR24     ,icBEAM)=39
      kytbl(kwR25     ,icBEAM)=40
      kytbl(kwR26     ,icBEAM)=41
      kytbl(kwR33     ,icBEAM)=42
      kytbl(kwR34     ,icBEAM)=43
      kytbl(kwR35     ,icBEAM)=44
      kytbl(kwR36     ,icBEAM)=45
      kytbl(kwR44     ,icBEAM)=46
      kytbl(kwR45     ,icBEAM)=47
      kytbl(kwR46     ,icBEAM)=48
      kytbl(kwR55     ,icBEAM)=49
      kytbl(kwR56     ,icBEAM)=50
      kytbl(kwR66     ,icBEAM)=51
      kytbl(kwCOUPLE     ,icBEAM)=52
      kytbl(kwMAX     ,icBEAM)=53
cc for PHSROT
      idummy=sethtb('phsrot  ',icDEF,icProt)
      kytbl(0,icProt)=sethtb('PHSROT  ',icDEF,icProt)
      kytbl(kwAX      ,icProt)=1
      kytbl(kwBX      ,icProt)=2
      kytbl(kwPX      ,icProt)=3
      kytbl(kwAY      ,icProt)=4
      kytbl(kwBY      ,icProt)=5
      kytbl(kwPY      ,icProt)=6
      kytbl(kwR1      ,icProt)=7
      kytbl(kwR2      ,icProt)=8
      kytbl(kwR3      ,icProt)=9
      kytbl(kwR4      ,icProt)=10
      kytbl(kwEX      ,icProt)=11
      kytbl(kwEPX     ,icProt)=12
      kytbl(kwEY      ,icProt)=13
      kytbl(kwEPY     ,icProt)=14
      kytbl(kwZX      ,icProt)=15
      kytbl(kwZPX     ,icProt)=16
      kytbl(kwZY      ,icProt)=17
      kytbl(kwZPY     ,icProt)=18
      kytbl(kwEMIX    ,icProt)=19
      kytbl(kwEMIY    ,icProt)=20
      kytbl(kwDP      ,icProt)=21
      kytbl(kwAZ      ,icProt)=22
      kytbl(kwDZ      ,icProt)=23
      kytbl(kwSIGZ    ,icProt)=24
      kytbl(kwPZ      ,icProt)=25
      kytbl(kwEMIZ    ,icProt)=26
      kytbl(kwBZ      ,icProt)=27
      kytbl(kwJDY     ,icProt)=28
      kytbl(kwJDPY    ,icProt)=29
      kytbl(kwD11     ,icProt)=30
      kytbl(kwD12     ,icProt)=31
      kytbl(kwD13     ,icProt)=32
      kytbl(kwD14     ,icProt)=33
      kytbl(kwD15     ,icProt)=34
      kytbl(kwD16     ,icProt)=35
      kytbl(kwD21     ,icProt)=36
      kytbl(kwD22     ,icProt)=37
      kytbl(kwD23     ,icProt)=38
      kytbl(kwD24     ,icProt)=39
      kytbl(kwD25     ,icProt)=40
      kytbl(kwD26     ,icProt)=41
      kytbl(kwD31     ,icProt)=42
      kytbl(kwD32     ,icProt)=43
      kytbl(kwD33     ,icProt)=44
      kytbl(kwD34     ,icProt)=45
      kytbl(kwD35     ,icProt)=46
      kytbl(kwD36     ,icProt)=47
      kytbl(kwD41     ,icProt)=48
      kytbl(kwD42     ,icProt)=49
      kytbl(kwD43     ,icProt)=50
      kytbl(kwD44     ,icProt)=51
      kytbl(kwD45     ,icProt)=52
      kytbl(kwD46     ,icProt)=53
      kytbl(kwD51     ,icProt)=54
      kytbl(kwD52     ,icProt)=55
      kytbl(kwD53     ,icProt)=56
      kytbl(kwD54     ,icProt)=57
      kytbl(kwD55     ,icProt)=58
      kytbl(kwD56     ,icProt)=59
      kytbl(kwD61     ,icProt)=60
      kytbl(kwD62     ,icProt)=61
      kytbl(kwD63     ,icProt)=62
      kytbl(kwD64     ,icProt)=63
      kytbl(kwD65     ,icProt)=64
      kytbl(kwD66     ,icProt)=65
      kytbl(kwB11     ,icProt)=66
      kytbl(kwB12     ,icProt)=67
      kytbl(kwB13     ,icProt)=68
      kytbl(kwB14     ,icProt)=69
      kytbl(kwB15     ,icProt)=70
      kytbl(kwB16     ,icProt)=71
      kytbl(kwB22     ,icProt)=72
      kytbl(kwB23     ,icProt)=73
      kytbl(kwB24     ,icProt)=74
      kytbl(kwB25     ,icProt)=75
      kytbl(kwB26     ,icProt)=76
      kytbl(kwB33     ,icProt)=77
      kytbl(kwB34     ,icProt)=78
      kytbl(kwB35     ,icProt)=79
      kytbl(kwB36     ,icProt)=80
      kytbl(kwB44     ,icProt)=81
      kytbl(kwB45     ,icProt)=82
      kytbl(kwB46     ,icProt)=83
      kytbl(kwB55     ,icProt)=84
      kytbl(kwB56     ,icProt)=85
      kytbl(kwB66     ,icProt)=86
      kytbl(kwCOUPLE     ,icProt)=87
      kytbl(kwMAX     ,icProt)=88
cc for MARK
      idummy=sethtb('mark    ',icDEF,icMark)
      kytbl(0,icMARK)=sethtb('MARK    ',icDEF,icMark)
      kytbl(kwAX      ,icMark)=1
      kytbl(kwBX      ,icMark)=2
      kytbl(kwPX      ,icMark)=3
      kytbl(kwAY      ,icMark)=4
      kytbl(kwBY      ,icMark)=5
      kytbl(kwPY      ,icMark)=6
      kytbl(kwEX      ,icMark)=7
      kytbl(kwEPX     ,icMark)=8
      kytbl(kwEY      ,icMark)=9
      kytbl(kwEPY     ,icMark)=10
      kytbl(kwR1      ,icMark)=11
      kytbl(kwR2      ,icMark)=12
      kytbl(kwR3      ,icMark)=13
      kytbl(kwR4      ,icMark)=14
      kytbl(kwDETR    ,icMark)=15
      kytbl(kwDX      ,icMark)=16
      kytbl(kwDPX     ,icMark)=17
      kytbl(kwDY      ,icMark)=18
      kytbl(kwDPY     ,icMark)=19
      kytbl(kwDZ      ,icMark)=20
      kytbl(kwDDP     ,icMark)=21
      kytbl(kwDP      ,icMark)=22
      kytbl(kwAZ      ,icMark)=23
      kytbl(kwOFFSET  ,icMark)=24
      kytbl(kwSIGZ    ,icMark)=25
      kytbl(kwGEO     ,icMark)=26
      kytbl(kwJDX     ,icMark)=27
      kytbl(kwJDPX    ,icMark)=28
      kytbl(kwJDY     ,icMark)=29
      kytbl(kwJDPY    ,icMark)=30
      kytbl(kwJDZ     ,icMark)=31
      kytbl(kwJDPZ    ,icMark)=32
      kytbl(kwEMIX    ,icMark)=33
      kytbl(kwEMIY    ,icMark)=34
      kytbl(kwCOUPLE  ,icMark)=35
      kytbl(kwMAX     ,icMark)=36
cc for apert
      idummy=sethtb('apert   ',icDEF,icAprt)
      kytbl(0,icAprt)=sethtb('APERT   ',icDEF,icAprt)
      kytbl(kwDX1     ,icAprt)=1
      kytbl(kwDX2     ,icAprt)=2
      kytbl(kwDY1     ,icAprt)=3
      kytbl(kwDY2     ,icAprt)=4
      kytbl(kwJDPX    ,icAprt)=5
      kytbl(kwJDPY    ,icAprt)=6
      kytbl(kwDP      ,icAprt)=7
      kytbl(kwCOUPLE  ,icAprt)=8
      kytbl(kwROT     ,icAprt)=9
      kytbl(kwAX      ,icAprt)=10
      kytbl(kwAY      ,icAprt)=11
      kytbl(kwDX      ,icAprt)=12
      kytbl(kwDY      ,icAprt)=13      
      kytbl(kwMAX     ,icAprt)=14
cc for mon
      idummy=sethtb('moni    ',icDEF,icMoni)
      kytbl(0,icMONI)=sethtb('MONI    ',icDEF,icMoni)
      kytbl(kwDX      ,icMONI)=1
      kytbl(kwDY      ,icMONI)=2
      kytbl(kwOFFSET  ,icMONI)=3
      kytbl(kwCOUPLE     ,icMONI)=4
      kytbl(kwROT     ,icMONI)=5
      kytbl(kwMAX     ,icMONI)=6
cc for spch
      idummy=sethtb('spch    ',icDEF,icSPCH)
      kytbl(0,icSPCH)=sethtb('SPCH    ',icDEF,icSPCH)
      kytbl(kwAX      ,icSPCH)=1
      kytbl(kwBX      ,icSPCH)=2
      kytbl(kwPX      ,icSPCH)=3
      kytbl(kwAY      ,icSPCH)=4
      kytbl(kwBY      ,icSPCH)=5
      kytbl(kwPY      ,icSPCH)=6
      kytbl(kwR1      ,icSPCH)=7
      kytbl(kwR2      ,icSPCH)=8
      kytbl(kwR3      ,icSPCH)=9
      kytbl(kwR4      ,icSPCH)=10
      kytbl(kwEX      ,icSPCH)=11
      kytbl(kwEPX     ,icSPCH)=12
      kytbl(kwEY      ,icSPCH)=13
      kytbl(kwEPY     ,icSPCH)=14
      kytbl(kwZX      ,icSPCH)=15
      kytbl(kwZPX     ,icSPCH)=16
      kytbl(kwZY      ,icSPCH)=17
      kytbl(kwZPY     ,icSPCH)=18
      kytbl(kwMAX     ,icSPCH)=19
c.....for debug
c     do 9999 i=1,kwMAX
c       print *,i,kytbl(i,0),pname(kytbl(i,0))
c9999 continue
c     do 9998 i=1,icQUAD
c       do 9998 j=0,kwMAX
c       print *,j,i,kytbl(j,i),pname(kytbl(0,i)),pname(kytbl(j,0))
c9998 continue
c.....end debug
c
       rlist(idval(sethtb8('%       ',icUNIT,ktcaloc(3))))=0.01d0
       rlist(idval(sethtb8('rad     ',icUNIT,ktcaloc(3))))=1.00d0
       rlist(idval(sethtb8('RAD     ',icUNIT,ktcaloc(3))))=1.00d0
       rlist(idval(sethtb8('mrad    ',icUNIT,ktcaloc(3))))=1.00d-3
       rlist(idval(sethtb8('MRAD    ',icUNIT,ktcaloc(3))))=1.00d-3
       rlist(idval(sethtb8('DEG     ',icUNIT,ktcaloc(3))))=pi/180.d0
       rlist(idval(sethtb8('deg     ',icUNIT,ktcaloc(3))))=pi/180.d0
       rlist(idval(sethtb8('M       ',icUNIT,ktcaloc(3))))=1.00d0
       rlist(idval(sethtb8('m       ',icUNIT,ktcaloc(3))))=1.00d0
       rlist(idval(sethtb8('cm      ',icUNIT,ktcaloc(3))))=1.00d-2
       rlist(idval(sethtb8('CM      ',icUNIT,ktcaloc(3))))=1.0d-2
       rlist(idval(sethtb8('mm      ',icUNIT,ktcaloc(3))))=1.0d-3
       rlist(idval(sethtb8('MM      ',icUNIT,ktcaloc(3))))=1.0d-3
       rlist(idval(sethtb8('T       ',icUNIT,ktcaloc(3))))=1.0d0
       rlist(idval(sethtb8('GAUSS   ',icUNIT,ktcaloc(3))))=1.0d-4
       rlist(idval(sethtb8('gauss   ',icUNIT,ktcaloc(3))))=1.0d-4
c      Standard energy unit wasJoule.(MKSA unit system)
c      Standard energy unit is Now eV.
       rlist(idval(sethtb8('JOULE   ',icUNIT,ktcaloc(3))))
     &                 = 1.0d0/elemch
       rlist(idval(sethtb8('eV      ',icUNIT,ktcaloc(3))))
     &                 = 1.0d0
       rlist(idval(sethtb8('KeV     ',icUNIT,ktcaloc(3))))
     &                 = 1.d3
       rlist(idval(sethtb8('MeV     ',icUNIT,ktcaloc(3))))
     &                 = 1.d6
       rlist(idval(sethtb8('GeV     ',icUNIT,ktcaloc(3))))
     &                 = 1.d9
       rlist(idval(sethtb8('EV      ',icUNIT,ktcaloc(3))))
     &                 = 1.0d0
       rlist(idval(sethtb8('KEV     ',icUNIT,ktcaloc(3))))
     &                 = 1.d3
       rlist(idval(sethtb8('MEV     ',icUNIT,ktcaloc(3))))
     &                 = 1.d6
       rlist(idval(sethtb8('GEV     ',icUNIT,ktcaloc(3))))
     &                 = 1.d9
c
       rlist(idval(sethtb8('V       ',icUNIT,ktcaloc(3))))
     &                 = 1.d0
       rlist(idval(sethtb8('KV      ',icUNIT,ktcaloc(3))))
     &                 = 1.d3
       rlist(idval(sethtb8('MV      ',icUNIT,ktcaloc(3))))
     &                 = 1.d6
       rlist(idval(sethtb8('GV      ',icUNIT,ktcaloc(3))))
     &                 = 1.d9
c
       rlist(idval(sethtb8('HZ      ',icUNIT,ktcaloc(3))))
     &                 = 1.d0
       rlist(idval(sethtb8('KHZ     ',icUNIT,ktcaloc(3))))
     &                 = 1.d3
       rlist(idval(sethtb8('MHZ     ',icUNIT,ktcaloc(3))))
     &                 = 1.d6
       rlist(idval(sethtb8('GHZ     ',icUNIT,ktcaloc(3))))
     &                 = 1.d9
c
       idummy=sethtb('normal  ',icRAND,1 )
       idummy=sethtb('NORMAL  ',icRAND,1 )
       idummy=sethtb('uniform ',icRAND,2 )
       idummy=sethtb('UNIFORM ',icRAND,2 )
c
       call defglb('STACKSIZ',icGLI,idummy)
       call IsetGL('STACKSIZ',2**21,idummy)
       call defglb('$SEED',icGLR,idummy)
       call RsetGL('$SEED',17.d0,idummy)
       call defglb('SEED',icGLR,idummy)
       call RsetGL('SEED',17.d0,idummy)
       call defglb('PI',icGLR,idummy)
       call RsetGL('PI',asin(1.d0)*2.d0,idummy)
       call defglb('MOMENTUM',icGLR,idummy)
       call RsetGL('MOMENTUM',30.d9,idummy)
       call defglb('PBUNCH',icGLR,idummy)
       call RsetGL('PBUNCH',1.d10,idummy)
       call defglb('NBUNCH',icGLR,idummy)
       call RsetGL('NBUNCH',1.d0,idummy)
       call defglb('EMITDIV',icGLR,idummy)
       call RsetGL('EMITDIV',1.d0,idummy)
       call defglb('EMITDIVB',icGLR,idummy)
       call RsetGL('EMITDIVB',1.d0,idummy)
       call defglb('EMITDIVQ',icGLR,idummy)
       call RsetGL('EMITDIVQ',1.d0,idummy)
       call defglb('EMITDIVS',icGLR,idummy)
       call RsetGL('EMITDIVS',1.d0,idummy)
       call defglb('FRINGDIV',icGLR,idummy)
       call RsetGL('FRINGDIV',1.d0,idummy)
       call defglb('MINCOUP',icGLR,idummy)
       call RsetGL('MINCOUP',1.d-2,idummy)
       call defglb('LOSSAMPL',icGLR,idummy)
       call RsetGL('LOSSAMPL',1.d0,idummy)
       call defglb('LOSSDZ',icGLR,idummy)
       call RsetGL('LOSSDZ',100.d0,idummy)
       call defglb('BBCUT',icGLR,idummy)
       call RsetGL('BBCUT',0.d0,idummy)
       call defglb('EMITX',icGLR,idummy)
       call RsetGL('EMITX',0.0d0,idummy)
       call defglb('EMITY',icGLR,idummy)
       call RsetGL('EMITY',0.0d0,idummy)
       call defglb('EMITZ',icGLR,idummy)
       call RsetGL('EMITZ',0.0d0,idummy)
       call defglb('EMITXE',icGLR,idummy)
       call RsetGL('EMITXE',0.0d0,idummy)
       call defglb('EMITYE',icGLR,idummy)
       call RsetGL('EMITYE',0.0d0,idummy)
       call defglb('EMITZE',icGLR,idummy)
       call RsetGL('EMITZE',0.0d0,idummy)
       call defglb('SIGE',icGLR,idummy)
       call RsetGL('SIGE',0.0d0,idummy)
       call defglb('SIGZ',icGLR,idummy)
       call RsetGL('SIGZ',0.0d0,idummy)
       call defglb('GCUT',icGLR,idummy)
       call RsetGL('GCUT',1.d35,idummy)
       call defglb('TDXI',icGLR,idummy)
       call RsetGL('TDXI',0.0d0,idummy)
       call defglb('TDYI',icGLR,idummy)
       call RsetGL('TDYI',0.0d0,idummy)
       call defglb('TDZI',icGLR,idummy)
       call RsetGL('TDZI',0.0d0,idummy)
       call defglb('$MASS$',icGLR,idummy)
       call RsetGL('$MASS$',elmass,idummy)
       call defglb('MASS',icGLR,idummy)
       call RsetGL('MASS',elmass,idummy)
       call defglb('NP',icGLR,idummy)
       call RsetGL('NP',0.d0,idummy)
       call defglb('TURNS',icGLR,idummy)
       call RsetGL('TURNS',0.d0,idummy)
       call defglb('CHARGE',icGLR,idummy)
       call RsetGL('CHARGE',0.d0,idummy)
c
       call defglb('OMEGA0',icGLR,idummy)
       call RsetGL('OMEGA0',0.0d0,idummy)
       call defglb('DTSYNCH',icGLR,idummy)
       call RsetGL('DTSYNCH',0.0d0,idummy)
       call defglb('PHICAV',icGLR,idummy)
       call RsetGL('PHICAV',0.0d0,idummy)
       call defglb('EFFVCRATIO',icGLR,idummy)
       call RsetGL('EFFVCRATIO',1.0d0,idummy)
       call defglb('HRF',icGLR,idummy)
       call RsetGL('HRF',5120.0D0,idummy)
       call defglb('FSHIFT',icGLR,idummy)
       call RsetGL('FSHIFT',0.0d0,idummy)
       call defglb('PSPAN',icGLR,idummy)
       call RsetGL('PSPAN',0.0d0,idummy)
       call defglb('ESHIFT',icGLR,idummy)
       call RsetGL('ESHIFT',0.0d0,idummy)
       call defglb('PHIS',icGLR,idummy)
       call RsetGL('PHIS',0.D00,idummy)
       call defglb('MAXORDER',icGLR,idummy)
       call RsetGL('MAXORDER',6.D00,idummy)
       call defglb('ADDDENSE',icGLR,idummy)
       call RsetGL('ADDDENSE',6.D00,idummy)
       call defglb('ADDTERMS',icGLR,idummy)
       call RsetGL('ADDTERMS',0.D00,idummy)
       call defglb('DJPLOT',icGLR,idummy)
       call RsetGL('DJPLOT',0.D00,idummy)
       call defglb('PHSPLOTS',icGLR,idummy)
       call RsetGL('PHSPLOTS',0.D00,idummy)
       call defglb('DAPWIDTH',icGLR,idummy)
       call RsetGL('DAPWIDTH',7.D00,idummy)
       call defglb('RADIAT',icGLI,idummy)
       call IsetGL('RADIAT',-1,idummy)
       call defglb('PSPACNX',icGLR,idummy)
       call RsetGL('PSPACNX',128.0D0,idummy)
       call defglb('PSPACNY',icGLR,idummy)
       call RsetGL('PSPACNY',128.0D0,idummy)
       call defglb('PSPACNZ',icGLR,idummy)
       call RsetGL('PSPACNZ',1.0D0,idummy)
       call defglb('PSPACDX',icGLR,idummy)
       call RsetGL('PSPACDX',1.0d-3,idummy)
       call defglb('PSPACDY',icGLR,idummy)
       call RsetGL('PSPACDY',1.0d-3,idummy)
       call defglb('PSPACDZ',icGLR,idummy)
       call RsetGL('PSPACDZ',1.0d-3,idummy)
       call defglb('PSPACNTURN',icGLR,idummy)
       call RsetGL('PSPACNTURN',1.0D0,idummy)
       call defglb('PSPACNTURNCALC',icGLR,idummy)
       call RsetGL('PSPACNTURNCALC',0.0D0,idummy)
c
c      call defglb('$RADIZZZ',icGLI,idummy)
c      call IsetGL('$RADIZZZ',-1,idummy)
c
       call defflg('RAD',FLAGOF)
       call defflg('RADCOD',FLAGOF)
       call defflg('PHOTONS',FLAGOF)
       call defflg('EMIT',FLAGON)
       call defflg('INTRA',FLAGOF)
       call defflg('POL'  ,FLAGOF)
       call defflg('COD',FLAGON)
       call defflg('RFSW',FLAGOF)
       call defflg('EMIOUT',FLAGOF)
       call defflg('DAPERT',FLAGOF)
       call defflg('FLUC',FLAGON)
       call defflg('CMPLOT',FLAGOF)
       call defflg('FOURIE',FLAGOF)
       call defflg('SMEAR',FLAGON)
       call defflg('GEOCAL',FLAGON)
       call defflg('DEBUG',FLAGOF)
       call defflg('ECHO',FLAGOF)
       call defflg('LOG',FLAGOF)
       call defflg('CTIME',FLAGOF)
       call defflg('INTRES',FLAGON)
       call defflg('HALFRES',FLAGON)
       call defflg('SUMRES',FLAGON)
       call defflg('DIFFRES',FLAGOF)
       call defglb('NPARA',icGLR,idummy)
       call RsetGL('NPARA',1.d0,idummy)
       msglvl=0
       return
       end
