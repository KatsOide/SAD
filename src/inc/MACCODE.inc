c     Don't confuse, Emacs. This is -*- fortran -*- mode!
c     internal code list
c$$$      integer icDRFT,icBEND,icQUAD,icSEXT,icOCTU,icDECA,icDODECA
c$$$      parameter(icDRFT=1,icBEND=2,icQUAD=4,icSEXT=6,
c$$$     $     icOCTU=8,icDECA=10,icDODECA=12)
c$$$      integer icSOL,icWIG,icUND
c$$$      parameter (icSOL=20,icWIG=19,icUND=18)
c$$$      integer icST,icMULT,icTEST
c$$$      parameter(icST=21,icMULT=22,icTEST=30)
c$$$      integer icCAVI,icTCAV,icMARK,icMONI,icAPRT
c$$$      parameter(icCAVI=31,icTCAV=32,icMARK=41,icMONI=42,icAPRT=43)
c$$$      integer icINS,icMAP
c$$$      parameter (icINS=34,icMAP=33)
c$$$      integer icCOORD
c$$$      parameter (icCOORD=35)
c$$$      integer icBEAM,icPROT
c$$$      parameter (icBEAM=36,icPROT=37)
c$$$      integer icSPCH
c$$$      parameter (icSPCH=38)
c$$$      integer icMXEL,icLINE,icNULL
c$$$      parameter(icMXEL=99,icLINE=100,icNULL=0)
c$$$      integer icCELL
c$$$      parameter (icCELL=icLINE)
c$$$      integer icRSVD,icDEF,icACT,icPROC
c$$$      integer icVAR,icKWRD,icUNIT,icRAND,icENV,icFLAG
c$$$      integer icGLI,icGLR,icGLL
c$$$      integer icGRAF,icPART
c$$$      parameter(icRSVD=258,icDEF=icRSVD+2,icACT=icDEF+2,icPROC=icACT+2)
c$$$      parameter(icVAR =icPROC+2, icKWRD=icVAR+2)
c$$$      parameter(icUNIT=icKWRD+2, icRAND=icUNIT+2)
c$$$      parameter(icENV =icRAND+2, icFLAG=icENV+2)
c$$$      parameter(icGLI =icFLAG+4,icGLL=icGLI+4,icGLR =icGLI+8)
c$$$      parameter(icGraf =icGLL+2,icPART=icGraf+2)
