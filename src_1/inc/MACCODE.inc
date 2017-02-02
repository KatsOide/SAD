c     Don't confuse, Emacs. This is -*- fortran -*- mode!
c     internal code list
      integer icDRFT,icBEND,icQUAD,icSEXT,icOCTU,icDECA,icDODECA
      parameter(icDRFT=1,icBEND=2,icQUAD=4,icSEXT=6,
     $     icOCTU=8,icDECA=10,icDODECA=12)
      integer icSOL,icWIG,icUND
      parameter (icSOL=20,icWIG=19,icUND=18)
      integer icST,icMULT,icTEST
      parameter(icST=21,icMULT=22,icTEST=30)
      integer icCAVI,icTCAV,icMARK,icMONI,icAPRT
      parameter(icCAVI=31,icTCAV=32,icMARK=41,icMONI=42,icAPRT=43)
      integer icINS,icMAP
      parameter (icINS=34,icMAP=33)
      integer icCOORD
      parameter (icCOORD=35)
      integer icBEAM,icPROT
      parameter (icBEAM=36,icPROT=37)
      integer icSPCH
      parameter (icSPCH=38)
      integer icMXEL,icLINE,icNULL
      parameter(icMXEL=99,icLINE=100,icNULL=0)
      integer icCELL
      parameter (icCELL=icLINE)
      integer icRSVD,icDEF,icACT,icPROC
      integer icVAR,icKWRD,icUNIT,icRAND,icENV,icFLAG
      integer icGLI,icGLR,icGLL
      integer icGRAF,icPART
      parameter(icRSVD=258,icDEF=icRSVD+2,icACT=icDEF+2,icPROC=icACT+2)
      parameter(icVAR =icPROC+2, icKWRD=icVAR+2)
      parameter(icUNIT=icKWRD+2, icRAND=icUNIT+2)
      parameter(icENV =icRAND+2, icFLAG=icENV+2)
      parameter(icGLI =icFLAG+4,icGLL=icGLI+4,icGLR =icGLI+8)
      parameter(icGraf =icGLL+2,icPART=icGraf+2)
