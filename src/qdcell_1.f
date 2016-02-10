      subroutine qdcell(latt,utwiss,itwissp,
     $     gammab,dtwiss,kk,ll,idp,hstab,vstab,iv,
     $     ctrans,iclast,nfam,nut,disp,nzcod)
      use tfstk
      use ffs
      use tffitcode
      implicit none
      integer*4 nfam,nut,kk,ll,idp,iv,k,l
      integer*4 latt(2,nlat),itwissp(nlat),iclast(-nfam:nfam)
      real*8 gammab(nlat),dtwiss(mfittry),
     $     utwiss(ntwissfun,-nfam:nfam,nut),ctrans(27,-nfam:nfam),
     $     dpsix,dpsiy,cosmux,cosmuy,sinmux,sinmuy,
     $     bxr,byr,trx,try,akx,aky,
     $     x11,x12,x21,x22,y11,y12,y21,y22,
     $     dx11,dx12,dx21,dx22,dy11,dy12,dy21,dy22,
     $     detr,r11,r12,r21,r22,c1,ddetr,
     $     dr11,dr12,dr21,dr22,
     $     s11,s12,s13,s14,s21,s22,s23,s24,
     $     s31,s32,s33,s34,s41,s42,s43,s44,
     $     sr11,sr12,sr13,sr14,sr21,sr22,sr23,sr24,
     $     sr31,sr32,sr33,sr34,sr41,sr42,sr43,sr44,
     $     tm11,tm12,tm13,tm14,tm21,tm22,tm23,tm24,
     $     tm31,tm32,tm33,tm34,tm41,tm42,tm43,tm44,
     $     dalphx,dbetax,dalphy,dbetay,
     $     dsr11,dsr12,dsr13,dsr14,dsr21,dsr22,dsr23,dsr24,
     $     dsr31,dsr32,dsr33,dsr34,dsr41,dsr42,dsr43,dsr44,
     $     dtm11,dtm12,dtm13,dtm14,dtm21,dtm22,dtm23,dtm24,
     $     dtm31,dtm32,dtm33,dtm34,dtm41,dtm42,dtm43,dtm44,
     $     dumm11,dumm12,dumm13,dumm14,
     $     dumm21,dumm22,dumm23,dumm24,
     $     ddumm11,ddumm12,ddumm13,ddumm14,
     $     ddumm21,ddumm22,ddumm23,ddumm24,
     $     ds11,ds12,ds13,ds14,ds21,ds22,ds23,ds24,
     $     ds31,ds32,ds33,ds34,ds41,ds42,ds43,ds44,
     $     da11,da12,da21,da22,da13,da23,
     $     db11,db12,db21,db22,db13,db23,
     $     dmux,dmuy,sinl,cosl,
     $     dr1,dr2,dr3,dr4,
     $     cod1,cod2,cod3,cod4,
     $     dcod1,dcod2,dcod3,dcod4,
     $     ddcod1,ddcod2,ddcod3,ddcod4,
     $     a11,a12,a21,a22,b11,b12,b21,b22,
     $     ep1,ep2,ep3,ep4,dc,dc1,detimx,detimy,xy,
     $     deta,ddeta,aa,daa,detb,ddetb,bb,dbb
      logical*4 cell0,hstab,vstab,disp,nzcod,htrx,htry,normal
c-deb
c      print *,'--------- enter qdcell --k,l=',k,l,'-------'
c          ------------------------------------
      k=itwissp(kk)
      l=itwissp(ll)
      htrx=.false.
      htry=.false.
      cell0=cell
      if(cell)then
        dpsix = utwiss(3,idp,nut) - utwiss(3,idp,1)
        dpsiy = utwiss(6,idp,nut) - utwiss(6,idp,1)
        cosmux = cos(dpsix)
        sinmux = sin(dpsix)
        cosmuy = cos(dpsiy)
        sinmuy = sin(dpsiy)
        bxr=sqrt(utwiss(2,idp,nut)/utwiss(2,idp,1))
        trx=(bxr+1.d0/bxr)*cosmux
     1      +(bxr*utwiss(1,idp,1)-utwiss(1,idp,nut)/bxr)*sinmux
        byr=sqrt(utwiss(5,idp,nut)/utwiss(5,idp,1))
        try=(byr+1.d0/byr)*cosmuy
     1      +(byr*utwiss(4,idp,1)-utwiss(4,idp,nut)/byr)*sinmuy
        htrx=abs(trx) .ge. 2.d0
        htry=abs(try) .ge. 2.d0
        if(htrx)then
          akx=(utwiss(1,idp,nut)-utwiss(1,idp,1))/utwiss(2,idp,1)
        endif
        if(htry)then
          aky=(utwiss(4,idp,nut)-utwiss(4,idp,1))/utwiss(5,idp,1)
        endif
c        cell=.not. (htrx .or. htry)
      else
        sinmux=0.d0
        cosmux=1.d0
        sinmuy=0.d0
        cosmuy=1.d0
        dpsix=0.d0
        dpsiy=0.d0
      endif
      call qdtwis(latt,utwiss,itwissp,
     $     gammab,dtwiss,ctrans,iclast,
     $     kk,ll,idp,iv,nfam,nut,disp,nzcod)
c          ------------------------------------
      if(.not. cell)then
        cell=cell0
        return
      endif
      cell=cell0
C====== Get Transfer Matrix tm**  =========================
      x11=cosmux+utwiss(1,idp,l)*sinmux
      x22=cosmux-utwiss(1,idp,l)*sinmux
      y11=cosmuy+utwiss(4,idp,l)*sinmuy
      y22=cosmuy-utwiss(4,idp,l)*sinmuy
      x12=utwiss(2,idp,l)*sinmux
      x21=-(1.d0+utwiss(1,idp,l)*utwiss(1,idp,l))*sinmux
     1                                         /utwiss(2,idp,l)
      y12=utwiss(5,idp,l)*sinmuy
      y21=-(1.d0+utwiss(4,idp,l)*utwiss(4,idp,l))*sinmuy
     1                                         /utwiss(5,idp,l)
C
C   ( transformation matrix at point l )
      r11 = utwiss(mfitr1,idp,l)
      r12 = utwiss(mfitr2,idp,l)
      r21 = utwiss(mfitr3,idp,l)
      r22 = utwiss(mfitr4,idp,l)
      detr=r11*r22-r12*r21
      c1 = sqrt(1.d0 - detr)
        normal=utwiss(mfitdetr,idp,l) .lt. 1.d0
        if(normal) then
         s11 = c1
         s12 = 0.d0
         s13 = -r22
         s14 =  r12
         s21 = 0.d0
         s22 = c1
         s23 =  r21
         s24 = -r11
         s31 =  r11
         s32 =  r12
         s33 = c1
         s34 = 0.d0
         s41 =  r21
         s42 =  r22
         s43 = 0.d0
         s44 = c1
         sr11 = c1
         sr12 = 0.d0
         sr13 =  r22
         sr14 = -r12
         sr21 = 0.d0
         sr22 = c1
         sr23 = -r21
         sr24 =  r11
         sr31 = -r11
         sr32 = -r12
         sr33 = c1
         sr34 = 0.d0
         sr41 = -r21
         sr42 = -r22
         sr43 = 0.d0
         sr44 = c1
            tm11 =  sr11*(x11*s11)
     1             +sr13*(y11*s31 + y12*s41)
     1             +sr14*(y21*s31 + y22*s41)
            tm12 =  sr11*(x12*s22)
     1             +sr13*(y11*s32 + y12*s42)
     1             +sr14*(y21*s32 + y22*s42)
            tm13 =  sr11*(x11*s13 + x12*s23)
     1             +sr13*(y11*s33)
     1             +sr14*(y21*s33)
            tm14 =  sr11*(x11*s14 + x12*s24)
     1             +sr13*(y12*s44)
     1             +sr14*(y22*s44)
            tm21 =  
     1             +sr22*(x21*s11)
     1             +sr23*(y11*s31 + y12*s41)
     1             +sr24*(y21*s31 + y22*s41)
            tm22 =  
     1             +sr22*(x22*s22)
     1             +sr23*(y11*s32 + y12*s42)
     1             +sr24*(y21*s32 + y22*s42)
            tm23 =  
     1             +sr22*(x21*s13 + x22*s23)
     1             +sr23*(y11*s33)
     1             +sr24*(y21*s33)
            tm24 =  
     1             +sr22*(x21*s14 + x22*s24)
     1             +sr23*(y12*s44)
     1             +sr24*(y22*s44)
            tm31 =  sr31*(x11*s11)
     1             +sr32*(x21*s11)
     1             +sr33*(y11*s31 + y12*s41)
            tm32 =  sr31*(x12*s22)
     1             +sr32*(x22*s22)
     1             +sr33*(y11*s32 + y12*s42)
            tm33 =  sr31*(x11*s13 + x12*s23)
     1             +sr32*(x21*s13 + x22*s23)
     1             +sr33*(y11*s33)
            tm34 =  sr31*(x11*s14 + x12*s24)
     1             +sr32*(x21*s14 + x22*s24)
     1             +sr33*(y12*s44)
            tm41 =  sr41*(x11*s11)
     1             +sr42*(x21*s11)
     1             +sr44*(y21*s31 + y22*s41)
            tm42 =  sr41*(x12*s22)
     1             +sr42*(x22*s22)
     1             +sr44*(y21*s32 + y22*s42)
            tm43 =  sr41*(x11*s13 + x12*s23)
     1             +sr42*(x21*s13 + x22*s23)
     1             +sr44*(y21*s33)
            tm44 =  sr41*(x11*s14 + x12*s24)
     1             +sr42*(x21*s14 + x22*s24)
     1             +sr44*(y22*s44)
        else
c         print *,'#### QDCELL  !!!case of detr>1'
         s11 = -r22
         s12 =  r12
         s13 = c1
         s14 = 0.d0
         s21 =  r21
         s22 = -r11
         s23 = 0.d0
         s24 = c1
         s31 = c1
         s32 = 0.d0
         s33 =  r11
         s34 =  r12
         s41 = 0.d0
         s42 = c1
         s43 =  r21
         s44 =  r22
         sr11 =  r22
         sr12 = -r12
         sr13 = c1
         sr14 = 0.d0
         sr21 = -r21
         sr22 =  r11
         sr23 = 0.d0
         sr24 = c1
         sr31 = c1
         sr32 = 0.d0
         sr33 = -r11
         sr34 = -r12
         sr41 = 0.d0
         sr42 = c1
         sr43 = -r21
         sr44 = -r22
            tm11 =  sr11*(x11*s11 + x12*s21)
     1             +sr12*(x21*s11 + x22*s21)
     1             +sr13*(y11*s31)
            tm12 =  sr11*(x11*s12 + x12*s22)
     1             +sr12*(x21*s12 + x22*s22)
     1             +sr13*(y12*s42)
            tm13 =  sr11*(x11*s13)
     1             +sr12*(x21*s13)
     1             +sr13*(y11*s33 + y12*s43)
            tm14 =  sr11*(x12*s24)
     1             +sr12*(x22*s24)
     1             +sr13*(y11*s34 + y12*s44)
            tm21 =  sr21*(x11*s11 + x12*s21)
     1             +sr22*(x21*s11 + x22*s21)
     1             +sr24*(y21*s31)
            tm22 =  sr21*(x11*s12 + x12*s22)
     1             +sr22*(x21*s12 + x22*s22)
     1             +sr24*(y22*s42)
            tm23 =  sr21*(x11*s13)
     1             +sr22*(x21*s13)
     1             +sr24*(y21*s33 + y22*s43)
            tm24 =  sr21*(x12*s24)
     1             +sr22*(x22*s24)
     1             +sr24*(y21*s34 + y22*s44)
            tm31 =  sr31*(x11*s11 + x12*s21)
     1             +sr33*(y11*s31)
     1             +sr34*(y21*s31)
            tm32 =  sr31*(x11*s12 + x12*s22)
     1             +sr33*(y12*s42)
     1             +sr34*(y22*s42)
            tm33 =  sr31*(x11*s13)
     1             +sr33*(y11*s33 + y12*s43)
     1             +sr34*(y21*s33 + y22*s43)
            tm34 =  sr31*(x12*s24)
     1             +sr33*(y11*s34 + y12*s44)
     1             +sr34*(y21*s34 + y22*s44)
            tm41 =  
     1             +sr42*(x21*s11 + x22*s21)
     1             +sr43*(y11*s31)
     1             +sr44*(y21*s31)
            tm42 =  
     1             +sr42*(x21*s12 + x22*s22)
     1             +sr43*(y12*s42)
     1             +sr44*(y22*s42)
            tm43 =  
     1             +sr42*(x21*s13)
     1             +sr43*(y11*s33 + y12*s43)
     1             +sr44*(y21*s33 + y22*s43)
            tm44 =  
     1             +sr42*(x22*s24)
     1             +sr43*(y11*s34 + y12*s44)
     1             +sr44*(y21*s34 + y22*s44)
        end if
C   ( 4*4 transfer matrix( assuming correct COD found ))
C
C         M  = s-inv. Md .s
C
C====== Get derivative of Transfer Matrix dtm**  ==========
C
C        dM  = dsf-inv. Md .s + sf-inv. dMd .s
C
      dalphx = (1.d0+utwiss(1,idp,l)*utwiss(1,idp,l))*dtwiss(1)
      dbetax = 0.5d0*dtwiss(2)
      dalphy = (1.d0+utwiss(4,idp,l)*utwiss(4,idp,l))*dtwiss(4)
      dbetay = 0.5d0*dtwiss(5)
      dx11=dbetax*x11
     1    +(-sinmux+utwiss(1,idp,l)*cosmux)*dtwiss(3)
      dx12=dbetax*x12
     1    +utwiss(2,idp,l)*cosmux*dtwiss(3)
      dx21=-x21*dbetax+(
     1      -(utwiss(1,idp,l)*sinmux+cosmux)*dalphx
     1      -(1.d0+utwiss(1,idp,l)**2)*cosmux*dtwiss(3))
     1    /utwiss(2,idp,l)
      dx22=-dbetax*x22
     1     -sinmux*dalphx
     1     -(sinmux+utwiss(1,idp,l)*cosmux)*dtwiss(3)
c---
      dy11=dbetay*y11
     1    +(-sinmuy+utwiss(4,idp,l)*cosmuy)*dtwiss(6)
      dy12=dbetay*y12
     1    +utwiss(5,idp,l)*cosmuy*dtwiss(6)
      dy21=-y21*dbetay+(
     1      -(utwiss(4,idp,l)*sinmuy+cosmuy)*dalphy
     1      -(1.d0+utwiss(4,idp,l)**2)*cosmuy*dtwiss(6))
     1    /utwiss(5,idp,l)
      dy22=-dbetay*y22
     1     -sinmuy*dalphy
     1     -(sinmuy+utwiss(4,idp,l)*cosmuy)*dtwiss(6)
c---
c     c^2+Det(R)=1   ---> dc = - d(Det(R))/(2c)
      dr11 = dtwiss(mfitr1)
      dr12 = dtwiss(mfitr2)
      dr21 = dtwiss(mfitr3)
      dr22 = dtwiss(mfitr4)
      ddetr = dr11*r22 + dr22*r11
     1       -dr12*r21 - dr21*r12
      dc1   = -0.5d0*ddetr/c1
        if(normal) then
         dsr11 = dc1
         dsr12 = 0.d0
         dsr13 =  dr22
         dsr14 = -dr12
         dsr21 = 0.d0
         dsr22 = dc1
         dsr23 = -dr21
         dsr24 =  dr11
         dsr31 = -dr11
         dsr32 = -dr12
         dsr33 = dc1
         dsr34 = 0.d0
         dsr41 = -dr21
         dsr42 = -dr22
         dsr43 = 0.d0
         dsr44 = dc1
      dtm11 =  dsr11*(x11*s11)+sr11*(dx11*s11)
     1     +dsr13*(y11*s31 + y12*s41)+sr13*(dy11*s31 + dy12*s41)
     1     +dsr14*(y21*s31 + y22*s41)+sr14*(dy21*s31 + dy22*s41)
      dtm12 =  dsr11*(x12*s22)+sr11*(dx12*s22)
     1     +dsr13*(y11*s32 + y12*s42)+sr13*(dy11*s32 + dy12*s42)
     1     +dsr14*(y21*s32 + y22*s42)+sr14*(dy21*s32 + dy22*s42)
      dtm13 =  dsr11*(x11*s13 + x12*s23)+sr11*(dx11*s13 + dx12*s23)
     1     +dsr13*(y11*s33)+sr13*(dy11*s33)
     1     +dsr14*(y21*s33)+sr14*(dy21*s33)
      dtm14 =  dsr11*(x11*s14 + x12*s24)+sr11*(dx11*s14 + dx12*s24)
     1     +dsr13*(y12*s44)+sr13*(dy12*s44)
     1     +dsr14*(y22*s44)+sr14*(dy22*s44)
      dtm21 =  
     1     +dsr22*(x21*s11)+sr22*(dx21*s11)
     1     +dsr23*(y11*s31 + y12*s41)+sr23*(dy11*s31 + dy12*s41)
     1     +dsr24*(y21*s31 + y22*s41)+sr24*(dy21*s31 + dy22*s41)
      dtm22 =  
     1     +dsr22*(x22*s22)+sr22*(dx22*s22)
     1     +dsr23*(y11*s32 + y12*s42)+sr23*(dy11*s32 + dy12*s42)
     1     +dsr24*(y21*s32 + y22*s42)+sr24*(dy21*s32 + dy22*s42)
      dtm23 =  
     1     +dsr22*(x21*s13 + x22*s23)+sr22*(dx21*s13 + dx22*s23)
     1     +dsr23*(y11*s33)+sr23*(dy11*s33 )
     1     +dsr24*(y21*s33)+sr24*(dy21*s33 )
      dtm24 =  
     1     +dsr22*(x21*s14 + x22*s24)+sr22*(dx21*s14 + dx22*s24)
     1     +dsr23*(y12*s44)+sr23*(dy12*s44)
     1     +dsr24*(y22*s44)+sr24*(dy22*s44)
      dtm31 =  dsr31*(x11*s11)+sr31*(dx11*s11)
     1     +dsr32*(x21*s11)+sr32*(dx21*s11)
     1     +dsr33*(y11*s31 + y12*s41)+sr33*(dy11*s31 + dy12*s41)
      dtm32 =  dsr31*(x12*s22)+sr31*(dx12*s22)
     1     +dsr32*(x22*s22)+sr32*(dx22*s22)
     1     +dsr33*(y11*s32 + y12*s42)+sr33*(dy11*s32 + dy12*s42)
      dtm33 =  dsr31*(x11*s13 + x12*s23)+sr31*(dx11*s13 + dx12*s23)
     1     +dsr32*(x21*s13 + x22*s23)+sr32*(dx21*s13 + dx22*s23)
     1     +dsr33*(y11*s33 )+sr33*(dy11*s33 )
      dtm34 =  dsr31*(x11*s14 + x12*s24)+sr31*(dx11*s14 + dx12*s24)
     1     +dsr32*(x21*s14 + x22*s24)+sr32*(dx21*s14 + dx22*s24)
     1     +dsr33*(y12*s44)+sr33*(dy12*s44)
      dtm41 =  dsr41*(x11*s11)+sr41*(dx11*s11)
     1     +dsr42*(x21*s11)+sr42*(dx21*s11)
     1     +dsr44*(y21*s31 + y22*s41)+sr44*(dy21*s31 + dy22*s41)
      dtm42 =  dsr41*(x12*s22)+sr41*(dx12*s22)
     1     +dsr42*(x22*s22)+sr42*(dx22*s22)
     1     +dsr44*(y21*s32 + y22*s42)+sr44*(dy21*s32 + dy22*s42)
      dtm43 =  dsr41*(x11*s13 + x12*s23)+sr41*(dx11*s13 + dx12*s23)
     1     +dsr42*(x21*s13 + x22*s23)+sr42*(dx21*s13 + dx22*s23)
     1     +dsr44*(y21*s33 )+sr44*(dy21*s33 )
      dtm44 =  dsr41*(x11*s14 + x12*s24)+sr41*(dx11*s14 + dx12*s24)
     1     +dsr42*(x21*s14 + x22*s24)+sr42*(dx21*s14 + dx22*s24)
     1     +dsr44*(y22*s44)+sr44*(dy22*s44)
        else
c         print *,'#### QDCELL  !!!case of detr>1'
         dsr11 =  dr22
         dsr12 = -dr12
         dsr13 = dc1
         dsr14 = 0.d0
         dsr21 = -dr21
         dsr22 =  dr11
         dsr23 = 0.d0
         dsr24 = dc1
         dsr31 = dc1
         dsr32 = 0.d0
         dsr33 = -dr11
         dsr34 = -dr12
         dsr41 = 0.d0
         dsr42 = dc1
         dsr43 = -dr21
         dsr44 = -dr22
      dtm11 =  dsr11*(x11*s11 + x12*s21)+sr11*(dx11*s11 + dx12*s21)
     1     +dsr12*(x21*s11 + x22*s21)+sr12*(dx21*s11 + dx22*s21)
     1     +dsr13*(y11*s31)+sr13*(dy11*s31)
      dtm12 =  dsr11*(x11*s12 + x12*s22)+sr11*(dx11*s12 + dx12*s22)
     1     +dsr12*(x21*s12 + x22*s22)+sr12*(dx21*s12 + dx22*s22)
     1     +dsr13*(y12*s42)+sr13*(dy12*s42)
      dtm13 =  dsr11*(x11*s13)+sr11*(dx11*s13)
     1     +dsr12*(x21*s13)+sr12*(dx21*s13)
     1     +dsr13*(y11*s33 + y12*s43)+sr13*(dy11*s33 + dy12*s43)
      dtm14 =  dsr11*(x12*s24)+sr11*(dx12*s24)
     1     +dsr12*(x22*s24)+sr12*(dx22*s24)
     1     +dsr13*(y11*s34 + y12*s44)+sr13*(dy11*s34 + dy12*s44)
      dtm21 =  dsr21*(x11*s11 + x12*s21)+sr21*(dx11*s11 + dx12*s21)
     1     +dsr22*(x21*s11 + x22*s21)+sr22*(dx21*s11 + dx22*s21)
     1     +dsr24*(y21*s31)+sr24*(dy21*s31)
      dtm22 =  dsr21*(x11*s12 + x12*s22)+sr21*(dx11*s12 + dx12*s22)
     1     +dsr22*(x21*s12 + x22*s22)+sr22*(dx21*s12 + dx22*s22)
     1     +dsr24*(y22*s42)+sr24*(dy22*s42)
      dtm23 =  dsr21*(x11*s13)+sr21*(dx11*s13)
     1     +dsr22*(x21*s13)+sr22*(dx21*s13)
     1     +dsr24*(y21*s33 + y22*s43)+sr24*(dy21*s33 + dy22*s43)
      dtm24 =  dsr21*(x12*s24)+sr21*(dx12*s24)
     1     +dsr22*(x21*s14 + x22*s24)+sr22*(dx21*s14 + dx22*s24)
     1     +dsr24*(y21*s34 + y22*s44)+sr24*(dy21*s34 + dy22*s44)
      dtm31 =  dsr31*(x11*s11 + x12*s21)+sr31*(dx11*s11 + dx12*s21)
     1     +dsr33*(y11*s31)+sr33*(dy11*s31)
     1     +dsr34*(y21*s31)+sr34*(dy21*s31)
      dtm32 =  dsr31*(x11*s12 + x12*s22)+sr31*(dx11*s12 + dx12*s22)
     1     +dsr33*(y12*s42)+sr33*(dy12*s42)
     1     +dsr34*(y22*s42)+sr34*(dy22*s42)
      dtm33 =  dsr31*(x11*s13)+sr31*(dx11*s13)
     1     +dsr33*(y11*s33 + y12*s43)+sr33*(dy11*s33 + dy12*s43)
     1     +dsr34*(y21*s33 + y22*s43)+sr34*(dy21*s33 + dy22*s43)
      dtm34 =  dsr31*(x12*s24)+sr31*(dx12*s24)
     1     +dsr33*(y11*s34 + y12*s44)+sr33*(dy11*s34 + dy12*s44)
     1     +dsr34*(y21*s34 + y22*s44)+sr34*(dy21*s34 + dy22*s44)
      dtm41 =  
     1     +dsr42*(x21*s11 + x22*s21)+sr42*(dx21*s11 + dx22*s21)
     1     +dsr43*(y11*s31)+sr43*(dy11*s31)
     1     +dsr44*(y21*s31)+sr44*(dy21*s31)
      dtm42 =  
     1     +dsr42*(x21*s12 + x22*s22)+sr42*(dx21*s12 + dx22*s22)
     1     +dsr43*(y12*s42)+sr43*(dy12*s42)
     1     +dsr44*(y22*s42)+sr44*(dy22*s42)
      dtm43 =  
     1     +dsr42*(x21*s13)+sr42*(dx21*s13)
     1     +dsr43*(y11*s33 + y12*s43)+sr43*(dy11*s33 + dy12*s43)
     1     +dsr44*(y21*s33 + y22*s43)+sr44*(dy21*s33 + dy22*s43)
      dtm44 =  
     1     +dsr42*(x22*s24)+sr42*(dx22*s24)
     1     +dsr43*(y11*s34 + y12*s44)+sr43*(dy11*s34 + dy12*s44)
     1     +dsr44*(y21*s34 + y22*s44)+sr44*(dy21*s34 + dy22*s44)
        end if
C====== Get derivative of transformation matrix dr*  ======
c          with periodic boundary condition
      if(normal)then
        xy=1.d0
      else
        xy=-1.d0
         dumm11=tm11
         dumm12=tm12
         dumm21=tm21
         dumm22=tm22
         tm11  =tm33
         tm12  =tm34
         tm21  =tm43
         tm22  =tm44
         tm33  =dumm11
         tm34  =dumm12
         tm43  =dumm21
         tm44  =dumm22
         dumm13=tm13
         dumm14=tm14
         dumm23=tm23
         dumm24=tm24
         tm13  =tm31
         tm14  =tm32
         tm23  =tm41
         tm24  =tm42
         tm31  =dumm13
         tm32  =dumm14
         tm41  =dumm23
         tm42  =dumm24
c        --------------------
         ddumm11=dtm11
         ddumm12=dtm12
         ddumm21=dtm21
         ddumm22=dtm22
         dtm11  =dtm33
         dtm12  =dtm34
         dtm21  =dtm43
         dtm22  =dtm44
         dtm33  =ddumm11
         dtm34  =ddumm12
         dtm43  =ddumm21
         dtm44  =ddumm22
         ddumm13=dtm13
         ddumm14=dtm14
         ddumm23=dtm23
         ddumm24=dtm24
         dtm13  =dtm31
         dtm14  =dtm32
         dtm23  =dtm41
         dtm24  =dtm42
         dtm31  =ddumm13
         dtm32  =ddumm14
         dtm41  =ddumm23
         dtm42  =ddumm24
         s11 = c1
         s12 = 0.d0
         s13 = -r22
         s14 =  r12
         s21 = 0.d0
         s22 = c1
         s23 =  r21
         s24 = -r11
         s31 =  r11
         s32 =  r12
         s33 = c1
         s34 = 0.d0
         s41 =  r21
         s42 =  r22
         s43 = 0.d0
         s44 = c1
      end if
      write(*,'(4(1p4g12.4/))')
c     $     tm11,tm12,tm13,tm14,
c     1     tm21,tm22,tm23,tm24,
c     1     tm31,tm32,tm33,tm34,
c     1     tm41,tm42,tm43,tm44,
     1     dtm11,dtm12,dtm13,dtm14,
     1     dtm21,dtm22,dtm23,dtm24,
     1     dtm31,dtm32,dtm33,dtm34,
     1     dtm41,dtm42,dtm43,dtm44
c          -------------------------------
      call qdmdia(tm11,tm12,tm13,tm14,
     1            tm21,tm22,tm23,tm24,
     1            tm31,tm32,tm33,tm34,
     1            tm41,tm42,tm43,tm44,
     1            dtm11,dtm12,dtm13,dtm14,
     1            dtm21,dtm22,dtm23,dtm24,
     1            dtm31,dtm32,dtm33,dtm34,
     1            dtm41,dtm42,dtm43,dtm44,
     1            dr1,dr2,dr3,dr4,r11)
c          -------------------------------
      dtwiss(mfitr1) = dr1
      dtwiss(mfitr2) = dr2
      dtwiss(mfitr3) = dr3
      dtwiss(mfitr4) = dr4
c      write(*,'(a,1p6g15.7)')'qdcell ',dr1,dr2,dr3,dr4,detr,ddetr
C====== Get derivative of transformation matrix ds**  =====
c          with periodic boundary condition
      dc   = -0.5d0*(dr1*r22+r11*dr4-dr2*r21-r12*dr3)/c1
      ds11 = dc
      ds12 = 0.d0
      ds13 = -dr4
      ds14 = dr2
      ds21 = 0.d0
      ds22 = dc
      ds23 = dr3
      ds24 = -dr1
      ds31 = dr1
      ds32 = dr2
      ds33 = dc
      ds34 = 0.d0
      ds41 = dr3
      ds42 = dr4
      ds43 = 0.d0
      ds44 = dc
C------dA,dB:derivative of diagonalized transfer matrix ------
         da11 = ds11*(tm11*s11-tm13*s31-tm14*s41)
     1        + ds13*(tm31*s11-tm33*s31-tm34*s41)
     1        + ds14*(tm41*s11-tm43*s31-tm44*s41)
     1        + s11*(dtm11*s11-dtm13*s31-dtm14*s41)
     1        + s13*(dtm31*s11-dtm33*s31-dtm34*s41)
     1        + s14*(dtm41*s11-dtm43*s31-dtm44*s41)
     1        + s11*(tm11*ds11-tm13*ds31-tm14*ds41)
     1        + s13*(tm31*ds11-tm33*ds31-tm34*ds41)
     1        + s14*(tm41*ds11-tm43*ds31-tm44*ds41)
         da12 = ds11*(tm12*s22-tm13*s32-tm14*s42)
     1        + ds13*(tm32*s22-tm33*s32-tm34*s42)
     1        + ds14*(tm42*s22-tm43*s32-tm44*s42)
     1        + s11*(dtm12*s22-dtm13*s32-dtm14*s42)
     1        + s13*(dtm32*s22-dtm33*s32-dtm34*s42)
     1        + s14*(dtm42*s22-dtm43*s32-dtm44*s42)
     1        + s11*(tm12*ds22-tm13*ds32-tm14*ds42)
     1        + s13*(tm32*ds22-tm33*ds32-tm34*ds42)
     1        + s14*(tm42*ds22-tm43*ds32-tm44*ds42)
         db21 = 
     1        + ds22*(tm21*s11+tm22*s21-tm23*s31-tm24*s41)
     1        + ds23*(tm31*s11+tm32*s21-tm33*s31-tm34*s41)
     1        + ds24*(tm41*s11+tm42*s21-tm43*s31-tm44*s41)
     1        + s22*(dtm21*s11+dtm22*s21-dtm23*s31-dtm24*s41)
     1        + s23*(dtm31*s11+dtm32*s21-dtm33*s31-dtm34*s41)
     1        + s24*(dtm41*s11+dtm42*s21-dtm43*s31-dtm44*s41)
     1        + s22*(tm21*ds11+tm22*ds21-tm23*ds31-tm24*ds41)
     1        + s23*(tm31*ds11+tm32*ds21-tm33*ds31-tm34*ds41)
     1        + s24*(tm41*ds11+tm42*ds21-tm43*ds31-tm44*ds41)
         da22 = 
     1        + ds22*(tm22*s22-tm23*s32-tm24*s42)
     1        + ds23*(tm32*s22-tm33*s32-tm34*s42)
     1        + ds24*(tm42*s22-tm43*s32-tm44*s42)
     1        + s22*(dtm22*s22-dtm23*s32-dtm24*s42)
     1        + s23*(dtm32*s22-dtm33*s32-dtm34*s42)
     1        + s24*(dtm42*s22-dtm43*s32-dtm44*s42)
     1        + s22*(tm22*ds22-tm23*ds32-tm24*ds42)
     1        + s23*(tm32*ds22-tm33*ds32-tm34*ds42)
     1        + s24*(tm42*ds22-tm43*ds32-tm44*ds42)
         db11 = ds31*(-tm11*s13-tm12*s23+tm13*s33)
     1        + ds32*(-tm21*s13-tm22*s23+tm23*s33)
     1        + ds33*(-tm31*s13-tm32*s23+tm33*s33)
     1        + s31*(-dtm11*s13-dtm12*s23+dtm13*s33)
     1        + s32*(-dtm21*s13-dtm22*s23+dtm23*s33)
     1        + s33*(-dtm31*s13-dtm32*s23+dtm33*s33)
     1        + s31*(-tm11*ds13-tm12*ds23+tm13*ds33)
     1        + s32*(-tm21*ds13-tm22*ds23+tm23*ds33)
     1        + s33*(-tm31*ds13-tm32*ds23+tm33*ds33)
         db12 = ds31*(-tm11*s14-tm12*s24+tm14*s44)
     1        + ds32*(-tm21*s14-tm22*s24+tm24*s44)
     1        + ds33*(-tm31*s14-tm32*s24+tm34*s44)
     1        + s31*(-dtm11*s14-dtm12*s24+dtm14*s44)
     1        + s32*(-dtm21*s14-dtm22*s24+dtm24*s44)
     1        + s33*(-dtm31*s14-dtm32*s24+dtm34*s44)
     1        + s31*(-tm11*ds14-tm12*ds24+tm14*ds44)
     1        + s32*(-tm21*ds14-tm22*ds24+tm24*ds44)
     1        + s33*(-tm31*ds14-tm32*ds24+tm34*ds44)
         db21 = ds41*(-tm11*s13-tm12*s23+tm13*s33)
     1        + ds42*(-tm21*s13-tm22*s23+tm23*s33)
     1        + ds44*(-tm41*s13-tm42*s23+tm43*s33)
     1        + s41*(-dtm11*s13-dtm12*s23+dtm13*s33)
     1        + s42*(-dtm21*s13-dtm22*s23+dtm23*s33)
     1        + s44*(-dtm41*s13-dtm42*s23+dtm43*s33)
     1        + s41*(-tm11*ds13-tm12*ds23+tm13*ds33)
     1        + s42*(-tm21*ds13-tm22*ds23+tm23*ds33)
     1        + s44*(-tm41*ds13-tm42*ds23+tm43*ds33)
         db22 = ds41*(-tm11*s14-tm12*s24+tm14*s44)
     1        + ds42*(-tm21*s14-tm22*s24+tm24*s44)
     1        + ds44*(-tm41*s14-tm42*s24+tm44*s44)
     1        + s41*(-dtm11*s14-dtm12*s24+dtm14*s44)
     1        + s42*(-dtm21*s14-dtm22*s24+dtm24*s44)
     1        + s44*(-dtm41*s14-dtm42*s24+dtm44*s44)
     1        + s41*(-tm11*ds14-tm12*ds24+tm14*ds44)
     1        + s42*(-tm21*ds14-tm22*ds24+tm24*ds44)
     1        + s44*(-tm41*ds14-tm42*ds24+tm44*ds44)
c        deta=tm11*tm22-tm12*tm21
c        ddeta=dtm11*tm22-dtm12*tm21
c     $       +tm11*dtm22-tm12*dtm21
c        aa=1.d0/sqrt(deta)
c        daa=-.5d0*ddeta/deta*aa
c        da11=dtm11*aa+tm11*daa
c        da12=dtm12*aa+tm12*daa
c        da21=dtm21*aa+tm21*daa
c        da22=dtm22*aa+tm22*daa
c        detb=tm33*tm44-tm34*tm43
c        ddetb=dtm33*tm44-dtm34*tm43
c     $       +tm33*dtm44-tm34*dtm43
c        bb=1.d0/sqrt(detb)
c        dbb=-.5d0*ddetb/detb*bb
c        db11=dtm33*bb+tm33*dbb
c        db12=dtm34*bb+tm34*dbb
c        db21=dtm43*bb+tm43*dbb
c        db22=dtm44*bb+tm44*dbb
c        write(*,'(a,2(1p8g15.7/))')'qdcell-a ',
c     $       x11,tm11*aa,x12,tm12*aa,
c     $       x21,tm21*aa,x22,tm22*aa,
c     $       y11,tm33*bb,y12,tm34*bb,
c     $       y21,tm43*bb,y22,tm44*bb
c-deb
c     print *,'da11,da12=',da11,da12
c     print *,'da21,da22=',da21,da22
C====== Get derivative of twiss parameters  ===============
C
C   Changed by K. O.  6/5/1992
C   Dtwiss(3) and dtwiss(6) below are still not perfect.
C   The change of R at l is not included yet.
C
C
       if(sinmux.ne.0.d0) then
         dmux=-0.5d0*(da11+da22)/sinmux
         dtwiss(1)=(.5d0*(da11-da22)
     $        -utwiss(1,idp,l)*cosmux*dmux)/sinmux
         dtwiss(2)=(da12-utwiss(2,idp,l)*cosmux*dmux)/x12
       else
         dmux=0.d0
       endif
       if(k .lt. l)then
         sinl=sin(dpsix-utwiss(3,idp,l))
         cosl=cos(dpsix-utwiss(3,idp,l))
         dtwiss(3)=dmux+sinl*((cosl+utwiss(1,idp,l)*sinl)*dtwiss(2)
     1        -sinl*dtwiss(1))
       else
         sinl=sin(utwiss(3,idp,l))
         cosl=cos(utwiss(3,idp,l))
         dtwiss(3)=-sinl*((cosl-utwiss(1,idp,l)*sinl)*dtwiss(2)
     1        +sinl*dtwiss(1))
       endif
       write(*,'(a,1p3g15.7)')'qdcell-axbx ',
     $      dtwiss(1),dtwiss(2)*utwiss(2,idp,l),dtwiss(3)
       dtwiss(1)=dtwiss(1)/(1.d0+utwiss(1,idp,l)*utwiss(1,idp,l))
       if(sinmuy.ne.0.d0) then
         dmuy=-0.5d0*(db11+db22)/sinmuy
         dtwiss(4)=(.5d0*(db11-db22)
     $        -utwiss(4,idp,l)*cosmuy*dmuy)/sinmuy
         dtwiss(5)=(db12-utwiss(5,idp,l)*cosmuy*dmuy)/y12
       else
         dmuy=0.d0
       end if
       if(k .lt. l)then
C     
C     A bug was found below on 7/1/1994. sin(dpsix- ....)
C     
         sinl=sin(dpsiy-utwiss(6,idp,l))
         cosl=cos(dpsiy-utwiss(6,idp,l))
         dtwiss(6)=dmuy+sinl
     1        *((cosl+utwiss(4,idp,l)*sinl)*dtwiss(5)
     1        -sinl*dtwiss(4))
       else
         sinl=sin(utwiss(6,idp,l))
         cosl=cos(utwiss(6,idp,l))
         dtwiss(6)=-sinl
     1        *((cosl-utwiss(4,idp,l)*sinl)*dtwiss(5)
     1        +sinl*dtwiss(4))
       endif
       write(*,'(a,1p3g15.7)')'qdcell-ayby ',
     $      dtwiss(4),dtwiss(5)*utwiss(5,idp,l),dtwiss(6)
       dtwiss(4)=dtwiss(4)/(1.d0+utwiss(4,idp,l)*utwiss(4,idp,l))
c      write(*,'(a,3i5,1p6g12.4)')'qdcell ',k,l,idp,
c     $     dtwiss(1),dtwiss(2),dtwiss(3),
c     $     dtwiss(4),dtwiss(5),dtwiss(6)
c      if(dtwiss(3) .eq. 0.d0 .or.
c     $     abs(dtwiss(3)) .gt. 1.d10)then
c      endif
c      write(*,*)'qdcell-y ',k,l,idp,dtwiss(4),dtwiss(5),dtwiss(6)
c------- dtwiss(1) = dalpha/(1+alpha^2)
c        dtwiss(2) = dbeta/beta
C====== Get derivative of orbit  ==========================
c..........going to 2*2 world
      cod1 = s11*utwiss(mfitdx,idp,l)+s12*utwiss(mfitdpx,idp,l)
     1     + s13*utwiss(mfitdy,idp,l)+s14*utwiss(mfitdpy,idp,l)
      cod2 = s21*utwiss(mfitdx,idp,l)+s22*utwiss(mfitdpx,idp,l)
     1     + s23*utwiss(mfitdy,idp,l)+s24*utwiss(mfitdpy,idp,l)
      cod3 = s31*utwiss(mfitdx,idp,l)+s32*utwiss(mfitdpx,idp,l)
     1     + s33*utwiss(mfitdy,idp,l)+s34*utwiss(mfitdpy,idp,l)
      cod4 = s41*utwiss(mfitdx,idp,l)+s42*utwiss(mfitdpx,idp,l)
     1     + s43*utwiss(mfitdy,idp,l)+s44*utwiss(mfitdpy,idp,l)
      ddcod1 = s11*dtwiss(mfitdx)+s12*dtwiss(mfitdpx)
     1     + s13*dtwiss(mfitdy)+s14*dtwiss(mfitdpy)
      ddcod2 = s21*dtwiss(mfitdx)+s22*dtwiss(mfitdpx)
     1     + s23*dtwiss(mfitdy)+s24*dtwiss(mfitdpy)
      ddcod3 = s31*dtwiss(mfitdx)+s32*dtwiss(mfitdpx)
     1     + s33*dtwiss(mfitdy)+s34*dtwiss(mfitdpy)
      ddcod4 = s41*dtwiss(mfitdx)+s42*dtwiss(mfitdpx)
     1     + s43*dtwiss(mfitdy)+s44*dtwiss(mfitdpy)
c..........dx = (I-M)^-1.(dM.x+delx)
      a11 = cosmux+utwiss(1,idp,l)*sinmux
      a12 = utwiss(2,idp,l)*sinmux
      a21 = -(1.d0+utwiss(1,idp,l)*utwiss(1,idp,l))*sinmux
     1                           /utwiss(2,idp,l)
      a22 = cosmux-utwiss(1,idp,l)*sinmux
      b11 = cosmuy+utwiss(4,idp,l)*sinmuy
      b12 = utwiss(5,idp,l)*sinmuy
      b21 = -(1.d0+utwiss(4,idp,l)*utwiss(4,idp,l))*sinmuy
     1                           /utwiss(5,idp,l)
      b22 = cosmuy-utwiss(4,idp,l)*sinmuy
      detimx = (1.d0-a11)*(1.d0-a22)-a12*a21
      detimy = (1.d0-b11)*(1.d0-b22)-b12*b21
      if(detimx.ne.0.d0) then
        dcod1 = ((1.d0-a22)*(da11*cod1+da12*cod2+ddcod1)
     1          +a12*(da21*cod1+da22*cod2+ddcod2))/detimx
        dcod2 = (a21*(da11*cod1+da12*cod2+ddcod1)
     1          +(1.d0-a11)*(da21*cod1+da22*cod2+ddcod2))/detimx
      else
c        write(*,*)'qdcell @ src/qdcell.f: ',
c     $        '2x2 horizontal matrix (I - A) is degenerated!',
c     $        '(FIXME)'
        dcod1=0.d0
        dcod2=0.d0
      end if
      if(detimy.ne.0.d0) then
        dcod3 = ((1.d0-b22)*(db11*cod3+db12*cod4+ddcod3)
     1          +b12*(db21*cod3+db22*cod4+ddcod4))/detimy
        dcod4 = (b21*(db11*cod3+db12*cod4+ddcod3)
     1          +(1.d0-b11)*(db21*cod3+db22*cod4+ddcod4))/detimy
      else
c        write(*,*)'qdcell @ src/qdcell.f: ',
c     $        '2x2 vetical matrix (I - B) is degenerated!',
c     $        '(FIXME)'
        dcod3=0.d0
        dcod4=0.d0
      end if
c..........going back to 4*4 world
      dtwiss(mfitdx)  = s11*dcod1+s12*dcod2-s13*dcod3-s14*dcod4
      dtwiss(mfitdpx) = s21*dcod1+s22*dcod2-s23*dcod3-s24*dcod4
      dtwiss(mfitdy)  =-s31*dcod1-s32*dcod2+s33*dcod3+s34*dcod4
      dtwiss(mfitdpy) =-s41*dcod1-s42*dcod2+s43*dcod3+s44*dcod4
C
C====== Get derivative of EX,EY,EPX,EPY(dispersion)  ======
c
c(Note) dispersion is defined in 2*2 world.
c..........getting a*3,b*3 and their derivative
      da13 = dtwiss(7)
      da23 = dtwiss(8)
      db13 = dtwiss(9)
      db23 = dtwiss(10)
      dtwiss(mfitpex) =da13
      dtwiss(mfitpepx)=da23
      dtwiss(mfitpey) =db13
      dtwiss(mfitpepy)=db23
      ep1=s11*utwiss(7,idp,l)+s12*utwiss(8,idp,l)
     1   -s13*utwiss(9,idp,l)-s14*utwiss(10,idp,l)
      ep2=s21*utwiss(7,idp,l)+s12*utwiss(8,idp,l)
     1   -s23*utwiss(9,idp,l)-s24*utwiss(10,idp,l)
      ep3=-s31*utwiss(7,idp,l)-s32*utwiss(8,idp,l)
     1    +s33*utwiss(9,idp,l)+s34*utwiss(10,idp,l)
      ep4=-s41*utwiss(7,idp,l)-s42*utwiss(8,idp,l)
     1    +s43*utwiss(9,idp,l)+s44*utwiss(10,idp,l)
      if( cosmux.ne.1.d0 ) then
       dtwiss(7) = (
     1               +0.5d0*(da13+a12*da23-a22*da13) )
     1             /(1.d0-cosmux)
     1     +ds11*ep1+ds12*ep2+ds13*ep3+ds14*ep4
       dtwiss(8) = (
     1               +0.5d0*(da23+a21*da13-a11*da23) )
     1             /(1.d0-cosmux)
     1     +ds21*ep1+ds22*ep2+ds23*ep3+ds24*ep4
      end if
      if( cosmuy.ne.1.d0 ) then
       dtwiss(9) = (
     1               +0.5d0*(db13+b12*db23-b22*db13) )
     1             /(1.d0-cosmuy)
     1     +ds31*ep1+ds32*ep2+ds33*ep3+ds34*ep4
       dtwiss(10)= (
     1               +0.5d0*(db23+b21*db13-b11*db23) )
     1             /(1.d0-cosmuy)
     1     +ds41*ep1+ds42*ep2+ds43*ep3+ds44*ep4
      end if
c-deb
c     print *,'dex,depx=',dtwiss(7),dtwiss(8)
c
      return
      end
