      subroutine qdmdia(tm11,tm12,tm13,tm14,
     1                  tm21,tm22,tm23,tm24,
     1                  tm31,tm32,tm33,tm34,
     1                  tm41,tm42,tm43,tm44,
     1                  dtm11,dtm12,dtm13,dtm14,
     1                  dtm21,dtm22,dtm23,dtm24,
     1                  dtm31,dtm32,dtm33,dtm34,
     1                  dtm41,dtm42,dtm43,dtm44,
     1                  dr1,dr2,dr3,dr4,r1)
c
c----This subroutine finds a derivative of transformation matrix---
c      which diagonalizes 4*4 transfer matrix.
c               M = S-1.Md.S
c         where
c               Md = (A, 0)
c                    (0, B)
c
c      tm** ; transfer matrix( input )
c      dtm**; derivative of transfer matrix( input )
c      dr*  ; derivative of transformation matrix( output )
c               S = (c, J.(r)t.J)
c                   (r,    c    )
c               dr* = d(r*)/dk
c---------------------------------------------------
      implicit none
      real*8 tm11,tm12,tm13,tm14,
     1     tm21,tm22,tm23,tm24,
     1     tm31,tm32,tm33,tm34,
     1     tm41,tm42,tm43,tm44,
     1     dtm11,dtm12,dtm13,dtm14,
     1     dtm21,dtm22,dtm23,dtm24,
     1     dtm31,dtm32,dtm33,dtm34,
     1     dtm41,dtm42,dtm43,dtm44,
     1     dr1,dr2,dr3,dr4,r1,
     $     trpq,dtrpq,trab2,dtrab2,c2,c,dc,da,aa,
     $     dni11,dni12,dni21,dni22,
     $     detn,ddetn,trqk,dtrqk,trqj,dtrqj,
     $     trqkj2,dtrqkj,aaa,bbb,daa,dbb,ccc,dcc,
     $     trql,dtrql,yy11,yy12,yy21,yy22,
     $     dx11,dx12,dx21,dx22,xx11,xx12,xx21,xx22,
     $     dy11,dy12,dy21,dy22,trp,p1,p2,p3,p4,
     $     dtrp,dp1,dp2,dp3,dp4
      real*8 , parameter :: eps=1.d-7
c
      trpq = (tm11 + tm22 - tm33 - tm44)
c.....trpq =(TrP-TrQ)
      dtrpq = (dtm11 + dtm22 - dtm33 - dtm44)
c.....dtrpq =d(TrP-TrQ)/dk
      trab2 = 4.d0*(
     $     (tm31+tm24)*(tm42+tm13)+(tm14-tm32)*(tm41-tm23))
     1     +trpq**2
c      trab2 = 4.d0*(tm31*tm42-tm32*tm41
c     $     +tm13*tm24-tm14*tm23
c     1     +tm31*tm13+tm32*tm23
c     1     +tm41*tm14+tm42*tm24)
c     1     +trpq*trpq
c.....trab2=(TrA-TrB)**2 should be positive or zero.
c
      dtrab2 = 4.d0*(
     $     (dtm31+dtm24)*(tm42+tm13)+(dtm14-dtm32)*(tm41-tm23)
     $     +(tm31+tm24)*(dtm42+dtm13)+(tm14-tm32)*(dtm41-dtm23))
     $     +2.d0*trpq*dtrpq
c      dtrab2 = 8.d0*(dtm31*tm42-dtm32*tm41
c     1             + tm31*dtm42-tm32*dtm41)
c     1        +4.d0*(dtm31*tm13+dtm32*tm23
c     1              +dtm41*tm14+dtm42*tm24
c     1              +tm31*dtm13+tm32*dtm23
c     1              +tm41*dtm14+tm42*dtm24)
c     1        +2.d0*trpq*dtrpq
c      write(*,*)'qdmdia ',8.d0*(dtm31*tm42-dtm32*tm41
c     1             + tm31*dtm42-tm32*dtm41),
c     1        +4.d0*(dtm31*tm13+dtm32*tm23
c     1              +dtm41*tm14+dtm42*tm24
c     1              +tm31*dtm13+tm32*dtm23
c     1              +tm41*dtm14+tm42*dtm24),
c     1        +2.d0*trpq*dtrpq
c.....dtrab2=d((TrA-TrB)**2)/dk
c
c
      if( trab2.lt.eps ) then
       if( abs(dtrab2).lt.eps ) go to 100
       trab2 = eps
      end if
c
c    +++++++++++++++++++++
c==> +  (TrA-TrB)^2 > 0  +
c    +++++++++++++++++++++
c-----Normal case
c     :::::::::::
c         print *,'----qdcell(normal case)-----'
c         print *,' ___tm=                    '
c         print *,      tm11,tm12,tm13,tm14
c         print *,      tm21,tm22,tm23,tm24
c         print *,      tm31,tm32,tm33,tm34
c         print *,      tm41,tm42,tm43,tm44
c         print *,' ___dtm=                    '
c         print *,      dtm11,dtm12,dtm13,dtm14
c         print *,      dtm21,dtm22,dtm23,dtm24
c         print *,      dtm31,dtm32,dtm33,dtm34
c         print *,      dtm41,dtm42,dtm43,dtm44
c         print *,' ___________________________'
          c2 = 0.5d0*abs(trpq)/sqrt(trab2)
          c  = sqrt(0.5d0+c2)
          dc= 0.5d0*c2*(dtrpq/trpq-.5d0*dtrab2/trab2)/c
          if( r1*(tm31+tm24) .le. 0.d0 ) then
            aa = -1.d0/c/sqrt(trab2)
          else
            aa =  1.d0/c/sqrt(trab2)
          end if
          da  = aa*(-dc/c -0.5d0*dtrab2/trab2 )
          dr1 = da*(tm31+tm24) + aa*(dtm31+dtm24)
          dr2 = da*(tm32-tm14) + aa*(dtm32-dtm14)
          dr3 = da*(tm41-tm23) + aa*(dtm41-dtm23)
          dr4 = da*(tm42+tm13) + aa*(dtm42+dtm13)
c          write(*,'(a,1p8g11.3)')'qdmdia ',
c     $         aa,da,c,dc,trab2,dtrab2,trpq,dtrpq
c         print *,' da,aa=',da,aa
c         print *,(tm31+tm24),(dtm31+dtm24)
c         print *,(tm32-tm14),(dtm32-dtm14)
c         print *,(tm41-tm23),(dtm41-dtm23)
c         print *,(tm42+tm13),(dtm42+dtm13)
c         print *,'c2,trpq=',c2,trpq
c         print *,dr1,dr2,dr3,dr4
      return
c
 100  continue
c    +++++++++++++++++++++
c==> +  (TrA-TrB)^2 = 0  +
c    +++++++++++++++++++++
c              This necessarily leads TrP = TrQ.
c
      detn = tm13*tm24-tm14*tm23
      ddetn= dtm13*tm24+tm13*dtm24-dtm14*tm23-tm14*dtm23
      if( abs(detn).lt.eps ) then
        if( abs(ddetn).lt.1.d0 ) go to 800
        detn = sign(eps,detn)
      end if
c         -----------------
c-------> +  det(N) =/ 0  +
c         -----------------
      trqk = tm33 - tm44
      dtrqk= dtm33 - dtm44
      if( abs(trqk).lt.eps ) then
        if( abs(dtrqk).lt.eps ) go to 400
        trqk = sign(eps,trqk)
      end if
c-----------Case of Tr(QK) =/ 0
      trqj = tm43 - tm34
      dtrqj= dtm43 - dtm34
      trqkj2= (trqk+trqj)*(trqk-trqj)
      dtrqkj= (dtrqk+dtrqj)*(trqk-trqj)+(trqk+trqj)*(dtrqk-dtrqj)
      if( abs(trqkj2).lt.eps ) then
        if( abs(dtrqkj).lt.eps ) go to 200
        trqkj2 = eps
      end if
c...........Case of (TrQK)**2 =/ (TrQJ)**2
              bbb = -detn/(trqj+trqk)
              dbb = -((dtrqj+dtrqk)*bbb+ddetn)/(trqj+trqk)
            go to 300
 200  continue
c...........Case of (TrQK)**2 = (TrQJ)**2
            if( abs(trqj).lt.eps ) trqj=sign(eps,trqj)
c           [ Tr(QJ) =/ 0 ]
              bbb = -0.5d0*detn/trqj
              dbb = -(dtrqj*bbb+0.5d0*ddetn)/trqj
 300  continue
            aaa  = (-detn-bbb*trqj )/trqk
            daa  = (-dtrqk*aaa-ddetn-dbb*trqj-bbb*dtrqj)/trqk
            ccc  = 0.d0
            dcc  = 0.d0
            go to 700
c-----------Case of Tr(QK) = 0
 400  continue
      trqj = tm43 - tm34
      dtrqj= dtm43 - dtm34
      if( abs(trqj).lt.eps ) then
        if( abs(dtrqj).lt.eps ) go to 500
        trqj = sign(eps,trqj)
      end if
c...........Case of TrQJ =/ 0
            aaa = detn/trqj
            bbb =-detn/trqj
            ccc = 0.d0
            daa = (-dtrqj*aaa+ddetn)/trqj
            dbb =-(dtrqj*bbb+ddetn)/trqj
            dcc = 0.d0
            go to 700
 500  continue
c...........Case of TrQJ = 0
      trql = tm34 + tm43
      dtrql= dtm34 + dtm43
      if( abs(trql).lt.eps ) then
        if( abs(dtrql).lt.eps ) go to 600
        trql = sign(eps,trql)
      end if
c           [ Tr(QL) =/ 0 ]
            aaa = 0.d0
            bbb = detn/trql
            ccc =-bbb
            daa = 0.d0
            dbb = (-dtrql*bbb+ddetn)/trql
            dcc =-dbb
            go to 700
 600  continue
c           [ Tr(QL) =  0 ]
      if( detn.gt.0.d0 ) then
c           ( Det(N) > 0 )
            aaa = 0.d0
            bbb = sqrt(detn/2.d0)
            ccc = 0.d0
            daa = 0.d0
            dbb = 0.5d0*ddetn/sqrt(2.d0*detn)
            dcc = 0.d0
      else
c           ( Det(N) < 0 )
            aaa = sqrt(-detn/2.d0)
            bbb = 0.d0
            ccc = 0.d0
            daa =-0.5d0*ddetn/sqrt(-2.d0*detn)
            dbb = 0.d0
            dcc = 0.d0
      end if
 700  continue
      xx11 =  aaa
      xx12 =  bbb + ccc
      xx21 = -bbb + ccc
      xx22 = -aaa
      dx11 =  daa
      dx12 =  dbb + dcc
      dx21 = -dbb + dcc
      dx22 = -daa
      dni11= -ddetn*tm24/detn+dtm24
      dni12=  ddetn*tm14/detn-dtm14
      dni21=  ddetn*tm23/detn-dtm23
      dni22= -ddetn*tm13/detn+dtm13
      dr1  = ( dx11*tm24-dx12*tm23 + xx11*dni11+xx12*dni21 )/detn
      dr2  = (-dx11*tm14+dx12*tm13 + xx11*dni12+xx12*dni22 )/detn
      dr3  = ( dx21*tm24-dx22*tm23 + xx21*dni11+xx22*dni21 )/detn
      dr4  = (-dx21*tm14+dx22*tm13 + xx21*dni12+xx22*dni22 )/detn
      return
 800  continue
c         -----------------
c-------> +  det(N) =  0  +
c         -----------------
      trqk = tm33 - tm44
      trqj = tm43 - tm34
      trql = tm43 + tm34
      dtrqk= dtm33 - dtm44
      dtrqj= dtm43 - dtm34
      dtrql= dtm43 + dtm34
      if( abs(trqk).lt.eps ) go to 1000
 900  continue
c-----------Case of Tr(QK) =/ 0
            aaa = 0.d0
            bbb = 0.d0
            ccc = 1.d0/trqk
            daa = 0.d0
            dbb = 0.d0
            dcc = -dtrqk*ccc/trqk
      go to 1300
1000  continue
      if( abs(trqj).lt.eps ) go to 1100
c-----------Case of Tr(QJ) =/ 0
            aaa = -1.d0/trqj
            bbb = 0.d0
            ccc = 0.d0
            daa = -dtrqj*aaa/trqj
            dbb = 0.d0
            dcc = 0.d0
      go to 1300
1100  continue
      if( abs(trql).lt.eps ) go to 1200
c-----------Case of Tr(QL) =/ 0
            aaa = 0.d0
            bbb = -1.d0/trql
            ccc = 0.d0
            daa = 0.d0
            dbb = -dtrql*bbb/trql
            dcc = 0.d0
      go to 1300
1200  continue
      if( abs(dtrqk).ge.eps ) then
        trqk = sign(eps,trqk)
        go to 900
      else if( abs(dtrqj).ge.eps ) then
        trqj = sign(eps,trqj)
        go to 1000
      else if( abs(dtrql).ge.eps ) then
        trql = sign(eps,trql)
        go to 1100
      else
        go to 1400
      end if
1300  continue
      yy11 = aaa + bbb
      yy12 = ccc
      yy21 = ccc
      yy22 = aaa - bbb
      dy11 = daa + dbb
      dy12 = dcc
      dy21 = dcc
      dy22 = daa - dbb
      dr1  = dy11*(-tm23)+dy12*(-tm24) + yy11*(-dtm23)+yy12*(-dtm24)
      dr2  = dy11*tm13   +dy12*tm14    + yy11*dtm13   +yy12*dtm14
      dr3  = dy21*(-tm23)+dy22*(-tm24) + yy21*(-dtm23)+yy22*(-dtm24)
      dr3  = dy21*tm13   +dy22*tm14    + yy21*dtm13   +yy22*dtm14
      return
1400  continue
c-----------Case of Tr(QJ) = Tr(QL) = Tr(QK) = 0
      trp = tm11 + tm22
      p1  = tm11 - 0.5d0*trp
      p2  = tm12
      p3  = tm21
      p4  = tm22 - 0.5d0*trp
      dtrp= dtm11 + dtm22
      dp1 = dtm11 - 0.5d0*dtrp
      dp2 = dtm12
      dp3 = dtm21
      dp4 = dtm22 - 0.5d0*dtrp
      if( abs(p3).lt.eps ) then
        if( abs(dp3).lt.eps ) go to 1500
        p3 = sign(eps,p3)
      end if
c...........Case of (P-1/2*TrP*I)3 =/ 0
      dr1 = 0.d0
      dr2 = -dp3*tm24/p3/p3+dtm24/p3
      dr3 = 0.d0
      dr4 = dp3*tm23/p3/p3-dtm23/p3
      return
1500  continue
      if( abs(p1).lt.eps ) then
        if( abs(dp1).lt.eps ) then
c...........Case of (P-1/2*TrP*I)3 = 0 and (P-1/2*TrP*I)1 =  0
            dr1 = 0.d0
            dr2 = 0.d0
            dr3 = 0.d0
            dr4 = 0.d0
            return
        end if
        p1 = sign(eps,p1)
      end if
c...........Case of (P-1/2*TrP*I)3 = 0 and (P-1/2*TrP*I)1 =/ 0
      dr1 = -tm24*dp1/p1/p1+dtm24/p1
      dr2 = 0.d0
      dr3 = tm23*dp1/p1/p1-dtm23/p1
      dr4 = 0.d0
      return
      end
