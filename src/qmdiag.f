      subroutine qmdiag(tm11,tm12,tm13,tm14,
     1     tm21,tm22,tm23,tm24,
     1     tm31,tm32,tm33,tm34,
     1     tm41,tm42,tm43,tm44,
     1     r1,r2,r3,r4,c,stab,nanq,lfno)
c     
c---- This subroutine find a transformation matrix---
c     which diagonalizes 4*4 transfer matrix.
c     M = S-1.Md.S
c     where
c     Md = (A, 0)
c     (0, B)
c     
c     tm** ; transfer matrix( input )
c     r*   ; transformation matrix( output )
c     S = (c, J.(r)t.J)
c     (r,    c    )
c---------------------------------------------------
      use tfstk, only: ktfenanq
      implicit none
      real*8 , intent(in)::tm11,tm12,tm13,tm14,
     1     tm21,tm22,tm23,tm24,
     1     tm31,tm32,tm33,tm34,
     1     tm41,tm42,tm43,tm44
      real*8 , intent(out)::r1,r2,r3,r4,c
      real*8 trpq,trab2,c2,aa,
     $     detn,trqk,trqj,aaa,bbb,ccc,p1,p2,p3,p4,trp,trqkj2,
     $     trql,xx11,xx12,xx21,xx22,yy11,yy12,yy21,yy22
      integer*4 , intent(inout)::lfno
      logical*4, intent(out):: stab,nanq
      stab=.true.
      nanq=.false.
c     
      trpq = tm11 + tm22 - tm33 - tm44
c.....trpq =(TrP-TrQ)
      trab2 = 4.d0*(
     $     (tm31+tm24)*(tm42+tm13)+(tm14-tm32)*(tm41-tm23))
     1     +trpq**2
c      trab2 = 8.d0*(tm31*tm42-tm32*tm41)
c     1     +4.d0*(tm31*tm13+tm32*tm23
c     1     +tm41*tm14+tm42*tm24)
c     1     +trpq*trpq
c.....trab2=(TrA-TrB)**2 should be positive or zero.
c     ::::::::::::::::
C--   deb
c      print *,'===qmdiag=== trab2 =',trab2
      if(ktfenanq(trab2))then
        write(lfno,'(a,1p5g15.7)')
     $       '***qmdiag---> Matrix not a number: ',trab2,
     $       tm11,tm22,tm33,tm44
        lfno=0
        stab=.false.
        nanq=.true.
        return
      endif
      if( trab2.ge.1.d-14 ) then
c     +++++++++++++++++++++
c==   > +  (TrA-TrB)^2 > 0  +
c     +++++++++++++++++++++
c     
c-----Normal case
c     :::::::::::
C--   deb
c     print *,'---> trab2>0'
        c2 = 0.5d0*(1.d0 + abs(trpq)/sqrt(trab2))
        c  = sqrt(c2)
        if( trpq.ge.0.d0 ) then
          aa = -1.d0/c/sqrt(trab2)
        else
          aa =  1.d0/c/sqrt(trab2)
        end if
        r1 = aa*(tm31 + tm24)
        r2 = aa*(tm32 - tm14)
        r3 = aa*(tm41 - tm23)
        r4 = aa*(tm42 + tm13)
      else if( abs(trab2).le.1.d-14 ) then
c     +++++++++++++++++++++
c==   > +  (TrA-TrB)^2 = 0  +
c     +++++++++++++++++++++
c     This necessarily leads TrP = TrQ.
c     
C--   deb
c     print *,'---> trpq=0 ,trab2=0'
        detn = tm13*tm24-tm14*tm23
        if( abs(detn).ge.1.d-14 ) then
c     -----------------
c-------> +  det(N) =/ 0  +
c     -----------------
          trqk = tm33 - tm44
          if( abs(trqk).ge.1.d-7 ) then
c-----------Case of Tr(QK) =/ 0
            trqj = tm43 - tm34
            trqkj2 = (trqk+trqj)*(trqk-trqj)
            if( abs(trqkj2).ge.1.d-14 ) then
c...........Case of (TrQK)**2 =/ (TrQJ)**2
              c   = 1.d0
              bbb = -detn/(trqj+trqk)
              ccc = 0.d0
            else
c...........Case of (TrQK)**2 = (TrQJ)**2
              if( abs(trqj).ge.1.d-7 ) then
c     [ Tr(QJ) =/ 0 ]
c     ..This condition must be always satisfied
c     because Tr(QK)=/0 and Tr(QK)**2=Tr(QJ)**2.
                c   = 1.d0
                bbb = -0.5d0*detn/trqj
                ccc = 0.d0
              else
                print *,'+++QMDIAG Logics is broken. Check tm**'
                print *,tm11,tm12,tm13,tm14
                print *,tm21,tm22,tm23,tm24
                print *,tm31,tm32,tm33,tm34
                print *,tm41,tm42,tm43,tm44
                print *,'trqk,trqj,trqkj2=',trqk,trqk,trqkj2
c     begin initialize for preventing compiler warning
                bbb = 0.d0
c     end   initialize for preventing compiler warning
              end if
            end if
            aaa  = ( detn*(1.d0-2.d0*c*c)/c - bbb*trqj )/trqk
          else
c-----------Case of Tr(QK) = 0
            trqj = tm43 - tm34
            if( abs(trqj).ge.1.d-7 ) then
c...........Case of TrQJ =/ 0
              c   = 1.d0
              aaa = detn/trqj
              bbb =-detn/trqj
              ccc = 0.d0
            else
c...........Case of TrQJ = 0
              trql = tm34 + tm43
              if( abs(trql).ge.1.d-7 ) then
c     [ Tr(QL) =/ 0 ]
                c   = 1.d0
                aaa = 0.d0
                bbb = detn/trql
                ccc =-detn/trql
              else
c     [ Tr(QL) =  0 ]
                c   = 1.d0/sqrt(2.d0)
                if( detn.gt.0.d0 ) then
c     ( Det(N) > 0 )
                  aaa = 0.d0
                  bbb = sqrt(detn/2.d0)
                  ccc = 0.d0
                else
c     ( Det(N) < 0 )
                  aaa = sqrt(-detn/2.d0)
                  bbb = 0.d0
                  ccc = 0.d0
                end if
              end if
            end if
          end if
          xx11 =  aaa
          xx12 =  bbb + ccc
          xx21 = -bbb + ccc
          xx22 = -aaa
          r1   = ( xx11*tm24-xx12*tm23)/detn
          r2   = (-xx11*tm14+xx12*tm13)/detn
          r3   = ( xx21*tm24-xx22*tm23)/detn
          r4   = (-xx21*tm14+xx22*tm13)/detn
        else
c     -----------------
c-------> +  det(N) =  0  +
c     -----------------
c     print *,'---detn=0,detn=',detn
          trqk = tm33 - tm44
          trqj = tm43 - tm34
          trql = tm43 + tm34
          if( abs(trqk).ge.1.d-7 ) then
c-----------Case of Tr(QK) =/ 0
            aaa = 0.d0
            bbb = 0.d0
            ccc = 1.d0/trqk
          else if( abs(trqj).ge.1.d-7 ) then
c-----------Case of Tr(QJ) =/ 0
            aaa = -1.d0/trqj
            bbb = 0.d0
            ccc = 0.d0
          else if( abs(trql).ge.1.d-7 ) then
c-----------Case of Tr(QL) =/ 0
            aaa = 0.d0
            bbb = -1.d0/trql
            ccc = 0.d0
          else
            go to 1000
          end if
          c    = 1.d0
          yy11 = aaa + bbb
          yy12 = ccc
          yy21 = ccc
          yy22 = aaa - bbb
          r1   = yy11*(-tm23) + yy12*(-tm24)
          r2   = yy11*tm13    + yy12*tm14
          r3   = yy21*(-tm23) + yy22*(-tm24)
          r4   = yy21*tm13    + yy22*tm14
          go to 2000
 1000     continue
c-----------Case of Tr(QJ) = Tr(QL) = Tr(QK) = 0
          trp = tm11 + tm22
          p1  = tm11 - 0.5d0*trp
          p2  = tm12
          p3  = tm21
          p4  = tm22 - 0.5d0*trp
          if( abs(p3).ge.1.d-7 ) then
c...........Case of (P-1/2*TrP*I)3 =/ 0
            c  = 1.d0
            r1 = 0.d0
            r2 = tm24/p3
            r3 = 0.d0
            r4 = -tm23/p3
          else if( abs(p1).ge.1.d-7 ) then
c...........Case of (P-1/2*TrP*I)3 = 0 and (P-1/2*TrP*I)1 =/ 0
            c  = 1.d0
            r1 = tm24/p1
            r2 = 0.d0
            r3 = -tm23/p1
            r4 = 0.d0
          else
c...........Case of (P-1/2*TrP*I)3 = 0 and (P-1/2*TrP*I)1 =  0
            c  = 1.d0/sqrt(2.d0)
            r1 = 1.d0
            r2 = 0.d0
            r3 = 0.d0
            r4 = c
          end if
 2000     continue
        end if
      else
c     ?????????????????????
c==   > ?  (TrA-TrB)^2 < 0  ?
c     ?????????????????????
c     
c     print *,'+++QMDIAG Strange transfer matrix !!!!!'
c     print *,tm11,tm12,tm13,tm14
c     print *,tm21,tm22,tm23,tm24
c     print *,tm31,tm32,tm33,tm34
c     print *,tm41,tm42,tm43,tm44
        if(lfno .gt. 0)then
          write(lfno,*)
     $   '***qmdiag---> Sum resonance: (TrA-TrB)^2 =',trab2
          lfno=0
        endif
        stab=.false.
      endif
      return
      end
