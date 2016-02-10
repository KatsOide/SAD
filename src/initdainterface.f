      subroutine initdainterface(isp1,kx,irtc)
      use tfstk
      implicit none
      integer*8 kx
      integer nvar,nord,lvec
      common /dainfo/nvar,nord,lvec
      integer*4 isp1,irtc,narg,itfmessage
      narg=isp-isp1           !number of arguments, always >=1 .
      if(narg .gt. 1)then
        irtc=itfmessage(9,'General::narg','"1"')
c        irtc=-1               !return without error message,
c                             !leave unevaluated.
        return
      endif
      if(ktfrealq(ktastk(isp)))then     !Check the type of arg
        nord=int(rtastk(isp))
        nvar=6
        call mkmachinef(nord)
        call initda(nord,lvec)     !If real, rtastk(isp1+i)has
c                                  !the i-th value.
        kx = ktfoper + mtfnull
        irtc=0                     !Setthe return code if no error.
        return
      else
        irtc=itfmessage(9,'General::wrongtype','"Real number"')
      endif
      return
      end
c
      subroutine execdapackage(isp1,kx,irtc)
      use tfstk
      implicit none
      integer*8 kx
      integer nvar,nord,lvec
      common /dainfo/nvar,nord,lvec
      integer*4 isp1,irtc,narg,itfmessage
      narg=isp-isp1           !number of arguments, always >=1 .
      if(narg .gt. 1)then
        irtc=itfmessage(9,'General::narg','"1"')
c        irtc=-1               !return without error message,
c                             !leave unevaluated.
        return
      endif
      if(ktfrealq(ktastk(isp)))then     !Check the type of arg
        nord=int(rtastk(isp))
        nvar=6
        call mkmachinef(nord)
        call execda(nord)     !If real, rtastk(isp1+i) 
c                                 has the i-th value.
        kx = ktfoper + mtfnull
        irtc=0                     !Setthe return code if no error.
        return
      else
        irtc=itfmessage(9,'General::wrongtype','"Real number"')
      endif
      return
      end

c
c =============================================================
c
      subroutine mkmachinef(nord)
      use tfstk
      use ffs
      use temw, only: rr=>r, rri=>ri
      use tffitcode
      implicit real*8 (a-h,o-z)
c
c  ilist(1,ilattp+1)=latt(1,1)
c  ilist(1,ifiele)=iele(1)
c  ilist(1,ifkpel)=kpele(1)
c
      data dx,dy,dteta/3*0./
      data il,ib,iq,is,icv/5*0/
c     aa=45.d0/datan(1.d0)
      integer ent_edge,exit_edge
      logical preele,sucele
      aa=1.
      lfile=11
      open(lfile,FILE='Machine.file',STATUS='UNKNOWN')
      write(lfile,*) 'Operator Koiso=(Nord=',nord,' );'
      write(lfile,*) 'Beam BEAM=(energy=',pgev*1E-9,
     &   ', N_particle=',pbunch,' );'
c      write(*,'(1p,6E12.4)') rr
      write(lfile,*) 'Linmap R='
      do 9 i=1,6
         if(i.eq.1) then
            write(lfile,*) '({{',rri(i,1),',',rri(i,2),',',rri(i,3),','
     &           ,rri(i,4),',',rri(i,5),',',rri(i,6),'},'
         else if(i.eq.6) then
            write(lfile,*) '{',rri(i,1),',',rri(i,2),',',rri(i,3),','
     &           ,rri(i,4),',',rri(i,5),',',rri(i,6),'}});'
         else
            write(lfile,*) '{',rri(i,1),',',rri(i,2),',',rri(i,3),','
     &           ,rri(i,4),',',rri(i,5),',',rri(i,6),'},'
         endif
 9    continue
c
      write(*,'(//A/)')
     &     '***** Create Machine file named "Machine.file" *****'
      write(*,*) ' nlat, nele ',nlat,nele
      do 10 i=1,nlat-1
c         j=latt(1,i)
         j=ilist(1,ilattp+i)
         id=idtype(j)
         rle=rlist(idval(j)+1)
c
c         if(id.ge.2.and.id.le.40) then
         jpre=0
         ipre=i-1
 1       continue
         if(ipre.le.0) then
            ipre=1
            jpre=ilist(1,ilattp+ipre)
            goto 2
         endif
         jpre=ilist(1,ilattp+ipre)
         if(idtype(jpre).gt.40) then
            ipre=ipre-1
            goto 1
         endif
 2       jsuc=0
         isuc=i+1
 3       continue
         if(isuc.gt.nlat-1) then
            isuc=nlat-1
            jsuc=ilist(1,ilattp+nlat-1)
            goto 4
         endif
         jsuc=ilist(1,ilattp+isuc)
         if(idtype(jsuc).gt.40) then
            isuc=isuc+1
            goto 3
         endif
 4       continue
         idpre=idtype(jpre)
         idsuc=idtype(jsuc)
         if(idpre.le.0.or.idpre.gt.40) then
            idpre=id
            write(*,*) 'idpre is set id which is',id,i,'-th element'
         endif
         if(idsuc.le.0.or.idsuc.gt.40) then
            idsuc=id
            write(*,*) 'idsuc is set id which is',id,i,'-th element'
         endif
c
c         if(i.lt.nlat) jsuc=ilist(1,ilattp+i+1)
c         write(6,*) i,j

         if(id.eq.1) then
            il=il+1
            write(lfile,*) 'Drift ',pname(j),'( L=',rle,');'
c
c -----  Bending magnet  ----------------------------------
c
         else if(id.eq.2) then
            ib=ib+1
            ang=rlist(ilist(2,ilattp+i)+2)
            e1=rlist(ilist(2,ilattp+i)+3)*ang*aa
            e2=rlist(ilist(2,ilattp+i)+4)*ang*aa
            teta=rlist(idval(j)+5)*aa
            dteta=rlist(ilist(2,ilattp+i)+5)*aa-teta
            dx=rlist(ilist(2,ilattp+i)+9)
            dy=rlist(ilist(2,ilattp+i)+10)
            disfrin=rlist(ilist(2,ilattp+i)+12)
            rhoi=ang/rle
            if(disfrin.eq.0.) then
            ent_edge=1
            exit_edge=1
            preele=idpre.ge.2.and.idpre.le.40
            sucele=idsuc.ge.2.and.idsuc.le.40
            if(preele.
     +              and.rlist(ilist(2,ilattp+ipre)).ne.0.) ent_edge=0
            if(sucele.
     +              and.rlist(ilist(2,ilattp+isuc)).ne.0.) exit_edge=0
            else
               ent_edge=0
               exit_edge=0
            endif
            write(lfile,*) 'Bend ',pname(j),'( L=',rle,', phi=',ang,
     +         ', e1=',e1,', e2=',e2,', ent_edge=',ent_edge,
     &         ', exit_edge=',exit_edge,' ) '
            if(dx.ne.0..or.dy.ne.0..or.teta.ne.0.) then
               write(lfile,*) '  dx=',dx,', dy=',dy,
     &              ', dtheta=',teta
            endif
            write(lfile,*) ';'
c
c -----  Quadrupole magnet  ----------------------------------
c
         else if(id.eq.4) then
            iq=iq+1
            rkl1=rlist(ilist(2,ilattp+i)+2)
            teta=rlist(idval(j)+4)*aa
            dx=rlist(ilist(2,ilattp+i)+5)
            dy=rlist(ilist(2,ilattp+i)+6)
            dteta=rlist(ilist(2,ilattp+i)+4)*aa-teta
            disfrin=rlist(ilist(2,ilattp+i)+9)
            e1=0.
            e2=0.
            if(disfrin.eq.0.) then
            ent_edge=1
            exit_edge=1
            preele=idpre.ge.2.and.idpre.le.40
            sucele=idsuc.ge.2.and.idsuc.le.40
            if(preele.
     +              and.rlist(ilist(2,ilattp+ipre)).ne.0.) ent_edge=0
            if(sucele.
     +              and.rlist(ilist(2,ilattp+isuc)).ne.0.) exit_edge=0
            else
               ent_edge=0
               exit_edge=0
            endif
            write(lfile,*) 'Quad ',pname(j),'( L=',rle,', K1=',rkl1,
     +          ', Ndiv=3, ent_edge=',ent_edge,', exit_edge=',
     +           exit_edge,' ) '
            if(dx.ne.0..or.dy.ne.0..or.teta.ne.0.) then
               write(lfile,*) '  dx=',dx,', dy=',dy,
     &              ', dtheta=',teta
            endif
            write(lfile,*) ';'
c
c -----  Sextupole magnet  ----------------------------------
c
         else if(id.eq.6) then
            is=is+1
            rkl2=rlist(ilist(2,ilattp+i)+2)
            teta=rlist(idval(j)+4)*aa
            dx=rlist(ilist(2,ilattp+i)+5)
            dy=rlist(ilist(2,ilattp+i)+6)*aa-teta
            dteta=rlist(ilist(2,ilattp+i)+4)
            e1=0.
            e2=0.
            if(rle.gt.0.) then
c              write(lfile,*) dx,dy,teta,dteta
               write(lfile,*) 'Sext ',pname(j),'(L=',rle,
     +              ', K2=',rkl2,' )'
            else
c              write(lfile,*) dx,dy,teta,dteta
               write(lfile,*) 'Thin ',pname(j),
     +             '(n_MP=2, K_n=',rkl2,' ) '
            endif
            if(dx.ne.0..or.dy.ne.0..or.teta.ne.0.) then
               write(lfile,*) '  dx=',dx,', dy=',dy,
     &              ', dtheta=',teta
            endif
            write(lfile,*) ';'
c
c -----  Octupole magnet  ----------------------------------
c
         else if(id.eq.8) then
            is=is+1
            rkl2=rlist(ilist(2,ilattp+i)+2)
            teta=rlist(idval(j)+4)*aa
            dx=rlist(ilist(2,ilattp+i)+5)
            dy=rlist(ilist(2,ilattp+i)+6)*aa-teta
            dteta=rlist(ilist(2,ilattp+i)+4)
            e1=0.
            e2=0.
            if(rle.gt.0.) then
               write(lfile,*) 'Oct ',pname(j),'(L=',rle,
     +              ', K3=',rkl2,' )'
            else
c              write(lfile,*) dx,dy,teta,dteta
               write(lfile,*) 'Thin ',pname(j),
     +             '(n_MP=3, K_n=',rkl2,' )'
            endif
            if(dx.ne.0..or.dy.ne.0..or.teta.ne.0.) then
               write(lfile,*) '  dx=',dx,', dy=',dy,
     &              ', dtheta=',teta
            endif
            write(lfile,*) ';'
c
c -----  Cavity  ----------------------------------
c
         else if(id.eq.31) then
            icv=icv+1
            rfv=rlist(ilist(2,ilattp+i)+2)
            rharm=rlist(ilist(2,ilattp+i)+3)
            freq=rlist(ilist(2,ilattp+i)+5)
            write(lfile,*) 'Cavity ',pname(j),'(volt=',rfv,
     +           ', freq=',freq,' );'
c
c -----  Beam-beam collision  ----------------------------------
c
         else if(id.eq.36) then
            i_pin0=ilist(2,ilattp+i)
            write(lfile,*) 'IP ',pname(j),
     &           '( N_particle=',rlist(i_pin0+29),
     &           ', N_slice=',rlist(i_pin0+28),
     &           ', x_angle=',rlist(i_pin0+21),','
            write(lfile,*) '      ',
     &           '  ax=',rlist(i_pin0+1),
     &           ', bx=',rlist(i_pin0+2),
     &           ', ay=',rlist(i_pin0+3),
     &           ', by=',rlist(i_pin0+4),','
            write(lfile,*) '      ',
     &           '  ex=',rlist(i_pin0+9),
     &           ', epx=',rlist(i_pin0+10),
     &           ', ey=',rlist(i_pin0+11),
     &           ', epy=',rlist(i_pin0+12),','
            write(lfile,*) '      ',
     &           '  zx=',rlist(i_pin0+13),
     &           ', zpx=',rlist(i_pin0+14),
     &           ', zy=',rlist(i_pin0+15),
     &           ', zpy=',rlist(i_pin0+16),','
            write(lfile,*) '      ',
     &           '  emx=',rlist(i_pin0+22),
     &           ', emy=',rlist(i_pin0+23),
     &           ', emz=',rlist(i_pin0+24)*rlist(i_pin0+27),
     &           ', bz=',rlist(i_pin0+27)/rlist(i_pin0+24),
     &           ' );'
c
c -----  Phase space rotator  ----------------------------------
c
         else if(id.eq.37) then
            i_pin0=ilist(2,ilattp+i)
            write(lfile,*) 'R_ph_rot ',pname(j),
     &           '( nux=',rlist(i_pin0+3),
     &           ', nuy=',rlist(i_pin0+6),
     &           ', nuz=',rlist(i_pin0+25),','
            write(lfile,*) '      ',
     &           '  ax=',rlist(i_pin0+1),
     &           ', bx=',rlist(i_pin0+2),
     &           ', ay=',rlist(i_pin0+4),
     &           ', by=',rlist(i_pin0+5),','
            write(lfile,*) '      ',
     &           '  ex=',rlist(i_pin0+11),
     &           ', epx=',rlist(i_pin0+12),
     &           ', ey=',rlist(i_pin0+13),
     &           ', epy=',rlist(i_pin0+14),','
            write(lfile,*) '      ',
     &           '  zx=',rlist(i_pin0+15),
     &           ', zpx=',rlist(i_pin0+16),
     &           ', zy=',rlist(i_pin0+17),
     &           ', zpy=',rlist(i_pin0+18),','
            write(lfile,*) '      ',
     &           '  az=',rlist(i_pin0+22),
     &           ', bz=',rlist(i_pin0+27),
     &           ' );'
c
c -----  Marker  ----------------------------------
c
         else if(id.eq.41) then
            imk=imk+1
c
c -----  Monitor  ----------------------------------
c
         else if(id.eq.42) then
            imon=imon+1
         else if(id.eq.0) then
            write(*,*) 'Element type=0   i=',i
         else
            write(lfile,'(A)') '//This element is not translated',id,
     &           pname(j)
         end if
 10   continue
c
      close(lfile)
      return
c
      end
