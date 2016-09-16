      subroutine ptol(word,latt,mult,master,
     1                twiss,gammab,idp,dp,istr,nstr,lfno)
      use tfstk
      use ffs
c     Tor Raubenheimer's analytical estimation (SLAC-PUB-4937)
      use tffitcode
      use ffs_pointer, only:idelc,idvalc,idtypec
      implicit real*8(a-h,o-z)
      logical over,full
      character*(*) word,vo(4)*8,autofg*8,line*79,name*8
      complex*16 ct,cts,cc
      integer*8 latt(nlat),ip
      dimension twiss(nlat,-ndim:ndim,ntwissfun),gammab(nlat),
     1          istr(nstra,4),mult(nlat)
      dimension twisss(ntwissfun),t(10),tt(10),ct(5),cts(5),cc(5),a(10)
      integer*4 master(nlat)
      data ale/0.05/,ndiv/1/,lpmax/50/
      include 'inc/common.inc'
c
      full=.false.
      lp=0
      call getwdl(word)
      if(word.eq.'FULL') then
        full=.true.
        call getwdl(word)
      endif
      do 10 i=1,10
        t(i)=0d0
        tt(i)=0d0
  10  continue
      do 20 i=1,5
        ct(i)=(0d0,0d0)
        cts(i)=(0d0,0d0)
        cc(i)=(0d0,0d0)
  20  continue
      name=' '
      line=' '
      ts2=0d0
      ts4=0d0
      psix=twiss(nlat,idp,3)-twiss(1,idp,3)
      psiy=twiss(nlat,idp,6)-twiss(1,idp,6)
      psix2=0.5d0*psix
      psiy2=0.5d0*psiy
      a(1)=1d0
      a(3)=1d0/4d0
      a(5)=1d0/8d0
      a(6)=1d0/4d0
      a(7)=1d0/16d0
      a(8)=1d0/32d0
      a(9)=1d0/32d0
      ndivt=0
      vs=0.d0
      vsb=0.d0
      if(full) then
        write(lfno,'(2a)')'           <-------- Ey**2/By ------->',
     1                          '<----- Emiy/Emix (LC) ---->'
        write(lfno,'(2a)')' Element   rot_q**2 ym_s**2  yc**2/By ',
     1                          'rot_q**2 ym_s**2  yc**2/By '
      endif
      do 100 j=0,nstr
        if(j.eq.0) then
          is=1
          if(nstr.eq.0) then
            if=nlat-1
          else
            if=max(1,istr(istr(1,2),1)-1)
          endif
        elseif(j.eq.nstr) then
          is=istr(istr(nstr,2),1)
          if=nlat-1
        else
          is=istr(istr(j,2),1)
          if=istr(istr(j+1,2),1)-1
        endif
        ida=0
        do 90 i=is,if
          if(idtypec(i).eq.icquad .or.
     1      idtypec(i).eq.icbend .or.
     1      idtypec(i).eq.icsext ) then
            id=idtypec(i)
            if(id.eq.icbend.and.
     $           rlist(idvalc(i)+8).eq.0d0)goto 90
            if(master(i).gt.0) then
              ndivt=max(ndivt,1)
              if(ida.eq.icbend .or. ida.eq.icquad) then
c10Oct93 vs was replaced with vv in the case of bend or quad.
c     (vsb=K1 of bend)
                if(ida.eq.icbend) then
                  vv=vsb
                else
                  vv=vs
                endif
                t(1)=t(1)+(vv*tt(1)/ndivt)**2
                t(2)=t(2)+vv*tt(2)/ndivt
                ct(1)=ct(1)+vv*cc(1)/ndivt
                t(6)=t(6)+(vv*tt(6)/ndivt)**2
                if(full) then
                  line(10:)=autofg((vv*tt(1)/ndivt)**2*a(1),'8.5')
                  line(37:)=autofg((vv*tt(6)/ndivt)**2*a(6),'8.5')
                endif
                tt(1)=0d0
                tt(2)=0d0
                cc(1)=(0d0,0d0)
                tt(6)=0d0
                ndivt=0
              elseif(ida.eq.icsext) then
                t(3)=t(3)+(vs*tt(3)/ndivt)**2
                t(4)=t(4)+vs*tt(4)/ndivt
                ct(2)=ct(2)+vs*cc(2)/ndivt
                t(7)=t(7)+(vs*tt(7)/ndivt)**2
                ct(3)=ct(3)+vs*cc(3)/ndivt
                ct(4)=ct(4)+vs*cc(4)/ndivt
                ct(5)=ct(5)+vs*cc(5)/ndivt
                if(full) then
                  line(19:)=autofg((vs*tt(3)/ndivt)**2*a(3),'8.5')
                  line(46:)=autofg((vs*tt(7)/ndivt)**2*a(7),'8.5')
                endif
                tt(3)=0d0
                tt(4)=0d0
                cc(2)=(0d0,0d0)
                tt(7)=0d0
                cc(3)=(0d0,0d0)
                cc(4)=(0d0,0d0)
                cc(5)=(0d0,0d0)
                ndivt=0
              endif
              if(full) then
                if(line(1:8).ne.' ') then
                  write(lfno,'(2a)')' ',line
                  lp=lp+1
                  if(lp.gt.lpmax) then
                    write(lfno,'(2a)')
     1               '           <-------- Ey**2/By ------->',
     1                          '<----- Emiy/Emix (LC) ---->'
                    write(lfno,'(2a)')
     1               ' Element   rot_q**2 ym_s**2  yc**2/By ',
     1                          'rot_q**2 ym_s**2  yc**2/By '
                    lp=0
                  endif
                endif
                call elname(i,name)
                line=name
              endif
              ida=idtypec(i)
            endif
            i1=i+1
            if(id.eq.icquad .or. id.eq.icbend) then
              tt(1)=tt(1)+0.5d0*(sqrt(twiss(i,idp,5))*twiss(i,idp,7)
     1                   +sqrt(twiss(i1,idp,5))*twiss(i1,idp,7))
              tt(2)=tt(2)-0.5d0*(twiss(i,idp,5)+twiss(i1,idp,5))
              cc(1)=cc(1)-0.5d0*(twiss(i,idp,5)*dcmplx(
     1               cos(2d0*twiss(i,idp,6)),sin(2d0*twiss(i,idp,6)))
     1                         +twiss(i1,idp,5)*dcmplx(
     1               cos(2d0*twiss(i1,idp,6)),sin(2d0*twiss(i1,idp,6))
     1                                                 ) )
              tt(6)=tt(6)+0.5d0*(sqrt(twiss(i,idp,2)*twiss(i,idp,5))
     1                         +sqrt(twiss(i1,idp,2)*twiss(i1,idp,5)))
            elseif(id.eq.icsext) then
              psi1=twiss(i,idp,3)+2d0*twiss(i,idp,6)
              psi2=twiss(i,idp,3)-2d0*twiss(i,idp,6)
              psi11=twiss(i1,idp,3)+2d0*twiss(i1,idp,6)
              psi21=twiss(i1,idp,3)-2d0*twiss(i1,idp,6)
              tt(3)=tt(3)+0.5d0*(sqrt(twiss(i,idp,5))*twiss(i,idp,7)
     1                         +sqrt(twiss(i1,idp,5))*twiss(i1,idp,7))
              tt(4)=tt(4)+0.5d0*(twiss(i,idp,5)*twiss(i,idp,7)
     1                         +twiss(i1,idp,5)*twiss(i1,idp,7))
              cc(2)=cc(2)+0.5d0*twiss(i,idp,5)*twiss(i,idp,7)*
     1         dcmplx(cos(2d0*twiss(i,idp,6)),sin(2d0*twiss(i,idp,6)))
     1                   +0.5d0*twiss(i1,idp,5)*twiss(i1,idp,7)*
     1         dcmplx(cos(2d0*twiss(i1,idp,6)),sin(2d0*twiss(i1,idp,6)))
              tt(7)=tt(7)+0.5d0*(sqrt(twiss(i,idp,2)*twiss(i,idp,5))
     1                         +sqrt(twiss(i1,idp,2)*twiss(i1,idp,5)))
              cc(3)=cc(3)+0.5d0*sqrt(twiss(i,idp,2))*twiss(i,idp,5)
     1              *dcmplx(cos(psi1),sin(psi1))
     1                   +0.5d0*sqrt(twiss(i1,idp,2))*twiss(i1,idp,5)
     1              *dcmplx(cos(psi11),sin(psi11))
              cc(4)=cc(4)+0.5d0*sqrt(twiss(i,idp,2))*twiss(i,idp,5)
     1              *dcmplx(cos(psi2),sin(psi2))
     1                   +0.5d0*sqrt(twiss(i1,idp,2))*twiss(i1,idp,5)
     1              *dcmplx(cos(psi21),sin(psi21))
              cc(5)=cc(5)+0.5d0*sqrt(twiss(i,idp,2))*twiss(i,idp,5)
     1              *dcmplx(cos(twiss(i,idp,3)),sin(twiss(i,idp,3)))
     1                   +0.5d0*sqrt(twiss(i1,idp,2))*twiss(i1,idp,5)
     1              *dcmplx(cos(twiss(i1,idp,3)),sin(twiss(i1,idp,3)))
            endif
            ip=latt(i)
c           ip=idvalc(i)
            als=rlist(ip+1)
            vs=rlist(ip+2)
            vsa=rlist(idvalc(i)+2)
            vsb=rlist(ip+8)
            e2s=rlist(ip+4)
            ndiv=min(1,int(als/ale)-1)
            ndivt=ndivt+ndiv+1
            do 60 k=1,ntwissfun
              twisss(k)=twiss(i+1,idp,k)
   60       continue
            do 70 n=1,ndiv
              rn=dble(n)/(ndiv+1)
              rlist(ip+1)=als*rn
              if(id .eq. icbend)then
                rlist(idvalc(i)+2)=vsa*rn
                rlist(ip+2)=vs*rn
                rlist(ip+4)=0.d0
                rlist(ip+8)=vsb*rn
              else
                rlist(ip+2)=vs*rn
              endif
              call qtwiss(twiss,idp,i,i+1,over)
              if(id.eq.icquad.or.id.eq.icbend) then
                tt(1)=tt(1)+sqrt(twiss(i1,idp,5))*twiss(i1,idp,7)
                tt(2)=tt(2)-twiss(i1,idp,5)
                cc(1)=cc(1)-twiss(i1,idp,5)*dcmplx(
     1             cos(2d0*twiss(i1,idp,6)),sin(2d0*twiss(i1,idp,6)))
                tt(6)=tt(6)+sqrt(twiss(i1,idp,2)*twiss(i1,idp,5))
              elseif(id.eq.icsext) then
                psi1=twiss(i1,idp,3)+2d0*twiss(i1,idp,6)
                psi2=twiss(i1,idp,3)-2d0*twiss(i1,idp,6)
                tt(3)=tt(3)+sqrt(twiss(i1,idp,5))*twiss(i1,idp,7)
                tt(4)=tt(4)+twiss(i1,idp,5)*twiss(i1,idp,7)
                cc(2)=cc(2)+twiss(i1,idp,5)*twiss(i1,idp,7)*dcmplx(
     1            cos(2d0*twiss(i1,idp,6)),sin(2d0*twiss(i1,idp,6)))
                tt(7)=tt(7)+sqrt(twiss(i1,idp,2)*twiss(i1,idp,5))
                cc(3)=cc(3)+sqrt(twiss(i1,idp,2))*twiss(i1,idp,5)
     1                *dcmplx(cos(psi1),sin(psi1))
                cc(4)=cc(4)+sqrt(twiss(i1,idp,2))*twiss(i1,idp,5)
     1                *dcmplx(cos(psi2),sin(psi2))
                cc(5)=cc(5)+sqrt(twiss(i1,idp,2))*twiss(i1,idp,5)
     1                *dcmplx(cos(twiss(i1,idp,3)),sin(twiss(i1,idp,3)))
              endif
   70       continue
            rlist(ip+1)=als
            if(id .eq. icbend)then
              rlist(ip+2)=vs
              rlist(idvalc(i)+2)=vsa
              rlist(ip+4)=e2s
              rlist(ip+8)=vsb
            else
              rlist(ip+2)=vs
            endif
            do 80 k=1,ntwissfun
              twiss(i+1,idp,k)=twisss(k)
  80        continue
          endif
   90   continue
        if(j.eq.0 .and. nstr.ne.0) then
          do 91 k=1,5
            cts(k)=ct(k)
   91     continue
          ts2=t(2)
          ts4=t(4)
        elseif(j.eq.nstr) then
          t(5)=t(5)+(ts2+ts4+t(2)+t(4))**2
     1                 +abs(cts(1)+cts(2)+ct(1)+ct(2))**2
          t(8)=t(8)+abs(cts(3)+ct(3))**2+abs(cts(5)+ct(5))**2
          t(9)=t(9)+abs(cts(4)+ct(4))**2+abs(cts(5)+ct(5))**2
          if(full) then
            sum1=((ts2+ts4+t(2)+t(4))**2
     1            +abs(cts(1)+cts(2)+ct(1)+ct(2))**2)*a(5)
            sum2=(abs(cts(3)+ct(3))**2 +abs(cts(5)+ct(5))**2)*a(8)
     1          +(abs(cts(4)+ct(4))**2 +abs(cts(5)+ct(5))**2)*a(9)
          endif
        else
          t(5)=t(5)+(t(2)+t(4))**2+abs(ct(1)+ct(2))**2
          t(8)=t(8)+abs(ct(3))**2+abs(ct(5))**2
          t(9)=t(9)+abs(ct(4))**2+abs(ct(5))**2
          if(full) then
            sum1=((t(2)+t(4))**2+abs(ct(1)+ct(2))**2)*a(5)
            sum2= (abs(ct(3))**2+abs(ct(5))**2)*a(8)
     1           +(abs(ct(4))**2+abs(ct(5))**2)*a(9)
          endif
        endif
        if(full) then
          if(j.ne.0 .or. nstr.eq.0) then
            if(line.ne.' ') then
              if(ida.eq.icbend .or. ida.eq.icquad ) then
                if(ida.eq.icbend) then
                  vv=vsb
                else
                  vv=vs
                endif
                t(1)=t(1)+(vv*tt(1)/ndivt)**2
                t(2)=t(2)+vv*tt(2)/ndivt
                ct(1)=ct(1)+vv*cc(1)/ndivt
                t(6)=t(6)+(vv*tt(6)/ndivt)**2
                line(10:)=autofg((vv*tt(1)/ndivt)**2*a(1),'8.5')
                line(37:)=autofg((vv*tt(6)/ndivt)**2*a(6),'8.5')
                write(lfno,'(2a)')' ',line
                line=' '
                tt(1)=0d0
                tt(2)=0d0
                cc(1)=(0d0,0d0)
                tt(6)=0d0
                ndivt=0
              elseif(ida.eq.icsext) then
                t(3)=t(3)+(vs*tt(3)/ndivt)**2
                t(4)=t(4)+vs*tt(4)/ndivt
                ct(2)=ct(2)+vs*cc(2)/ndivt
                t(7)=t(7)+(vs*tt(7)/ndivt)**2
                ct(3)=ct(3)+vs*cc(3)/ndivt
                ct(4)=ct(4)+vs*cc(4)/ndivt
                ct(5)=ct(5)+vs*cc(5)/ndivt
                line(19:)=autofg((vs*tt(3)/ndivt)**2*a(3),'8.5')
                line(46:)=autofg((vs*tt(7)/ndivt)**2*a(7),'8.5')
                write(lfno,'(2a)')' ',line
                line=' '
                tt(3)=0d0
                tt(4)=0d0
                cc(2)=(0d0,0d0)
                tt(7)=0d0
                cc(3)=(0d0,0d0)
                cc(4)=(0d0,0d0)
                cc(5)=(0d0,0d0)
                ndivt=0
              endif
            endif
            call elname(if+1,name)
            line=name
            line(28:)=autofg(sum1,'8.5')
            line(55:)=autofg(sum2,'8.5')
            write(lfno,'(2a)')' ',line
            line=' '
          endif
        endif
        do 92 k=1,5
          ct(k)=(0d0,0d0)
   92   continue
        t(2)=0d0
        t(4)=0d0
  100 continue
c
c ===== Dispersion =====
c ----- Emittance produced by dispersion
      write(lfno,'(a)')' <<< Emittance produced by disersion >>>'
      write(lfno,'(a)')' Emiy = 2Je/Jy(dE/E)^2 Ey**2/By.'
      write(lfno,'(a)')' Ey**2/By :'
c ----- dispersion due to Quad rotation
      vo(1)=autofg(t(1),'8.5')
      write(lfno,'(a)')' --- due to Quad rotation ---'
      write(lfno,'(3a)')' Ey**2/By=',vo(1),
     1                        ' /(2 Sin(Pi nu_y)**2) dtheta(quad)**2'
c ----- dispersion due to Sext misalign
      vo(2)=autofg(t(3)/4d0,'8.5')
      write(lfno,'(a)')' --- due to Sext misalignment ---'
      write(lfno,'(3a)')' Ey**2/By=',vo(2),
     1                        ' /(2 Sin(Pi nu_y)**2) dely(sext)**2'
c ----- dispersion due to Steerings
      t4=pi2/max(1,nstr)
      vo(3)=autofg(t4,'8.5')
      write(lfno,'(a)')' --- due to Correction Dipoles ---'
      write(lfno,'(3a)')' Ey**2/By=',vo(3),' xi_y dy**2/By/Tan(Pi nu_y)'
c ----- dispersion due to Orbit
      vo(4)=autofg(t(5)/8d0,'8.5')
      write(lfno,'(a)')' --- due to Orbit ---'
      write(lfno,'(3a)')' Ey**2/By= (dy**2/By) +',
     1                         vo(4),' /(2 Sin(Pi nu_y)**2 dy**2/By.'
c ----- evaluation for the present tune
      vo(1)=autofg(psix/pi2,'8.5')
      vo(2)=autofg(psiy/pi2,'8.5')
      write(lfno,'(5a)')'   [ For current tunes : nu_x ',vo(1),' nu_y ',
     1             vo(2),' ]'
      write(lfno,'(a)')' Emiy = 2Je/Jy(dE/E)^2 * ('
      vo(1)=autofg(t(1)/2d0/sin(psiy2)**2,'8.5')
      vo(2)=autofg(t(3)/8d0/sin(psiy2)**2,'8.5')
      vo(3)=autofg(t4/tan(psiy2),'8.5')
      vo(4)=autofg(t(5)/16d0/sin(psiy2)**2,'8.5')
      write(lfno,'(5a)')'       ',vo(1),' dtheta(quad)**2 + ',
     1                       vo(2),' dely(sext)**2 + '
      write(lfno,'(5a)')'       (',vo(3),'xi_y + 1) dy**2/By + ',
     1                       vo(4),' dy**2/By. )'
c ===== Linear Coupling =====
c ----- Emittance produced by Linear coupling
      write(lfno,'(a)')' <<< Emittance produced by linear coupling>>>'
      vo(1)=autofg(t(6)/4d0,'8.5')
      vo(2)=autofg(t(7)/16d0,'8.5')
      write(lfno,'(a)')
     1   ' --- Emittance due to Quad rotation Sext misalignment ---'
      write(lfno,'(2a)')' Emiy/Emix = tauy/taux *',
     1   ' (1/Sin(Pi(nu_x+nu_y))**2 + 1/Sin(Pi(nu_x-nu_y))**2)'
      write(lfno,'(5a)')
     1   '     *(',vo(1),' dtheta(quad)**2 + ',vo(2), ' dely(sext)**2 )'
      write(lfno,'(a)')' --- Emittance due to orbit ---'
      write(lfno,'(a)')' Emiy/Emix = tauy/taux * dy**2/By'
      vo(3)=autofg(t(8)/32d0,'8.5')
      vo(4)=autofg(t(9)/32d0,'8.5')
      write(lfno,'(5a)')
     1   '     *(',vo(3),' /Sin(Pi(nu_x+nu_y))**2 + ',
     1             vo(4),' /Sin(Pi(nu_x-nu_y))**2 ).'
c ----- evaluation for the present tune
      vo(1)=autofg(psix/pi2,'8.5')
      vo(2)=autofg(psiy/pi2,'8.5')
      write(lfno,'(5a)')'   [ For current tunes : nu_x ',vo(1),
     $     ' nu_y ',vo(2),' ]'
      vo(1)=autofg(t(6)/4d0*
     1        (1d0/sin(psix2+psiy2)**2+1d0/sin(psix2-psiy2)**2),'8.5')
      vo(2)=autofg(t(7)/16d0*
     1        (1d0/sin(psix2+psiy2)**2+1d0/sin(psix2-psiy2)**2),'8.5')
      vo(3)=autofg(t(8)/32d0/sin(psix2+psiy2)**2+
     1             t(9)/32d0/sin(psix2-psiy2)**2,'8.5')
      write(lfno,'(5a)')' Emiy/Emix = tauy/taux *(',
     1           vo(1),' dtheta(quad)**2 + ',vo(2),' dely(sect)**2 +'
      write(lfno,'(3a)')'                        ',vo(3),' dy**2/By ).'
      return
      end
