      subroutine temitf(plot,lfno)
      use ffs_pointer
      use tmacro
      use tffitcode
      use ffs_flag, only:trpt,wspac,intra
      use temw, only:tfinibeam,tfetwiss,iaez,beamplt
      use maccbk, only:i00
      implicit none
      integer*4 nparam,ntitle
      parameter (nparam=59,ntitle=26)
      integer*4 lfno,i,in,lfnos,in1,lu,j,indexs,lene
      real*8 scale(nparam),cod(6),dps,dpsa,ddl,rgetgl1,v
      real*8 sx(-2:2),sy(-2:2),stbl(4,nparam)
      real*8 trans(6,12),beam(21,2),ctrb(21,21),param(nparam,-2:2)
      character*24 title(nparam),unit
      character*78 buff
      character*11 autofg
      logical*4 plot,stab
      external trim
      real*8 r(6,13),rsave(6,13)
      data (title(i),i=1,15)/
     1           'COD x,mm,0.001          ',
     1           '    px/p0,mrad,0.001    ',
     1           '    y,mm,0.001          ',
     1           '    py/p0,mrad,0.001    ',
     1           '    z,mm,0.001          ',
     1           '    dp/p0, ,1.          ',
     1           'nuX, ,1.                ',
     1           'nuY, ,1.                ',
     1           'nuZ, ,1.                ',
     1           'U0,MeV,1E6              ',
     1           'Vc,MV ,1E6              ',
     1           'dz,mm ,0.001            ',
     1           'alpha, ,1.              ',
     1           'dl,mm,0.001             ',
     1           'dV/P0, ,1.              '/
      data (title(i),i=16,26)/
     1           'T0/tauX, ,1.            ',
     1           'T0/tauY, ,1.            ',
     1           'T0/tauZ, ,1.            ',
     1           'DX, ,1.                 ',
     1           'DY, ,1.                 ',
     1           'DZ, ,1.                 ',
     1           'EMIX,m,1.               ',
     1           'EMIY,m,1.               ',
     1           'EMIZ,m,1.               ',
     1           'SIGe, ,1.               ',
     1           'SIGz,mm,0.001           '/
      do 10 i=1,ntitle
        in=indexs(title(i),',',1)
        in=indexs(title(i),',',in+1)
        read(title(i)(in+1:),*)scale(i)
10    continue
      codin=twiss(1,0,mfitdx:mfitddp)
c      tw=tfetwiss(ri,codin,normali)
c      write(*,'(a,1p6g15.7)')'temitf-etwiss ',
c     $       tw(mfitax:mfitny)/[1d0,1d0,m_2pi,1d0,1d0,m_2pi]
      beamin=0.d0
      if(trpt)then
        beamin=tfinibeam(1)
      endif
      beamplt=wspac .or. intra
      cod=codin
      call temit(trans,cod,beam,ctrb,
     1     .not. trpt,iaez,plot,param(1,0),stab,lfno)
      dps=rgetgl1('PSPAN')
      dpsa=dps*.5d0
      if(dps .gt. 0.d0)then
        rsave=r
        lfnos=lfno
        lfno=0
        circ=c/omega0*pi2;
        ddl=dps*param(13,0)*.5d0*circ
        do 110 i=-2,2
          if(i .ne. 0)then
            dleng=param(14,0)+i*ddl
            call rsetgl1('FSHIFT',-dleng/circ)
            call tsetdvfs
            cod=codin
            call temit(trans,cod,beam,ctrb,
     1           .true.,iaez,.false.,param(1,i),stab,lfno)
            if(.not. stab .and. lfnos .ne. 0)then
              write(lfnos,9101)'Unstable at "dp/p0" =',i*dps*.5d0
9101          format(1x,a,f10.6)
            endif
          endif
          sx(i)=i*dps*.5d0
110     continue
        lfno=lfnos
        dleng=param(14,0)
        call rsetgl1('FSHIFT',-dleng/circ)
        call tsetdvfs
        trf0=param(12,0)
        codin=param(1:6,0)
        r=rsave
        if(lfno .ne. 0)then
          write(lfno,9001)-dps,-.5d0*dps,0.d0,dps*.5d0,dps
 9001     format(' "dp/p0"   ',5f12.6)
        endif
        do 120 i=1,ntitle
          in=indexs(title(i),',',1)
          buff(1:12)=title(i)(1:in-1)
          in1=indexs(title(i),',',in+1)
          unit=title(i)(in+1:in1-1)
          call trim(unit)
          lu=max(lene(unit),1)
          buff(73:78)=' '
          buff(79-lu:78)=unit
          do 130 j=0,4
            v=param(i,j-2)/scale(i)
            buff(j*12+13:(j+2)*12)=' '//autofg(v,'11.8')
            sy(j-2)=v
130       continue
          write(lfno,'(a)')buff
          stbl(1,i)=((sy(1)-sy(-1))/.75d0-(sy(2)-sy(-2))/6.d0)/dps
          stbl(2,i)=((sy(1)+sy(-1))/.75d0-(sy(2)+sy(-2))/12.d0
     1               -2.5d0*sy(0))/dpsa**2
          stbl(3,i)=(-sy(1)+sy(-1)+(sy(2)-sy(-2))*.5d0)/dpsa**3
          stbl(4,i)=(-4.d0*(sy(1)+sy(-1))+sy(2)+sy(-2)+6.d0*sy(0))
     1              /dpsa**4
120     continue
        if(lfno .ne. 0)then
          write(lfno,*)
          write(lfno,'(2a)')
     1'               f           f''          f''''         f''''''',
     1'        f'''''''''
          do 210 i=1,ntitle
            in=indexs(title(i),',',1)
            buff(1:12)=title(i)(1:in-1)
            buff(13:24)=' '//autofg(param(i,0)/scale(i),'11.8')
            buff(25:36)=' '//autofg(stbl(1,i),'11.8')
            buff(37:48)=' '//autofg(stbl(2,i),'11.8')
            buff(49:60)=' '//autofg(stbl(3,i),'11.8')
            buff(61:72)=' '//autofg(stbl(4,i),'11.8')
            in1=indexs(title(i),',',in+1)
            unit=title(i)(in+1:in1-1)
            call trim(unit)
            lu=max(lene(unit),1)
            buff(73:78)=' '
            buff(79-lu:78)=unit
            write(lfno,'(a)')buff(1:78)
 210      continue
        endif
      endif
      return
      end
