      subroutine twsdrw(latt,pos,icomp,word,wordp,lfnd,
     1                 twiss,idp,imon,emon,nmon,
     1                 title,case,exist)
      use tfstk
      use ffs
      use tffitcode
      use ffs_pointer, only:idelc,idtypec,pnamec
      implicit real*8 (a-h,o-z)
      type (ffs_bound) fbound
      parameter (nkey=35,nstyle=8)
      integer*8 latt(nlat),it,it1,jp
      real*8 pos(nlat)
      integer*4 icomp(nlat),imon(*)
      real*8 twiss(nlat,-ndim:ndim,ntwissfun),emon(*)
      real*8 ymax(nkey,2),ymin(nkey,2),gmin(2),gmax(2),fctr(2)
      integer*4 icat(nkey),ipw(2,nkey),kvar(nkey),
     1     line(nkey,2),lno(2),mp(nkey,2)
      real*8 sc(2,nkey),wc(nkey+9)
      character*8 keywrd(nkey),
     1     chart*80,charc*80,ct(nkey),cc(nkey),ctitl*16,ccase*16,
     1     outs*160,patt
      character*20 style(0:nstyle-1)
      character*(*) word,wordp,title,case
cslac character*40 dat
      character*30 dat
      logical abbrev,exist,fin,err,lat,tmatch,only,monly,monly1
      data keywrd /'AX      ','BX      ','NX      ','AY      ',
     1             'BY      ','NY      ','EX      ','EPX     ',
     1             'EY      ','EPY     ','R1      ','R2      ',
     1             'R3      ','R4      ','DX      ','DPX     ',
     1             'DY      ','DPY     ','PEX     ','PEPX    ',
     1             'PEY     ','PEPY    ','DETR    ','DMX     ',
     1             'DPMX    ','DMY     ','DPMY    ','SIGX    ',
     1             'SIGPX   ','SIGY    ','SIGPY   ','ROTATE  ',
     1             'DOX     ','DOY     ','MOMENTUM'/
      data icat   / 0         ,99        ,0         ,0
     1             ,99        ,0         ,1         ,0
     1             ,1         ,0         ,0         ,1
     1             ,-1        ,0         ,1         ,0
     1             ,1         ,0         ,1         ,0
     1             ,1         ,0         ,0         ,1
     1             ,0         ,1         ,0         ,1
     1             ,0         ,1         ,0         ,0
     1             ,1         ,1         ,9/
      data style  /'10                  ','10                  ',
     1             '.2 .06              ','.2 .06 .06 .06      ',
     1             '.2 .06 .12  .06     ','10                  ',
     1             '10                  ','10                  '/
      data (ct(i),i=1,4)  /'A0X1    ','2B06O0X1','N0X1    ','A0Y1    '/
      data (cc(i),i=1,4)  /'GXLX    ','MGUUDXLX','GXLX    ','GXLX    '/
      data (ct(i),i=5,8)  /'2B06O0Y1','N0Y1    ','H0X1    ','HA0X1   '/
      data (cc(i),i=5,8)  /'MGUUDXLX','GXLX    ','GXLX    ','GPXLX   '/
      data (ct(i),i=9,12) /'H0Y1    ','HA0Y1   ','R1      ','R2      '/
      data (cc(i),i=9,12) /'GXLX    ','GPXLX   ','        ','        '/
      data (ct(i),i=13,16)/'R3      ','R4      ','DX      ','DXA     '/
      data (cc(i),i=13,16)/'        ','        ','LL      ','LLP     '/
      data (ct(i),i=17,20)/'DY      ','DYA     ','H0PX1   ','HA0PX1  '/
      data (cc(i),i=17,20)/'LL      ','LLP     ','GXLLX   ','GPXLLX  '/
      data (ct(i),i=21,24)/'H0PY1   ','HA0PY1  ','DETR    ','DX      '/
      data (cc(i),i=21,24)/'GXLLX   ','GPXLLX  ','LLL     ','L       '/
      data (ct(i),i=25,28)/'DXP     ','DY      ','DYA     ','S0X1    '/
      data (cc(i),i=25,28)/'L P     ','L       ','L P     ','GXLX    '/
      data (ct(i),i=29,32)/'S0XA1   ','S0Y1    ','S0YA1   ','ROTATION'/
      data (cc(i),i=29,32)/'GXLPX   ','GXLX    ','GXLPX   ',' LLLLLLL'/
      data (ct(i),i=33,35)/'DX0M1   ','DY0M1   ','P       '/
      data (cc(i),i=33,35)/'LLXLX   ','LLXLX   ','L       '/
      data wc  /8.,   14.,   8.0,   8.,    14.,   8.0,
c               alfax sqrtbx psix   alfay  sqrtby psiy
     1          8.5,  10.0,  8.5,   10.0,  8.0,   9.5,
c               etax  etapx  etay   etapy  r1     r2
     1          9.5,  9.5,   8.5,   10.5,  8.5,   10.5,
c               r3    r4     dx     dxp    dy     dyp
     1          11.5, 13.5,  11.5,  13.5,  17.5,
c               etapx etapxp etapy  etapyp detr
     1          9.0,  11.0,  9.0,   11.0,
c               dX    dXp    dY     dYp
     1          8,    9.0,   8.,    9.,    33.5,
c               sigxp sigxp  sigy   sigyp  rotation
     1          15.2, 15.2,  5,
c               dox   doy    p
     1          3.0,  2.5,   18.5,  13.0,  20.5,  20.0,    32.0,
c               space  '   (sqrtm)  (m)    (mm)   (m-1)   (10^-3m^-1)
     1          32.,     16.4/
c               (GeV/c) (10^-n)
c     begin initialize for preventing compiler warning
      monly1=.false.
c     end   initialize for preventing compiler warning
      lat=.false.
      ls=1
      le=nlat
      only=.false.
      monly=.false.
      call getwdl2(word,wordp)
      ls=ielm(wordp,exist)
      if(exist) then
        call getwdl2(word,wordp)
        le=ielm(wordp,exist)
        if(exist) then
          call getwdl2(word,wordp)
          do 920 i=1,nlat-1
            if( tmatch(pname(ilist(2,latt(i))),wordp) ) then
              only=.true.
              patt=wordp
              call getwdl2(word,wordp)
              goto 921
            endif
920       continue
921       if(abbrev(word,'MON_ITOR','_')) then
            monly=.true.
            call getwdl2(word,wordp)
          endif
        else
          le=nlat
        endif
      else
        ls=1
        do 926 i=1,nlat-1
          if( tmatch(pname(ilist(2,latt(i))),wordp) ) then
            only=.true.
            patt=wordp
            call getwdl2(word,wordp)
            goto 927
          endif
926     continue
927     if(abbrev(word,'MON_ITOR','_')) then
          monly=.true.
          call getwdl2(word,wordp)
        endif
      endif
      call tffsbound(fbound)
      if(fbound%lb .ne. 1 .or. fbound%le .ne. nlat)then
        ls1=ls
        ls=max(min(ls1,le),fbound%lb)
        le=min(max(ls1,le),fbound%le+1)
      endif
      altotal=pos(nlat)-pos(1)
      alen=pos(le)-pos(ls)
      if(alen .le. 0d0)then
        alen=pos(nlat)+alen
      endif
      alen1=pos(ls)
      alen2=pos(ls)+alen
      ale=alen/100.d0
      np=mod(le-ls+1+nlat,nlat)
      if(np.eq.0) np=nlat
      istrt=ls
      if(ls.gt.le)then
        istop=nlat-1
      else
        istop=le-1
      endif
      fin=.false.
911   do 910 i=istrt,istop
        id=idtypec(i)
        if(id .le. icdodeca .or. id .eq. icmult .or.
     $       id .eq. iccavi .or. id .eq. ictcav)then
          ndiv=int(max(
     $         pos(i+1)-pos(i),rlist(latt(i)+kytbl(kwL,id)))/ale)
          np=np+ndiv+1
        endif
910   continue
      if(ls.gt.le.and..not.fin) then
        istrt=1
        istop=le-1
        fin=.true.
        goto 911
      endif
      call drwkwd(word,wordp,keywrd,icat,ipw,sc,nkey,lfnd,exist,err)
      if(err) return
      if(exist) call getwdl2(word,wordp)
      if(word.eq.'LAT') then
        lat=.true.
        call getwdl2(word,wordp)
      endif
      npage=0
      do 10 i=1,nkey
        npage=max(npage,ipw(1,i))
10    continue
      call fdate1(dat)
      if(npage.gt.1)then
        write(lfnd,'(a,i2,3a)')
     1             ' (   SAD/FFS  DRAW ',npage,' pages ',dat,' )'
      else
        write(lfnd,'(a,i2,3a)')
     1             ' (   SAD/FFS  DRAW ',npage,' page ',dat,' )'
      endif
      do 1000 ipage=1,npage
        write(lfnd,*)
     $       'NEWFRAME;SET CARD 80;SET FONT DUPLEX;SET TITLE SIZE -3'
        write(lfnd,*)'SET WINDOW X 3.3 12.1 Y 2.5  8.35'
        write(lfnd,*)'TITLE 8.5 9.1 CENTER '' '''
        write(lfnd,*)'MORE  ''',title(1:min(59,lene(title))),''''
        write(lfnd,*)'CASE '' ',case(1:min(59,lene(case))) ,''''
        write(lfnd,*)'SET TITLE SIZE -1.6'
cslac   write(lfnd,*)'TITLE 7 8.5 ''',dat,''''
        write(lfnd,*)'TITLE 8 8.5 ''',dat,  ''''
        nvar=0
        nwin=0
        do 20 i=1,nkey
          if(ipw(1,i).eq.ipage)then
            nvar=nvar+1
            nwin=max(nwin,ipw(2,i))
          endif
20      continue
        if(nvar.eq.0)goto 1000
        winl=2.5
        winr=11
        dwin=(8.3-2.5)/nwin
        do 1010 iwin=1,nwin
          winb=8.3-dwin*iwin
          wint=winb+dwin
          tsize=-3d0/sqrt(dble(nwin))
          write(lfnd,'(a,f6.2)')' SET TITLE SIZE',tsize
          km=0
          do 40 i=1,nkey
            if(ipw(1,i).eq.ipage)then
              if(ipw(2,i).eq.iwin) then
                km=km+1
                kvar(km)=i
              endif
            endif
40        continue
          it=ktaloc(np*2*km)
          lno(1)=0
          lno(2)=0
          left=icat(kvar(1))
          do 110 k=1,km
            if(icat(kvar(k)).eq.left) then
              lno(1)=lno(1)+1
              line(lno(1),1)=kvar(k)
            else
              lno(2)=lno(2)+1
              line(lno(2),2)=kvar(k)
            endif
110       continue
          it1=it
          do 120 l=1,2
            gmin(l)= 1d20
            gmax(l)=-1d20
            do 122 j=1,lno(l)
              if(line(j,l).eq.33 .or. line(j,l).eq.34) then
                monly1=monly
                monly=.true.
              endif
              call tdrwdt(line(j,l),rlist(it1),mp(j,l),
     1             ymin(j,l),ymax(j,l),ls,le,
     $             fbound%fb,fbound%fe,ale,np*km,
     $             only,monly,patt,
     1             latt,twiss,idp,pos,imon,emon,nmon)
              if(line(j,l).eq.33 .or. line(j,l).eq.34) then
                monly=monly1
              endif
              it1=it1+mp(j,l)*2
              if(icat(line(j,l)).eq.99) then
                ymin(j,l)=sqrt(ymin(j,l))
                ymax(j,l)=sqrt(ymax(j,l))
              elseif(line(j,l).eq.3.or.line(j,l).eq.6) then
                ymin(j,l)=ymin(j,l)/pi2
                ymax(j,l)=ymax(j,l)/pi2
              endif
              gmin(l)=min(gmin(l),sc(2,line(j,l)))
              gmax(l)=max(gmax(l),sc(1,line(j,l)))
122         continue
            if(gmin(l).ge.1d20) then
              call minmax(ymin(1,l),gmin(l),xx,lno(l))
            endif
            if(gmax(l).le.-1d20) then
              call minmax(ymax(1,l),xx,gmax(l),lno(l))
            endif
            if(gmin(l).eq.gmax(l)) then
              gmin(l)=gmin(l)-1d-4
              gmax(l)=gmax(l)+1d-4
            endif
            if(abs(gmin(l)) .gt. 0.01d0*abs(gmax(l)))then
              gmin(l)=1.05d0*gmin(l)-.05d0*gmax(l)
            else
              gmin(l)=0.d0
            endif
            if(abs(gmax(l)) .gt. 0.01d0*abs(gmin(l)))then
              gmax(l)=1.05d0*gmax(l)-.05d0*gmin(l)
            else
              gmax(l)=0.d0
            endif
120       continue
          do 130 l=1,2
            if(lno(l).eq.0)goto 130
            j=line(1,l)
            if(gmax(l)-gmin(l).lt.0.1d-3) then
              fctr(l)=1d6
              if(icat(j).eq.1)then
                ctitl=' (MM)'
                ccase='  GL '
                spaces=wc(nkey+1)+wc(nkey+5)
              elseif(icat(j).eq.-1)then
                ctitl=' (102-63M2-13)'
                ccase='    X  XLX  X'
                spaces=wc(nkey+1)+wc(nkey+7)
              elseif(icat(j).eq.99)then
                fctr(l)=1d0
                ctitl=' (2M062O1)'
                ccase='  MLUUUDU '
                spaces=wc(nkey+1)+wc(nkey+3)
              elseif(icat(j).eq.9)then
                fctr(l)=1d-3
                ctitl=' (KEV/C)'
                ccase='  LL  L '
                spaces=wc(nkey+1)+wc(nkey+8)
              else
                ctitl=' (102-63)'
                ccase='    X  X'
                spaces=wc(nkey+1)+wc(nkey+9)
              endif
            elseif(gmax(l)-gmin(l).lt.0.2) then
              fctr(l)=1d3
              if(icat(j).eq.1)then
                ctitl=' (MM)'
                ccase='  LL '
                spaces=wc(nkey+1)+wc(nkey+5)
              elseif(icat(j).eq.-1)then
                ctitl=' (102-33M2-13)'
                ccase='    X  XLX  X'
                spaces=wc(nkey+1)+wc(nkey+7)
              elseif(icat(j).eq.99)then
                fctr(l)=1d0
                ctitl=' (2M062O1)'
                ccase='  MLUUUDU '
                spaces=wc(nkey+1)+wc(nkey+3)
              elseif(icat(j).eq.9)then
                fctr(l)=1d-6
                ctitl=' (MEV/C)'
                ccase='   L  L '
                spaces=wc(nkey+1)+wc(nkey+8)
              else
                ctitl=' (102-33)'
                ccase='    X  X'
                spaces=wc(nkey+1)+wc(nkey+9)
              endif
            else
              fctr(l)=1d0
              if(icat(j).eq.1)then
                ctitl=' (M)'
                ccase='  L '
                spaces=wc(nkey+1)+wc(nkey+4)
              elseif(icat(j).eq.-1)then
                ctitl=' (M2-13)'
                ccase='  LX  X '
                spaces=wc(nkey+1)+wc(nkey+6)
              elseif(icat(j).eq.99)then
                fctr(l)=1d0
                ctitl=' (2M062O1)'
                ccase='  MLUUUDU '
                spaces=wc(nkey+1)+wc(nkey+3)
              elseif(icat(j).eq.9)then
                fctr(l)=1d-9
                ctitl=' (GEV/C)'
                ccase='   L  L '
                spaces=wc(nkey+1)+wc(nkey+8)
              else
                ctitl=' '
                ccase=' '
                spaces=0d0
              endif
            endif
            ip=1
            do 132 k=1,lno(l)
              if(k.eq.1) then
                chart(ip:ip)=''''
                charc(ip:ip)=''''
                ip=ip+1
              else
                chart(ip:ip)=','
                charc(ip:ip)=' '
                ip=ip+1
              spaces=spaces+wc(nkey+2)
              endif
              ln=lene(ct(line(k,l)))
              chart(ip:ip+ln-1)=ct(line(k,l))(1:ln)
              charc(ip:ip+ln-1)=cc(line(k,l))(1:ln)
              ip=ip+ln
              spaces=spaces+wc(line(k,l))
132         continue
            if(icat(j).ne.0 .or. icat(j).eq.0.and.fctr(l).ne.1d0)then
c              ???lene ignores last ')'
              ln=lene(ctitl)+1
              chart(ip:ip+ln-1)=ctitl(1:ln)
              ln=lene(ccase)+1
              charc(ip:ip+ln-1)=ccase(1:ln)
              ip=ip+ln
            endif
            chart(ip:ip)=''''
            charc(ip:ip)=''''
c           ip=ip+1
            if(l.eq.1)then
              xtitl=winl-1.4d0
              atitl=90d0
            else
              xtitl=winr+1.4d0
              atitl=-90d0
            endif
            if(nwin.eq.1) then
              xtitl=xtitl+sign(1d0,dble(l)-1.5d0)*0.4d0
            endif
            ytitl=winb+dwin/2d0
            write(outs,'(a,2f6.2,a,f4.0,a,f4.1,1x,a)')
     1      'TITLE',xtitl,ytitl,' CENTER ANGLE ',atitl,
     1              ' SPACE=',sngl(spaces/5d0),chart(1:ip)
            ln=lene(outs)
            if(ln.gt.71) then
              chart(1:70)=outs(1:70)
              chart(71:71)=''''
              write(lfnd,'(a)') chart(1:71)
              write(lfnd,*) 'MORE ''',outs(71:ln)
            else
              write(lfnd,*) outs(1:ln)
            endif
            write(lfnd,*)'CASE ',charc(1:ip)
130       continue
          write(lfnd,'(a,f6.3,1x,a)')
     1       ' SET OUTLINE ON;SET TICKS SIZE',sngl(0.1/dble(nwin)),'ON'
          write(lfnd,*)'SET WINDOW X ',sngl(winl),sngl(winr)
          write(lfnd,*)'SET WINDOW Y ',sngl(winb),sngl(wint)
          if(lno(2).ne.0) then
            write(lfnd,*)'SET TICKS RIGHT OFF'
          endif
c         if(iwin.ne.nwin) then
c           write(lfnd,*)'SET TICKS BOTTOM OFF'
c         endif
          write(lfnd,'(a,1p,2g15.7,a,2g15.7)')
     1      'SET LIMIT X ',sngl(alen1),sngl(alen2),
     $         ' Y ',sngl(gmin(1)*fctr(1)),sngl(gmax(1)*fctr(1))
          write(lfnd,'(a,f6.2,a)')
     1      ' SET LABELS SIZE=',sngl(-2.4d0/sqrt(dble(nwin))),
     1                          ' ALL OFF LEFT ON;PLOT AXES'
          write(lfnd,*)'SET TICKS ALL OFF;SET OUTLINE ALL OFF'
          if(iwin.eq.nwin) then
          write(lfnd,'(a,f6.2,a)')
     1      ' SET LABELS SIZE=',sngl(-2.0d0/sqrt(dble(nwin))),
     1                        ' ALL OFF BOTTOM ON;PLOT AXES'
          endif
          if(lno(2).ne.0) then
            write(lfnd,'(a,1p,2g15.7,a,2g15.7)')
     1      'SET LIMIT X ',sngl(alen1),sngl(alen2),
     $           ' Y ',sngl(gmin(2)*fctr(2)),sngl(gmax(2)*fctr(2))
            write(lfnd,*)
     $           'SET LABELS SIZE=',sngl(-2.4/sqrt(dble(nwin))),
     1           ' ALL OFF RIGHT ON'
            write(lfnd,*)
     1      'SET TICKS ALL OFF RIGHT ON;PLOT AXES'
          endif
          write(lfnd,*)
     1    'SET LABELS ALL OFF;SET TICKS ALL OFF;SET INTENSITY 1'
          jp=it
          do 150 l=1,2
            do 152 k=1,lno(l)
              write(lfnd,'(a,1p,2g15.7,a,2g15.7)')
     1             'SET LIMIT X ',sngl(alen1),sngl(alen2),
     $             ' Y ',sngl(gmin(l)),sngl(gmax(l))
              call tdinit(lfnd,'JOIN 1',style(
     1              mod( (l-1)*lno(max(1,l-1)) +k,nstyle) ) )
              if(icat(line(k,l)).eq.99) then
                do 160 j=1,mp(k,l)
                  call tdput(rlist(jp),sqrt(rlist(jp+1)))
                  jp=jp+2
160             continue
              elseif(line(k,l).eq.3.or.line(k,l).eq.6) then
                write(lfnd,*)sngl(alen1),' 0;',sngl(alen2),
     $               ' 0;JOIN 0 DOTS'
                do 162 j=1,mp(k,l)
                  call tdput(rlist(jp),rlist(jp+1)/pi2)
                  jp=jp+2
162             continue
              else
                write(lfnd,*)sngl(alen1),' 0;',
     $               sngl(alen2),' 0;JOIN 0 DOTS'
                do 170 j=1,mp(k,l)
                  call tdput(rlist(jp),rlist(jp+1))
                  jp=jp+2
170             continue
              endif
              call tdterm
152         continue
150       continue
          call tfree(int8(it))
1010    continue
1000  continue
      if(lat) then
        it=ktaloc(nlat)
        call tdlat(latt,ls,le,pos,rlist(it),
     1             icomp,' ',exist,lfnd)
        call tfree(int8(it))
      elseif(word.ne.' ') then
        istrt=ls
        if(ls.gt.le) then
          istop=nlat-1
        else
          istop=le-1
        endif
        fin=.false.
1111    do 1110 i=istrt,istop
          if( tmatch(pnamec(i),wordp) ) then
            goto 1100
          endif
1110    continue
        if(le.lt.ls.and..not.fin)then
          istrt=1
          istop=le-1
          fin=.true.
          goto 1111
        else
          exist=.false.
          return
        endif
1100    continue
        it=ktaloc(nlat)
        call tdlat(latt,ls,le,pos,rlist(it),
     1                            icomp,wordp,exist,lfnd)
        call tfree(int8(it))
      endif
      return
      end

      subroutine minmax(a,x,y,n)
      implicit real*8 (a-h,o-z)
      dimension a(n)
      x0=1.d35
      x1=-1.d35
      do 10 i=1,n
        x0=min(x0,a(i))
        x1=max(x1,a(i))
10    continue
      x=x0
      y=x1
      return
      end
