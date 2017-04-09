      subroutine pstati(word,latt,twiss,imon,emon,nmon,istr,nstr,
     1                   title,case,lfno)
      use tfstk
      use ffs
      use tffitcode
      use tfcsi, only:cssetp
      implicit real*8(a-h,o-z)
      parameter (ndl=16,nde=8,nd=20,lio=50)
      integer*8 ix,isa,isb,isc,is
      logical abbrev,exist,start,ghisto,lod,append
      character*(*) word,title,case
      character*256 tfconvstr
      dimension latt(*),twiss(*),imon(*),emon(*),istr(*)
     z,         xst(4,4),est(4,4)
      common /pstatis/na,nb,nc,isa,isb,isc,start
      exist=.true.
      if(word.eq.'SUM+') then
        if(start .and. max(nmon,nstr).ne.0) then
          if(na.eq.0) then
            isa=ktaloc(nd*ndl)
          elseif(mod(na,nd).eq.0) then
            call palocx(isa,ndl*na,ndl*(na+nd))
          endif
          na=na+1
          ix=ktaloc(max(nmon,nstr))
          call pstati1(latt,twiss,imon,emon,nmon,istr,nstr,xst,
     1         rlist(ix))
          call tfree(ix)
          call tmov(xst,rlist(isa+(na-1)*ndl),ndl)
        endif
      elseif(word.eq.'SUM-') then
        call getwdl(word)
        if(word.ne.'E') then
          exist=.false.
          if(na.eq.0) goto 999
          na=na-1
          if(na.eq.0) then
            call tfree(isa)
            goto 999
          elseif(mod(na,nd).eq.0) then
            call palocx(isa,ndl*(na+1),ndl*na)
          endif
        else
          if(nb.eq.0) goto 100
          nb=nb-1
          if(nb.eq.0) then
            call tfree(isb)
            goto 100
          elseif(mod(nb,nd).eq.0) then
            call palocx(isb,nde*(nb+1),nde*nb)
          endif
  100     if(nc.eq.0) goto 999
          nc=nc-1
          if(nc.eq.0) then
            call tfree(isc)
            goto 999
          elseif(mod(nc,nd).eq.0) then
            call palocx(isc,nde*(nc+1),nde*nc)
          endif
        endif
      elseif(word.eq.'SUM') then
        call getwdl(word)
        if(word.eq.'START') then
          start=.true.
        elseif(word.eq.'STOP') then
          start=.false.
        else
          if(word.eq.'OUT' .or. word.eq.'APPEND') then
            append=word.eq.'APPEND'
            itype=itfpeeko(ia,x,next)
            lfni1=int(x+.5d0)
            if(itype.eq.1) then
              call cssetp(next)
              write(word,'(''ftn'',i2.2)') lfni1
            elseif(itype.eq.101) then
              call cssetp(next)
              word=tfconvstr(101,ia,x,ncc,'*')
              if(word.eq.' ') then
                call permes('?Missing filename for SUM','.',' ',lfno)
                call getwdl(word)
                return
              endif
              call texpfn(word)
            endif
            inquire(unit=lio,opened=lod,iostat=ios)
            if(lod) close(lio)
            if(append) then
              open (lio,file=word,iostat=ios,access='APPEND',
     $             status='UNKNOWN',err=1991)
            else
              open (lio,file=word,iostat=ios,status='UNKNOWN',err=1991)
            endif
            go to 1992
 1991       call permes(' File open error in SUM',' ',' ',lfno)
            return
 1992       continue
            if(ios .ne. 0) then
              call permes(' Cannot open',' ',word,lfno)
              return
            endif
            exist=.true.
            ghisto=.true.
          else
            exist=.false.
            ghisto=.false.
          endif
          if(na.ne.0) then
            call pstati2(rlist(isa),xst,est,na,lfno)
          endif
          if(nb.ne.0) then
            write(lfno,1001) nb
 1001       format(1X,'Ndata=',I3,T19,'mean',T33,'sdev',T48,'min',T64,
     &           'max',t78,'cl=68%',t94,'cl=87%',t110,'cl=95%')
            call pstati3(rlist(isb),xst,est,nb,ghisto,title,case,lio,
     $           lfno)
          endif
          if(nc.ne.0) then
            write(lfno,1002) nc
 1002       format(1X,'Ndata=',I3,'  with intrabeam scattering')
            call pstati3(rlist(isc),xst,est,nc,ghisto,title,case,lio,
     $           lfno)
          endif
          if(ghisto .and. lfni1.ne.6) close(lio)
        endif
      elseif(abbrev(word,'SUMC_LR','_')) then
        if(na.ne.0) then
          call tfree(is)
          na=0
        endif
        if(nb.ne.0) then
          call tfree(isb)
          nb=0
        endif
        if(nc.ne.0) then
          call tfree(isc)
          nc=0
        endif
        start=.false.
      endif
999   if(exist) then
        call getwdl(word)
      endif
      return
      end
