      subroutine tfaprt(latt,mult,lfno,word,
     1                  twiss,gammab)
      use tfstk
      use ffs
      use tffitcode
      implicit none
      integer*4 m,n,lfno,i,i1,it,nl,ia,ifany,lene
      real*8 a,al,bpole,ak,bpole1,brho1,s,sigx,sigy,
     $     sigx0,sigx1,sigy0,sigy1,sigxp,sigyp,theta,v,
     $     getva,sigxy0,sigxy1
      integer*4 latt(2,nlat),mult(nlat)
      real*8 twiss(nlat,-ndim:ndim,ntwissfun),gammab(nlat)
      real*8 beam1(21),beam2(21)
      character*(*) word
      character*8 patt,name,aout
      character*8 autofg
      character*73 buff
      logical*4 match,begin,temat,exist
      external trim
      ia(m,n)=((m+n+abs(m-n))**2+2*(m+n)-6*abs(m-n))/8
c     begin initialize for preventing compiler warning
      v=0.d0
      ak=0.d0
c     end   initialize for preventing compiler warning
      begin=.true.
      nl=0
      bpole=1.d0
 3    bpole1=getva(exist)
      if(exist)then
        bpole=bpole1
      endif
      call getwdl(word)
1     patt=word
      if(patt .eq. ' ' .and. begin)then
        patt='*'
      endif
      match=.false.
2     if(patt .eq. ' ')then
        go to 11
      endif
      do 10 i=1,nlat-1
        it=idtype(latt(1,i))
        if(it .eq. 1)then
          go to 10
        endif
        if(temat(latt,i,mult,name,patt))then
          if(mod(nl,65) .eq. 0)then
            begin=.false.
            aout=autofg(bpole,'S6.3')
            call trim(aout)
            write(lfno,'(a)')
     1      'SigXent SigYent CoXYent Element  SigXext SigYext CoXYext'
     1            //' SigXpk  SigYpk  a@'//aout(1:lene(aout))//'T'
            write(lfno,'(a)')
     1      ' micron  micron                   micron  micron        '
     1            //'  micron  micron   mm'
          endif
          if(.not. match)then
            call getwdl(word)
            match=.true.
          endif
          if(it .gt. 1 .and. it .le. 6)then
            v=rlist(latt(2,i)+2)
            al=rlist(latt(2,i)+1)
            if(al .ne. 0.d0)then
              ak=abs(v/al)
            else
              ak=0.d0
            endif
          elseif(it .eq. icmult)then
            v=rlist(latt(2,i)+13)
            al=rlist(latt(2,i)+1)
            if(al .ne. 0.d0)then
              ak=abs(v/al)
            else
              ak=0.d0
            endif
          endif
          if(it .eq. 2)then
            theta=rlist(latt(2,i)+5)
          elseif(it .eq. 4 .or. it .eq. 6)then
            theta=rlist(latt(2,i)+4)
          else
            theta=0.d0
          endif
          brho1=brho*gammab(i)/gammab(1)
          i1=i+1
          call tfbeam(twiss,gammab,i ,theta,beam1)
          call tfbeam(twiss,gammab,i1,theta,beam2)
          sigx0=sqrt(beam1(ia(1,1)))
          sigy0=sqrt(beam1(ia(3,3)))
          sigxy0=beam1(ia(1,3))/sigx0/sigy0
          sigx1=sqrt(beam2(ia(1,1)))
          sigy1=sqrt(beam2(ia(3,3)))
          sigxy1=beam2(ia(1,3))/sigx1/sigy1
          sigx=0.d0
          sigy=0.d0
          if(it .eq. 4 .or. it .eq. icmult)then
            if(v .gt. 0.d0)then
              if(beam1(ia(1,2)) .gt. 0.d0
     1           .and. beam2(ia(1,2)) .lt. 0.d0)then
                s=ak*beam1(ia(1,1))+beam1(ia(2,2))
                sigx=(s+sqrt(s**2
     1                 -4.d0*ak*(beam1(ia(1,1))*beam1(ia(2,2))
     1                          -beam1(ia(2,1))**2) ))/2.d0/ak
              endif
            elseif(v .lt. 0.d0)then
              if(beam1(ia(3,4)) .gt. 0.d0
     1           .and. beam2(ia(3,4)) .lt. 0.d0)then
                s=ak*beam1(ia(3,3))+beam1(ia(4,4))
                sigy=(s+sqrt(s**2
     1                 -4.d0*ak*(beam1(ia(3,3))*beam1(ia(4,4))
     1                          -beam1(ia(4,3))**2) ))/2.d0/ak
              endif
            endif
            sigx=sqrt(sigx)
            sigy=sqrt(sigy)
            if(ak .ne. 0.d0)then
              a=bpole/ak/brho1
            else
              a=bpole
            endif
            aout=autofg(a*1.d3,'6.3')
          elseif(it .eq. 6)then
c           b=ak*max(sigx0,sigy0,sigx1,sigy1)**2*.5d0*brho1
            if(ak .ne. 0.d0)then
              a=sqrt(2.d0*bpole/ak/brho1)
            else
              a=1.d0
            endif
            aout=autofg(a*1.d3,'6.3')
          else
            aout=' '
          endif
          sigxp=max(sigx0,sigx1,sigx)
          sigyp=max(sigy0,sigy1,sigy)
          nl=nl+1
          buff( 1: 8)=autofg(sigx0*1.d6,'8.5')
          buff( 9:16)=autofg(sigy0*1.d6,'8.5')
          buff(17:24)=autofg(sigxy0    ,'8.5')
          buff(25:33)=' '//name
          buff(34:41)=autofg(sigx1*1.d6,'8.5')
          buff(42:49)=autofg(sigy1*1.d6,'8.5')
          buff(50:57)=autofg(sigxy1    ,'8.5')
          buff(58:65)=autofg(sigxp*1.d6,'8.5')
          buff(66:73)=autofg(sigyp*1.d6,'8.5')
          write(lfno,9001)buff,aout
9001      format(2a)
        endif
10    continue
11    if(patt .eq. '*')then
        return
      endif
      if(match)then
        go to 1
      else
        if(ifany(word,'*%',1) .eq. 0)then
          if(begin)then
            match=.true.
            patt='*'
            go to 2
          endif
          return
        else
          go to 3
        endif
      endif
      end
