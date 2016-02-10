      subroutine mstack2(iop,ax,latt,twiss,istr,nstr,imon,emon,
     1                   nmon,id)
      use tfstk
      use ffs
      use tffitcode
      implicit real*8(a-h,o-z)
      parameter (kstack=3)
      logical addxerr,addyerr
      dimension latt(2,nlat),twiss(nlat,-ndim:ndim,ntwissfun)
      dimension id(*)
      dimension istr(nstra,4),imon(nmona,4),emon(nmona,4)
      real*8 tgauss
      external pack
      common /mcstack/icoma,ipnt(kstack),iistck(kstack),
     1                mstkmx(kstack)
      include 'inc/common.inc'
c
      icom=icoma
      ip=ipnt(icom)
c
      if(icom.eq.1) then
        if(iop.eq.1) then
          ip=ip+1
          id(ip)=italoc(1+(3*nstr+1)/2)
          ilist(1,id(ip))=nstr
          do 10 i=1,nstr
            rlist(id(ip)+i)=rlist(latt(2,istr(istr(i,2),1))+11)
   10     continue
          call pmovi(istr(1,2),rlist(id(ip)+1+nstr),nstr)
        elseif(iop.eq.2) then
          no=ilist(1,id(ip))
          do 11 i=1,no
            j=ilist(mod(i-1,2)+1,id(ip)+no+1+(i-1)/2)
            if(istr(j,3).eq.0) then
              rlist(latt(2,istr(j,1))+11)=rlist(id(ip)+i)
            endif
   11     continue
          call tfree(int8(id(ip)))
          ip=ip-1
        elseif(iop.eq.3)then
          no=ilist(1,id(ip))
          ll=1
          do 13 i=1,nstr
            j=istr(i,2)
            ls=ll
            do 12 l=ls,no
              j1=ilist(mod(l-1,2)+1,id(ip)+no+1+(l-1)/2)
              if(j.eq.j1)then
                rlist(latt(2,istr(j,1))+11)=rlist(id(ip)+l)
     1                                     +rlist(latt(2,istr(j,1))+11)
                ll=l
                goto 13
              endif
   12       continue
            istr(j,3)=1
   13     continue
          nstrold=nstr
          call pack(istr(1,2),istr(1,3),nstr,nstrold)
          if(nstr.ne.nstrold) then
            nstract=nstr
            write(*,'(2(A,I4))')
     1           '  Correction dipole available :',nstr,' in ',nstrold
          endif
          call tfree(int8(id(ip)))
          ip=ip-1
        elseif(iop.eq.4)then
          no=ilist(1,id(ip))
          ll=1
          do 15 i=1,nstr
            j=istr(i,2)
            ls=ll
            do 14 l=ls,no
              j1=ilist(mod(l-1,2)+1,id(ip)+no+1+(l-1)/2)
              if(j.eq.j1)then
                rlist(latt(2,istr(j,1))+11)=rlist(id(ip)+l)
     1                                     -rlist(latt(2,istr(j,1))+11)
                ll=l
                goto 15
              endif
   14       continue
            istr(j,3)=1
   15     continue
          nstrold=nstr
          call pack(istr(1,2),istr(1,3),nstr,nstrold)
          if(nstr.ne.nstrold) then
            nstract=nstr
            write(*,'(2(A,I4))')
     1           '  Correction dipole available :',nstr,' in ',nstrold
          endif
          call tfree(int8(id(ip)))
          ip=ip-1
        elseif(iop.eq.5 .or. iop.eq.6 .or. iop.eq.7)then
          do 16 i=1,nstr
            j=istr(i,2)
            rlist(latt(2,istr(j,1))+11)=ax*rlist(latt(2,istr(j,1))+11)
   16     continue
        elseif(iop.eq.8)then
          itemp=italoc(nstr)
          do 17 i=1,nstr
            rlist(itemp-1+i)=rlist(latt(2,istr(istr(i,2),1))+11)
   17     continue
          no=ilist(1,id(ip))
          do 18 i=1,no
            j=ilist(mod(i-1,2)+1,id(ip)+1+(i-1)/2)
            if(istr(j,3).eq.0) then
              rlist(latt(2,istr(j,1))+11)=rlist(id(ip)+i)
            endif
   18     continue
          call tfree(int8(id(ip)))
          id(ip)=italoc(1+(3*nstr+1)/2)
          ilist(1,id(ip))=nstr
          call tmov(rlist(itemp),rlist(id(ip)+1),nstr)
          call pmovi(istr(1,2),rlist(id(ip)+1+nstr),nstr)
          call tfree(int8(itemp))
        endif
      elseif(icom.eq.2) then
        if(iop.eq.1) then
          ip=ip+1
          id(ip)=italoc(4*nlat)
          call tmov(twiss(1,0,15),rlist(id(ip)),nlat)
          call tmov(twiss(1,0,16),rlist(id(ip)+nlat),nlat)
          call tmov(twiss(1,0,17),rlist(id(ip)+2*nlat),nlat)
          call tmov(twiss(1,0,18),rlist(id(ip)+3*nlat),nlat)
        elseif(iop.eq.2)then
          call tmov(rlist(id(ip)),twiss(1,0,15),nlat)
          call tmov(rlist(id(ip)+nlat),twiss(1,0,16),nlat)
          call tmov(rlist(id(ip)+2*nlat),twiss(1,0,17),nlat)
          call tmov(rlist(id(ip)+3*nlat),twiss(1,0,18),nlat)
          call tfree(int8(id(ip)))
          ip=ip-1
        elseif(iop.eq.3)then
          call padd(rlist(id(ip)),twiss(1,0,15),nlat)
          call padd(rlist(id(ip)+nlat),twiss(1,0,16),nlat)
          call padd(rlist(id(ip)+2*nlat),twiss(1,0,17),nlat)
          call padd(rlist(id(ip)+3*nlat),twiss(1,0,18),nlat)
          call tfree(int8(id(ip)))
          ip=ip-1
        elseif(iop.eq.4)then
          call psub(rlist(id(ip)),twiss(1,0,15),nlat)
          call psub(rlist(id(ip)+nlat),twiss(1,0,16),nlat)
          call psub(rlist(id(ip)+2*nlat),twiss(1,0,17),nlat)
          call psub(rlist(id(ip)+3*nlat),twiss(1,0,18),nlat)
          call tfree(int8(id(ip)))
          ip=ip-1
        elseif(iop.eq.5 .or. iop.eq.6 .or. iop.eq.7)then
          call ptimes(ax,twiss(1,0,15),nlat)
          call ptimes(ax,twiss(1,0,16),nlat)
          call ptimes(ax,twiss(1,0,17),nlat)
          call ptimes(ax,twiss(1,0,18),nlat)
        elseif(iop.eq.8)then
          itemp=italoc(nlat)
          call tmov(twiss(1,0,15),rlist(itemp),nlat)
          call tmov(rlist(id(ip)),twiss(1,0,15),nlat)
          call tmov(twiss(1,0,16),rlist(itemp),nlat)
          call tmov(rlist(id(ip)+nlat),twiss(1,0,16),nlat)
          call tmov(twiss(1,0,17),rlist(itemp),nlat)
          call tmov(rlist(id(ip)+2*nlat),twiss(1,0,17),nlat)
          call tmov(twiss(1,0,18),rlist(itemp),nlat)
          call tmov(rlist(id(ip)+3*nlat),twiss(1,0,18),nlat)
          call tfree(int8(itemp))
        endif
      elseif(icom.eq.3) then
        if(iop.eq.1)then
          ip=ip+1
          id(ip)=italoc(1+(5*nmon+1)/2)
          ilist(1,id(ip))=nmon
          if(simulate) then
            do 20 i=1,nmon
              j=imon(imon(i,2),1)
              nq=imon(imon(i,2),4)
              rlist(id(ip)+i)=twiss(j,0,15)-twiss(j,ndim,15)
     1             -rlist(latt(2,nq)+5)+emon(imon(i,2),1)
              rlist(id(ip)+nmon+i)=twiss(j,0,17)-twiss(j,ndim,17)
     1             -rlist(latt(2,nq)+6)+emon(imon(i,2),2)
              ilist(mod(i-1,2)+1,id(ip)+2*nmon+1+(i-1)/2)=imon(i,2)
 20         continue
          else
            do i=1,nmon
              j=imon(imon(i,2),1)
              rlist(id(ip)+i)=twiss(j,0,15)
              rlist(id(ip)+nmon+i)=twiss(j,0,17)
              ilist(mod(i-1,2)+1,id(ip)+2*nmon+1+(i-1)/2)=imon(i,2)
            enddo
          endif
        elseif(iop.eq.2)then
          no=ilist(1,id(ip))
          do 21 i=1,no
            k=ilist(mod(i-1,2)+1,id(ip)+2*no+1+(i-1)/2)
            if(imon(k,3).eq.0)then
              j=imon(k,1)
              twiss(j,0,15)=rlist(id(ip)+i)
              twiss(j,0,17)=rlist(id(ip)+no+i)
            endif
   21     continue
          call tfree(int8(id(ip)))
          ip=ip-1
        elseif(iop.eq.3)then
          addxerr=simulate.and.errval(2).ne.0d0
          addyerr=simulate.and.errval(4).ne.0d0
          no=ilist(1,id(ip))
          ll=1
          do 23 i=1,nmon
            j=imon(i,2)
            ls=ll
            do 22 l=ls,no
              j1=ilist(mod(l-1,2)+1,id(ip)+2*no+1+(l-1)/2)
              if(j.eq.j1) then
                k=imon(j,1)
                nq=imon(j,4)
                twiss(k,0,15)=rlist(id(ip)+l)
     1                       +twiss(k,0,15)-twiss(k,ndim,15)
     1                         -rlist(latt(2,nq)+5)+emon(j,1)
c             ++  added 08Sep93 ++
                if(addxerr) twiss(k,0,15)=
     $               twiss(k,0,15)+errval(2)*tgauss()
c             ++
                twiss(k,0,17)=rlist(id(ip)+no+l)
     1                       +twiss(k,0,17)-twiss(k,ndim,17)
     1                         -rlist(latt(2,nq)+6)+emon(j,2)
c             ++  added 08Sep93 ++
                if(addyerr) twiss(k,0,17)=
     $               twiss(k,0,17)+errval(4)*tgauss()
c             ++
                ll=l
                goto 23
              endif
   22       continue
            imon(j,3)=1
   23     continue
          nmonold=nmon
          call pack(imon(1,2),imon(1,3),nmon,nmonold)
          if(nmon.ne.nmonold) then
            nmonact=nmon
            write(*,'(2(A,I4))')
     $           '  BPM available: ',nmon,' in ',nmonold
          endif
          call tfree(int8(id(ip)))
          ip=ip-1
        elseif(iop.eq.4)then
          addxerr=simulate.and.errval(2).ne.0d0
          addyerr=simulate.and.errval(4).ne.0d0
          no=ilist(1,id(ip))
          ll=1
          do 25 i=1,nmon
            j=imon(i,2)
            ls=ll
            do 24 l=ls,no
              j1=ilist(mod(l-1,2)+1,id(ip)+2*no+1+(l-1)/2)
              if(j.eq.j1) then
                k=imon(j,1)
                nq=imon(j,4)
                twiss(k,0,15)=rlist(id(ip)+l)
     1                       -(twiss(k,0,15)-twiss(k,ndim,15)
     1                         -rlist(latt(2,nq)+5)+emon(j,1))
c             ++  added 08Sep93 ++
                if(addxerr) twiss(k,0,15)=
     $               twiss(k,0,15)-errval(2)*tgauss()
c             ++
                twiss(k,0,17)=rlist(id(ip)+no+l)
     1                       -(twiss(k,0,17)-twiss(k,ndim,17)
     1                         -rlist(latt(2,nq)+6)+emon(j,2))
c             ++  added 08Sep93 ++
                if(addyerr) twiss(k,0,17)=
     $               twiss(k,0,17)-errval(4)*tgauss()
c             ++
                ll=l
                goto 25
              endif
   24       continue
            imon(j,3)=1
   25     continue
          nmonold=nmon
          call pack(imon(1,2),imon(1,3),nmon,nmonold)
          if(nmon.ne.nmonold) then
            nmonact=nmon
            write(*,'(2(A,I4))')
     $           '  BPM available: ',nmon,' in ',nmonold
          endif
          call tfree(int8(id(ip)))
          ip=ip-1
        elseif(iop.eq.5 .or. iop.eq.6 .or. iop.eq.7)then
          do 26 i=1,nmon
            j=imon(imon(i,2),1)
            twiss(j,0,15)=ax*twiss(j,0,15)
            twiss(j,0,17)=ax*twiss(j,0,17)
   26     continue
        elseif(iop.eq.8)then
          itemp=italoc(2*nmon)
          do 27 i=1,nmon
            j=imon(imon(i,2),1)
            nq=imon(imon(i,2),4)
            rlist(itemp-1+i)=twiss(j,0,15)-twiss(j,ndim,15)
     1                     -rlist(latt(2,nq)+5)+emon(imon(i,2),1)
            rlist(itemp-1+nmon+i)=twiss(j,0,17)-twiss(j,ndim,17)
     1                     -rlist(latt(2,nq)+6)+emon(imon(i,2),2)
   27     continue
          no=ilist(1,id(ip))
          do 28 i=1,no
            j=ilist(mod(i-1,2)+1,id(ip)+1+2*no+(i-1)/2)
            if(imon(j,3).eq.0) then
              twiss(imon(j,1),0,15)=rlist(id(ip)+1+i)
              twiss(imon(j,1),0,17)=rlist(id(ip)+1+no+i)
            endif
   28     continue
        endif
      endif
      ipnt(icom)=ip
      return
      end
