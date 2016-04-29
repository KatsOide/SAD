      subroutine tfinit(latt,ival,klp,vlim,
     1                  iele,iele1,iele2,ibzl,
     1                  errk,couple,mult)
      use tfstk
      use ffs
      use tffitcode
      implicit none
      integer*4 i,l,j,ikx,lele,ib,ibz,ibznext,ibzb,k,ibg,ibb
      integer*4 latt(2,nlat),klp(nele),ival(nele),ibzl(3,nlat)
      integer*4 iele(nlat),iele1(nlat),iele2(nlat),mult(nlat)
      real*8 errk(2,nlat),couple(nlat),vlim(nele,2),v
      do i=1,nele
        klp(i)=mult(i)
        ival(i)=0
      enddo
      do 10 l=1,nlat-1
        ikx=iele1(l)
        ilist(2,latt(2,l))=0
        couple(l)=1.d0
        errk(1,l)=1.d0
        errk(2,l)=0.d0
        if(ikx .gt. 0)then
          lele=idtype(latt(1,klp(ikx)))
          vlim(ikx,1)=-1.d10
          vlim(ikx,2)=1.d10
          go to (110,120,10,140,10,160,10,160,10,160,10,160),lele
          go to 210
110       ival(ikx)=1
          vlim(ikx,1)=.1d0
          go to 200
120       ival(ikx)=2
          go to 200
140       ival(ikx)=2
          go to 200
160       ival(ikx)=2
          go to 200
210       if(lele .eq. icSOL)then
            ival(ikx)=2
          elseif(lele .eq. icMULT)then
            ival(ikx)=kytbl(kwK1,icMULT)
          elseif(lele .eq. icCAVI)then
            ival(ikx)=2
          elseif(lele .eq. 32)then
            ival(ikx)=2
          elseif(lele .eq. 41)then
            ival(ikx)=0
          elseif(lele .eq. icMONI)then
            ival(ikx)=0
          elseif(lele .eq. 43)then
            ival(ikx)=0
          elseif(lele .eq. 35)then
            ival(ikx)=0
          elseif(lele .eq. 34)then
            ival(ikx)=0
          elseif(lele .eq. 36)then
            ival(ikx)=0
          elseif(lele .eq. 37)then
            ival(ikx)=0
          else
            ival(ikx)=0
            go to 10
          endif
 200      if(ival(ikx) .gt. 0)then
            v=rlist(idval(latt(1,klp(ikx)))+ival(ikx))
            if(v .ne. 0.d0)then
              errk(1,l)=rlist(latt(2,l)+ival(ikx))/v
            endif
          endif
        endif
10    continue
      mult(nlat)=0
      iele(nlat)=nlat
      iele1(nlat)=0
      call tfinimult(latt,1,mult,iele,iele1,klp)
      ib=1
      ibz=0
      ibzb=0
      ibznext=0
      do i=1,nlat
        ibzl(3,i)=0
      enddo
      do i=1,nlat-1
        if(idtype(latt(1,i)) .eq. icSOL)then
          if(ibz .ne. 0 .and.
     $         rlist(idval(latt(1,i))+kytbl(kwBND,icSOL))
     $         .ne. 0.d0)then
            ibznext=0
            if(ibzb .ne. 0)then
              if(rlist(idval(latt(1,i))+kytbl(kwGEO,icSOL))
     $             .ne. 0.d0)then
                ibg=i
                ibb=ibzb
              else
                ibg=ibzb
                ibb=i
              endif
              do k=ibzb,i
                ibzl(2,k)=ibg
                ibzl(3,k)=ibb
              enddo
            endif
            ibzb=0
          else
            if(ibz .eq. 0)then
              ibzb=i
            endif
            if(rlist(latt(2,i)+ilist(1,latt(2,i))) .gt. 0.d0)then
              ibz=i
              ibznext=i
            else
              ibznext=0
              do j=i+1,nlat-1
                if(idtype(latt(1,j)) .eq. icSOL)then
                  ibz=j
                  ibznext=j
                  exit
                endif
              enddo
            endif
            ibzl(2,i)=ibzb
          endif
        else
          ibzl(2,i)=ibzb
        endif
        ibzl(1,i)=ibz
        ibz=ibznext
      enddo
      if(ibz .ne. 0)then
        Write(*,*)'Missing end of solenoid: ',ibz
        call forcesf()
      endif
      ibzl(1,nlat)=0
      ibzl(2,nlat)=0
      return
      end

      real*8 function tfbzs(i,ibz)
      use tfstk
      use ffs
      use tffitcode
      implicit none
      integer*4 i,ibz,lp
      real*8 bzthre
      parameter (bzthre=1.d-20)
      ibz=ilist(i*3-2,ifibzl)
      if(ibz .gt. 0)then
        lp=ilist(2,ibz+ilattp)
        tfbzs=charge*(rlist(lp+kytbl(kwBZ,icSOL))
     $       +rlist(lp+kytbl(kwDBZ,icSOL)))
     $       *rlist(lp+ilist(1,lp))
     $       /(amass*rlist(ifgamm+i-1)/c)
        if(abs(tfbzs) .lt. bzthre)then
          tfbzs=0.d0
        endif
      else
        tfbzs=0.d0
      endif
c      write(*,*)'tfbzs ',i,ibz,brho,
c     $     rlist(ifgamm+i-1),rlist(ifgamm),tfbzs
      return
      end

      subroutine tfbndsol(i,ibg,ibb)
      use tfstk
      use ffs
      use tffitcode
      implicit none
      integer*4 i,ibg,ibb
      ibg=ilist(i*3-1,ifibzl)
      ibb=ilist(i*3  ,ifibzl)
      return
      end

      logical*4 function tfinsol(i)
      use tfstk
      use ffs
      use tffitcode
      implicit none
      integer*4 i,ibz,lp,ld
      ibz=ilist(i*3-2,ifibzl)
      tfinsol=.false.
      if(ibz .ne. 0)then
        if(ibz .lt. i)then
          tfinsol=.true.
        elseif(ibz .gt. i)then
          ld=ilist(1,i+ilattp)
          if(idtype(ld) .ne. icSOL)then
            tfinsol=.true.
          elseif(rlist(idval(ld)+kytbl(kwBND,icSOL)) .eq. 0.d0)then
            tfinsol=.true.
          endif            
        else
          lp=ilist(2,i+ilattp)
          if(rlist(lp+ilist(1,lp)) .lt. 0.d0)then
            tfinsol=.true.
          endif
        endif
      endif
      return
      end

      subroutine tffsrenumber(latt,mult,iele,iele1,iele2,klp,lfno)
      use tfstk
      use ffs
      use tffitcode
      implicit none
      integer*4 latt(2,nlat),mult(nlat),
     $     iele(nlat),iele1(nlat),iele2(nlat),klp(nele),
     $     next,ielm,lfno,i0
      character*(MAXPNAME+16) name
      logical*4 exist
      call peekwdp(name,next)
      i0=ielm(latt,name,1,mult,exist)
      if(exist)then
        call cssetp(next)
        call tfinimult(latt,i0,mult,iele,iele1,iele2,klp)
      else
        call termes(lfno,
     $       'Missing origin component for RENUM_BER',' ')
      endif
      return
      end

      subroutine tfinimult(latt,i0,mult,iele,iele1,klp)
      use tfstk
      use ffs
      use tffitcode
      implicit none
      integer*4 latt(2,nlat),i0,mult(nlat),im(nele),
     $     iele(nlat),iele1(nlat),klp(nele),i,
     $     ii,iie,k,ltyp,idx
      do i=1,nlat
        mult(i)=0
      enddo
      do i=1,nele
        klp(i)=0
        im(i)=0
c        ilist(i-1,im)=0
      enddo
      do i=1,nlat-1
        ii=mod(i+nlat+i0-3,nlat-1)+1
        iie=iele1(ii)
        k=klp(iie)
        if(k .eq. 0)then
          klp(iie)=ii
          iele(ii)=ii
        else
          if(mult(k) .eq. 0)then
            mult(k)=1
            ltyp=idtype(latt(1,ii))
            if(ltyp .gt. icNULL .and. ltyp .lt. icMXEL)then
              idx=kytbl(kwINDX,ltyp)
              if(idx .ne. 0)then
                mult(k)=max(1,int(rlist(idval(latt(1,ii))+idx)))
              endif
            endif
            im(iie)=mult(k)
c            ilist(iie-1,im)=mult(k)
          endif
          im(iie)=im(iie)+1
          mult(ii)=im(iie)
c          ilist(iie-1,im)=ilist(iie-1,im)+1
c          mult(ii)=ilist(iie-1,im)          
          iele(ii)=k
        endif
      enddo
      return
      end
