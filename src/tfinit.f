      subroutine tfinit
      use kyparam
      use tfstk
      use ffs
      use ffs_pointer
      use tffitcode
      implicit none
      type (sad_comp), pointer :: cmps
      integer*4 i,l,j,ikx,lele,ib,ibz,ibznext,ibzb,k,ibg,ibb
      integer*8 idv
      real*8 v
      evarini=.true.
      nelvx(1:nele)%klp=mult(1:nele)
      nelvx(1:nele)%ival=0
      do l=1,nlat-1
        ikx=iele1(l)
        couple(l)=1.d0
        errk(1,l)=1.d0
        errk(2,l)=0.d0
        if(ikx .gt. 0)then
          lele=idtypec(nelvx(ikx)%klp)
          nelvx(ikx)%vlim=(/-1.d10,1.d10/)
c          go to (110,120,10,140,10,160,10,160,10,160,10,160),lele
c          go to 210
          select case (lele)
          case (icDRFT)
            nelvx(ikx)%ival=1
            nelvx(ikx)%vlim(1)=.1d0
          case (icBEND,icQUAD,icSEXT,icOCTU,icDECA,icDODECA,
     $           icSOL,icCAVI,icTCAV)
            nelvx(ikx)%ival=2
          case (icMULT)
            nelvx(ikx)%ival=ky_K1_MULT
          case (icMARK)
            nelvx(ikx)%ival=0
            idv=idvalc(l)
            twiss(l,0,1:ntwissfun)=rlist(idv+1:idv+ntwissfun)
            cycle
          case default
            nelvx(ikx)%ival=0
            cycle
          end select
          if(nelvx(ikx)%ival .gt. 0)then
            call loc_comp(idvalc(nelvx(ikx)%klp),cmps)
            v=cmps%value(nelvx(ikx)%ival)
c            v=rlist(idvalc(klp(ikx))+ival(ikx))
            if(v .ne. 0.d0)then
              errk(1,l)=tfvalvar(l,nelvx(ikx)%ival)/v
c              errk(1,l)=rlist(latt(l)+ival(ikx))/v
            endif
          endif
        endif
      enddo
      mult(nlat)=0
      icomp(nlat)=nlat
      iele1(nlat)=0
      call tfinimult(1)
      ib=1
      ibz=0
      ibzb=0
      ibznext=0
c      do i=1,nlat
        ibzl(3,1:nlat)=0
c      enddo
      do i=1,nlat-1
        if(idtypec(i) .eq. icSOL)then
          if(ibz .ne. 0 .and.
     $         rlist(idvalc(i)+ky_BND_SOL)
     $         .ne. 0.d0)then
            ibznext=0
            if(ibzb .ne. 0)then
              if(rlist(idvalc(i)+ky_GEO_SOL)
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
            if(direlc(i) .gt. 0.d0)then
              ibz=i
              ibznext=i
            else
              ibznext=0
              do j=i+1,nlat-1
                if(idtypec(j) .eq. icSOL)then
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
        call abort
      endif
      ibzl(1,nlat)=0
      ibzl(2,nlat)=0
      return
      end

      real*8 function tfbzs(i,ibz)
      use kyparam
      use tfstk
      use ffs
      use sad_main, only:sad_comp
      use ffs_pointer, only:direlc,compelc
      implicit none
      integer*4 ,intent(in):: i
      integer*4 ,intent(out):: ibz
      type (sad_comp), pointer ::cmp
      real*8 ,parameter ::bzthre=1.d-20
      ibz=ilist(i*3-2,ifibzl)
      if(ibz .gt. 0)then
        call compelc(ibz,cmp)
        tfbzs=charge*(cmp%value(ky_BZ_SOL)
     $       +cmp%value(ky_DBZ_SOL))
     $       *direlc(ibz)
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
      use ffs_pointer
      implicit none
      integer*4 ,intent(in):: i
      integer*4 ,intent(out):: ibg,ibb
      ibg=ibzl(2,i)
      ibb=ibzl(3,i)
      return
      end

      logical*4 function tfinsol(i)
      use kyparam
      use tfstk
      use ffs
      use tffitcode
      use ffs_pointer, only:direlc,idelc
      implicit none
      integer*4 ,intent(in):: i
      integer*4 ibz,ld
      ibz=ilist(i*3-2,ifibzl)
      tfinsol=.false.
      if(ibz .ne. 0)then
        if(ibz .lt. i)then
          tfinsol=.true.
        else
          ld=idelc(i)
          if(ibz .gt. i)then
            if(idtype(ld) .ne. icSOL)then
              tfinsol=.true.
            elseif(rlist(idval(ld)+ky_BND_SOL) .eq. 0.d0)then
              tfinsol=.true.
            endif            
          elseif(rlist(idval(ld)+ky_BND_SOL) .eq. 0.d0)then
            tfinsol=.true.
          elseif(direlc(i) .lt. 0.d0)then
            tfinsol=.true.
          endif
        endif
      endif
      return
      end

      subroutine tffsrenumber(lfno)
      use tfstk
      use ffs
      use ffs_pointer
      use tffitcode
      use tfcsi,only:ipoint
      implicit none
      integer*4 ,intent(in):: lfno
      integer*4 next,ielm,i0
      character*(MAXPNAME+16) name
      logical*4 exist
      call peekwdp(name,next)
      i0=ielm(name,exist)
      if(exist)then
        ipoint=next
        call tfinimult(i0)
        evarini=.true.
      else
        call termes(lfno,
     $       'Missing origin component for RENUM_BER',' ')
      endif
      return
      end

      subroutine tfinimult(i0)
      use tfstk
      use ffs
      use ffs_pointer
      use tffitcode
      implicit none
      integer*4 ,intent(in):: i0
      integer*4 im(nele),i,ii,iie,k,ltyp,idx
      mult(1:nlat)=0
      nelvx(1:nele)%klp=0
      im=0
      do i=1,nlat-1
        ii=mod(i+nlat+i0-3,nlat-1)+1
        iie=iele1(ii)
        k=nelvx(iie)%klp
        if(k .eq. 0)then
          nelvx(iie)%klp=ii
          icomp(ii)=ii
        else
          if(mult(k) .eq. 0)then
            mult(k)=1
            ltyp=idtypec(ii)
            if(ltyp .gt. icNULL .and. ltyp .lt. icMXEL)then
              idx=kytbl(kwINDX,ltyp)
              if(idx .ne. 0)then
                mult(k)=max(1,
     $               int(rlist(idvalc(ii)+idx)))
              endif
            endif
            im(iie)=mult(k)
c            ilist(iie-1,im)=mult(k)
          endif
          im(iie)=im(iie)+1
          mult(ii)=im(iie)
c          ilist(iie-1,im)=ilist(iie-1,im)+1
c          mult(ii)=ilist(iie-1,im)          
          icomp(ii)=k
        endif
      enddo
      return
      end
