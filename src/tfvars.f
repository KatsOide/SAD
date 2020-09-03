      subroutine tfvars(nvar,kx,irtc,ret,lfno)
      use tfstk
      use ffs
      use ffs_pointer
      use tffitcode
      use tfcsi,only:ipoint
      use ffs_seg
      use eeval
      implicit none
      type (sad_descriptor) kx,kxr
      type (sad_dlist), pointer :: kli,fvl
      type (sad_rlist), pointer :: kla
      type (sad_comp), pointer :: cmps
      integer*4 lfno,i,iv,kk,level, itfuplevel,itfdownlevel,
     $     isp1,next,ifany
      integer*4 nvar,irtc,lenw
      real*8 x3,vmin,vmax,coup
      integer*4 k
      logical*4 ret,exist,tmatch
      character*15 autofg,v1,v2,v3,v4,v5,v6
      character*(MAXPNAME) key,tfkwrd
      character*(MAXPNAME+16) name,ncoup,ele
      integer*8  ifv,ifvh,ifnvar,ifvkey
      save ifv,ifvh,ifnvar,ifvkey
      data ifv /0/
      if(ifv .eq. 0)then
        ifv=ktadaloc(0,3)
        ifvh=ktfsymbolz('VariableRange',13)
        ifnvar=ktsalocb(0,'        ',MAXPNAME+16)
        ifvkey=ktsalocb(0,'        ',MAXPNAME)
        klist(ifv)=ktfsymbol+ktfcopy1(ifvh)
        klist(ifv+1)=ktfstring+ifnvar
        klist(ifv+2)=ktfstring+ifvkey
      endif
      write(lfno,'(a)')
     $     '!Variable     Keyword    Now          '//
     $  '!   Previous    Saved         Minimum      Maximum '//
     $     '      Couple      Coefficient'
      isp1=isp
      call peekwd(ele,next)
      if(ele .eq. ' ')then
        ele='*'
      endif
 1    exist=nvar .le. 0
      do i=1,nvar
        iv=nvevx(i)%ivarele
        if(nvevx(i)%ivcomp .eq. 0)then
          kk=nelvx(iv)%klp
        else
          kk=nvevx(i)%ivcomp
        endif
        k=idelc(kk)
        if(nvevx(i)%ivcomp .eq. 0)then
          name=pname(k)
        else
          call elnameK(nvevx(i)%ivcomp,name)
        endif
        if(tmatch(name,ele))then
          exist=.true.
          ipoint=next
          v1=autofg(nvevx(i)%valvar,'15.12')
          v2=autofg(nvevx(i)%valvar2,'12.9')
          call loc_comp(idval(k),cmps)
          x3=tfvcmp(cmps,nvevx(i)%ivvar)
c          x3=rlist(idval(k)+ivvar(i))
          v3=autofg(x3,'12.9')
          key=tfkwrd(idtype(k),nvevx(i)%ivvar)
          if(nvevx(i)%ivvar .eq. nelvx(iv)%ival)then
            vmin=nelvx(iv)%vlim(1)
            vmax=nelvx(iv)%vlim(2)
            call elname1(icomp(kk),ncoup,nvevx(i)%ivcomp .ne. 0)
            coup=couple(kk)
          else
            vmin=-1.d99
            vmax=1.d99
            ncoup='        '
            coup=1.d0
          endif
          ilist(1,ifnvar)=lenw(name)
          call tfpadstr(name,ifnvar+1,ilist(1,ifnvar))
          ilist(1,ifvkey)=lenw(key)
          call tfpadstr(key,ifvkey+1,ilist(1,ifvkey))
          rlist(ifv+3)=nvevx(i)%valvar
          call tclrfpe
          level=itfuplevel()
          call descr_sad(dlist(ifv),fvl)
          kxr=tfleval(fvl,.true.,irtc)
          if(irtc .ne. 0)then
            level=itfdownlevel()
            if(ierrorprint .ne. 0)then
              call tfaddmessage(' ',2,6)
            endif
            call termes(6,'Error in VariableRange ',
     $           pname(k)(:lpname(k))//' '//key(1:lenw(key)))
            kx=kxadalocnull(-1,nvar)
            return
          elseif(tfreallistq(kxr,kla))then
            if(kla%nl .eq. 2)then
              vmin=max(vmin,kla%rbody(1))
              vmax=min(vmax,kla%rbody(2))
            endif
          endif
          level=itfdownlevel()
          if(ret)then
            isp=isp+1
            dtastk(isp)=kxadaloc(-1,9,kli)
            kli%dbody(1)=kxsalocb(0,name,lenw(name))
            kli%dbody(2)=kxsalocb(0,key,lenw(key))
            kli%rbody(3)=nvevx(i)%valvar
            kli%rbody(4)=nvevx(i)%valvar2
            kli%rbody(5)=x3
            kli%rbody(6)=vmin
            kli%rbody(7)=vmax
            kli%dbody(8)%k=ktfstring+ktsalocb(0,ncoup,lenw(ncoup))
            kli%rbody(9)=coup
          endif
c     Note:
c     * POSITION of `name'  element: ivcomp(i) == 0 ? kk : ivcomp(i)
c     * POSITION of `ncoup' element: icomp(kk) if ivvar(i) .eq. ival(iv)
c     * icomp(kk) != 0 for 1 =< kk =< nlat
          if(nvevx(i)%ivvar .eq. nelvx(iv)%ival)then
            if((nvevx(i)%ivcomp .eq. 0 .and. icomp(kk) .eq. kk)
     $            .or. (icomp(kk) .eq. nvevx(i)%ivcomp))then
              ncoup='<--'
            endif
          endif
          v4=autofg(vmin,'12.8')
          v5=autofg(vmax,'12.8')
          v6=autofg(coup,'12.8')
          write(lfno,9001)name(1:12),key(1:8),v1,v2(1:12),v3(1:12),
     $         v4(1:12),v5(1:12),ncoup(1:12),v6(1:12)
 9001     format(1x,a,1x,a,a,' !',a,a,1x,a,2x,a,1x,a,a)
        endif
      enddo
      if(.not. exist)then
        if(ifany(ele,'*%{',1) .le. 0)then
          ele='*'
          go to 1
        endif
      endif
      if(ret)then
        kx=kxmakelist(isp1)
        isp=isp1
        irtc=0
      endif
      return
      end
