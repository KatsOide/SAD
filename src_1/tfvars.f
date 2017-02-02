      subroutine tfvars(nvar,kx,irtc,ret,lfno)
      use tfstk
      use ffs
      use ffs_pointer
      use tffitcode
      implicit none
      type (sad_descriptor) kx
      type (sad_list), pointer :: kli,kla
      integer*8 kxr,isp1
      integer*4 lfno,i,iv,kk,level, itfuplevel,itfdownlevel,
     $     next,ifany
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
        iv=ivarele(i)
        if(ivcomp(i) .eq. 0)then
          kk=klp(iv)
        else
          kk=ivcomp(i)
        endif
        k=ilist(2,latt(kk))
        if(ivcomp(i) .eq. 0)then
          name=pname(k)
        else
          call elnameK(ivcomp(i),name)
        endif
        if(tmatch(name,ele))then
          exist=.true.
          call cssetp(next)
          v1=autofg(valvar2(i,1),'15.12')
          v2=autofg(valvar2(i,2),'12.9')
          x3=rlist(idval(k)+ivvar(i))
          v3=autofg(x3,'12.9')
          key=tfkwrd(idtype(k),ivvar(i))
          if(ivvar(i) .eq. ival(iv))then
            vmin=vlim(iv,1)
            vmax=vlim(iv,2)
            call elname1(iele(kk),ncoup,ivcomp(i) .ne. 0)
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
          rlist(ifv+3)=valvar2(i,1)
          call tclrfpe
          level=itfuplevel()
          call tfleval(klist(ifv-3),kxr,.true.,irtc)
          if(irtc .ne. 0)then
            level=itfdownlevel()
            if(ierrorprint .ne. 0)then
              call tfaddmessage(' ',2,6)
            endif
            call termes(6,'Error in VariableRange ',
     $           pname(k)(1:lenw(pname(k)))//' '//key(1:lenw(key)))
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
            kli%rbody(3)=valvar2(i,1)
            kli%rbody(4)=valvar2(i,2)
            kli%rbody(5)=x3
            kli%rbody(6)=vmin
            kli%rbody(7)=vmax
            kli%body(8)=ktfstring+ktsalocb(0,ncoup,lenw(ncoup))
            kli%rbody(9)=coup
          endif
c     Note:
c     * POSITION of `name'  element: ivcomp(i) == 0 ? kk : ivcomp(i)
c     * POSITION of `ncoup' element: iele(kk) if ivvar(i) .eq. ival(iv)
c     * iele(kk) != 0 for 1 =< kk =< nlat
          if(ivvar(i) .eq. ival(iv))then
            if((ivcomp(i) .eq. 0 .and. iele(kk) .eq. kk)
     $            .or. (iele(kk) .eq. ivcomp(i)))then
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
