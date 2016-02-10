      subroutine  mstack(word,latt,twiss,istr,nstr,imon,emon,nmon)
      use ffs
      use tffitcode
      implicit real*8(a-h,o-z)
      parameter (kstack=3)
      logical abbrev,exist,number
      character*(*) word
      dimension latt(2,nlat),twiss(nlat,-ndim:ndim,ntwissfun)
      common /mcstack/icoma,ipnt(kstack),iistck(kstack),mstkmx(kstack)
      dimension istr(*),imon(*),emon(*)
c     begin initialize for preventing compiler warning
      number=.false.
c     end   initialize for preventing compiler warning
      if(word.eq.'PUSH') then
        iop=1
      elseif(word.eq.'DROP') then
        iop=2
      elseif(word.eq.'ADD') then
        iop=3
      elseif(word.eq.'SUB') then
        iop=4
      elseif(abbrev(word,'FAC_TOR','_'))then
        iop=5
      elseif(abbrev(word,'DIV_IDE','_')) then
        iop=6
      elseif(abbrev(word,'NEG_ATE','_')) then
        iop=7
        ax=-1d0
      elseif(word.eq.'SWAP') then
        iop=8
      else
        return
      endif
      if(iop.eq.5.or.iop.eq.6) then
        ax=getva(exist)
        if(exist .and. iop.eq.6) ax=1d0/ax
        number=exist
      endif
      call getwdl(word)
      exist=.true.
      if(word.eq.'C')then
        icom=1
      elseif(word.eq.'O') then
        icom=2
      elseif(word.eq.'M') then
        icom=3
      else
        if(icoma.eq.0) then
          if(iop.eq.5 .or. iop.eq.6) then
            print *,' No default object: select C or O'
          else
            print *,' No default stack object: select C, O, or M'
          endif
          return
        endif
        icom=icoma
        exist=.false.
      endif
      icoma=icom
      if(.not.number .and. (iop.eq.5.or.iop.eq.6)) then
        ax=getva(exist)
        if(.not.exist) then
          call getwdl(word)
          return
        endif
        if(iop.eq.6) ax=1d0/ax
      endif
      if(exist) call getwdl(word)
c93/10/9
c     if(icom.eq.3 .and. (iop.eq.5.or.iop.eq.6)) then
c       print *,'  Unexpected object (M): select C or O'
c       return
c     endif
      call mstack1(iop,ax,latt,twiss,istr,nstr,imon,emon,nmon)
!     print *,'stack level:',ipnt(1),ipnt(2),ipnt(3)
      return
      end
