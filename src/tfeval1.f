      module eval1
      use tfstk
      
      contains
      subroutine tfflagordef(isp1,kx,irtc)
      type (sad_descriptor) ,intent(out):: kx
      type (sad_symbol), pointer :: sym
      type (sad_string), pointer :: str
      integer*4 ,intent(out):: irtc
      integer*4 ,intent(in):: isp1
      integer*4 nc
      real*8 vx,fflogi
      logical*4 exist
      character*8 name
      if(isp /= isp1+1)then
        call tfdefinition(isp1,kx,irtc)
        return
      endif
      if(ktfsymbolq(dtastk(isp),sym))then
        call sym_symstr(sym,str)
        nc=min(8,str%nch)
        name(1:nc)=str%str(1:nc)
        call capita(name(1:nc))
        vx=fflogi(name(1:nc),exist)
c        write(*,*)'flagordef-vx ',exist,vx,'"'//name(1:nc)//'"'
        if(exist)then
          kx=dfromr(vx)
          irtc=0
        else
          call tfdefinition(isp1,kx,irtc)
        endif
      else
        call tfdefinition(isp1,kx,irtc)
      endif
      return
      end

      function tfeval1to(k1,k2,iopc,old,irtc) result(kx)
      use eeval
      implicit none
      type (sad_descriptor) ,intent(in):: k1,k2
      type (sad_descriptor) kx,kv,kr,ku,ks,tfeval1,tfset1
      type (sad_dlist), pointer :: kl1
      integer*4 ,intent(in):: iopc
      integer*4 ,intent(out):: irtc
      logical*4 ,intent(in):: old
      if(ktflistq(k1,kl1))then
        kv=tfleval(kl1,.true.,irtc)
        if(irtc /= 0)then
          return
        endif
      elseif(ktfsymbolq(k1))then
        kv=tfsyeval(k1,irtc)
        if(irtc /= 0)then
          return
        endif
      elseif(ktfpatq(k1))then
        kv=tfpateval(k1,irtc)
        if(irtc /= 0)then
          return
        endif
      else
        kv=k1
      endif
      if(iopc .eq. mtfaddto)then
        kr=tfeval1(kv,k2,mtfplus,irtc)
      elseif(iopc .eq. mtftimesby)then
        kr=tfeval1(kv,k2,mtfmult,irtc)
      elseif(iopc .eq. mtfsubtractfrom)then
        ku=tfeval1(sad_descr(-1.d0),k2,mtfmult,irtc)
        if(irtc /= 0)then
          return
        endif
        kr=tfeval1(kv,ku,mtfplus,irtc)
      else
        ku=tfeval1(k2,sad_descr(-1.d0),mtfpower,irtc)
        if(irtc /= 0)then
          return
        endif
        kr=tfeval1(kv,ku,mtfmult,irtc)
      endif
      if(irtc /= 0)then
        return
      endif
      ks=tfeevaldef(k1,irtc)
      if(irtc /= 0)then
        return
      endif
      kx=tfset1(ks,kr,mtfset,irtc)
      if(old)then
        kx=kv
      endif
      return
      end

      end module eval1

      logical*4 function tfgetstoredp(ks0,kp,def,irtc)
      use eval1
      implicit none
      type (sad_descriptor) ,intent(in):: ks0
      type (sad_descriptor) ks
      type (sad_dlist), pointer :: lists
      type (sad_symdef), pointer :: symd
      integer*8 kp
      integer*4 ,intent(out):: irtc
      integer*4 itfmessageexp,itfmessage
      logical*4 ,intent(out):: def
      logical*4 ev
      def=.false.
      irtc=0
      tfgetstoredp=.false.
      ks=ks0
      do while(ktflistq(ks,lists))
        def=.true.
        ks=lists%head
      enddo
      if(ktfsymbolqdef(ks%k,symd))then
        if(symd%sym%override .eq. 0)then
          irtc=itfmessageexp(999,'General::invset',ks0)
          return
        endif
        if(ktfprotectedqo(symd%sym))then
          irtc=itfmessage(999,'General::protect','""')
          return
        endif
        if(def)then
          call loc_sad(ktfaddrd(ks0),lists)
          call tfgetdefargp(lists,ktfaddr(ks),kp,ev,irtc)
          if(irtc /= 0)then
            return
          endif
        else
          kp=ksad_loc(symd%value%k)
        endif
      else
        return
      endif
      tfgetstoredp=.true.
      return
      end function 


      recursive function tfeval1(k1,k2,iopc1,irtc) result(kx)
      use eval1
      use complex,only:tfcmplx
      implicit none
      type (sad_descriptor) ,intent(in):: k1,k2
      type (sad_descriptor) kx,tfefunrefu,tfeexpr,tfearray
      integer*4 ,intent(in):: iopc1
      integer*4 ,intent(out):: irtc
      integer*4 isp1
      select case (iopc1)
      case (mtfrevpower)
        kx=tfeval1(k2,k1,mtfpower,irtc)
      case (mtfcomp)
        irtc=0
        kx=k2
      case (mtffun:mtfruledelayed,mtfpattest,mtfslot,mtfslotseq,
     $       mtfalt,mtfrepeated,mtfrepeatednull)
        kx=tfeexpr(k1,k2,iopc1)
        irtc=0
      case (mtfnull:mtfdiv,mtfpower,mtfgreater:mtfless,mtfand:mtfnot,
     $       mtfleftbra:mtfrightbrace,
     $       mtfcomplex:mtfcomma,mtfdot,mtfend)
        if(tfnumberq(k2) .and.
     $     (tfnumberq(k1) .or. iopc1 .eq. mtfnot))then
          kx=tfcmplx(k1,k2,iopc1,irtc)
          return
        endif
        if(tflistq(k2) .or. tflistq(k1))then
          kx=tfearray(k1,k2,iopc1,irtc)
        else
          kx=tfeexpr(k1,k2,iopc1)
          irtc=0
        endif
      case default
        isp=isp+1
        isp1=isp
        ktastk(isp)=ktfoper+iopc1
        isp=isp+1
        dtastk(isp)=k1
        isp=isp+1
        dtastk(isp)=k2
        kx=tfefunrefu(isp1,irtc)
        isp=isp1-1
      end select
      return
      end
