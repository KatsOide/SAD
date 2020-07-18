      subroutine tfdelete(def,del,unset)
      use tfstk
      use tfcode
      use efun
      implicit none
      type (sad_symdef) ,intent(inout):: def
      type (sad_symdef), pointer :: def1
      type (sad_defhash), pointer :: dhash
      type (sad_descriptor) kx
      integer*8 ka1,ka10,kp1,kadi,kadi0,kp0
      integer*4 i,kk,irtc,isp0
      logical*4 ,intent(in):: del,unset
      if(unset)then
        if(def%upval .ne. 0)then
          isp0=isp
          isp=isp+1
          ktastk(isp)=ktfoper+mtfunset
          isp=isp+1
          dtastk(isp)=sad_descr(def%sym)
          kx=tfefunref(isp0+1,.true.,irtc)
          if(irtc .ne. 0 .and. ierrorprint .ne. 0)then
            call tfreseterror
          endif
          isp=isp0
        endif
      endif
      ka1=def%upval
      def%upval=0
      do kk=1,2
        do while(ka1 .ne. 0)
          call loc_defhash(ka1,dhash)
          ka10=dhash%next
          if(dhash%gen .eq. maxgeneration)then
            do i=0,dhash%nhash
              kadi=dhash%dhash(i)%k
              do while(kadi .ne. 0)
                kadi0=klist(kadi)
                call tfcleardaloc(kadi)
                call tfree(kadi)
                kadi=kadi0
              enddo
            enddo
          else
            call tfcleardaloc(ka1)
          endif
          call tfree(ka1)
          ka1=ka10
        enddo
        ka1=def%downval
        def%downval=0
      enddo
      call tflocald(def%value)
      if(del)then
        kp0=def%prev
        kp1=def%next
        if(kp1 .ne. 0)then
          call loc1_symdef(kp1,def1)
          def1%prev=kp0
          if(max(0,def1%sym%gen) .eq. max(0,def%sym%gen) .and.
     $         def1%sym%override .eq. 0)then
            def1%sym%override=-2
          endif
        endif
        klist(kp0)=kp1
        def%sym%override=0
        def%sym%attr=def%len-6
        def%len=6
        call tfree(sad_loc(def%next))
        call tflocal1(sad_loc(def%sym%loc))
      else
        def%value=dtfcopy1(sad_descr(def%sym))
      endif
      return
      end

      integer*8 function ktfcopy(k)
      use tfstk, ktfc=>ktfcopy
      implicit none
      integer*8 ,intent(in):: k
      ktfcopy=ktfc(k)
      return
      end
