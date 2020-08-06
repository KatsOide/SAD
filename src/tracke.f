      subroutine tracke(latt,l20,sv,np,cmd,name,lfno)
      use tfstk
      use ffs
      use tffitcode
      implicit none
      integer*8 latt(nlat)
      integer*4 l20,lfno,np
      character*(*) cmd,name
      real*8 result(6,7),sv(5)
      save result
      call tracker(latt,l20,sv,result,np,cmd,name,lfno)
      return
      end

      subroutine tracker(latt,l20,sv,ss,np,cmd,name,lfno)
      use tfstk
      use ffs
      use tffitcode
      use track_tt ,ikptbl=>itt1,ix=>itt2,ix1=>itt3,
     $     l1=>itt4,l2s=>itt5
c
c     method   : function
c     STANDBY  : set up intitiol distribution
c     TRACK    : track particles from l1 to l2 with initial condition
c     CONT     : track particles from last l2 to new l2
c     PASS     : 
c     RESET    : reset distribution and free memory
c     SAVE     : record current position and distribution for TRACK/CONT
c     POS      : get position and angel (NY)
c     DISP     :                        (NY)
c
      implicit none
      integer*8 kx
      integer*8 latt(nlat)
      integer*4 lfno,l20,np,nl0,l2
      integer*4 irtc
      real*8 sv(5),sa(6),ss(6,7),es
      save sa
      character*(*) cmd,name
      logical*4 normal
c
c     This routine is called from MEA_SURE/DRAW/TRC_OD command @ tffsa.f
c     (Subroutine call flow at 2008/01/25)
c     Under tffsa() subroutine, SADScript interpreter is already initialized.
c
c     Store current random seed into NISTACK$FIXSEED stack
      call tfevalb('NISTACK$FIXSEED@Push[]',kx,irtc)
      l2=l20
      if(l20 .lt. nlat)then
        if(idtype(ilist(2,latt(l20))) .eq. icMARK)then
          l2=l20+1
        elseif(l20 .gt. 1 .and.
     $         idtype(ilist(2,latt(l20-1))) .eq. icMARK)then
          l2=l20-1
        endif
      endif
      if(cmd .eq. 'STANDBY')then
        novfl=0
        ix=ktaloc(np0*8)
        ix1=ktaloc(np0*8)
        ikptbl=ktaloc(np0*3)
        call tkptblini(ilist(1,ikptbl))
        nl0=nlat
        call ttinit(latt,
     1        rlist(ix      ),rlist(ix+np0  ),
     1        rlist(ix+np0*2),rlist(ix+np0*3),
     1        rlist(ix+np0*4),rlist(ix+np0*5),
     1        rlist(ix+np0*6))
        call tmov(rlist(ix),rlist(ix1),np0*8)
        l1=1
        l2s=1
      elseif(cmd .eq. 'TRACK' .or. cmd .eq. 'CONT'
     1       .or. cmd .eq. 'PASS')then
        np=np0
        if(cmd .ne. 'CONT')then
          call tmov(rlist(ix1),rlist(ix),np0*8)
        else
          l1=l2s
        endif
        nl0=nlat
        l2s=l2
        if(l2 .gt. l1)then
          call tturn0(np,l1,l2s,
     1          rlist(ix      ),rlist(ix+np0  ),
     1          rlist(ix+np0*2),rlist(ix+np0*3),
     1          rlist(ix+np0*4),rlist(ix+np0*5),
     1          rlist(ix+np0*6),rlist(ix+np0*7),ilist(1,ikptbl),1,
     $         normal)
        endif
        if(cmd .ne. 'PASS')then
          if(np .le. 1)then
            sv(1)=1.d20
            sv(2)=1.d20
            sv(3)=1.d20
            sv(4)=1.d20
            sv(5)=1.d20
          else
            call ttstat(np,
     1          rlist(ix      ),rlist(ix+np0  ),
     1          rlist(ix+np0*2),rlist(ix+np0*3),
     1          rlist(ix+np0*4),rlist(ix+np0*5),
     1          rlist(ix+np0*6),0.d0,
     1          name,sa,ss,es,.false.,.true.,lfno)
c$$$            call ttstat(np,
c$$$     1          rlist(ix      ),rlist(ix+np0  ),
c$$$     1          rlist(ix+np0*2),rlist(ix+np0*3),
c$$$     1          rlist(ix+np0*4),rlist(ix+np0*5),
c$$$     1          rlist(ix+np0*6),rlist(ilist(2,ifwakep+4)),
c$$$     1          name,
c$$$     1          sa,ss,es,.false.,.true.,lfno)
            sv(1)=sa(1)
            sv(2)=sa(3)
            sv(3)=ss(1,1)
            sv(4)=ss(3,3)
            sv(5)=ss(1,3)
            ss(:,7)=sa
c$$$            if(twake .or. lwake)then
c$$$              nb=ilist(1,ifwakep)
c$$$              sb=rlist(ifwakep+1)
c$$$              ns=ilist(1,ifwakep+2)
c$$$              write(lfno,9082)nb,sb,ns,pbunch
c$$$9082          format(' Bunches =',i4,
c$$$     1               ' Spacing =',f8.5,' m',
c$$$     1               ' Slices =',i3,
c$$$     1               ' Particles/bunch =',1pg12.4)
c$$$              if(bunchsta)then
c$$$                npb=np0/nb
c$$$                do 8010 i=1,nb
c$$$                  write(lfno,9081)i
c$$$9081              format(' Bunch #',i4,':')
c$$$                  ioff=ix+(i-1)*npb
c$$$                  call ttstat(npb,
c$$$     1                rlist(ioff      ),rlist(ioff+np0  ),
c$$$     1                rlist(ioff+np0*2),rlist(ioff+np0*3),
c$$$     1                rlist(ioff+np0*4),rlist(ioff+np0*5),
c$$$     1                rlist(ioff+np0*6),
c$$$     1                rlist(ilist(2,ifwakep+4)+(i-1)*npb),
c$$$     1                name,
c$$$     1                sa,ss,es,.false.,.false.,lfno)
c$$$8010            continue
c$$$              endif
c$$$            endif
          endif
        endif
c     Following lines are added by N. Yamamoto Apr. 25, '93
      elseif(cmd .eq. 'POS')then
          sv(1)=sa(1)
          sv(2)=sa(2)
          sv(3)=sa(3)
          sv(4)=sa(4)
      elseif(cmd .eq. 'DISP')then
          if(ss(6,6) .ne. 0.d0 ) then
            sv(1)=ss(1,6)/ss(6,6)
            sv(2)=ss(2,6)/ss(6,6)
            sv(3)=ss(3,6)/ss(6,6)
            sv(4)=ss(4,6)/ss(6,6)
          end if
          sv(5)=ss(6,6)
c     end of lines added by N. Yamamoto
      elseif(cmd .eq. 'SAVE')then
        call tmov(rlist(ix),rlist(ix1),np0*8)
        l1=l2s
      else
        nl0=nlat
        call tltrm(ilist(1,ikptbl))
        call tfree(ix)
        call tfree(ix1)
        call tfree(ikptbl)
      endif
c     Restore prior random seed from NISTACK$FIXSEED stack if needed
      if(fseed)then
        call tfevalb('NISTACK$FIXSEED@Pop[]',kx,irtc)
      else
        call tfevalb('NISTACK$FIXSEED@Discard[]',kx,irtc)
      endif
      call tclrfpe
      return
      end
