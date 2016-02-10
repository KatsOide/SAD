c$Header: /SAD/cvsroot/oldsad/src/tfmemcheck.f,v 1.61 2010/08/21 01:54:58 oide Exp $
      subroutine tfmemcheck(isp1,itx,iax,vx,irtc)
      implicit none
      include 'inc/TFCBK.inc'
      include 'inc/TFCODE.inc'
      include 'inc/TFSTK.inc'
      include 'inc/TFMEM.inc'
      integer*4 isp1,itx,iax,irtc,i,j,k,it,ia,ig,nc,i1,k1,
     $     iad,l,iad1,k0,itavaloc,jj,
     $     nf,mf,itadaloc,itadvaloc,iax5,iaxi,isp0,m,
     $     itfmessage,iadi,iadi1,ii,ih
      integer*8 ip1,ip0,ip,ix,jx
      real*8 vx,xnmem,xnnet
      character*16 tfgetstr,name
      logical*4 tfcheckelement,freel,check
      if(isp1+1 .eq. isp .and. itastk(1,isp) .eq. ntfreal)then
        check=.true.
        freel=vstk(ivstkoffset+isp) .eq. 2.d0
      else
        check=.false.
        freel=.false.
      endif
      if(check)then
        do m=0,nsymhash
          j=itfcontext+m
          i=ilist(2,j)
          do while(i .ne. j)
            k=ilist(2,i-1)-4
            k0=i
            do while(k .gt. 0)
              it=ilist(1,k+3)
              ia=ilist(2,k+3)
              ig=ilist(2,k+8)
              if(ilist(2,k+2)-4 .ne. k0)then
                write(*,*)
     $      'Inconsistent previous definition generation pointer:',
     $               ilist(2,k+2)-4,
     $               ' <> previous name point ',k0
                go to 9000
              endif
              do l=1,2
                iad=ilist(l,k+1)
                iad1=k+1
                do while(iad .ne. 0)
                  if(iad .lt. 0)then
                    write(*,*)
     $                   'Negative definition pointer: direction =',l,
     $                   ' pointer =',iad,' previous pointer =',iad1
                    go to 9000
                  endif
                  if(ilist(1,iad+1) .eq. maxgeneration)then
                    if(ilist(3-l,iad) .ne. iad1)then
                      write(*,*)
     $               'Inconsistent definition pointer: direction =',
     $                     l,' in dispatch table, ',
     $                     ' pointer =',iad,' previous pointer =',iad1,
     $                     ' <> stored at the pointer ',ilist(3-l,iad)
                      go to 9000
                    endif
                    do jj=1,2
                      do ii=iad+2,iad+ilist(2,iad+1)+2
                        ih=ii-iad-2
                        iadi1=ii
                        iadi=ilist(jj,ii)
                        do while(iadi .ne. 0)
                          if(iadi .lt. 0)then
                            write(*,*)
     $                 'Negative definition pointer: direction =',l,
     $                  ' pointer =',iad,' previous pointer =',iadi,
     $                           ' hash =',ih
                            go to 9000
                          endif
                          if(ilist(3-l,iadi) .ne. iadi1)then
                            write(*,*)
     $                   'Inconsistent definition pointer: direction =',
     $                           l,
     $                   ' pointer =',iadi1,' previous pointer =',iadi1,
     $                   ' <> stored at the pointer ',ilist(3-l,iadi),
     $                   ' hash =',ih
                            go to 9000
                          endif
                          if(.not. tfcheckelement(
     $                         ntflist,ilist(2,iadi+1),.false.))then
                            write(*,*)
     $                   'Error in argument of definition: direction =',
     $                           l,
     $                  ' location of argument =',ilist(2,iadi+1),
     $                  ' pointer =',iadi, ' previous pointer =',iadi1,
     $                           ' hash =',ih
                            go to 9000
                          endif
                          if(.not. tfcheckelement(ilist(1,iadi+2),
     $                         ilist(2,iadi+2),.false.))then
                            write(*,*)
     $                'Error in value of definition: direction =',l,
     $                ' pointer =',iadi, ' previous pointer =',iadi1,
     $                    ' hash =',ih
                            go to 9000
                          endif
                          iadi1=iadi
                          iadi=ilist(l,iadi)
                        enddo
                      enddo
                    enddo
                  else
                    if(ilist(3-l,iad) .ne. iad1)then
                      write(*,*)
     $             'Inconsistent definition pointer: direction =',
     $                     l,
     $             ' pointer =',iad,' previous pointer =',iad1,
     $               ' <> stored at the pointer ',ilist(3-l,iad)
                      go to 9000
                    endif
                    if(.not. tfcheckelement(
     $                   ntflist,ilist(2,iad+1),.false.))then
                      write(*,*)
     $              'Error in argument of definition: direction =',
     $                     l,
     $              ' location of argument =',ilist(2,iad+1),
     $                ' pointer =',iad, ' previous pointer =',iad1
                      go to 9000
                    endif
                    if(.not. tfcheckelement(ilist(1,iad+2),
     $                   ilist(2,iad+2),.false.))then
                      write(*,*)
     $              'Error in value of definition: direction =',l,
     $              ' pointer =',iad, ' previous pointer =',iad1
                      go to 9000
                    endif
                    if(.not. tfcheckelement(ilist(1,iad+4),
     $                   ilist(2,iad+4),.false.))then
                      write(*,*)
     $          'Error in value of original definition: direction =',
     $                     l,
     $              ' pointer =',iad, ' previous pointer =',iad1
                      go to 9000
                    endif
                  endif
                  iad1=iad
                  iad=ilist(l,iad)
                enddo
              enddo
              if(ilist(1,k+2) .lt. 0)then
                write(*,*)'Negative attribute: ',ilist(1,k+2)
                go to 9000
              endif
              if(.not. tfcheckelement(it,ia,.false.))then
                go to 9000
              endif
              if(ilist(1,k-1) .lt. 10 .or.
     $             ilist(1,k-1) .gt. 13)then
                write(*,*)'Inconsistent size of definition table:',
     $               ilist(1,k-1)
                go to 9000
              endif
              k1=ilist(2,k-1)-4
              if(k1 .lt. -4)then
                write(*,*)
     $               'Negative pointer to next generation table:',k1
                go to 9000
              endif
              if(ilist(2,k+7) .ne. i+4)then
                write(*,*)'Inconsistent pointer to name table: ',
     $               ilist(2,k+7),' name table =',i
                go to 9000
              endif
              k0=k
              k=k1
            enddo
            i1=ilist(2,i)
            if(i1 .le. 0)then
              write(*,*)'Zero or negative pointer to next name:',i1
              go to 9000
            endif
            i=i1
          enddo
        enddo
        do l=1,levele
          j=itflocal+l
          i1=j
          i=ilist(2,j)
          do while(i .ne. j)
            if(ilist(1,i) .ne. i1)then
              write(*,*)
     $       'Inconsistent pointer to previous temporary element:',
     $             ilist(1,i)
              go to 9100
            endif
            if(.not. tfcheckelement(ilist(1,i+1),i+2,.true.))then
           call tfdebugprint(ilist(1,i+1),i+2,0.d0,'MemoryCheck[]',3)
              go to 9100
            endif
            if(ilist(2,i) .le. 0)then
              write(*,*)
     $             'Wrong pointer to next temporary element:',i
              go to 9100
            endif
            i1=i
            i=ilist(2,i)
          enddo
        enddo
      endif
      nf=0
      mf=0
      isp0=isp
      do ix=0,1
        do jx=0,mhash
          ip1=ich+jx*4+ix*2
          ip0=ip1
          ip=klist(ip1)
          do while(ip .ne. ip0)
            if(klist(ip+1) .ne. ip1)then
              write(*,*)
     $           'Free area (hash-list) is inconsistent: area =',ip,
     $           ' previous area =',ip1,' stored previous area =',
     $           klist(ip+1),' index =',ix
              go to 9200
            endif
            ip1=ip
            ip=klist(ip1)
          enddo
        enddo
      enddo
      do ix=0,nindex
        ip0=icp+ix*2
        ip=klist(ip0)
        do while(ip .ne. ip0)
          ip1=klist(ip)
          if(klist(ip1+1) .ne. ip)then
            write(*,*)
     $           'Free area (size-list) is inconsistent: area =',ip,
     $           ' next area =',ip1,' previous(next) area =',
     $           klist(ip1+1),' index =',ix
            go to 9200
          endif
          nf=nf+1
          if(ix .lt. nindex)then
            m=ix+1
            if(m .gt. 6)then
              if(klist(ip+6) .ne. m)then
                write(*,*)
     $           'Free area size is inconsistent: area =',ip,
     $           ' index =',m,' stored size =',klist(ip+6)
                go to 9200
              endif
            endif
          else
            m=klist(ip+6)
          endif
          mf=mf+m
          if(freel)then
            isp=isp+1
            itastk(1,isp)=ip-2
            itastk(2,isp)=m
          endif
          ip=ip1
        enddo
      enddo
      xnnet=nnet
      xnmem=nmem
      if(freel)then
        iax=itadvaloc(-1,5)
        iax5=itadaloc(0,nf)
        call tfsetlist(ntflist,iax5,0.d0,iax,5)
        do i=1,nf
          iaxi=itavaloc(0,2)
          rlist(iaxi+2)=itastk(1,isp0+i)
          rlist(iaxi+3)=itastk(2,isp0+i)
          call tfsetlist(ntflist,iaxi,0.d0,iax5,i)
        enddo
        isp=isp0
      else
        iax=itavaloc(-1,4)
      endif
      call tfsetlist(ntfreal,0,xnnet,iax,1)
      call tfsetlist(ntfreal,0,xnmem,iax,2)
      call tfsetlist(ntfreal,0,dble(nf),iax,3)
      call tfsetlist(ntfreal,0,mf+xnnet-xnmem,iax,4)
      itx=ntflist
      irtc=0
      return
 9000 name=tfgetstr(i+4,nc)
      write(*,*)'Name = ',name(1:nc),
     $     ', hash =',m,
     $     ', location =',k,
     $     ', parent location =',i,
     $     ', element type =',it,
     $     ', element pointer =',ia,
     $     ', generation =',ig,
     $     ', stack pointer =',isp,
     $     ', current generation =',lgeneration
      irtc=itfmessage(999,'General::memory',' ')
      itx=ntfoper
      iax=mtfnull
      return
 9100 write(*,*)'Level = ',l,
     $     ', root =',j,
     $     ', temporary element =',i,
     $     ', previous =',ilist(1,i),
     $     ', next =',ilist(2,i),
     $     ', type =',ilist(1,i+1),
     $     ', size =',ilist(1,i-1)
 9200 irtc=itfmessage(999,'General::memory',' ')
      itx=ntfoper
      iax=mtfnull
      return
      end

      logical*4 function tfcheckelement(it,ia,temp)
      implicit none
      include 'inc/TFCBK.inc'
      include 'inc/TFCODE.inc'
      integer*4 it,ia,i,m,iad,iav,nc,nw,lg,nrecmax
      parameter (nrecmax=1024)
      character*2 ord
      integer*4 nrec,kpat,ip
      logical*4 temp,tfcheckelementrc
      save nrec
      data nrec/0/
      tfcheckelement=.false.
      if(it .lt. 0)then
        write(*,*)'Negative type:',it
        return
      endif
      if(it .eq. ntfreal)then
        tfcheckelement=.true.
        return
      endif
      if(ia .lt. 0)then
        write(*,*)'Negative pointer:'
        return
      endif
      if(it .eq. ntfoper)then
        if(ia .lt. 0 .or. ia .gt. 2048)then
          write(*,*)'Illegal function number: ',ia
          return
        endif
      else
        if(ia .eq. 0)then
          write(*,*)'Zero pointer:'
          return
        endif
        if(it .eq. ntfsymbol)then
          if(ilist(2,ia-3) .ne. 0)then
            if(ilist(1,ia-1) .ne. ntfsymbol)then
              write(*,*)'Inconsistent element type for definition:',
     $             ilist(1,ia-1)
              return
            endif
          else
            if(ilist(1,ia-1) .eq. ntfstring)then
              if(ilist(2,ia-1) .ne. ia)then
      write(*,*)'Inconsistent element type (string as symbol):',
     $               ilist(1,ia-1),ia,ilist(2,ia-1)
                return
              endif
            elseif(ilist(1,ia-1) .ne. ntfsymbol)then
              write(*,*)'Inconsistent element type (symbol):',
     $             ilist(1,ia-1)
              return
            endif
          endif
        else
          if(ilist(1,ia-1) .ne. it)then
            write(*,*)'Inconsistent element type:',ilist(1,ia-1),it
            return
          endif
        endif
        if(ilist(1,ia-2) .gt. -1 .and. .not. temp)then
          write(*,*)'Wrong degeneration count:',ilist(1,ia-2)
          return
        endif
        if(it .eq. ntflist)then
          m=ilist(2,ia-1)
          if(m .lt. 0)then
            write(*,*)'Negative length of list:',m
            return
          endif
          nrec=nrec+1
          if(nrec .lt. nrecmax)then
            if(.not. tfcheckelementrc(
     $           ilist(1,ia),ilist(2,ia),.false.))then
              write(*,*)'Error in head of list: type =',ilist(1,ia),
     $             ', pointer =',ilist(2,ia)
              return
            endif
          endif
          iad=ilist(1,ia+1)
          if(iad .lt. 0)then
            write(*,*)'Negative pointer for list descriptor:',iad
            return
          endif
          iav=ilist(2,ia+1)
          if(iad .gt. 0)then
            if(nrec .lt. nrecmax)then
              do i=1,m
                if(.not. tfcheckelementrc(ilist(1,iad+i),
     $               ilist(2,iad+i),.false.))then
                  write(*,*)'Error in',i,ord(i),' element: type =',
     $                 ilist(1,iad+i),', pointer =',ilist(2,iad+i),
     $                 'length =',m,' List pointer =',ia
                  call tfdebugprint(ilist(1,ia),ilist(2,ia),0.d0,
     $                 ' head:',3)
                  if(i .ne. 1)then
                    call tfdebugprint(ilist(1,iad+1),ilist(2,iad+1),
     $                   rlist(iav+1),' 1st element: ',3)
                  endif
                  return
                endif
              enddo
            endif
          endif
          if(iav .lt. 0)then
            write(*,*)'Negative pointer for real elements:',iav
            return
          endif
          nrec=nrec-1
        elseif(it .eq. ntfstring)then
          nc=ilist(1,ia)
          nw=ilist(1,ia-3)-4
          if(nw .lt. nc/8+1)then
            write(*,*)'Inconsistent length of string:',
     $           nc,' size of allocated memory:',nw
            return
          endif
        elseif(it .eq. ntfsymbol)then
          nw=ilist(1,ia-3)
          if((nw .lt. 4 .and. nw .ne. 0) .or. nw .gt. 9)then
            if(it .ne. ntfsymbol .or. ilist(2,ia-1) .ne. ia)then
              write(*,*)'Inconsistent length of symbol: ',nw
              return
            endif
          endif
          lg=ilist(2,ia)
          if(lg .lt. -3)then
            write(*,*)'Illegal generation of symbol:',lg
            return
          endif
          ip=ilist(2,ia-1)
          if(ip .le. 0)then
            write(*,*)'Illegal name table of symbol:',ip
            return
          endif
        elseif(it .eq. ntfpat)then
          nw=ilist(1,ia-3)
          if(nw .lt. 13 .or. nw .gt. 16)then
            write(*,*)
     $           'Inconsistent size of pattern: length =',nw
            return
          endif
          kpat=ilist(1,ia)
          if(kpat .lt. -3)then
            write(*,*)'Illegal type of pattern: ',kpat
            return
          elseif(kpat .ge. 0)then
            if(.not. tfcheckelementrc(kpat,ilist(2,ia),.false.))then
              write(*,*)'Error in pattern expression: type =',kpat,
     $             ' pointer =',ilist(2,ia)
              return
            endif
          elseif(ilist(2,ia) .ne. 0)then
            if(.not. tfcheckelementrc(
     $           ntfsymbol,ilist(2,ia),.false.))then
              write(*,*)'Error in pattern head: type =',kpat,
     $             ' pointer =',ilist(2,ia)
              return
            endif
          endif
        else
          write(*,*)'Undefined type'
          return
        endif
      endif
      tfcheckelement=.true.
      return
      end

      logical*4 function tfcheckelementrc(it,ia,temp)
      implicit none
      integer*4 it,ia
      logical*4 temp,tfcheckelement
      tfcheckelementrc=tfcheckelement(it,ia,temp)
      return
      end

      subroutine tfreecheck(tag,it,ia,v)
      implicit none
      integer*4 it,ia,irtc
      real*8 v
      character*(*) tag
      call tfreecheck1(tag,it,ia,v,irtc)
      return
      end

      subroutine tfreecheck1(tag,it,ia,v,irtc)
      implicit none
      include 'inc/TFCBK.inc'
      include 'inc/TFMEM.inc'
      integer*8 ix,jx,ip1,ip,nf,mf,ip0,m
      integer*4 it,ia,irtc
      real*8 v
      character*(*) tag
      irtc=0
      nf=0
      mf=0
      do ix=0,1
        do jx=0,mhash
          ip1=ich+jx*4+ix*2
          ip0=ip1
          ip=klist(ip1)
          do while(ip .ne. ip0)
            if(klist(ip+1) .ne. ip1)then
              write(*,*)
     $           'Free area (hash-list) is inconsistent: area =',ip,
     $           ' previous area =',ip1,' stored previous area =',
     $           klist(ip+1),' index =',ix
              go to 9200
            endif
            ip1=ip
            ip=klist(ip1)
          enddo
        enddo
      enddo
      do ix=0,nindex
        ip0=icp+ix*2
        ip=klist(ip0)
        do while(ip .ne. ip0)
          ip1=klist(ip)
          if(klist(ip1+1) .ne. ip)then
            write(*,*)
     $           'Free area (size-list) is inconsistent: area =',ip,
     $           ' next area =',ip1,' previous(next) area =',
     $           klist(ip1+1),' index =',ix
            go to 9200
          endif
          nf=nf+1
          if(ix .lt. nindex)then
            m=ix+1
            if(m .gt. 6)then
              if(klist(ip+6) .ne. m)then
                write(*,*)
     $           'Free area size is inconsistent: area =',ip,
     $           ' index =',m,' stored size =',klist(ip+6)
                go to 9200
              endif
            endif
          else
            m=klist(ip+6)
          endif
          mf=mf+m
          ip=ip1
        enddo
      enddo
      return
 9200 if(it .ge. 0)then
        call tfdebugprint(it,ia,v,tag,3)
      endif
      irtc=-1
      return
      end

      integer*4 function memuse()
      implicit none
      include 'inc/TFMEM.inc'
      memuse=nnet
      return
      end

      subroutine tfgarbagecollect(isp1,itx,iax,vx,irtc)
      implicit none
      include 'inc/TFCBK.inc'
      include 'inc/TFCODE.inc'
      integer*4 isp1,itx,iax,irtc
      real*8 vx
      return
      end

      subroutine tfmemcheckprint(tag)
      implicit none
      include 'inc/TFCBK.inc'
      include 'inc/TFCODE.inc'
      include 'inc/TFSTK.inc'
      integer*4 isp0,it1,ia1,irtc
      real*8 v1
      character*(*) tag
      isp0=isp
      call tfmemcheck(isp,it1,ia1,v1,irtc)
      if(irtc .eq. 0)then
        write(*,*)tag
      else
        call tflocal1(irtc)
      endif
      isp=isp0
      return
      end
