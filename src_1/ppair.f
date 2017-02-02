      subroutine ppair(latt,l1,l2,ntyp,ipair)
      use tfstk
      use ffs
      use tffitcode
      character*(MAXPNAME) name
      integer*8 latt(nlat)
      dimension ipair(2,*),ns(2),nf(2)
      goto 1
      entry ppair1(latt,l1,l2,ntyp,npair)
      is=0
      npair=0
      if(l1.le.l2) then
        ns(1)=l1
        ns(2)=nlat-1
        nf(1)=l2
        nf(2)=1
      else
        ns(1)=l1
        ns(2)=1
        nf(1)=nlat-1
        nf(2)=l2
      endif
      do 10 k=1,2
        do 12 l=ns(k),nf(k)
          if(idtype(ilist(2,latt(l))).eq.ntyp) then
            if(l.eq.1.or.idtype(ilist(2,latt(l-1))).ne.ntyp) then
              is=is+1
              if(mod(is,2).eq.0) then
                if(name.eq.pname(ilist(2,latt(l)))) then
                  npair=is/2
                else
                  is=is-1
                  name=pname(ilist(2,latt(l)))
                endif
              else
                name=pname(ilist(2,latt(l)))
              endif
            endif
          endif
   12   continue
   10 continue
      return
    1 continue
      is=0
      ll=0
      if(l1.le.l2) then
        ns(1)=l1
        ns(2)=nlat-1
        nf(1)=l2
        nf(2)=1
      else
        ns(1)=l1
        ns(2)=1
        nf(1)=nlat-1
        nf(2)=l2
      endif
      do 20 k=1,2
        do 20 l=ns(k),nf(k)
          if(idtype(ilist(2,latt(l))).eq.ntyp) then
            if(l.eq.1.or.idtype(ilist(2,latt(l-1))).ne.ntyp) then
              is=is+1
              if(mod(is,2).eq.0) then
                if(name.eq.pname(ilist(2,latt(l)))) then
                  ip=is/2
                  ipair(1,ip)=ll
                  ipair(2,ip)=l
                else
                  is=is-1
                  name=pname(ilist(2,latt(l)))
                  ll=l
                endif
              else
                name=pname(ilist(2,latt(l)))
                ll=l
              endif
            endif
          endif
   20 continue
      return
      end
