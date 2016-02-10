      subroutine tqr(a,w,eig,ibtab,vx,n,ndim)
      implicit none
      integer*4 itmax
      parameter (itmax=30)
      real*8 vmax,vmin,threj,alpha,decth
      parameter (vmax=1.d10,vmin=1.d0/vmax,
     $     decth=1.d-16,
     $     threj=1.d-8,alpha=1.d-6)
      integer*4 ibp,ib,ie,i,j,ie1,iter,is,
     $     is1,is2,i1,i2,i3,n,ndim,jm,k,j1,ii
      integer*4 ibtab(n)
      real*8 a(n,n),w(ndim,n),eig(2,n),vx(n)
      real*8 px(2),py(2)
      real*8 anorm,c,s,p,q,v1,v2,am,am1,w1,w2,w3,aa,ww,
     $     dec,sqrd,r,da,pk,pl,u,v,pm,pn,a1,a2,x,y,ee,sa,
     $     r1,r2,s1,aa1
      logical*4 paired,jordan,isnan
      complex*16 ca,cc,cu,cr,ck,cm
      equivalence (ck,px),(cm,py)
      equivalence (px,pk),(px(2),pl),(py,pm),(py(2),pn)
c     begin initialize for preventing compiler warning
      is2=0
      w1=0.d0
      w2=0.d0
      w3=0.d0
c     end   initialize for preventing compiler warning
      jordan=.false.
      ibp=1
      ibtab(1)=1
      ib=1
      ie=n
      do 10 i=1,ie-2
        a(i+2,i)=0.d0
        a(min(i+3,ie),i)=0.d0
10    continue
      anorm=0.d0
      do 1020 i=1,n
        do 1030 j=max(1,i-1),n
          anorm=anorm+abs(vx(i)/vx(j)*a(i,j))
1030    continue
1020  continue
      anorm=anorm*2/n**2
1     if(ie .le. ib+1)then
        if(ib .eq. 1)then
          go to 2000
        else
          ie=ib-1
          ibp=ibp-1
          ib=ibtab(ibp)
          go to 1
        endif
      endif
      ie1=ie-1
      iter=0
2     continue
c       write(*,'(1P4G15.7)')((vx(i)/vx(j)*a(i,j),j=1,4),i=1,4)
      do 310 i=ie-1,ib,-1
        i1=i+1
        s=max(anorm,abs(a(i1,i1))+abs(a(i,i)))
        if(abs(vx(i1)/vx(i)*a(i1,i))+s .eq. s)then
          a(i1,i)=0.d0
          if(i1 .eq. ie)then
            ie=ie-1
            go to 1
          endif
          ibp=ibp+1
          ibtab(ibp)=ib
          ib=i1
          go to 1
        endif
310   continue
      do 4001 i=ib,ie
        if(abs(vx(i)) .lt. vmin  .or. abs(vx(i)) .gt. vmax)then
          do 4010 j=max(i-1,ib),n
            a(i,j)=a(i,j)*vx(i)
4010      continue
          do 4020 j=1,min(ie,i+1)
            a(j,i)=a(j,i)/vx(i)
4020      continue
          do 4030 j=1,n
            w(j,i)=w(j,i)/vx(i)
4030      continue
          vx(i)=1.d0
        endif
4001  continue
      do 4510 is=ie-2,ib,-1
        is1=is+1
        is2=is1+1
        a1=a(ie,ie)-a(is,is)
        a2=a(ie1,ie1)-a(is,is)
        w1=((a1*a2-a(ie1,ie)*a(ie,ie1))/a(is1,is)+a(is,is1))
     1     *vx(is)/vx(is1)
        if(w1 .eq. 0.d0)then
          w1=-alpha*a(ie1,ie)*a(ie,ie1)/a(is1,is)*vx(is)/vx(is1)
        endif
        w2=a(is1,is1)-a(is,is)-a1-a2
        w3=a(is2,is1)*vx(is2)/vx(is1)
        if(is .eq. ib)then
          go to 4511
        endif
        u=abs(a(is,is-1)*vx(is)/vx(is-1))*(abs(w2)+abs(w3))
        v=abs(w1)*(abs(a(is-1,is-1))+abs(a(is,is))+abs(a(is1,is1)))
        if(u+v .eq. v)then
          a(is1,is-1)=0.d0
          a(is2,is-1)=0.d0
          go to 4511
        endif
4510  continue
4511  continue
      a(is2,is)=0.d0
      a(min(is2+1,ie),is)=0.d0
      am =abs(dcmplx(w1,w2))
      am1=abs(dcmplx(am,w3))
      if(is .gt. ib)then
        if(am1 .ne. 0.d0)then
          a(is,is-1)=a(is,is-1)*w1/am1*vx(is)
        endif
      endif
      do 210 i=is-1,ie-2
        i1=i+1
        i2=i1+1
        i3=i2+1
        if(i1 .gt. is)then
          w1=vx(i1)*a(i1,i)
          w2=vx(i2)*a(i2,i)
          a(i2,i)=0.d0
          am =abs(dcmplx(w1,w2))
          if(i3 .le. ie)then
            w3=vx(i3)*a(i3,i)
            a(i3,i)=0.d0
            am1=abs(dcmplx(am,w3))
          else
            w3=0.d0
            am1=0.d0
          endif
        endif
        if(am .ne. 0.d0)then
          c=w1/am
          s=w2/am
          if(abs(c) .gt. abs(s))then
            v1=vx(i1)/c
            v2=c*vx(i2)
            p=s*vx(i1)/v2
            q=s*vx(i2)/v1
            vx(i1)=v1
            if(abs(v2) .gt. vmin .and. abs(v2) .lt. vmax)then
              vx(i2)=v2
              do 220 j=i1,n
                a(i2,j)=a(i2,j)-p*a(i1,j)
                a(i1,j)=a(i1,j)+q*a(i2,j)
 220          continue
              do 230 j=1,min(ie,i3)
                a(j,i1)=a(j,i1)+p*a(j,i2)
                a(j,i2)=a(j,i2)-q*a(j,i1)
 230          continue
              do 231 j=1,n
                w(j,i1)=w(j,i1)+p*w(j,i2)
                w(j,i2)=w(j,i2)-q*w(j,i1)
 231          continue
            else
              vx(i2)=1.d0
              do j=i1,n
                aa=a(i2,j)-p*a(i1,j)
                a(i1,j)=a(i1,j)+q*aa
                a(i2,j)=v2*aa
              enddo
              do j=1,min(ie,i3)
                a(j,i1)=a(j,i1)+p*a(j,i2)
                a(j,i2)=(a(j,i2)-q*a(j,i1))/v2
              enddo
              do j=1,n
                w(j,i1)=w(j,i1)+p*w(j,i2)
                w(j,i2)=(w(j,i2)-q*w(j,i1))/v2
              enddo
            endif
          else
            v1=vx(i2)/s
            v2=vx(i1)*s
            if(i1 .gt. is)then
              a(i1,i)=w2/vx(i2)
            endif
            p=c*vx(i2)/v2
            q=c*vx(i1)/v1
            vx(i1)=v1
            if(abs(v2) .gt. vmin .and. abs(v2) .lt. vmax)then
              vx(i2)=v2
              do 221 j=i1,n
                aa=a(i2,j)
                a(i2,j)=p*aa-a(i1,j)
                a(i1,j)=aa-q*a(i2,j)
 221          continue
              do 232 j=1,min(ie,i3)
                aa=a(j,i1)
                a(j,i1)= p*aa+a(j,i2)
                a(j,i2)= q*a(j,i1)-aa
 232          continue
              do 233 j=1,n
                ww=w(j,i1)
                w(j,i1)= p*ww+w(j,i2)
                w(j,i2)= q*w(j,i1)-ww
 233          continue
            else
              vx(i2)=1.d0
              do j=i1,n
                aa=a(i2,j)
                aa1=p*aa-a(i1,j)
                a(i1,j)=aa-q*aa1
                a(i2,j)=aa1*v2
              enddo
              do j=1,min(ie,i3)
                aa=a(j,i1)
                a(j,i1)= p*aa+a(j,i2)
                a(j,i2)= (q*a(j,i1)-aa)/v2
              enddo
              do j=1,n
                ww=w(j,i1)
                w(j,i1)= p*ww+w(j,i2)
                w(j,i2)= (q*w(j,i1)-ww)/v2
              enddo
            endif
          endif
        endif
        if(am1 .ne. 0.d0)then
          c=am/am1
          s=w3/am1
          if(abs(c) .gt. abs(s))then
            v1=vx(i1)/c
            v2=c*vx(i3)
            p=s*vx(i1)/v2
            q=s*vx(i3)/v1
            vx(i1)=v1
            vx(i3)=v2
            do 240 j=i1,n
              a(i3,j)=a(i3,j)-p*a(i1,j)
              a(i1,j)=a(i1,j)+q*a(i3,j)
240         continue
            do 250 j=1,min(ie,i3+1)
              a(j,i1)=a(j,i1)+p*a(j,i3)
              a(j,i3)=a(j,i3)-q*a(j,i1)
250         continue
            do 251 j=1,n
              w(j,i1)=w(j,i1)+p*w(j,i3)
              w(j,i3)=w(j,i3)-q*w(j,i1)
251         continue
          else
            v1=vx(i3)/s
            v2=vx(i1)*s
            if(i1 .gt. is)then
              a(i1,i)=w3/vx(i3)
            endif
            p=c*vx(i3)/v2
            q=c*vx(i1)/v1
            vx(i1)=v1
            vx(i3)=v2
            do 241 j=i1,n
              aa=a(i3,j)
              a(i3,j)=p*aa-a(i1,j)
              a(i1,j)=aa-q*a(i3,j)
241         continue
            do 252 j=1,min(ie,i3+1)
              aa=a(j,i1)
              a(j,i1)= p*aa+a(j,i3)
              a(j,i3)= q*a(j,i1)-aa
252         continue
            do 253 j=1,n
              ww=w(j,i1)
              w(j,i1)= p*ww+w(j,i3)
              w(j,i3)= q*w(j,i1)-ww
253         continue
          endif
c          if(abs(vx(i1)) .gt. vmax .or. abs(vx(i1)) .lt. vmin .or.
c     $         abs(vx(i2)) .gt. vmax .or. abs(vx(i2)) .lt. vmin)then
c            write(*,*)'tqs-10 ',i1,vx(i1),vx(i2),v1,v2,c,s,
c     $           vx(8437),a(8437,8435)
c          endif
        endif
c       write(*,'(1P4G15.7)')((vx(ii)/vx(j)*a(ii,j),j=1,4),ii=1,4)
210   continue
      if(is .gt. ib)then
        a(is,is-1)=a(is,is-1)/vx(is)
      endif
      iter=iter+1
      if(iter .gt. itmax)then
        write(*,*)' TEIGEN convergence failed. Range =',ib,ie
        write(*,*)
     1  '        Lower right corner =',
     1  a(ie-1,ie-1),vx(ie)/vx(ie-1)*a(ie,ie-1),a(ie,ie)
        a(ie,ie-1)=0.d0
        if(isnan(a(ie,ie)))then
          a(ie,ie)=0.d0
          a(ie-1,ie-1)=0.d0
        endif
        if(isnan(vx(ie)))then
          vx(ie)=1.d0
        endif
        if(isnan(vx(ie-1)))then
          vx(ie-1)=1.d0
        endif
        iter=0
        ie=ie-1
        go to 1
      endif
      go to 2
2000  do 4101 i=1,n
        do 4110 j=max(i-1,1),n
          a(i,j)=a(i,j)*vx(i)
4110    continue
        do 4120 j=1,min(n,i+1)
          a(j,i)=a(j,i)/vx(i)
4120    continue
        do 4130 j=1,n
          w(j,i)=w(j,i)/vx(i)
4130    continue
4101  continue
      paired=.false.
      do 2010 i=1,n-1
        if(paired)then
          paired=.false.
          go to 2010
        endif
        i1=i+1
        if(a(i1,i) .ne. 0.d0)then
          if(abs(a(i1,i)) .gt. abs(a(i,i1)))then
            do 2020 j=i,n
              aa=a(i,j)
              a(i,j)=a(i1,j)
              a(i1,j)=aa
2020        continue
            do 2030 j=1,i1
              aa=a(j,i)
              a(j,i)=a(j,i1)
              a(j,i1)=aa
2030        continue
            do 2040 j=1,n
              ww=w(j,i)
              w(j,i)=w(j,i1)
              w(j,i1)=ww
2040        continue
          endif
          if(a(i1,i) .eq. 0.d0)then
            eig(1,i)=a(i,i)
            eig(2,i)=0.d0
            go to 2010
          endif
          s=a(i,i)-a(i1,i1)
          dec=s**2+4.d0*a(i1,i)*a(i,i1)
          s1=a(i,i)**2+a(i1,i1)**2
          if(dec .ge. -decth*s1)then
            sqrd=sqrt(abs(dec))
            p=2.d0*a(i1,i)/(s+sign(sqrd,s))
            do 2050 j=i,n
              a(i1,j)=a(i1,j)-p*a(i,j)
2050        continue
            do 2060 j=1,i1
              a(j,i)=a(j,i)+p*a(j,i1)
2060        continue
            do 2070 j=1,n
              w(j,i)=w(j,i)+p*w(j,i1)
2070        continue
            eig(1,i)=a(i,i)
            eig(1,i1)=a(i1,i1)
            eig(2,i)=0.d0
            eig(2,i1)=0.d0
            a(i1,i)=0.d0
          else
            p=-s/a(i,i1)*.5d0
            r=sign(
     1      sqrt(abs((a(i1,i)-p*(a(i,i)-a(i1,i1)+p*a(i,i1)))/a(i,i1))),
     1      a(i,i1))
            do 2080 j=i,n
              a(i1,j)=a(i1,j)-p*a(i,j)
              a(i ,j)=a(i ,j)*r
2080        continue
            do 2090 j=1,i1
              a(j,i)=(a(j,i)+p*a(j,i1))/r
2090        continue
            do 2100 j=1,n
              w(j,i)=(w(j,i)+p*w(j,i1))/r
2100        continue
            eig(1,i )=a(i ,i )
            eig(1,i1)=a(i1,i1)
            eig(2,i )=a(i ,i1)
            eig(2,i1)=a(i1,i )
          endif
          paired=.true.
        else
          eig(1,i)=a(i,i)
          eig(2,i)=0.d0
        endif
2010  continue
      if(.not. paired)then
        eig(1,n)=a(n,n)
        eig(2,n)=0.d0
      endif
      do 2110 i=2,n
        if(a(i,i-1) .ne. 0.d0)then
          go to 2110
        endif
        i1=i+1
        if(eig(2,i) .eq. 0.d0)then
          jm=0
          do 2120 j=i-1,1,-1
            if(a(j+1,j) .ne. 0.d0)then
              go to 2120
            endif
            if(eig(2,j) .eq. 0.d0)then
              da=a(j,j)-a(i,i)
              sa=abs(a(j,j))+abs(a(i,i))
              if(abs(da) .gt. threj*sa)then
                p=a(j,i)/da
                do 2130 k=i,n
                  a(j,k)=a(j,k)+p*a(i,k)
 2130           continue
                do 2140 k=1,j
                  a(k,i)=a(k,i)-p*a(k,j)
 2140           continue
                do 2150 k=1,n
                  w(k,i)=w(k,i)-p*w(k,j)
 2150           continue
              else
                jordan=.true.
              endif
            else
              j1=j-1
              u=a(i,i)-a(j1,j1)
              v=a(j1,j)
              r=u**2+v**2
              pk=(a(j1,i)*u+a(j ,i)*v)/r
              pl=(a(j ,i)*u-a(j1,i)*v)/r
              do 2160 k=i,n
                a(j1,k)=a(j1,k)-pk*a(i,k)
                a(j ,k)=a(j ,k)-pl*a(i,k)
2160          continue
              do 2170 k=1,j
                a(k,i)=a(k,i)+pk*a(k,j1)+pl*a(k,j)
2170          continue
              do 2180 k=1,n
                w(k,i)=w(k,i)+pk*w(k,j1)+pl*w(k,j)
2180          continue
            endif
2120      continue
        else
          do 2210 j=i-1,1,-1
            if(a(j+1,j) .ne. 0.d0)then
              go to 2210
            endif
            if(eig(2,j) .eq. 0.d0)then
              x=a(i,i)-a(j,j)
              y=a(i,i1)
              r=x**2+y**2
              pk=(a(j,i )*x+a(j,i1)*y)/r
              pl=(a(j,i1)*x-a(j,i )*y)/r
              do 2220 k=i,n
                a(j,k)=a(j,k)-pk*a(i,k)-pl*a(i1,k)
2220          continue
              do 2230 k=1,j
                a(k,i )=a(k,i )+pk*a(k,j)
                a(k,i1)=a(k,i1)+pl*a(k,j)
2230          continue
              do 2240 k=1,n
                w(k,i )=w(k,i )+pk*w(k,j)
                w(k,i1)=w(k,i1)+pl*w(k,j)
2240          continue
            else
              j1=j-1
              cu=dcmplx(a(j1,j1)-a(i,i),-a(i,i1))
              v =a(j1,j)
              cr=cu**2+v**2
              if(cr .ne. (0.d0,0.d0))then
                ca=dcmplx(a(j1,i),a(j1,i1))
                cc=dcmplx(a(j ,i),a(j ,i1))
                ck=(cu*ca-v*cc)/cr
                cm=(cu*cc+v*ca)/cr
                do 2250 k=i,n
                  a(j1,k)=a(j1,k)+pk*a(i,k)+pl*a(i1,k)
                  a(j,k )=a(j ,k)+pm*a(i,k)+pn*a(i1,k)
2250            continue
                do 2260 k=1,j
                  a(k,i )=a(k,i )-pk*a(k,j1)-pm*a(k,j)
                  a(k,i1)=a(k,i1)-pl*a(k,j1)-pn*a(k,j)
2260            continue
                do 2270 k=1,n
                  w(k,i )=w(k,i )-pk*w(k,j1)-pm*w(k,j)
                  w(k,i1)=w(k,i1)-pl*w(k,j1)-pn*w(k,j)
2270            continue
              endif
            endif
2210      continue
        endif
2110  continue
      if(jordan)then
        do 2510 i=n,2,-1
          if(eig(2,i) .eq. 0.d0)then
            jm=0
            do 2520 j=i-1,1,-1
              if(eig(1,j) .eq. eig(1,i) .and. eig(2,j) .eq. 0.d0)then
                r=a(j,i)
                if(abs(r) .gt. 1.d-30)then
                  if(jm .eq. 0)then
                    r1=sqrt(abs(r))
                    r2=r/r1
                    do 2530 k=j,n
                      a(j,k)=a(j,k)/r1
2530                continue
                    do 2540 k=1,j
                      a(k,j)=a(k,j)*r1
2540                continue
                    do 2550 k=1,n
                      w(k,j)=w(k,j)*r1
2550                continue
                    do k=1,i
                      a(k,i)=a(k,i)/r2
                    enddo
                    do k=i,n
                      a(i,k)=a(i,k)*r2
                    enddo
                    do k=1,n
                      w(k,i)=w(k,i)/r2
                    enddo
                    jm=j
                  else
                    do 2560 k=jm,n
                      a(j,k)=a(j,k)-r*a(jm,k)
2560                continue
                    do 2570 k=1,j
                      a(k,jm)=a(k,jm)+r*a(k,j)
2570                continue
                    do 2580 k=1,n
                      w(k,jm)=w(k,jm)+r*w(k,j)
2580                continue
                  endif
                else
                  a(j,i)=0.d0
                endif
              endif
2520        continue
          endif
2510    continue
      endif
c     do 1000 k=1,n
c       write(*,'(1X,:1P10G13.6)')(a(k,j),j=1,n)
c1000  continue
c     do 1010 k=1,n
c       write(*,'(1X,:1P10G13.6)')(w(k,j),j=1,n)
c1010  continue
c     write(*,*)
      if(mod(n,2) .eq. 0)then
        ii=0
        do 3010 i=1,n
          if(eig(2,i) .gt. 0.d0 .and. ii .ne. 0)then
            do 3020 j=1,n
              ww=w(j,i-1)
              w(j,i-1)=w(j,i  )
              w(j,i  )=w(j,i+1)
              w(j,i+1)=ww
3020        continue
            ee=eig(1,i-1)
            eig(1,i-1)=eig(1,i  )
            eig(1,i  )=eig(1,i+1)
            eig(1,i+1)=ee
            ee=eig(2,i-1)
            eig(2,i-1)=eig(2,i  )
            eig(2,i  )=eig(2,i+1)
            eig(2,i+1)=ee
          endif
          ii=1-ii
3010    continue
      endif
      return
      end
