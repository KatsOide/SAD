      subroutine tcftr(ca,n,conj)
      implicit none
      integer*4 idim
      parameter (idim=30)
      integer*4 n,m(idim),i,ln,j,l,ka,kb,k
      complex*16 ca(n),cex(idim),ch,cw,cw1
      logical*4 conj
      data m/
     $     1,2,4,8,16,
     $     32,64,128,256,512,
     $     1024,2048,4096,8192,16384,
     $     32768,65536,131072,262144,524288,
     $     1048576,2097152,4194304,8388608,16777216,
     $     33554432,67108864,134217728,268435456,536870912/
      data cex/
     $     (-1.d0                 ,0.d0),
     $     (0.d0                  ,1.d0),
     $     (0.707106781186547524d0,0.707106781186547524d0),
     $     (0.923879532511286756d0,0.382683432365089772d0),
     $     (0.980785280403230449d0,0.195090322016128268d0),
     $     (0.995184726672196886d0,0.098017140329560602d0),
     $     (0.998795456205172393d0,0.0490676743274180143d0),
     $     (0.99969881869620422d0 ,0.024541228522912288d0),
     $     (0.999924701839144541d0,0.01227153828571992608d0),
     $     (0.999981175282601143d0,0.00613588464915447536d0),
     $     (0.999995293809576172d0,0.00306795676296597627d0),
     $     (0.99999882345170191d0 ,0.001533980186284765612d0),
     $     (0.999999705862882219d0,0.000766990318742704527d0),
     $     (0.999999926465717851d0,0.000383495187571395589d0),
     $     (0.999999981616429294d0,0.000191747597310703307d0),
     $     (0.999999995404107313d0,0.0000958737990959773459d0),
     $     (0.999999998851026828d0,0.0000479368996030668845d0),
     $     (0.999999999712756707d0,0.0000239684498084182187d0),
     $     (0.999999999928189177d0,0.00001198422490506970642d0),
     $     (0.999999999982047294d0,5.99211245264242784d-6),
     $     (0.999999999995511824d0,2.99605622633466075d-6),
     $     (0.999999999998877956d0,1.498028113169011229d-6),
     $     (0.999999999999719489d0,7.49014056584715721d-7),
     $     (0.999999999999929872d0,3.74507028292384124d-7),
     $     (0.999999999999982468d0,1.87253514146195345d-7),
     $     (0.999999999999995617d0,9.36267570730980828d-8),
     $     (0.999999999999998904d0,4.68133785365490927d-8),
     $     (0.999999999999999726d0,2.34066892682745528d-8),
     $     (0.999999999999999932d0,1.170334463413727718d-8),
     $     (0.999999999999999983d0,5.85167231706863869d-9)
     $     /
      do 10 i=1,idim-1
        if(n .eq. m(i+1))then
          ln=i
          go to 1
        endif
 10   continue
      n=0
      return
 1    j=1
      do 20 i=1,n-1
        if( j .gt. i )then
*     VOPTION NOFVAL
          ch=ca(j)
          ca(j)=ca(i)
          ca(i)=ch
        endif
        l=ln
 2      if(j .le. m(l))then
          j=j+m(l)
        else
          j=j-m(l)
          l=l-1
          go to 2
        endif
 20   continue
      do 110 i=1,ln
        cw=(1.d0,0.d0)
        if(conj)then
          cw1=conjg(cex(i))
        else
          cw1=cex(i)
        endif
        do 120 j=1,m(i)
*     VOPTION NOFVAL
          do 130 k=0,n-1,m(i+1)
            ka=k+j
            kb=ka+m(i)
            ch=cw*ca(kb)
            ca(kb)=ca(ka)-ch
            ca(ka)=ca(ka)+ch
 130      continue
          cw=cw*cw1
 120    continue
 110  continue
      return
      end

      subroutine trftr(a,m,conj)
      use iso_c_binding
      implicit none
      integer*4 m,i,i2,j2
      real*8 , target::a(0:m-1)
      real*8 pi,w,ci,si,ai2,ai21,aj2,aj21,a1
      parameter (pi=3.14159265358979324d0)
      logical*4 conj
      complex*16 , pointer ::ca(:)
      call c_f_pointer(c_loc(a),ca,[m/2])
      call tcftr(ca,m/2,conj)
      if(conj)then
        w=-2.d0*pi/m
      else
        w= 2.d0*pi/m
      endif
      do i=1,m/4
        ci=cos(w*i)
        si=sin(w*i)
        i2  =i*2
        j2  =m-i2
        ai2 =a(i2)
        ai21=a(i2+1)
        aj2 =a(j2)
        aj21=a(j2+1)
        a(i2  )=.5d0*(ai2 +aj2 +ci*(ai21+aj21)+si*(ai2-aj2))
        a(j2  )=.5d0*(ai2 +aj2 -ci*(ai21+aj21)-si*(ai2-aj2))
        a(i2+1)=.5d0*(ai21-aj21+si*(ai21+aj21)-ci*(ai2-aj2))
        a(j2+1)=.5d0*(aj21-ai21+si*(ai21+aj21)-ci*(ai2-aj2))
      enddo
      a1=a(1)
      a(1)=a(0)-a1
      a(0)=a(0)+a1
      return
      end

      subroutine tftrr(a,m,conj)
c
c Inversion of trftr
c tftrr(trftr(f)) = m*f
c
      use iso_c_binding
      implicit none
      integer*4 m,i,i2,j2
      real*8 ,target:: a(0:m-1)
      real*8 pi,w,ci,si,ai2,ai21,aj2,aj21,a1
      parameter (pi=3.14159265358979324d0)
      logical*4 conj
      complex*16 ,pointer :: ca(:)
      if(conj)then
        w= 2.d0*pi/m
      else
        w=-2.d0*pi/m
      endif
      do i=1,m/4
        ci=cos(w*i)
        si=sin(w*i)
        i2  =i*2
        j2  =m-i2
        ai2 =a(i2)
        ai21=a(i2+1)
        aj2 =a(j2)
        aj21=a(j2+1)
        a(i2  )=ai2 +aj2 -ci*(ai21+aj21)+si*(ai2-aj2)
        a(j2  )=ai2 +aj2 +ci*(ai21+aj21)-si*(ai2-aj2)
        a(i2+1)=ai21-aj21+si*(ai21+aj21)+ci*(ai2-aj2)
        a(j2+1)=aj21-ai21+si*(ai21+aj21)+ci*(ai2-aj2)
      enddo
      a1=a(1)
      a(1)=a(0)-a1
      a(0)=a(0)+a1
      call c_f_pointer(c_loc(a),ca,[m/2])
      call tcftr(ca,m/2,conj)
      return
      end
