      module gammaf
      use mathfun
      use macmath

      integer*4 ,parameter,private ::ibmax=100,nstg=100
      real*8 ,parameter,private ::hg=0.2d0
      complex*16 , parameter,private ::one=(1.d0,1.d-18)
      real*8 ,parameter,private:: berbf(0:ibmax)=[
     $     1.d0,-0.5d0,0.0833333333333333333d0,0.d0,
     -     -0.00138888888888888889d0,0.d0,
     -     0.0000330687830687830688d0,0.d0,
     -     -8.2671957671957672d-7,0.d0,
     -     2.0876756987868099d-8,0.d0,
     -     -5.28419013868749318d-10,0.d0,
     -     1.33825365306846788d-11,0.d0,
     -     -3.38968029632258287d-13,0.d0,
     -     8.58606205627784456d-15,0.d0,
     -     -2.17486869855806187d-16,0.d0,
     -     5.50900282836022952d-18,0.d0,
     -     -1.39544646858125233d-19,0.d0,
     -     3.53470703962946747d-21,0.d0,
     -     -8.95351742703754685d-23,0.d0,
     -     2.26795245233768306d-24,0.d0,
     $     -5.74479066887220245d-26,0.d0,
     $     1.4551724756148649d-27,0.d0,
     $     -3.68599494066531018d-29,0.d0,
     $     9.33673425709504467d-31,0.d0,
     $     -2.36502241570062993d-32,0.d0,
     $     5.9906717624821343d-34,0.d0,
     $     -1.51745488446829026d-35,0.d0,
     $     3.84375812545418823d-37,0.d0,
     $     -9.73635307264669104d-39,0.d0,
     $     2.46624704420068096d-40,0.d0,
     $     -6.24707674182074369d-42,0.d0,
     $     1.58240302446449143d-43,0.d0,
     $     -4.00827368594893597d-45,0.d0,
     $     1.01530758555695563d-46,0.d0,
     $     -2.57180415824187175d-48,0.d0,
     $     6.51445603523381493d-50,0.d0,
     $     -1.65013099068965246d-51,0.d0,
     $     4.17983062853947589d-53,0.d0,
     $     -1.05876346677029088d-54,0.d0,
     -     2.6818791912607706661d-56,0.d0,
     -     -6.7932793511074212095d-58,0.d0,
     -     1.7207577616681404905d-59,0.d0,
     -     -4.3587303293488938434d-61,0.d0,
     -     1.1040792903684666751d-62,0.d0,
     -     -2.7966655133781345072d-64,0.d0,
     -     7.0840365016794701985d-66,0.d0,
     -     -1.7944074082892240666d-67,0.d0,
     -     4.5452870636110961071d-69,0.d0,
     -     -1.1513346631982051813d-70,0.d0,
     -     2.9163647710923613547d-72,0.d0,
     -     -7.3872382634973375626d-74,0.d0,
     -     1.8712093117637953062d-75,0.d0,
     -     -4.7398285577617994055d-77,0.d0,
     -     1.200612599335450652d-78,0.d0,
     -     -3.041187241514292383d-80]
      real*8 ,parameter,private :: stg(0:nstg)=[
     $     0.577215664901532861d0, -0.0728158454836767249d0,
     $     -0.00969036319287231848d0, 0.00205383442030334587d0, 
     $     0.00232537006546730006d0, 0.000793323817301062702d0, 
     $     -0.000238769345430199610d0, -0.000527289567057751046d0, 
     $     -0.000352123353803039510d0, -0.0000343947744180880482d0, 
     $     0.000205332814909064795d0, 0.000270184439543903527d0, 
     $     0.000167272912105140193d0, -0.0000274638066037601589d0, 
     $     -0.000209209262059299946d0, -0.000283468655320241447d0, 
     $     -0.000199696858308969775d0, 0.0000262770371099183367d0, 
     $     0.000307368408149252827d0, 0.000503605453047355629d0, 
     $     0.000466343561511559449d0, 0.000104437769756000116d0, 
     $     -0.000541599582203997702d0, -0.00124396209040824578d0, 
     $     -0.00158851127890356156d0, -0.00107459195273848882d0,
     $     0.000656803518637154432d0, 0.00347783691361853821d0, 
     $     0.00640006853170062946d0, 0.00737115177047223913d0, 
     $     0.00355772885557316095d0, -0.00751332599781522893d0, 
     $     -0.0257037291084204018d0, -0.0451067341080802199d0, 
     $     -0.0511269280215084644d0, -0.0203730436038613127d0, 
     $     0.0724821588168113337d0, 0.236026382274301503d0,
     $     0.428963446384809153d0, 0.517921842692923719d0,
     $     0.248721559394615465d0,-0.719574846901300351d0,
     -     -2.63879492733573454d0,-5.26493031235502383d0,
     -     -7.18874588950352728d0,-5.07234458991637249d0,
     -     6.60991560909696581d0,34.0397749821587482d0,
     -     78.682479763242585d0,125.844387631978469d0,
     -     126.823602651322717d0,-19.1969118730278558d0,
     -     -463.188923026716811d0,-1340.65914437689219d0,
     -     -2572.45474040443552d0,-3457.14120864538995d0,
     -     -2055.2758162319743d0,5372.28221320319128d0,
     -     24019.3893776069882d0,57424.3192969640755d0,
     -     98543.2545901460421d0,111670.957814941079d0,
     -     5333.66521050076434d0,-390972.687313396396d0,
     -     -1.30318071253251981d6,-2.84507655260861212d6,
     -     -4.5405266097377241d6,-4.3419051390015162d6,
     -     2.87156694597246051d6,2.66049085466867709d7,
     -     7.93216631192990606d7,1.66215134046825436d8,
     -     2.55153258308238993d8,2.12655631691854037d8,
     -     -2.98767089431166184d8,-1.9194874277328024d9,
     -     -5.51557425812922027d9,-1.1483450987926256d10,
     -     -1.75701522777772639d10,
     -     -1.39610214580125182d10,
     -     2.51634410107906368d10,1.51058510830039652d11,
     -     4.37904431188825412d11,9.31706846884597532d11,
     -     1.47209981910894167d12,1.25904496781436343d12,
     -     -1.95881022472885916d12,
     -     -1.29515454993526074d13,
     -     -3.92971538821761472d13,
     -     -8.75304473973060639d13,
     -     -1.47161049421415044d14,
     -     -1.47797700074868832d14,1.294632141268515d14,
     -     1.18856297595292327d15,3.9206846269612364d15,
     -     9.341670850346208d15,1.70752474522235763d16,
     -     2.07179835426089506d16,
     -     -2.85430784969050842d15,
     -     -1.12584810772062443d17,
     -     -4.25340157170802696d17]
      real*8 ,parameter,private ::alogsqrt2pi=log(m_sqrt2pi),
     $     alogpi=1.1447298858494001741d0,
     $     ath1=0.5d0,ath2=0.5d0,cth=59.d0,ztlim=18.d0,
     $     dz0=-.5d0*log(m_2pi),zt0=-0.5d0,zt2=m_pi**2/6.d0,
     $     mm_pi=-m_pi,epso=5.d-17**2/4.d0,
     $     ggamma=4.7421875d0,
     $     gc0=0.99999999999999709182d0,
     $     gc1=57.156235665862923517d0,
     $     gc2=-59.597960355475491248d0,
     $     gc3=14.136097974741747174d0,
     $     gc4=-0.49191381609762019978d0,
     $     gc5=.33994649984811888699d-4,
     $     gc6=.46523628927048575665d-4,
     $     gc7=-.98374475304879564677d-4,
     $     gc8=.15808870322491248884d-3,
     $     gc9=-.21026444172410488319d-3,
     $     gc10=.21743961811521264320d-3,
     $     gc11=-.16431810653676389022d-3,
     $     gc12=.84418223983852743293d-4,
     $     gc13=-.26190838401581408670d-4,
     $     gc14=.36899182659531622704d-5,

c$$$     $     ggamma=9.d0,
c$$$     $     gc0=1.000000000000000174663d0,
c$$$     $     gc1=5716.400188274341379136d0,
c$$$     $     gc2=-14815.30426768413909044d0,
c$$$     $     gc3=14291.49277657478554025d0,
c$$$     $     gc4=-6348.160217641458813289d0,
c$$$     $     gc5=1301.608286058321874105d0,
c$$$     $     gc6=-108.1767053514369634679d0,
c$$$     $     gc7=2.605696505611755827729d0,
c$$$     $     gc8=-0.7423452510201416151527d-2,
c$$$     $     gc9=0.5384136432509564062961d-7,
c$$$     $     gc10=-0.4023533141268236372067d-8

c$$$     $     ggamma=7.d0,
c$$$     $     gc0=0.99999999999980993227684700473478d0,
c$$$     $     gc1=676.520368121885098567009190444019d0,
c$$$     $     gc2=-1259.13921672240287047156078755283d0,
c$$$     $     gc3=771.3234287776530788486528258894d0,
c$$$     $     gc4=-176.61502916214059906584551354d0,
c$$$     $     gc5=12.507343278686904814458936853d0,
c$$$     $     gc6=-0.13857109526572011689554707d0,
c$$$     $     gc7=9.984369578019570859563d-6,
c$$$     $     gc8=1.50563273514931155834d-7,

c$$$     $     gamma=5,
c$$$     $     c0=0.999999999999997524d0,
c$$$     $      c1=76.1800917309326077d0,
c$$$     $      c2=-86.5053204008552d0,
c$$$     $      c3=24.01409906379226027d0,
c$$$     $      c4=-1.231743354454618365d0,
c$$$     $      c5=0.001216872118636531519d0,
c$$$     $      c6=-0.00001408915744128554778d0,
c$$$     $      c7=4.00491935010387864d-6,
c$$$     $      c8=-4.93961492826482964d-7,
     $     algg=m_sqrt2pi/exp(ggamma)
      integer*4 ,parameter,private :: nolog=8,
     $     nogam=800,nocgam=2000,
     $     nopg=48,nozt=16,nozt2=30,nogam2=3000,nopl=500
      real*8 ,parameter,private ::sconf=2.d0**55,xmp=30.d0,xmth=0.96d0,
     $     xmth1=1.1d0,xmeps=1.d-6,
     $     dpmg1=-1.4132139976024971836d0,
     $     ddpmg1=-0.10668467023117844470d0,
     $     dpzmg1=1.1217084897192276354d0,
     $     dpmg2=-0.295789424393121511d0,
     $     ddpmg2=1.58338508513914499d0,
     $     dpzmg2=0.853223513334086317d0,
     $     epsba=1.d-10,epsba1=1.d-7
      complex*16 ,parameter,private :: xs1=dcmplx(0.5d0,sqrt(.75d0)),
     $     xs2=dcmplx(0.5d0,-sqrt(.75d0)),czero=(0.d0,0.d0),
     $     cone=(1.d0,0.d0)
      real*8 ,parameter,private :: epsx=0.3d0
      integer*4 ,parameter,private :: itmaxg=67
      real*8 ,parameter,private :: zimth=1.4d0,
     $     chgth=8.d0,veryl=1.d250,arth=0.8d0,vlim=1.d30,
     $     slimpl=1.d-4,slims=4.5d0,lgzlim=1.4d0,plith=4.d0,
     $     plrth=4.d0,pllzth=1.5d0*m_pi,simth=4.d0
      integer*4 ,parameter,private ::kg2max=20
      complex*16 ,private::cpgn,cpgz
      real*8 ,private::pgnc,pgn,pgz,ag(itmaxg)

      contains

      subroutine aginit
      implicit none
      integer*4 i
      real*8 an
      if(ag(1) == 0.d0)then
        an=-1.d0
        do i=1,itmaxg
          an=an+2.d0
          ag(i)=exp(-(an*hg)**2)
        enddo
      endif
c      write(*,*)'agint ',itmaxg,ag(itmaxg)
      return
      end subroutine

      real*8 recursive  function bernb(n) result(f1)
      implicit none
      real*8 ,intent(in):: n
      if(dble(n) .ne. anint(n))then
        f1=1.d0/0.d0
      else
        f1=bernbf(nint(n))*factorial(n)
      endif
      return
      end

      real*8 recursive pure function bernbf(n) result(f1)
      implicit none
      integer*4 ,intent(in):: n
      real*8 u,rn1,k,k1
      integer*4 n1,i
      if(n < 0)then
        f1=0.d0
      elseif(n <= ibmax)then
        f1=berbf(n)
      elseif(mod(n,2) /= 0)then
        f1=0.d0
      else
        n1=n/2
        rn1=dble(n1+1)
        u=pochh(dble(n),1.d0-rn1)
        f1=rn1*bernbf(n1)*u
        k=0.d0
        do i=1,n1-1
          k1=k+1.d0
          u=u*(rn1-1.d0+k)*(rn1-k)/k1
          f1=f1+(rn1+k)*bernbf(n1+i)*u
          k=k1
        enddo
        f1=-f1/rn1/(n+1)
      endif
      return
      end function

      complex*16 pure function berpol(n,x) result(f1)
      implicit none
      integer*4 ,intent(in):: n
      complex*16 ,intent(in):: x
      complex*16 u,df,x2
      integer*4 k,no
      if(x == czero)then
        f1=berbf(n)*factorial(dble(n))
      elseif(imag(x) .eq. 0.d0)then
        f1=dcmplx(berpolr(n,dble(x)),0.d0)
      else
        select case (n)
        case (0)
          f1=berbf(0)
        case(1)
          f1=berbf(1)+berbf(0)*x
        case default
          no=8
          if(mod(n,2) == 0)then
            if(n == 2)then
              f1=2.d0*(berbf(2)+x*(berbf(1)+x*berbf(0)/2.d0))
            else
              x2=x**2
              u=factorial(dble(n))
c              write(*,'(a,1p10g12.4)')'bp-0',n,u
              f1=u*bernbf(n)
              do k=2,n-2,2
                u=u*x2/dble(k*(k-1))
                df=bernbf(n-k)*u
                f1=f1+df
                no=no+8
                if(abs((abs(dble(df))+abs(imag(df)))/f1)**2
     $               <= no*epso)then
c                  u=u*x2**((n-2-k)/2)/pochh(dble(k+1),dble(n-2-k))
                  return
                endif
              enddo
              f1=f1+u*x/dble(n-1)*(berbf(1)+berbf(0)*x/dble(n))
c              write(*,'(a,1p10g12.4)')':   ',u,f1
            endif
          else
            if(n == 3)then
              f1=6.d0*x*(berbf(2)+x*(berbf(1)/2.d0+x*berbf(0)/6.d0))
            else
              x2=x**2
              u=factorial(dble(n))
              f1=u*bernbf(n-1)*x
              u=u*x
              do k=3,n-2,2
                u=u*x2/dble(k*(k-1))
                df=bernbf(n-k)*u
                f1=f1+df
                no=no+8
                if(abs((abs(dble(df))+abs(imag(df)))/f1)**2
     $               <= no*epso)then
c                  u=u*x2**((n-2-k)/2)/pochh(dble(k+1),dble(n-2-k))
                  return
                endif
              enddo
              f1=f1+u*x/dble(n-1)*(berbf(1)+berbf(0)*x/dble(n))
            endif
          endif
        end select
      endif
      return
      end

      real*8 pure function berpolr(n,x) result(f1)
      implicit none
      integer*4 ,intent(in):: n
      real*8 ,intent(in):: x
      real*8 u,df
      integer*4 k,no
      if(x == 0.d0)then
        f1=berbf(n)*factorial(dble(n))
      else
        select case (n)
        case (0)
          f1=berbf(0)
        case(1)
          f1=berbf(1)+berbf(0)*x
        case default
          no=8
          if(mod(n,2) == 0)then
            if(n == 2)then
              f1=2.d0*(berbf(2)+x*(berbf(1)+x*berbf(0)/2.d0))
c     write(*,'(a,1p10g12.4)')'bp ',n,x,f1
            else
              u=factorial(dble(n))
              f1=u*bernbf(n)
              do k=2,n-2,2
                u=u*x**2/dble(k*(k-1))
                df=bernbf(n-k)*u
                f1=f1+df
                no=no+8
                if((df/f1)**2 <= no*epso)then
                  return
                endif
              enddo
              f1=f1+u*x/dble(n-1)*(berbf(1)+berbf(0)*x/dble(n))
            endif
          else
            if(n == 3)then
              f1=6.d0*x*(berbf(2)+x*(berbf(1)/2.d0+x*berbf(0)/6.d0))
            else
              u=factorial(dble(n))
              f1=u*bernbf(n-1)*x
              u=u*x
              do k=3,n-2,2
                u=u*x**2/dble(k*(k-1))
                df=bernbf(n-k)*u
                f1=f1+df
                no=no+8
                if((df/f1)**2 <= no*epso)then
                  return
                endif
              enddo
              f1=f1+u*x/dble(n-1)*(berbf(1)+berbf(0)*x/dble(n))
            endif
          endif
        end select
      endif
      return
      end

      real*8 pure function stgn(n) result(f1)
      implicit none
      integer*4 ,intent(in):: n
      real*8 ,parameter :: ln2=log(2.d0)
      real*8 ak,df
      integer*4 k,no
      if(n < 0)then
        f1=0.d0
      elseif(n <= nstg)then
        f1=stg(n)
      else
        f1=-berpolr(n+1,0.d0)
        k=2
        ak=2.d0
        no=10
        do
          df=(berpolr(n+1,log(ak)/ln2)/ak
     $         -berpolr(n+1,log(ak+1.d0)/ln2)/(ak+1.d0))
          f1=f1+df
          no=no+nolog*2+4
c          write(*,'(a,i7,1p10g12.4)')'stgn ',n,ak,f1,df
          if((df/f1)**2 <= no*epso)then
            f1=f1*ln2**n/(n+1)
            exit
          endif
          k=k+2
          ak=ak+2.d0
        enddo
      endif
      return
      end function

      real*8 pure function rgamma(x) result(f1)
      implicit none
      real*8 ,intent(in):: x
      f1=gamma(x)
      return
      end function

      real*8 pure function factorial(x)
      implicit none
c     Including pi = m_pi
      real*8 ,intent(in):: x
      integer*4 ix
      integer*4 ,parameter ::lt=23
      real*8 ,parameter :: fac(0:lt)=[
     $1d0,
     $1d0,
     $2d0,
     $6d0,
     $24d0,
     $120d0,
     $720d0,
     $5040d0,
     $40320d0,
     $362880d0,
     $3628800d0,
     $39916800d0,
     $479001600d0,
     $6227020800d0,
     $87178291200d0,
     $1307674368000d0,
     $20922789888000d0,
     $355687428096000d0,
     $6402373705728000d0,
     $121645100408832000d0,
     $2432902008176640000d0,
     $51090942171709440000d0,
     $1124000727777607680000d0,
     $25852016738884976640000d0]
      if(x < 0.d0)then
        factorial=exp(-log_gamma(1.d0-x))*m_pi*x/sinp(x)
      else
        ix=nint(x)
        if(ix <= lt .and. x-ix == 0.d0)then
          factorial=fac(ix)
        else
          factorial=exp(log_gamma(x+1.d0))
        endif
      endif
      return
      end

      real*8 pure function aloggamma(x) result(f)
      implicit none
      real*8 ,intent(in):: x
      f=log_gamma(x)
      return
      end function

      real*8 pure function aloggamma1(x) result(f)
      implicit none
      real*8 ,intent(in):: x
      f=log_gamma(x+1.d0)
      return
      end function
c$$$c Lanczos' formula
c$$$c$$$      use macmath
c$$$      implicit none
c$$$c     Including pi = m_pi
c$$$      real*8 ,intent(in):: x
c$$$      real*8 x1
c$$$      if(x < 0.d0)then
c$$$        x1=-x
c$$$        aloggamma1=log(pi*x/sinp(x))-
c$$$     $       ((x1+.5d0)*log(x1+ggamma+.5d0)-(x1+ggamma+.5d0)
c$$$     $       +alogsqrt2pi+log(
c$$$     $       gc0+gc1/(1.d0+x1)+gc2/(2.d0+x1)+gc3/(3.d0+x1)
c$$$     $       +gc4/(4.d0+x1)+gc5/(5.d0+x1)+gc6/(6.d0+x1)
c$$$     $       +gc7/(7.d0+x1)+gc8/(8.d0+x1)))
c$$$c     $       +gc9/(9.d0+x1)+gc10/(10.d0+x1))
c$$$      else
c$$$        aloggamma1=(x+.5d0)*log(x+ggamma+.5d0)-(x+ggamma+.5d0)
c$$$     $       +alogsqrt2pi+log(
c$$$     $       gc0+gc1/(1.d0+x)+gc2/(2.d0+x)+gc3/(3.d0+x)
c$$$     $       +gc4/(4.d0+x)+gc5/(5.d0+x)+gc6/(6.d0+x)
c$$$     $       +gc7/(7.d0+x)+gc8/(8.d0+x))
c$$$c     $       +gc9/(9.d0+x)+gc10/(10.d0+x))
c$$$      endif
c$$$      return
c$$$      end

      complex*16 recursive pure function clg2(z) result(f1)
c Lanczos' formula
      implicit none
      dimension f1(2)
      complex*16 ,intent(in):: z
      complex*16 c1(2),c2(2)
      if(abs(imag(z)) > lgzlim)then
        c1=clg2(.5d0*z)
        c2=clg2(.5d0*z+.5d0)
        f1(1)=c1(1)+c2(1)+(z-1.d0)*m_ln2
        f1(2)=c1(2)*c2(2)/m_sqrtpi
      else
        f1(1)=(z-.5d0)*log((z+ggamma-.5d0)/m_e)
        f1(2)=algg*(
     $       gc0+gc1/z      +gc2/(1.d0+z)
     $       +gc3/(2.d0+z)  +gc4/(3.d0+z)
     $       +gc5/(4.d0+z)  +gc6/(5.d0+z)
     $       +gc7/(6.d0+z)  +gc8/(7.d0+z)
     $       +gc9/(8.d0+z)  +gc10/(9.d0+z)
     $       +gc11/(10.d0+z)+gc12/(11.d0+z)
     $       +gc13/(12.d0+z)+gc14/(13.d0+z)
     $       )
      endif
      return
      end function

      complex*16 recursive pure function clg1(z) result(f1)
c Lanczos' formula
      implicit none
      dimension f1(2)
      complex*16 ,intent(in):: z
      complex*16 c1(2),c2(2)
      if(abs(imag(z)) > lgzlim)then
        c1=clg2(.5d0*z+.5d0)
        c2=clg2(.5d0*z+1.d0)
        f1(1)=c1(1)+c2(1)+z*m_ln2
        f1(2)=c1(2)*c2(2)/m_sqrtpi
      else
        f1(1)=(z+.5d0)*log((z+ggamma+.5d0)/m_e)
        f1(2)=algg*(
     $       gc0+gc1/(1.d0+z)+gc2/(2.d0+z)
     $       +gc3/(3.d0+z)  +gc4/(4.d0+z)
     $       +gc5/(5.d0+z)  +gc6/(5.d0+z)
     $       +gc7/(7.d0+z)  +gc8/(8.d0+z)
     $       +gc9/(8.d0+z)  +gc10/(9.d0+z)
     $       +gc11/(10.d0+z)+gc12/(11.d0+z)
     $       +gc13/(12.d0+z)+gc14/(13.d0+z)
     $       )
      endif
      return
      end function

      complex*16 pure function cgamma(s) result(f1)
      implicit none
      complex*16 ,intent(in):: s
      complex*16 c1(2)
      if(imag(s) == 0.d0)then
        f1=dcmplx(gamma(dble(s)),0.d0)
      elseif(dble(s) < 0.d0)then
        f1=-m_pi*cgammai(-s)/s/csinp(s)
      else
        c1=clg2(s)
        f1=exp(c1(1))*c1(2)
      endif
      return
      end function

      real*8 pure function gammai(x)
      implicit none
      real*8 ,intent(in):: x
      if(anint(x) == x .and. x <= 0.d0)then
        gammai=0.d0
      else
        gammai=1.d0/gamma(x)
      endif
      return
      end function

      complex*16 pure function cgammai(x)
      implicit none
      complex*16 ,intent(in):: x
      if(imag(x) == 0.d0)then
        cgammai=dcmplx(gammai(dble(x)),0.d0)
      else
        cgammai=1.d0/cgamma(x)
      endif
      return
      end function

      complex*16 recursive pure function cloggamma(z) result(f1)
c Lanczos' formula
      implicit none
c     Including pi = m_pi
      complex*16 ,intent(in):: z
      complex*16 z1,c1(2)
      real*8 x1
      if(dble(z) < 0.d0)then
        if(imag(z) == 0.d0)then
          x1=1.d0-dble(z)
          f1=dcmplx(-log_gamma(x1)+alogpi-log(abs(sinp(x1))),
     $         -m_pi*floor(x1))
        else
          z1=zeroim(1.d0-z)
          f1=-cloggamma(z1)+alogpi-log(csinp(z1))
     $         +dcmplx(0.d0,m_2pi*floor(.5d0*(dble(z)+.5d0)))
        endif
      else
        c1=clg2(z)
        f1=c1(1)+log(c1(2))
      endif
      return
      end

      complex*16 recursive pure function cloggamma1(z) result(f1)
c Lanczos' formula
      implicit none
c     Including pi = m_pi
      complex*16 ,intent(in):: z
      complex*16 z1,c1(2)
      real*8 x1
      if(dble(z) < 0.d0)then
        if(imag(z) == 0.d0)then
          x1=-dble(z)
          f1=dcmplx(-log_gamma(x1)+alogpi-log(abs(sinp(x1))),
     $         -m_pi*floor(x1))
        else
          z1=zeroim(-z)
          f1=-cloggamma(z1)+alogpi-log(csinp(z1))
     $         +dcmplx(0.d0,m_2pi*floor(.5d0*(dble(z)+.5d0)))
        endif
      else
        c1=clg1(z)
        f1=c1(1)+log(c1(2))
      endif
      return
      end

      complex*16 pure function cfactorial(x)
      implicit none
      complex*16 ,intent(in):: x
      cfactorial=exp(cloggamma1(x))
      return
      end function

      complex*16 recursive pure function cpochh(a,n) result(f)
      implicit none
      complex*16 ,intent(in):: a,n
      complex*16 an
      if(imag(a) == 0.d0 .and. imag(n) == 0.d0)then
        f=dcmplx(pochh(dble(a),dble(n)),0.d0)
      elseif(n == czero)then
        f=cone
      elseif(n == cone)then
        f=a
      else
        an=a+n
        if(imag(an) == 0.d0 .and.
     $       anint(dble(an)) == dble(an) .and. dble(an) <= 0.d0
     $       .and. (anint(dble(a)) /= dble(a)
     $       .or. dble(a) < 0.d0))then
          f=czero
        else
          f=exp(cloggamma(an)-cloggamma(a))
        endif
      endif
      return
      end function

      real*8 recursive pure function pochh(a,n) result(f)
      implicit none
      real*8 ,intent(in):: a,n
      real*8 ,parameter ::xth=100.d0
      real*8 an
      if(n == 0.d0)then
        f=1.d0
      elseif(n == 1.d0)then
        f=a
      elseif(a == 0.d0)then
        if(n > 0.d0)then
          f=gamma(n)
        else
          f=gammai(1.d0-n)
        endif
      else
        an=a+n
        if(an > 0.d0 .and. a > 0.d0)then
          f=exp(log_gamma(an)-log_gamma(a))
        elseif(an /= anint(an) .or. a /= anint(a))then
          f=gamma(an)*gammai(a)
        elseif(anint(n) == n)then
          if(n > 0.d0)then
            if(an > 0.d0)then
              f=0.d0
            else
              f=(-1.d0)**nint(n)*pochh(1.d0-an,n)
            endif
          else
            f=(-1.d0)**nint(n)/pochh(1.d0-a,-n)
          endif
        else
          f=0.d0
        endif
      endif
      return
      end function

      complex*16 function cbeta2(a,b) result(f1)
      implicit none
      complex*16 ,intent(in):: a,b
      if(imag(a) == 0.d0 .and. imag(b) == 0.d0 .and.
     $     dble(a) > 0.d0 .and. dble(b) > 0.d0)then
        f1=dcmplx(exp(log_gamma(dble(a))+log_gamma(dble(b))
     $       -log_gamma(dble(a)+dble(b))),
     $       0.d0)
      else
        f1=exp(cloggamma(a)+cloggamma(b)-cloggamma(a+b))
      endif
      return
      end function

      real*8 recursive function zeta(s) result(f1)
      implicit none
      real*8 ,intent(in):: s
      real*8 n,df,sh,sh1,s1,ps
      integer*4 no
      if(s == 0.d0)then
        f1=zt0
        return
      endif
        s1=1.d0-s
        sh=.5d0*s
        sh1=.5d0-sh
        f1=-m_pi**sh/s/s1
        ps=m_pi**(s-.5d0)
        n=1.d0
        no=nolog*4
        do
          df=n**(-s)*gamma2(sh,m_pi*n**2)
     $         +ps*n**(-s1)*gamma2(sh1,m_pi*n**2)
          f1=f1+df
          no=no+nogam2*2
c          write(*,'(a,1p10g12.4)')'zt ',n,s,df,f1*gammai(sh)
          if((df/f1)**2 <= no*epso)then
            f1=f1*gammai(sh)
            exit
          endif
          n=n+1.d0
        enddo
      return
      end function

      real*8 recursive function dzeta(x) result(f1)
      implicit none
      real*8 ,intent(in):: x
      real*8 u,df,x1,x2,d,xd
      real*8 n
      integer*4 k,no
      if(x == 0.d0)then
        f1=dz0
        return
      endif
      x1=1.d0-x
      if(abs(x1) < ztlim)then
        u=-1.d0
        f1=-1.d0/x1**2
        no=0
        do k=1,nstg
          if(k .gt. nstg)then
            df=stgn(k)*u
            no=no+40
          else
            df=stg(k)*u
          endif
          f1=f1+df
          no=no+1
          if(abs(df/f1)**2 <= no*epso)then
            exit
          endif
          u=u*x1/k
        enddo
      elseif(abs(x+0.5d0) > cth)then
        f1=0.d0
      elseif(x < 0.d0)then
        f1=m_2pi**(-x1)*gamma(x1)*
     $       ((m_pi*cosp(.5d0*x) +
     $       2.d0*(polygamma(x1) + log(m_2_pi))*
     $       sinp(.5d0*x))*zeta(x1)
     $       - 2.d0*sinp(.5d0*x)*dzeta(x1))
      else
        f1=0.d0
        n=1.d0
        x2=2.d0**(-x)
        d=1.d0/(1.d0-x2)
        xd=x2*log(2.d0)*d
        no=nolog*2+1
        do
          df=-n**(-x)*(xd+log(n))*d
          f1=f1+df
          no=no+nolog*2+4
          if(abs(df/f1)**2 <= no*epso)then
            exit
          endif
          n=n+2.d0
        enddo
      endif
      return
      end function

      complex*16 recursive  function dczeta(x) result(f1)
      implicit none
      complex*16 ,intent(in):: x
      complex*16 u,df,x1,x2,d,xd
      real*8 n
      integer*4 k,no
      x1=zeroim(1.d0-x)
      if(abs(x1) < ztlim)then
        u=-1.d0
        f1=-1.d0/x1**2
        no=nolog
        do k=1,nstg
          if(k .gt. nstg)then
            df=stgn(k)*u
            no=no+40
          else
            df=stg(k)*u
          endif
          f1=f1+df
          no=no+1
          if(abs((abs(dble(df))+abs(imag(df)))/f1)**2
     $         <= no*epso)then
            exit
          endif
          u=u*x1/k
        enddo
      elseif(abs(dble(x)+0.5d0) > cth)then
        f1=0.d0
      elseif(dble(x) < 0.d0)then
c        f1=m_2pi**(-x1)*cgamma(x1)*
        f1=exp(-x1*log(m_2pi)+cloggamma(x1))*
     $       ((m_pi*ccosp(.5d0*x) +
     $       2.d0*(cpolygamma(x1) + log(m_2_pi))*
     $       csinp(.5d0*x))*czeta(x1)
     $       - 2.d0*csinp(.5d0*x)*dczeta(x1))
      else
        f1=0.d0
        n=1.d0
        x2=2.d0**(-x)
        d=1.d0/(1.d0-x2)
        xd=x2*log(2.d0)*d
        no=nolog*2+4
        do
          df=-n**(-x)*(xd+log(n))*d
          f1=f1+df
          no=no+nolog*2+4
          if(abs((abs(dble(df))+abs(imag(df)))/f1)**2
     $         <= no*epso)then
            exit
          endif
          n=n+2.d0
        enddo
      endif
      return
      end function

      complex*16 recursive function czeta(s) result(f1)
      implicit none
      complex*16 ,intent(in):: s
      complex*16 df,sh,sh1,s1,ps
      real*8 n
c      real*8 ,parameter :: simth=(1.d0+sqrt(2.d0))*.5d0*8.d0
      integer*4 no
      if(s == czero)then
        f1=(-0.5d0,0.d0)
      elseif(abs(imag(s)) > simth)then
        f1=czetab(s)
      else
        s1=1.d0-s
        sh=.5d0*s
        sh1=.5d0*s1
        f1=-m_pi**sh/s/s1
        ps=m_pi**(s-.5d0)
        n=1.d0
        no=nolog*4
        do
          df=n**(-s)*cgamma2r(sh,m_pi*n**2)
     $         +ps*n**(-s1)*cgamma2r(sh1,m_pi*n**2)
          f1=f1+df
          no=no+(nogam2+nolog)*2
          if(abs((abs(dble(df))+abs(imag(df)))/f1)**2
     $         <= no*epso)then
            f1=f1*cgammai(sh)
            return
          endif
          n=n+1.d0
        enddo
      endif
      end function

      complex*16 recursive function czetab(s) result(f1)
      implicit none
      complex*16 ,intent(in):: s
      complex*16 u,df,ns
      real*8 k,k1,adf,adf0,n
      integer*4 no,i
      if(dble(s) < 0.5d0)then
        f1=cgamma(1.d0-s)*m_2pi**s/m_pi*csinp(.5d0*s)*czetab(1.d0-s)
      else
        n=ceiling(abs(s))+1.d0
        f1=0.d0
        do i=1,int(n)
          f1=f1+1.d0/dble(i)**s
        enddo
        ns=n/n**s
        u=s*ns/n**2
        f1=f1+ns*(1.d0/(s-1.d0)-.5d0/n)+u*berbf(2)
        k=2.d0
        adf0=veryl
        no=(int(n)+1)*nolog+10
        do
          k1=k+2.d0
          u=u*(s+k)*(s+k-1.d0)/n**2
          df=u*bernbf(int(k1))
          adf=abs(dble(df))+abs(imag(df))
          if(adf > adf0 .and. k > 2.d0*n)then
c            write(*,'(a,1p10g12.4)')'czb-n ',k,n,u,df,f1
            exit
          endif
          adf0=adf
          f1=f1+df
          no=no+14
          if(adf**2/abs(f1)**2 <= no*epso)then
c            write(*,'(a,1p10g12.4)')'czb-c ',k,n,u,df,f1
            exit
          endif
          k=k1
        enddo
      endif
      return
      end function

      complex*16 recursive function czetas(s) result(f1)
      implicit none
      complex*16 ,intent(in):: s
      complex*16 u,df,s1
      real*8 n,adf,adf0
      integer*4 no
      s1=1.d0-s
      f1=-1.d0/s1+stg(0)
      u=cone
      n=1.d0
      no=6
      adf0=veryl
      do
        u=u*s1/n
        if(n <= dble(nstg))then
          df=u*stg(int(n))
          no=no+6
        else
          df=u*stgn(int(n))
          no=no+40
        endif
        adf=abs(dble(df))+abs(imag(df))
        if(adf > adf0 .and. n > 40.d0)then
          exit
        endif
        adf0=adf
        f1=f1+df
        if(adf**2/abs(f1)**2 <= no*epso)then
          exit
        endif
        n=n+1.d0
      enddo
      return
      end function

      real*8 recursive function zeta2(s,a) result(f1)
      implicit none
      real*8 ,intent(in):: s,a
      real*8 n,ak,n1
      integer*4 i,nt
      if(s == 0.d0)then
        f1=.5d0-a
        return
      elseif(s == 1.d0)then
        f1=1.d0/0.d0
        return
      elseif(anint(s) == s .and. s < 0.d0 .and. a >= 0.d0)then
        f1=berpolr(1-nint(s),a)/(s-1.d0)
        return
      elseif(a == 0.d0 .or. a == 1.d0)then
        f1=zeta(s)
        return
      endif
      n=floor(dble(a))
      nt=int(n)
      n1=dble(nt)
      if(nt > 0)then
        if(a /= n1)then
          f1=zeta2(s,a-n1)-abs(a-n1)**(-s)
        else
          f1=zeta(s)
        endif
        ak=a-n1+1.d0
        do i=1,nt-1
          if(ak /= 0.d0)then
            f1=f1-abs(ak)**(-s)
          endif
          ak=ak+1.d0
        enddo
      elseif(nt < 0)then
        f1=zeta2(s,a-n1)+abs(a)**(-s)
        ak=a+1.d0
        do i=1,-nt-1
          if(ak /= 0.d0)then
            f1=f1+abs(ak)**(-s)
          endif
          ak=ak+1.d0
        enddo
      else
        f1=zconv(s,a)
c        write(*,'(a,1p10g12.4)')'z2 ',s,a,f1
      endif
      return
      end function

      complex*16 recursive function czeta2(s,a) result(f1)
      implicit none
      complex*16 ,intent(in):: s,a
      integer*4 ,parameter :: nmax=1000
      complex*16 ak
      real*8 n,n1
      integer*4 i,nt
      if(s == czero)then
        f1=.5d0-a
c        write(*,'(a,1p10g12.4)')'cz2-0 ',s,a,f1
        return
      elseif(s == cone)then
        f1=1.d0/czero
c        write(*,'(a,1p10g12.4)')'cz2-1 ',s,a,f1
        return
      elseif(imag(s) == 0.d0 .and. imag(a) == 0.d0)then
        f1=dcmplx(zeta2(dble(s),dble(a)),0.d0)
c        write(*,'(a,1p10g12.4)')'cz2-r ',s,a,f1
        return
      elseif(imag(s) == 0.d0 .and. anint(dble(s)) == dble(s)
     $       .and. dble(s) < 0.d0 .and. dble(a) >= 0.d0)then
        f1=berpol(1-nint(dble(s)),a)/(dble(s)-1.d0)
c        write(*,'(a,1p10g12.4)')'cz2-bn ',s,a,f1
        return
      elseif(a == cone .or. a == czero)then
        f1=czeta(s)
c        write(*,'(a,1p10g12.4)')'cz2-01 ',s,a,f1
        return
      endif
      n=ceiling(dble(a))-1.d0
      nt=int(n)
      n1=dble(nt)
      if(nt > 0)then
        f1=czeta2(s,a-n1)
        ak=zeroim(a-n1)
        do i=1,nt
          if(ak /= czero)then
            f1=f1-(ak*ak)**(-s*.5d0)
          endif
          ak=ak+1.d0
        enddo
c        write(*,'(a,1p10g12.4)')'cz2-nt>0 ',dble(nt),s,a,f1
      elseif(nt < 0)then
        f1=czeta2(s,a-n1)
        ak=zeroim(a)
        do i=1,-nt
          if(ak /= czero)then
            f1=f1+(ak*ak)**(-s*.5d0)
          endif
          ak=ak+1.d0
        enddo
c        write(*,'(a,1p10g12.4)')'cz2-nt<0 ',dble(nt),s,a,f1
c      elseif(imag(a) < 0.d0)then
c        f1=conjz(czeta2(conjz(s),conjz(a)))
c        write(*,'(a,1p10g12.4)')'cz2-cj ',s,a,f1
      elseif(abs(imag(a)) > zimth)then
        f1=(czeta2(s,.5d0*a)+czeta2(s,.5d0*a+.5d0))*2.d0**(-s)
c        write(*,'(a,1p10g12.4)')'cz2-2a ',s,a,f1
c        write(*,'(a,1p10g12.4)')'cz2-vz ',s,a,f1
      elseif(dble(a) >= .5d0)then
        f1=czconvb(s,a)
      else
        f1=czconv(s,zeroim(a))
      endif
      return
      end function

      complex*16 function czconvb(s,a) result(f1)
      implicit none
      complex*16 ,intent(in):: s,a
      complex*16 u,df
      real*8 k,k1
      integer*4 no
      real*8 ,parameter :: lm=m_pi
      f1=cgamma2(s,lm*a)*a**(-s)
      k=1.d0
      no=nocgam
      do
        df=cgamma2(s,lm*(a+k))*zeroim(a+k)**(-s)
        no=no+nocgam+nolog
        f1=f1+df
        if(((abs(dble(df))+abs(imag(df)))/abs(f1))**2 <= no*epso)then
          exit
        endif
        k=k+1.d0
      enddo
      u=lm**(s-1.d0)
      f1=f1+u*berpol(0,a)/(s-1.d0)
      k=0.d0
      no=no+30
      do
        k1=k+1.d0
        u=-u*lm/k1
        df=u*berpol(int(k1),a)/(s+k)
        f1=f1+df
c        write(*,'(a,1p10g12.4)')'czcb-k ',k,s,a,df,f1
        no=no+100
        if(((abs(dble(df))+abs(imag(df)))/abs(f1))**2 <= no*epso)then
          f1=f1*cgammai(s)
          exit
        endif
        k=k1
      enddo
      return
      end function

      complex*16 function czconv(s,a) result(f1)
      implicit none
      complex*16 ,intent(in):: s,a
      complex*16 df,sh,sh1,sh2,sh3,a2,an,an1,an2,an3,cgi,cgi1,c
      real*8 n,adf,adf0
      integer*4 no
      real*8 ,parameter :: lm=m_pi,lm1=m_pi**2/lm
      a2=a**2
      sh=.5d0*s
      sh1=sh+.5d0
      sh2=zeroim(.5d0-sh)
      cgi=cgammai(sh)
      cgi1=cgammai(sh1)
      f1=(cgamma2(sh,lm*a2)*cgi+cgamma2(sh1,lm*a2)*cgi1)
     $     /a2**sh+2.d0*m_pi**sh/(s-1.d0)*cgi
      if(dble(a) <= 0.d0 .or. dble(a) > 1.d0)then
        write(*,'(a,1p10g12.4)')'czc-00 ',a,sh,sh1,lm*a2,f1
      endif
      n=1.d0
      no=(nogam2+nolog)*2
      adf0=veryl
      do
        an=a+n
        if(an /= czero)then
          an2=an**2
          df=(cgamma2(sh,lm*an2)*cgi+cgamma2(sh1,lm*an2)*cgi1)
     $     /an2**sh
          no=no+nogam2*2+nolog
c          write(*,'(a,1p10g12.4)')'czc-cg  ',n,a,df,f1
        else
          df=czero
        endif
        an1=a-n
        if(an1 /= czero)then
          an3=an1**2
          df=df+(cgamma2(sh,lm*an3)*cgi+cgamma2(sh1,lm*an3)*cgi1
     $         *merge(1.d0,-1.d0,dble(an1) >= 0.d0))/an3**sh
          no=no+nogam2*2+nolog
        endif
c        write(*,'(a,1p10g12.4)')'czc-cg1 ',n,a,df,f1
        adf=abs(dble(df))+abs(imag(df))
        if(adf > adf0)then
c          write(*,'(a,1p5g21.10)')'czc-nc ',n,a,f1,df
          exit
        endif
        adf0=adf
        f1=f1+df
        if((adf/abs(f1))**2 <= no*epso)then
c          write(*,'(a,1p10g12.4)')'czc-cv ',n,a,f1,df
          exit
        endif
        n=n+1.d0
      enddo
      f1=f1*.5d0
      sh3=1.d0-sh
      n=1.d0
      c=m_pi**(s-0.5d0)
      adf0=veryl
      do
c        write(*,'(a,1p10g12.4)')'czc-10 ',n,sh2,sh3,lm1*n**2
        df=c*(cgamma2(sh2,dcmplx(lm1*n**2,0.d0))*cgi*ccosp(2.d0*n*a)
     $       +cgamma2(sh3,dcmplx(lm1*n**2,0.d0))*cgi1*csinp(2.d0*n*a))
     $       *n**(s-1.d0)
        adf=abs(dble(df))+abs(imag(df))
        if(adf > adf0)then
          exit
        endif
        adf0=adf
        no=no+nogam2*2+nolog+30
        f1=f1+df
c        write(*,'(a,1p10g12.4)')'czc1 ',n,s,a,df
c        write(*,'(1p5g21.10)')f1
        if((adf/abs(f1))**2 <= no*epso)then
          exit
        endif
        n=n+1.d0
      enddo
      return
      end function

      real*8 function zconv(s,a) result(f1)
      implicit none
      real*8 ,intent(in):: s,a
      real*8 df,sh,sh1,sh2,sh3,a2,an,an2,cgi,cgi1,c
      real*8 n,adf,adf0
      integer*4 no
      real*8 ,parameter :: lm=m_pi,lm1=m_pi**2/lm
      a2=a**2
      sh=.5d0*s
      sh1=sh+.5d0
      sh2=.5d0-sh
      cgi=gammai(sh)
      cgi1=gammai(sh1)
      f1=(gamma2(.5d0*s,lm*a2)*cgi+gamma2(sh1,lm*a2)*cgi1)
     $     /a2**sh+2.d0*lm**sh/(s-1.d0)*cgi
      n=1.d0
      no=nogam2
      adf0=veryl
      do
        an=a+n
        an2=an**2
        df=(gamma2(sh,lm*an2)*cgi+gamma2(sh1,lm*an2)*cgi1)
     $     /an2**sh
        an=a-n
        an2=an**2
        df=df+(gamma2(sh,lm*an2)*cgi-gamma2(sh1,lm*an2)*cgi1)
     $       /an2**sh
        adf=abs(df)
        if(adf > adf0)then
          exit
        endif
        adf0=adf
        no=no+nogam2*2
        f1=f1+df
c        write(*,'(a,1p10g12.4)')'zc ',n,f1,df
        if(adf**2 <= no*f1**2*epso)then
          exit
        endif
        n=n+1.d0
      enddo
      f1=f1*.5d0
      sh3=1.d0-sh
      n=1.d0
      c=m_pi**(s-0.5d0)
      adf0=veryl
      do
        df=c*(gamma2(sh2,lm1*n**2)*cgi*cosp(2.d0*n*a)
     $       +gamma2(sh3,lm1*n**2)*cgi1*sinp(2.d0*n*a))
     $       *n**(s-1.d0)
        no=no+(nogam2+nolog)
        adf=abs(df)
        if(adf > adf0)then
          exit
        endif
        adf0=adf
        f1=f1+df
c        write(*,'(a,1p10g12.4)')'zc1 ',n,f1,df,gamma2(sh2,lm1*n**2),
c     $       gamma2(sh3,lm1*n**2)
        if(adf**2 <= no*f1**2*epso)then
          exit
        endif
        n=n+1.d0
      enddo
      return
      end function

      complex*16 recursive function chzeta2(s,a) result(f1)
      implicit none
      complex*16 ,intent(in):: s,a
      complex*16 eim,af,a1
      real*8 fan
      a1=zeroim(a)
      if(dble(a1) >= 0.d0)then
        if(imag(s) == 0.d0 .and. imag(a) == 0.d0)then
          f1=zeta2(dble(s),dble(a))
c          write(*,'(a,1p10g12.4)')'chz2-r ',s,a1,f1
        else
          f1=czeta2(s,a1)
        endif
      elseif(imag(a1) >= 0.d0)then
        fan=floor(-dble(a))
        eim=cexpp(zeroim(dcmplx(imag(s),-dble(s))))
        if(a1+fan == 0.d0)then
          f1=eim*czeta2(s,a1)+czeta(s)*(1.d0-eim)
        else
          af=a1+fan
          f1=eim*czeta2(s,a1)
     $         +(czeta2(s,a1+fan+1.d0)
     $         +(1.d0+floor(dble(a1))+fan)*zeroim(af*af)**(-s/2.d0))
     $         *(1.d0-eim)
        endif
c        write(*,'(a,1p10g12.4)')'chz2-p ',s,a,fan,eim,f1
      else
        f1=conjz(chzeta2(conjz(s),conjz(a)))
c        write(*,'(a,1p10g12.4)')'chz2-n ',s,a,f1
      endif
      return
      end

      complex*16 function dchzeta2(s,a) result(f1)
      implicit none
      complex*16 ,intent(in):: s,a
      complex*16 eim,deim,af2,gp,af1
      real*8 fan,g
      if(imag(s) == 0.d0 .and. imag(a) == 0.d0)then
        f1=dhzeta2(dble(s),dble(a))
      elseif(dble(a) >= 0.d0)then
        f1=dczeta2(s,a)
      elseif(imag(a) >= 0.d0)then
        fan=floor(-dble(a))
        eim=cexpp(dcmplx(imag(s),-dble(s)))
        deim=(0.d0,mm_pi)*eim
c        write(*,'(a,1p10g12.6)')'dchz2-0 ',s,a,a+fan
        if(a+fan == czero)then
          f1=deim*czeta2(s,a)+eim*dczeta2(s,a)
     $         +dczeta(s)*(1.d0-eim)
     $         -czeta(s)*deim
        else
          g=1.d0+floor(dble(a))+fan
          af1=a+fan+1.d0
          if(g == 0.d0)then
            f1=deim*czeta2(s,a)+eim*dczeta2(s,a)
     $           +dczeta2(s,af1)*(1.d0-eim)
     $           -czeta2(s,af1)*deim
          else
            af2=(a+fan)**2
            gp=g/af2**(s/2.d0)
            f1=deim*czeta2(s,a)+eim*dczeta2(s,a)
     $           +(dczeta2(s,af1)-gp*(.5d0*log(af2)))*(1.d0-eim)
     $           -(czeta2(s,af1)+gp)*deim
          endif
        endif
c        write(*,'(a,1p10g12.6)')'dchz2-posima ',s,a,f1
      else
        fan=floor(-dble(a))
        af1=a+fan+1.d0
        eim=cexpp(dcmplx(imag(s),-dble(s)))
        deim=(0.d0,mm_pi)*eim
        f1=deim*czeta2(s,a)+eim*dczeta2(s,a)
     $       +dczeta2(s,af1)*(1.d0-eim)-czeta2(s,af1)*deim
c        write(*,'(a,1p10g12.6)')'dchz2-negima ',s,a,f1
      endif
      return
      end

      complex*16 function dhzeta2(s,a) result(f1)
      implicit none
      real*8 ,intent(in):: s,a
      complex*16 eim,deim
      real*8 fan,g,gp,af2
      if(a >= 0.d0)then
        f1=dcmplx(dzeta2(s,a),0.d0)
      else
        fan=floor(-a)
        eim=dcmplx(cosp(s),sinp(s))
        deim=(0.d0,mm_pi)*eim
c        write(*,'(a,1p10g12.6)')'dhz2-0 ',s,a,a+fan
        if(a+fan == 0.d0)then
          f1=deim*zeta2(s,a)+eim*dzeta2(s,a)
     $         +dzeta(s)*(1.d0-eim)
     $         -zeta(s)*deim
        else
          g=1.d0+floor(a)+fan
          if(g == 0.d0)then
            f1=deim*zeta2(s,a)+eim*dzeta2(s,a)
     $           +dzeta2(s,a+fan+1.d0)*(1.d0-eim)
     $           -zeta2(s,a+fan+1.d0)*deim
          else
            af2=(a+fan)**2
            gp=g/af2**(s/2.d0)
            f1=deim*zeta2(s,a)+eim*dzeta2(s,a)
     $           +(dzeta2(s,a+fan+1.d0)
     $           -gp*(.5d0*log(af2)))*(1.d0-eim)
     $           -(zeta2(s,a+fan+1.d0)+gp)*deim
          endif
        endif
      endif
      return
      end

      complex*16 recursive function dczeta2(s,a) result(f1)
      implicit none
      complex*16 ,intent(in):: s,a
      complex*16 u,f,a1,u1,du,du1,sn
      real*8 n
      if(imag(s) == 0.d0 .and. imag(a) == 0.d0)then
        f1=dcmplx(dzeta2(dble(s),dble(a)),0.d0)
        return
      endif
      if(dble(a) < ath1)then
        f1=dczeta2(s,a+1.d0)-0.5d0*(a**2)**(-s/2.d0)*log(a**2)
        return
      endif
      a1=1.d0-a
      if(abs(a1) >= ath2)then
        f1=(dczeta2(s,.5d0*a)+dczeta2(s,.5d0*a+.5d0)
     $       -(czeta2(s,.5d0*a)+czeta2(s,.5d0*a+.5d0))*log(2.d0)
     $       )/2.d0**s
        return
      endif
      u=1.d0
      du=0.d0
      f=dczeta(s)
      n=1.d0
      do
        sn=s+n
        if(sn == cone)then
          u1=czero
          du1=u/n*a1
          u=czero
          du=du1/(n+1.d0)*a1
c          f1=f+m_euler*du1+zeta(2.d0)*du
          f1=f+a1**n*(-dzsc(n)+m_euler/n)+zt2*du
        elseif(sn == czero)then
          u1=-u/n*a1
          du1=(u-du)/n*a1
          u=czero
          du=u1/(n+1.d0)*a1
c          f1=f+dz0*u1+zt0*du1+m_euler*du
          f1=f-0.5d0*log(m_2pi)*u1-0.5d0*du1
     $         -a1**(n+1.d0)*(-dzsc(n+1.d0)+m_euler/(n+1.d0))
        else
          u1=u*(sn-1.d0)/n*a1
          du1=(u+du*(sn-1.d0))/n*a1
          u=u1*sn/(n+1.d0)*a1
          du=(u1+du1*sn)/(n+1.d0)*a1
          f1=f+dczeta(sn)*u1+dczeta(sn+1.d0)*u
     $         +czeta(sn)*du1+czeta(sn+1.d0)*du
        endif
        if(f1 == f)then
c          write(*,'(a,1p10g12.4)')'dcz2 ',s,a,sn,f1
          exit
        endif
        f=f1
        n=n+2.d0
      enddo
      return
      end function

      real*8  function dzsc(n) result(f1)
      implicit none
      real*8 ,intent(in):: n
      real*8 s,fk1
      integer*4 i
      if(n <= 1.d0)then
        f1=0.d0
      else
        s=0.d0
        fk1=1.d0
        do i=1,int(n)-1
          s=s*i+fk1
          fk1=fk1*i
        enddo
        fk1=fk1*n
        f1=s/fk1
      endif
c      write(*,'(a,1p8g15.7)')'dzsc ',n,s,fk1,f1
      return
      end function

      real*8 recursive function dzeta2(s,a) result(f1)
      use tfstk, only:ktfenanq
      implicit none
      real*8 ,intent(in):: s,a
      real*8 u,f,a1,u1,du,du1,sn
      real*8 n
      if(a == 0.d0)then
        f1=dzeta(s)
        return
      elseif(a < ath1)then
        f1=dzeta2(s,a+1.d0)-0.5d0*(a**2)**(-s/2.d0)*log(a**2)
        return
      endif
      a1=1.d0-a
      if(abs(a1) >= ath2)then
        f1=(dzeta2(s,.5d0*a)+dzeta2(s,.5d0*a+.5d0)
     $       -(zeta2(s,.5d0*a)+zeta2(s,.5d0*a+.5d0))*log(2.d0)
     $       )/2.d0**s
        return
      endif
      u=1.d0
      du=0.d0
      f=dzeta(s)
      n=1.d0
      if(s == anint(s) .and. s <= 0)then
        do
          sn=s+n
          if(sn == 1.d0)then
            u1=0.d0
            du1=u/n*a1
            u=0.d0
            du=du1/(n+1.d0)*a1
            f1=f+a1**n*(-dzsc(n)+m_euler/n)+zt2*du
          elseif(sn == 0.d0)then
            u1=-u/n*a1
            du1=(u-du)/n*a1
            u=0.d0
            du=u1/(n+1.d0)*a1
            f1=f-0.5d0*log(m_2pi)*u1-0.5d0*du1
     $           -a1**(n+1.d0)*(-dzsc(n+1.d0)+m_euler/(n+1.d0))
          else
            u1=u*(sn-1.d0)/n*a1
            du1=(u+du*(sn-1.d0))/n*a1
            u=u1*sn/(n+1.d0)*a1
            du=(u1+du1*sn)/(n+1.d0)*a1
            f1=f+dzeta(sn)*u1+dzeta(sn+1.d0)*u
     $           +zeta(sn)*du1+zeta(sn+1.d0)*du
          endif
          if(f1 == f)then
            exit
          endif
          f=f1
          n=n+2.d0
        enddo
      else
        do
          sn=s+n
          u1=u*(sn-1.d0)/n*a1
          du1=(u+du*(sn-1.d0))/n*a1
          u=u1*sn/(n+1.d0)*a1
          du=(u1+du1*sn)/(n+1.d0)*a1
          f1=f+dzeta(sn)*u1+dzeta(sn+1.d0)*u
     $         +zeta(sn)*du1+zeta(sn+1.d0)*du
c          write(*,'(a,1p10g12.4)')'dz2-n ',du,du1,u,u1,f1
          if(f1 == f)then
            exit
          endif
          f=f1
          n=n+2.d0
        enddo
      endif
      return
      end function

      complex*16  function zetads(s,a) result(f1)
      implicit none
      complex*16 ,intent(in):: s,a
      integer*4 ,parameter ::nmax=100
      complex*16 u,f,s1,akp(0:nmax)
      real*8 n
      integer*4 k,km
      n=0.d0
      s1=s-1.d0
      f=0.d0
      km=0
      akp(0)=a**(-s1)
      do
        u=1.d0
        f1=akp(0)
        do k=1,int(n)
          if(k > km)then
            akp(k)=(a+k)**(-s1)
            km=k
          endif
          u=-u/dble(k)*dble(n-k+1)
          f1=f1+u*akp(k)
        enddo
        f1=f+f1/(n+1.d0)
c        write(*,'(a,1p8g15.7)')'zetads ',n,f1/s1,akp(km)
        if(f == f1 .or. km >= nmax)then
          f1=f1/s1
          return
        endif
        f=f1
        n=n+1.d0
      enddo
      end function

      complex*16  function zetabk(z,a) result(f1)
      implicit none
      complex*16 ,intent(in):: z,a
      complex*16 u,f,z1
      integer*4 k
      z1=z-1.d0
      u=1.d0/z1
      f=u+.5d0/a
      do k=2,ibmax,2
        u=u*(z1+k-3)*(z1+k-2)/a**2
        f1=f+berbf(k)*u
c        write(*,'(a,i5,1p8g15.7)')'zetabk ',k,z1,berbf(k),u,f1
        if(f1 == f)then
          exit
        endif
        f=f1
      enddo
      f1=f1/a**z1
      end function

      real*8 pure recursive function polygamma(x1) result(f1)
      use tfstk, only:ktfenanq
      implicit none
      real*8 ,intent(in):: x1
      real*8 x
      if(x1 <= 0.d0 .and. anint(x1) == x1)then
        f1=-1.d0/0.d0
        return
      endif
      x=x1-1.d0
      if(x < -1.d0)then
        f1=-m_pi*cosp(x)/sinp(x)+polygamma(-x)
        f1=merge(1.d0/0.d0,f1,ktfenanq(f1))
      else
        f1=-1.d0+(.5d0+x)/(0.5d0+ggamma+x)-(
     $       gc1/(1.d0+x)**2   +gc2/(2.d0+x)**2
     $       +gc3/(3.d0+x)**2  +gc4/(4.d0+x)**2
     $       +gc5/(5.d0+x)**2  +gc6/(6.d0+x)**2
     $       +gc7/(7.d0+x)**2  +gc8/(8.d0+x)**2
     $       +gc9/(9.d0+x)**2  +gc10/(10.d0+x)**2
     $       +gc11/(11.d0+x)**2+gc12/(12.d0+x)**2
     $       +gc13/(13.d0+x)**2+gc14/(14.d0+x)**2
     $       )/(
     $       gc0+gc1/(1.d0+x)+gc2/(2.d0+x)
     $       +gc3/(3.d0+x)  +gc4/(4.d0+x)
     $       +gc5/(5.d0+x)  +gc6/(6.d0+x)
     $       +gc7/(7.d0+x)  +gc8/(8.d0+x)
     $       +gc9/(9.d0+x)  +gc10/(10.d0+x)
     $       +gc11/(11.d0+x)+gc12/(12.d0+x)
     $       +gc13/(13.d0+x)+gc14/(14.d0+x)
     $       )+log(0.5d0+ggamma+x)
c     $       +gc9/(9.d0+x)+gc10/(10.d0+x))
      endif
      return
      end function

      complex*16 pure recursive function cpolygamma(x1) result(f1)
      use tfstk, only:ktfenanq
      implicit none
      complex*16 ,intent(in):: x1
      complex*16 x
      if(imag(x1) == 0.d0)then
        f1=dcmplx(polygamma(dble(x1)),0.d0)
        return
      endif
      x=zeroim(x1-1.d0)
      if(dble(x) < -1.d0)then
        f1=-m_pi*ccosp(x)/csinp(x)+cpolygamma(-x)
        f1=merge(dcmplx(1.d0/0.d0,0.d0),f1,
     $       ktfenanq(dble(f1)) .or. ktfenanq(imag(f1)))
      else
        f1=-1.d0+(.5d0+x)/(0.5d0+ggamma+x)-(
     $       gc1/(1.d0+x)**2   +gc2/(2.d0+x)**2
     $       +gc3/(3.d0+x)**2  +gc4/(4.d0+x)**2
     $       +gc5/(5.d0+x)**2  +gc6/(6.d0+x)**2
     $       +gc7/(7.d0+x)**2  +gc8/(8.d0+x)**2
     $       +gc9/(9.d0+x)**2  +gc10/(10.d0+x)**2
     $       +gc11/(11.d0+x)**2+gc12/(12.d0+x)**2
     $       +gc13/(13.d0+x)**2+gc14/(14.d0+x)**2
     $       )/(
     $       gc0+gc1/(1.d0+x)+gc2/(2.d0+x)
     $       +gc3/(3.d0+x)  +gc4/(4.d0+x)
     $       +gc5/(5.d0+x)  +gc6/(6.d0+x)
     $       +gc7/(7.d0+x)  +gc8/(8.d0+x)
     $       +gc9/(9.d0+x)  +gc10/(10.d0+x)
     $       +gc11/(11.d0+x)+gc12/(12.d0+x)
     $       +gc13/(13.d0+x)+gc14/(14.d0+x)
     $       )+log(0.5d0+ggamma+x)
      endif
      return
      end function

      complex*16 function chpolygamma2(n,x) result(f1)
      implicit none
      complex*16 ,intent(in):: n,x
      if(imag(n) == 0.d0)then
        if(dble(n) == 0.d0)then
          f1=cpolygamma(x)
          return
        elseif(dble(n) == anint(dble(n))
     $         .and. dble(n) >= 0.d0)then
          if(x == cone)then
            f1=-(-1.d0)**nint(dble(n))
     $           *gamma(dble(n)+1.d0)*czeta(n+1.d0)
          else
            f1=-(-1.d0)**nint(dble(n))
     $           *gamma(dble(n)+1.d0)*chzeta2(n+1.d0,x)
          endif
          return
        endif
      endif
      if(x == cone)then
        f1=(dczeta(n+1.d0)+(m_euler+cpolygamma(-n))*czeta(n+1.d0))
     $       *cgammai(-n)
      else
        f1=(dchzeta2(n+1.d0,x)
     $       +(m_euler+cpolygamma(-n))*chzeta2(n+1.d0,x))
     $       *cgammai(-n)
      endif
      return
      end function

      complex*16 function cgpolygamma2(n,x) result(f1)
      implicit none
      complex*16 ,intent(in):: n,x
      if(imag(n) == 0.d0)then
        if(dble(n) == 0.d0)then
          f1=cpolygamma(x)
          return
        elseif(imag(x) == 0.d0)then
          f1=dcmplx(gpolygamma2(dble(n),dble(x)),0.d0)
          return
        elseif(dble(n) == anint(dble(n))
     $         .and. dble(n) >= 0.d0)then
          if(x == cone)then
            f1=-(-1.d0)**nint(dble(n))
     $           *gamma(dble(n)+1.d0)*czeta(n+1.d0)
          else
            f1=-(-1.d0)**nint(dble(n))
     $           *gamma(dble(n)+1.d0)*czeta2(n+1.d0,x)
          endif
          return
        endif
      endif
      if(x == cone)then
        f1=(dczeta(n+1.d0)+(m_euler+cpolygamma(-n))*czeta(n+1.d0))
     $       *cgammai(-n)
      else
        f1=(dczeta2(n+1.d0,x)+(m_euler+cpolygamma(-n))*czeta2(n+1.d0,x))
     $       *cgammai(-n)
      endif
      return
      end function

      real*8 function gpolygamma2(n,x) result(f1)
      implicit none
      real*8 ,intent(in):: n,x
      if(n == 0.d0)then
        f1=polygamma(x)
      elseif(n == anint(n) .and. n > 0.d0)then
        if(x == 1.d0)then
          f1=-(-1.d0)**nint(n)*gamma(n+1.d0)*zeta(n+1.d0)
        else
          f1=-(-1.d0)**nint(n)*gamma(n+1.d0)*zeta2(n+1.d0,x)
        endif
      elseif(x == 1.d0)then
        f1=(dzeta(n+1.d0)+(m_euler+polygamma(-n))*zeta(n+1.d0))
     $       *gammai(-n)
      else
c        write(*,'(a,1p10g12.4)')'gp2 ',n,x,dzeta2(n+1.d0,x),
c     $       polygamma(-n),zeta2(n+1.d0,x)
        f1=(dzeta2(n+1.d0,x)+(m_euler+polygamma(-n))*zeta2(n+1.d0,x))
     $       *gammai(-n)
      endif
      return
      end function

      real*8 pure function ferf(x)
      implicit none
      real*8 , intent(in)::x
c      ferf=sign(gammap(.5d0,x**2),x)
      ferf=erf(x)
      return
      end

      real*8 pure function ferfc(x)
      implicit none
      real*8 ,intent(in):: x
      ferfc=erfc(x)
c      real*8 x,gammaq,gammap
c      if(x < 0.d0)then
c        erfc=1.d0+gammap(.5d0,x**2)
c      else
c        erfc=gammaq(.5d0,x**2)
c      endif
      return
      end

      real*8 pure function productlog(x)
      implicit none
c     Including m_e(Napier's constant: Exp[1])
      real*8 ,intent(in):: x
      real*8 w,f1,d
      real*8 ,parameter ::eps=3.d-16,en=m_e
      if(x <= -1.d0/en)then
        productlog=-1.d0
        return
      endif
      if(x == 0.d0)then
        productlog=0.d0
        return
      endif
      if(x < -0.25d0)then
        w=sqrt1(2.d0*en*x+1.d0)
c        w=sqrt(2.d0*(en*x+1.d0))-1.d0
      elseif(x < 2.d0)then
        w=.5d0*sqrt1(4.d0*x)
c        w=.5d0*(sqrt(4.d0*x+1.d0)-1.d0)
      else
        w=log(x/log(x/log(x)))
      endif
      f1=x*exp(-w)
      d=(f1-w)/(f1+1.d0)
      do while(abs(d) > eps*abs(w))
        w=w+d
        f1=x*exp(-w)
        d=(f1-w)/(f1+1.d0)
      enddo
      productlog=w+d
      return
      end

      complex*16 pure function cproductlog(z)
      implicit none
c     Including m_e(Napier's constant: Exp[1])
      complex*16 ,intent(in):: z
      complex*16 w,f1,d
      real*8 ,parameter ::eps=3.d-16,en=m_e
      if(imag(z) == 0.d0)then
        if(dble(z) >= -1.d0/en)then
          cproductlog=productlog(dble(z))
          return
        endif
      endif
      if(abs(z+1.d0/en) < 0.2d0)then
        w=sqrt(2.d0*(en*z+1.d0))-1.d0
      elseif(abs(z) < 3.d0)then
        w=.5d0*(sqrt(4.d0*z+1.d0)-1.d0)
      else
        w=log(z/log(z/log(z)))
      endif
      f1=z*exp(-w)
      d=(f1-w)/(f1+1.d0)
      do while(abs(d) > abs(w)*eps)
        w=w+d
        f1=z*exp(-w)
        d=(f1-w)/(f1+1.d0)
      enddo
      cproductlog=w+d
      return
      end

      complex*16 function cpolylogs(s,z,f0,n) result(f)
      implicit none
      complex*16 ,intent(in):: s,z,f0
      integer*4 ,intent(in):: n
      complex*16 v,df
      integer*4 i,no
      f=f0
      v=z
      no=0
      do i=2,n
        v=v*z
        df=v/dble(i)**s
        f=f+df
        no=no+nolog+4
        if(abs((abs(dble(df))+abs(imag(df)))/f)**2
     $       <= no*epso)then
          exit
        endif
      enddo
c      write(*,'(a,1p10g12.4)')'cpls ',dble(i),s,z,v,df,f
      return
      end function 

      complex*16 recursive function cpolylog(s,z) result(f)
      implicit none
      complex*16 ,intent(in):: s,z
      complex*16 lz,s1,es
      real*8 la
      real*8 ,parameter :: lm=m_pi,plzth=0.5d0,epspl=log(epso*nopl)/2.d0
      if(z == czero)then
        f=czero
        return
      elseif(z == cone)then
        if(dble(s) > 1.d0)then
          f=czeta(s)
        elseif(dble(s) == 0.d0)then
          f=0.d0/0.d0
        else
          f=1.d0/0.d0
        endif
        return
      endif
      la=log(abs(z))
      if(11.d0*max(0.d0,la) < epspl+dble(s)*log(12.d0))then
        f=cpolylogs(s,z,z,11)
c      elseif(abs(z) <= .5d0)then
c        f=cpolylogs(s,z,z,10000)
c      elseif(abs(z) > 2.d0)then
c        f=(1.d0,m_2pi)**s*cgammai(s)
c     $       *chzeta2(1.d0-s,.5d0-log(zeroim(-1.d0/z))/(0.d0,m_2pi))
c     $       -(-1.d0)**s*cpolylog(s,1.d0/z)
      elseif(dble(s) >= 0.d0 .and. s == anint(dble(s)))then
        f=cpolylog1i(int(dble(s)),z)
      elseif(abs(imag(s)) > plith)then
        f=cpolylogl(s,z)
      else
        es=cexpp(dcmplx(imag(s),-dble(s)))
        s1=1.d0-s
        if(imag(z) > 0.d0 .or. dble(z) <= -1.d0
     $       .or. imag(z) == 0.d0 .and.  abs(z) < 1.d0)then
          lz=log(z)/(0.d0,m_2pi)
        else
          lz=0.5d0+log(zeroim(-z))/(0.d0,m_2pi)
        endif
        f=exp((0.d0,m_pi_2)*(s+1.d0)-s1*log(m_2pi)+cloggamma(s1))
     $       *(es*chzeta2(s1,lz)-chzeta2(s1,1.d0-lz))
c        f=exp((0.d0,m_pi_2)*(s+1.d0)-s1*log(m_2pi)+cloggamma(s1))
c     $       *cplzconv(s1,lz)
      endif
      return
      end function

      complex*16 recursive function cplzconv(s,a) result(f1)
      implicit none
      complex*16 ,intent(in):: s,a
      complex*16 df,sh,sh1,sh2,sh3,a2,an,an1,an2,an3,
     $     cgi,cgi1,c,a1,a21,es
      real*8 n,adf,adf0
      integer*4 no
      real*8 ,parameter :: lm=m_pi,lm1=m_pi**2/lm
      es=cexpp(dcmplx(-imag(s),dble(s)-1.d0))
c     es=cexpp(dcmplx(imag(s),-dble(s)))
      if(imag(a) > zimth .or.
     $     dble(a) <= 0.d0 .or. dble(a) >= 1.d0)then
        f1=es*chzeta2(s,a)-chzeta2(s,1.d0-a)
        write(*,'(a,1p10g12.4)')'cplz-a ',s,a,es,f1
      else          
        a1=1.d0-a
        a2=a**2
        a21=a1**2
        sh=.5d0*s
        sh1=sh+.5d0
        sh2=zeroim(.5d0-sh)
        cgi=cgammai(sh)
        cgi1=cgammai(sh1)
        f1=es*(cgamma2(sh,lm*a2)*cgi+cgamma2(sh1,lm*a2)*cgi1)
     $       /a2**sh
     $       -(cgamma2(sh,lm*a21)*cgi+cgamma2(sh1,lm*a21)*cgi1)
     $       /a21**sh
     $       +2.d0*(es-1.d0)*m_pi**sh/(s-1.d0)*cgi
        n=1.d0
        no=(nogam2+nolog)*4
        adf0=veryl
        do
          an=a+n
          if(an /= czero)then
            an2=an**2
            df=es*(cgamma2(sh,lm*an2)*cgi+cgamma2(sh1,lm*an2)*cgi1)
     $           /an2**sh
            no=no+nogam2*2+nolog
          else
            df=czero
          endif
          an=a1+n
          if(an /= czero)then
            an2=an**2
            df=df-(cgamma2(sh,lm*an2)*cgi+cgamma2(sh1,lm*an2)*cgi1)
     $           /an2**sh
            no=no+nogam2*2+nolog
          endif
          an1=a-n
          if(an1 /= czero)then
            an3=an1**2
c     write(*,'(a,1p10g12.4)')'czc-cg ',n,sh,an3,cgi,
c     $         cgamma2(sh,lm*an3)*cgi
            df=df+es*(cgamma2(sh,lm*an3)*cgi-cgamma2(sh1,lm*an3)*cgi1)
     $           /an3**sh
            no=no+nogam2*2+nolog
          endif
          an1=a1-n
          if(an1 /= czero)then
            an3=an1**2
c     write(*,'(a,1p10g12.4)')'czc-cg ',n,sh,an3,cgi,
c     $         cgamma2(sh,lm*an3)*cgi
            df=df-(cgamma2(sh,lm*an3)*cgi-cgamma2(sh1,lm*an3)*cgi1)
     $           /an3**sh
            no=no+nogam2*2+nolog
          endif
          adf=abs(dble(df))+abs(imag(df))
          if(adf > adf0)then
            exit
          endif
          adf0=adf
          f1=f1+df
          if((adf/abs(f1))**2 <= no*epso)then
c     write(*,'(a,1p5g21.10)')'czc-cv ',n,f1,df
            exit
          endif
          n=n+1.d0
        enddo
        f1=f1*.5d0
        sh3=1.d0-sh
        n=1.d0
        c=m_pi**(s-0.5d0)
        adf0=veryl
        do
c     write(*,'(a,1p10g12.4)')'czc-10 ',n,sh2,sh3,lm1*n**2
          df=c*(cgamma2(sh2,dcmplx(lm1*n**2,0.d0))
     $         *(es*ccosp(2.d0*n*a)-ccosp(2.d0*n*a1))*cgi
     $         +cgamma2(sh3,dcmplx(lm1*n**2,0.d0))
     $         *(es*csinp(2.d0*n*a)-csinp(2.d0*n*a1))*cgi1)
     $         *n**(s-1.d0)
          adf=abs(dble(df))+abs(imag(df))
          if(adf > adf0)then
            exit
          endif
          adf0=adf
          no=no+nogam2*2+nolog+30
          f1=f1+df
          if((adf/abs(f1))**2 <= no*epso)then
c     write(*,'(a,1p10g12.4)')'czc1 ',n,s,a,df
c     write(*,'(1p5g21.10)')f1
            exit
          endif
          n=n+1.d0
        enddo
      endif
      return
      end function

      complex*16 recursive function cpolyloga(s,z) result(f)
      implicit none
      complex*16 ,intent(in):: s,z
      complex*16 lz,u,df
      real*8 k,adf,adf0
      integer*4 no
      lz=log(zeroim(-z))
      u=-1.d0*lz**s*cgammai(s+1.d0)
      f=u*berbf(0)
      k=2.d0
      no=nolog+nogam
      adf0=veryl
      do
        u=-u*m_pi**2/lz**2*(s-k)*(s-k+1.d0)
        df=(1.d0-2.d0**(1.d0-k))*bernbf(int(k))*u
        adf=abs(dble(df))+abs(imag(df))
        if(adf > adf0)then
          exit
        endif
        adf0=adf
        no=no+nolog+16
        if((adf/abs(f))**2 <= no*epso)then
          exit
        endif
        k=k+2.d0
      enddo
      return
      end function

      complex*16 recursive function cpolylogl(s,z) result(f)
      implicit none
      complex*16 ,intent(in):: s,z
      complex*16 lz,v,df,s1,z1
      real*8 k
      integer*4 no
      lz=log(zeroim(z))
      if(abs(lz) > pllzth)then
        z1=sqrt(z)
        f=(cpolylog(s,z1)+cpolylog(s,-z1))/2.d0**s1
      else
        s1=1.d0-s
        f=cgamma(s1)/(-lz)**s1+czeta(s)
        no=nocgam+nolog+nozt
        v=cone
        k=1.d0
        do
          v=v*lz/k
          df=czeta(s-k)*v
          k=k+1.d0
          v=v*lz/k
          df=df+czeta(s-k)*v
          f=f+df
          no=no+nozt*2+14
          if(((abs(dble(df))+abs(imag(df)))/abs(f))**2
     $         <= no*epso)then
c            write(*,'(a,1p10g12.4)')'cpll ',k,s,z,df,f
            exit
          endif
          k=k+1.d0
        enddo
      endif
      return
      end function

      complex*16 recursive function cpolylog1i(s,z) result(f)
      implicit none
      complex*16 ,intent(in):: z
      integer*4 ,intent(in):: s
      complex*16 u,df,al,es,cs1
      real*8 n
      integer*4 no,m,s1
      real*8 ,parameter :: lm=m_pi
      s1=1-s
      cs1=dcmplx(dble(s1),0.d0)
      es=(0.d0,1.d0)**s1
      al=zeroim(0.5d0-log(zeroim(-z))/(0.d0,m_2pi))
      f=es*cgamma2(cs1,lm*(1.d0-al))/zeroim(1.d0-al)**s1
      n=2.d0
      no=nogam2+nolog*3
      do
        df=es*cgamma2(cs1,lm*(n-al))/zeroim(n-al)**s1
        f=f+df
        no=no+nogam2+nolog
        if((abs(dble(df))+abs(imag(df))/abs(f))**2 <= no*epso)then
          exit
        endif
        n=n+1.d0
      enddo
c      write(*,'(a,1p10g12.4)')'cpl-s0 ',s,z,f
      f=f+cgamma2(cs1,lm*al)/es/zeroim(al)**s1
c      write(*,'(a,1p10g12.4)')'cpl-s01 ',cgamma2(s1,lm*al),f
      n=1.d0
      no=nogam2+nolog
      do
        df=cgamma2(cs1,lm*(n+al))/es/zeroim(n+al)**s1
        f=f+df
        no=no+nogam2+nolog
        if(((abs(dble(df))+abs(imag(df)))/abs(f))**2 <= no*epso)then
          exit
        endif
        n=n+1.d0
      enddo
c      write(*,'(a,1p10g12.4)')'cpl-s1 ',dble(n),s,z,f,df
      u=-m_pi/lm**s
      if(s == 0)then
        f=f+u*berpol(0,al)
      else
        f=f+u*berpol(0,al)*sinp(.5d0*s)/(m_pi_2*s)
      endif
      no=50
      do m=1,ibmax,2
        u=-u*(0.d0,lm)/dble(m)
        if(s == m)then
          df=u*berpol(m,al)
        else
          df=u*berpol(m,al)
     $         *sinp(.5d0*dble(m-s))/(m_pi_2*dble(m-s))
        endif
        u=-u*(0.d0,lm)/dble(m+1)
        if(s == m+1)then
          df=df+u*berpol(m+1,al)
        else
          df=df+u*berpol(m+1,al)
     $         *sinp(.5d0*dble(m+s1))/(m_pi_2*(dble(m+s1)))
        endif
        f=f+df
        no=no+100
        if(((abs(dble(df))+abs(imag(df)))/abs(f))**2 <= no*epso)then
          exit
        endif
      enddo
      f=f/(m_2pi**s1)
      return
      end function

      complex*16 recursive function cpolylog1(s,z) result(f)
      implicit none
      complex*16 ,intent(in):: s,z
      complex*16 u,df,s1,al,es
      real*8 n
      integer*4 no,m
      real*8 ,parameter :: lm=m_pi
      s1=1.d0-s
      es=(0.d0,1.d0)**s1
      al=zeroim(0.5d0-log(zeroim(-z))/(0.d0,m_2pi))
      f=es*cgamma2(s1,lm*(1.d0-al))/zeroim(1.d0-al)**s1
      n=2.d0
      no=nogam2+nolog*3
      do
        df=es*cgamma2(s1,lm*(n-al))/zeroim(n-al)**s1
        f=f+df
        no=no+nogam2+nolog
        if((abs(dble(df))+abs(imag(df))/abs(f))**2 <= no*epso)then
          exit
        endif
        n=n+1.d0
      enddo
c      write(*,'(a,1p10g12.4)')'cpl-s0 ',s,z,f
      f=f+cgamma2(s1,lm*al)/es/zeroim(al)**s1
c      write(*,'(a,1p10g12.4)')'cpl-s01 ',cgamma2(s1,lm*al),f
      n=1.d0
      no=nogam2+nolog
      do
        df=cgamma2(s1,lm*(n+al))/es/zeroim(n+al)**s1
        f=f+df
        no=no+nogam2+nolog
        if(((abs(dble(df))+abs(imag(df)))/abs(f))**2 <= no*epso)then
          exit
        endif
        n=n+1.d0
      enddo
c      write(*,'(a,1p10g12.4)')'cpl-s1 ',dble(n),s,z,f,df
      u=-m_pi/lm**s
      if(s == czero)then
        f=f+u*berpol(0,al)
      else
        f=f+u*berpol(0,al)*csinp(.5d0*s)/(m_pi_2*s)
      endif
      no=50
      do m=1,ibmax,2
        u=-u*(0.d0,lm)/dble(m)
        if(s == dcmplx(dble(m),0.d0))then
          df=u*berpol(m,al)
        else
          df=u*berpol(m,al)
     $         *csinp(.5d0*(dble(m)-s))/(m_pi_2*(dble(m)-s))
        endif
        u=-u*(0.d0,lm)/dble(m+1)
        if(s == dcmplx(dble(m+1),0.d0))then
          df=df+u*berpol(m+1,al)
        else
          df=df+u*berpol(m+1,al)
     $         *csinp(.5d0*(dble(m+1)-s))/(m_pi_2*(dble(m+1)-s))
        endif
        f=f+df
        no=no+100
        if(((abs(dble(df))+abs(imag(df)))/abs(f))**2 <= no*epso)then
          exit
        endif
      enddo
      f=f/(m_2pi**s1)
      return
      end function

      complex*16 recursive function cpolylogz(s,z) result(f)
      implicit none
      complex*16 ,intent(in):: s,z
      real*8 la
      real*8 ,parameter :: epspl=log(epso*nopl)/2.d0
      if(z == czero)then
        f=czero
        return
      elseif(z == cone)then
        if(dble(s) > 1.d0)then
          f=czeta(s)-cone
        elseif(dble(s) == 0.d0)then
          f=0.d0/0.d0
        else
          f=1.d0/0.d0
        endif
        return
      endif
      la=log(abs(z))
      if(14.d0*max(0.d0,abs(la)) < epspl+dble(s)*log(16.d0/2.d0))then
        f=cpolylogs(s,z,czero,15)
      else
        f=cpolylog(s,z)-z
      endif
      return
      end function

      complex*16 function cplconv(s,z) result(f1)
      implicit none
      complex*16 ,intent(in):: z,s
      integer*4 ,parameter :: nmax=1000
      real*8 ,parameter ::nf=1.d5
      complex*16 u,z1,vn(0:nmax),df,df1,df0
      real*8 n,adf,adf0,v,v1,k,k1,n1,az
      integer*4 nv,no,i,j,nj
      u=-z/zeroim(1.d0-z)
      z1=u
      az=abs(imag(log(z1**2)))/m_2_pi
      do nj=1,4
        if(abs(nj*az-anint(nj*az)) <= 0.25d0)then
          exit
        endif
      enddo
      nj=min(1,nj)
      vn(0)=-cone
      nv=0
      f1=-z1
      no=nolog+8
      n=1.d0
      adf0=veryl
c      write(*,'(a,1p10g12.4)')'chlconv-0 ',z,s,a,u,f1
      do
        df=czero
        do j=1,nj
          n1=n+1.d0
          v=-1.d0
          v1=-1.d0
          k=0.d0
          df0=-cone
          df1=-cone
          do i=1,min(nint(n1),nmax)
            k1=k+1.d0
            v=-v*(n-k)/k1
            v1=-v1*(n1-k)/k1
            if(i > nv)then
              vn(i)=k1**(-s)
              nv=i
            endif
            df0=df0+v*vn(i)
            df1=df1+v1*vn(i)
            k=k1
          enddo
          u=u*z1
          df=df+u*(df0+z1*df1)
          u=u*z1
          n=n1+1.d0
        enddo
        no=no+nint(n)*(nolog*2+16)*nj+12
        adf=abs(dble(df))+abs(imag(df))
        if(adf > adf0)then
c          if(adf**2 > no*nf*abs(f1)**2*epso)then
c            f1=(veryl,0.d0)
c          endif
          exit
        endif
        adf0=adf
        f1=f1+df
c        write(*,'(a,1p10g12.4)')'cplc ',n,z,f1,df
        if(adf**2 <= no*abs(f1)**2*epso)then
          exit
        endif
      enddo
      return
      end function

      complex*16 recursive function cplconvg(s,z) result(f1)
      implicit none
      complex*16 ,intent(in):: z,s
      complex*16 df,sh,sh1,sh2,sh3,cgi,cgi1,c,u,un,un1,s1
      real*8 n,adf,adf0
      integer*4 no
      real*8 ,parameter :: lm=m_pi,lm1=m_pi**2/lm
      s1=1.d0-s
      sh=.5d0*s
      sh1=sh+.5d0
      sh2=zeroim(.5d0-sh)
      cgi=cgammai(sh)
      cgi1=cgammai(sh1)
c      f1=(cgamma2(sh,lm*a2)*cgi+cgamma2(sh1,lm*a2)*cgi1)/a2**sh
      f1=-lm**sh*cgi/s
      n=1.d0
      no=nogam2*2
      adf0=veryl
      u=cone
      do
        u=u*z
        df=(cgamma2(sh,dcmplx(lm*n**2,0.d0))*cgi*(u+1.d0/u)
     $       +cgamma2(sh1,dcmplx(lm*n**2,0.d0))*cgi1*(u-1.d0/u))
     $       /n**s
        adf=abs(dble(df))+abs(imag(df))
        if(adf > adf0)then
          exit
        endif
        adf0=adf
        no=no+nogam2*2
        f1=f1+df
        if(adf**2 <= no*abs(f1)**2*epso)then
          exit
        endif
        n=n+1.d0
      enddo
      f1=f1*.5d0
      sh3=zeroim(1.d0-sh)
      n=1.d0
      adf0=veryl
      c=.5d0*m_pi**(s-0.5d0)
      un=log(zeroim(z))/(0.d0,m_2pi)
      f1=f1+c*(cgamma2(sh2,lm1*un**2)*cgi
     $     +cgamma2(sh3,lm1*un**2)*cgi1*dcmplx(0.d0,
     $     merge(1.d0,merge(1.d0,-1.d0,dble(un)==0.d0),
     $     dble(un)>0.d0)))
     $     /un**s1
      un1=un
      do
        un=un+1.d0
        df=(cgamma2(sh2,lm1*un**2)*cgi
     $       +cgamma2(sh3,lm1*un**2)*cgi1*dcmplx(0.d0,
     $       merge(1.d0,merge(1.d0,-1.d0,dble(un)==0.d0),
     $       dble(un)>0.d0)))
     $       /un**s1
        un1=un1-1.d0
        df=c*(df+cgamma2(sh2,lm1*un1**2)*cgi
     $       +cgamma2(sh3,lm1*un1**2)*cgi1*dcmplx(0.d0,
     $       merge(1.d0,merge(1.d0,-1.d0,dble(un1)==0.d0),
     $       dble(un1)>0.d0))
     $       /un1**s1)
        adf=abs(dble(df))+abs(imag(df))
        if(adf > adf0)then
          exit
        endif
        adf0=adf
        f1=f1+df
c      if(abs(imag(a)) < 1.d0 .and. abs(dble(z)) < 2.d0
c     $     .and. abs(imag(a)) > 0.39d0)then
c        write(*,'(a,1p10g12.4)')'chlcg-22  ',z,a,un1,f1,df
c      endif
        no=no+(nocgam+nolog)*4
        if(adf**2 <= no*abs(f1)**2*epso)then
          exit
        endif
      enddo
      return
      end function

      complex*16  function loga(z) result(f)
      implicit none
      complex*16 ,intent(in):: z
      f=log(zeroim(z))
      if(imag(f) < 0.d0)then
        f=f+(0.d0,m_2pi)
      endif
      return
      end function

      complex*16 recursive function chlerchn(z,s,n) result(f1)
      implicit none
      integer*4 ,parameter :: nmax=1000
      complex*16 ,intent(in):: z,s
      integer*4 ,intent(in):: n
      complex*16 u
      real*8 ak
      integer*4 i
      if(n == 0)then
        f1=cpolylog(s,z)
      elseif(n == 1)then
        f1=cpolylog(s,z)/z
      elseif(imag(z) < 0.d0)then
        f1=conjz(chlerchn(conjz(z),conjz(s),n))
      elseif(n > 0)then
        f1=cpolylog(s,z)
        u=1.d0
        ak=1.d0
        do i=1,n-1
          u=u*z
          f1=f1-u*ak**(-s)
          ak=ak+1.d0
        enddo
        f1=f1*z**(-n)
      else
        f1=z**(-n)*cpolylog(s,z)+dble(n)**(-s)
        u=1.d0
        ak=1.d0
        do i=1,-n-1
          u=u*z
          f1=f1+u*ak**(-s)
          ak=ak+1.d0
        enddo
      endif
      return
      end function

      complex*16 recursive function chlerch(z,s,a) result(f1)
      implicit none
      complex*16 ,intent(in):: z,s,a
      complex*16 df,u,z1,ak
      real*8 m,m1,n,n1,az
      integer*4 no,nt,i,j,nj
      real*8 ,parameter ::mmax=50.d0
      if(z == czero)then
        f1=a**(-s)
      elseif(z == cone)then
        f1=chzeta2(s,a)
      elseif(imag(a) == 0.d0 .and. anint(dble(a)) == dble(a))then
        f1=chlerchn(z,s,nint(dble(a)))
      elseif(imag(z) < 0.d0)then
        f1=conjz(chlerch(conjz(z),conjz(s),conjz(a)))
      elseif(abs(a) > arth)then
        n=anint(dble(a)+0.1d0)
        nt=int(n)
        n1=dble(nt)
        f1=czero
        if(nt > 0)then
          if(a /= dcmplx(n1,0.d0))then
            f1=chlerch(z,s,a-n1)-zeroim(a-n1)**(-s)
          endif
          u=1.d0
          ak=zeroim(a-n1+1.d0)
          do i=1,nt-1
            u=u*z
            if(ak /= czero)then
              f1=f1-u*ak**(-s)
            endif
            ak=ak+1.d0
          enddo
          f1=f1*z**(-nt)
        elseif(nt < 0)then
          f1=z**(-nt)*chlerch(z,s,a-n1)+a**(-s)
          u=1.d0
          ak=zeroim(a+1.d0)
          do i=1,-nt-1
            u=u*z
            if(ak /= czero)then
              f1=f1+u*ak**(-s)
            endif
            ak=ak+1.d0
          enddo
        else
          z1=zeroim(z**2)
          f1=(chlerch(z1,s,.5d0*a)+z*chlerch(z1,s,.5d0*a+.5d0))
     $         *2.d0**(-s)
        endif
      else
        az=abs(imag(log(a)))/m_2pi
        do nj=1,4
          if(abs(nj*2.d0*az-anint(nj*2.d0*az)) <= 0.25d0)then
            exit
          endif
        enddo
        nj=min(4,nj)
        u=1.d0
        f1=a**(-s)+z*(1.d0+a)**(-s)+cpolylogz(s,z)
        m=0.d0
        no=nopl+nolog*2
        do
          df=0.d0
          do j=1,nj
            m1=m+1.d0
            u=-u*(m+s)*a/m1
            df=df+u*cpolylogz(s+m1,z)
            m=m1
          enddo
          f1=f1+df
          no=no+(nopl+14)*nj+1
c          write(*,'(a,1p10g12.4)')'chl ',m,s,u,df,f1
          if((abs(dble(df))+abs(imag(df)))**2
     $         <= no*abs(f1)**2*epso)then
            exit
          endif
        enddo
      endif
      return
      end function

      complex*16 function chlconv(z,s,a) result(f1)
      implicit none
      complex*16 ,intent(in):: z,s,a
      integer*4 ,parameter :: nmax=1000
      real*8 ,parameter ::nf=1.d5
      complex*16 u,z1,vn(0:nmax),df,df1,df0
      real*8 n,adf,adf0,v,v1,k,k1,n1,az
      integer*4 nv,no,i,j,nj
      u=1.d0/zeroim(1.d0-z)
      z1=-z*u
      az=abs(imag(log(z1**2)))/m_2_pi
      do nj=1,4
        if(abs(nj*az-anint(nj*az)) <= 0.25d0)then
          exit
        endif
      enddo
      nj=min(4,nj)
      if(a ==czero)then
        vn(0)=czero
      else
        vn(0)=zeroim(a)**(-s)
      endif
      nv=0
      f1=u*vn(0)
      no=nolog+8
      n=1.d0
      adf0=veryl
c      write(*,'(a,1p10g12.4)')'chlconv-0 ',z,s,a,u,f1
      do
        df=czero
        do j=1,nj
          n1=n+1.d0
          v=1.d0
          v1=1.d0
          k=0.d0
          df0=vn(0)
          df1=vn(0)
          do i=1,min(nint(n1),nmax)
            k1=k+1.d0
            v=-v*(n-k)/k1
            v1=-v1*(n1-k)/k1
            if(i > nv)then
              if(a+k1 /= czero)then
                vn(i)=zeroim(a+k1)**(-s)
              else
                vn(i)=czero
              endif
              nv=i
            endif
            df0=df0+v*vn(i)
            df1=df1+v1*vn(i)
            k=k1
          enddo
          u=u*z1
          df=df+u*(df0+z1*df1)
          u=u*z1
          n=n1+1.d0
        enddo
        no=no+nint(n)*(nolog*2+16)*nj+12
        adf=abs(dble(df))+abs(imag(df))
        if(adf > adf0)then
c          if(adf**2 > no*nf*abs(f1)**2*epso)then
c            f1=(veryl,0.d0)
c          endif
          exit
        endif
        adf0=adf
        f1=f1+df
c        write(*,'(a,1p10g12.4)')'hlp-n ',n,z,a,f1,df
        if(adf**2 <= no*abs(f1)**2*epso)then
          exit
        endif
      enddo
      return
      end function

      logical*4  function omegaa(z,a) result(f1)
      implicit none
      complex*16 ,intent(in):: z,a
      if(dble(a) > 0.d0)then
        f1=(imag(z) /= 0.d0 .or. dble(z) < 1.d0)
      else
        f1=abs(z) < 1.d0
      endif
      return
      end function

      complex*16 recursive function chlconvg(z,s,a) result(f1)
      implicit none
      complex*16 ,intent(in):: z,s,a
      complex*16 df,sh,sh1,sh2,sh3,a2,an,an2,an3,cgi,cgi1,c,
     $     u,un,un1,es,eu,f0
      real*8 n,adf,adf0
      integer*4 no
      real*8 ,parameter :: lm=m_pi,lm1=m_pi**2/lm
      a2=a**2
      sh=.5d0*s
      sh1=sh+.5d0
      sh2=zeroim(.5d0-sh)
      cgi=cgammai(sh)
      cgi1=cgammai(sh1)
      f1=(cgamma2(sh,lm*a2)*cgi+cgamma2(sh1,lm*a2)*cgi1)/a2**sh
      write(*,'(a,1p10g12.4)')'chlcg ',z,s,a,f1
      n=1.d0
      no=nogam2*2
      adf0=veryl
      u=1.d0
      do
        u=u*z
        an=zeroim(a+n)
        if(an /= czero)then
          an2=an**2
          df=(cgamma2(sh,lm*an2)*cgi+cgamma2(sh1,lm*an2)*cgi1)
     $     /an2**sh*u
        else
          df=czero
        endif
        an=zeroim(a-n)
        if(an /= czero)then
          an3=an**2
          df=df+(cgamma2(sh,lm*an3)*cgi-cgamma2(sh1,lm*an3)*cgi1)
     $         /an3**sh/u
        endif
        adf=abs(dble(df))+abs(imag(df))
        if(adf > adf0)then
          exit
        endif
        adf0=adf
        no=no+nogam2*4
        f1=f1+df
        if(adf**2 <= no*abs(f1)**2*epso)then
          exit
        endif
        n=n+1.d0
      enddo
      f1=f1*.5d0
      f0=f1
      sh3=zeroim(1.d0-sh)
      n=1.d0
      eu=1.d0
      adf0=veryl
      c=.5d0*m_pi**(s-0.5d0)/z**a
c      if(imag(z) == 0.d0 .and. dble(z) < 0.d0)then
c        un=log(abs(z))/(0.d0,m_2pi)+0.5d0
c      else
        un=log(zeroim(z))/(0.d0,m_2pi)
c      endif
      es=cexpp(2.d0*dcmplx(imag(a),-dble(a)))
      f1=f1+c*(cgamma2(sh2,lm1*un**2)*cgi
     $     +cgamma2(sh3,lm1*un**2)*cgi1*dcmplx(0.d0,
     $     merge(1.d0,merge(1.d0,-1.d0,dble(un)==0.d0),
     $     dble(un)>0.d0)))
     $     /(un*un)**sh2
      un1=un
c      if(abs(imag(a)) < 1.d0 .and. abs(dble(z)) < 2.d0
c     $     .and. abs(imag(a)) > 0.398d0)then
c        write(*,'(a,1p10g12.4)')'chlcg-20  ',z,a,un,un**2,f1
c        write(*,'(a,1p10g12.4)')':         ',f0,c
c     $       cgamma2(sh2,lm1*un**2),cgamma2(sh3,lm1*un**2)
c      endif
      do
        un=un+1.d0
        eu=eu*es
        df=eu*(cgamma2(sh2,lm1*un**2)*cgi
     $       +cgamma2(sh3,lm1*un**2)*cgi1*dcmplx(0.d0,
     $       merge(1.d0,merge(1.d0,-1.d0,dble(un)==0.d0),
     $       dble(un)>0.d0)))
     $       /(un*un)**sh2
c      if(abs(imag(a)) < 1.d0 .and. abs(dble(z)) < 2.d0
c     $     .and. abs(imag(a)) > 0.37d0)then
c        write(*,'(a,1p10g12.4)')'chlcg-21 ',z,a,un,df
c      endif
        un1=un1-1.d0
        df=c*(df+cgamma2(sh2,lm1*un1**2)*cgi
     $       +cgamma2(sh3,lm1*un1**2)*cgi1*dcmplx(0.d0,
     $       merge(1.d0,merge(1.d0,-1.d0,dble(un1)==0.d0),
     $       dble(un1)>0.d0))
     $       /(un1*un1)**sh2/eu)
        adf=abs(dble(df))+abs(imag(df))
        if(adf > adf0)then
          exit
        endif
        adf0=adf
        f1=f1+df
c      if(abs(imag(a)) < 1.d0 .and. abs(dble(z)) < 2.d0
c     $     .and. abs(imag(a)) > 0.39d0)then
c        write(*,'(a,1p10g12.4)')'chlcg-22  ',z,a,un1,f1,df
c      endif
        no=no+(nocgam+nolog)*4
        if(adf**2 <= no*abs(f1)**2*epso)then
          exit
        endif
      enddo
      return
      end function

      complex*16 recursive function clerch(z,s,a) result(f1)
      implicit none
      complex*16 ,intent(in):: z,s,a
      complex*16 es,afa1
      real*8 fa1
      integer*4 nfa1
      real*8 ,parameter :: lm=m_pi,lm1=m_pi
      if(z == czero)then
        f1=(a*a)**(-.5d0*s)
      elseif(z == cone)then
        f1=czeta2(s,a)
      elseif(dble(a) >= 0.d0)then
        f1=chlerch(z,s,a)
      elseif(imag(a) < 0.d0)then
        f1=conjz(clerch(conjz(z),conjz(s),conjz(a)))
      else
        es=cexpp(dcmplx(-imag(s),dble(s)))
        fa1=floor(-dble(a))
        nfa1=int(fa1)
        afa1=a+fa1
        if(afa1 /= czero)then
          f1=es*chlerch(z,s,a)+(1.d0-es)*z**nfa1
     $         *(z*chlerch(z,s,afa1+1.d0)-
     $         +(1.d0+fa1+floor(dble(a)))/(afa1**2)**(s/2))
        else
          f1=es*chlerch(z,s,a)+(1.d0-es)*z**nfa1
     $         *z*chlerch(z,s,afa1+1.d0)
        endif
      endif
      return
      end function

      complex*16 recursive function chg(a,b,c,x,reg) result(f)
      implicit none
      complex*16,intent(in):: a,b,c,x
      logical*4 ,intent(in):: reg
      real*8 ax,ax1
      if(dble(b) < dble(a))then
        f=chg(b,a,c,x,reg)
        return
      elseif(imag(a) == 0.d0)then
        if(imag(b) == 0.d0 .and. imag(c) == 0.d0
     $       .and. imag(x) == 0.d0 .and. dble(x) <= 1.d0)then
          f=dcmplx(hgrr(dble(a),dble(b),dble(c),dble(x),reg),0.d0)
          return
        elseif(anint(dble(a)) == dble(a) .and. dble(a) <= 0.d0)then
          f=chgp(dble(a),b,c,x,reg)
          return
        endif
      endif
      if(abs(x-xs1) < epsx .or. abs(x-xs2) < epsx)then
        f=chg1s(a,b,c,x,reg)
      elseif(dble(x) <= 0.5d0)then
        ax=abs(x)
        if(ax .ge. 1.d0)then
          f=chg1(a,b,c,x,reg)
c          write(*,'(a,1p10g12.4)')'chg-1 ',a,b,c,x,f
        else
          ax1=abs(x-1.d0)
          if(ax1 > 1.d0)then
            f=zeroim(1.d0-x)**(-a)*chg(a,c-b,c,zeroim(x/(x-1.d0)),reg)
c            write(*,'(a,1p10g12.4)')'chg-2 ',a,b,c,x,f
          else
            f=chg3(a,b,c,x,reg)
c            write(*,'(a,1p10g12.4)')'chg-3 ',a,b,c,x,f
          endif
        endif
      else
        ax1=abs(x-1.d0)
        if(ax1 > 1.d0)then
          f=chg6(a,b,c,x,reg)
c          write(*,'(a,1p10g12.4)')'chg-6 ',a,b,c,x,f
        else
          ax=abs(x)
          if(ax > 1.d0)then
            f=chg5(a,b,c,x,reg)
c            write(*,'(a,1p10g12.4)')'chg-5 ',a,b,c,x,f
          else
            f=chg4(a,b,c,x,reg)
c            write(*,'(a,1p10g12.4)')'chg-4 ',a,b,c,x,f
          endif
        endif
      endif
      return
      end function

      real*8 recursive function hgrr(a,b,c,x,reg) result(f)
      implicit none
      real*8 ,intent(in):: a,b,c,x
      logical*4 ,intent(in):: reg
      real*8 ,parameter ::bth=2.d0**32
      if(b < a)then
        f=hgrr(b,a,c,x,reg)
        return
      endif
      if(anint(a) == a .and. a <= 0.d0)then
        f=hgrp(a,b,c,x,reg)
      elseif(x <= -1.d0)then
        f=hgrr1(a,b,c,x,reg)
      elseif(x < 0.d0)then
        if(b < bth .or. x <= -0.5d0)then
          f=(1.d0-x)**(-a)*hgrr(a,c-b,c,x/(x-1.d0),reg)
        else
          f=hgrr3(a,b,c,x,reg)
        endif
      elseif(x == 0.d0)then
        if(reg)then
          f=gammai(c)
        else
          f=1.d0
        endif
      elseif(x <= 0.5d0)then
        f=hgrr3(a,b,c,x,reg)
      elseif(x <= 1.d0)then
        f=hgrr4(a,b,c,x,reg)
      else
        f=0.d0
      endif
      return
      end function

      complex*16 function chgp(a,b,c,x,reg) result(f)
      implicit none
      complex*16 ,intent(in):: b,c,x
      real*8 ,intent(in):: a
      logical*4 ,intent(in):: reg
      complex*16 g0,g1
      real*8 k
      integer*4 i
      if(reg)then
        f=cgammai(c)
      else
        f=cone
      endif
      if(a == 0.d0)then
        return
      endif
      g0=f
      f=f-f*b/c*x
      if(a == -1.d0)then
        return
      endif
      g1=f
      k=-1.d0
      do i=2,-nint(a)
        f=(k*(x-1.d0)*g0-(c-2.d0*k+(k-b)*x)*g1)/(k-c)
        g0=g1
        g1=f
        k=k-1.d0
      enddo
      return
      end function

      real*8 function hgrp(a,b,c,x,reg) result(f)
      implicit none
      real*8 ,intent(in):: a,b,c,x
      real*8 g0,g1,k
      integer*4 i
      logical*4 ,intent(in):: reg
      if(reg)then
        f=gammai(c)
      else
        f=1.d0
      endif
      if(a == 0.d0)then
        return
      endif
      g0=f
      f=f-f*b/c*x
      if(a == -1.d0)then
        return
      endif
      g1=f
      k=-1.d0
      do i=2,-nint(a)
        f=(k*(x-1.d0)*g0-(c-2.d0*k+(k-b)*x)*g1)/(k-c)
        g0=g1
        g1=f
        k=k-1.d0
      enddo
      return
      end function

      complex*16 function chg1s(a,b,c,x,reg) result(f1)
      implicit none
      complex*16 ,intent(in):: a,b,c,x
      logical*4 ,intent(in):: reg
      complex*16 f,u,x1,df,b2c,g,g0,g1
      real*8 k,k1
      integer*4 no
      u=zeroim(1.d0-.5d0*x)**(-a)
      if(reg)then
        g0=cgammai(c)
        no=nocgam
      else
        g0=1.d0
        no=20
      endif
      g1=g0*(1.d0-2.d0*b/c)
      x1=x/(x-2.d0)
      f=u*g0
      u=u*a*x1
      f=f+u*g1
      k=1.d0
      b2c=b*2.d0-c
      do
        k1=k+1.d0
        u=u*(a+k)/k1*x1
        g=(g0*k-b2c*g1)/(c+k)
        df=u*g
        f1=f+df
        no=no+10
c        write(*,'(a,i10,1p10g12.4)')'chg1s ',no,k1,u,f1,df
        if((abs(dble(df))+abs(imag(df)))**2
     $       <= no*abs(f1)**2*epso)then
          return
        endif
        f=f1
        k=k1
        g0=g1
        g1=g
      enddo
      end function

      complex*16 function chg1(a,b,c,x,reg) result(f1)
      implicit none
      complex*16 ,intent(in):: a,b,c,x
      logical*4 ,intent(in):: reg
      complex*16 sba,x1,df,u,lx,ba,b1,lc
      real*8 m,k,k1
      integer*4 i,no
      x1=zeroim(1.d0/(1.d0-x))
      ba=b-a
      m=anint(dble(ba))
      if(dble(ba) /= m .or. imag(ba) /= 0.d0)then
        sba=csinp(ba)/m_pi
        if(abs(sba) > epsba1)then
          lx=log(x1)
          if(reg)then
            f1=(chg(a,c-b,1.d0-ba,x1,.true.)
     $           *exp(a*lx-cloggamma(b)-cloggamma(c-a))
     $           -chg(b,c-a,ba+1.d0,x1,.true.)
     $           *exp(b*lx-cloggamma(a)-cloggamma(c-b))
     $           )/sba
          else
            lc=cloggamma(c)
            f1=(chg(a,c-b,1.d0-ba,x1,.true.)
     $           *exp(a*lx-cloggamma(b)-cloggamma(c-a)+lc)
     $           -chg(b,c-a,ba+1.d0,x1,.true.)
     $           *exp(b*lx-cloggamma(a)-cloggamma(c-b)+lc)
     $           )/sba
          endif
          return
        endif
      endif
      b1=a+m
      if(m /= 0.d0)then
        f1=exp(log_gamma(m)-cloggamma(b1)-cloggamma(c-a))
        u=f1
        do i=1,int(m)-1
          u=-u*(a+(i-1))*(c-b1+(i-1))/i/(m-i)*x1
          f1=f1+u
        enddo
        no=int(m)*4+nocgam*3
      else
        f1=czero
        no=0
      endif
      lx=-log(x1)
      u=exp(m*log(zeroim(-x1))-log_gamma(m+1.d0)
     $     -cloggamma(a)-cloggamma(c-b1))
c      u=zeroim(-x1)**m*gammai(m+1.d0)*cgammai(a)*cgammai(c-b1)
      f1=f1+u*(lx-m_euler+polygamma(m+1.d0)
     $     -cpolygamma(b1)-cpolygamma(c-a))
c      write(*,'(a,1p10g12.4)')'chg1 ',x,u,f1
      no=no+nocgam*3+nopg*3
      k=1.d0
      do
        k1=k+1.d0
        u=u*(b1+k-1.d0)*(c-a+k-1.d0)/k/(k+m)*x1
        df=u*(lx+polygamma(k1)+polygamma(m+k1)
     $       -cpolygamma(b1+k)-cpolygamma(c-a+k))
        no=no+nopg*4
        f1=f1+df
        if((abs(dble(df))+abs(imag(df)))**2
     $       <= no*abs(f1)**2*epso)then
          f1=f1*x1**a
          f1=f1+(ba-m)*(cpolygamma(b1)*(1.d0-f1)
     $         +cpolygamma(b1+1.d0)*a*b1/c*x)
          if(.not. reg)then
            f1=f1*cgamma(c)
          endif
          exit
        endif
        k=k1
      enddo
      return
      end function

      real*8 function hgrr1(a,b,c,x,reg) result(f1)
      implicit none
      real*8 ,intent(in):: a,b,c,x
      logical*4 ,intent(in):: reg
      real*8 sba,x1,m,k,k1,df,u,lx,ba,b1
      integer*4 i,no
      x1=1.d0/(1.d0-x)
      ba=b-a
      m=anint(ba)
      if(ba /= m)then
        sba=sinp(ba)/m_pi
        if(abs(sba) > epsba)then
          if(reg)then
            f1=(hgrr(a,c-b,1.d0-ba,x1,.true.)
     $           *x1**a*gammai(b)*gammai(c-a)
     $           -hgrr(b,c-a,ba+1.d0,x1,.true.)
     $           *x1**b*gammai(a)*gammai(c-b))/sba
          else
            f1=(hgrr(a,c-b,1.d0-ba,x1,.true.)
     $           *x1**a*gammai(b)*pochh(c-a,a)
     $           -hgrr(b,c-a,ba+1.d0,x1,.true.)
     $           *x1**b*gammai(a)*pochh(c-b,b))/sba
          endif
          return
        endif
      endif
      b1=a+m
      if(m /= 0.d0)then
        f1=pochh(b1,-a)*gammai(c-a)
c        f1=exp(log_gamma(m)-log_gamma(b1)-log_gamma(c-a))
        u=f1
        do i=1,int(m)-1
          u=-u*(a+(i-1))*(c-b1+(i-1))/i/(m-i)*x1
          f1=f1+u
        enddo
        no=int(m)*5+nogam
      else
        f1=czero
        no=0.d0
      endif
      lx=-log(x1)
      u=(-x1)**nint(m)*gammai(m+1.d0)*gammai(a)*gammai(c-b1)
c      u=(-x1)**nint(m)*exp(-log_gamma(m+1.d0)-log_gamma(a)
c     $     -log_gamma(c-b1))
      f1=f1+u*(lx-m_euler+polygamma(m+1.d0)
     $     -polygamma(b1)-polygamma(c-a))
      no=no+nogam*3+nopg*3
      k=1.d0
      do
        k1=k+1.d0
        u=u*(b1+k-1.d0)*(c-a+k-1.d0)/k/(k+m)*x1
        df=u*(lx+polygamma(k1)+polygamma(m+k1)
     $       -polygamma(b1+k)-polygamma(c-a+k))
        f1=f1+df
        no=no+nopg*4
        if(df**2 <= no*f1**2*epso)then
          f1=f1*x1**a
          f1=f1+(ba-m)*polygamma(b1)*(1.d0-f1)
          if(.not. reg)then
            f1=f1*gamma(c)
          endif
          return
        endif
        k=k1
      enddo
      return
      end function

      complex*16 function chg3(a,b,c,x,reg) result(f1)
      implicit none
      complex*16 ,intent(in):: a,b,c,x
      logical*4 ,intent(in):: reg
      complex*16 u
      real*8 s,s1
      integer*4 no
      if(imag(c) == 0.d0 .and.
     $     anint(dble(c)) == dble(c) .and. dble(c) <= 0.d0)then
        s=-dble(c)+1.d0
        if(reg)then
          u=exp(cloggamma(a+s)-cloggamma(a)+cloggamma(b+s)-cloggamma(b)
     $         -log_gamma(s+1.d0)+s*log(zeroim(x))-cloggamma(c))
          no=nocgam*6
        else
          u=exp(cloggamma(a+s)-cloggamma(a)+cloggamma(b+s)-cloggamma(b)
     $         -log_gamma(s+1.d0)+s*log(zeroim(x)))
          no=nocgam*5
        endif
c        u=cpochh(a,dcmplx(s,0.d0))*cpochh(b,dcmplx(s,0.d0))
c     $       *gammai(s+1.d0)*zeroim(x)**s
      else
        if(reg)then
          u=cgammai(c)
          no=nocgam
        else
          u=1.d0
          no=0
        endif
        s=0.d0
      endif
      f1=u
      do
        s1=s+1.d0
        u=u*(a+s)*(b+s)/(c+s)/s1*x
        f1=f1+u
        no=no+9
        if(abs(u)**2 <= no*abs(f1)**2*epso)then
          exit
        endif
        s=s1
      enddo
      end function

      real*8  function hgrr3(a,b,c,x,reg) result(f1)
      implicit none
      real*8 ,intent(in):: a,b,c,x
      logical*4 ,intent(in):: reg
      real*8 s,s1,u,x1
      real*8 ,parameter ::xthb=-20.d0
      integer*4 no
      x1=x*b
      if(anint(c) == c .and. c <= 0.d0)then
        s=-c+1.d0
        if(reg)then
          u=x**s*pochh(a,s)*pochh(b,s)*gammai(s+1.d0)*gammai(c)
c          u=exp(log_gamma(a+s)-log_gamma(a)+log_gamma(b+s)-log_gamma(b)
c     $         -log_gamma(s+1.d0)+s*log(x)-log_gamma(c))
          no=nogam*4+nolog
        else
          u=x**s*pochh(a,s)*pochh(b,s)*gammai(s+1.d0)
c          u=exp(log_gamma(a+s)-log_gamma(a)+log_gamma(b+s)-log_gamma(b)
c     $         -log_gamma(s+1.d0)+s*log(x))
          no=nogam*3+nolog
        endif
c        u=pochh(a,s)*pochh(b,s)*gammai(s+1.d0)*x**s
      else
        if(reg)then
          u=gammai(c)
          no=nogam
        else
          u=1.d0
          no=0
        endif
        s=0.d0
      endif
      f1=u
      do
        s1=s+1.d0
        u=u*(a+s)*(b+s)/(c+s)/s1*x
        f1=f1+u
        no=no+4
        if(u**2 <= no*f1**2*epso)then
          return
        endif
        s=s1
      enddo
      end function

      complex*16 recursive function chg4(a,b,c,x,reg) result(f1)
      implicit none
      complex*16 ,intent(in):: a,b,c,x
      logical*4 ,intent(in):: reg
      complex*16 cab,scab,x1,lx,u,df,c1,lc
      real*8 m,k,k1
      integer*4 i,no
      x1=zeroim(1.d0-x)
      cab=c-a-b
      m=anint(dble(cab))
      if(m /= dble(cab) .or. imag(cab) /= 0.d0)then
        scab=csinp(cab)/m_pi
        if(abs(scab) > epsba1)then
          lx=log(x1)
          if(reg)then
            f1=(chg(a,b,1.d0-cab,x1,.true.)
     $           *exp(-cloggamma(c-a)-cloggamma(c-b))
     $           -chg(c-a,c-b,cab+1.d0,x1,.true.)
     $           *exp(lx*cab-cloggamma(a)-cloggamma(b)))/scab
          else
            lc=cloggamma(c)
            f1=(chg(a,b,1.d0-cab,x1,.true.)
     $           *exp(-cloggamma(c-a)-cloggamma(c-b)+lc)
     $           -chg(c-a,c-b,cab+1.d0,x1,.true.)
     $           *exp(lx*cab-cloggamma(a)-cloggamma(b)+lc))/scab
          endif
          return
        endif
      endif
      if(dble(m) < 0.d0)then
        f1=x1**m*chg4(a+m,b+m,c,x,reg)
        return
      endif
      c1=a+b+m
      if(m /= 0.d0)then
c        u=cgammai(a+m)*cgammai(b+m)*gamma(m)
c        u=cgammai(a+m)*cpochh(b+m,-b)
        u=exp(-cloggamma(a+m)-cloggamma(b+m)+log_gamma(m))
        f1=u
        do i=1,int(m)-1
          u=-u*(a+(i-1))*(b+(i-1))/i/(m-i)*x1
          f1=f1+u
        enddo
        no=int(m)*3+nocgam*3
      else
        f1=czero
        no=0
      endif
      if(x1 /= czero)then
        lx=log(x1)
c        u=-zeroim(-x1)**m*cgammai(a)*cgammai(b)*gammai(m+1.d0)
        u=-exp(log(zeroim(-x1))*m-cloggamma(a)-cloggamma(b)
     $       -log_gamma(m+1.d0))
        f1=f1+u*(lx+m_euler-polygamma(m+1.d0)
     $       +cpolygamma(a+m)+cpolygamma(b+m))
        no=no+nocgam*3+nopg*3
        k=0.d0
        do
          k1=k+1.d0
          u=u*(a+m+k)*(b+m+k)/k1/(k1+m)*x1
          df=u*(lx-polygamma(k1+1.d0)-polygamma(m+k1+1.d0)
     $         +cpolygamma(a+m+k1)+cpolygamma(b+m+k1))
          f1=f1+df
          no=no+nopg*4
          if((abs(dble(df))+abs(imag(df)))**2
     $         <= no*abs(f1)**2*epso)then
            f1=f1-(cab-m)*(cpolygamma(c1)*(1.d0-f1)
     $           +cpolygamma(c1+1.d0)*a*b/c1*x)
            if(.not. reg)then
              f1=f1*cgamma(c)
            endif
            exit
          endif
          k=k1
        enddo
      endif
      return
      end function

      real*8 recursive function hgrr4(a,b,c,x,reg) result(f1)
      implicit none
      real*8 ,intent(in):: a,b,c,x
      logical*4 ,intent(in):: reg
      real*8 cab,scab,x1,m,k,k1,lx,u,df,c1
      integer*4 i,no
      x1=1.d0-x
      cab=c-a-b
      m=anint(cab)
      if(m /= cab)then
        scab=sinp(cab)/m_pi
        if(abs(scab) > epsba)then
          if(reg)then
            f1=(hgrr(a,b,1.d0-cab,x1,.true.)
     $           *gammai(c-a)*gammai(c-b)
     $           -hgrr(c-a,c-b,cab+1.d0,x1,.true.)
     $           *gammai(a)*gammai(b)*x1**cab)/scab
          else
            f1=(hgrr(a,b,1.d0-cab,x1,.true.)
     $           *gammai(c-a)*pochh(c-b,b)
     $           -hgrr(c-a,c-b,cab+1.d0,x1,.true.)
     $           *gammai(a)*pochh(b,c-b)*x1**cab)/scab
          endif
          return
        endif
      endif
      if(m < 0.d0)then
        f1=x1**m*hgrr4(a+m,b+m,c,x,reg)
        return
      endif
      c1=a+b+m
      if(m /= 0.d0)then
        u=gammai(a+m)*pochh(b+m,-b)
        f1=u
        do i=1,int(m)-1
          u=-u*(a+(i-1))*(b+(i-1))/i/(m-i)*x1
          f1=f1+u
        enddo
        no=int(m)*3+nogam*3
      else
        f1=0.d0
        no=0
      endif
      if(x1 /= 0.d0)then
        lx=log(x1)
        u=-(-x1)**m*gammai(a)*gammai(b)*gammai(m+1.d0)
        f1=f1+u*(lx+m_euler-polygamma(m+1.d0)
     $       +polygamma(a+m)+polygamma(b+m))
        no=no+nogam*3+nopg*3
        k=0.d0
        do
          k1=k+1.d0
          u=u*(a+m+k)*(b+m+k)/k1/(k1+m)*x1
          df=u*(lx-polygamma(k1+1.d0)-polygamma(m+k1+1.d0)
     $         +polygamma(a+m+k1)+polygamma(b+m+k1))
          no=no+nopg*4
          if(df**2 <= no*f1**2*epso)then
            f1=f1-(cab-m)*(polygamma(c1)*(1.d0-f1)
     $           +polygamma(c1+1.d0)*a*b/c1*x)
            if(.not. reg)then
              f1=f1*gamma(c)
            endif
            exit
          endif
          k=k1
        enddo
      endif
      return
      end function

      complex*16 recursive function chg5(a,b,c,x,reg) result(f1)
      implicit none
      complex*16 ,intent(in):: a,b,c,x
      logical*4 ,intent(in):: reg
      complex*16 cab,scab,x1,u,df,c1
      complex*16 clx,lc,lx,lx1
      real*8 k,k1,m
      integer*4 i,no
      x1=1.d0-1.d0/x
      cab=c-a-b
      m=anint(dble(cab))
      if(m /= cab .or. imag(cab) /= 0.d0)then
        scab=csinp(cab)/m_pi
        if(abs(scab) >= epsba1)then
          lx=log(x)
          lx1=log(zeroim(1.d0-x))
          if(reg)then
            f1=( chg(a,a-c+1.d0,1.d0-cab,x1,.true.)
     $           *exp(-a*lx-cloggamma(c-a)-cloggamma(c-b))
     $           -chg(c-a,1.d0-a,1.d0+cab,x1,.true.)
     $           *exp((a-c)*lx+lx1*cab-cloggamma(a)-cloggamma(b)))/scab
          else
            lc=cloggamma(c)
            f1=( chg(a,a-c+1.d0,1.d0-cab,x1,.true.)
     $           *exp(-a*lx-cloggamma(c-a)-cloggamma(c-b)+lc)
     $           -chg(c-a,1.d0-a,1.d0+cab,x1,.true.)
     $           *exp((a-c)*lx+lx1*cab-cloggamma(a)-cloggamma(b)+lc)
     $           )/scab
          endif
          return
        endif
      endif
      if(m < 0.d0)then
        f1=zeroim(1.d0-x)**m*chg5(a+m,b+m,c,x,reg)
        return
      endif
      c1=m+a+b
      if(m /= 0.d0)then
        u=cpochh(b+m,-b)*cgammai(a+m)
        f1=u
        do i=1,int(m)-1
          u=u*(a+(i-1))*(b+m-i)/i/(m-i)*x1
          f1=f1+u
        enddo
        no=int(m)*4+nocgam*2
      else
        f1=czero
        no=0
      endif
      clx=log(zeroim(-x1))
      u=-x1**m*gammai(m+1.d0)*cgammai(a)
      if(imag(b) == 0.d0 .and. dble(b) <= 0.d0 .and.
     $     dble(b) == anint(dble(b)))then
        f1=f1+u*(-1.d0)**b*cgamma(1.d0-b)
      else
        f1=f1+u*(clx+m_euler-polygamma(m+1.d0)
     $       +cpolygamma(a+m)+cpolygamma(b))*cgammai(b)
        no=no+nopg*4
      endif
      no=no+nocgam*3+nopg*3
      k=0.d0
      do
        k1=k+1.d0
        u=-u*(a+m+k)/k1/(k1+m)*x1
        if(imag(b) == 0.d0 .and. dble(b) <= k1 .and.
     $       dble(b-k1) == anint(dble(b-k1)))then
          df=u*(-1.d0)**(k-b)*cgamma(k1-b+1.d0)
        else
          df=u*(clx-polygamma(k1+1.d0)-polygamma(m+k1+1.d0)
     $         +cpolygamma(a+m+k1)+cpolygamma(b-k1))
     $         *cgammai(b-k1)
          no=no+nopg*4
        endif
        no=no+nocgam
        f1=f1+df
        if((abs(dble(df))+abs(imag(df)))**2
     $       <= no*abs(f1)**2*epso)then
          f1=f1*zeroim(x)**(-a)
          f1=f1-(cab-m)*(cpolygamma(c1)*(1.d0-f1)
     $         +cpolygamma(c1+1.d0)*a*b/c1*x)
          if(.not. reg)then
            f1=f1*cgamma(c)
          endif
          exit
        endif
        k=k1
      enddo
      return
      end function

      complex*16 function chg6(a,b,c,x,reg) result(f1)
      implicit none
      complex*16 ,intent(in):: a,b,c,x
      logical*4 ,intent(in):: reg
      complex*16 ba,sba,x1,u,df,d,clx,b1,lc,lx
      real*8 k,k1,m
      integer*4 i,no
      x1=1.d0/x
      ba=b-a
      m=anint(dble(ba))
      if(imag(ba) /= 0.d0 .or. m /= dble(ba))then
        sba=csinp(ba)/m_pi
        if(abs(sba) > epsba1)then
          lx=log(zeroim(-x))
          if(reg)then
            f1=(chg(a,a-c+1.d0,1.d0-ba,x1,.true.)
     $           *exp(-a*lx-cloggamma(b)-cloggamma(c-a))
     $           -chg(b,b-c+1.d0,ba+1.d0,x1,.true.)
     $           *exp(-b*lx-cloggamma(a)-cloggamma(c-b)))/sba
          else
            lc=cloggamma(c)
            f1=(chg(a,a-c+1.d0,1.d0-ba,x1,.true.)
     $           *exp(-a*lx-cloggamma(b)-cloggamma(c-a)+lc)
     $           -chg(b,b-c+1.d0,ba+1.d0,x1,.true.)
     $           *exp(-b*lx-cloggamma(a)-cloggamma(c-b)+lc))/sba
          endif
          return
        endif
      endif
      b1=a+m
      if(m /= 0.d0)then
        u=exp(log_gamma(m)-cloggamma(b1)-cloggamma(c-a))
c        u=cpochh(b1,m-b1)*cgammai(c-a)
        f1=u
        do i=1,int(m)-1
          u=u*(a+(i-1))*(c-a-i)/i/(m-i)*x1
          f1=f1+u
        enddo
        no=int(m)*4+nocgam*2
      else
        f1=czero
        no=0
      endif
      clx=log(zeroim(-x))
c      u=x**(-m)*cgammai(a)*gammai(m+1.d0)
      u=exp(-m*log(x)-cloggamma(a)-log_gamma(m+1.d0))
      d=c-b1
      if(imag(d) == 0.d0 .and.
     $     anint(dble(d)) == dble(d) .and. dble(d) <= 0.d0)then
        f1=f1+u*(-1.d0)**dble(d)*gamma(1.d0-dble(d))
        no=no+nocgam*3
      else
        f1=f1+u*(clx-m_euler+polygamma(m+1.d0)
     $       -cpolygamma(b1)-cpolygamma(d))*cgammai(d)
        no=no+nocgam*3+nopg*3
      endif
      k=0.d0
      do
        k1=k+1.d0
        d=c-b1-k1
        u=-u*(b1+k)*x1/k1/(k1+m)
        if(imag(d) == 0.d0 .and.
     $       anint(dble(d)) == dble(d) .and. dble(d) <= 0.d0)then
          df=u*(-1.d0)**dble(d)*gamma(1.d0-dble(d))
          no=no+nogam
        else
          df=u*(clx+polygamma(k1+1.d0)+polygamma(m+k1+1.d0)
     $         -cpolygamma(b1+k1)-cpolygamma(d))*cgammai(d)
          no=no+nogam+nopg*3
        endif
        f1=f1+df
        if((abs(dble(df))+abs(imag(df)))**2
     $       <= no*abs(f1)**2*epso)then
          f1=f1*zeroim(-x)**(-a)
          f1=f1+(ba-m)*(cpolygamma(b1)*(1.d0-f1)
     $         +cpolygamma(b1+1.d0)*a*b1/c*x)
          if(.not. reg)then
            f1=f1*cgamma(c)
          endif
          exit
        endif
        k=k1
      enddo
      return
      end function

      complex*16 recursive function confhg0(c,x0) result(f1)
      implicit none
      complex*16 ,intent(in):: c,x0
      complex*16 f,u,x
      real*8 s,s1
      integer*4 no
      x=zeroim(x0)
      if(imag(c) == 0.d0 .and. imag(x) == 0.d0)then
        f1=confhgrr0(dble(c),dble(x))
      else
        if(imag(c) == 0.d0 .and.
     $       anint(dble(c)) == dble(c) .and. dble(c) <= 0.d0)then
          s=-dble(c)+1.d0
          u=exp(log(x)*s-log_gamma(s+1.d0))
c          u=gammai(s+1.d0)*x**s
        else
          u=cgammai(c)
          s=0.d0
        endif
        no=nocgam
        f=u
        do
          s1=s+1.d0
          u=u/(c+s)/s1*x
          f1=f+u
          no=no+5
          if(abs(u)**2 <= no*abs(f1)**2*epso)then
            exit
          endif
          f=f1
          s=s1
        enddo
      endif
      return
      end function

      real*8 recursive function confhgrr0(c,x) result(f1)
      implicit none
      real*8 ,intent(in):: c,x
      real*8 f,s,s1,u
      integer*4 no
      if(anint(c) == c .and. c <= 0.d0)then
        s=-c+1.d0
c        u=gammai(s+1.d0)*x**s
        u=exp(s*log(x)-log_gamma(s+1.d0))
      else
        u=gammai(c)
        s=0.d0
      endif
      no=nogam
      f=u
      do
        s1=s+1.d0
        u=u/(c+s)/s1*x
        f1=f+u
        if(u**2 <= no*f1**2*epso)then
          return
        endif
        f=f1
        s=s1
      enddo
      return
      end function

      complex*16 recursive function confhg1(a,c,x,reg) result(f1)
      implicit none
      complex*16 ,intent(in):: a,c,x
      complex*16 u,ach(0:2),cch(0:0),ac,xa,xac,lx
      real*8 s,s1
      logical*4 ,intent(in) , optional:: reg
      logical*4 reg1
      integer*4 no
      if(present(reg))then
        reg1=reg
      else
        reg1=.true.
      endif
      if(imag(c) == 0.d0)then
        if(imag(a) == 0.d0 .and. imag(x) == 0.d0)then
          f1=confhgrr1(dble(a),dble(c),dble(x),reg1)
          return
        elseif(dble(c) == anint(dble(c)) .and. dble(c) <= 0.d0)then
          if(reg1)then
            f1=cpochh(a,1.d0-c)*x**(1.d0-c)
     $           *confhg1(a-c+1.d0,2.d0-c,x,.true.)
          else
            f1=1.d0/czero
          endif
          return
        endif
      endif
      if(a == c)then
        if(reg1)then
          f1=exp(x-cloggamma(c))
        else
          f1=exp(x)
        endif
        return
      elseif(abs(x) > chgth+max(abs(a),abs(c-a)))then
c        write(*,'(a,1p10g12.4)')'chg1-th ',a,c,x
        do
          if(imag(x) == 0.d0)then
            if(dble(x) >= 0.d0)then
              if(imag(a) == 0.d0 .and. dble(a) == anint(dble(a)))then
                xa=dcmplx((-dble(x))**(-nint(dble(a))),0.d0)
              else
                exit
              endif
            else
              if(imag(a) == 0.d0)then
                xa=dcmplx((-dble(x))**(-dble(a)),0.d0)
              else
                exit
              endif
            endif
          else
            lx=log(-x)
            if(dble(x) >= 0.d0 .and. imag(x) > 0.d0)then
              lx=lx+(0.d0,m_2pi)
            endif
            xa=exp(-a*lx)
          endif
          ac=a-c
          if(imag(x) == 0.d0)then
            if(dble(x) < 0.d0)then
              if(imag(ac) == 0.d0 .and. dble(ac) == anint(dble(ac)))then
                xac=dcmplx((dble(x))**nint(dble(ac)),0.d0)
              else
                exit
              endif
            else
              if(imag(ac) == 0.d0)then
                xac=dcmplx((dble(x))**(dble(ac)),0.d0)
              else
                exit
              endif
            endif
          else
            lx=log(x)
            if(dble(x) < 0.d0 .and. imag(x) > 0.d0)then
              lx=lx-(0.d0,m_2pi)
            endif
            xac=exp(ac*lx)
          endif
          ach(0)=czero
          cch(0)=czero
          ach(1)=a
          ach(2)=ac+1.d0
          if(reg)then
            f1=cgammai(c-a)*xa*chgpq(ach,cch,zeroim(-1.d0/x),.false.)
          else
            f1=cpochh(c-a,a)*xa*chgpq(ach,cch,zeroim(-1.d0/x),.false.)
          endif
          ach(1)=-ac
          ach(2)=1.d0-a
          if(reg)then
            f1=f1+cgammai(a)*xac*exp(x)*chgpq(ach,cch,1.d0/x,.false.)
          else
            f1=f1+cpochh(a,c-a)*xac*exp(x)*chgpq(ach,cch,1.d0/x,.false.)
          endif
c          write(*,'(a,1p10g12.4)')'chg1-pq ',a,c,x,f1
          return
        enddo
      endif
      if(dble(x) .ge. 0.d0)then
        if(reg1)then
          u=cgammai(c)
          no=nocgam
        else
          u=1.d0
          no=0
        endif
        s=0.d0
        f1=u
        if(a == cone)then
          do
            u=u/(c+s)*x
            f1=f1+u
            no=no+8
            if(abs(u)**2 <= no*abs(f1)**2*epso)then
              exit
            endif
            s=s+1.d0
          enddo
        else
          do
            s1=s+1.d0
            if(a+s == czero)then
              exit
            else
              u=u*(a+s)/(c+s)/s1*x
            endif
            f1=f1+u
            no=no+14
            if(abs(u)**2 <= no*abs(f1)**2*epso)then
              exit
            endif
            s=s1
          enddo
        endif
c        write(*,'(a,1p10g12.4)')'chg1-nt ',a,c,x,f1
      else
        f1=exp(x)*confhg1(c-a,c,-x,reg1)
      endif
      return
      end function

      real*8 recursive function confhgrr1(a,c,x,reg) result(f1)
      implicit none
      real*8 ,intent(in):: a,c,x
      real*8 s,s1,u,ac,xa,xac
      complex*16 ach(0:2),cch(0:0)
      integer*4 no
      logical*4 ,intent(in), optional:: reg
      logical*4 reg1
      if(present(reg))then
        reg1=reg
      else
        reg1=.true.
      endif
      if(c == anint(c) .and. c <= 0.d0)then
        if(reg1)then
          f1=pochh(a,1.d0-c)*x**(1.d0-c)
     $         *confhgrr1(a-c+1.d0,2.d0-c,x,.true.)
        else
          f1=1.d0/0.d0
        endif
        return
      elseif(a == c)then
        if(reg)then
          f1=exp(x)*gammai(c)
        else
          f1=exp(x)
        endif
        return
      elseif(abs(x) > chgth)then
        do
          if(x >= 0.d0)then
            if(a == anint(a))then
              xa=(-x)**(-nint(a))
            else
              exit
            endif
          else
            xa=(-x)**(-a)
          endif
          ac=a-c
          if(x < 0.d0)then
            if(ac == anint(ac))then
              xac=x**nint(ac)
            else
              exit
            endif
          else
            xac=x**ac
          endif
          ach(0)=cone
          cch(0)=cone
          ach(1)=dcmplx(a,0.d0)
          ach(2)=dcmplx(ac+1.d0,0.d0)
          if(reg)then
            f1=gammai(c-a)*xa*dble(hgpq(ach,cch,-1.d0/x,.false.))
          else
            f1=pochh(c-a,a)*xa*dble(hgpq(ach,cch,-1.d0/x,.false.))
          endif
c          write(*,'(a,1p10g12.4)')'hg1-1 ',a,c,xa,f1
          ach(1)=dcmplx(-ac,0.d0)
          ach(2)=dcmplx(1.d0-a,0.d0)
          if(reg)then
            f1=f1+gammai(a)
     $           *xac*exp(x)*dble(hgpq(ach,cch,1.d0/x,.false.))
          else
            f1=f1+pochh(a,c-a)
     $           *xac*exp(x)*dble(hgpq(ach,cch,1.d0/x,.false.))
          endif
c          write(*,'(a,1p10g12.4)')'hg1-2 ',a,c,xac,f1
          return
        enddo
      endif
      if(x >= 0.d0)then
        if(reg1)then
          u=gammai(c)
          no=nogam
        else
          u=1.d0
          no=0
        endif
        s=0.d0
        f1=u
        if(a == 1.d0)then
          do
            s1=s+1.d0
            if(s1 == 0.d0)then
              exit
            else
              u=u/(c+s)*x
            endif
            f1=f1+u
            no=no+2
c            write(*,'(a,i5,1p10g12.4)')'hg1-2 ',s,a,c,x,u,f1
            if(u**2 <= no*abs(f1)**2*epso)then
              exit
            endif
            s=s1
          enddo
        else
          do
            if(a+s == 0.d0)then
              exit
            else
              s1=s+1.d0
              u=u*(a+s)/(c+s)/s1*x
            endif
            f1=f1+u
            no=no+3
            if(u**2 <= no*f1**2*epso)then
              exit
            endif
            s=s1
          enddo
        endif
      else
        f1=exp(x)*confhgrr1(c-a,c,-x,reg1)
      endif
      return
      end function

      recursive function tfhg(isp1,irtc) result(kx)
      use tfstk
      implicit none
      type (sad_descriptor) kx
      integer*4 ,intent(in):: isp1
      integer*4 ,intent(out):: irtc
      integer*4 iconf
      logical*4 reg
      reg=.false.
      iconf=0
      if(isp == isp1+5)then
        reg=.true.
        if(ktastk(isp1+3) == dxnullo%k)then
          iconf=2
        elseif(ktastk(isp1+4) == dxnullo%k)then
          iconf=1
        endif
      elseif(isp == isp1+2)then
        iconf=2
      elseif(isp == isp1+3)then
        iconf=1
      elseif(isp /= isp1+4)then
        go to 9000
      endif
      kx=kxhg(dtastk(isp1+1),dtastk(isp1+2),dtastk(isp1+3-iconf),
     $     dtastk(isp1+4-iconf),reg,iconf,irtc)
      return
 9000 irtc=-1
      kx=dxnullo
      return
      end function

      recursive function kxhg(ka,kb,kc,k,reg,iconf,irtc) result(kx)
      use tfstk
      implicit none
      type (sad_descriptor) kx
      type (sad_descriptor) ,intent(in):: k,ka,kb,kc
      type (sad_dlist),pointer::kxl,kl
      type (sad_rlist),pointer::klv
      integer*4 ,intent(out):: irtc
      complex*16 ca,cb,cc,cx
      integer*4 ,intent(in):: iconf
      integer*4 i
      logical*4 ,intent(in):: reg
      logical*4 d
      irtc=0
      if(.not. tfnumberq(kc,cc) .or.
     $     (iconf /= 2 .and. .not. tfnumberq(ka,ca))
     $       .or. (iconf == 0 .and. .not. tfnumberq(kb,cb)))then
        kx=dxnullo
        irtc=-1
        return
      endif
      if(tfnumberq(k,cx))then
        kx=kxhgc(ca,cb,cc,cx,iconf,reg)
      elseif(tfreallistq(k,klv))then
        kx=kxadaloc(-1,klv%nl,kxl)
        d=.false.
        do i=1,klv%nl
          kxl%dbody(i)=dtfcopyd(kxhgc(ca,cb,cc,
     $         dcmplx(klv%rbody(i),0.d0),iconf,reg),d)
        enddo
        if(.not. d)then
          kxl%attr=ior(kxl%attr,lnonreallist)-lnonreallist
        endif
      elseif(tflistq(k,kl))then
        kx=kxadaloc(-1,kl%nl,kxl)
        d=.false.
        do i=1,kl%nl
          kxl%dbody(i)=dtfcopyd(kxhg(ka,kb,kc,kl%dbody(i),
     $         reg,iconf,irtc),d)
          if(irtc /= 0)then
            kxl%dbody(i:kl%nl)=dxnullo
            exit
          endif
        enddo
        if(.not. d)then
          kxl%attr=ior(kxl%attr,lnonreallist)-lnonreallist
        endif
      else
        kx=dxnullo
        irtc=-1
      endif
      return
      end function 

      recursive function kxhgc(a,b,c,x,iconf,reg) result(kx)
      use tfstk
      implicit none
      type (sad_descriptor) kx
      complex*16 ,intent(in):: a,b,c,x
      integer*4 ,intent(in):: iconf
      complex*16 cx
      logical*4 ,intent(in):: reg
      select case (iconf)
      case (0)
        cx=chg(a,b,c,x,reg)
      case (1)
        cx=confhg1(a,c,x,reg)
      case (2)
        cx=confhg0(c,x)
        if(.not. reg)then
          cx=cx*cgamma(c)
        endif
      case (3)
        cx=chgu(a,c,x)
      case (4)
        cx=chlerch(x,a,c)
      case (5)
        cx=clerch(x,a,c)
      case default
        kx=dxnullo
        return
      end select
      kx=kxcalocc(-1,cx)
      return
      end function

      function kxhgpq(isp1,reg,irtc) result(kx)
      use tfstk
      implicit none
      type (sad_descriptor) kx
      type (sad_dlist),pointer ::kla,klb
      integer*4 ,intent(in):: isp1
      integer*4 ,intent(out):: irtc
      logical*4 ,intent(in):: reg
      complex*16 ,allocatable::a(:),b(:)
      integer*4 na,nb,m
      logical*4 cmplm,realm,veca,vecb
      irtc=-1
      kx=dxnullo
      call tfmatrixmaybeq(dtastk(isp1+1),cmplm,realm,veca,na,m,kla)
      if(m /= 0)then
        return
      endif
      call tfmatrixmaybeq(dtastk(isp1+2),cmplm,realm,vecb,nb,m,klb)
      if(m /= 0)then
        return
      endif
      allocate(a(0:na))
      call tfl2cm(kla,a(1:na),na,0,.true.,irtc)
      if(irtc /= 0)then
        return
      endif
      a(0)=merge(cone,czero,veca)
      allocate(b(0:nb))
      call tfl2cm(klb,b(1:nb),nb,0,.true.,irtc)
      if(irtc /= 0)then
        return
      endif
      b(0)=merge(cone,czero,vecb)
      kx=kxhgpqa(a,b,dtastk(isp),reg,irtc)
      return
      end function

      recursive function kxhgpqa(a,b,k,reg,irtc) result(kx)
      use tfstk
      implicit none
      type (sad_descriptor) kx
      type (sad_descriptor) ,intent(in):: k
      type (sad_dlist) ,pointer ::kl
      type (sad_dlist) ,pointer ::kxl
      complex*16 ,intent(in):: a(0:),b(0:)
      logical*4 ,intent(in):: reg
      integer*4 ,intent(out):: irtc
      complex*16 z
      integer*4 i
      logical*4 d
      irtc=0
      if(tfnumberq(k,z))then
        kx=kxcalocc(-1,chgpq(a,b,z,reg))
      elseif(tflistq(k,kl))then
        kx=kxadaloc(-1,kl%nl,kxl)
        d=.false.
        do i=1,kl%nl
          kxl%dbody(i)=dtfcopyd(kxhgpqa(a,b,kl%dbody(i),reg,irtc),d)
          if(irtc /= 0)then
            kxl%dbody(i:kl%nl)=dxnullo
            return
          endif
        enddo
        if(.not. d)then
          kxl%attr=ior(kxl%attr,lnonreallist)-lnonreallist
        endif
      endif
      return
      end function

      complex*16  function chgpq(a,b,z,reg) result(f)
      implicit none
      complex*16 ,intent(in):: z,a(0:),b(0:)
      logical*4 ,intent(in):: reg
      integer*4 na,nb
      complex*16 u
      real*8 k,au0,au
      integer*4 i,no
      if(imag(z) == 0.d0 .and.
     $     a(0) /= czero .and. b(0) /= czero)then
        f=hgpq(a,b,dble(z),reg)
        return
      endif
      na=size(a)-1
      nb=size(b)-1
c      write(*,'(a,2i5,1p10g12.4)')'chgpq ',na,nb,b(1)
      if(nb == 0 .and. na == 0)then
        f=exp(z)
        return
      elseif(nb == 1)then
        select case (na)
        case (0)
          if(reg)then
            f=confhg0(b(1),z)
          else
            f=confhg0(b(1),z)*cgammai(b(1))
          endif
          return
        case (1)
          f=confhg1(a(1),b(1),z,reg)
          return
        case (2)
          f=chg(a(1),a(2),b(1),z,reg)
          return
        end select
      endif
      u=czero
      if(reg)then
        do i=1,nb
          u=u+cloggamma(b(i))
        enddo
        no=nb*nocgam
      endif
      u=exp(-u)
      f=u
      k=0.d0
      no=0
      au0=veryl
      main: do
        do i=1,na
          u=u*(a(i)+k)
          if(u == czero)then
            return
          endif
        enddo
        do i=1,nb
          u=u/(b(i)+k)
        enddo
        k=k+1.d0
        u=u*z/k
        au=abs(u)
        if(au > au0)then
c          write(*,'(a,1p10g12.4)')'hgpq-nc ',k,au,au0,u,f
          return
        endif
        au0=au
        no=no+2*(na+nb+1)
        f=f+u
        if(au**2 <= no*abs(f)**2*epso)then
c          write(*,'(a,1p10g12.4)')'hgpq-cv ',k,au,au0,u,f
          return
        endif
      enddo main
      end function 

      complex*16  function hgpq(a,b,z,reg) result(f)
      implicit none
      real*8 ,intent(in):: z
      complex*16 ,intent(in):: a(0:),b(0:)
      logical*4 ,intent(in):: reg
      integer*4 na,nb
      real*8 u,k,f1,au,au0
      integer*4 i,no
      na=size(a)-1
      nb=size(b)-1
      if(nb == 0 .and. na == 0)then
        f=exp(z)
        return
      elseif(nb == 1)then
        select case (na)
        case (0)
          if(reg)then
            f=confhgrr0(dble(b(1)),z)*gammai(dble(b(1)))
          else
            f=confhgrr0(dble(b(1)),z)
          endif
          return
        case (1)
          if(reg)then
            f=confhgrr1(dble(a(1)),dble(b(1)),z)*gammai(dble(b(1)))
          else
            f=confhgrr1(dble(a(1)),dble(b(1)),z)
          endif
          return
        case (2)
          f=chg(a(1),a(2),b(1),dcmplx(z,0.d0),reg)
          return
        end select
      endif
      u=1.d0
      if(reg)then
        do i=1,nb
          u=u*gamma(dble(b(i)))
        enddo
        no=nb*nogam
      endif
      f1=u
      k=0.d0
      no=0
      au0=veryl
      main: do
        do i=1,na
          u=u*(dble(a(i))+k)
          if(u == 0.d0)then
c            write(*,'(a,3i5,1p10g12.4)')'hgpq-u0 ',i,na,nb,k,a(i),u
            f=f1
            return
          endif
        enddo
        do i=1,nb
          u=u/(dble(b(i))+k)
        enddo
        k=k+1.d0
        u=u*z/k
        au=abs(u)
        if(au > au0)then
          f=f1
          return
        endif
        au0=au
        no=no+2*(na+nb+1)
        f1=f1+u
c        write(*,'(a,3i5,1p10g12.4)')'hgpq-0 ',i,na,nb,z,u,f
        if(u**2 <= no*f1**2*epso)then
          f=f1
          return
        endif
      enddo main
      end function 

      real*8  function pgin(t) result(f1)
      implicit none
      real*8 ,intent(in):: t
      f1=polygamma(pgz*t+1.d0)*(1.d0-t)**pgn
      return
      end

      real*8  function pgin1(t) result(f1)
      implicit none
      real*8 ,intent(in):: t
      f1=gpolygamma2(2.d0,pgz*t+1.d0)*(1.d0-t)**pgn
      return
      end

      real*8  function pgip(t) result(f1)
      implicit none
      real*8 ,intent(in):: t
      f1=gpolygamma2(pgnc,pgz*t+1.d0)*(1.d0-t)**pgn
      return
      end

      complex*16  function cpgin(t) result(f1)
      implicit none
      real*8 ,intent(in):: t
      f1=cpolygamma(cpgz*t+1.d0)*(1.d0-t)**cpgn
      return
      end

      complex*16  function cpgin1(t) result(f1)
      implicit none
      real*8 ,intent(in):: t
      f1=cgpolygamma2((2.d0,0.d0),cpgz*t+1.d0)*(1.d0-t)**cpgn
      return
      end

      complex*16  function cpgip(t) result(f1)
      implicit none
      real*8 ,intent(in):: t
      f1=cgpolygamma2(dcmplx(pgnc,0.d0),cpgz*t+1.d0)*
     $     (1.d0-t)**cpgn
      return
      end

      complex*16 function cpolygamma2(n,x0) result(f1)
      use tfstk, only:ktfenanq
      implicit none
      complex*16 ,intent(in):: n,x0
      complex*16 u,f,lx,df,x
      complex*16 ,external:: cbint
      integer*4 i
      integer*4 ,parameter :: nmax=1000
      complex*16 xk,xn,xk1,sabcp,cab,f10,xkk,xkk1,j1n,sabc
      complex*16 ,save :: ns=dcmplx(-veryl,0.d0),cgi(-2:nmax),
     $     chk1(0:nmax),chk3(0:nmax),chk2(0:nmax),
     $     chk4(0:nmax),chk5(0:nmax),chk6(0:nmax)
      real*8 j,j1,axk,m,k,dfa
      integer*4 ,save :: im1,im2,im3,im4,im5,im6
      integer*4 jj,icg,no
      x=zeroim(x0)
      if(x == czero)then
        if(dble(n) .ge. -1.d0)then
          f1=dcmplx(1.d0/0.d0,0.d0)
        else
          f1=czero
        endif
        return
      elseif(imag(n) == 0.d0)then
        if(dble(n) == 0.d0)then
          f1=cpolygamma(x)
          return
        elseif(imag(x) == 0.d0 .and. dble(x) .ge. 0.d0)then
          f1=dcmplx(polygamma2(dble(n),dble(x)),0.d0)
          return
        endif
        if(dble(n) == anint(dble(n)) .and.
     $         dble(n) > 0.d0)then
          f1=(-1.d0)**(dble(n)+1.d0)*factorial(dble(n))
     $         *chzeta2(n+1.d0,x)
          return
        endif
      endif
      if(abs(x) < xmth)then
        u=x*cgammai(2.d0-n)
        f=(m_euler-log(x)+cpolygamma(-n))*cgammai(-n)/x
     $       -m_euler*cgammai(1.d0-n)+zt2*u
        k=2.d0
        no=nocgam
        do
          u=-u*x*k/(k-n)
          df=zeta(k+1.d0)*u
          f=f+df
          if((abs(dble(df))+abs(imag(df)))**2
     $         <= no*abs(f)**2*epso)then
            f1=f/x**n
            exit
          endif
          k=k+1.d0
          no=no+nozt+10
        enddo
      elseif(abs(n+1.d0)+abs(x-1.d0) < xmeps)then
        f1=(n+1.d0)*(dpmg1+ddpmg1*.5d0*(n+1.d0)+
     $       dpzmg1*(x-1.d0))
     $       +(x-1.d0)*(-m_euler+(x-1.d0)*.5d0*zt2)
      elseif(abs(n+1.d0)+abs(x-2.d0) < xmeps)then
        f1=(n+1.d0)*(dpmg2+ddpmg2*.5d0*(n+1.d0)+
     $       dpzmg2*(x-2.d0))
     $       +(x-2.d0)*((1.d0-m_euler)+(x-1.d0)*.5d0*(zt2-1.d0))
c$$$      elseif(abs(x) > xmp)then
c$$$        u=-cgammai(-n)/x**(n+1.d0)
c$$$        f=cgpolygamma2(n,x)+dz0*u
c$$$        k=2.d0
c$$$        do
c$$$          u=-u*(n+k)/x
c$$$          f1=f+u*gpolygamma2(-k,1.d0)
c$$$          if(f1 == f)then
c$$$            exit
c$$$          endif
c$$$          f=f1
c$$$          k=k+1.d0
c$$$        enddo
      else
        if(ns /= n)then
          ns=n
          im1=0
          im2=0
          im3=0
          im4=0
          im5=0
          im6=0
          icg=1
          cgi(-2)=cgammai(2.d0-n)
          cgi(-1)=cgammai(1.d0-n)
          cgi(0)=cgammai(-n)
          cgi(1)=cgammai(-1.d0-n)
          chk1(0)=-m_euler-cpolygamma(1.d0-n)
          chk2(0)=cgi(-2)
          chk3(0)=chk2(0)
          chk4(0)=cgammai(2.d0+n)
          chk5(0)=chk4(0)
          chk6(0)=(-m_euler-cpolygamma(-n))
        endif
        f=(n*log(x)-m_euler*(x+n)-n*cpolygamma(-n))*cgi(-1)/x
     $       +cgi(-2)*zt2*x
        no=nogam*2+nopg*2
        k=1.d0
        f1=f
        do
          f10=f
          dfa=abs(f)**2*epso
          xk=zeroim(-x/k)
          xkk=xk/k
          axk=abs(xk)
c          write(*,'(a,1p10g12.4)')'cp2-k ',k,xk,n,f
c          if(.true.)then
c            f=f+(cgi(-2)-chg(cone,(2.d0,0.d0),
c     $           (2.d0,0.d0)-n,xk))*xkk
          if(abs(xk-xs1) < epsx .or. abs(xk-xs2) < epsx)then
            f=f+(cgi(-2)-chg1s(cone,(2.d0,0.d0),
     $           (2.d0,0.d0)-n,xk,.true.))*xkk
          elseif(dble(xk) <= 0.5d0)then
            if(axk .ge. 1.d0)then
              xk1=zeroim(1.d0/(1.d0-xk))
              xkk1=xkk*xk1
              lx=-log(xk1)
              u=-cgi(0)*xkk1*xk1
              f=f-cgi(-1)*xkk1-u*(lx+chk1(0))+cgi(-2)*xkk
              j=1.d0
              do jj=1,nmax
                j1=j+1.d0
                u=u*(j-n)/j*xk1
                if(jj > im1)then
                  chk1(jj)=polygamma(j1)-cpolygamma(j1-n)
                  im1=jj
                  no=no+nopg*2
                endif
                df=-u*(lx+chk1(jj))
                f=f+df
                if((abs(dble(df))+abs(imag(df)))**2
     $               <= no*dfa)then
                  no=no+16*3+4+jj*(3*3+1)
c                  write(*,'(a,1p10g12.4)')'pgm-1 ',xk,f1,f10
                  exit
                endif
                j=j1
              enddo
            elseif(abs(xk-1.d0) > 1.d0)then
c$$$              f=(one-x)**(-a)*chg(a,c-b,c,x/(x-1.d0))
              xk1=zeroim(xk/(xk-1.d0))
              xkk1=xk1/k
              f=f+cgi(-2)*(xkk1+xkk)
              xn=xk1*xkk1
              do i=1,nmax
                if(i > im2)then
                  chk2(i)=chk2(i-1)*((i-1)-n)/((i+1)-n)
                  im2=i
                endif
                df=chk2(i)*xn
                f=f+df
                if((abs(dble(df))+abs(imag(df)))**2
     $               <= no*dfa)then
c                  write(*,'(a,1p10g12.4)')'pgm-2 ',xk,f1,f10
                  no=no+6*3+i*(6+1)
                  exit
                endif
                xn=xn*xk1
              enddo
            else
              xn=-xk*xkk
              do i=1,nmax
                if(i > im3)then
                  chk3(i)=chk3(i-1)*(i+1)/((1+i)-n)
                  im3=i
                endif
                df=chk3(i)*xn
                f=f+df
                if((abs(dble(df))+abs(imag(df)))**2 <= no*dfa)then
                  no=no+3*3+i*(3+1)
c                  if(k <= 3.d0 .or. mod(k,100000.d0) == 0.d0)then
c                    write(*,'(a,1p10g12.4)')'pgm-3 ',k,xk,f10,f-f10
c                  endif
                  exit
                endif
                xn=xn*xk
              enddo
            endif
          else
c            if(abs(xk-1.d0) .ge. 1.d0)then
c              f=f+(cgi(-2)-chg(cone,(2.d0,0.d0),
c     $             (2.d0,0.d0)-n,xk))*xkk
c              write(*,'(a,1p10g12.4)')'pgmh-3 ',xk,f,f-f10
            if(abs(xk-1.d0) .ge. 1.d0)then
              xk1=1.d0/xk
              xkk1=xk1*xkk
              lx=log(zeroim(-xk))
              u=xk1*xkk1
              f=f+cgi(-1)*xkk1+u*(lx+chk6(0))*cgi(0)
     $             +cgi(-2)*xkk
              j=0.d0
              icg=0
              do jj=1,nmax
                j1=j+1.d0
                u=-u*xk1/j1
                j1n=-j1-n
                if(jj > icg)then
                  cgi(jj)=cgammai(j1n)
                  icg=jj
                  no=no+nocgam
                endif
                if(jj > im6)then
                  chk6(jj)=polygamma(j1+1.d0)-cpolygamma(j1n)
                  im6=jj
                  no=no+nopg*2
                endif
                if(imag(j1n) == 0.d0 .and.
     $               dble(j1n) == anint(dble(j1n)) .and.
     $               dble(j1n) <= 0.d0)then
                  df=-u*(-1.d0)**(1.d0-dble(j1n))*gamma(1.d0-dble(j1n))
                  no=no+nogam
                else
                  df=u*(lx+chk6(jj))*cgi(jj)
                endif
                f=f+df
                if((abs(dble(df))+abs(imag(df)))**2 <= no*dfa)then
                  no=no+19*3+jj*(4*3+2)
c                  write(*,'(a,i8,1p10g12.4)')'cpg-6 ',no,j1+n,xk,df,f
                  exit
                endif
                j=j1
              enddo
c            elseif(axk > 1.d0)then
c              f=f+(cgi(-2)-chg(cone,(2.d0,0.d0),
c     $             (2.d0,0.d0)-n,xk))*xkk
c              write(*,'(a,1p10g12.4)')'pgmh-5 ',xk,f,f-f10
            elseif(axk > 1.d0)then
              xk1=zeroim(1.d0-1.d0/xk)
              cab=-1.d0-n
              m=dble(cab)
              f=f+cgi(-2)*xkk
              if(imag(n) == 0.d0 .and. m == anint(m)
     $             .or. csinp(cab) == 0.d0)then
                if(m /= 0.d0)then
c                  u=-gamma(m)*gammai(m+2.d0)*gammai(m+1.d0)/k
                  u=-gammai(m+2.d0)/m/k
                  f=f+u
                  do i=1,int(m)-1
                    u=u*(3.d0+m-i)/(m-i)*xk1
                    f=f+u
                  enddo
                  no=no+nogam+int(m)*3
                endif
                lx=log(zeroim(-xk1))
                u=xk1**m*gammai(m+1.d0)/k
c                write(*,'(a,1p10g12.4)')'cpg2-n50 ',xk,m,u,f-f10
                f=f+u*((lx+1.d0)-xk1*(lx-1.d0))
                u=u*xk1
                j=2.d0
                no=no+nogam
c                write(*,'(a,1p10g12.4)')'cpg2-n51 ',xk,m,u,f-f10
                do
                  u=u/j*max(1.d0,j-2.d0)*xk1
                  f=f-u
                  if(abs(u)**2 <= no*dfa)then
                    no=no+20*3+7+int(j)*(3*3+1)
c                    write(*,'(a,1p10g12.4)')'cpg2-n5 ',xk,m,j,u,f1-f10
                    exit
                  endif
                  j=j+1.d0
                enddo
              else
                sabcp=m_pi/csinp(cab)
                xkk1=-cgi(-1)*cgi(0)/k*sabcp
c     f1=( chg(a,a-c+1.d0,1.d0-cab,x1)
c     $             *cgammai(c-a)*cgammai(c-b)*x**(-a)
c     $             -chg(c-a,1.d0-a,1.d0+cab,x1)
c     $             *(one-x)**cab*x**(a-c)
c     $             *cgammai(a)*cgammai(b))/sabc*m_pi
                f=f+chk5(0)*xkk1
     $               +cgi(0)*zeroim(1.d0-xk)**cab*xk**(n-1.d0)*sabcp*xkk
                xn=xk1*xkk1
                do i=1,nmax
                  if(i > im5)then
                    chk5(i)=chk5(i-1)*(n+i-1)/(i+1+n)
                    im5=i
                  endif
                  df=chk5(i)*xn
                  f=f+df
                  if((abs(dble(df))+abs(imag(df)))**2 <= no*dfa)then
c                    write(*,'(a,i8,1p10g12.4)')'pgm-5 ',no,xk,f,df
                    no=no+(23*3+2)+i*(3+1)
                    exit
                  endif
                  xn=xn*xk1
                enddo
              endif
c            elseif(.true.)then
c              f=f+(cgi(-2)-chg(cone,(2.d0,0.d0),
c     $             (2.d0,0.d0)-n,xk))*xkk
c              write(*,'(a,1p10g12.4)')'pgmh-4 ',xk,f,f-f10
            else
              xk1=zeroim(1.d0-xk)
              cab=-1.d0-n
              m=dble(cab)
              f=f+cgi(-2)*xkk
              sabc=csinp(cab)
              if(imag(n) == 0.d0 .and. m == anint(m)
     $             .or.  sabc == 0.d0)then
                if(m /= 0.d0)then
c                   u=-exp(-log_gamma(m+1.d0)-log_gamma(m+2.d0))*xkk
                  u=-gammai(1.d0+m)*gammai(2.d0+m)*xkk
                  f=f+u
                  do i=1,int(m)-1
                    u=-u*(1.d0+i)/(m-i)*xk1
                    f=f+u
                  enddo
                  no=no+nogam*2
                endif
                if(xk1 /= czero)then
                  lx=log(xk1)
                  u=zeroim(-xk1)**m*gammai(m+1.d0)*xkk
                  f=f+u*(lx+m_euler+polygamma(2.d0+m))
                  j=0.d0
                  no=no+nogam+nopg
                  do
                    j1=j+1.d0
                    u=u*(j1+m+1.d0)/j1*xk1
                    df=u*(lx-polygamma(j1+1.d0)+polygamma(2.d0+m+j1))
                    f=f+df
                    no=no+nopg*2
                    if((abs(dble(df))+abs(imag(df)))**2 <= no*dfa)then
c                      write(*,'(a,1p10g12.4)')'cpg2-n4 ',
c     $                     xk,m,j1,f1,f1-f10
                      exit
                    endif
                    j=j1
                  enddo
                endif
c              elseif(xk == cone)then
c                f1=f+dcmplx(1.d0/0.d0,0.d0)
c                return
              else
                if(xk1 == 0.d0 .and.
     $               (imag(cab) /= 0.d0 .or. dble(cab) < 0.d0))then
                  f1=dcmplx(-1.d0/0.d0,-1.d0/0.d0)
                  return
                endif
                sabcp=m_pi/sabc
c$$$  f1=(chg(a,b,1.d0-cab,x1)
c$$$  $             *cgammai(c-a)*cgammai(c-b)
c$$$  $             -chg(c-a,c-b,cab+1.d0,x1)*x1**cab
c$$$  $             *cgammai(a)*cgammai(b))/sabc*m_pi
                xkk1=-cgi(-1)*cgi(0)*xkk*sabcp
                f=f+chk4(0)*xkk1
                xn=xk1*xkk1
                do i=1,nmax
                  if(i > im4)then
                    chk4(i)=chk4(i-1)*(i+1)/(i+1+n)
                    im4=i
                  endif
                  u=chk4(i)*xn
                  f=f+u
                  if(abs(u)**2 <= no*dfa)then
                    no=no+25+i*4
                    exit
                  endif
                  xn=xn*xk1
                enddo
                u=cgi(0)*xk1**cab*xkk*sabcp
c                write(*,'(a,1p10g12.4)')'pgm-41 ',xk,xk1,cab,sabcp,u
                f=f+u
                j=0.d0
                do
                  j1=j+1.d0
                  u=u*(j1-n)/j1*xk1
                  f=f+u
                  if(abs(u)**2 <= no*dfa)then
                    no=no+(3*4+1)+int(j1)*(3*3+1)
c                    write(*,'(a,i8,1p10g12.4)')'pgm-4 ',no,xk,f-f10
                    exit
                  endif
                  j=j1
                enddo
              endif
            endif
          endif
          if(abs(f-f10)**2 <= no*dfa)then
            f1=f/x**n
            exit
          endif
          k=k+1.d0
        enddo
c      real*8 ,parameter ::eps=1.d-9
c$$$      elseif(dble(n) > 0.d0)then
c$$$        cpgz=x
c$$$        in=ceiling(dble(n))
c$$$        pgnc=in+2
c$$$        cpgn=pgnc-n-1.d0
c$$$        u=-cgammai(1.d0-n)/x**n
c$$$        f1=((cbint(cpgip,0.d0,1.d0,eps,eps)*x
c$$$     $       +gpolygamma2(pgnc-1.d0,1.d0))*x/cpgn
c$$$     $       +gpolygamma2(pgnc-2.d0,1.d0))/(cpgn-1.d0)
c$$$     $       *cgammai(cpgn-1.d0)*x**(cpgn-1.d0)
c$$$     $       -((log(x)-cpolygamma(-n)-m_euler)*n/x-m_euler)*u
c$$$c        write(*,'(a,1p10g12.4)')'cp2-p ',pgnc,f1,
c$$$c     $      ( cbint(pgip,0.d0,1.d0,1.d-2,1.d-11)*x
c$$$c     $       +gpolygamma2(pgnc,1.d0))/pgn
c$$$        do i=1,in-1
c$$$          u=-u*x*i/(i-n)
c$$$          f1=f1+zeta(i+1.d0)*u
c$$$        enddo
c$$$      elseif(dble(n) .ge. -2.d0)then
c$$$        cpgz=x
c$$$        cpgn=1.d0-n
c$$$        f1=(((cbint(cpgin1,0.d0,1.d0,eps,eps)*x
c$$$     $       +m_pi**2/6.d0)*x/cpgn
c$$$     $       -m_euler)/(-n)
c$$$     $       -(log(x)-cpolygamma(-n)-m_euler)/x)/x**n*cgammai(-n)
c$$$      else
c$$$        cpgz=x
c$$$        cpgn=-1.d0-n
c$$$        f1=(cbint(cpgin,0.d0,1.d0,eps,eps)
c$$$     $       -(log(x)-cpolygamma(-n)-m_euler)/x)/x**n*cgammai(-n)
      endif
      return
      end function

      real*8 recursive function polygamma2(n,x) result(f1)
      use tfstk, only:ktfenanq
      implicit none
      real*8 ,intent(in):: n,x
      integer*4 ,parameter :: nmax=1000
      real*8 ,parameter ::eps=1.d-9
      real*8 f,xk,xn,u,xk1,lx,f10,xkk,xkk1,df
      real*8 ,save :: ns=dcmplx(-veryl,0.d0),cgi(-2:nmax),
     $     chk1(0:nmax),chk2(0:nmax)
      real*8 k,j,j1,axk,dfa
      integer*4 ,save :: im1,im2
      integer*4 i,jj,icg,no
      if(x == 0.d0)then
        if(n .ge. -1.d0)then
          f1=1.d0/0.d0
        else
          f1=0.d0
        endif
      elseif(n == anint(n) .and. n .ge. 0.d0)then
        f1=gpolygamma2(n,x)
      elseif(abs(x) < xmth)then
        u=gammai(2.d0-n)*x
        f=(m_euler-log(x)+polygamma(-n))*gammai(-n)/x
     $       -m_euler*gammai(1.d0-n)+zt2*u
        no=nogam*2
        k=2.d0
        do
          u=-u*x*k/(k-n)
          df=zeta(k+1.d0)*u
          f=f+df
          if(df**2 <= no*f**2*epso)then
            f1=f/x**n
            exit
          endif
          k=k+1.d0
          no=no+nozt
        enddo
      elseif(abs(n+1.d0)+abs(x-1.d0) < xmeps)then
        f1=(n+1.d0)*(dpmg1+ddpmg1*.5d0*(n+1.d0)+
     $       dpzmg1*(x-1.d0))
     $       +(x-1.d0)*(-m_euler+(x-1.d0)*.5d0*zt2)
      elseif(abs(n+1.d0)+abs(x-2.d0) < xmeps)then
        f1=(n+1.d0)*(dpmg2+ddpmg2*.5d0*(n+1.d0)+
     $       dpzmg2*(x-2.d0))
     $       +(x-2.d0)*(1.d0-m_euler+(x-2.d0)*.5d0*(-1.d0+zt2))
c$$$      elseif(abs(x) > xmp)then
c$$$        u=-gammai(-n)/x**(n+1.d0)
c$$$        f=gpolygamma2(n,x)+dz0*u
c$$$        k=2.d0
c$$$        do
c$$$          u=-u*(n+k)/x
c$$$          f1=f+u*gpolygamma2(-k,1.d0)
c$$$          if(f1 == f)then
c$$$            exit
c$$$          endif
c$$$          f=f1
c$$$          k=k+1.d0
c$$$        enddo
      else
        if(ns /= n)then
          ns=n
          im1=0
          im2=0
          icg=1
          cgi(-2)=gammai(2.d0-n)
          cgi(-1)=gammai(1.d0-n)
          cgi(0)=gammai(-n)
          cgi(1)=gammai(-1.d0-n)
          chk1(0)=-m_euler-polygamma(1.d0-n)
          chk2(0)=cgi(-2)
        endif
        f=(n*log(x)-m_euler*(x+n)-n*polygamma(-n))*cgi(-1)/x
     $       +cgi(-2)*zt2*x
c        f=(m_euler-log(x)+polygamma(-n))*cgi(0)/x
c     $       -m_euler*cgi(-1)+cgi(-2)*zt2*x
        k=1.d0
        no=nogam*2
        do
          f10=f
          dfa=f**2*epso
          xk=-x/k
          xkk=xk/k
          axk=abs(xk)
          if(xk <= 0.5d0)then
            if(axk .ge. 1.d0)then
              xk1=1.d0/(1.d0-xk)
              xkk1=xkk*xk1
              lx=-log(xk1)
              u=-cgi(0)*xkk1*xk1
              f=f-cgi(-1)*xkk1-u*(lx+chk1(0))+cgi(-2)*xkk
              j=1.d0
              do jj=1,nmax
                j1=j+1.d0
                u=u*(j-n)/j*xk1
                if(jj > im1)then
                  chk1(jj)=polygamma(j1)-polygamma(j1-n)
                  im1=jj
                endif
                df=-u*(lx+chk1(jj))
                f=f+df
                if(df**2 <= no*dfa)then
                  no=no+16*1+4+jj*(3*1+1)
c                  write(*,'(a,i5,1p10g12.4)')'pgm-1 ',jj,xk,f10-f1
                  exit
                endif
                j=j1
              enddo
            elseif(abs(xk-1.d0) > 1.d0)then
c$$$              f=(one-x)**(-a)*chg(a,c-b,c,x/(x-1.d0))
              xk1=xk/(xk-1.d0)
              xkk1=xk1/k
              f=f+cgi(-2)*(xkk1+xkk)
              xn=xk1*xkk1
              do i=1,nmax
                if(i > im2)then
                  chk2(i)=chk2(i-1)*((i-1)-n)/((i+1)-n)
                  im2=i
                endif
                df=chk2(i)*xn
                f=f+df
                if(df**2 <= no*dfa)then
c                  write(*,'(a,1p10g12.4)')'pgm-2 ',xk,f1,f10
                  no=no+6*1+i*(6+1)
                  exit
                endif
                xn=xn*xk1
              enddo
            else
              go to 9000
            endif
          else
            go to 9000
          endif
          if(abs(f-f10)**2 <= no*dfa)then
            f1=f/x**n
            exit
          endif
          k=k+1.d0
        enddo
      endif
      return
 9000 write(*,'(a,1p8g15.7)')
     $     'polygamma2-implementation error ',n,x
      f1=0.d0
      return
      end function

      complex*16 function cgammaq(a,x) result(f1)
      implicit none
      complex*16 ,intent(in):: a,x
      if(a == czero .and. x == czero)then
        f1=czero
      else
        f1=cgammai(a)*cgamma2(a,x)
c        f1=1.d0-confhg1(a,a+1.d0,-x)*zeroim(x)**a
      endif
      return
      end

      real*8 function gammaq(a,x) result(f1)
      implicit none
      real*8 ,intent(in):: a,x
      if(a == 0.d0 .and. x == 0.d0)then
        f1=0.d0
      else
        f1=gammai(a)*gamma2(a,x)
c        f1=1.d0-confhgrr1(a,a+1.d0,-x)*x**a
      endif
      return
      end

      complex*16   function cgammap(a,x) result(f1)
      implicit none
      complex*16 ,intent(in):: a,x
      if(a == czero .and. x == czero)then
        f1=cone
      else
        f1=1.d0-cgammaq(a,x)
c        f1=confhg1(a,a+1.d0,-x)*zeroim(x)**a
      endif
      return
      end

      real*8   function gammap(a,x) result(f1)
      implicit none
      real*8  ,intent(in):: a,x
      if(a == 0.d0 .and. x == 0.d0)then
        f1=1.d0
      else
        f1=1.d0-gammaq(a,x)
      endif
      return
      end

      complex*16 recursive  function cgamma2(s,z) result(f1)
      use tfstk
      implicit none
      complex*16 ,intent(in):: s,z
      complex*16 lg0,lg1,lg2,u,df
      real*8 n,n1,adf0,adf
      integer*4 no,j,k
      real*8 ,parameter ::czlim=1.d-50
      if(abs(z) < czlim)then
        if(dble(s) > 0.d0)then
          f1=cgamma(s)
        else
          f1=1.d0/czero
        endif
        return
      endif
      lg0=1.d0
      lg1=1.d0-s+z
      main: do
        if(lg1 == czero)then
          exit
        endif
        u=exp(log(z)*s-z)
c     u=zeroim(z)**s*exp(-z)
        f1=u/lg0/lg1
        if(ktfenanq(dble(f1)+imag(f1)))then
          exit
        endif
        n=1.d0
        no=nolog
        adf0=veryl
        do k=1,kg2max
          df=0.d0
          do j=1,2
            n1=n+1.d0
            lg2=((z-s+n+n1)*lg1-(n-s)*lg0)/n1
            if(lg2 == czero)then
              exit main
            endif
            lg0=lg1
            lg1=lg2
            u=u*(n-s)/n1
            df=df+u/lg0/lg1
            n=n1
          enddo
          adf=abs(dble(df))+abs(imag(df))
c     write(*,'(a,1p10g12.4)')'cgm2-r ',n,s,z,df,f1
          if(adf > adf0 .or. ktfenanq(dble(f1)+imag(f1)))then
            exit main
          endif
          adf0=adf
          f1=f1+df
          no=no+50
          if((adf/abs(f1))**2 <= no*epso)then
            return
          endif
          n=n1
        enddo
        exit
      enddo main
c      write(*,'(a,1p10g12.4)')'cgm2-nc ',s,z,df,f1
      f1=exp(-z)*chguaa(1.d0-s,z,.false.)
c      write(*,'(a,1p10g12.4)')'cgm2-h ',s,z,f1
      return
      end function

      complex*16 recursive  function cgamma2r(s,z) result(f1)
      use tfstk
      implicit none
      complex*16 ,intent(in):: s
      real*8 ,intent(in):: z
      complex*16 lg0,lg1,lg2,u,df
      real*8 n,n1,adf0,adf
      integer*4 no,j,k
      real*8 ,parameter ::czlim=1.d-50
      if(abs(z) < czlim)then
        if(dble(s) > 0.d0)then
          f1=cgamma(s)
        else
          f1=1.d0/czero
        endif
        return
      endif
      lg0=1.d0
      lg1=1.d0-s+z
      main: do
        if(lg1 == czero)then
          exit
        endif
        u=exp(log(z)*s-z)
c     u=zeroim(z)**s*exp(-z)
        f1=u/lg0/lg1
        if(ktfenanq(dble(f1)+imag(f1)))then
          exit
        endif
        n=1.d0
        no=nolog
        adf0=veryl
        do k=1,kg2max
          df=0.d0
          do j=1,4
            n1=n+1.d0
            lg2=((z-s+n+n1)*lg1-(n-s)*lg0)/n1
            if(lg2 == czero)then
              exit main
            endif
            lg0=lg1
            lg1=lg2
            u=u*(n-s)/n1
            df=df+u/lg0/lg1
            n=n1
          enddo
          adf=abs(dble(df))+abs(imag(df))
c     write(*,'(a,1p10g12.4)')'cgm2-r ',n,s,z,df,f1
          if(adf > adf0 .or. ktfenanq(dble(f1)+imag(f1)))then
            exit main
          endif
          adf0=adf
          f1=f1+df
          no=no+100
          if((adf/abs(f1))**2 <= no*epso)then
            return
          endif
          n=n1
        enddo
        exit
      enddo main
c      write(*,'(a,1p10g12.4)')'cgm2-nc ',s,z,df,f1
      f1=exp(-z)*chguaa(1.d0-s,dcmplx(z,0.d0),.false.)
c      write(*,'(a,1p10g12.4)')'cgm2-h ',s,z,f1
      return
      end function

      real*8 recursive function gamma2(s,z) result(f1)
      implicit none
      real*8 ,intent(in):: s,z
      real*8 lg0,lg1,lg2,u,df,n,n1,adf0,adf
      integer*4 no,j,k
      if(z == 0.d0)then
        if(s > 0.d0)then
          f1=gamma(s)
        else
          f1=1.d0/0.d0
        endif
        return
      endif
      lg0=1.d0
      lg1=1.d0-s+z
      main: do
        if(lg1 == 0.d0)then
          exit
        endif
        u=exp(s*log(abs(z))-z)
        f1=u/lg0/lg1
        n=1.d0
        no=nolog*2
        adf0=veryl
        do k=1,kg2max
          df=0.d0
          do j=1,4
            n1=n+1.d0
            lg2=((z-s+n+n1)*lg1-(n-s)*lg0)/n1
            if(lg2 == 0.d0)then
              exit main
            endif
            lg0=lg1
            lg1=lg2
            u=u*(n-s)/n1
            df=df+u/lg0/lg1
            n=n1
          enddo
          adf=abs(df)
          if(adf > adf0)then
            exit main
          endif
          adf0=adf
          f1=f1+df
          no=no+50
          if((df/f1)**2 <= no*epso)then
c            write(*,'(a,1p10g12.4)')'gm2-c ',n,s,z,df,f1
            return
          endif
          n=n1
        enddo
        exit
      enddo main
      f1=exp(-z)*hguaa(1.d0-s,z,.false.)
c      write(*,'(a,1p10g12.4)')'gm2-h ',s,z,f1
      return
      end function

      complex*16 recursive  function cgamma3(a,x1,x2) result(f1)
      implicit none
      complex*16 ,intent(in):: a,x1,x2
      f1=cgamma2(a,x1)-cgamma2(a,x2)
c      f1=exp(-x1)*chguaa(1.d0-a,x1)-exp(-x2)*chguaa(1.d0-a,x2)
      return
      end function

      real*8 recursive   function gamma3(a,x1,x2) result(f1)
      implicit none
      real*8 ,intent(in):: a,x1,x2
      f1=gamma2(a,x1)-gamma2(a,x2)
c      f1=exp(-x1)*hguaa(1.d0-a,x1)-exp(-x2)*hguaa(1.d0-a,x2)
      return
      end function

      complex*16 recursive function chgu(a,b,x0) result(f1)
      implicit none
      complex*16 ,intent(in):: a,b,x0
      complex*16 lx,u,ab1,x,df
      real*8 k,k1,rb,ra,rab1
      integer*4 n,i,no
      if(a == b)then
        f1=chguaa(a,x0,.true.)
        return
      endif
      ra=dble(a)
      rb=dble(b)
      x=zeroim(x0)
      if(imag(a) == 0.d0 .and. imag(b) == 0.d0
     $     .and. imag(x) == 0.d0 .and. dble(x) .ge. 0.d0)then
        f1=dcmplx(hgu(ra,rb,dble(x)),0.d0)
        return
      endif
      ab1=a-b+1.d0
      rab1=dble(ab1)
      if(x == czero)then
        if(rb < 1.d0)then
          f1=cpochh(ab1,-a)
        else
          f1=1.d0/czero
        endif
      elseif(imag(a) == 0.d0 .and.
     $       ra == anint(ra) .and. ra <= 0.d0)then
        n=nint(-ra)
        u=x**n
        k=0.d0
        f1=u
        do i=1,n
          k1=k+1.d0
          u=u*(ra+k)/k1*(rb-ra-k1)/x
          f1=f1+u
          k=k1
        enddo
      elseif(imag(ab1) == 0.d0 .and.
     $       rab1 == anint(rab1) .and. rab1 <= 0.d0)then
c        f1=(-1.d0)**nint(-rab1)/x**(b-1.d0)*cpochh(2.d0-b,-ab1)
c     $       *confhg1(ab1,2.d0-b,x,.false.)
        f1=(-1.d0)**nint(-rab1)*exp(-(b-1.d0)*log(x)-cloggamma(2.d0-b)
     $       +cloggamma(1.d0-a))
     $       *confhg1(ab1,2.d0-b,x,.false.)
      elseif(imag(b) == 0.d0 .and. anint(rb) == rb)then
        n=nint(rb-1.d0)
        if(n <= -1.d0)then
          f1=x**(-n)*chgu(ab1,2.d0-b,x)
        else
          lx=log(x)
          if(n > 0)then
            u=cpochh(2.d0-a,b-2.d0)*cgammai(a)*gammai(rb-1.d0)/x
c            u=exp((1.d0-rb)*lx+log_gamma(rb-1.d0)-cloggamma(a))
            f1=u
c            write(*,'(a,1p10g12.4)')'chgu   ',a,x,u,
c     $           cpochh(2.d0-a,b-2.d0),f1
            k=1.d0
            do i=2,n
              k1=k+1.d0
              u=u*k*(rb-k1)/(k1-a)/x
              f1=f1+u
c              write(*,'(a,1p10g12.4)')'chgu ',a,b,x,u,f1
              k=k1
            enddo
          else
            f1=czero
          endif
c          u=-(-1.d0)**n*cgammai(a-rb+1.d0)*gammai(rb)
          u=-(-1.d0)**n*exp(-cloggamma(ab1)-log_gamma(rb))
          f1=f1+u*(lx+cpolygamma(a)+m_euler-polygamma(rb))
          k=0.d0
          no=nolog*2+nocgam*3+n*5
          do
            k1=k+1.d0
            u=u*(a+k)/k1/(rb+k)*x
            df=u*(lx+cpolygamma(a+k1)-polygamma(k1+1.d0)
     $           -polygamma(rb+k1))
            f1=f1+df
            no=no+nocgam*3
            if((abs(dble(df))+abs(imag(df)))**2
     $           <= no*abs(f1)**2*epso)then
              exit
            endif
            k=k1
          enddo
        endif
      else
        f1=(cgammai(ab1)*confhg1(a,b,x,.true.)
     $       -exp((1.d0-b)*log(x)-cloggamma(a))
     $       *confhg1(ab1,2.d0-b,x,.true.))*m_pi/csinp(b)
      endif
      return
      end function

      complex*16 recursive function chguaa(a,x,gm) result(f1)
      implicit none
      complex*16 ,intent(in):: a,x
      logical*4 ,intent(in):: gm
      complex*16 lx,u,ab1,df
      real*8 k,k1,ra,rab1
      integer*4 n,i,no
      ra=dble(a)
      if(x == czero)then
        if(dble(a) < 1.d0)then
          f1=cgamma(1.d0-a)
        else
          f1=1.d0/czero
        endif
        return
      endif
      if(imag(a) == 0.d0
     $     .and. imag(x) == 0.d0 .and. dble(x) .ge. 0.d0)then
        f1=dcmplx(hguaa(ra,dble(x),.true.),0.d0)
        return
      endif
      ab1=1.d0
      rab1=1.d0
      if(imag(a) == 0.d0 .and. ra == anint(ra))then
        if(ra <= 0.d0)then
          n=nint(-ra)
          u=x**n
          k=0.d0
          f1=u
          do i=1,n
            k1=k+1.d0
            u=-u*(ra+k)/x
            f1=f1+u
            k=k1
          enddo
        else
          n=nint(ra-1.d0)
          if(n <= -1.d0)then
            f1=x**(-n)*chgu(cone,2.d0-a,x)
          else
            lx=log(x)
            if(n > 0)then
              u=pochh(2.d0-ra,ra-2.d0)*gammai(ra)*gammai(ra-1.d0)/x
              f1=u
c              write(*,'(a,1p10g12.4)')'chguaa ',a,x,u,
c     $             pochh(2.d0-ra,ra-2.d0),f1
              k=1.d0
              do i=2,n
                u=-u*k/x
                f1=f1+u
                k=k+1.d0
              enddo
            else
              f1=czero
            endif
            u=-(-1.d0)**n*gammai(ra)
            f1=f1+u*(lx+m_euler)
            k=0.d0
            no=nolog*2+nocgam+n*3
            do
              k1=k+1.d0
              u=u/k1*x
              df=u*(lx-polygamma(k1+1.d0))
              f1=f1+df
              no=no+nocgam
              if((abs(dble(df))+abs(imag(df)))**2
     $             <= no*abs(f1)**2*epso)then
                exit
              endif
              k=k1
            enddo
          endif
        endif
      elseif(gm)then
        f1=exp(x)*cgamma2(1.d0-a,x)
      else
        f1=cgammai(a)*(exp(x)
     $       -x**(1.d0-a)*confhg1(cone,2.d0-a,x,.true.))*m_pi/csinp(a)
      endif
      return
      end function

      real*8 recursive function hgu(a,b,x) result(f1)
      implicit none
      real*8 ,intent(in):: a,b,x
      real*8 lx,u,k,k1,ab1,df
      integer*4 n,i,no
      ab1=a-b+1.d0
      if(x == 0.d0)then
        if(b < 1.d0)then
          f1=pochh(ab1,-a)
        else
          f1=1.d0/0.d0
        endif
      elseif(a == anint(a) .and. a <= 0.d0)then
        n=nint(-a)
        u=x**n
        k=0.d0
        f1=u
        do i=1,n
          k1=k+1.d0
          u=u*(a+k)/k1*(b-a-k1)/x
          f1=f1+u
          k=k1
        enddo
      elseif(ab1 == anint(ab1) .and. ab1 <= 0.d0)then
        f1=(-1.d0)**nint(-ab1)*pochh(2.d0-b,-ab1)*
     $       confhgrr1(ab1,2.d0-b,x,.false.)/x**(b-1.d0)
c        write(*,'(a,1p10g12.4)')'hgu-ab1 ',ab1,f1,
c     $       confhgrr1(ab1,2.d0-b,x,.false.),pochh(2.d0-b,-ab1)
      elseif(b == anint(b))then
        n=nint(b-1.d0)
        if(n .ge. 0)then
          if(n == 0)then
            f1=0.d0
            no=0
          else
c     u=exp(-log_gamma(a)+log_gamma(b-a)-log_gamma(2.d0-a)
c     $         -log_gamma(b-1.d0))/x
            u=pochh(2.d0-a,b-2.d0)*gammai(a)*gammai(b-1.d0)/x
            f1=u
            k=1.d0
            do i=2,n
              k1=k+1.d0
              u=u*k*(b-k1)/(k1-a)/x
              f1=f1+u
              k=k1
            enddo
            no=nogam*4+(n-1)*3
          endif
          lx=log(x)
c     u=-(-1.d0)**n*exp(-log_gamma(ab1)-log_gamma(b))
          u=-(-1.d0)**n*gammai(ab1)*gammai(b)
          f1=f1+u*(lx+polygamma(a)+m_euler-polygamma(b))
          k=0.d0
          no=no+nogam*2+nolog+4
          do
            k1=k+1.d0
            u=u*(a+k)/k1/(b+k)*x
            df=u*(lx+polygamma(a+k1)-polygamma(k1+1.d0)
     $           -polygamma(b+k1))
            f1=f1+df
            no=no+nogam*6
            if(df**2 <= no*f1**2*epso)then
              exit
            endif
            k=k1
          enddo
        else
          f1=x**(-n)*hgu(ab1,2.d0-b,x)
        endif
      else
        f1=(gammai(ab1)*confhgrr1(a,b,x,.true.)
     $       -x**(1.d0-b)*gammai(a)
     $       *confhgrr1(ab1,2.d0-b,x,.true.))*m_pi/sinp(b)
      endif
      return
      end function

      real*8 recursive function hguaa(a,x,gm) result(f1)
      use tfstk
      implicit none
      real*8 ,intent(in):: a,x
      logical*4 ,intent(in):: gm
      real*8 lx,u,k,k1,df
      integer*4 n,i,no
      if(x == 0.d0)then
        if(a < 1.d0)then
          f1=gamma(1.d0-a)
        else
          f1=1.d0/0.d0
        endif
      elseif(a == anint(a))then
        if(a <= 0.d0)then
          n=nint(-a)
          u=x**n
          k=0.d0
          f1=u
          do i=1,n
            k1=k+1.d0
            u=-u*(a+k)/x
            f1=f1+u
            k=k1
          enddo
        else
          n=nint(a-1.d0)
          if(n .ge. 0)then
            if(n == 0)then
              f1=0.d0
              no=0
            else
c     u=exp(-log_gamma(a)+log_gamma(b-a)-log_gamma(2.d0-a)
c     $         -log_gamma(b-1.d0))/x
              u=pochh(2.d0-a,a-2.d0)*gammai(a)*gammai(a-1.d0)/x
              f1=u
              k=1.d0
              do i=2,n
                k1=k+1.d0
                u=-u*k/x
                f1=f1+u
                k=k1
              enddo
              no=nogam*4+(n-1)*3
            endif
            lx=log(x)
            no=no+nolog
c     u=-(-1.d0)**n*exp(-log_gamma(ab1)-log_gamma(b))
c            u=-(-1.d0)**n*gammai(a)
            u=merge(-1.d0,1.d0,mod(n,2)==0)*gammai(a)
            f1=f1+u*(lx+m_euler)
            k=0.d0
            do
              k1=k+1.d0
              u=u/k1*x
              df=u*(lx-polygamma(k1+1.d0))
              f1=f1+df
              no=no+nopg+5
              if((df/f1)**2 <= no*epso .or. ktfenanq(f1))then
                exit
              endif
              k=k1
            enddo
          else
            f1=x**(-n)*hgu(1.d0,2.d0-a,x)
          endif
        endif
      elseif(gm)then
        f1=exp(x)*gamma2(1.d0-a,x)
      else
        f1=gammai(a)*(exp(x)-x**(1.d0-a)
     $       *confhgrr1(1.d0,2.d0-a,x,.true.))*m_pi/sinp(a)
c        write(*,'(a,1p10g12.4)')'hguaa ',a,x,
c     $       confhgrr1(1.d0,2.d0-a,x,.true.)
      endif
      return
      end function

      complex*16 function cerf(z) result(f1)
      implicit none
      complex*16 ,intent(in):: z
      f1=2.d0*z/m_sqrtpi
     $     *confhg1((.5d0,0.d0),(1.5d0,0.d0),zeroim(-z**2),.false.)
      return
      end function

      complex*16 function cerfc(z) result(f1)
      implicit none
      complex*16 ,intent(in):: z
      f1=1.d0-2.d0*z/m_sqrtpi*confhg1((.5d0,0.d0),(1.5d0,0.d0),
     $     zeroim(-z**2),.false.)
      return
      end function

      complex*16 recursive function cinverseerf(x) result(f1)
      implicit none
      complex*16 ,intent(in):: x
      complex*16 df,f,x2
      real*8, parameter:: c=sqrt(m_pi_2),eps=5.d-17**2,
     $     dimax=2.d0,drmax=0.5d0,dmin=0.1d0
      integer*4 ,parameter :: imax=1000
      real*8 adf,adf0,fact
      integer*4 no,i
      if(dble(x) < 0.d0)then
        f1=-cinverseerf(zeroim(-x))
c        write(*,'(a,1p10g12.4)')'cierf-1 ',x,f1
        return
      elseif(sign(1.d0,imag(x)) < 0.d0)then
        f1=conjz(cinverseerf(conjz(x)))
c        write(*,'(a,1p10g12.4)')'cierf-2 ',x,f1
        return
      endif
      If(abs(x-1.d0) + abs(x-1.8d0) < 1.2d0)then
        f1=log(sqrt(2.d0/m_pi)/(1.d0-x))
        f1=sqrt(zeroim(f1-.5d0*log(2.d0*f1)))
        no=nolog*2
      elseif(abs(x) < 1.d0)then
        x2=m_pi*x**2
        f1=.5d0*m_sqrtpi*x*(1.d0+x2*(1.d0/12.d0
     $       +x2*(7.d0/480.d0+x2*(127.d0/40320.d0
     $       +x2*4369.d0/5806080.d0))))
        no=20+nolog
      else
        f1=sqrt(-4.d0/m_pi*log(0.5d0-(0.d0,1.d0)*x))
        no=nolog+10
      endif
      f=f1
      adf0=veryl
      fact=1.d0
      do i=1,imax
        df=(x-cerf(f))*c*exp(f**2)
        adf=(abs(dble(df))+abs(imag(df)))**2
        fact=max(1.d0,sqrt(adf/adf0))
        if(fact > 1.d4)then
          exit
        endif
        f1=f+df/fact
c        if(fact /= 1.d0)then
c          write(*,'(a,1p10g12.4)')'cierf-i ',x,fact,f1,df
c        endif
        if(adf/fact**2 <= no*abs(f1)**2*eps)then
          exit
        endif
        f=f1
        adf0=adf*4.d0
        no=no+40
      enddo
c      write(*,'(a,1p10g12.4)')'cierf-3 ',x,f1
      return
      end

      real*8  recursive function inverseerf(x) result(f1)
      implicit none
      real*8 ,intent(in):: x
      real*8 df,x2
      integer*4 no,i
      real*8, parameter:: c=sqrt(m_pi_2),eps=5.d-17
      integer*4 ,parameter :: imax=1000
      if(x < 0.d0)then
        f1=-inverseerf(-x)
        return
      endif
      if(abs(x) == 1.d0)then
        f1=x/0.d0
        return
      elseif(x > 0.8d0)then
        f1=log(sqrt(2.d0/m_pi)/(1.d0-x))
        f1=sqrt(f1-.5d0*log(2.d0*f1))
        no=nolog*2+10
      else
        x2=m_pi*x**2
        f1=.5d0*m_sqrtpi*x*(1.d0+x2*(1.d0/12.d0
     $       +x2*(7.d0/480.d0+x2*(127.d0/40320.d0
     $       +x2*4369.d0/5806080.d0))))
        no=11
      endif
      do i=1,imax
        df=(x-erf(f1))*c*exp(f1**2)
        f1=f1+df
        if(df**2 <= no*f1**2*eps)then
          exit
        endif
        no=no+nolog
      enddo
      return
      end

      complex*16 function cbesj(n,z) result(f)
      implicit none
      complex*16 ,intent(in):: n,z
      f=zeroim(.5d0*z)**n*confhg0(n+1.d0,-.25d0*z**2)
      return
      end

      complex*16 function cbesi(n,z) result(f)
      implicit none
      complex*16 ,intent(in):: n,z
      f=zeroim(.5d0*z)**n*confhg0(n+1.d0,.25d0*z**2)
      return
      end

      complex*16 function cbesy(n,z0) result(f)
      implicit none
      complex*16 ,intent(in):: n,z0
      complex*16 z
      z=zeroim(z0)
      f=((.5d0*z)**n*ccosp(n)*confhg0(n+1.d0,-.25d0*z**2)
     $     -(2.d0/z)**n*confhg0(1.d0-n,-.25d0*z**2))/csinp(n)
      return
      end

      complex*16 function cbesk(n,z0) result(f)
      implicit none
      complex*16 ,intent(in):: n,z0
      complex*16 z
      z=zeroim(z0)
      f=.5d0*m_pi/csinp(n)*((2.d0/z)**n*confhg0(1.d0-n,.25d0*z**2)
     $     -(.5d0*z)**n*confhg0(1.d0+n,.25d0*z**2))
      return
      end
      
      real*8 function gamma0(x) result(f1)
      implicit none
c     Including Euler's gamma(euler)
      real*8 ,intent(in):: x
      if(x == 0.d0)then
        f1=1.d0/0.d0
      elseif(x < 1.d0)then
        f1=gamma0ser(x)-log(x)-m_euler
      else
        f1=gamma0cf(x)
      endif
      return
      end

      complex*16 function cgamma0(x) result(f1)
      implicit none
c     Including Euler's gamma(euler)
      complex*16 ,intent(in):: x
      real*8 ,parameter :: xth=10.d0
      if(imag(x) == 0.d0 .and. dble(x) >= 0.d0)then
        f1=dcmplx(gamma0(dble(x)),0.d0)
      elseif(abs(x) < xth)then
        f1=cgamma0ser(x)-log(x)-m_euler
      else
        f1=cgamma0cf(x)
        if(imag(x) == 0.d0)then
          f1=f1+dcmplx(0.d0,-m_pi)
        endif
      endif
      return
      end

      real*8 function gamma0log(x)
      implicit none
c     Including Euler's gamma(euler)
      real*8 ,intent(in):: x
      if(x < 0.d0)then
        gamma0log=x/0.d0
        return
      endif
      if(x < 1.d0)then
        gamma0log=gamma0ser(x)
      else
        gamma0log=gamma0cf(x)+log(x)+euler
      endif
      return
      end

      real*8 function gammaser(a,x)
      implicit none
      integer*4 i
      real*8 ,intent(in):: a,x
      real*8 gln,ap,sum,del
      integer*4 ,parameter::itmax=300
      real*8 ,parameter ::eps=1.d-13
      if(x < 0.d0)then
        gammaser=x/0.d0
        return
      elseif(x == 0.d0)then
        gammaser=0.d0
        return
      elseif(a == 0.d0)then
        gammaser=1.d0
        return
      endif
      ap=a
      sum=1.d0
      del=sum
      do i=1,itmax
        ap=ap+1.d0
        del=del*x/ap
        sum=sum+del
        if(abs(del) < abs(sum)*eps)then
          go to 1
        endif
      enddo
 1    gln=log_gamma(a)
      gammaser=sum*exp(-x+a*log(x)-gln)/a
      return
      end
      
      real*8 function gammacf(a,x) result(f1)
      implicit none
      real*8 ,intent(in):: a,x
      integer*4 i
      real*8 a1,gln,b,c,d,h,an,del
      integer*4 ,parameter::itmax=300
      real*8 ,parameter ::eps=1.d-13,fpmin=1.d-30
      if(a == 0.d0)then
        f1=0.d0
        return
      endif
      a1=a-1.d0
      b=x-a1
      c=1.d0/fpmin
      d=1.d0/b
      h=d
      do i=1,itmax
        an=-i*(i-a)
        b=b+2.d0
        d=an*d+b
        if(abs(d) < fpmin)then
          d=fpmin
        endif
        c=b+an/c
        if(abs(c) < fpmin)then
          c=fpmin
        endif
        d=1.d0/d
        del=d*c
        h=h*del
        if(abs(del-1.d0) < eps)then
          exit
        endif
      enddo
      gln=log_gamma(a1+1.d0)
      f1=exp(-x+a*log(x)-gln)*h
      return
      end

      complex*16 function cgamma0ser(x) result(f1)
      implicit none
      complex*16 ,intent(in):: x
      integer*4 i
      complex*16 sum,del
      real*8 ap
      integer*4 ,parameter ::itmax=300
      real*8 ,parameter ::eps=1.d-13
      ap=1.d0
      sum=x
      del=x
      do i=1,itmax
        del=-del*x*ap
        ap=ap+1.d0
        del=del/ap**2
        sum=sum+del
        if(abs(del) < abs(sum)*eps)then
          exit
        endif
      enddo
      f1=sum
      return
      end
      
      real*8 function gamma0ser(x)
      implicit none
      real*8 ,intent(in):: x
      integer*4 i
      real*8 ap,sum,del
      integer*4 ,parameter ::itmax=300
      real*8 ,parameter ::eps=1.d-13
      if(x <= 0.d0)then
        gamma0ser=x/0.d0
        return
      endif
      ap=1.d0
      sum=x
      del=x
      do i=1,itmax
        del=-del*x*ap
        ap=ap+1.d0
        del=del/ap**2
        sum=sum+del
        if(abs(del) < abs(sum)*eps)then
          go to 1
        endif
      enddo
 1    gamma0ser=sum
      return
      end
      
      real*8 function gamma0cf(x)
      implicit none
      integer*4 i
      real*8 ,intent(in):: x
      real*8 a1,b,c,d,h,an,del
      integer*4 ,parameter ::itmax=300
      real*8 ,parameter ::eps=1.d-13,fpmin=1.d-30
      a1=-1.d0
      b=x-a1
      c=1.d0/fpmin
      d=1.d0/b
      h=d
      do i=1,itmax
        an=-i**2
        b=b+2.d0
        d=an*d+b
        if(abs(d) < fpmin)then
          d=fpmin
        endif
        c=b+an/c
        if(abs(c) < fpmin)then
          c=fpmin
        endif
        d=1.d0/d
        del=d*c
        h=h*del
        if(abs(del-1.d0) < eps)then
          go to 1
        endif
      enddo
 1    gamma0cf=exp(-x)*h
      return
      end

      complex*16 function cgamma0cf(x) result(f1)
      implicit none
      integer*4 i
      complex*16 ,intent(in):: x
      complex*16 b,c,d,h,del
      real*8 an
      integer*4 ,parameter ::itmax=4000
      real*8 ,parameter ::eps=1.d-13,fpmin=1.d-30
      b=x+1.d0
      c=1.d0/fpmin
      d=1.d0/b
      h=d
      do i=1,itmax
        an=-i**2
        b=b+2.d0
        d=an*d+b
        if(abs(d) < fpmin)then
          d=fpmin
        endif
        c=b+an/c
        if(abs(c) < fpmin)then
          c=fpmin
        endif
        d=1.d0/d
        del=d*c
        h=h*del
c        write(*,'(a,1p10g12.4)')'cg0cf ',x,del,h
        if(abs(del-1.d0) < eps)then
          exit
        endif
      enddo
      f1=exp(-x)*h
      return
      end

      complex*16 function cerfs(z)
      implicit none
      complex*16 ,intent(in):: z
c     Including m_2_sqrtpi:	2 / Sqrt[Pi]
      real*8 ,parameter ::r=m_2_sqrtpi
      complex*16 z2
      z2=z**2
      cerfs=r*z*(1.d0-z2*(1.d0/3.d0-z2*(1.d0/10.d0-
     $     z2*(1.d0/42.d0-z2*(1.d0/216.d0-z2/1320.d0)))))
      return
      end

      complex*16 function cerfcd(z)
      implicit none
c     Including m_2_sqrtpi:	2 / Pi
      complex*16 ,intent(in):: z
      real*8 an1,an2,an0
      real*8 ,parameter ::r=m_2_pi,eps=5.d-17
      complex*16 cs,cs1,ca,ca2,cd,z1,z2
      integer*4 i,no
      no=0
      if(imag(z) > 0.d0)then
        z1=dcmplx(imag(z),dble(z))
        an0=anint(dble(z1)/2.d0/hg)*2.d0
        z2=z1-an0*hg
        an1=an0+1.d0
        an2=an0-1.d0
        ca=exp(2.d0*hg*z2)
        ca2=ca**2
        cs=ag(1)*(ca/an1+1.d0/an2/ca)
        do i=2,itmaxg
          an1=an1+2.d0
          an2=an2-2.d0
          ca=ca*ca2
          cd=ag(i)*(ca/an1+1.d0/an2/ca)
          cs1=cs+cd
          if(abs(cd)**2 <= no*abs(cs1)**2*eps)then
            exit
          endif
          cs=cs1
          no=no+9
        enddo
        cs=exp(z1**2-z2**2)*cs
        cerfcd=dcmplx(1.d0-r*imag(cs),-r*dble(cs))
      else
        z1=dcmplx(-imag(z),dble(z))
        an0=anint(dble(z1)/2.d0/hg)*2.d0
        z2=z1-an0*hg
        an1=an0+1.d0
        an2=an0-1.d0
        ca=exp(2.d0*hg*z2)
        ca2=ca**2
        cs=ag(1)*(ca/an1+1.d0/an2/ca)
        do i=2,itmaxg
          an1=an1+2.d0
          an2=an2-2.d0
          ca=ca*ca2
          cd=ag(i)*(ca/an1+1.d0/an2/ca)
          cs1=cs+cd
          if(abs(cd)**2 <= no*abs(cs1)**2*eps)then
            exit
          endif
          cs=cs1
          no=no+9
        enddo
c        write(*,*)'cerfcd convergence error'
        cs=exp(z1**2-z2**2)*cs
        cerfcd=dcmplx(1.d0-r*imag(cs),r*dble(cs))
      endif
      return
      end

      complex*16 function cerfcf(z)
      implicit none
c     Including m_2_sqrtpi:	2 / Sqrt[Pi]
      complex*16 ,intent(in):: z
      real*8 ,parameter:: fpmin=1.d-300,eps=1.d-15,r=m_2_sqrtpi
      integer*4 i
      integer*4 ,parameter::itmax=1000
      real*8 an,a
      complex*16 cb,cd,ch,cdel,cc,z2
      z2=z**2
      cb=2.d0*z2+1.d0
      cc=1.d0/fpmin
      cd=1.d0/cb
      ch=cd
      an=-1.d0
      do i=1,itmax
        an=an+2.d0
        a=-an*(an+1.d0)
        cb=cb+4.d0
        cd=1.d0/(a*cd+cb)
        cc=cb+a/cc
        cdel=cc*cd
        ch=ch*cdel
        if(abs(cdel-1.d0) < eps)then
          go to 2
        endif
      enddo
c      write(*,*)'erfc convergence error'
 2    cerfcf=ch*z*r*exp(-z2)
      return
      end

      end module
