/*************************************************************
*                                                            *
* Differential Algebra Library at KEK.                       *
*            C++ version 2            93/06/15               *
*                                     K.Ohmi                 *
*************************************************************/
#define Nresist 6
#define N_var_max 100      //
#define L_poly_max 256     //    affect only polynomial expression
#define N_poly_max 100     //
#define N_dim_max 30       // affect only asin table
#include <iostream>
using std::istream;
using std::ostream;
using std::cout;
using std::ios;
using std::endl;
#include <cstdio>
using std::sscanf;
#include <cstdlib>
using std::atof;
using std::atoi;
using std::exit;
#include <cstring>
using std::strchr;
using std::strcmp;
using std::strcpy;
using std::strcspn;
using std::strlen;
using std::strspn;
using std::strstr;
using std::strtok;
using std::memcpy;
using std::memset;
#include <cctype>
using std::isxdigit;
#include <cmath>
using std::fabs;
using std::sqrt;
using std::pow;
using std::sin;
using std::cos;
using std::asin;
using std::atan;
using std::sinh;
using std::cosh;
using std::exp;
using std::log;

#include <dacpp.h>

void addr_tbl(void);
int poly_ana(char*,int*,double*);
void mon_ana(char,char*,int*,double*);

void p_da_init(int,int,const char*[]);
void p_da_init(int,int);

int Nvar,Nord,Lvec,Ncv;
int *d1,*d2,*Id1,*NMS,*NM1,*ixref,*iv1,*iv2; 
int Nvr2,Naddr,MaxC,Lvec2;
static double *z0,*z1,*z2,*z3,*z4,*z5;

static const char *var_name[N_var_max];
// constructor
da::da(void) { v=new double[L_vec];}
//
// copy constructor
//
da::da(const da& x)
{
   v=new double[Lvec];
   for(int i=0;i<Lvec;i++) v[i]=x.v[i];
}

 

 
void da_init(char* s,int No,int Npo=0)
{
  const char *sep=" ,";
  char *p;
  int i,ialoc;
  int Npvar=0;

  p=strtok(s,sep);
  var_name[0]=p;
  i=1;
  Nvar=0;
  while((p=strtok(NULL,sep))!=NULL) {
     if(*p!='&') {
        if(i<N_var_max) var_name[i]=p;
        i++;
     }
     else{
	Ncv=i;
     }
     if(!strcmp(p,"&&")) {
        if(i%2==1) {
          cout << "Number of variable is odd. Added 1 dummy variable @ \n";
          var_name[i]="@";
	  i++;
	}
	Nvar=i;
     }
  }
  
  if(Nvar==0) {
     if(i%2==1) {
       cout << "Number of variable is odd. Added 1 dummy variable @ \n";
       var_name[i]="@";
       i++;
     }
     Nvar=i;
     Npvar=0;
  }
  else{
     if(i%2==1) {
       cout <<"Number of variable is odd. Added 1 dummy variable @@ \n";
       var_name[i]="@@";
       i++;
     }
     Npvar=i-Nvar;
  }
  int Nvtot=i;
  cout << " Number of variables  = " << Ncv << " (cvar)  " 
       << Nvar << " (var)  " << Npvar << " (pvar)  " 
	<< Nvtot << " (tot)\n ";
  
  if(Nvar>=N_var_max) 
     cout << " variable name can not be stored in var_name[]\n";
  cout << " variables   ";
  for(i=0;i<Nvar;i++) cout << var_name[i];
  cout << '\n';
  
  Nord=No;
  Nvr2=Nvar/2;
  Naddr=1;  for(i=0;i<Nvr2;i++) Naddr=Naddr*(Nord+1);
  MaxC=Nord; for(i=0;i<Nvr2-1;i++) MaxC=MaxC*(Nord+1);
  Lvec=icomb(Nord+Nvar,Nvar);
  Lvec2=icomb(Nord+Nvr2,Nvr2);
  N_ord=Nord;
  N_var=Nvar;
  L_vec=Lvec;
  N_cv=Ncv;

  ialoc=2*(Naddr+1)+Lvec2+2*(Nord+1)+2*Nvar;

  d1=new int[ialoc];
  cout << "*** " << ialoc << "*" << sizeof(int) 
	<< "  Byte allocation for DA table \n";
  memset(d1,0,ialoc*sizeof(int));
/*	Make address table	*/
  addr_tbl();

  z0=new double[Nresist*Lvec];
  z1=z0+Lvec;  z2=z1+Lvec;  z3=z2+Lvec;  z4=z3+Lvec;  z5=z4+Lvec;
  cout << "*** " << Nresist << '*' << sizeof(double) 
	<< "   Byte allocation for DA resistor \n";
  
  //p_da init
  if(Npo>0) p_da_init(Npvar,Npo,var_name+Nvar);

}

int poly_ana(const char* s_in,int* Iadr,double* Coef)
{

   const char *sep="+-";
   char *p,s[L_poly_max],op,op_f='+';
   int iadr;
   double coef;

   int i=0,n,m,ls,lsin;
   if((lsin=strlen(s_in))>L_poly_max) {
     cout << " Input polynomial is too long : Max " << L_poly_max
	<< '<' << lsin <<'\n';
     exit(1);
   }
   for(i=0;i<lsin;i++) {
     if(strchr("(){}",*(s_in+i))!=NULL) {
        cout << " Donot use (){} to express polynomial : " << s_in 
	<< '\n';
        exit(1);
     }
   }
   strcpy(s,s_in);  
   ls=strlen(s);

   p=s;
   i=0;
   while(1) {
      n=strcspn(p,sep);
      op=*(p+n);
      *(p+n)='\0';
      if((m=strspn(p," "))<n) {
	 p+=m;
	 mon_ana(op_f,p,&iadr,&coef);
	 *(Iadr+i)=iadr;
	 *(Coef+i)=coef;
	 i++;
         p+=(n-m+1);
	 if(i>=N_poly_max) {
	    cout << "  No. of monomials is too many : Max = " 
		<< N_poly_max << '\n';
	    exit(1);
	 }
         if(p>=s+ls) break;
      } else p+=(n+1);
      op_f=op;
   }
   return i;

   /* term_ana(i,operator,term);*/
}
void mon_ana(char sign,char* term,int* iadr,double* coef)
{
   int *iv;
   char op='*',op_f=' ';
   int n,m,nop,lt,i,lv,jmem=0,hit=0;
   char *p,*var;

   iv=iv1;
   memset(iv,0,sizeof(int)*Nvar);

   if(sign=='+') *coef=1.;
   if(sign=='-') *coef=-1.;
   
   p=term;
   lt=strlen(term);

   while(term+lt>p) {
     n=strcspn(p," *^");
     nop=strspn(p+n," *^");
     op=' ';
     for(i=0;i<nop;i++) {
       if(*(p+n+i)=='*') {op='*'; break;}
       if(*(p+n+i)=='^') {op='^'; break;}
       op='*';
       if((p+n+i)==term+lt-1) op=' ';
     }
     *(p+n)='\0';
     m=strspn(p," ");
     var=p+m;
     lv=strcspn(var," ");
     *(var+lv)='\0';
     //   cout<< var << op;
     if(op_f=='^') iv[jmem]+=(atoi(p)-1); 
     else {
       for(int j=0;j<Nvar;j++) {
         if(!strcmp(var,var_name[j])) {
	   hit=1;
	   jmem=j;
	   iv[jmem]++;
	   break;
         }
       }
       op_f=op;
       if(!hit) {
         if(strcspn(var,"0123456789.E ")) {
 	   cout << "Input polynomial is  not matched : " << var << '\n';
	   exit(1);
         }
         *coef=*coef*atof(var);
       }
     }
     p+=(n+nop);
   }
   int C1M=ccal(Nord,Nvr2,iv1);
   int C2M=ccal(Nord,Nvr2,iv2);
   *iadr=*(d1+C1M)+*(d2+C2M);
}   

void InitializeDifferentialAlgebra(int No,int Nv,int Npo,int Npv)
{
  int i,ialoc;
  int Npvar=0;

  if(Nv%2==0) Nvar=Nv; 
  else {
     cout <<"Number of variable is odd. Added 1 dummy variable \n";
     Nvar=Nv+1;}
  if(Npv%2==0) Npvar=Npv; 
  else {
     cout <<"Number of external variable is odd. Added 1 dummy variable \n";
     Npvar=Npv+1;}

  int Nvtot=Nvar+Npvar;
  cout << " Number of variables  = " << Ncv << " (cvar)  " 
       << Nvar << " (var)  " << Npvar << " (pvar)  " 
	<< Nvtot << " (tot)\n ";
    
  Nord=No;
  Nvr2=Nvar/2;
  Ncv=Nvar;
  Naddr=1;  for(i=0;i<Nvr2;i++) Naddr=Naddr*(Nord+1);
  MaxC=Nord; for(i=0;i<Nvr2-1;i++) MaxC=MaxC*(Nord+1);
  Lvec=icomb(Nord+Nvar,Nvar);
  Lvec2=icomb(Nord+Nvr2,Nvr2);
  N_ord=Nord;
  N_var=Nvar;
  L_vec=Lvec;
  N_cv=Ncv;

  ialoc=2*(Naddr+1)+Lvec2+2*(Nord+1)+2*Nvar;

  d1=new int[ialoc];
  cout << "*** " << ialoc << "*" << sizeof(int) 
	<< "  Byte allocation for DA table \n";
  memset(d1,0,ialoc*sizeof(int));
/*	Make address table	*/
  addr_tbl();

  z0=new double[Nresist*Lvec];
  z1=z0+Lvec;  z2=z1+Lvec;  z3=z2+Lvec;  z4=z3+Lvec;  z5=z4+Lvec;
  cout << "*** " << Nresist << '*' << sizeof(double) 
	<< "   Byte allocation for DA resistor \n";
  
  //p_da init
  if(Npo>0) p_da_init(Npvar,Npo);
  cout.flush();
}
void addr_tbl(void)
{
/*
**************************************************************
*                                                            *
* Differential Algebra Library at KEK.                       *
*            sad version 1            92/07/19               *
*                                     K.Ohmi                 *
**************************************************************
*/

  int C1,D10,D20;
  int Nc,II,Nth,ICP;
  int i,j;

  d2=d1+Naddr+1;
  Id1=d2+Naddr+1;
  NMS=Id1+Lvec2;
  NM1=NMS+Nord+1;
  ixref=NM1+Nord+1;
  iv1=ixref+Nvar;
  iv2=iv1+Nvr2;

  //int* iv1=new int[Nvr2];

  cout << " \n \n Differential Algebra Library \n";
  cout << " Order and number of variables = " 
	<< Nord << "   " << Nvar << '\n';
  cout << " Number of canonical variables = " << Ncv << '\n';
  cout << " Total address is " << Naddr
	<< " :  Length of DA vector = " << Lvec << '\n';
  Nvr2=Nvar/2;

  for(i=0;i<Naddr;i++) {
    ICP=Naddr-1-i;
    ivcalp(ICP,Nord,Nvr2,iv1);

    Nc=0;
    for(j=0;j<Nvr2;j++) Nc+=*(iv1+j);

    if(Nc<=Nord) {

      C1=ccal(Nord,Nvr2,iv1);

      *(NM1+Nc)=(*(NM1+Nc))+1;
      if(C1<=MaxC) *(d1+C1)=*(NM1+Nc)-1;
/*      if(C1<=MaxC) *(d1+C1)=*(NM1+Nc); fortran version */
/*      cout.width(4);
	cout<< " Nc, iv1   ";
	cout.width(4); cout << Nc ;
        cout.width(4); cout << *iv1;
	cout.width(4); cout << *(iv1+1) << '\n';*/
    }
  }

  *NMS=1;
  //cout << "\n\n    Order      NM1       NMS  \n";
  for(i=0;i<Nord;i++) {
    *(NMS+i+1)=(*(NMS+i))+(*(NM1+i+1));
    //cout.width(10);
    //cout << i+1;
    //cout.width(10);
    //cout << *(NM1+i+1);
    //cout.width(10);
    //cout << *(NMS+i+1) << '\n';
  }


  if(Nord<4 && Nvar==4) 
    cout<<"\n        iv(1)     iv(2)       c1         d1        d2 \n";

  for(i=0;i<Naddr;i++) {
    ICP=Naddr-i-1;
    ivcalp(ICP,Nord,Nvr2,iv1);

    Nc=0;
    for(j=0;j<Nvr2;j++) Nc+=*(iv1+j);

    if(Nc<=Nord) {

      C1=ccal(Nord,Nvr2,iv1);

      D10=0;
      if(Nc>0) D10=*(NMS+Nc-1);
      D20=0;
      for(j=0;j<Nc;j++) D20+=*(NM1+j)**(NMS+Nord-j);

      if(C1<=MaxC) {

        *(d2+C1)=(*(d1+C1))**(NMS+Nord-Nc)+D20;
        *(d1+C1)=*(d1+C1)+D10;

/*c err FIXED 12/DEC/91 */

        II=*(d1+C1);

	if(II<=*(NMS+Nord)) {
	  *(Id1+II)=C1;
	}

	if(Nord<=4 && Nvar==4) {
	  cout.width(10); cout << *iv1 ;
	  cout.width(10); cout << *(iv1+1) ;
	  cout.width(10); cout << C1 ;
	  cout.width(10); cout << *(d1+C1) ;
	  cout.width(10); cout << *(d2+C1) << '\n';}
	      
      }
    }
  }

  Nth=0;
  //cout<< "\n\n      Order  No. till the order \n";
  for(i=0;i<=Nord;i++){
    for(j=0;j<=i;j++) Nth+=*(NM1+j)**(NM1+i-j);
    //cout.width(10); cout << i ;
    //cout.width(10); cout << Nth << '\n' ;
  }
  //cout <<"\n Dinv\n";
  //for(i=0;i<Lvec2;i++) {
     //cout.width(10); cout << i ;
     //cout.width(10); cout << *(Id1+i) << '\n' ;
  //}
   int icx=1;
   for(i=0;i<Nvr2;i++) {
      if(i>0) icx=icx*(Nord+1); else icx=1;
      *(ixref+i)=*(d1+icx);
      *(ixref+i+Nvr2)=*d1+(*(d2+icx));
   }
   //cout<< "\n  1st derivative reference \n";
   //for(i=0;i<Nvar;i++) {
   // cout.width(10); cout << i ;
   // cout.width(10); cout << *(ixref+i) << '\n' ;
   //}
   //delete [] iv1;
}

//
//
//  *******   DA Library   ***************************
//
//
da& da::operator=(const char* s)
{
   int Iadr[N_poly_max];
   double Coef[N_poly_max];

   int N_mon=poly_ana(s,Iadr,Coef);
   memset(v,0,sizeof(double)*Lvec);
   for(int i=0;i<N_mon;i++) *(v+Iadr[i])=Coef[i];
   return *this;
}

double da::coef(const char* s_in)
{
   int iadr;
   double coef;
   int n=strlen(s_in);
   char s[L_poly_max];
   for(int i=0;i<n;i++) {
     if(strchr("+-(){}",*(s_in+i))!=NULL) {
       cout<< " Input error ( remove +-(){} from " << s_in << ") \n";
       exit(1);
     }
   }
   strcpy(s,s_in);
   mon_ana('+',s,&iadr,&coef);
   return *(v+iadr);
}

double da::lin_da(int i)
{
   return(v[*(ixref+i)]);
}

void da::dBase(int i, double f=1.)
{
   for(int j=0;j<L_vec;j++) *(v+j)=0.;
   v[*(ixref+i)]=f;
}
void da::dBase(double* f)
{
   for(int j=0;j<L_vec;j++) *(v+j)=0.;
   for(int i=0;i<Nvar;i++) v[*(ixref+i)]=*(f+i);
}
void da::Dcin(int i,double f)
{
   v[*(ixref+i)]=f;
}


void da_mul(double* x,double* y,double* z)
{

  double xi;
  int C1M,C2M,C1,C2,iD2;
  int i,j,i1A,i2A,N2A,N1A,NA;
  int IMN,i1B0,i2B0,i1B,i2B,N2B;
  int No,N1A0;
  int is=0;

  memset(z0,0,sizeof(double)*Lvec);
  No=Nord;
  if(x[0]==0.) { is=1;}
  if(y[0]==0.) { No=Nord-1; }

/*   NA : order of x=x1*x2  ,  N1A : order of x1 , N2A : order of x2 
     i1A : */

  int i1A0=0;
  int i2A0=0;
  for(N2A=0;N2A<=No;N2A++) {
    if(is) { N1A0=1; } else {N1A0=0;}
    for(i2A=0;i2A<*(NM1+N2A);i2A++) {   
      /* loop up to number of N2A-th order*/
      N1A=N1A0;
      C2M=*(Id1+i2A+i2A0);
      for(i1A=is;i1A<*(NMS+No-N2A);i1A++){
	i=i1A+i1A0;
	xi=*(x+i);
	C1M=*(Id1+i1A);

	if(i1A==*(NMS+N1A)) N1A++;
	NA=N2A+N1A;

	i1B0=0;  i2B0=0;
	for(N2B=0;N2B<=Nord-NA;N2B++) {
	  for(i2B=0;i2B<*(NM1+N2B);i2B++) {
	    C2=C2M+*(Id1+i2B+i2B0);
	    iD2=*(d2+C2);
	    for(i1B=0;i1B<*(NMS+Nord-NA-N2B);i1B++) {
	      j=i1B+i1B0;
	      C1=C1M+*(Id1+i1B);
	      IMN=*(d1+C1)+iD2;
	      *(z0+IMN)+=xi*(*(y+j));
	    }
	    i1B0+=*(NMS+Nord-N2B);
	  }
	  i2B0+=*(NM1+N2B);
	}
      }
      i1A0+=*(NMS+Nord-N2A);
    }
    i2A0+=*(NM1+N2A);
    is=0;
  }

  memcpy(z,z0,sizeof(double)*Lvec);
}

da mul(const da& x,const da& y,int No,int Mx=0,int My=0)
{
   da z0;
  double xi;
  int C1M,C2M,C1,C2,iD2;
  int i,j,i1A,i2A,N2A,N1A,NA;
  int IMN,i1B0,i2B0,i1B,i2B,N2B;
  int Nox,is,N1A0;

  memset(z0.v,0,sizeof(double)*Lvec);
  Nox=No-My;
  //  Noy=No-Mx;
/*   NA : order of x=x1*x2  ,  N1A : order of x1 , N2A : order of x2 
     i1A : */

  
  int i1A0=0;
  int i2A0=0;
  for(N2A=0;N2A<=Nox;N2A++) {
    if(N2A>=Mx) { N1A0=0; is=0; }
    else {N1A0=Mx-N2A; is=*(NMS+N1A0-1);}

    for(i2A=0;i2A<*(NM1+N2A);i2A++) {   
      /* loop up to number of N2A-th order*/
       N1A=N1A0;
       C2M=*(Id1+i2A+i2A0);
       for(i1A=is;i1A<*(NMS+Nox-N2A);i1A++){
          i=i1A+i1A0;
          xi=x.v[i];
          C1M=*(Id1+i1A);

	if(i1A==*(NMS+N1A)) N1A++;
	NA=N2A+N1A;

	i1B0=0;  i2B0=0;
	for(N2B=0;N2B<=No-NA;N2B++) {
	  for(i2B=0;i2B<*(NM1+N2B);i2B++) {
	    C2=C2M+*(Id1+i2B+i2B0);
	    iD2=*(d2+C2);
	    for(i1B=0;i1B<*(NMS+No-NA-N2B);i1B++) {
	      j=i1B+i1B0;
	      C1=C1M+*(Id1+i1B);
	      IMN=*(d1+C1)+iD2;
	      z0.v[IMN]+=xi*y.v[j];
	    }
	    i1B0+=*(NMS+Nord-N2B);
	  }
	  i2B0+=*(NM1+N2B);
	}
      }
      i1A0+=*(NMS+Nord-N2A);
    }
    i2A0+=*(NM1+N2A);
  }
  return z0;
}

// y*x where x is No-th order polynomial
da muln(const da& y,const da& x,int No)
{
   da z0;
  double xi;
  int C1M,C2M,C1,C2,iD2;
  int i,j,i1A,i2A,N2A,N1A,NA;
  int IMN,i1B0,i2B0,i1B,i2B,N2B;
  int is,N1A0;

  memset(z0.v,0,sizeof(double)*Lvec);

/*   NA : order of x=x1*x2  ,  N1A : order of x1 , N2A : order of x2 
     i1A : */

  
  int i1A0=0;
  int i2A0=0;
  for(N2A=0;N2A<=No;N2A++) {
    N1A0=No-N2A; 
    if(N1A0) is=*(NMS+N1A0-1); else is=0;

    for(i2A=0;i2A<*(NM1+N2A);i2A++) {   
      /* loop up to number of N2A-th order*/
       N1A=N1A0;
       C2M=*(Id1+i2A+i2A0);
       for(i1A=is;i1A<*(NMS+No-N2A);i1A++){
          i=i1A+i1A0;
          xi=x.v[i];
          C1M=*(Id1+i1A);

	 if(i1A<*(NMS+N1A)) {	//N1A++;
	 NA=N2A+N1A;

	 i1B0=0;  i2B0=0;
	 for(N2B=0;N2B<=Nord-NA;N2B++) {
	  for(i2B=0;i2B<*(NM1+N2B);i2B++) {
	    C2=C2M+*(Id1+i2B+i2B0);
	    iD2=*(d2+C2);
	    for(i1B=0;i1B<*(NMS+Nord-NA-N2B);i1B++) {
	      j=i1B+i1B0;
	      C1=C1M+*(Id1+i1B);
	      IMN=*(d1+C1)+iD2;
	      z0.v[IMN]+=xi*y.v[j];
	    }
	    i1B0+=*(NMS+Nord-N2B);
	  }
	  i2B0+=*(NM1+N2B);
	 }
	 }
       }
       i1A0+=*(NMS+Nord-N2A);
    }
    i2A0+=*(NM1+N2A);
  }
  return z0;
}


da operator/(double c,const da& x)
{
  int i;
  da z1,z3,z;
  if(x.v[0]==0.) {
    cout<< " DA inverse is not usable because of x(1)=0 \n";
    exit(1);
  }
  double y1=c/(x.v[0]);

  z3=x/(-x.v[0]);
  z3.v[0]=0.;

  //     z3  : Expanding variable

  z1=z3;
  *z1.v=1.;
  for(i=2;i<=Nord;i++) {
     z1=mul(z3,z1,i,1,0);
     z1.v[0]=1.;
  }
  for(i=0;i<Lvec;i++) z.v[i]=y1*z1.v[i];
  return z;
}


da dpow(const da& x,int n,int M =0)
{
   da z;
   int No_cal,Mx;
   Mx=M;
   if(Mx==0 && x.v[0]==0.) Mx=1;
   if(n!=0) {
      z=x;
      if(Mx>0) {
	 if(Nord<n*Mx) {memset(z.v,0,sizeof(double)*Lvec); 
	    return z; }
         for(int i=1;i<n;i++) {
	    No_cal=Nord-(n-i-1)*Mx;
	    z=mul(z,x,No_cal,i,Mx);
	 }
      }
      else{
	 for(int i=1;i<n;i++) z.mul(z,x);
      }
   }
   else z=1.;
   return z;
}

da pow(const da& x,double k)
{

  da z1,z3,z4;

  if(x.v[0]==0.) {
    cout<< " DA pow is not usable because of x(0)=0 \n";
    exit(1);
  }

   // Bug ??? can not set break point hereafter but result is OK

  z3=x/(x.v[0]);
  z3.v[0]=0.;

  //     z3  : Expanding variable

  double f=(k-(double)(Nord-1))/(double)(Nord);

  if(N_ord>1) {
     z4.mul(z3,f);
     *z4.v=1.;
     for(int i=1;i<(Nord-1);i++) {
       f=(k-(double)(Nord-i-1))/(double)(Nord-i);
       z4=mul(z3,z4,i+1,1,0);
       z4.mul(z4,f);
       z4.v[0]=1.;
     }
  }
  else{ z4=1.; }
  z4.mul(z3,z4);
  z4.mul(z4,k);
  z4.v[0]=1.;

  double y1=pow(x.v[0],k);
  for(int i=0;i<Lvec;i++) z4.v[i]=z4.v[i]*y1;
  return z4;
}   


da sin(const da& x)
{
   da z1,z2,zc,zs,z;
   double f;
   int Mx=0;
   int No,Ne,i;
   
  double y1s=sin(x.v[0]);
  double y1c=cos(x.v[0]);
  if(x.v[0]==0.) Mx=1;
  z1=x;
  z1.v[0]=0.;

  // z1 Expanding variable  z2=z1^2
  // zc cos like da
  // zs sin like da
  z2.mul(z1,z1);

  if(Nord==0) {zc=1.; zs=0.; }
  else if(Nord==1) {zc=1.; zs=z1; }
  else if(Nord==2) {zc=1.-z2*0.5; zs=z1;}
  else{
    if(Nord%2) {No=Nord; Ne=Nord-1;} else {No=Nord-1; Ne=Nord;}
  //
  // sin like da
    f=-1./(double)(No*(No-1));
    zs.mul(z2,f);
    *zs.v=1.;
    for(i=2;i<No-2;i+=2) {
       f=-1./(double)((No-i)*(No-i-1));
       zs=mul(z2,zs,i+2,2,0); 
       zs.mul(zs,f);
       *zs.v=1.;
    }
    zs=mul(z1,zs,Nord,1,0);
  //
  //cos like da
    if(!Mx) {
      f=-1./(double)(Ne*(Ne-1));
      zc.mul(z2,f);
      *zc.v=1.;
      for(i=2;i<=Ne-2;i+=2) {
        f=-1./(double)((Ne-i)*(Ne-i-1));
        zc=mul(z2,zc,i+2,2,0); 
        zc.mul(zc,f);
        *zc.v=1.;
      }
    }
  }

  for(i=0;i<Lvec;i++) z.v[i]=zs.v[i]*y1c;
  if(!Mx) for(i=0;i<Lvec;i++) z.v[i]+=(zc.v[i]*y1s);
  return z;
}



da cos(const da& x)
{
   da z1,z2,zc,zs,z;
   double f;
   int Mx=0;
   int No,Ne,i;
   
  double y1s=sin(x.v[0]);
  double y1c=cos(x.v[0]);
  if(x.v[0]==0.) Mx=1;
  z1=x;
  z1.v[0]=0.;

  // z1 Expanding variable  z2=z1^2
  // zc cos like da
  // zs sin like da
  z2.mul(z1,z1);

  if(Nord==0) {zc=1.; zs=0.; }
  else if(Nord==1) {zc=1.; zs=z1;  }
  else if(Nord==2) {zc=1.-z2*0.5; zs=z1;}
  else{
    if(Nord%2) {No=Nord; Ne=Nord-1;} else {No=Nord-1; Ne=Nord;}
  //
  // sin like da
    if(!Mx){
    f=-1./(double)(No*(No-1));
    zs.mul(z2,f);
    *zs.v=1.;
    for(i=2;i<No-2;i+=2) {
       f=-1./(double)((No-i)*(No-i-1));
       zs=mul(z2,zs,i+2,2,0); 
       zs.mul(zs,f);
       *zs.v=1.;
    }
    zs=mul(z1,zs,Nord,1,0);
    }
  //
  //cos like da
//    if(!Mx) {
      f=-1./(double)(Ne*(Ne-1));
      zc.mul(z2,f);
      *zc.v=1.;
      for(i=2;i<=Ne-2;i+=2) {
        f=-1./(double)((Ne-i)*(Ne-i-1));
        zc=mul(z2,zc,i+2,2,0); 
        zc.mul(zc,f);
        *zc.v=1.;
      }
  //}
  }

  
  for(i=0;i<Lvec;i++) z.v[i]=zc.v[i]*y1c;
  if(!Mx) for(i=0;i<Lvec;i++) z.v[i]-=(zs.v[i]*y1s);
  return z;
}


da sinh(const da& x)
{
  da z1,z2,zc,zs,z;
  double f;
  int Mx=0;
  int No,Ne,i;

  double y1s=sinh(x.v[0]);
  double y1c=cosh(x.v[0]);
  if(x.v[0]==0.) Mx=1;
  z1=x;
  z1.v[0]=0.;

  // z1 Expanding variable  z2=z1^2
  // zc cos like da
  // zs sin like da
  z2.mul(z1,z1);

  if(Nord==0) {zc=1.; zs=0.; }
  else if(Nord==1) {zc=1.; zs=z1;  }
  else if(Nord==2) {zc=1.+z2*0.5; zs=z1;}
  else{
    if(Nord%2) {No=Nord; Ne=Nord-1;} else {No=Nord-1; Ne=Nord;}
  //
  // sin like da
    f=1./(double)(No*(No-1));
    zs.mul(z2,f);
    *zs.v=1.;
    for(i=2;i<No-2;i+=2) {
       f=1./(double)((No-i)*(No-i-1));
       zs=mul(z2,zs,i+2,2,0); 
       zs.mul(zs,f);
       *zs.v=1.;
    }
    zs=mul(z1,zs,Nord,1,0);
  //
  //cos like da
    if(!Mx) {
      f=1./(double)(Ne*(Ne-1));
      zc.mul(z2,f);
      *zc.v=1.;
      for(i=2;i<=Ne-2;i+=2) {
        f=1./(double)((Ne-i)*(Ne-i-1));
        zc=mul(z2,zc,i+2,2,0); 
        zc.mul(zc,f);
        *zc.v=1.;
      }
    }
  }

  for(i=0;i<Lvec;i++) z.v[i]=zc.v[i]*y1s;
  if(!Mx) for(i=0;i<Lvec;i++) z.v[i]+=(zs.v[i]*y1c);
  return z;
}

da cosh(const da& x)
{
   da z1,z2,zc,zs,z;
   double f;
   int Mx=0;
   int No,Ne,i;
   
  double y1s=sinh(x.v[0]);
  double y1c=cosh(x.v[0]);
  if(x.v[0]==0.) Mx=1;
  z1=x;
  z1.v[0]=0.;

  // z1 Expanding variable  z2=z1^2
  // zc cos like da
  // zs sin like da
  z2.mul(z1,z1);

  if(Nord==0) {zc=1.; zs=0.; }
  else if(Nord==1) {zc=1.; zs=z1;  }
  else if(Nord==2) {zc=1.+z2*0.5; zs=z1;}
  else{
    if(Nord%2) {No=Nord; Ne=Nord-1;} else {No=Nord-1; Ne=Nord;}
  //
  // sin like da
    if(!Mx){
    f=1./(double)(No*(No-1));
    zs.mul(z2,f);
    *zs.v=1.;
    for(i=2;i<No-2;i+=2) {
       f=1./(double)((No-i)*(No-i-1));
       zs=mul(z2,zs,i+2,2,0); 
       zs.mul(zs,f);
       *zs.v=1.;
    }
    zs=mul(z1,zs,Nord,1,0);
    }
  //
  //cos like da
      f=1./(double)(Ne*(Ne-1));
      zc.mul(z2,f);
      *zc.v=1.;
      for(i=2;i<=Ne-2;i+=2) {
        f=1./(double)((Ne-i)*(Ne-i-1));
        zc=mul(z2,zc,i+2,2,0); 
        zc.mul(zc,f);
        *zc.v=1.;
      }
  //}
  }

  for(i=0;i<Lvec;i++) z.v[i]=zc.v[i]*y1c;
  if(!Mx) for(i=0;i<Lvec;i++) z.v[i]+=(zs.v[i]*y1s);
  return z;
}



da exp(const da& x)
{
  int i;
  da z3,z4,z; 

  double y1=exp(x.v[0]);
  z3=x;
  z3.v[0]=0.;
  
  //     z3  : Expanding variable

  z4=z3/((double)(Nord));
  *z4.v=1.;
  for(i=1;i<Nord;i++) {
     z4=mul(z3,z4,i+1,1,0);
     z4=z4/((double)(Nord-i));
     *z4.v=1.;
  }
  for(i=0;i<Lvec;i++) {
     z.v[i]=z4.v[i]*y1;
  }
  return z;
}




da log(const da& x)
{
  int i;
  da z3,z4,z;
  if(x.v[0]<=0.) {
    cout <<" DA log is not usable because of x(1)<=0 \n";
    exit(1);
  }
  double y1=log(x.v[0]);

  z3=x/(x.v[0]);
  z3.v[0]=0.;

  //     z3  : Expanding variable
  if(Nord>1){
    double f=-(double)(Nord-1)/(double)Nord;
    z4.mul(z3,f);
    z4.v[0]=1.;
    for(i=1;i<Nord-1;i++){
      f=-(double)(Nord-i-1)/(double)(Nord-i);
      z4=mul(z3,z4,i+1,1,0);
      z4.mul(z4,f);
      z4.v[0]=1.;
    }
  }
  else {
    z4=1.;
  }
  z=z4*z3;
  z.v[0]=y1;
  return z;
}

da sqrt(const da& x)
{
  int i;
  da z1,z4,z3,z;
  if(x.v[0]<=0.) {
    cout <<" DA sqrt is not usable because of x(1)<=0 \n";
    exit(1);
  }

  double y1=sqrt(x.v[0]);

  z3=x/(x.v[0]);
  z3.v[0]=0.;

  //     z3  : Expanding variable

  double f=-(double)(2*Nord-3)/(double)(2*Nord);

  if(Nord>1) {
     z4.mul(z3,f);
     *z4.v=1.;
     for(i=1;i<(Nord-1);i++) {
       f=-(double)(2*(Nord-i)-3)/(double)(2*(Nord-i));
       z4=mul(z3,z4,i+1,1,0);
       z4.mul(z4,f);
       z4.v[0]=1.;
     }
  }
  else{ z4=1.; }
  z4.mul(z3,z4);
  z4.mul(z4,0.5);
  z4.v[0]=1.;
  for(i=0;i<Lvec;i++) z.v[i]=z4.v[i]*y1;
  return z;
}

da asin(const da& x)
{
  int i;
   double c[N_dim_max];
   
   // Make taylor coefficient table
   if(Nord>=N_dim_max) {
    cout<< " Nord is larger than max value (asin) "
	<< Nord << '>' << N_dim_max;
    exit(1);
  }
  double cr;
  for(int n=1;n<=N_dim_max;n++) {
   double cn=0.;
   double sign;
   for(int r=0;r<=n-1;r++) {
     if(r%2) sign=-1.; else sign=1.;
     int m=2*r-1;
     int k1=m;
     if(m<1) k1=1; else{
     while(m>=3){ 
        m=m-2;
        k1=k1*m;
     }}
     m=2*n-2*r-3;
     int k2=m;
     if(m<1) k2=1; else{
     while(m>=3){ 
        m=m-2;
        k2=k2*m;
     }}
     double x1=pow((1.+x.v[0]),-0.5-r)*pow((1.-x.v[0]),0.5-n+r);
     if(n!=1) 
         cr=sign*icomb(n-1,r)*k1*k2*x1;
     else cr=sign*k1*k2*x1;
     cn+=cr;
   }
   double c1=1.; for(i=0;i<n-1;i++) { c1=c1*0.5;}
   double nki;
   nki=(double)n;
   c[0]=asin(x.v[0]);
   for(i=1;i<n;i++) {
      nki=nki*(double)i; }
   c[n]=cn/nki*c1;
   if(fabs(c[n])<1.e-15) c[n]=0.;
//   cout << " c[" << n << "]= " << c[n] << '\n';
//  cout.flush();
  }
  // main of asin da
  da z4;
  da z1=x;
  
  z1.v[0]=0.;
  //   Closed orbit 1.e-10 m
  if(fabs(x.v[0])>1.e-10 && c[Nord]!=0. && c[Nord-1]!=0.){
    if(Nord>1) {
     double f=c[Nord]/c[Nord-1];
     z4.mul(z1,f);
     *z4.v=1.;       // z4=1+c_n/c_{n-1} z1
     for(i=1;i<(Nord-1);i++) {
       f=c[Nord-i]/c[Nord-i-1];
       z4=mul(z1,z4,i+1,1,0);
       z4.mul(z4,f);
       z4.v[0]=1.;
     }
    }
    else{ z4=1.; }
    z4.mul(z1,z4);
    z4.mul(z4,c[1]);
    z4.v[0]=c[0];        // c_0+c_1*z*z4
  }
  else {
     int No;
    if(Nord%2) {No=Nord;} else {No=Nord-1;}
    if(Nord!=1) {
    da z2;
    z2.mul(z1,z1);
    double f=c[No]/c[No-2];
    z4.mul(z2,f);
    *z4.v=1.;
    for(i=2;i<No-2;i+=2) {
       f=c[No-i]/c[No-i-2];
       z4=mul(z2,z4,i+1,1,0);
       z4.mul(z4,f);
       z4.v[0]=1.;
    }
    z4.mul(z1,z4);
    } else {z4=z1;}
    z4.mul(z4,c[1]);

  }
//   cout << "z4" << z4;
  return z4;
}

da atan(const da& x)
{
  da z1,z2,z4;
  double f;
  int No,i;

  z1=x;
  z1.v[0]=0.;
  if(fabs(x.v[0])>1.e-10) {
    z1=z1/(1.+x.v[0]*x.v[0]+x.v[0]*z1);
  }

  if(Nord%2) {No=Nord;} else {No=Nord-1;}
  if(Nord!=1) {
    da z2;
    z2.mul(z1,z1);
    f=-double(No-2)/double(No);
    z4.mul(z2,f);
    *z4.v=1.;
    for(i=2;i<No-2;i+=2) {
      f=-(double)(No-i-2)/(double)(No-i);
      z4=mul(z2,z4,i+2,2,0);
      z4.mul(z4,f);
      z4.v[0]=1.;
    }
    z4.mul(z1,z4);
  } 
  else {z4=z1;}
//   cout << "z4" << z4;
  if(fabs(x.v[0])>1.e-10) {
    z4+=atan(x.v[0]);
  }
  return z4;
}

da erf(const da& x)
{
   da xf,z;
   xf=0.;
   if(vget(x,0)!=0.) {
      cout << " x.v[0] must be zero in this version \n";
   }

   z=x*x;

   int No2=(Nord-1)/2;
  //   z  : Expanding variable

  xf=-(double)(2*No2-1)*z/((double)((2*No2+1)*No2));
  xf.ci(0,1.);
  for(int i=1;i<No2;i++) {
     xf=mul(z,xf,i,1,0);
     xf=-(double)(2*(No2-i)-1)*z/((double)((2*(No2-i)+1)*(No2-i)));
     xf.ci(0,1.);
  }
  return xf*x;
}

void da::daprint(const char* s)
{
  int C1M,C2M;
  int i,j,i1A,i2A,N2A,N1A,NA;
  double xi;
  
  cout << endl << endl << s << endl << "  daprint"
       << " Nord=" << Nord
       << " Nvar=" << Nvar
       << " Lvec=" << Lvec
       << endl << endl;
  int i1A0=0;
  int i2A0=0;
  for(N2A=0;N2A<=Nord;N2A++) {
    for(i2A=0;i2A<*(NM1+N2A);i2A++) {   
      /* loop up to number of N2A-th order*/
      N1A=0;
      C2M=*(Id1+i2A+i2A0);
      ivcal(C2M,Nord,Nvr2,iv2);
      for(i1A=0;i1A<*(NMS+Nord-N2A);i1A++){
	i=i1A+i1A0;
	xi=*(v+i);
	if(i1A==*(NMS+N1A)) N1A++;
	if(xi!=0.){
	  C1M=*(Id1+i1A);
	  ivcal(C1M,Nord,Nvr2,iv1);


	  NA=N2A+N1A;

	  {
	    int previous_width = cout.width(5);
	    cout << i; cout.width(previous_width);
	    cout << " : " << NA << " " << xi << " ";
	  }
	  for(j=0;j<Nvr2;j++) cout << " " << *(iv1+j) << " " << *(iv2+j);
	  cout << endl;
	}
      }
      i1A0+=*(NMS+Nord-N2A);
    }
    i2A0+=*(NM1+N2A);
  } 
}

da diag_part(const da& x)
{
  int C1M,C2M;
  int i,j,i1A,i2A,N2A,N1A,prsw;
  double xi;
  da y;
  
  y=0.;
  
  int i1A0=0;
  int i2A0=0;
  for(N2A=0;N2A<=Nord;N2A++) {
    for(i2A=0;i2A<*(NM1+N2A);i2A++) {   
      /* loop up to number of N2A-th order*/
      N1A=0;
      C2M=*(Id1+i2A+i2A0);
      ivcal(C2M,Nord,Nvr2,iv2);
      for(i1A=0;i1A<*(NMS+Nord-N2A);i1A++){
	i=i1A+i1A0;
	xi=x.v[i];
	if(i1A==*(NMS+N1A)) N1A++;
	if(xi!=0.){
	  C1M=*(Id1+i1A);
	  ivcal(C1M,Nord,Nvr2,iv1);

	  prsw=1;
	  for(j=0;j<Nvr2;j++) {
	     if(*(iv1+j)!=*(iv2+j)) prsw=0;
	  }
	  if(prsw) y.v[i]=xi;
	}
      }
      i1A0+=*(NMS+Nord-N2A);
    }
    i2A0+=*(NM1+N2A);
  } 
  return y;
}

void da::daprint(int max_ord)
{
  int C1M,C2M;
  int i,j,i1A,i2A,N2A,N1A,NA;
  double xi;
  
  cout << "{{" << Nord << ", " << Nvar << "}";
  int i1A0=0;
  int i2A0=0;
  for(N2A=0;N2A<=Nord;N2A++) {
    for(i2A=0;i2A<*(NM1+N2A);i2A++) {   
      /* loop up to number of N2A-th order*/
      N1A=0;
      C2M=*(Id1+i2A+i2A0);
      ivcal(C2M,Nord,Nvr2,iv2);
      for(i1A=0;i1A<*(NMS+Nord-N2A);i1A++){
	i=i1A+i1A0;
	xi=*(v+i);
	if(i1A==*(NMS+N1A)) N1A++;
	if(xi!=0.){
	  C1M=*(Id1+i1A);
	  ivcal(C1M,Nord,Nvr2,iv1);

	  NA=N2A+N1A;
	  if(NA<=max_ord) {
	    cout.setf(ios::hex);
	    cout << ",\n{" << NA;
	    for(j=0;j<Nvr2;j++) 
	      cout << *(iv1+j) << *(iv2+j);
	    cout.unsetf(ios::dec);
	    cout << ", ";
	    cout.width(15); cout << xi << '}';
	  }
	}
      }
      i1A0+=*(NMS+Nord-N2A);
    }
    i2A0+=*(NM1+N2A);
  } 
  cout << "};\n";
}

void da::NFPrint(ostream& s)
{
  int C1M,C2M;
  int i,j,i1A,i2A,N2A,N1A,NA;
  double xi;
  
  s.precision(15);
  s << "{{" << Nord << ", " << Nvar << "}";
  int i1A0=0;
  int i2A0=0;
  for(N2A=0;N2A<=Nord;N2A++) {
    for(i2A=0;i2A<*(NM1+N2A);i2A++) {   
      /* loop up to number of N2A-th order*/
      N1A=0;
      C2M=*(Id1+i2A+i2A0);
      ivcal(C2M,Nord,Nvr2,iv2);
      for(i1A=0;i1A<*(NMS+Nord-N2A);i1A++){
	i=i1A+i1A0;
	xi=v[i];
	if(i1A==*(NMS+N1A)) N1A++;
	if(xi!=0.){
	  C1M=*(Id1+i1A);
	  ivcal(C1M,Nord,Nvr2,iv1);

	  NA=N2A+N1A;
	  for(j=0;j<Nvr2;j++) {
	    if((*(iv1+j))%2 || (*(iv2+j))%2) goto NOWRITE;
	  }
	  s.setf(ios::hex);
	  s << ",\n{{" << NA;
	  for(j=0;j<Nvr2;j++) 
		s << ',' << *(iv1+j) << ',' << *(iv2+j);
	  s.unsetf(ios::dec);
	  s << "}, ";
          s << xi << '}';
	NOWRITE:;
	}
      }
      i1A0+=*(NMS+Nord-N2A);
    }
    i2A0+=*(NM1+N2A);
  }
  s << "};\n";
  s.precision(6);

}

ostream& operator<<(ostream& s,da& x)
{
  int C1M,C2M;
  int i,j,i1A,i2A,N2A,N1A,NA;
  double xi;
  
  s.precision(15);
  s << "{{" << Nord << ", " << Nvar << "}";
  int i1A0=0;
  int i2A0=0;
  for(N2A=0;N2A<=Nord;N2A++) {
    for(i2A=0;i2A<*(NM1+N2A);i2A++) {   
      /* loop up to number of N2A-th order*/
      N1A=0;
      C2M=*(Id1+i2A+i2A0);
      ivcal(C2M,Nord,Nvr2,iv2);
      for(i1A=0;i1A<*(NMS+Nord-N2A);i1A++){
	i=i1A+i1A0;
	xi=x.v[i];
	if(i1A==*(NMS+N1A)) N1A++;
	if(xi!=0.){
	  C1M=*(Id1+i1A);
	  ivcal(C1M,Nord,Nvr2,iv1);


	  NA=N2A+N1A;
	  s.setf(ios::hex);
	  s << ",\n{{" << NA;
	  for(j=0;j<Nvr2;j++) 
		s << ',' << *(iv1+j) << ',' << *(iv2+j);
	  s.unsetf(ios::dec);
	  s << "}, ";
          s << xi << '}';
	}
      }
      i1A0+=*(NMS+Nord-N2A);
    }
    i2A0+=*(NM1+N2A);
  }
  s << "};\n";
  s.precision(6);
  return s; 
}

double get_parm(char*,const char*);

istream& operator>>(istream& s,da& x)
{
  int C1M,C2M,Iaddr,readerr;
  char buf[256],*p,dlim;
  double xvec;

  int i,No,Nv,lbuf;

  readerr=0;
  s.get(buf,255);	if(!s) readerr=1;
  s.get(dlim);		if(!s) readerr=1;
  lbuf=strlen(buf);
  if(s.eof() || lbuf==0 || readerr) {
    cout << " DA read error\n"; exit(1);
  }
  p=strstr(buf,"{{")+2;
  sscanf(p,"%d %d",&No,&Nv);
  if(No<0 || No>30 || Nv<0 || Nv>30) {
    cout << " DA read error  Nord=" << No << "  Nvar=" << Nv << '\n';
    exit(1);
  }

  readerr=0;
  while(1) {
    s.get(buf,255);		if(!s) readerr=1;
    s.get(dlim);		if(!s) readerr=1;
    lbuf=strlen(buf);
    if(s.eof() || lbuf==0 || readerr) break;

    p=strchr(buf,'{')+1;
    while(*p==' ') p++;
    i=0;
    while(*p!=' ') {
      if(!isxdigit(*p)) cout << " DAread error  Not xdigit " << p;
      if(i%2) iv2[i/2]=*p-'0'; else iv1[i/2]=*p-'0';
      i++; p++;
    }
    if(i%2) {
      cout << " number of iv table is odd.\n";
      exit(1);
    }
    xvec=atof(p);
    C1M=ccal(Nord,Nvr2,iv1);
    C2M=ccal(Nord,Nvr2,iv2);
    Iaddr=*(d1+C1M)+*(d2+C2M);
    x.v[Iaddr]=xvec;
  }
  return s;
}

void pprint(da& x)
{
  int C1M,C2M;
  int i,j,i1A,i2A,N2A,N1A,ivj;
  double xi;
  
  cout << "\n polynomial print \n";
  int i1A0=0;
  int i2A0=0;
  for(N2A=0;N2A<=Nord;N2A++) {
    for(i2A=0;i2A<*(NM1+N2A);i2A++) {   
      /* loop up to number of N2A-th order*/
      N1A=0;
      C2M=*(Id1+i2A+i2A0);
      ivcal(C2M,Nord,Nvr2,iv2);
      for(i1A=0;i1A<*(NMS+Nord-N2A);i1A++){
	i=i1A+i1A0;
	xi=*(x.v+i);
	if(i1A==*(NMS+N1A)) N1A++;
	if(xi!=0.){
	  C1M=*(Id1+i1A);
	  ivcal(C1M,Nord,Nvr2,iv1);

	  cout.setf(ios::showpos);
	  cout<< xi; cout.unsetf(ios::showpos);
	  for(j=0;j<Nvar;j++) {
	     ivj=*(iv1+j);
	     if(ivj==1) cout<< "*" << var_name[j];
             if(ivj>1) cout<< "*" << var_name[j] << "^" << ivj;
	  }
	}
      }
      i1A0+=*(NMS+Nord-N2A);
    }
    i2A0+=*(NM1+N2A);
  } 
  cout << endl;
}

da dif(const da& x,const int k)
{
  int C1M,C2M;
  int i,i1A,i2A,N2A,N1A;
  int kk,nodk,Iaddr;
  da z;

  if(k>=Nvar) {
    cout<< " k(" << k << ") >= Nvar(" << Nvar << ") :  terminated \n";
    exit(1);
  }

  memset(z.v,0,sizeof(double)*Lvec);

  int i1A0=0;
  int i2A0=0;
  for(N2A=0;N2A<=Nord;N2A++) {
    for(i2A=0;i2A<*(NM1+N2A);i2A++) {   
      /* loop up to number of N2A-th order*/
      N1A=0;
      C2M=*(Id1+i2A+i2A0);
      ivcal(C2M,Nord,Nvr2,iv2);

/****** difinition of k    x:0 y:1 ... px:Nvr2 py:Nvr2+1.... */
      nodk=0;
      if(k>=Nvr2) {
	kk=k-Nvr2;
	nodk=*(iv2+kk);
	if(nodk>0) {
	  (*(iv2+kk))--;
	  C2M=ccal(Nord,Nvr2,iv2);
	}
      }

      for(i1A=0;i1A<*(NMS+Nord-N2A);i1A++){
	i=i1A+i1A0;
	if(i1A==*(NMS+N1A)) N1A++;

	C1M=*(Id1+i1A);
	ivcal(C1M,Nord,Nvr2,iv1);

        if(k<Nvr2) {
          kk=k;
          nodk=*(iv1+kk);
          if(nodk>0) {
            (*(iv1+kk))--;
            C1M=ccal(Nord,Nvr2,iv1);
          }
	}

        if(nodk>0) {
          Iaddr=*(d1+C1M)+*(d2+C2M);
          z.v[Iaddr]=x.v[i]*nodk;
	}
      }
      i1A0+=*(NMS+Nord-N2A);
    }
    i2A0+=*(NM1+N2A);
  }
  return z;
}
      

da itg(const da& x,const int k)
{
  da z0;
  int C1M,C2M;
  int i,i1A,i2A,N2A,N1A,NA;
  int kk,nodk,Iaddr;

  if(k>=Nvar) {
    cout<< " k(" << k << ") >= Nvar( " << Nvar << ") :  terminated \n";
    exit(1);
  }

  memset(z0.v,0,sizeof(double)*Lvec);

  int i1A0=0;
  int i2A0=0;
  for(N2A=0;N2A<=Nord;N2A++) {
    for(i2A=0;i2A<*(NM1+N2A);i2A++) {   
      /* loop up to number of N2A-th order*/
      N1A=0;
      C2M=*(Id1+i2A+i2A0);
      ivcal(C2M,Nord,Nvr2,iv2);

/****** difinition of k    x:0 y:1 ... px:Nvr2 py:Nvr2+1.... */
      nodk=0;
      if(k>=Nvr2) {
	kk=k-Nvr2;
	nodk=*(iv2+kk);
	if(nodk<Nord) {
	  (*(iv2+kk))++;
	  C2M=ccal(Nord,Nvr2,iv2);
	}
      }

      for(i1A=0;i1A<*(NMS+Nord-N2A);i1A++){
	i=i1A+i1A0;
	if(i1A==*(NMS+N1A)) N1A++;
        NA=N1A+N2A;

	C1M=*(Id1+i1A);
	ivcal(C1M,Nord,Nvr2,iv1);

        if(k<Nvr2) {
          kk=k;
          nodk=*(iv1+kk);
          if(nodk<Nord) {
            (*(iv1+kk))++;
            C1M=ccal(Nord,Nvr2,iv1);
          }
         }

        if(NA<Nord) {
          Iaddr=*(d1+C1M)+*(d2+C2M);
          *(z0.v+Iaddr)=*(x.v+i)/(nodk+1.);
	}
      }
      i1A0+=*(NMS+Nord-N2A);
    }
    i2A0+=*(NM1+N2A);
  }
  return z0;
}


da poi(const da& f,const da& g)
{
  da z1,z2,z3,z4;
  z4=0.;
  for(int i=0;i<Nvr2;i++){
    z1=dif(f,i);
    z2=dif(g,i+Nvr2);
    z3.mul(z1,z2);
    z4+=z3;

    z1=dif(f,i+Nvr2);
    z2=dif(g,i);
    z3.mul(z1,z2);
    z4-=z3;
  }
  return z4;
}

da line_itg(const da& x)
{
   da z;
   int C2M;
   int i,i1A,i2A,N2A,N1A,NA;
   z=0.;
   int i1A0=0;
   int i2A0=0;
   for(N2A=0;N2A<=Nord;N2A++) {
      for(i2A=0;i2A<*(NM1+N2A);i2A++) {   
      /* loop up to number of N2A-th order*/
         N1A=0;
         C2M=*(Id1+i2A+i2A0);
         ivcal(C2M,Nord,Nvr2,iv2);
         for(i1A=0;i1A<*(NMS+Nord-N2A);i1A++){
	    i=i1A+i1A0;
	    if(i1A==*(NMS+N1A)) N1A++;
	    if(x.v[i]!=0.){
	   //	  C1M=*(Id1+i1A);
	   //	  ivcal(C1M,Nord,Nvr2,iv1);
	       NA=N2A+N1A;
	       z.v[i]=x.v[i]/((double)(NA+1));
	    }
         }
         i1A0+=*(NMS+Nord-N2A);
      }
      i2A0+=*(NM1+N2A);
  } 
  return z;
}


// da::msk  get only k-th order of da x

void da::msk(const da& x,const int k)
{
  int C2M;
  int i,i1A,i2A,N2A,N1A,NA;
  double xi;

  memset(z0,0,sizeof(double)*Lvec);
  int i1A0=0;
  int i2A0=0;
  for(N2A=0;N2A<=Nord;N2A++) {
    for(i2A=0;i2A<*(NM1+N2A);i2A++) {   
      /* loop up to number of N2A-th order*/
      N1A=0;
      C2M=*(Id1+i2A+i2A0);
      ivcal(C2M,Nord,Nvr2,iv2);
      for(i1A=0;i1A<*(NMS+Nord-N2A);i1A++){
	i=i1A+i1A0;
	xi=*(x.v+i);
	if(i1A==*(NMS+N1A)) N1A++;
	if(xi!=0.){
	   //	  C1M=*(Id1+i1A);
	   //	  ivcal(C1M,Nord,Nvr2,iv1);
	  NA=N2A+N1A;
	  if(NA==k) *(z0+i)=xi;
	}
      }
      i1A0+=*(NMS+Nord-N2A);
    }
    i2A0+=*(NM1+N2A);
  } 
  memcpy(v,z0,sizeof(double)*Lvec);
}


/* *****************************************

      subroutine da_map(x,y,z)
      include 'dah'

      real*8 x(*),y(*),z

      z=0.d0
      nvl2=Nval/2
      do 10 i=1,NV
      ic1=id1inv+i-1)
      ic2=ilist(id2inv+i-1)

     ivcal(ic1,Nord,nvl2,iv1);
     ivcal(ic2,Nord,nvl2,iv2);

      yi=1.d0
      do 20 j=1,nvl2
        do 30 iv=1,ilist(iv1-1+j)
   30   yi=yi*y(j)
        do 31 iv=1,ilist(iv2-1+j)
   31   yi=yi*y(nvl2+j)
   20 continue

      z=z+x(i)*yi
   10 continue
      end

************************************************************ */

void da::varmsk(const da& x,const int k)
{
  int C2M;
  int i,i1A,i2A,N2A,N1A,NA;
  double xi;

  memset(z0,0,sizeof(double)*Lvec);
  int i1A0=0;
  int i2A0=0;
  for(N2A=0;N2A<=Nord;N2A++) {
    for(i2A=0;i2A<*(NM1+N2A);i2A++) {   
      /* loop up to number of N2A-th order*/
      N1A=0;
      C2M=*(Id1+i2A+i2A0);
      ivcal(C2M,Nord,Nvr2,iv2);
      for(i1A=0;i1A<*(NMS+Nord-N2A);i1A++){
	i=i1A+i1A0;
	xi=*(x.v+i);
	if(i1A==*(NMS+N1A)) N1A++;
	if(xi!=0.){
	   //	  C1M=*(Id1+i1A);
	   //	  ivcal(C1M,Nord,Nvr2,iv1);
	  NA=N2A+N1A;
	  if(NA==k) *(z0+i)=xi;
	}
      }
      i1A0+=*(NMS+Nord-N2A);
    }
    i2A0+=*(NM1+N2A);
  } 
  memcpy(v,z0,sizeof(double)*Lvec);
}


// Return min order of da.

int min_ord(const da& f)
{

  int C2M;
  int i,i1A,i2A,N2A,N1A,NA;
  int Nmin=Nord;
  int i1A0=0;
  int i2A0=0;
  for(N2A=0;N2A<=Nord;N2A++) {
    for(i2A=0;i2A<*(NM1+N2A);i2A++) {   
      /* loop up to number of N2A-th order*/
      N1A=0;
      C2M=*(Id1+i2A+i2A0);
      ivcal(C2M,Nord,Nvr2,iv2);
      for(i1A=0;i1A<*(NMS+Nord-N2A);i1A++){
	i=i1A+i1A0;
	if(i1A==*(NMS+N1A)) N1A++;
	if(f.v[i]!=0.){

	  NA=N2A+N1A;
	  if(NA<Nmin) Nmin=NA;
	}
      }
      i1A0+=*(NMS+Nord-N2A);
    }
    i2A0+=*(NM1+N2A);
  }
  return Nmin;
}
   
/*************************************************************

      subroutine da_tdo(F)
      include 'dah'
      real*8 f(*)

      open(11,FILE='td.out')
      write(11,'( a)') ' TITLE TOP '' Lie Hamiltonian |(RFn)| '''
      write(11,'( a)') ' SET ORDER X Y '
      write(11,'( a)') ' SET SCALE  Y LOG '

      nvl2=Nval/2
      do 9 i=1,NV
      if(f(i)==0.) goto 9
      ic1=ilist(id1inv+i-1)
      ic2=ilist(id2inv+i-1)

     ivcal(ic1,Nord,nvl2,iv1);
     ivcal(ic2,Nord,nvl2,iv2);
      write(11,100) float(i),dabs(F(i);,
     +    ilist(iv1),ilist(iv1+1),ilist(iv2),ilist(iv2+1)
 9    continue
 100  format(f5.1,E15.6,4I3)
      write(11,'( a)') ' PLOT'
      write(11,'( a)') ' JOIN 1'
      close(11)
      end


*************************************************************/
#include "concate.cc"
#include "oldmul.cc"
#include "perturb.cc"
