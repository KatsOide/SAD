#include <iostream>
using std::cout;
using std::endl;
#include <cstdlib>
using std::exit;

// This file defines following 4 variable.
// Do not include <dacpp.h> to avoid  muti-define.
int N_ord,N_var,L_vec,N_cv;
int Np_ord,Np_var,Lp_vec;

void ivcalp(int c1,int n,int v,int* iv)
{
  int ci,c,i;

  ci=c1;
  for(i=0;i<v;i++) {
    c=ci/(n+1);
    *(iv+v-i-1)=ci-c*(n+1);
    ci=c;
  }
}


void ivcal(int c1,int n,int v,int* iv)
{
  int ci,c,i;

  ci=c1;
  for(i=0;i<v;i++){
    c=ci/(n+1);
    *(iv+i)=ci-c*(n+1);
    ci=c;
  }
}


int ccal(int n,int v2,int* iv)
{
  int c1=0,i,nn,j;

  for(i=0;i<v2;i++){
    nn=*(iv+i);
    if(i>0) { for(j=0;j<i;j++) nn=nn*(n+1);}
      c1+=nn;
  }
  return(c1);
}



int icomb(int n,int m)
{

  double a,b,c;
  int ic,i;

  if(n<m) {
     cout << " Not suported 1st arg<2nd arg in Comb " << n << " " << m << endl;
     exit(1);
  }
  a=n;
  b=m;
  if(m>0){
  for(i=1;i<m;i++) {
    a=a*(n-i);
    b=b*(m-i);
  }
  c=a/b; } else c=1.;
  ic=(int)(c+0.5);
  //  cout << " Combination check  " << c << endl;
  return(ic);
}
