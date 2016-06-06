#include <iostream>
using std::cout;
#include <cstring>
using std::strstr;

#include <element.h>
#include <ph_rot.h>
// Phase space rotator (Test element)


Ph_rot::Ph_rot(char* s) : Element(s) 
{
  if(strstr(s,"{{")!=NULL && strstr(s,"}}")!=NULL) {
    ReadMathList(s,*Linmap);
  }
  else{
    nux=get_parm(s,"nux");
    nuy=get_parm(s,"nuy");
    nuz=get_parm(s,"nuz");
    mux=nux*pi2; muy=nuy*pi2; muz=nuz*pi2; 
    twiss[0]=get_parm(s,"ax"); twiss[1]=get_parm(s,"bx");
    twiss[2]=get_parm(s,"ay"); twiss[3]=get_parm(s,"by");
    twiss[4]=get_parm(s,"az"); twiss[5]=get_parm(s,"bz");
    twiss[6]=get_parm(s,"r1"); twiss[7]=get_parm(s,"r2");
    twiss[8]=get_parm(s,"r3"); twiss[9]=get_parm(s,"r4");
    twiss[10]=get_parm(s,"ex"); twiss[11]=get_parm(s,"epx");
    twiss[12]=get_parm(s,"ey"); twiss[13]=get_parm(s,"epy");
    twiss[14]=get_parm(s,"zx"); twiss[15]=get_parm(s,"zpx");
    twiss[16]=get_parm(s,"zy"); twiss[17]=get_parm(s,"zpy");
    if(twiss[1]==0.) twiss[1]=1.;
    if(twiss[3]==0.) twiss[3]=1.;
    if(twiss[5]==0.) twiss[5]=1.;

    Linmap=new matrix(6,6);

    matrix NtoP(6,6),PtoN(6,6),U(6,6);

    NtoP=set_NtoP(twiss);
    PtoN=SInverse(NtoP);
    U=set_U(nux,nuy,nuz);

    *Linmap=(matrix)NtoP*U*PtoN;

  }
}

Ph_rot::Ph_rot(double* K) : Element(K) 
{ 
  int i;
  for(i=0;i<18;i++) {twiss[i]=K[i+1];}
  nux=K[19]; nuy=K[20]; nuz=K[21];
  mux=nux*pi2; muy=nuy*pi2; muz=nuz*pi2;
  if(twiss[1]==0.) {twiss[1]=1.; cout << "bx set to 1";}
  if(twiss[3]==0.) {twiss[3]=1.; cout << "by set to 1";}
  if(twiss[5]==0.) {twiss[5]=1.; cout << "bz set to 1";}

  Linmap=new matrix(6,6);
  matrix NtoP(6,6),PtoN(6,6),U(6,6);

  NtoP=set_NtoP(twiss);
  PtoN=SInverse(NtoP);
  U=set_U(nux,nuy,nuz);
  
  *Linmap=(matrix)NtoP*U*PtoN;
}

void Ph_rot::Mapping(map_double& x)
{
  double x0,x1,x2,x3,x4,x5;

  x0=(*Linmap)[0][0]*x[0]+(*Linmap)[0][1]*x[3]
    +(*Linmap)[0][2]*x[1]+(*Linmap)[0][3]*x[4]
    +(*Linmap)[0][4]*x[2]+(*Linmap)[0][5]*x[5];
  x3=(*Linmap)[1][0]*x[0]+(*Linmap)[1][1]*x[3]
    +(*Linmap)[1][2]*x[1]+(*Linmap)[1][3]*x[4]
    +(*Linmap)[1][4]*x[2]+(*Linmap)[1][5]*x[5];
  x1=(*Linmap)[2][0]*x[0]+(*Linmap)[2][1]*x[3]
    +(*Linmap)[2][2]*x[1]+(*Linmap)[2][3]*x[4]
    +(*Linmap)[2][4]*x[2]+(*Linmap)[2][5]*x[5];
  x4=(*Linmap)[3][0]*x[0]+(*Linmap)[3][1]*x[3]
    +(*Linmap)[3][2]*x[1]+(*Linmap)[3][3]*x[4]
    +(*Linmap)[3][4]*x[2]+(*Linmap)[3][5]*x[5];
  x2=(*Linmap)[4][0]*x[0]+(*Linmap)[4][1]*x[3]
    +(*Linmap)[4][2]*x[1]+(*Linmap)[4][3]*x[4]
    +(*Linmap)[4][4]*x[2]+(*Linmap)[4][5]*x[5];
  x5=(*Linmap)[5][0]*x[0]+(*Linmap)[5][1]*x[3]
    +(*Linmap)[5][2]*x[1]+(*Linmap)[5][3]*x[4]
    +(*Linmap)[5][4]*x[2]+(*Linmap)[5][5]*x[5];

  x[0]=x0; x[3]=x3;
  x[1]=x1; x[4]=x4;
  x[2]=x2; x[5]=x5;
}

void Ph_rot::Mapping(map_da& x)
{
  da x0,x1,x2,x3,x4,x5;

  x0=(*Linmap)[0][0]*x[0]+(*Linmap)[0][1]*x[3]
    +(*Linmap)[0][2]*x[1]+(*Linmap)[0][3]*x[4]
    +(*Linmap)[0][4]*x[2]+(*Linmap)[0][5]*x[5];
  x3=(*Linmap)[1][0]*x[0]+(*Linmap)[1][1]*x[3]
    +(*Linmap)[1][2]*x[1]+(*Linmap)[1][3]*x[4]
    +(*Linmap)[1][4]*x[2]+(*Linmap)[1][5]*x[5];
  x1=(*Linmap)[2][0]*x[0]+(*Linmap)[2][1]*x[3]
    +(*Linmap)[2][2]*x[1]+(*Linmap)[2][3]*x[4]
    +(*Linmap)[2][4]*x[2]+(*Linmap)[2][5]*x[5];
  x4=(*Linmap)[3][0]*x[0]+(*Linmap)[3][1]*x[3]
    +(*Linmap)[3][2]*x[1]+(*Linmap)[3][3]*x[4]
    +(*Linmap)[3][4]*x[2]+(*Linmap)[3][5]*x[5];
  x2=(*Linmap)[4][0]*x[0]+(*Linmap)[4][1]*x[3]
    +(*Linmap)[4][2]*x[1]+(*Linmap)[4][3]*x[4]
    +(*Linmap)[4][4]*x[2]+(*Linmap)[4][5]*x[5];
  x5=(*Linmap)[5][0]*x[0]+(*Linmap)[5][1]*x[3]
    +(*Linmap)[5][2]*x[1]+(*Linmap)[5][3]*x[4]
    +(*Linmap)[5][4]*x[2]+(*Linmap)[5][5]*x[5];

  x[0]=x0; x[3]=x3;
  x[1]=x1; x[4]=x4;
  x[2]=x2; x[5]=x5;
}

void Ph_rot::Mapping(map_p_da& x)
{
  p_da x0,x1,x2,x3,x4,x5;

  x0=(*Linmap)[0][0]*x[0]+(*Linmap)[0][1]*x[3]
    +(*Linmap)[0][2]*x[1]+(*Linmap)[0][3]*x[4]
    +(*Linmap)[0][4]*x[2]+(*Linmap)[0][5]*x[5];
  x3=(*Linmap)[1][0]*x[0]+(*Linmap)[1][1]*x[3]
    +(*Linmap)[1][2]*x[1]+(*Linmap)[1][3]*x[4]
    +(*Linmap)[1][4]*x[2]+(*Linmap)[1][5]*x[5];
  x1=(*Linmap)[2][0]*x[0]+(*Linmap)[2][1]*x[3]
    +(*Linmap)[2][2]*x[1]+(*Linmap)[2][3]*x[4]
    +(*Linmap)[2][4]*x[2]+(*Linmap)[2][5]*x[5];
  x4=(*Linmap)[3][0]*x[0]+(*Linmap)[3][1]*x[3]
    +(*Linmap)[3][2]*x[1]+(*Linmap)[3][3]*x[4]
    +(*Linmap)[3][4]*x[2]+(*Linmap)[3][5]*x[5];
  x2=(*Linmap)[4][0]*x[0]+(*Linmap)[4][1]*x[3]
    +(*Linmap)[4][2]*x[1]+(*Linmap)[4][3]*x[4]
    +(*Linmap)[4][4]*x[2]+(*Linmap)[4][5]*x[5];
  x5=(*Linmap)[5][0]*x[0]+(*Linmap)[5][1]*x[3]
    +(*Linmap)[5][2]*x[1]+(*Linmap)[5][3]*x[4]
    +(*Linmap)[5][4]*x[2]+(*Linmap)[5][5]*x[5];
   
  x[0]=x0; x[3]=x3;
  x[1]=x1; x[4]=x4;
  x[2]=x2; x[5]=x5;
}

void Ph_rot::Mapping(pBeam& x)
{
  int i;
  double x0,x1,x2,x3,x4,x5;

  for(i=0;i<x.np;i++) {
    x0=(*Linmap)[0][0]*x.x[i]+(*Linmap)[0][1]*x.px[i]
      +(*Linmap)[0][2]*x.y[i]+(*Linmap)[0][3]*x.py[i]
      +(*Linmap)[0][4]*x.z[i]+(*Linmap)[0][5]*x.pz[i];
    x3=(*Linmap)[1][0]*x.x[i]+(*Linmap)[1][1]*x.px[i]
      +(*Linmap)[1][2]*x.y[i]+(*Linmap)[1][3]*x.py[i]
      +(*Linmap)[1][4]*x.z[i]+(*Linmap)[1][5]*x.pz[i];
    x1=(*Linmap)[2][0]*x.x[i]+(*Linmap)[2][1]*x.px[i]
      +(*Linmap)[2][2]*x.y[i]+(*Linmap)[2][3]*x.py[i]
      +(*Linmap)[2][4]*x.z[i]+(*Linmap)[2][5]*x.pz[i];
    x4=(*Linmap)[3][0]*x.x[i]+(*Linmap)[3][1]*x.px[i]
      +(*Linmap)[3][2]*x.y[i]+(*Linmap)[3][3]*x.py[i]
      +(*Linmap)[3][4]*x.z[i]+(*Linmap)[3][5]*x.pz[i];
    x2=(*Linmap)[4][0]*x.x[i]+(*Linmap)[4][1]*x.px[i]
      +(*Linmap)[4][2]*x.y[i]+(*Linmap)[4][3]*x.py[i]
      +(*Linmap)[4][4]*x.z[i]+(*Linmap)[4][5]*x.pz[i];
    x5=(*Linmap)[5][0]*x.x[i]+(*Linmap)[5][1]*x.px[i]
      +(*Linmap)[5][2]*x.y[i]+(*Linmap)[5][3]*x.py[i]
      +(*Linmap)[5][4]*x.z[i]+(*Linmap)[5][5]*x.pz[i];
    x.x[i]=x0; x.px[i]=x3;
    x.y[i]=x1; x.py[i]=x4;
    x.z[i]=x2; x.pz[i]=x5;
  }
}


