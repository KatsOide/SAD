#include <iostream>
using std::istream;
using std::ostream;
#include <cstdlib>
using std::exit;
#include <cstring>
using std::strcmp;
using std::strstr;
#include <cmath>
using std::sqrt;

#include <matrix.h>
#include <element.h>
#include <track.h>

pBeam::pBeam(char* str,int N,double rn)
{
  int i;
  N_particle=N;
  ndump=N;
//  np=N;
  x=new double[6*N_particle];
  y=x+N_particle;
  z=y+N_particle;
  px=z+N_particle;
  py=px+N_particle;
  pz=py+N_particle;
  sv=new int[N_particle];
  if(strstr(str,"Init")!=NULL) {
    for(i=0;i<N_particle*6;i++) x[i]=0.;
    for(i=0;i<N_particle;i++) sv[i]=i+1;
    np=N_particle;
  }
  rn_particle=rn;
}


pBeam::pBeam(int N,double rn)
{
  N_particle=N;
  ndump=N;
//  np=N;
  x=new double[6*N_particle];
  y=x+N_particle;
  z=y+N_particle;
  px=z+N_particle;
  py=px+N_particle;
  pz=py+N_particle;
  sv=new int[N_particle];
  rn_particle=rn;
}


pBeam::pBeam(const pBeam& pb) 
{
  int i;
  N_particle=pb.N_particle;
  ndump=pb.ndump;
  np=pb.np;
  x=new double[6*N_particle];
  y=x+N_particle;
  z=y+N_particle;
  px=z+N_particle;
  py=px+N_particle;
  pz=py+N_particle;
  sv=new int[N_particle];
  rn_particle=pb.rn_particle;
  for(i=0;i<N_particle*6;i++) x[i]=pb.x[i]; 
  for(i=0;i<N_particle;i++) sv[i]=pb.sv[i];
}
      
pBeam& pBeam::operator=(const pBeam& pb) 
{
  int i;
  if(N_particle!=pb.N_particle) {
    cout << "N_particle is inconsistent"
      << N_particle << " " << pb.N_particle << '\n';
    exit(1);
  }
  if(rn_particle!=pb.rn_particle) {
    cout << "rn_particle is inconsistent"
      << rn_particle << " " << pb.rn_particle << '\n';
  }
  np=pb.np;
  ndump=pb.ndump;
  for(i=0;i<N_particle*6;i++) x[i]=pb.x[i]; 
  for(i=0;i<N_particle;i++) sv[i]=pb.sv[i];
  return *this;
}
      
pBeam& pBeam::operator=(double f)
{
  int i;
  np=N_particle;
  for(i=0;i<N_particle*6;i++) x[i]=f;
  for(i=0;i<N_particle;i++) sv[i]=i+1;
  return *this;
}

map_double pBeam::particle(int i)
{
  map_double xx;
  for(int j=0;j<6;j++) xx[j]=x[i*6+j];
  xx.setsv(sv[i]);
  return xx;
}


ostream& operator<<(ostream& s,pBeam& x)
{
  s.precision(5);
  if(x.ndump>x.N_particle) x.ndump=x.N_particle;
  for(int i=0;i<x.ndump;i++) {
    s.width(11); s << x.x[i] << ' ';
    s.width(11); s << x.px[i] << ' '; 
    s.width(11); s << x.y[i] << ' ';
    s.width(11); s << x.py[i] << ' ';
    s.width(11); s << x.z[i] << ' ';
    s.width(11); s << x.pz[i] << ' ';
    s.width(4);  s << x.sv[i] << '\n';
  }
  s.precision(0);
  return s;
}


istream& operator>>(istream& s,pBeam& x)
{
  int sv,i;
  for(i=0;i<x.N_particle;i++) {
    s >> x.x[i] >> x.px[i] >> x.y[i] >> x.py[i] >> 
      x.z[i] >> x.pz[i] >> sv;
    if(sv>0) x.sv[i]=i+1; else x.sv[i]=-i-1;
    if(!s) break;
  }
  x.np=i;
  for(i=x.np;i<x.N_particle;i++) {
    x.x[i]=0.; x.px[i]=0.;
    x.y[i]=0.; x.py[i]=0.; 
    x.z[i]=0.; x.pz[i]=0.;  x.sv[i]=-i-1;
  }
  cout << x.np << " particles were read (Allocated=" << x.N_particle << ")\n";
  cout.flush();
  return s;
}

void MathListOut(pBeam& x)
{
//  cout.precision(10);

  cout << "{";
  for(int i=0;i<x.N_particle;i++) {
    cout << '{' << x.x[i] << ',' << x.px[i] << ',' 
      << x.y[i] << ',' << x.py[i] << ','
      << x.z[i] << ',' << x.pz[i] << ',' << x.sv[i];
    if(i!=x.N_particle-1) cout << "}\n";
    else cout << "}};\n";
  }
}


void is_survive(pBeam& z)
{
  pBeam x=z;
  int ii=0,ik=x.np;
  for(int i=0;i<x.np;i++) {
    if(x.x[i]*x.x[i]+x.y[i]*x.y[i]<0.01) {
      z.x[ii]=x.x[i];
      z.px[ii]=x.px[i];
      z.y[ii]=x.y[i];
      z.py[ii]=x.py[i];
      z.z[ii]=x.z[i];
      z.pz[ii]=x.pz[i];
      z.sv[ii]=x.sv[i];
      ii++;
    }
    else {
      if(x.sv[i]>0) {
	ik--;
	z.x[ik]=x.x[i];
	z.px[ik]=x.px[i];
	z.y[ik]=x.y[i];
	z.py[ik]=x.py[i];
	z.z[ik]=x.z[i];
	z.pz[ik]=x.pz[i];
	z.sv[ik]=-x.sv[i];
      }
    }
  }
  z.np=ii;
//  cout << " ii = " << ii << "  ik=" << ik; cout.flush();
}

float rgauss();
float ran1();
double tgauss();
double tran();
void Normalize(matrix&, double*);

void pBeam::Initialize(const char* s,const matrix& Benv)
{
  double eig[12],amp[6];
  matrix R(6,6),Ri(6,6);
  matrix tmp=Benv;
  int n,i,j;

  teigen(Benv,R,eig);
  Normalize(R,eig);
  cout << "\n Beam envelope matrix \n" << tmp;
  cout << "\n Eigen matrix diagonalizing Benv matrix\n" << R;
  cout << "\n Eigen values of the Benv matrix (sqrt)\n";
  for(i=0;i<6;i++) {cout.width(12);  cout << eig[2*i] << ' ';} cout << '\n';
//  Ri=Transpose(R);
//  cout << Ri*Benv*R;
  cout.flush();

  if(!strcmp(s,"Gaussian")) {
    for(n=0;n<N_particle;n++) {
      for(i=0;i<6;i++) amp[i]=sqrt(eig[2*i])*tgauss();
      x[n]=0.; px[n]=0.; y[n]=0.; py[n]=0.; z[n]=0.; pz[n]=0.;
      for(j=0;j<6;j++) {
	x[n]+=R[0][j]*amp[j];
	px[n]+=R[1][j]*amp[j];
	y[n]+=R[2][j]*amp[j];
	py[n]+=R[3][j]*amp[j];
	z[n]+=R[4][j]*amp[j];
	pz[n]+=R[5][j]*amp[j];
      }
      sv[n]=n+1;
    }
  }
  else if(!strcmp(s,"Uniform")) {
    for(n=0;n<N_particle;n++) {
      for(i=0;i<6;i++) amp[i]=sqrt(eig[2*i])*(2.*tran()-1.);
      x[n]=0.; px[n]=0.; y[n]=0.; py[n]=0.; z[n]=0.; pz[n]=0.;
      for(j=0;j<6;j++) {
	x[n]+=R[0][j]*amp[j];
	px[n]+=R[1][j]*amp[j];
	y[n]+=R[2][j]*amp[j];
	py[n]+=R[3][j]*amp[j];
	z[n]+=R[4][j]*amp[j];
	pz[n]+=R[5][j]*amp[j];
      }
      sv[n]=n+1;
    }
  }
  np=N_particle;
}

// !!!!!!!!! Forbidden !!!!!!!!!!!!!!
//void pBeam::Initialize(const char* s,const Beam& BEAM)

/*
void pBeam::ReadFile(const char* s)
{
  ifstream trin(s,ios::in);
  trin >> 
  trin.close();
}
*/
map_double CenterofMass(const pBeam& x)
{
  int n;
  map_double y;
  y=0.;
  for(n=0;n<x.np;n++) {
    y[0]+=x.x[n];
    y[1]+=x.y[n];
    y[2]+=x.z[n];
    y[3]+=x.px[n];
    y[4]+=x.py[n];
    y[5]+=x.pz[n];
  }
  for(n=0;n<6;n++) y[n]=y[n]/(double)x.np;
  return y;
}
  
void BeamSizeMonitor(const pBeam& x,double* beam)
{
  int j;
  for(j=0;j<27;j++) beam[j]=0.;
  for(j=0;j<x.np;j++) {
    beam[0]+=x.x[j];
    beam[1]+=x.y[j];
    beam[2]+=x.z[j];
    beam[3]+=x.px[j];
    beam[4]+=x.py[j];
    beam[5]+=x.pz[j];
    beam[6]+=x.x[j]*x.x[j];
    beam[7]+=x.x[j]*x.px[j];
    beam[8]+=x.x[j]*x.y[j];
    beam[9]+=x.x[j]*x.py[j];
    beam[10]+=x.x[j]*x.z[j];
    beam[11]+=x.x[j]*x.pz[j];
    beam[12]+=x.px[j]*x.px[j];
    beam[13]+=x.px[j]*x.y[j];
    beam[14]+=x.px[j]*x.py[j];
    beam[15]+=x.px[j]*x.z[j];
    beam[16]+=x.px[j]*x.pz[j];
    beam[17]+=x.y[j]*x.y[j];
    beam[18]+=x.y[j]*x.py[j];
    beam[19]+=x.y[j]*x.z[j];
    beam[20]+=x.y[j]*x.pz[j];
    beam[21]+=x.py[j]*x.py[j];
    beam[22]+=x.py[j]*x.z[j];
    beam[23]+=x.py[j]*x.pz[j];
    beam[24]+=x.z[j]*x.z[j];
    beam[25]+=x.z[j]*x.pz[j];
    beam[26]+=x.pz[j]*x.pz[j];
  }
  for(j=0;j<27;j++) {
    beam[j]=beam[j]/x.np;
  }
}

void BeamSizeMonitor(const pBeam& x,double* beam,matrix& Benv)
{
  if(getNr(Benv)!=6 || getNc(Benv)!=6) {
    cout << "Benv should be 6x6 matrix  @BeamSizeMonitor\n";
    exit(1);
  }
  int j,i;
  for(j=0;j<6;j++) beam[j]=0.;
  Benv=0;
  for(j=0;j<x.np;j++) {
    beam[0]+=x.x[j];
    beam[1]+=x.px[j];
    beam[2]+=x.y[j];
    beam[3]+=x.py[j];
    beam[4]+=x.z[j];
    beam[5]+=x.pz[j];
    Benv[0][0]+=x.x[j]*x.x[j];
    Benv[0][1]+=x.x[j]*x.px[j];
    Benv[0][2]+=x.x[j]*x.y[j];
    Benv[0][3]+=x.x[j]*x.py[j];
    Benv[0][4]+=x.x[j]*x.z[j];
    Benv[0][5]+=x.x[j]*x.pz[j];
    Benv[1][1]+=x.px[j]*x.px[j];
    Benv[1][2]+=x.px[j]*x.y[j];
    Benv[1][3]+=x.px[j]*x.py[j];
    Benv[1][4]+=x.px[j]*x.z[j];
    Benv[1][5]+=x.px[j]*x.pz[j];
    Benv[2][2]+=x.y[j]*x.y[j];
    Benv[2][3]+=x.y[j]*x.py[j];
    Benv[2][4]+=x.y[j]*x.z[j];
    Benv[2][5]+=x.y[j]*x.pz[j];
    Benv[3][3]+=x.py[j]*x.py[j];
    Benv[3][4]+=x.py[j]*x.z[j];
    Benv[3][5]+=x.py[j]*x.pz[j];
    Benv[4][4]+=x.z[j]*x.z[j];
    Benv[4][5]+=x.z[j]*x.pz[j];
    Benv[5][5]+=x.pz[j]*x.pz[j];
  }
  for(j=0;j<6;j++) {
    beam[j]=beam[j]/x.np;
    for(i=0;i<6;i++) {
      Benv[j][i]/=x.np;
    }
  } 
  for(j=0;j<6;j++) {
    for(i=j;i<6;i++) {
      Benv[j][i]-=beam[j]*beam[i];
    }
  }
  Benv[1][0]=Benv[0][1];
  Benv[2][0]=Benv[0][2];
  Benv[3][0]=Benv[0][3];
  Benv[4][0]=Benv[0][4];
  Benv[5][0]=Benv[0][5];
  Benv[2][1]=Benv[1][2];
  Benv[3][1]=Benv[1][3];
  Benv[4][1]=Benv[1][4];
  Benv[5][1]=Benv[1][5];
  Benv[3][2]=Benv[2][3];
  Benv[4][2]=Benv[2][4];
  Benv[5][2]=Benv[2][5];
  Benv[4][3]=Benv[3][4];
  Benv[5][3]=Benv[3][5];
  Benv[5][4]=Benv[4][5];
}


