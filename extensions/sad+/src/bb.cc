#include <iostream>
using std::cout;
#include <cstdlib>
using std::exit;
#include <cmath>
using std::fabs;
using std::sqrt;
using std::pow;
using std::sin;
using std::cos;
using std::tan;
using std::exp;

#include <Complex.h>
#include <phys_const.h>
#include <map_double.h>
#include <c_da.h>
#include <map_da.h>
#include <c_p_da.h>
#include <map_p_da.h>
#include <element.h>
#include <bb.h>

#define KO_SLICE 0
#define max_n_slice 50

Complex cerf(Complex);
double erf_fl(double);
matrix set_Cros(double);
void m_zslice(matrix&,matrix&,double*);
void Normalize(matrix&,double*);
double benv_norm_axis(matrix&,double*);
//extern void teigen(double*,double*,int ,int );
void KOslicing(int,double sigz,double* z_slice);
void KHslicing(int,double sigz,double* z_slice);
double SliceBenv(matrix*,matrix*,double*);
void BenvCalc(double emx,double emy,double emz,double* twiss,
	 matrix* Beam_envelope);

// Interaction Pointbeam


IP::IP(char *s) : Element(s)
{
  int i;
  double twiss[21];
  lmon=0;
  x_angle=get_parm(s,"x_angle"); 
  twiss[0]=get_parm(s,"ax");  twiss[1]=get_parm(s,"bx");
  twiss[10]=get_parm(s,"ex"); twiss[11]=get_parm(s,"epx");
  twiss[2]=get_parm(s,"ay");  twiss[3]=get_parm(s,"by");
  twiss[12]=get_parm(s,"ey"); twiss[13]=get_parm(s,"epy");
  twiss[4]=get_parm(s,"az");  twiss[5]=get_parm(s,"bz");
  twiss[14]=get_parm(s,"zx"); twiss[15]=get_parm(s,"zpx");
  twiss[16]=get_parm(s,"zy"); twiss[17]=get_parm(s,"zpy");
  twiss[6]=get_parm(s,"r1");  twiss[7]=get_parm(s,"r2");
  twiss[8]=get_parm(s,"r3");  twiss[9]=get_parm(s,"r4");
  twiss[18]=(1.+twiss[0]*twiss[0])/twiss[1];  //gx
  twiss[19]=(1.+twiss[2]*twiss[2])/twiss[3];  //gy
  twiss[20]=(1.+twiss[4]*twiss[4])/twiss[5];  //gz

  double emx,emy,emz;
  emx=get_parm(s,"emx");  emy=get_parm(s,"emy");  emz=get_parm(s,"emz");
  
  cosx=cos(x_angle);
  sinx=sin(x_angle);
  tanx=tan(x_angle);

  cod[0]=get_parm(s,"dx");
  cod[3]=get_parm(s,"dpx");
  cod[1]=get_parm(s,"dy");
  cod[4]=get_parm(s,"dpy");
  cod[2]=get_parm(s,"dz");
  cod[5]=get_parm(s,"dpz");
  if(x_angle!=0.) {
    cod[3]=cod[3]/cosx;
    cod[4]=cod[4]/cosx;
    cod[2]=cod[2]/cosx;
  }

  N_particle=get_parm(s,"N_particle");

   // Beam envelope

  Beam_envelope=new matrix(N_cv,N_cv);

  BenvCalc(emx,emy,emz,twiss,Beam_envelope);
// cout << "Beam envelope matrix (Crossing)\n" << *Beam_envelope;

   // Get 5x5 matrix
  Benv_slice=new matrix(5,5);
  tilt_angle= SliceBenv(Beam_envelope,Benv_slice,bcen_fac);
  sin_t=sin(tilt_angle);
  cos_t=cos(tilt_angle);

// -------------------------------------------------------------------
// Beam slicing
// -------------------------------------------------------------------

  n_slice=(int)get_parm(s,"N_slice");
  np_slice=new double[n_slice];
  z_slice=new double[n_slice];

  sigz=sqrt((*Beam_envelope)[4][4]);

  if(KO_SLICE) {
    KOslicing(n_slice,sigz,z_slice);}
  else{
    KHslicing(n_slice,sigz,z_slice);}   
   
  for(i=0;i<n_slice;i++) {
    np_slice[i]=N_particle/(double)n_slice;
  }
  double k=N_particle*r_e*m_e/2./pi/strong_beam_energy;
  sigx=sqrt((*Beam_envelope)[0][0]);
  sigy=sqrt((*Beam_envelope)[2][2]);
  xi_x=k*twiss[1]/(sigx*(sigx+sigy));
  xi_y=k*twiss[3]/(sigy*(sigx+sigy));

  luminosity=0.;
  dLum=0.;
}

IP::IP(double *K) : Element(K)
{
  Beam_envelope=new matrix(N_cv,N_cv);
  Benv_headon=new matrix(N_cv,N_cv);
  Benv_slice=new matrix(5,5);
//  n_slice=(int)K[4];
  np_slice=new double[max_n_slice];
  z_slice=new double[max_n_slice];

  this->IP_reset(K);
}

void IP::IP_reset(double *K)
{
  int i,j;
  length=0.;
  lmon=0;

  N_particle=K[0];
  x_angle=K[1];
  cosx=cos(x_angle);
  sinx=sin(x_angle);
  tanx=tan(x_angle);

  strong_beam_energy=K[2];
  weak_beam_energy=K[3];


  for(i=0;i<6;i++) cod[i]=K[i+6];
  if(x_angle!=0.) {
    cod[3]=cod[3]/cosx;
    cod[4]=cod[4]/cosx;
    cod[2]=cod[2]/cosx;
  }


   // Beam envelope


  int ij=0;
  for(i=0;i<6;i++) {
    for(j=i;j<6;j++) {
      (*Beam_envelope)[i][j]=K[12+ij];
      (*Beam_envelope)[j][i]=K[12+ij];
      ij++;
    }
  }
//  cout << "Beam envelope matrix \n" << *Beam_envelope;
//  cout.flush();
  // Crossing angle
  if(x_angle!=0.) {
    matrix mCros(N_cv,N_cv);
    mCros=set_Cros(x_angle);
    //     cout << "Crossing matrix \n" << mCros;
    *Benv_headon=((matrix)mCros)*(*Beam_envelope)*Transpose(mCros);
  }
//  cout << "Beam envelope matrix (Crossing)\n" << *Benv_headon;

   // Get 5x5 matrix
  tilt_angle= SliceBenv(Benv_headon,Benv_slice,bcen_fac);
  sin_t=sin(tilt_angle);
  cos_t=cos(tilt_angle);

// -------------------------------------------------------------------
// Beam slicing
// -------------------------------------------------------------------

  n_slice=(int)K[4];
  sigz=sqrt((*Benv_headon)[4][4]);

  if(KO_SLICE) {
    KOslicing(n_slice,sigz,z_slice);}
  else{
    KHslicing(n_slice,sigz,z_slice);}   

  for(i=0;i<n_slice;i++) {
    np_slice[i]=N_particle/(double)n_slice;
  }

  luminosity=0.;
  dLum=0.;
}



void BenvCalc(double emx,double emy,double emz,double* twiss,
	 matrix* Beam_envelope)
{
  matrix NtoP(N_cv,N_cv);
  NtoP=set_NtoP(twiss);
//  cout << "Normal to Physcal coord. matrix \n" << NtoP;

  matrix Emit(6,6);
  double emits[6];
  emits[0]=emx;emits[1]=emx;
  emits[2]=emy;emits[3]=emy;
  emits[4]=emz;emits[5]=emz;
  Emit.DiagonalMatrix(emits);
  *Beam_envelope=(matrix)NtoP*Emit*Transpose(NtoP);

  cout << "Beam envelope matrix \n" << *Beam_envelope;
}

double SliceBenv(matrix* Benv_headon,matrix* Benv_slice,double* bcen_fac)
{
  int i,j;
  double tilt_angle;

  matrix Beinv(N_cv,N_cv);
  Beinv=Inverse(*Benv_headon);
// cout << "Beinv matrix \n" << Beinv << "\ncheck\n" << Beinv*(*Benv_headon);
  matrix A(5,5),B(5,5);
  double v6[6],eig[10],eigr[5];
  m_zslice(Beinv,B,v6);
//   cout << "decomposed matrix\n" << A << '\n';
//   for(i=0;i<6;i++) cout << v6[i];
  teigen(B,A,eig);

  Normalize(A,eig);

  for(i=0;i<5;i++) {
    if(eig[2*i]<0.) {
      cout << "IP constructor : 5dim. beam envelope has zero or" << 
	" negative eigen values";
      exit(1);
    }
    eigr[i]=1./eig[2*i];
  }
  (*Benv_slice).DiagonalMatrix(eigr);
  *Benv_slice=(matrix)A*(*Benv_slice)*Transpose(A);
  for(i=0;i<5;i++) {
    bcen_fac[i]=0.;
    for(j=0;j<5;j++) bcen_fac[i]-=(*Benv_slice)[i][j]*v6[j];
  }
  tilt_angle=benv_norm_axis(*Benv_slice,bcen_fac);
  return tilt_angle;
}


void KOslicing(int n_slice,double sigz,double* z_slice)
{
  double z,f;
  int i,i0,n;
  n=(n_slice-1)/2;
//  cout << " KO slicing\n";
  if(n_slice%2) i0=0; else i0=1;
  for(i=i0;i<n_slice;i+=2) {
    z=0.;
    while(fabs(f=erf_fl(z)-(double)(i)/(double)n_slice)>1e-7) {
      z-=f*sqr_pi/2/exp(z*z);
    }
    z=z*sqrt(2.);
    if(n_slice%2) {
      z_slice[n+i/2]=-z*sigz;
      if(i!=0) z_slice[n-i/2]=z*sigz;
    }
    else {
      z_slice[n+i/2+1]=-z*sigz;
      z_slice[n-i/2]=z*sigz;
    }
  } 
}

void KHslicing(int n_slice,double sigz,double* z_slice)
{
  double z,f;
  int i,i0,n;
//  cout << " KH slicing\n";
  if(n_slice%2) i0=1; else i0=2;
  n=(n_slice-1)/2;
  double z_ba=0.,zc;
  for(i=i0;i<n_slice;i+=2) {
    z=0.;
    while(fabs(f=erf_fl(z)-(double)(i)/(double)n_slice)>1e-5) {
      z-=f*sqr_pi/2/exp(z*z);
    }
    z=z*sqrt(2.);
    if(i==2) z_ba=0.;
    zc=(exp(-z_ba*z_ba*0.5)-exp(-z*z*0.5))/sqr_pi/sqrt(2.);
    if(i==1) zc=0.;
    z_ba=z;
    if(n_slice%2) {
      z_slice[n+i/2]=-zc*sigz*n_slice;
      if(i!=1) z_slice[n-i/2]=zc*sigz*n_slice;
    }
    else {
      z_slice[n+i/2]=-zc*sigz*n_slice;
      z_slice[n-i/2+1]=zc*sigz*n_slice;
    }
  }
  z_slice[0]=exp(-z_ba*z_ba*0.5)/sqr_pi/sqrt(2.)*sigz*n_slice;
  z_slice[n_slice-1]=-z_slice[0];
}


IP::~IP(void)
{
   delete [] z_slice;
   delete [] np_slice;
   delete  Benv_slice;
   delete  Benv_headon;
   delete [] Beam_envelope;
}

void IP::print(void) 
{
  int i;
  cout << "#######  STRONG BEAM  #################################\n"
    << "N_particle of strong beam    " << N_particle
    << "\nEnergy of strong beam        " << strong_beam_energy
    << "\nEnergy of weak beam          " << weak_beam_energy
    << "\nCrossing angle               " << x_angle << "\n";
  cout << "Beam envelope matrix                \n" << *(Beam_envelope);
  cout << "Beam envelope matrix (Head-on frame)\n" << *(Benv_headon);
  cout << "On tilted frame   angle= " << tilt_angle << '\n' << *Benv_slice;
  for(i=0;i<5;i++) {cout.width(12); cout << bcen_fac[i] << ' ';}
  cout << '\n';
  cout << " \n Number of slice " << n_slice 
    << "\n       np_slice  z.slice " << '\n';
  for(i=0;i<n_slice;i++) { cout.width(5); cout << i+1;
    cout.width(10); cout << np_slice[i] << "  ";
    cout.width(10); cout << z_slice[i] << '\n';
  }
  cout << "\nsigma_x= " << sqrt((*Beam_envelope)[0][0]) 
    << ",  sigma_y= " << sqrt((*Beam_envelope)[2][2])
      << "\nsigma_z= " << sqrt((*Beam_envelope)[4][4]) 
	<< ",  sigma_e= " << sqrt((*Beam_envelope)[5][5]) << "\n";
//  cout << "\n Beam beam parameters \n";
//  cout << "xi_x= " << xi_x << "  xi_y= " << xi_y;
}


// ########################################################################
//         Mapping
// ########################################################################

void IP::Mapping(map_double& x)
{
  double sqrt_pi=sqrt(pi);
  double sqrt_2=sqrt(2.);
  double c0;  //double in all case!
//  following variables change for other types of mapping.
  Complex z1,z2,w1,w2,fz;
  double pxy2,pz,E,h,xx;
  double S,S2,Sigx,Sigy,Sigx_y,sgx,sgy,sqrti_Sig,c1,exy;
  double Unx,Uny,xy_fxy,dSigx_ds,dSigy_ds,g;
   
// %%%%% crossing angle correction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//  cout.setf(ios::scientific,ios::floatfield);
//  cout.precision(15);
//  cout << "Before BB " << x;
  if(x_angle!=0.) {
    pxy2=x[3]*x[3]+x[4]*x[4];
    E=1+x[5];
    pz=sqrt(E*E-pxy2);
    h=pxy2/(E+pz);
    x[5]=x[5]-tanx*x[3]+(tanx*tanx)*h;
    x[3]=(x[3]-tanx*h)/cosx;
    x[4]=x[4]/cosx;
    
    pxy2=x[3]*x[3]+x[4]*x[4];
    E=1+x[5];
    pz=pow(E*E-pxy2,-0.5); // Note!! This pz is 1/pz.
    //      h=h/(cosx*cosx);
    xx=x[0];
    x[0]=tanx*x[2]+(1.+sinx*x[3]*pz)*xx;
    x[1]=x[1]+sinx*x[4]*xx*pz;
    x[2]=x[2]/cosx-(sinx/cosx/cosx)*h*pz*xx;
    //cout.setf(ios::scientific,ios::floatfield);
    //cout.precision(15);
    //cout << "after lorentz boost\n" << x1;
  }
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  int icod=0,deform=0,i;
  for(i=0;i<5;i++) {if(cod[i]!=0.) icod=1;}
  for(i=0;i<4;i++) {if(bcen_fac[i]!=0.) deform=1;}

// cod correction
  if(icod) {
    x[0]-=cod[0];
    x[1]-=cod[1];
    x[2]+=cod[2];
    x[3]+=cod[3];
    x[4]+=cod[4];
  }
// tilt angle correction
  if(tilt_angle!=0.) {
    xx=cos_t*x[0]+sin_t*x[1];
    x[1]=cos_t*x[1]-sin_t*x[0];
    x[0]=xx;
    xx=cos_t*x[3]+sin_t*x[4];
    x[4]=cos_t*x[4]-sin_t*x[3];
    x[3]=xx;
  }
//
  for(i=0;i<n_slice;i++) {
    pxy2=x[3]*x[3]+x[4]*x[4];
    S=(x[2]-z_slice[i])*0.5;	// S=(z-z*)/2
    {
      // exp(-:D:)x		where D=(px^2+py^2)*S(z,z*)/2
      x[0]+=x[3]*S;		// X=x+px*S
      x[1]+=x[4]*S;		// Y=y+py*S
      x[5]-=pxy2*0.25;		// PZ=pz-(px^2+py^2)/4
    }
// strong beam off center from deformation.
    if(deform) {
      x[0]-=z_slice[i]*bcen_fac[0]; 
      x[1]-=z_slice[i]*bcen_fac[2]; 
      x[3]-=z_slice[i]*bcen_fac[1]; 
      x[4]-=z_slice[i]*bcen_fac[3]; 
    }
      
    // Beam-beam interaction  Body

    { // beam size
      S2=S*S;
      Sigx=(*Benv_slice)[0][0]-2.*(*Benv_slice)[0][1]*S+
	(*Benv_slice)[1][1]*S2;
      Sigy=(*Benv_slice)[2][2]-2.*(*Benv_slice)[2][3]*S+
	(*Benv_slice)[3][3]*S2;
      sgx=sqrt(Sigx);
      sgy=sqrt(Sigy);
      Sigx_y=Sigx-Sigy;
      sqrti_Sig=pow(Sigx_y,-0.5)/sqrt_2;
      
      c0=2.*np_slice[i]*r_e*m_e/weak_beam_energy;
      c1=(c0*sqrt_pi)*sqrti_Sig;
    }
    { // Px,Py kick
	 //cout << " c0 " << c0 << " (x,y) " << z1 << " S " << S <<'\n';
	 //cout << " sgx sgy  " << sgx << ' ' << sgy << '\n';
	 //cout << "Sigx_y " << Sigx_y << " sqrti_Sig " << sqrti_Sig << '\n';

      z1=Complex(x[0],x[1]);
      exy=exp(-x[0]*x[0]*0.5/Sigx-x[1]*x[1]*0.5/Sigy);
      z1=z1*sqrti_Sig;
      z2=Complex(sgy*x[0]/sgx,sgx*x[1]/sgy);
      z2=z2*sqrti_Sig;
      if(x[1]>0.) {
	w1=cerf(z1);
	w2=cerf(z2);
      }
      else {
	z1=conj(z1);
	w1=-conj(cerf(z1));
	z2=conj(z2);
	w2=-conj(cerf(z2));
      }
      //cout << " z1 and cerf(z1) " << z1 << ' ' << w1 << '\n';
      //cout << " z2 and cerf(z2) " << z2 << ' ' << w2 << '\n' 
      // << " c1 " << c1 << " exp(-x^2/Sig..) " << exy <<'\n';
      fz=c1*(w1-exy*w2);
      //      fx=imag(fz); fy=real(fz);
      x[3]-=imag(fz);
      x[4]-=real(fz);
    }
    {  // Pz kick
      dSigx_ds=(*Benv_slice)[1][1]*S-(*Benv_slice)[0][1];
      dSigy_ds=(*Benv_slice)[3][3]*S-(*Benv_slice)[2][3];
      xy_fxy=x[0]*imag(fz)+x[1]*real(fz);
      Unx=xy_fxy+c0*(sgy*exy/sgx-1.);
      Unx=Unx/Sigx_y*(-0.5);
      Uny=xy_fxy+c0*(sgx*exy/sgy-1.);
      Uny=Uny/Sigx_y*0.5;
//      g=0.5*(dSigx_ds*Unx+dSigy_ds*Uny);
      g=dSigx_ds*Unx+dSigy_ds*Uny;
      
      x[5]-=g;
    }
// strong beam off center from deformation.
    if(deform) {
      x[0]+=z_slice[i]*bcen_fac[0]; 
      x[1]+=z_slice[i]*bcen_fac[2]; 
      x[3]+=z_slice[i]*bcen_fac[1]; 
      x[4]+=z_slice[i]*bcen_fac[3]; 
    }

    //cout << " fx fy g " << imag(fz) << " " << real(fz) << " " << g << '\n';
    //cout.flush();
    {  // exp(:D:)x
      pxy2=x[3]*x[3]+x[4]*x[4];
      x[0]-=x[3]*S;
      x[1]-=x[4]*S;
      x[5]+=pxy2*0.25;
    }
  }
//
//  End slice loop
//

// tilt angle correction
  if(tilt_angle!=0.) {
    xx=cos_t*x[0]-sin_t*x[1];
    x[1]=cos_t*x[1]+sin_t*x[0];
    x[0]=xx;
    xx=cos_t*x[3]-sin_t*x[4];
    x[4]=cos_t*x[4]+sin_t*x[3];
    x[3]=xx;
  }

// cod correction
  if(icod) {
    x[0]+=cod[0];
    x[1]+=cod[1];
    x[2]-=cod[2];
    x[3]-=cod[3];
    x[4]-=cod[4];
  }
//
// %%%%% crossing angle correction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if(x_angle!=0.) {
    pxy2=x[3]*x[3]+x[4]*x[4];
    E=1+x[5];
    pz=sqrt(E*E-pxy2);
    h=pxy2/(E+pz);
    pz=1./pz;
    x[0]=(x[0]-sinx*x[2])/(1.+(sinx*(x[3]+sinx*h))*pz);
    x[1]=x[1]-sinx*x[4]*pz*x[0];
    x[2]=x[2]*cosx+(sinx*cosx)*h*pz*x[0];
    h=h*(cosx*cosx);
    xx=x[3];
    x[3]=xx*cosx+tanx*h;
    x[4]=x[4]*cosx;
    x[5]=x[5]+sinx*xx;
  }
//  cout.setf(ios::scientific,ios::floatfield);
//  cout.precision(15);
//  cout << "After BB " << x;
//  cout.flush();
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

}

void IP::Mapping(map_da& x)
{
  double sqrt_pi=sqrt(pi);
  double sqrt_2=sqrt(2.);
  double c0;  //double in all case!
//  following variables change for other types of mapping.
  c_da z1,w1,w2,fz,w3;
  da pxy2,pz,E,xx,h;
  da S,S2,Sigx,Sigy,Sigx_y,sgx,sgy,sqrti_Sig,c1,exy;
  da Unx,Uny,xy_fxy,dSigx_ds,dSigy_ds,g;
// %%%%% crossing angle correction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if(x_angle!=0.) {
    pxy2=x[3]*x[3]+x[4]*x[4];
    E=1+x[5];
    pz=sqrt(E*E-pxy2);
    h=pxy2/(E+pz);
    x[5]=x[5]-tanx*x[3]+(tanx*tanx)*h;
    x[3]=(x[3]-tanx*h)/cosx;
    x[4]=x[4]/cosx;

    pxy2=x[3]*x[3]+x[4]*x[4];
    E=1+x[5];
    pz=pow(E*E-pxy2,-0.5); // Note!! This pz is 1/pz.
    //      h=h/(cosx*cosx);
    xx=x[0];
    x[0]=tanx*x[2]+(1.+sinx*x[3]*pz)*xx;
    x[1]=x[1]+sinx*x[4]*xx*pz;
    x[2]=x[2]/cosx-(sinx/cosx/cosx)*h*pz*xx;
  }
  //   cout << x;
  //   is_symplectic(x);
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  int icod=0,deform=0,i;
  for(i=0;i<5;i++) {if(cod[i]!=0.) icod=1;}
  for(i=0;i<4;i++) {if(bcen_fac[i]!=0.) deform=1;}

// cod correction
  if(icod) {
    x[0]-=cod[0];
    x[1]-=cod[1];
    x[2]+=cod[2];
    x[3]+=cod[3];
    x[4]+=cod[4];
  }
// tilt angle correction
  if(tilt_angle!=0.) {
    xx=cos_t*x[0]+sin_t*x[1];
    x[1]=cos_t*x[1]-sin_t*x[0];
    x[0]=xx;
    xx=cos_t*x[3]+sin_t*x[4];
    x[4]=cos_t*x[4]-sin_t*x[3];
    x[3]=xx;
  }
//
 
  for(i=0;i<n_slice;i++) {
    pxy2=x[3]*x[3]+x[4]*x[4];
    S=(x[2]-z_slice[i])*0.5;	// S=(z-z*)/2
    {
      // exp(-:D:)x		where D=(px^2+py^2)*S(z,z*)/2
      x[0]+=x[3]*S;		// X=x+px*S
      x[1]+=x[4]*S;		// Y=y+py*S
      x[5]-=pxy2*0.25;		// PZ=pz-(px^2+py^2)/4
    }
// strong beam off center from deformation.
    if(deform) {
      x[0]-=z_slice[i]*bcen_fac[0]; 
      x[1]-=z_slice[i]*bcen_fac[2]; 
      x[3]-=z_slice[i]*bcen_fac[1]; 
      x[4]-=z_slice[i]*bcen_fac[3]; 
    }
      
    //is_symplectic(x);
    // Beam-beam interaction  Body
    //cout << " map before 4dim BB ----------------\n" << x;

      
    { // beam size
      //x1.dBase();
      S=(x[2]-z_slice[i])*0.5;	// S=(z-z*)/2
      S2=S*S;
      Sigx=(*Benv_slice)[0][0]-2.*(*Benv_slice)[0][1]*S+
	(*Benv_slice[1][1])*S2;
      Sigy=(*Benv_slice)[2][2]-2.*(*Benv_slice)[2][3]*S+
	(*Benv_slice)[3][3]*S2;
      //cout << "=========== TEST x1=unit map at bb.cc ==============\n";
      Sigx_y=Sigx-Sigy;
      sgx=sqrt(Sigx);
      sgy=sqrt(Sigy);
      sqrti_Sig=1./sqrt(Sigx_y)/sqrt_2;
      
      c0=2.*np_slice[i]*r_e*m_e/weak_beam_energy;
      c1=(c0*sqrt_pi)*sqrti_Sig;
    }
    { // Px,Py kick
      z1=c_da(x[0],x[1]);	// should be change when other types map.
      z1=z1*sqrti_Sig;
      //w1=cerf(z1);
      //w1=cerf_odd(z1);
      //cout << "z1 \n" << z1 ;
      exy=exp(-x[0]*x[0]*0.5/Sigx-x[1]*x[1]*0.5/Sigy);
      w3=exp(-z1*z1);
      z1=z1*Complex(0.,-1.);
      w1=erf(z1);
				// should be change when other types map.
      z1=c_da(sgy*x[0]/sgx,sgx*x[1]/sgy);	
      z1=Complex(0.,-1.)*z1*sqrti_Sig;
      //w2=cerf(z1);
      //w2=cerf_odd(z1);
      w2=erf(z1);
      //cout << "z2 w2\n" << z1 << w2; cout.flush();
      //fz=c1*(w1-exy*w2);
      fz=c1*w3*(w2-w1);
      //      fx=imag(fz); fy=real(fz);
      x[3]-=imag(fz);
      x[4]-=real(fz);
      //cout << "c1" << c1;
      //cout << "\n----- fz -----------------------------------------\n" << fz;
    }
    {  // Pz kick
      dSigx_ds=(*Benv_slice)[1][1]*S-(*Benv_slice)[0][1];
      dSigy_ds=(*Benv_slice)[3][3]*S-(*Benv_slice)[2][3];
      xy_fxy=x[0]*imag(fz)+x[1]*real(fz);
      Unx=xy_fxy+c0*(sgy*exy/sgx-1.);
      Unx=Unx/Sigx_y*(-0.5);
      Uny=xy_fxy+c0*(sgx*exy/sgy-1.);
      Uny=Uny/Sigx_y*0.5;
      g=dSigx_ds*Unx+dSigy_ds*Uny;
      //cout << "xy_fxy" << xy_fxy;
      //cout << "dSigx_ds  dSigy_ds" << dSigx_ds << dSigy_ds;
      //cout << "Unx  Uny " << Unx << Uny;
      //cout << "\n------ g ------------------------------------------\n" << g;
      
      x[5]-=g;
      //cout <<x1;
      //is_symplectic(x1); cout.flush();
      //x1=x1*x2;
      //cout << "concatenation 4dimBB and zmap\n\n" << x1; 
      //is_symplectic(x1); cout.flush();
      
      
    }
// strong beam off center from deformation.
    if(deform) {
      x[0]+=z_slice[i]*bcen_fac[0]; 
      x[1]+=z_slice[i]*bcen_fac[2]; 
      x[3]+=z_slice[i]*bcen_fac[1]; 
      x[4]+=z_slice[i]*bcen_fac[3]; 
    }
    {  // exp(:D:)x
      pxy2=x[3]*x[3]+x[4]*x[4];
      x[0]-=x[3]*S;
      x[1]-=x[4]*S;
      x[5]+=pxy2*0.25;
    }
  }
//
//  End slice loop
//
   
// tilt angle correction
  if(tilt_angle!=0.) {
    xx=cos_t*x[0]-sin_t*x[1];
    x[1]=cos_t*x[1]+sin_t*x[0];
    x[0]=xx;
    xx=cos_t*x[3]-sin_t*x[4];
    x[4]=cos_t*x[4]+sin_t*x[3];
    x[3]=xx;
  }

// cod correction
  if(icod) {
    x[0]+=cod[0];
    x[1]+=cod[1];
    x[2]-=cod[2];
    x[3]-=cod[3];
    x[4]-=cod[4];
  }
//
// %%%%% crossing angle correction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //   x1.dBase();
  if(x_angle!=0.) {
    pxy2=x[3]*x[3]+x[4]*x[4];
    E=1+x[5];
    pz=sqrt(E*E-pxy2);
    h=pxy2/(E+pz);
    pz=1./pz;
    x[0]=(x[0]-sinx*x[2])/(1.+(sinx*(x[3]+sinx*h))*pz);
    x[1]=x[1]-sinx*x[4]*pz*x[0];
    x[2]=x[2]*cosx+(sinx*cosx)*h*pz*x[0];
    h=h*(cosx*cosx);
    xx=x[3];
    x[3]=xx*cosx+tanx*h;
    x[4]=x[4]*cosx;
    x[5]=x[5]+sinx*xx;
  }
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
}


void IP::Mapping(map_p_da& x)
{
  double sqrt_pi=sqrt(pi);
  double sqrt_2=sqrt(2.);
  double c0;  //double in all case!
//  following variables change for other types of mapping.
  c_p_da z1,w1,w2,fz;
  p_da pxy2,E,pz,h,xx;
  p_da S,S2,Sigx,Sigy,Sigx_y,sgx,sgy,sqrti_Sig,c1,exy;
  p_da Unx,Uny,xy_fxy,dSigx_ds,dSigy_ds,g;
      
// %%%%% crossing angle correction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if(x_angle!=0.) {
    pxy2=x[3]*x[3]+x[4]*x[4];
    E=1+x[5];
    pz=sqrt(E*E-pxy2);
    h=pxy2/(E+pz);
    x[5]=x[5]-tanx*x[3]+(tanx*tanx)*h;
    x[3]=(x[3]-tanx*h)/cosx;
    x[4]=x[4]/cosx;

    pxy2=x[3]*x[3]+x[4]*x[4];
    E=1+x[5];
    pz=pow(E*E-pxy2,-0.5); // Note!! This pz is 1/pz.
    //      h=h/(cosx*cosx);
    xx=x[0];
    x[0]=tanx*x[2]+(1.+sinx*x[3]*pz)*xx;
    x[1]=x[1]+sinx*x[4]*xx*pz;
    x[2]=x[2]/cosx-(sinx/cosx/cosx)*h*pz*xx;
  }
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  int icod=0,deform=0,i;
  for(i=0;i<5;i++) {if(cod[i]!=0.) icod=1;}
  for(i=0;i<4;i++) {if(bcen_fac[i]!=0.) deform=1;}

// cod correction
  if(icod) {
    x[0]-=cod[0];
    x[1]-=cod[1];
    x[2]+=cod[2];
    x[3]+=cod[3];
    x[4]+=cod[4];
  }
// tilt angle correction
  if(tilt_angle!=0.) {
    xx=cos_t*x[0]+sin_t*x[1];
    x[1]=cos_t*x[1]-sin_t*x[0];
    x[0]=xx;
    xx=cos_t*x[3]+sin_t*x[4];
    x[4]=cos_t*x[4]-sin_t*x[3];
    x[3]=xx;
  }
//
   
  for(i=0;i<n_slice;i++) {
    pxy2=x[3]*x[3]+x[4]*x[4];
    S=(x[2]-z_slice[i])*0.5;	// S=(z-z*)/2
    {
      // exp(-:D:)x		where D=(px^2+py^2)*S(z,z*)/2
      x[0]+=x[3]*S;		// X=x+px*S
      x[1]+=x[4]*S;		// Y=y+py*S
      x[5]-=pxy2*0.25;		// PZ=pz-(px^2+py^2)/4
    }
// strong beam off center from deformation.
    if(deform) {
      x[0]-=z_slice[i]*bcen_fac[0]; 
      x[1]-=z_slice[i]*bcen_fac[2]; 
      x[3]-=z_slice[i]*bcen_fac[1]; 
      x[4]-=z_slice[i]*bcen_fac[3]; 
    }
      
    // Beam-beam interaction  Body

    { // beam size
      S2=S*S;
      Sigx=(*Benv_slice)[0][0]-2.*(*Benv_slice)[0][1]*S+
	(*Benv_slice)[1][1]*S2;
      Sigy=(*Benv_slice)[2][2]-2.*(*Benv_slice)[2][3]*S+
	(*Benv_slice)[3][3]*S2;
      Sigx_y=Sigx-Sigy;
      sgx=sqrt(Sigx);
      sgy=sqrt(Sigy);
      sqrti_Sig=pow(Sigx_y,-0.5)/sqrt_2;
      
      c0=2.*np_slice[i]*r_e*m_e/weak_beam_energy;
      c1=(c0*sqrt_pi)*sqrti_Sig;
    }
    { // Px,Py kick
      z1=c_p_da(x[0],x[1]);
      z1=z1*sqrti_Sig;
      w1=cerf(z1);
      exy=exp(-x[0]*x[0]*0.5/Sigx-x[1]*x[1]*0.5/Sigy);
      z1=c_p_da(sgy*x[0]/sgx,sgx*x[1]/sgy);
      z1=z1*sqrti_Sig;
      w2=cerf(z1);
      fz=c1*(w1-exy*w2);
      //      fx=imag(fz); fy=real(fz);
      x[3]-=imag(fz);
      x[4]-=real(fz);
    }
    {  // Pz kick
      dSigx_ds=(*Benv_slice)[1][1]*S-(*Benv_slice)[0][1];
      dSigy_ds=(*Benv_slice)[3][3]*S-(*Benv_slice)[2][3];
      xy_fxy=x[0]*imag(fz)+x[1]*real(fz);
      Unx=xy_fxy+c0*(sgy*exy/sgx-1.);
      Unx=Unx/Sigx_y*(-0.5);
      Uny=xy_fxy+c0*(sgx*exy/sgy-1.);
      Uny=Uny/Sigx_y*0.5;
      g=dSigx_ds*Unx+dSigy_ds*Uny;
      
      x[5]-=g;
    }
// strong beam off center from deformation.
    if(deform) {
      x[0]+=z_slice[i]*bcen_fac[0]; 
      x[1]+=z_slice[i]*bcen_fac[2]; 
      x[3]+=z_slice[i]*bcen_fac[1]; 
      x[4]+=z_slice[i]*bcen_fac[3]; 
    }
    {  // exp(:D:)x
      pxy2=x[3]*x[3]+x[4]*x[4];
      x[0]-=x[3]*S;
      x[1]-=x[4]*S;
      x[5]+=pxy2*0.25;
    }
  }
//
//  End slice loop
//
// tilt angle correction
  if(tilt_angle!=0.) {
    xx=cos_t*x[0]-sin_t*x[1];
    x[1]=cos_t*x[1]+sin_t*x[0];
    x[0]=xx;
    xx=cos_t*x[3]-sin_t*x[4];
    x[4]=cos_t*x[4]+sin_t*x[3];
    x[3]=xx;
  }
//

// cod correction
  if(icod) {
    x[0]+=cod[0];
    x[1]+=cod[1];
    x[2]-=cod[2];
    x[3]-=cod[3];
    x[4]-=cod[4];
  }
   
// %%%%% crossing angle correction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if(x_angle!=0.) {
    pxy2=x[3]*x[3]+x[4]*x[4];
    E=1+x[5];
    pz=sqrt(E*E-pxy2);
    h=pxy2/(E+pz);
    pz=1./pz;
    x[0]=(x[0]-sinx*x[2])/(1.+(sinx*(x[3]+sinx*h))*pz);
    x[1]=x[1]-sinx*x[4]*pz*x[0];
    x[2]=x[2]*cosx+(sinx*cosx)*h*pz*x[0];
    h=h*(cosx*cosx);
    xx=x[3];
    x[3]=xx*cosx+tanx*h;
    x[4]=x[4]*cosx;
    x[5]=x[5]+sinx*xx;
  }
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

}



// pBeam

#define MX 40
#define MY 31

static Complex W_[4][MX][MY];

void Seterf()
{
  Complex Z_,CWERF_;
  int I,J,NX,NY;
  double X,Y,DX,DY;

  Complex I_(0.,1.),C_1(-1.5,1.5),C_2(0.5,-0.5);
  Complex C_3(0,-2),C_4(1.5,0.5),C_5(-1.5,0.5);
  DX=0.1;
  DY=0.1;
  NX=39;
  NY=30;
  for(I=0;I<=NX;I++) {
    X=DX*I;
    for(J=0;J<=NY;J++) {
      Y=DY*J;
      Z_=Complex(X,Y);
      W_[0][I][J]=cerf(Z_);
    }
  }
  for(I=0;I<NX;I++) {
    for(J=0;J<NY;J++) {
      W_[1][I][J]=C_1*W_[0][I][J]+C_2*W_[0][I+1][J+1]+W_[0][I][J+1]
	-I_*W_[0][I+1][J];
      W_[2][I][J]=C_3*W_[0][I][J]+C_4*W_[0][I+1][J]+I_*W_[0][I+1][J+1]
	+C_5*W_[0][I][J+1];
      W_[3][I][J]=C_2*(W_[0][I][J+1]-W_[0][I+1][J]
		       +I_*(W_[0][I][J]-W_[0][I+1][J+1]));
    }
  }
}

void IP::Mapping(pBeam& x)
{
  int j;
  double sqrt_pi=sqrt(pi);
  double cn0;  //double in all case!
//  following variables change for other types of mapping.
  double fx0,fy0,d1x,d1y,asp,cdu;
  double pxy2,pz,E,h,xx;
  double S,S2,sig11x,sig11y,Sigx_y,sqrsig11i,cne,a2,a2x,a2y,hxy;
  double dux,duy,fxy,sig12x,sig12y,gxy;
// Tab -----------------------------------------------------
  double X,Y,Z2I,Z2R;
  int IX,IY;
  Complex DZ_,Z_1;
  Complex WERF1,WERF2;
  const double C1=0.4613135E0,C2=0.1901635E0,C3=0.09999216E0,C4=1.7844927E0;
  const double C5=2.883894E-3,C6=5.5253437E0,C7=0.5124242E0, C8=0.2752551E0;
  const double C9=0.05176536E0,C10=2.724745E0;
// Tab -----------------------------------------------------
  dLum=0.;
// %%%%% crossing angle correction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//  cout.setf(ios::scientific,ios::floatfield);
//  cout.precision(15);
//  cout << "Before BB " << x1;
  if(x_angle!=0.) {
    for(j=0;j<x.np;j++) {
      pxy2=x.px[j]*x.px[j]+x.py[j]*x.py[j];
      E=1+x.pz[j];
      pz=sqrt(E*E-pxy2);
      h=pxy2/(E+pz);
      x.pz[j]=x.pz[j]-tanx*x.x[j]+(tanx*tanx)*h;
      x.px[j]=(x.px[j]-tanx*h)/cosx;
      x.py[j]=x.py[j]/cosx;

      pxy2=x.px[j]*x.px[j]+x.py[j]*x.py[j];
      E=1+x.pz[j];
      pz=pow(E*E-pxy2,-0.5); // Note!! This pz is 1/pz.
      //      h=h/(cosx*cosx);
      xx=x.x[j];
      x.x[j]=tanx*x.z[j]+(1.+sinx*x.px[j]*pz)*xx;
      x.y[j]=x.y[j]+sinx*x.py[j]*xx*pz;
      x.z[j]=x.z[j]/cosx-(sinx/cosx/cosx)*h*pz*xx;
      //cout.setf(ios::scientific,ios::floatfield);
      //cout.precision(15);
      //cout << "after lorentz boost\n" << x1;
    }
  }
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  int icod=0,deform=0,i;
  for(i=0;i<5;i++) {if(cod[i]!=0.) icod=1;}
  for(i=0;i<4;i++) {if(bcen_fac[i]!=0.) deform=1;}

// cod correction
  if(icod) {
    for(j=0;j<x.np;j++) {
      x.x[j]-=cod[0];
      x.y[j]-=cod[1];
      x.z[j]+=cod[2];
      x.px[j]+=cod[3];
      x.py[j]+=cod[4];
    }
  }
// tilt angle correction
  if(tilt_angle!=0.) {
    for(j=0;j<x.np;j++) {
      xx=cos_t*x.x[j]+sin_t*x.y[j];
      x.y[j]=cos_t*x.y[j]-sin_t*x.x[j];
      x.x[j]=xx;
      xx=cos_t*x.px[j]+sin_t*x.py[j];
      x.py[j]=cos_t*x.py[j]-sin_t*x.px[j];
      x.px[j]=xx;
    }
  }


//

  for(i=0;i<n_slice;i++) {
    for(j=0;j<x.np;j++) {
      pxy2=x.px[j]*x.px[j]+x.py[j]*x.py[j];
      S=(x.z[j]-z_slice[i])*0.5;	// S=(z-z*)/2
      {
      // exp(-:D:)x		where D=(px^2+py^2)*S(z,z*)/2
	x.x[j]+=x.px[j]*S;		// X=x+px*S
	x.y[j]+=x.py[j]*S;		// Y=y+py*S
	x.pz[j]-=pxy2*0.25;		// PZ=pz-(px^2+py^2)/4
      }
// strong beam off center from deformation.
      if(deform) {
	x.x[j]-=z_slice[i]*bcen_fac[0]; 
	x.y[j]-=z_slice[i]*bcen_fac[2]; 
	x.px[j]-=z_slice[i]*bcen_fac[1]; 
	x.py[j]-=z_slice[i]*bcen_fac[3]; 
      }
      
      // Beam-beam interaction  Body

      { // beam size
	S2=S*S;
	sig11x=(*Benv_slice)[0][0]-2.*(*Benv_slice)[0][1]*S+
	  (*Benv_slice)[1][1]*S2;
	sig11y=(*Benv_slice)[2][2]-2.*(*Benv_slice)[2][3]*S+
	  (*Benv_slice)[3][3]*S2;
	sig12x=(*Benv_slice)[1][1]*S-(*Benv_slice)[0][1];
	sig12y=(*Benv_slice)[3][3]*S-(*Benv_slice)[2][3];
	Sigx_y=sig11x-sig11y;
	asp=sqrt(sig11y/sig11x);
	sqrsig11i=1./sqrt(2.*fabs(Sigx_y));
      
	cn0=2.*np_slice[i]*r_e*m_e/weak_beam_energy;
	cne=(cn0*sqrt_pi)*sqrsig11i;
      }
      { // Px,Py kick
	 //cout << " cn0 " << cn0 << " (x,y) " << z1 << " S " << S <<'\n';
	 //cout << " sig11x sig11y  " << sig11x << ' ' << sig11y << '\n';
	 //cout << "Sigx_y " << Sigx_y << " sqrsig11i " << sqrsig11i << '\n';

	d1x=fabs(x.x[j])*sqrsig11i;
	d1y=fabs(x.y[j])*sqrsig11i;
	X=d1x;
	Y=d1y;
	if(X<3.8999 && Y<2.9999) {
	  X=10.0*X;
	  Y=10.0*Y;
	  IX=int(X);
	  IY=int(Y);
	  DZ_=Complex(X-(double)(IX),Y-(double)(IY));
	  WERF1=((DZ_*W_[3][IX][IY]+W_[2][IX][IY])*DZ_
		 +W_[1][IX][IY])*DZ_+W_[0][IX][IY];
	}
	else if(X<6.0 && Y<6.0) {
	  Z2R=X*X-Y*Y;
	  Z2I=2.0*X*Y;
	  Z_1=Complex(-Y,X);
	  WERF1=Z_1*(C1/Complex(Z2R-C2,Z2I)+C3/Complex(Z2R-C4,Z2I)
		     +C5/Complex(Z2R-C6,Z2I));
	}
	else {
	  Z2R=X*X-Y*Y;
	  Z2I=2.0*X*Y;
	  Z_1=Complex(-Y,X);
	  WERF1=Z_1*(C7/Complex(Z2R-C8,Z2I)+C9/Complex(Z2R-C10,Z2I));
	}

	hxy=-x.x[j]*x.x[j]*0.5/sig11x-x.y[j]*x.y[j]*0.5/sig11y;

	X=d1x*asp;
	Y=d1y/asp;
	if(X<3.8999 && Y<2.9999) {
	  X=10.0*X;
	  Y=10.0*Y;
	  IX=int(X);
	  IY=int(Y);
	  DZ_=Complex(X-(double)(IX),Y-(double)(IY));
	  WERF2=((DZ_*W_[3][IX][IY]+W_[2][IX][IY])*DZ_
		 +W_[1][IX][IY])*DZ_+W_[0][IX][IY];
	}
	else if(X<6.0 && Y<6.0) {
	  Z2R=X*X-Y*Y;
	  Z2I=2.0*X*Y;
	  Z_1=Complex(-Y,X);
	  WERF2=Z_1*(C1/Complex(Z2R-C2,Z2I)+C3/Complex(Z2R-C4,Z2I)
		     +C5/Complex(Z2R-C6,Z2I));
	}
	else {
	  Z2R=X*X-Y*Y;
	  Z2I=2.0*X*Y;
	  Z_1=Complex(-Y,X);
	  WERF2=Z_1*(C7/Complex(Z2R-C8,Z2I)+C9/Complex(Z2R-C10,Z2I));
	}

	  
	a2=exp(hxy);
// %%%%%% dLuminosity  %%%%%%%%%%%%%%%%%%%
	dLum+=a2/pi2/sqrt(sig11x*sig11y);
//
	a2x=a2*imag(WERF2);
	a2y=a2*real(WERF2);
	fx0=cne*(imag(WERF1)-a2x);
	fy0=cne*(real(WERF1)-a2y);
	if(x.x[j]<0.) fx0=-fx0;
	if(x.y[j]<0.) fy0=-fy0;
        x.px[j]-=fx0;
        x.py[j]-=fy0;

      //cout << " z1 and cerf(z1) " << z1 << ' ' << w1 << '\n';
      //cout << " z2 and cerf(z2) " << z2 << ' ' << w2 << '\n' 
      // << " c1 " << c1 << " exp(-x^2/Sig..) " << a2 <<'\n';
      }
      {  // Pz kick
	cdu=-0.5/Sigx_y;
	fxy=x.x[j]*fx0+x.y[j]*fy0;
	dux=cdu*(fxy+cn0*(a2*asp-1.));
	duy=-cdu*(fxy+cn0*(a2/asp-1.));
	gxy=sig12x*dux+sig12y*duy;
      
	x.pz[j]-=gxy;
      }
// strong beam off center from deformation.
      if(deform) {
	x.x[j]+=z_slice[i]*bcen_fac[0]; 
	x.y[j]+=z_slice[i]*bcen_fac[2]; 
	x.px[j]+=z_slice[i]*bcen_fac[1]; 
	x.py[j]+=z_slice[i]*bcen_fac[3]; 
      }

      //cout << " fx fy g " << fx0 << " " << fy0 << " " << g << '\n';
      //cout.flush();
      {  // exp(:D:)x
	pxy2=x.px[j]*x.px[j]+x.py[j]*x.py[j];
	x.x[j]-=x.px[j]*S;
	x.y[j]-=x.py[j]*S;
	x.pz[j]+=pxy2*0.25;
      }
    }
  }
//
//  End slice loop
//

// tilt angle correction
  if(tilt_angle!=0.) {
    for(j=0;j<x.np;j++) {
      xx=cos_t*x.x[j]-sin_t*x.y[j];
      x.y[j]=cos_t*x.y[j]+sin_t*x.x[j];
      x.x[j]=xx;
      xx=cos_t*x.px[j]-sin_t*x.py[j];
      x.py[j]=cos_t*x.py[j]+sin_t*x.px[j];
      x.px[j]=xx;
    }
  }

// cod correction
  if(icod) {
    for(j=0;j<x.np;j++) {
      x.x[j]+=cod[0];
      x.y[j]+=cod[1];
      x.z[j]-=cod[2];
      x.px[j]-=cod[3];
      x.py[j]-=cod[4];
    }
  }
//
// %%%%% crossing angle correction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if(x_angle!=0.) {
    for(j=0;j<x.np;j++) {
      pxy2=x.px[j]*x.px[j]+x.py[j]*x.py[j];
      E=1+x.pz[j];
      pz=sqrt(E*E-pxy2);
      h=pxy2/(E+pz);
      pz=1./pz;
      x.x[j]=(x.x[j]-sinx*x.z[j])/(1.+(sinx*(x.px[j]+sinx*h))*pz);
      x.y[j]=x.y[j]-sinx*x.py[j]*pz*x.x[j];
      x.z[j]=x.z[j]*cosx+(sinx*cosx)*h*pz*x.x[j];
      h=h*(cosx*cosx);
      xx=x.px[j];
      x.px[j]=xx*cosx+tanx*h;
      x.py[j]=x.py[j]*cosx;
      x.pz[j]=x.pz[j]+sinx*xx;
    }
  }
//
// %%%%  Luminosity  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
  dLum=dLum/n_slice/x.np;
  if(lmon) luminosity+=dLum;

//  cout.setf(ios::scientific,ios::floatfield);
//  cout.precision(15);
//  cout << "After BB " << x1;
//  cout.flush();
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

}

