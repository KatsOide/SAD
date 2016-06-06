#include <fstream>
using std::ofstream;
using std::ios;
#include <cstdlib>
using std::atoi;
using std::exit;

#include <sadplus.h>

int main(int argc,char* argv[])
{
  int i;

  EMIT SAD;
  Beam BEAM;
  Accelerator TMP;
  int Noda=0;

  if(argc<1) exit(1);
  switch(argc) { // usage: $0 [[source file] order]
  case 2:
    Noda=atoi(argv[1]);
  case 1:
    OperatingParameterSet("fort.11",BEAM,SAD,TMP);
    break;

  default:
    Noda=atoi(argv[2]);
    OperatingParameterSet(argv[1],BEAM,SAD,TMP);
  }
  if(Noda<1) Noda=4;

//  Print(TMP);
//  cout << TMP;
  cout << "Line Length = " << Length(TMP) << '\n';

  InitializeDifferentialAlgebra(Noda);

  //  da x1;
  //x1.dBase(3,1.);
  //cout << x1;
//  exit(0);
  map_da x;
  x.dBase();

  TMP.Mapping(x);
  cout << x;
//  is_symplectic(x);
//  lin_map T;
//  is_symplectic(x);
//  T=lin_da(x);
//  PrintList((matrix)T);
  exit(0);
//  T.EigenSystem();
//  PrintTwiss(T);
  cout << SAD;
  da z;
//  z=fac_drg_type1(x);
//  z=fac_drg_type2(x);
//  z=fac_drg_type3(x);
  z=can_perturbation(x);

  ofstream daout("dap.out",ios::out);
  z.NFPrint(cout);
  z.NFPrint(daout);
  daout.close();
  exit(0);
//-----------------------------------------------------------------
//  map_double
//-----------------------------------------------------------------
  int nturn=5;
  map_double xd0,xd;
  xd0=0.;
  xd0[0]=1.1856E-5; 
  xd0[3]=7.6777E-6;
  xd0[1]=8.1925E-7;
  xd0[4]=1.6980E-5;
  xd0[2]=6.1386E-6;
  xd0[5]=3.5212E-4;
  ofstream trout("track.out",ios::out);
  for(i=0;i<nturn;i++) {
    xd=xd0;
    TMP.Mapping(xd);
//    cout << i << "th -turn\n" <<xd;
//    SAD.SynchrotronRadiation(xd);
//    cout <<xd <<'\n';
    trout << xd;
    cout << "start\n" << xd0 << "End \n" << xd << "\n";
    xd0=(xd+xd0)*0.5;
  }
  cout << "Final results\n" << xd;
  trout.close();
  exit(0);
//----------------------------------------------------------------
// pBeam
//----------------------------------------------------------------
  int np=100;

  pBeam xx(np);

//  matrix Benv=BeamEnvelope(BEAM);

  xx.Initialize("Gaussian",*SAD.Beam_envelope);
//  xx=0.;
//  xx[0][0]=0.0003;
//  xx[1][0]=0.00003;
//  for(i=0;i<np;i++) {
//    xx[0][i]=0.0001*(i+1);
//    xx[1][i]=0.00001*(i+1);
//  }
//  ofstream trout("track.out",ios::out);
  //   cout << xx; cout.flush();
//  trout << xx;
  for(i=0;i<nturn;i++) {
    if(i==0) TMP.LuminosityMonitorStart();
    TMP.Mapping(xx);
    SAD.SynchrotronRadiation(xx);
//    trout << xx;
  }
  cout << "Turn" << i <<"\n" << xx;
  cout << "Luminosity=" << TMP.Luminosity()/nturn << '\n';
//  trout.close();
  exit(0);
//  is_symplectic(x);
//  exit(0);
//  da z;
//  z=fac_drg_type1(x);
//  z=fac_drg_type2(x);
//  z=fac_drg_type3(x);
//  z=can_perturbation(x);

//  ofstream daout("dap.out",ios::out);
//  cout << Normal_expression(z);
//  daout.close();
}



