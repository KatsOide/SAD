//---------------------------------------------------------------------------------------------
// BBBREM — Monte Carlo simulation of radiative Bhabha scattering in the very forward direction
// R. Kleiss, H. Burkhardt http://dx.doi.org/10.1016/0010-4655(94)90085-X
// C++ version
//---------------------------------------------------------------------------------------------

#include <iostream>
#include <fstream>
#include <iomanip>

#include <string>
#include <sstream>

#include "RunFlags.h"
#include "BBBrem.h"

int main(int argc, char *argv[])
{
  std::vector<std::string> ArgvDefault={"45.6","0.01","10000","0","0"}; // default input parameters
  const char* cflags="bhtuvV";
  RunFlags Flag(argc,argv,cflags,ArgvDefault);

  double ebeam=std::stod(Flag.ArgvVec[1]);
  double k0   =std::stod(Flag.ArgvVec[2]);
  unsigned int nevent=std::stod(Flag.ArgvVec[3]); // double to string to also allow for 1.e5 events
  double seed  =std::stod(Flag.ArgvVec[4]);
  double cutoff=std::stod(Flag.ArgvVec[5]);

  // Flag.verbose=3; // extra testing
  if (Flag.verbose>2 && nevent>100 ) nevent=10;

  double roots=2*ebeam;
  bool Make4Vec= Flag.b||Flag.t; // generation of 4 vectors only needed when writing events

  if (Flag.verbose>1)
  {
    std::cout << __FILE__ << " line " << __LINE__ << " " << __FUNCTION__
    << " ebeam=" << ebeam
    << " k0=" << k0
    << " nevent="   << nevent
    << " cutoff="   << cutoff
    << " Make4Vec=" << Make4Vec
    << " verbose="  << Flag.verbose << '\n';
  }

  double Unweighted_up_to=100;
  bool WriteInvariants=true;

  if (Flag.verbose) std::cout << "BBBbrem start initialization" << '\n';

  BBBrem bbbrem(roots, k0, seed, cutoff, Flag.t, Flag.b, Flag.u, Unweighted_up_to, WriteInvariants, Flag.verbose); // construct BBBrem with initialization

  if (Flag.verbose) std::cout << "BBBbrem done with initialization. Start generation for ebeam=" << ebeam << " roots=" << roots << " GeV" << '\n';

  for(unsigned int k=0; k<nevent; ++k) // event loop
  {
    if (Flag.u) bbbrem.generate_UnWeighted();
    else bbbrem.generate(); // generate an event with weights  const char* fname="bbrems_ntuple.out";
  }

  if (Flag.verbose) std::cout << "BBBbrem done with generation" << '\n';

  bbbrem.finish();
}
