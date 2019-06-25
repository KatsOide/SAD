//---------------------------------------------------------------------------------------------
// BBBREM — Monte Carlo simulation of radiative Bhabha scattering in the very forward direction
// R. Kleiss, H. Burkhardt http://dx.doi.org/10.1016/0010-4655(94)90085-X
// C++ version
//---------------------------------------------------------------------------------------------

#ifndef BBBrem_h
#define BBBrem_h 1

#include <random>

class BBBrem {

public:

  BBBrem(double roots=91.2, double k0=0.01, unsigned int seed=0, double cutoff=0, bool write_evts_fl=false, bool write_evts_binary_fl=false,
         bool UnWeighted=false, double Unweighted_up_to=100,  bool WriteInvariants=true, unsigned int verbose=0);
  ~BBBrem()=default; // default destructor
  void generate_UnWeighted();
  void generate();
  void finish();
  void read_evts_binary(const std::string fname="bbrems.MyBinNtuple"); // example, how to read back optionally generated binary

  unsigned int Verbose() { return verbose; } // get
  void Verbose(unsigned int verbose) { this->verbose=verbose; } // set

protected:
  void setup_write_evts(int prec=6); // prec is number of digits for text output
  void        write_evts();  // text -- useful to see results for few events but slower than generation
  void write_evts_binary();  // binary output with minimal text header

  double s,t,d1,d2,s1; // invariants

  double sigapp,wmax,wnil,wneg,w0,w1,w2,sumwcut; // results for cross section calculation and statistics at end of generation

  double p1[4],p2[4],q1[4],q2[4],qk[4],weight,weight_remaining; // in and outgoing four vectors and weight
  double eb,pb,k0; // beam energy and momentum, from roots and me  and photon cutoff

  double me2s,logme2s,rin2pb; // constants at given roots : me2/s and log

  double z0,a1,ac; // k0/(1-k0) and log
  double sqrtmt,cutoff; // sqrt(-t) and optional cut in sqrtmt

  unsigned int seed;
  std::mt19937_64 mt; // Mersenne Twister https://en.wikipedia.org/wiki/Mersenne_Twister engine http://www.cplusplus.com/reference/random/mt19937_64
  std::uniform_real_distribution<double> RanD; // flat double in default range 0 to 1

  bool Make4Vec, write_evts_fl, write_evts_binary_fl, UnWeighted, WriteInvariants;
  double Unweighted_up_to;
  std::string fname,fname_bin;
  std::ofstream my_out_bin,my_out;

  unsigned int verbose;
};
#endif
