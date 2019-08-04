//---------------------------------------------------------------------------------------------
// BBBREM — Monte Carlo simulation of radiative Bhabha scattering in the very forward direction
// R. Kleiss, H. Burkhardt http://dx.doi.org/10.1016/0010-4655(94)90085-X
// C++ version
//---------------------------------------------------------------------------------------------

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>

#include "BBBrem.h"

// physical constants http://pdg.lbl.gov/2016/reviews/rpp2016-rev-phys-constants.pdf
double const alpha = 0.0072973525664;
double const me    = 0.5109989461e-03 , me2=me*me; // electron mass me and me^2
double const tomb  = 3.893793656e05 /1e6; // hbarc^2

/*
 // physical constants 1994 values
 double const alpha = 1./137.036;
 double const me    = 0.51099906e-03 , me2=me*me; // electron mass me and me^2
 double const tomb  = 3.8937966e05 /1e6; // hbarc^2
 */

double inline FourVecProd(const double p1[],const double p2[])
{
  return p1[0]*p2[0] -p1[1]*p2[1] -p1[2]*p2[2] -p1[3]*p2[3];
}

BBBrem::BBBrem(double roots, double k0, unsigned int seed, double cutoff, bool write_evts_fl, bool write_evts_binary_fl, bool UnWeighted, double Unweighted_up_to,
               bool WriteInvariants, unsigned int verbose) // constructor
: wmax(0), wnil(0), wneg(0), w0(0), w1(0), w2(0), sumwcut(0), // initialize counters
k0(k0), cutoff(cutoff), seed(seed), write_evts_fl(write_evts_fl), write_evts_binary_fl(write_evts_binary_fl),
UnWeighted(UnWeighted), WriteInvariants(WriteInvariants), Unweighted_up_to(Unweighted_up_to), verbose(verbose)  // initialize parameters
{
  
  mt.seed(seed);
  
  if (verbose>1) std::cout << __FILE__ << " line " << __LINE__ << " " << __PRETTY_FUNCTION__ << " constructor start"
  << " seed=" << this->seed << " verbose=" << verbose << '\n';
  
  std::cout
  << " ******************************************************" << '\n'
  << " *** bbbrem : beam-beam bremsstrahlung              ***" << '\n'
  << " *** radiative Bhabha scattering                    ***" << '\n'
  << " *** physics as published in 1994                   ***" << '\n'
  << " *** http://dx.doi.org/10.1016/0010-4655(94)90085-X ***" << '\n'
  << " *** authors R. Kleiss and H. Burkhardt             ***" << '\n'
  << " *** C++ version by H.B., 23/03/2017                ***" << '\n'
  << " ******************************************************" << '\n';
  
  // derived constants
  std::cout << "total energy                " << roots  << " GeV" << '\n';
  std::cout << "minimum photon energy       " << this->k0    << "  * beam energy" << '\n';
  if (this->cutoff>0)   std::cout << "optional cutoff in sqrt(-t) " << this->cutoff << " GeV" << '\n';
  if (this->UnWeighted) std::cout << "generated without weight" << '\n';
  
  if (k0<0 || k0>1)
  {
    std::cout << "*** Error*** wrong value for k0" << '\n';
    exit(1);
  }
  
  s = roots*roots;
  me2s = me2/s;
  logme2s = -log(me2s);
  z0 = k0/(1-k0);
  
  // approximate total cross section
  a1 = log((1+z0)/z0);
  double a2 = (log(1+z0))/z0;
  ac = a1/(a1+a2);
  sigapp = 8*pow(alpha,3)/me2*(-log(me2s))*(a1+a2)*tomb;
  std::cout << "approximate cross section " << sigapp << " millibarn" << '\n';
  
  // the initial-state momenta
  eb = roots*0.5;
  pb = sqrt(eb*eb-me2);
  rin2pb = 0.5/pb;
  
  p1[0] = eb;
  p1[1] = 0;
  p1[2] = 0;
  p1[3] = -pb;
  
  q1[0] = eb;
  q1[1] = 0;
  q1[2] = 0;
  q1[3] = pb;
  
  //-- extra setup if writing selected
  Make4Vec=write_evts_fl || write_evts_binary_fl;
  if (Make4Vec) setup_write_evts();
  
  if (verbose>1) std::cout << __FILE__ << " line " << __LINE__ << " " << __PRETTY_FUNCTION__ << " constructor end" << '\n';
};

void BBBrem::generate()
{
  // generate z
  double z;
  if (RanD(mt)<ac)
  {
    double temp0=RanD(mt);
    z = 1/(temp0*(exp(a1*RanD(mt))-1));
    if (verbose>2) std::cout << '\n' << "1 ac=" << ac << " temp0=RanD(mt)=" << temp0 << " z=" << z << '\n';
  }
  else
  {
    z = z0/RanD(mt);
    if (verbose>2) std::cout << '\n' << "2 ac=" << ac << " z = z0/RanD(mt)=" << z << '\n';
  }
  
  // bounds on t
  double y = me2s*z;
  double q0 = eb*y;
  double temp1 = pb*pb-eb*q0;
  double temp2 = temp1*temp1-me2*q0*q0;
  
  double sy=0,tmin=0,tmax=0;
  double ww2=0;
  // exit if temp2<0 (very very rare): put weight to 0
  // the `else' clause extends to the end of the routine
  if (temp2<0)
  {
    std::cout << " y too large: delta_t^2 =" << temp2 << '\n';
    weight = 0;
  }
  else
  {
    tmin = -2*(temp1+sqrt(temp2)); // in absolute big number, nearly -s
    tmax = me2*s*y*y/tmin; // in absolute very small number  roughly -z^2 m^6 / s^2
    
    // generate t
    sy = s*y;
    ww2 = sy+me2;
    double temp_1 = sy+tmax; // local to this part
    double rlamx = sqrt(temp_1*temp_1-4*ww2*tmax);
    if (temp_1<0) temp_1 = rlamx-temp_1;
    else  temp_1 = -4*ww2*tmax/(rlamx+temp_1);
    
    double b;
    do
    {
      b = exp(RanD(mt)*log(1+2*sy/temp_1));
      t = -b*z*z*me2/((b-1)*(b*z+b-1));
      if (t<tmin) std::cout << " t = " << t << " tmin=" << tmin << '\n';
    } while (t <= tmin);
  }
  
  // generate cgam
  double rlam = sqrt((sy-t)*(sy-t)-4*me2*t);
  double eps = 4*me2*ww2/(rlam*(rlam+ww2+me2-t));
  double rl = log((2+eps)/eps);
  double vgam = eps*(exp(RanD(mt)*rl)-1);
  double cgam = 1-vgam;
  double sgam = sqrt(vgam*(2-vgam));
  
  // generate azimuthal angles
  double const twopi = 2*M_PI; // mathematical constant
  double phi  = twopi*RanD(mt);
  double phig = twopi*RanD(mt);
  
  // construct momentum transfer q(mu) = q1(mu)-q2(mu)
  double ql = (2*eb*q0-t)*rin2pb;
  double qt = sqrt((tmax-t)*(t-tmin))*rin2pb;
  double q[4];
  q[0] = q0;
  q[1] = qt*sin(phi);
  q[2] = qt*cos(phi);
  q[3] = ql;
  
  if (Make4Vec) for(unsigned int i=0; i<4; ++i) q2[i]=q1[i]-q[i]; // construct momentum of outgoing positron in lab frame
  
  // find euler angles of p1(mu) in cm frame
  double r0 = eb+q[0];
  double w = sqrt(ww2);
  double rin2w = 0.5/w;
  double rinr0w = 1/(r0+w);
  double eta = -(sy+2*w*q[0]+t)*rin2w*rinr0w;
  
  double phot[4];
  phot[1] = -q[1]*(1+eta);
  phot[2] = -q[2]*(1+eta);
  phot[3] = pb*eta-ql*(1+eta);
  double phatl = rlam*rin2w; // long
  double phatt = hypot(phot[1],phot[2]); // trans
  double sfhat = phot[1]/phatt; // sin phi component
  double cfhat = phot[2]/phatt; // cos phi component
  double sthat = phatt/phatl;   // sin theta component
  double vthat;
  if (phot[3]>0) vthat = sthat*sthat/(1-sqrt(1-sthat*sthat));
  else          vthat = sthat*sthat/(1+sqrt(1-sthat*sthat));
  
  // rotate using these euler angles to get the qk direction in the cm
  double temp3 = sgam*sin(phig);
  double cfg = cos(phig);
  double temp4 = (vthat-1)*sgam*cfg+sthat*cgam;
  double qkhat[4];
  qkhat[0] = sy*rin2w;
  qkhat[1] = qkhat[0]*( cfhat*temp3+sfhat*temp4);
  qkhat[2] = qkhat[0]*(-sfhat*temp3+cfhat*temp4);
  qkhat[3] = qkhat[0]*(vthat+vgam-vthat*vgam-sthat*sgam*cfg-1);
  
  { // boost the photon momentum to the lab frame
    double temp5 = pb*qkhat[3];
    double temp6;
    if (temp5>0) temp6 = (me2*qkhat[0]*qkhat[0] +pb*pb*(qkhat[1]*qkhat[1]+qkhat[2]*qkhat[2])) /(eb*qkhat[0]+temp5);
    else temp6 = eb*qkhat[0]-temp5;
    qk[0] = (temp6+ FourVecProd(qkhat,q))/w;
  }
  if (Make4Vec)
  {
    double temp7 = (qk[0]+qkhat[0])*rinr0w;
    qk[1] = qkhat[1]+q[1]*temp7;
    qk[2] = qkhat[2]+q[2]*temp7;
    qk[3] = qkhat[3]+temp7*(-pb+ql);
    if (verbose>2)
    {
      std::cout << "y=" << y << " rlam=" << rlam << " t=" << t << " eps=" << eps
      << " vgam=" << vgam
      << " vthat=" << vthat
      << " sgam=" << sgam
      << " sfg=" << sin(phig)
      << " temp1=" << temp1
      << " temp1=" << temp1
      << " temp2=" << temp2
      << " temp3=" << temp3
      << " temp4=" << temp4
      //<< " temp5=" << temp5
      //<< " temp6=" << temp6
      << " temp7=" << temp7
      << " tmin=" << tmin
      << " tmax=" << tmax
      << '\n';
      std::cout << "qkhat=";
      for(unsigned int i=1; i<4; ++i) std::cout << std::setw(12) << qkhat[i] << ' '; std::cout << std::setw(12) << qkhat[0];
      std::cout << '\n';
    }
  }
  
  if (Make4Vec)
  {
    // construct p2 by momentum conservation
    for(unsigned int i=0; i<4; ++i) p2[i]=-q2[i]-qk[i];
    p2[0]+=2*eb;
  }
  
  // impose cut on the photon energy: qk[0]>eb*k0
  double c1,c2;
  if (qk[0]<eb*k0) weight = 0;
  else // the event is now accepted: compute matrix element and weight
  {
    c1 = log(1+z)/log((2+eps)/eps); // compute fudge factor c1
    { // compute fudge factor c2
      double temp8 = sy-tmax;
      double vnumx = sqrt(temp8*temp8-4*me2*tmax)+temp8;
      double temp9 = sy+tmax;
      double vdenx;
      if (temp9<0) vdenx = sqrt(temp9*temp9-4*ww2*tmax)-temp9;
      else vdenx = -4*ww2*tmax/(sqrt(temp9*temp9-4*ww2*tmax)+temp9);
      double temp10 = sy-tmin;
      double vnumn = sqrt(temp10*temp10-4*me2*tmin)+temp10;
      double temp11 = sy+tmin;
      double vdenn;
      if (temp11<0) vdenn = sqrt(temp11*temp11-4*ww2*tmin)-temp11;
      else vdenn = -4*ww2*tmin/(sqrt(temp11*temp11-4*ww2*tmin)+temp11);
      c2 = 2*logme2s/log((vnumx*vdenn)/(vdenx*vnumn));
      if (verbose>2)
      {
        std::cout << " sy=" << sy << " tmin=" << tmin << " tmax=" << tmax << " ww2=" << ww2
        << " vnumx=" << vnumx
        << " vnumn=" << vnumn
        << " vdenx=" << vdenx
        << " vdenn=" << vdenn
        << " temp8=" << temp8
        << " temp9=" << temp9
        << '\n';
      }
    }
    
    double zz;
    { // calculate zz
      // compute vector (small) r in cm frame, and (big) z
      double rlabl = (t-2*me2*y)*rin2pb;
      double rhat[4];
      rhat[0] = -(2*me2*y+(1-y)*t)*rin2w;
      double etar = rhat[0]*rinr0w;
      rhat[1] = -q[1]*(1+etar);
      rhat[2] = -q[2]*(1+etar);
      rhat[3] = rlabl+(pb-ql)*etar;
      if (verbose>2)
      {
        for(unsigned int i=1; i<4; ++i) std::cout << std::setw(12) << rhat[i] << ' '; std::cout << std::setw(12) << rhat[0];
        std::cout << '\n';
      }
      zz = s*( FourVecProd(rhat,qkhat) );
    }
    
    // the other invariants
    s1 = 4*eb*(eb-qk[0]);
    d1 = sy*rlam*(eps+vgam)*rin2w*rin2w;
    d2 = 0.5*sy;
    
    { // now that all invariantes are know, calculate the weight
      // the exact matrix element
      double rind1  = 1/d1;
      double rind12 = rind1*rind1;
      double rind2  = 1/d2;
      double rind22 = rind2*rind2;
      
      // Kleiss-Burkhardt approx cross section multiplied by t**2
      double temp20 = s +t-2*d2;
      double temp21 = s1+t+2*d1;
      double aa0 = (s*s+s1*s1+temp20*temp20+temp21*temp21)*rind1*rind2*(-t);
      double aa1 =  -4*me2*zz*zz*rind12*rind22;
      double aa2 =  -8*me2*(d1*d1+d2*d2)*rind1*rind2;
      double aa3 =  16*me2*me2*(d1-d2)*zz*rind12*rind22;
      double aa4 = -16*me2*me2*me2*(d1-d2)*(d1-d2)*rind12*rind22;
      double mex = aa0+aa1+aa2+aa3+aa4; // matrix element exact
      
      // the approximate matrix element without c1,2, multiplied by t**2
      double map = 4*s*s*rind1*rind2*(-t)*c1*c2;
      
      if (verbose>2)
      {
        std::cout
        << " temp20=" << temp20
        << " temp21=" << temp21
        << " zz=" << zz
        << " aa0=" << aa0
        << " aa1=" << aa1
        << " aa2=" << aa2
        << " aa3=" << aa3
        << " aa4=" << aa4 << '\n';
      }
      if (verbose>1)
      {
        if (verbose>2 || mex/map>10)
        std::cout
        << " c1=" << std::setw(10) << c1
        << " c2=" << std::setw(10) << c2
        << " mex=" << std::setw(12) << mex << " map=" << std::setw(12) << map << " mex/map=" << std::setw(10) << mex/map << '\n';
      }
      weight = mex/map*sigapp; // the weight
    }
    // the weight is now defined for both accepted and rejected events
    
  }
  
  // bookkeeping on the weight
  if (weight > wmax) wmax = weight;
  if (weight == 0) ++wnil;
  if (weight  < 0) ++wneg;
  ++w0;
  w1+=weight;
  w2+=weight*weight;
  sqrtmt=sqrt(-t);
  if (sqrtmt>cutoff) sumwcut+=weight;
  
  if (!UnWeighted)
  {
    if (write_evts_fl) write_evts();
    if (write_evts_binary_fl) write_evts_binary();
  }
  if (verbose>2) std::cout << __FILE__ << " line " << __LINE__ << " " << __PRETTY_FUNCTION__ << " end" << '\n';
}

void BBBrem::generate_UnWeighted()
{
  do
  {
    generate();
  } while (weight <= Unweighted_up_to*sigapp*RanD(mt)); // reject, to get weight 1 to Unweighted_up_to
  
  if (weight>Unweighted_up_to*sigapp) weight_remaining=weight/(Unweighted_up_to*sigapp); // calculate remaining weight, if MaxWeight not high enough
  else weight_remaining=1;
  
  if (write_evts_fl) write_evts();
  if (write_evts_binary_fl) write_evts_binary();
}

void BBBrem::finish() // analyze weight distribution, and print results
{
  double smc = w1/w0;
  double ser = sqrt(w2-w1*w1/w0)/w0;
  std::cout << "******** cross section evaluation ********" << '\n';
  std::cout << "last random number        " << std::setw(8)<< RanD(mt) << '\n';
  std::cout << "total # of events         " << std::setw(8)<< w0 << '\n';
  std::cout << "sum of weights            " << std::setw(8)<< w1 << '\n';
  std::cout << "sum of (weight**2)s       " << std::setw(8)<< w2 << '\n';
  if (UnWeighted) std::cout << "max weight occurred       " << std::setw(8)<< wmax/Unweighted_up_to << "    /sigapp=" << wmax/(Unweighted_up_to*sigapp) << '\n';
  else std::cout << "max weight occurred       " << std::setw(8)<< wmax << "    /sigapp=" << wmax/sigapp << '\n';
  std::cout << "# of zero weights         " << std::setw(8)<< wnil << "   fraction=" << wnil/w0 << '\n';
  std::cout << "# of negative weights     " << std::setw(8)<< wneg << '\n';
  std::cout << "computed  cross section   " << std::setw(8) << smc << "  millibarn" << '\n';
  if (cutoff>0)
  std::cout << "with cutoff               " << std::setw(8) << smc*sumwcut/w1 << "  millibarn" << '\n';
  std::cout << "                  +/-     " << std::setw(8) << ser << "  millibarn" << '\n';
  if (cutoff>0)
  std::cout << "fraction above cutoff     " << std::setw(8) << sumwcut/w1 << '\n';
  std::cout << "******************************************" << '\n';
  if (verbose && write_evts_fl)        std::cout << "text  output in  " << fname << '\n';
  if (verbose && write_evts_binary_fl) std::cout << "binary output in " << fname_bin << '\n';
}

void BBBrem::setup_write_evts(int prec)
{
  // set up optional output files
  fname    ="bbrems_ntuple.out";
  fname_bin="bbrems.MyBinNtuple";
  std::vector<std::string> varnames ={"q2_1","q2_2","q2_3","q2_4","p2_1","p2_2","p2_3","p2_4","qk_1","qk_2","qk_3","qk_4","weight","sqrtmt"};
  std::vector<std::string> varnames2={"t","d1","d2","s1"};
  if (WriteInvariants) for(auto varname2 : varnames2) varnames.push_back(varname2); // add to end
  if (write_evts_fl)
  {
    if (verbose>1) std::cout << "ascii output will be written to file " << fname << '\n';
    my_out.open(fname,std::ios::out);
    int wid=7+prec;
    for(auto varname : varnames) my_out << std::setw(wid) << varname;
    my_out << '\n';
    if (verbose>1) std::cout << __FILE__ << " line " << __LINE__ << " " << __PRETTY_FUNCTION__
    << " prec=" << prec << " number of digits written to text output file " << fname << '\n';
  }
  if (write_evts_binary_fl)
  {
    my_out_bin.open(fname_bin,std::ios::out | std::ios::trunc | std::ios::binary);
    my_out_bin  << "q2_1 q2_2 q2_3 q2_4 p2_1 p2_2 p2_3 p2_4 qk_1 qk_2 qk_3 qk_4 weight sqrtmt ";
    if (WriteInvariants) my_out_bin <<  "t d1 d2 s1 ";
    my_out_bin << std::endl;
    int16_t ShortOne=1; // BigEndian ShortOne = "0000 0001" LittleEndian "0001 0000"
    my_out_bin.write((char*) &ShortOne,sizeof(ShortOne) ); // write this directly as binary - to allow easy check in read if ByteSwap is needed
  }
}

void BBBrem::write_evts() // text output
{
  if (weight>0)
  {
    // write energy [0] component last, to be more similar to the fortran output
    my_out << " " << std::setprecision(6);
    for(unsigned int i=1; i<4; ++i) my_out << std::setw(12) << q2[i] << ' '; my_out << std::setw(12) << q2[0] << ' ';
    for(unsigned int i=1; i<4; ++i) my_out << std::setw(12) << p2[i] << ' '; my_out << std::setw(12) << p2[0] << ' ';
    for(unsigned int i=1; i<4; ++i) my_out << std::setw(12) << qk[i] << ' '; my_out << std::setw(12) << qk[0] << ' ';
    if (UnWeighted)  my_out << std::setw(12) << weight_remaining << ' '; else my_out << std::setw(12) << weight/sigapp << ' ';
    my_out << std::setw(12) << sqrtmt << ' ';
    if (WriteInvariants) my_out
    << std::setw(12) << t  << ' '
    << std::setw(12) << d1  << ' '
    << std::setw(12) << d2  << ' '
    << std::setw(12) << s1  << ' ';
    my_out << '\n';
  }
}

void BBBrem::write_evts_binary()
{
  const unsigned int WordLength=sizeof(double);
  if (weight>0)
  {
    double data[18];
    for(unsigned int i=0; i<3; ++i) data[i]  =q2[i+1]; data[ 3]=q2[0];
    for(unsigned int i=0; i<3; ++i) data[i+4]=p2[i+1]; data[ 7]=p2[0];
    for(unsigned int i=0; i<3; ++i) data[i+8]=qk[i+1]; data[11]=qk[0];
    if (UnWeighted) data[12]=weight_remaining; else data[12]=weight/sigapp;
    data[13]=sqrtmt;
    if (WriteInvariants)
    {
      data[14]=t;
      data[15]=d1;
      data[16]=d2;
      data[17]=s1;
      my_out_bin.write((char*) &data[0],18*WordLength);
    }
    else my_out_bin.write((char*) &data[0],14*WordLength);
  }
}

#include <sstream>
void BBBrem::read_evts_binary(const std::string fname) // example of how to read back the binary output
{
  std::ifstream BinFileIn(fname,std::ios::in | std::ios::binary);
  if (BinFileIn.is_open())
  { // ok, file is open
    if (verbose>1) std::cout << __FILE__ << " line " << __LINE__ << " " << __PRETTY_FUNCTION__  << " fname=" << fname << " open for read" << '\n';
  }
  else
  {
    std::cout << "read_evts_binary, *** Warning *** cannot open fname=" << fname << '\n';
    return;
  }
  std::string line;
  getline(BinFileIn,line);  // Stroustrup3  51, 598  http://www.cplusplus.com/reference/string/string/getline/
  if (verbose) std::cout << "header line:" << '\n';
  std::cout << line << '\n';
  std::istringstream istr(line);
  unsigned int nvar=0;
  while(istr)
  {
    std::string varnam;
    istr >> varnam;
    if (istr) ++nvar;
  }
  int16_t ShortOne=0;
  BinFileIn.read((char*) &ShortOne,sizeof(ShortOne) );
  if (ShortOne == 1)
  {
    if (verbose>1) std::cout << "read_evts_binary control word read OK. nvar=" << nvar << " now start to read binary data" << '\n';
  }
  else
  {
    std::cout << "read_evts_binary *** warning *** control read failed -- stop reading" << '\n';
    return;
  }
  const unsigned int Nvar=nvar,WordLength=8;
  double data[Nvar];
  size_t RecordLength=Nvar*WordLength;
  unsigned int nevts=0;
  while(BinFileIn)
  {
    BinFileIn.read((char*) data,RecordLength);
    ++nevts;
    if (nevts>10) break; // only read first events
    for(auto i=0;i<Nvar;++i) std::cout << std::setw(13) << data[i];
    std::cout << '\n';
  }
  if (verbose>1) std::cout << __FILE__ << " line " << __LINE__ << " " << __PRETTY_FUNCTION__  << " fname=" << fname << " done reading nevts=" << nevts << '\n';
}
