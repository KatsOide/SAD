// RunFlags.C   command line utility interface, H. Burkhardt

#include <iostream>
#include <iomanip>
#include <sstream>

#include "RunFlags.h"

RunFlags::RunFlags(int argc,char* const argv[],const char* cflags,const std::vector<std::string>& ArgvDefault)  // constructor
: verbose(0),b(false),h(false),t(false),u(false),v(false),V(false),errflg(false),ArgvDefault(ArgvDefault) // initialize
{
  argv0=argv[0];
  int ic;
  while( (ic=getopt(argc,argv,cflags)) != EOF ) // decode options
  {
    switch(ic)
    {
      case 'b':
        b=true;
        break;
      case 'h':
        h=true;
        break;
      case 't':
        t=true;
        break;
      case 'u':
        u=true;
        break;
      case 'v':
        verbose=1;
        v=true;
        break;
      case 'V':
        verbose=2;
        V=true;
        break;
      case '?':
        errflg=true;
    }
  }
  unsigned int nArg=std::max((unsigned int)argc-optind+1,(unsigned int)ArgvDefault.size()+1);
  ArgvVec.resize(nArg);
  ArgvVec[0]=argv0; // program name, always there
  for(auto i=0; i<ArgvDefault.size(); ++i) ArgvVec[i+1]=ArgvDefault[i];
  for(auto i=1; (i+optind-1)<argc; ++i) ArgvVec[i]=argv[i+optind-1];
  if (errflg)
  {
    std::cerr << "error with parameters" << '\n';
    ShowUsage();
    exit(2);
  }

  if (verbose>1)
  {
    std::cout << __PRETTY_FUNCTION__ << " in " << __FILE__ << " line " << std::setw(4) << __LINE__ << " " << '\n';
    for(int i=0; i<argc; ++i) std::cout << "arg" << i << "=" << argv[i] << '\n';
    std::cout << "argc=" << argc << " optind=" << optind << " ArgvVec.size()=" << ArgvVec.size() << " ArgvDefault.size()=" << ArgvDefault.size() << '\n';
    for(unsigned int i=0;i<ArgvVec.size();++i) std::cout << "ArgvVec[" << i << "]=" << ArgvVec[i] << '\n';
    std::cout << '\n';
    std::cout << Show() << '\n';
  }
  if (h) // show help and exit
  {
    ShowUsage();
    exit(0);
  }
}

std::string RunFlags::Show() const
{
  std::ostringstream OutStr;
  OutStr << "option settings are " << '\n'
  << "  b (write binary ntuple)=" << b << '\n'
  << "  h (help)=" << h << '\n'
  << "  t (write formatted text ntuple)=" << t << '\n'
  << "  u (generate (mostly) unweighted)" << u << '\n'
  << "  v (verbose)=" << v << '\n'
  << "  V (more verbose)=" << V << '\n'
  ;
  return OutStr.str();
}

void RunFlags::ShowUsage() const
{
  std::cout << "Help for " << argv0 << " to run BBBrem as unix command line tool" << '\n'
  << BOLD << "SYNOPSIS" << NORMAL << '\n'
  << "  " << argv0 << " [OPTION] [Eb] [rk] [nevent] [seed] [cutoff]" << '\n'
  << "   where                                                                             default " << '\n'
  << "      Eb is the beam energy in GeV                                                   " << ArgvVec[1] << '\n'
  << "      rk the relative photon cutoff  - a number between 0 and 1, typically 0.01      " << ArgvVec[2] << '\n'
  << "      nevent, # events to be generated                                               " << ArgvVec[3] << '\n'
  << "      seed used for the Mersenne Twister random generator                            " << ArgvVec[4] << '\n'
  << "      cutoff, generate only sqrt(-t) > cutoff [GeV]                                  " << ArgvVec[5] << '\n'
  << BOLD << "OPTIONS" << NORMAL << '\n'
  << "    -b   write binary ntuple" << '\n'
  << "    -h   just show this help" << '\n'
  << "    -t   write text ntuple" << '\n'
  << "    -u   generate (mostly) unweighted" << '\n'
  << "    -v   verbose (turn debug on)" << '\n'
  << "    -V   more verbose" << '\n'
  << BOLD << "EXAMPLEs" << NORMAL << '\n'
  << argv0 << " -v                            ; echo 'verbose run with default values'" << '\n'
  << argv0 << " -vt 45.6 0.01 1000 0 6.e-11   ; echo 'short run with cutoff, writing 1000 events to text file'" << std::endl
  ;
}
