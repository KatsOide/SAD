// RunFlags.h   command line utility interface, H. Burkhardt

#ifndef RunFlags_h
#define RunFlags_h 1

#include <vector>
#include <getopt.h>

#ifndef EscapeSequences_h
static const char* const NORMAL="\033[0m";
static const char* const BOLD="\033[1m";
#endif
#define EscapeSequences_h 1

class RunFlags
{
public:
  unsigned int verbose;
  bool b,h,t,u,v,V,errflg;
  std::vector<std::string> ArgvVec; // arguments, without option flags
  RunFlags(){;} // defult constructor
  RunFlags(int argc,char* const argv[],const char* cflags,const std::vector<std::string>& ArgvDefault=std::vector<std::string>());
  void ShowUsage() const;
  std::string Show() const;
private: // used in ShowUsage  -  to be provided with each main calling RunFlags
  const char* argv0;  // Name of calling program
  std::vector<std::string> ArgvDefault;
};

#endif
