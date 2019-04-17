#ifndef _EMHome_H_
#define _EMHome_H_

#include <math.h>
#include <mpi.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <numeric>
#include <queue>
#include <random>
#include <set>
#include <sstream>
#include <string>
#include <tuple> 
#include <unordered_set>
#include <unordered_map>
#include "armadillo"
#include "gbDef.h"
#include "gbCnf.h"
#include "gbUtl.h"

using std::pair; 
using std::cerr;
using std::cout;
using std::endl;
using std::ifstream;
using std::map;
using std::move;
using std::ofstream;
using std::set;
using std::setprecision;
using std::string;
using std::stringstream;
using std::to_string;
using std::unordered_map;
using std::vector;

class EMHome 
{
private:
  vector<Config> cnfs;
  
  class gbCnf;
  Config c0; //reference
  Config c1; //working on
  int NI;
  double engy; //individual energy
  double boltzEngy;
  double sumEngy;

public:
  int me, nProcs;

  unordered_map<string, double> dparams;
  unordered_map<string, int> iparams;
  unordered_map<string, string> sparams;

  EMHome(int argc, char* argv[]);
  ~EMHome();

  void parseArgs(int argc, char* argv[]);
  void initParam();
  void readParam();
  /*align all the configurations with a fixed thickness on each side of GB*/
  void runAlign(gbCnf&, double);
  /*gather energies from all the structures
   *and calculate probability */
  void getProb(gbCnf&, double, double);
  void estimateMean(gbCnf&);
};

#include "EMAlign.inl"
#endif
