#ifndef pfDefines
#define pfDefines

#include <functional>
#include <numeric>
#include <utility>
#include <vector>
#include "Atom.h"

using std::cout;
using std::endl;
using std::set;
using std::string;
using std::vector;
using std::pair;

#define EXT 0
#define DIM 3
#define BOPORD 8
#define MAXLEN 1024
#define D2A 57.2957795131  // (180. / PI)
#define PI 3.14159265359
#define INVPI 0.31830988618
#define PFROOT 0
#define EVA3_GPA 160.21766208
#define SQRT3 1.73205080757
#define MAXATOMS 2000

// nearest neighbors
// 1st 8 ; 2nd 6; 3rd 12; 4th 24; 5th 8
enum { N1 = 8, N2 = 14, N3 = 26, N4 = 50, N5 = 58 };
enum { XX = 0, YY = 1, ZZ = 2, XY = 3, YZ = 4, ZX = 5 };
enum { X = 0, Y = 1, Z = 2 };

class Config {
public:
  int natoms, ntypes;
  double engy, aveord;
  double oldEngy;
  double hEngy; // energy defined by convex hull

  vector<double> cell;    // lox, loy, loz, hix, hiy, hiz, xy xz yz
  vector<double> center;  // center postion
  vector<double> length;  // length of three edges
  vector<double> bvx, tvx, bvy, tvy, bvz, tvz;
  vector<Atom> atoms;
  vector<int> bondbn;
  vector<int> qlbn;
  pair<int, int> atomNum;
  double conc;
  double QE; //QuasiEntropy

  bool operator<(const Config &b) const { return this->engy < b.engy; }
  Config()
      : natoms(0),
        ntypes(0),
        engy(0.0),
        oldEngy(0.0),
        hEngy(0.0),
        QE(0.0),
        conc(0.0),
        aveord(0.0),
        cell(9),
        center(3),
        length(3),
        bvx(3),
        tvx(3),
        bvy(3),
        tvy(3),
        bvz(3),
        tvz(3){};
  Config(int n)
      : natoms(0),
        ntypes(0),
        engy(0.0),
        oldEngy(0.0),
        hEngy(0.0),
        QE(0.0),
        conc(0.0),
        aveord(0.0),
        cell(9),
        center(3),
        length(3),
        bvx(3),
        tvx(3),
        bvy(3),
        tvy(3),
        bvz(3),
        tvz(3){};
  ~Config(){};
friend class EMHome;
};

#endif
