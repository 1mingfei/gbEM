#ifndef _GB_CNF_H
#define _GB_CNF_H

#include "armadillo"
#include "EMHome.h"
using arma::mat;
using arma::vec;
using std::vector;

class EMHome::gbCnf 
{
  EMHome& hm;
  vector<Config>& cnfs;
  unordered_map<string, string>& sparams;

public:
  double rcut;
  int NI;  // number of individuals per generations
  gbCnf(EMHome& x, double rc = 6.0)
      : hm(x), cnfs(x.cnfs), sparams(x.sparams), rcut(rc), NI(x.NI) {}

  //gbInCnf.cpp
  Config readLmpData(const string& fname);

  // output
  //gbOutCnf.cpp
  void writeLmpData(Config&, string);
  void writeLmpDataDebug(Config&, string);

  //alignment
  //EMAlign.cpp
  Config chopConfig(Config&, double, double);
  double alignInPlane(const Config&, Config&, vector<Atom>&, vector<Atom>&);
  int getExpdParam(const Config&, const double);
  vector<Atom> expandCellZ(const Config&, const int);
  double calDist(vector<double>, const vector<Atom>&, int, int);
  void getNBL(Config&, double);
  /*return GB Y location value, and atm stores atoms in GB level bin*/
  double getGBLoc(Config&, vector<Atom>&);
  
  //gbInCnf.cpp
  void cnvVec2Mat(const vector<double>&, Config& c);
  void cnvMat2Vec(Config&);
  vector<double> cnvVecXY2VecAng(const vector<double>& v);

  //gbBox.cpp
  void initBox(Config&);
  void wrapAtomPos(Config&);
  void cnvprl2pst(Config&);
};


/**************************************************
 * check if atoms in the range
 **************************************************/
inline bool isIntheRange(const vector<double>& pst, const vector<double>& ctr,
                         const vector<double>& rng) {
  return std::abs(pst[X] - ctr[X]) <= rng[X] &&
         std::abs(pst[Y] - ctr[Y]) <= rng[Y] &&
         std::abs(pst[Z] - ctr[Z]) <= rng[Z];
}

/**************************************************
 * check if relative position is in the box
 **************************************************/
inline bool isInside(const vector<double>& pos) {
  return pos[0] <= 1 && pos[0] >= 0 && pos[1] >= 0 && pos[1] <= 1 &&
         pos[2] >= 0 && pos[2] <= 1;
}
#endif
