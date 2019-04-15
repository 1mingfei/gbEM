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

  // input
  Config readLmpData(const string& fname);

  // output
  void writeLmpData(Config&, string);
  void writeLmpDataDebug(Config&, string);

  Config chopConfig(Config&, double, double);
  int getExpdParam(const Config&, const double);
  vector<Atom> expandCellZ(const Config&, const int);
  double calDist(vector<double>, const vector<Atom>&, int, int);
  void getNBL(Config&, double);
  double getGBLoc(Config&);

  void cnvVec2Mat(const vector<double>&, Config& c);
  void cnvMat2Vec(Config&);
  vector<double> cnvVecXY2VecAng(const vector<double>& v);
  void cnvprl2pst(Config&);

  void initBox(Config&);
  void initDomShift(Config&);
  void initDom(Config&);
  void wrapAtomPos(Config&);
  void wrapAtomPosLocal(Config&);
  void wrapAtomPos(vector<double>& pst, const vector<double>& box,
                   const vector<double>& len);
  vector<double> calculateRelativePst(
    Config& cnf, const vector<double>& pst);
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
