#ifndef _ATOM_
#define _ATOM_

using std::vector;
//atom with common neigh analysis results
class Atom {
private:
public:
  int id, tp, CN, CNA;
  //int posStd, tpMean, tpStd;
  double pst[3], prl[3];
  vector<int> NBL;

  Atom() : id(0), tp(1) {};
  Atom(int n) : id(n), tp(1) {};
  Atom(int _id, double x, double y, double z) : id(_id), tp(1) {
    pst[0] = x, pst[1] = y, pst[2] = z;
  };
  Atom(int _id, int _tp, double x, double y, double z):
    id(_id), tp(_tp) {
    pst[0] = x, pst[1] = y, pst[2] = z;
  }
  ~Atom(){};

  bool operator<(const Atom &b) const { return id < b.id; }
};

#endif
