/*
 * Author: 1mingfei 
 * Date:   2019-04-14
 * Purpose: inline functions for EMAlign
 * self-explained
 * To-do:  lexicographic sort may need a grace tolerance
 */

/*calculate mean of a vector with predefined prob*/
inline double getProbMean(vector<double> const & A,
                          vector<double> const & prob) {
  assert(A.size() == prob.size());
  double sum = 0.0;
  for (int i = 0; i < A.size(); ++i)
    sum += A[i] * prob[i];
  return sum;
}

/*calculate std deviation of a vector with predefined prob*/
inline double getStdDevProb(vector<double> const & A,
                            vector<double> const & prob,
                            double const mean) {
  assert(A.size() == prob.size());
  double sum = 0.0;
  for (int i = 0; i < A.size(); ++i)
    sum += ((A[i] - mean) * (A[i] - mean) * prob[i]);
  double std = sqrt(sum); 
  return std;
}

/*calculate mean val of a vector*/
inline double meanV(vector<double> const & A) {
  return std::accumulate(A.begin(), A.end(), 0.0) / A.size();
}

/*calculate std deviation of a vector*/
inline double stddev(vector<double> const & A) {
  double mean = meanV(A);
  double sq_sum = std::inner_product(A.begin(), A.end(), A.begin(), 
    0.0, [](double const & x, double const & y) {return x + y;},
    [mean](double const & x, double const & y) {return (x - mean)*(y - mean);});
  return sq_sum / (A.size() - 1);
}

/*Sort points lexicographically*/
inline void sortAtomLexi(vector<Atom>& atmList) {
  sort(atmList.begin(), atmList.end(),
    [](const Atom& a, const Atom& b) -> bool
    {return a.pst[X] < b.pst[X] || (a.pst[X] == b.pst[X] && a.pst[Z] < b.pst[Z]);});
}

inline vector<double> calDisp(Atom& a1, Atom& a2) {
  vector<double> res(3);
  res[X] = a1.pst[X] - a2.pst[X];
  res[Y] = a1.pst[Y] - a2.pst[Y];
  res[Z] = a1.pst[Z] - a2.pst[Z];
  return res;
}

inline void wrapAtom(Atom& atm, vector<double> length) {
  /*wrap back atoms in x and z direction*/
  int tmp = std::floor(atm.pst[X]/length[X]);
  if (tmp) {
    atm.pst[X] -= tmp*length[X];
  }
  tmp = std::floor(atm.pst[Z]/length[Z]);
  if (tmp) {
    atm.pst[Z] -= tmp*length[Z];
  }
}

/*calculate in-plane misfit score function
 *Note1: this only calculate in plane atoms*/
inline double calInPlaneScore(vector<Atom> list0, vector<Atom> list1) {
  int minSize = std::min(list0.size(), list1.size());
  double sum = 0.0;
  for (int i = 0; i < minSize; ++i) {
    vector<double> disp = calDisp(list0[i], list1[i]);
    double factor;
    if (list0[i].tp == list1[i].tp) {
      factor = 1.0;
    } else {
      factor = 10.0;
    }
    sum += factor*std::sqrt(disp[X]*disp[X] + disp[Z]*disp[Z]);
  }
  sum /= double(minSize);
  return sum;
}
