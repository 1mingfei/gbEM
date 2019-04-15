#include "gbCnf.h"

/*calculate std deviation of a vector*/
inline double stddev(std::vector<double> const & A)
{
  double mean = std::accumulate(A.begin(), A.end(), 0.0) / A.size();
  double sq_sum = std::inner_product(A.begin(), A.end(), A.begin(), 
    0.0, [](double const & x, double const & y) {return x + y;},
    [mean](double const & x, double const & y) {return (x - mean)*(y - mean);});
  return sq_sum / ( A.size() - 1 );
}

/*Sort points lexicographically*/
inline void sortAtomLexi(vector<Atom>& atmList)
{
  sort(atmList.begin(), atmList.end(),
    [](const Atom& a, const Atom& b) -> bool
    {return a.pst[X] < b.pst[X] || (a.pst[X] == b.pst[X] && a.pst[Z] < b.pst[Z]);});
}

inline vector<double> calDisp(Atom& a1, Atom& a2)
{
  vector<double> res(3);
  res[X] = a1.pst[X] - a2.pst[X];
  res[Y] = a1.pst[Y] - a2.pst[Y];
  res[Z] = a1.pst[Z] - a2.pst[Z];
  return res;
}

inline void wrapAtom(Atom& atm, vector<double> length)
{
  /*wrap back atoms*/
  int tmp = std::floor(atm.pst[X]/length[X]);
  if (tmp)
    atm.pst[X] -= tmp*length[X];
  tmp = std::floor(atm.pst[Z]/length[Z]);
  if (tmp)
    atm.pst[Z] -= tmp*length[Z];
}

/*calculate in-plane misfit score function
 *Note1: this only calculate in plane atoms*/
inline double calInPlaneScore(std::vector<Atom> list0, std::vector<Atom> list1)
{
  int minSize = std::min(list0.size(), list1.size());
  double sum = 0.0;
  for (int i = 0; i < minSize; ++i)
  {
    std::vector<double> disp = calDisp(list0[i], list1[i]);
    double factor;
    if ((list0[i].tp == 4) && (list1[i].tp == 4))
      factor = 1.0;
    else if ((list0[i].tp == 1) && (list1[i].tp == 1))
      factor = 1.0;
    else if ((list0[i].tp == 1) && (list1[i].tp == 4))
      factor = 10.0;
    else if ((list0[i].tp == 4) && (list1[i].tp == 1))
      factor = 10.0;
    else
      factor = 1.0;
    sum += factor*(disp[X]*disp[X] + disp[Z]*disp[Z]);
  }
  sum /= double(minSize);
  return(sum);
}

void EMHome::runAlign(gbCnf& cnfModifier, double halfThick)
{
  c0 = std::move(cnfModifier.readLmpData(sparams["refFile"]));
  cnfModifier.getNBL(c0, 3.8);
  vector<Atom> c0GBAtom;
  double loc = cnfModifier.getGBLoc(c0, c0GBAtom);
  //c0 = cnfModifier.chopConfig(c0, loc, halfThick);
  c0 = cnfModifier.chopConfig(c0, loc, 1.0);
  cnfModifier.writeLmpDataDebug(c0, "20.txt");
  MPI_Barrier(MPI_COMM_WORLD);
  for (int i = 0; i < NI; ++i)
  {
    if (i%nProcs != me) continue;
    c1 = std::move(cnfModifier.readLmpData("final." + to_string(i) + ".txt"));
    cnfModifier.getNBL(c1, 3.8);
    vector<Atom> c1GBAtom;
    double loc = cnfModifier.getGBLoc(c1, c1GBAtom);
    //c1 = cnfModifier.chopConfig(c1, loc, halfThick);
    c1 = cnfModifier.chopConfig(c1, loc, 1.0);
    double bestScore = cnfModifier.alignInPlane(c0, c1, c0GBAtom, c1GBAtom);
    std::cout << i << " best score: " << bestScore << std::endl;
    cnfModifier.writeLmpDataDebug(c1, to_string(i) + ".txt");
  }
}

/*shift c1 x and z coord in order to be aligned with c0
 * return best score */
double EMHome::gbCnf::alignInPlane(const Config& c0, Config& c1,\
                                 vector<Atom>& c0GB, vector<Atom>& c1GB)
{
  sortAtomLexi(c0GB);
  sortAtomLexi(c1GB);
  vector<double> disp(3);
  double bestScore = std::numeric_limits<double>::max();
  int bestIndex = c0GB.size();
  vector<Atom> c1Copy = c1GB;
  for (int i = 0; i < c0GB.size(); ++i)
  {
    c1GB = c1Copy;
    disp = calDisp(c0GB[i], c1GB[0]);
    for (int j = 0; j < c1GB.size(); ++j)
    {
      /*apply displacement to c1*/
      c1GB[j].pst[X] += disp[X];
      c1GB[j].pst[Y] += 0.0; 
      c1GB[j].pst[Z] += disp[Z];

      /*wrap back atoms*/
      wrapAtom(c1GB[j], c1.length);

      /*calculate InPlane Misfit score*/
      double currScore = calInPlaneScore(c0GB, c1GB);
      //double currScore = calInPlaneScore(c0.atoms, c1.atoms);

      if (currScore < bestScore)
      {
        bestIndex = i;
        bestScore = currScore;
      }
    }
  }

  /*update all c1*/
  disp = calDisp(c0GB[bestIndex], c1Copy[0]);
  for (Atom& atm : c1.atoms)
  {
    atm.pst[X] += disp[X];
    atm.pst[Z] += disp[Z];
    wrapAtom(atm, c1.length);
  }
  return bestScore;
}

/*chop hlf distance away from center and move ctr to cell ctr*/
Config EMHome::gbCnf::chopConfig(Config& inCnf, double ctr, double hlf)
{
  Config outCnf;
  wrapAtomPos(inCnf);
  vector<Atom> tmpAtoms;
  for (auto& atm : inCnf.atoms)
  {
    if ((atm.pst[Y] >= (ctr - hlf)) && (atm.pst[Y] <= (ctr + hlf)))
    {
      atm.pst[Y] -= (ctr - inCnf.length[Y]/2.0);
      tmpAtoms.push_back(atm);
    }
  }

  outCnf.atoms = tmpAtoms;
  outCnf.natoms = tmpAtoms.size();
  outCnf.ntypes = inCnf.ntypes;
  outCnf.cell = inCnf.cell;
  outCnf.length = inCnf.length;
  outCnf.oldEngy = inCnf.oldEngy; //heir the energy
  outCnf.engy = inCnf.engy; //heir the energy

  return outCnf;
}

void EMHome::gbCnf::getNBL(Config& cnf, double Rcut = 3.8)
{
  int factor = getExpdParam(cnf, Rcut);
  vector<double> tmpLength;
  tmpLength = cnf.length;
  tmpLength[Z] *= double(factor);
  vector<Atom> tmpAtoms = expandCellZ(cnf, factor);
  for (int i = 0; i < cnf.atoms.size(); ++i)
  {
    vector<int> res;
    //for (int j = 0; j < cnf.atoms.size(); ++j)
    for (int j = 0; j < tmpAtoms.size(); ++j)
    {
      //double dist = calDist(cnf.length, cnf.atoms, i, j);
      double dist = calDist(tmpLength, tmpAtoms, i, j);
      if ((dist <= Rcut) && (j % cnf.atoms.size() - i != 0))
        res.push_back(j);
    }
    cnf.atoms[i].NBL = std::move(res);
    cnf.atoms[i].CN = cnf.atoms[i].NBL.size();
  }
}

int EMHome::gbCnf::getExpdParam(const Config& cnf, const double Rcut = 3.8)
{
  if (cnf.length[Z] > 2.0*Rcut)
    return 1;
  else
    //return 2;
    return(static_cast<int>((2.0*Rcut/cnf.length[Z])+1));
}

/* expand cell in +/- Z direction */
vector<Atom> EMHome::gbCnf::expandCellZ(const Config& cnf, const int factor)
{
  vector<Atom> res;
  int initSize = cnf.atoms.size();
  for (int i = 0; i < initSize; ++i)
  {
    res.push_back(cnf.atoms[i]);
  }
  for (int i = initSize; i < initSize*factor; ++i)
  {
    Atom atm = cnf.atoms[i%initSize];
    res.push_back(Atom(i, atm.tp, atm.pst[X], atm.pst[Y],
           atm.pst[Z] + cnf.length[Z]*int(i/initSize)));
  }
  return res;
}

double EMHome::gbCnf::calDist(vector<double> length, const vector<Atom>& atoms,\
    int i, int j)
{
  double xi = atoms[i].pst[X];
  double xj = atoms[j].pst[X];
  double yi = atoms[i].pst[Y];
  double yj = atoms[j].pst[Y];
  double zi = atoms[i].pst[Z];
  double zj = atoms[j].pst[Z];
  double a, b, c;

  if (xj - xi >= 0.5 * length[X])
    a = (xi - xj + length[X]);
  else if (xj - xi <  -0.5 * length[X])
    a = (xi - xj - length[X]); 
  else
    a = xi - xj;

  b = yi - yj;

  if (zj - zi >= 0.5 * length[Z])
    c = (zi - zj + length[Z]);
  else if (zj - zi <  -0.5 * length[Z])
    c = (zi - zj - length[Z]); 
  else
    c = zi - zj;

  double dist = sqrt(a*a + b*b + c*c);
  return dist;
}

/*return GB Y location value, and atm stores atoms in GB level bin*/
double EMHome::gbCnf::getGBLoc(Config& cnf, vector<Atom>& atm)
{
  int nBins = 20;
  vector<vector<double>> Score(nBins); //sum CN
  vector<double> Count(nBins, 0); //count how many atoms in each bin
  vector<double> Loc(nBins, 0); //in Y dimension
  vector<vector<Atom>> atomList(nBins);
  for (Atom& atm : cnf.atoms)
  {
    for (int i = 0; i < nBins; ++i)
    {
      if (((65.0 + static_cast<double>(i)) < atm.pst[Y]) &&\
          (atm.pst[Y] <= (65.0 + static_cast<double>(i+1))))
      {
        Score[i].push_back(atm.CN);
        Count[i] += 1.0;
        Loc[i] += atm.pst[Y];
        atomList[i].push_back(atm);
      }
    }
  }
  double bestScore = 0.0, bestLoc = 0.0;
  for (int i = 0; i < nBins; ++i)
  {
    if (Count[i] > 0.01)
    {
      //Score[i] /= Count[i];
      Loc[i] /= Count[i];
      //if (std::abs(Score[i] - 12.0) > bestScore)
      if (stddev(Score[i]) > bestScore)
      {
        bestScore = stddev(Score[i]);
        bestLoc = Loc[i];
        atm = atomList[i];
      }
    }
    //std::cout << i << " " << stddev(Score[i]) << " " << Loc[i] << std::endl;
  }
  return bestLoc;
}
