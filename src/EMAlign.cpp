#include "gbCnf.h"
inline double stddev(std::vector<double> const & func)
{
    double mean = std::accumulate(func.begin(), func.end(), 0.0) / func.size();
    double sq_sum = std::inner_product(func.begin(), func.end(), func.begin(), 0.0,
        [](double const & x, double const & y) { return x + y; },
        [mean](double const & x, double const & y) { return (x - mean)*(y - mean); });
    return sq_sum / ( func.size() - 1 );
}


void EMHome::runAlign(gbCnf& cnfModifier, double halfThick)
{
  c0 = std::move(cnfModifier.readLmpData(sparams["refFile"]));
  MPI_Barrier(MPI_COMM_WORLD);
  for (int i = 0; i < NI; ++i)
  {
    if (i%nProcs != me) continue;
    c1 = std::move(cnfModifier.readLmpData("final." + to_string(i) + ".txt"));
    Config c2;
    cnfModifier.getNBL(c1, 3.8);
    double loc = cnfModifier.getGBLoc(c1);
    c2 = cnfModifier.chopConfig(c1, loc, 10.0);
    cnfModifier.writeLmpData(c2, to_string(i) + ".txt");
  }
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
  std::cout << factor << std::endl;
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

double EMHome::gbCnf::getGBLoc(Config& cnf)
{
  int nBins = 20;
  vector<vector<double>> Score(nBins); //sum CN
  vector<double> Count(nBins, 0); //count how many atoms in each bin
  vector<double> Loc(nBins, 0); //in Y dimension
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
      }
    }
    //std::cout << i << " " << stddev(Score[i]) << " " << Loc[i] << std::endl;
  }
  return bestLoc;
}
