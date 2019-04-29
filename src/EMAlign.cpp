/*
 * Author: 1mingfei 
 * Date:   2019-04-17
 * Purpose: functions for EMHome class
 * align geometry
 * self-explained
 */

#include "gbCnf.h"
#include "userDefFnx.inl"

void EMHome::runAlign(gbCnf& cnfModifier, double halfThick)
{
  if (me == 0) 
    cout << "processing geometry alignment\n";
  c0 = move(cnfModifier.readLmpData(sparams["refFile"]));
  cnfModifier.getNBL(c0, dparams["Rcut"]);
  vector<Atom> c0GBAtom;
  double loc = 0.0;
  if (sparams["AlignFnx"] == "meanScore")
    loc = cnfModifier.getGBLoc(c0, c0GBAtom, meanScore);
  else if (sparams["AlignFnx"] == "stdScore")
    loc = cnfModifier.getGBLoc(c0, c0GBAtom, stdScore);

  c0 = cnfModifier.chopConfig(c0, loc, halfThick);
  if (me==0)
    cnfModifier.writeLmpDataDebug(c0, to_string(NI) + ".txt");
  for (int i = 0; i < NI; ++i)
  {
    if (i%nProcs != me) continue;
    c1 = move(cnfModifier.readLmpData("final." + to_string(i) + ".txt"));
    cnfModifier.getNBL(c1, dparams["Rcut"]);
    vector<Atom> c1GBAtom;
    double loc = 0.0;
    if (sparams["AlignFnx"] == "meanScore")
      loc = cnfModifier.getGBLoc(c1, c1GBAtom, meanScore);
    else if (sparams["AlignFnx"] == "stdScore")
      loc = cnfModifier.getGBLoc(c1, c1GBAtom, stdScore);

    c1 = cnfModifier.chopConfig(c1, loc, halfThick);
    double bestScore = cnfModifier.alignInPlane(c0, c1, c0GBAtom, c1GBAtom);
    cnfModifier.writeLmpDataDebug(c1, to_string(i) + ".txt");

    cout << "processed config " << i << " on processor " << me 
         << " out of total " << nProcs  << "processors, best score: " 
         << bestScore << endl;

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
      c1GB[j].pst[Z] += disp[Z];
      /*wrap back atoms*/
      wrapAtom(c1GB[j], c1.length);
    }
    sortAtomLexi(c1GB);

    /*calculate InPlane Misfit score*/
    double currScore = calInPlaneScore(c0GB, c1GB);

    if (currScore < bestScore)
    {
      bestIndex = i;
      bestScore = currScore;
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
    for (int j = 0; j < tmpAtoms.size(); ++j)
    {
      double dist = calDist(tmpLength, tmpAtoms[i], tmpAtoms[j]);
      if ((dist <= Rcut) && (j % cnf.atoms.size() - i != 0))
        res.push_back(j);
    }
    cnf.atoms[i].NBL = move(res);
    cnf.atoms[i].CN = cnf.atoms[i].NBL.size();
  }
}

int EMHome::gbCnf::getExpdParam(const Config& cnf, const double Rcut = 3.8)
{
  if (cnf.length[Z] > 2.0*Rcut)
    return 1;
  else
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

/*calculate distance between one atom in configuration and one from ref*/
double EMHome::gbCnf::calDist(const vector<double> length, const Atom& atm1,\
                              const Atom& atm2)
{
  double xi = atm1.pst[X];
  double xj = atm2.pst[X];
  double yi = atm1.pst[Y];
  double yj = atm2.pst[Y];
  double zi = atm1.pst[Z];
  double zj = atm2.pst[Z];
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
double EMHome::gbCnf::getGBLoc(Config& cnf, vector<Atom>& atm,
                       int (*f)(const vector<vector<double>>&, const vector<double>&))
{
  int nBins = 20;
  double bestLoc = 0.0;
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
  //get averaged location
  for (int i = 0; i < nBins; ++i)
    if (Count[i] > 0.01)
      Loc[i] /= Count[i];

  /*
  double bestScore = 0.0; //1210TB
  //double bestScore = std::numeric_limits<double>::max(); //FCC100
  for (int i = 0; i < nBins; ++i)
  {
    if (Count[i] > 0.01)
    {
      if (stddev(Score[i]) > bestScore) //1210TB
      //if (meanV(Score[i]) < bestScore) //FCC100
      {
        bestScore = stddev(Score[i]); //1210TB
        //bestScore = meanV(Score[i]); //FCC100
        bestLoc = Loc[i];
        atm = atomList[i];
      }
    }
  }
  */

  int index = (*f)(Score, Loc);
  bestLoc = Loc[index];
  atm = atomList[index];
  return bestLoc;
}

