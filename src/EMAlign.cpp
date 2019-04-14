#include "gbCnf.h"

void EMHome::runAlign(gbCnf& cnfModifier, double halfThick)
{
  c0 = std::move(cnfModifier.readLmpData(sparams["refFile"]));
  MPI_Barrier(MPI_COMM_WORLD);
  for (int i = 0; i < NI; ++i)
  {
    if (i%nProcs != me) continue;
    c1 = std::move(cnfModifier.readLmpData("final." + to_string(i) + ".txt"));
    Config c2;
    c2 = cnfModifier.chopConfig(c1, 80.0, 10.0);
    cnfModifier.writeLmpData(c2, to_string(me) + ".txt");
  }

}

Config EMHome::gbCnf::chopConfig(Config& inCnf, double ctr, double hlf)
{
  Config outCnf;
  wrapAtomPos(inCnf);
  vector<Atom> tmpAtoms;
  for (auto& atm : inCnf.atoms)
  {
    if ((atm.pst[Y] >= (ctr - hlf)) && (atm.pst[Y] <= (ctr + hlf)))
    {
      tmpAtoms.push_back(atm);
    }
  }

  cout << tmpAtoms.size() << endl;
  outCnf.atoms = tmpAtoms;
  outCnf.natoms = tmpAtoms.size();
  outCnf.ntypes = inCnf.ntypes;
  outCnf.length = inCnf.length;
  outCnf.oldEngy = inCnf.oldEngy; //heir the energy
  outCnf.engy = inCnf.engy; //heir the energy

  return outCnf;
}
