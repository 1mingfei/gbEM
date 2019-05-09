/*
 * Author: 1mingfei 
 * Date:   2019-04-14
 * Purpose: inline functions for EMHome
 * self-explained
 */

#include "gbCnf.h"
#include "EMHome.h"

EMHome::EMHome(int argc, char* argv[]) {
  MPI_Comm_size(MPI_COMM_WORLD, &nProcs);
  MPI_Comm_rank(MPI_COMM_WORLD, &me);

  iparams["Nconfigs"] = 4;
  dparams["halfThick"] = 10.0;
  dparams["T"] = 300.0;
  sparams["refFile"] = "lmp.init";
  sparams["AlignFnx"] = "meanScore";
  dparams["Rcut"] = 3.8;
  sparams["elem0"] = "Mg";
  sparams["elem1"] = "Zn";
  parseArgs(argc, argv);
  initParam();
  NI = iparams["Nconfigs"];
  double T = dparams["T"];
  sumEngy = 0.0;
  double halfThick = dparams["halfThick"];
  double N = double(iparams["N"]);
  gbCnf cnfModifier(*this);

  runAlign(cnfModifier, halfThick);

  cnfs.assign(NI, c0);
  getProb(cnfModifier, T, N);

  estimateMean(cnfModifier);
}

EMHome::~EMHome() {}
