/*
 * Author: 1mingfei 
 * Date:   2019-04-17
 * Purpose: functions for EMHome class 
 * estimate mean color
 * To-do: use a k-d Tree to better searching 1NN
 */

#include "gbCnf.h"

/*gather energies from all the structures
 *and calculate probability */
void EMHome::getProb(gbCnf& cnfModifier, double T, double N)
{
  MPI_Barrier(MPI_COMM_WORLD);
  if (me == 0)
    std::cout << "processing Boltzmann probability\n";
  
  int quotient = NI / nProcs;
  int remainder = NI % nProcs;
  int nCycle = remainder ? (quotient + 1) : quotient;
  /*calculate total boltzmann energy Z*/
  for (int j = 0; j < nCycle; ++j)
  {
    double* subEngys = new double[nProcs];
    for (int i = (j * nProcs); i < ((j + 1) * nProcs); ++i)
    {
      if (i % nProcs != me) continue;
      boltzEngy = 0.0; //important to initialize here
      if (i >= NI) continue;

      c1 = std::move(cnfModifier.readLmpData(to_string(i) + ".txt"));
      engy = c1.engy / N;
      boltzEngy = exp(-engy/KB/T);
      cout << "read in config " << i << " on processor " << me 
           << " out of total " << nProcs  << " processors, energy: "
           << engy << " BE " << boltzEngy << endl;
      /*cnf:  oldEngy -- normalized energy
       *      engy    -- Boltzmann energy */
      c1.oldEngy = engy;
      c1.engy = boltzEngy;
      cnfs[i] = move(c1);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Gather(&boltzEngy, 1, MPI_DOUBLE, subEngys, 1, MPI_DOUBLE, 0,
               MPI_COMM_WORLD);

    if (me == 0)
    {
      cout << "subset energies: ";
      for (int i = 0; i < nProcs; ++i)
        cout << subEngys[i] << "   ";
      cout << endl;
      for (int i = 0; i < nProcs; ++i)
        sumEngy += subEngys[i];
      cout << "running total: " << sumEngy << endl;
    }
    delete [] subEngys;
  }
  MPI_Bcast(&sumEngy, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  cout << "rank " << me << " sum energy: " << sumEngy << endl;
  /*calculate total boltzmann energy Z*/
  for (int i = 0; i < NI; ++i)
  {
    if (i % nProcs != me) continue;
    cnfs[i].oldEngy = cnfs[i].engy/sumEngy;
    cnfModifier.writeLmpDataDebug(cnfs[i], to_string(i) + ".txt");
  }
  /*Note: by now all the structures are stored in cnfs, with oldEngy represent
   *probability. */
}

/*calculate which is the nearest point in reference system and average pos 
 *then calculate type mean and std. */
void EMHome::estimateMean(gbCnf& cnfModifier)
{
  int nAtom = c0.atoms.size();
  for (int i = 0; i < nAtom; ++i)
  {
    if (i % nProcs != me) continue;

    vector<double> dists;
    vector<double> posX;
    vector<double> posY;
    vector<double> posZ;
    vector<double> type;
    vector<double> prob;

    for (int j = 0; j < NI; ++j)
    {
      /*get which atom in cnfs[j] is closest to atom[i] in c0. */
      double dist;
      int closest = getNNID(c0, i, cnfs[j], dist);
      Atom tmpAtom = cnfs[j].atoms[closest];
      
      dists.push_back(dist);
      posX.push_back(tmpAtom.pst[X]);
      posY.push_back(tmpAtom.pst[Y]);
      posZ.push_back(tmpAtom.pst[Z]);
      double tmpType = (tmpAtom.tp == 1) ? 1.0 : 2.0;
      type.push_back(tmpType);
      prob = cnfs[j].oldEngy;
    }

    /*get some stats*/
    double meanDist = getProbMean(dists, prob);
    double meanX = getProbMean(posX, prob);
    double meanY = getProbMean(posY, prob);
    double meanZ = getProbMean(posZ, prob);
    double meanType = getProbMean(type, prob);
    c0.atoms[i].tpMean = meanType;

    c0.atoms[i].posStd = getStdDevProb(dists, prob, meanDist);
    c0.atoms[i].tpStd = getStdDevProb(type, prob, meanType);

    cout << "on rank " << me << " atomID " c0.atoms[i].id 
         << " meanX: " << meanX << " meanY: " << meanY << " meanZ: " << meanZ
         << " posStd: " << c0.atoms[i].posStd 
         << " typeMean: " << c0.atoms[i].tpMean 
         << " typeStd: " << c0.atoms[i].tpStd << endl;
  }

}


