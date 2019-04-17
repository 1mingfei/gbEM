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

      c0 = std::move(cnfModifier.readLmpData(to_string(i) + ".txt"));
      engy = c0.engy / N;
      boltzEngy = exp(-engy/KB/T);
      cout << "read in config " << i << " on processor " << me 
           << " out of total " << nProcs  << " processors, energy: "
           << engy << " BE " << boltzEngy << endl;
      /* cnf:  oldEngy -- normalized energy
       *       engy    -- Boltzmann energy */
      c0.oldEngy = engy;
      c0.engy = boltzEngy;
      cnfs[i] = move(c0);
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

}
