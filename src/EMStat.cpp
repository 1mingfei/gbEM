#include "gbCnf.h"

/*gather energies from all the structures
 *and calculate probability */
void EMHome::getProb(gbCnf& cnfModifier, double T, double N)
{
  MPI_Barrier(MPI_COMM_WORLD);
  if (me == 0) std::cout << "processing Boltzmann probability\n";
  for (int i = 0; i < NI; ++i)
  {
    if (i%nProcs != me) continue;
    c0 = std::move(cnfModifier.readLmpData(to_string(i) + ".txt"));
    engy = c0.engy / N;
    boltzEngy = exp(-engy/KB/T);
    cout << "read in config " << i << " on processor " << me 
         << " out of total " << nProcs  << " processors, energy: "
         << engy << " BE " << boltzEngy << endl;
  }

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Gather(&boltzEngy, 1, MPI_DOUBLE, engys, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  if (me == 0)
  {
    for (int i = 0; i < NI; ++i)
      cout << engys[i] << "   ";
    cout << endl;
    for (int i = 0; i < NI; ++i)
      sumEngy += engys[i];
  }
  MPI_Bcast(&sumEngy, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  cout << "rank " << me << " sum energy: " << sumEngy << endl;

}
