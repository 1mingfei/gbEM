#include "gbCnf.h"
#include "EMHome.h"

typedef std::chrono::high_resolution_clock Clock;

int main(int argc, char* argv[]) {
  MPI_Init(&argc, &argv);
  int me;
  MPI_Comm_rank(MPI_COMM_WORLD, &me);
  auto t1 = Clock::now();
  EMHome* em = new EMHome(argc, argv);
  if (em) delete em;
  auto t2 = Clock::now();
  MPI_Barrier(MPI_COMM_WORLD);
  if (me == 0) {
    cout << "Delta t2-t1: "
         << std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count()
         << " seconds" << endl;
  }

  MPI_Finalize();
  return (0);
}
