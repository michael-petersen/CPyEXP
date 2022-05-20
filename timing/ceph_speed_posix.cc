#include <algorithm>
#include <iostream>
#include <sstream>
#include <fstream>
#include <chrono>
#include <random>
#include <vector>
#include <limits>

#include <mpi.h>
const int MegaByte = 1024 * 1024;
const int FileSize = 1024 * MegaByte;

int main(int argc, char **argv)
{
  // Initialize MPI
  MPI_Init(&argc, &argv);

  // Preliminaries: get communicator rank and size
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  int bufsize = FileSize/size;
  int nints = bufsize/sizeof(int);

  // Make a buffer of ints
  std::vector<int> buffer(nints);

  // Fill it with random values
  std::random_device rd;
  std::mt19937 mte(rd());

  const int maxint = std::numeric_limits<int>::max();
  std::uniform_int_distribution<int> dist(-maxint, maxint);

  std::generate(buffer.begin(), buffer.end(), [&] () { return dist(mte); });


  // Okay now start timing

  using unit_t = std::chrono::nanoseconds;
  auto start_time = std::chrono::steady_clock::now();

  // Create and open the file
  //
  std::ostringstream ostr;
  ostr << "test.file." << rank;
  std::ofstream out(ostr.str());

  // Write the blob
  out.write((const char *)buffer.data(), nints*sizeof(int));

  // Done
  out.close();

  // Get elapsed time
  //
  auto end_time = std::chrono::steady_clock::now();

  auto duration = 1.0e-9 *
    std::chrono::duration_cast<unit_t>(end_time - start_time).count();

  if (rank==0) {
    std::cout << "Time to complete: " << duration << std::endl;
    std::cout << "Rate: " << static_cast<double>(FileSize)/MegaByte/duration
	      << " MB/s" << std::endl;
  }

  MPI_Finalize();
  return 0;
}
