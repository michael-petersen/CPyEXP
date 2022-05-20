#include <algorithm>
#include <iostream>
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
  std::vector<int> buffer(bufsize);

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
  const std::string filename("test.file");
  MPI_File fh;
  MPI_File_open(MPI_COMM_WORLD, filename.c_str(),
		MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);

  // Position the file pointer
  //
  MPI_File_seek(fh, rank * bufsize, MPI_SEEK_SET);

  // Do the write
  //
  MPI_Status status;
  MPI_File_write(fh, buffer.data(), nints, MPI_INT, &status);

  int count;
  MPI_Get_count( &status, MPI_INT, &count );
  if (count != nints) {
    std::cout << "Wrong count [" << rank << "]: " << count  << std::endl;
  }


  // Done
  //
  MPI_File_close(&fh);

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
