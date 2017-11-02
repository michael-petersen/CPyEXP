/*
Working compile string:


// specify all libraries and includes
clang++ -Ofast  -v -I/usr/local/boost_1_61_0/include -L/usr/local/boost_1_61_0/lib  -L/usr/local/lib -I/usr/local/include -lboost_program_options halo_utils.c -o halo_utils

boost needs .cpp

*/
#include <boost/program_options.hpp>

#include <fstream>
#include <string>
#include <array>
#include <cmath>
#include <iostream>

using namespace std;



namespace po = boost::program_options;

int main(int argc, char** argv)
{

  double omega, ratio, dt;
  std::string model, cache;

  po::options_description desc(
"This routine reads in the input spherical models to generate basic quantities.\n\nCommand line options");

  desc.add_options()
    ("model,m",      po::value<std::string>(&model)->default_value("SLGridSph.model"), "input sph model")
    ("cache,c",     po::value<std::string>(&cache)->default_value(".sl_grid_sph"), "input sph cache")
    ("help,h",       "this help message");

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);
  
  if (vm.count("help")) {
    std::cout << desc << std::endl;
    return 0;
  }

  /*
    Structure of cache file

   (f, dtype=np.uint32,count=4)

  */

  streampos size;
  char * memblock;

  ifstream file (model.c_str(), ios::in|ios::binary|ios::ate);
  if (file.is_open())
  {
    size = file.tellg();
    memblock = new char [size];
    file.seekg (0, ios::beg);
    file.read (memblock, size);
    file.close();

    cout << "the entire file content is in memory\n";

    delete[] memblock;
  }
  
  else cout << "Unable to open file";
  return 0;
}

 
