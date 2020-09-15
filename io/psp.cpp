/*
Working compile string:


// specify all libraries and includes
clang++ -Ofast  -v -I/usr/local/boost_1_61_0/include -L/usr/local/boost_1_61_0/lib  -L/usr/local/lib -I/usr/local/include -lboost_program_options psp.cpp -o psp

boost needs .cpp

//http://www.cplusplus.com/doc/tutorial/files/

*/
#include <boost/program_options.hpp>

#include <fstream>
#include <string>
#include <array>
#include <cmath>
#include <iostream>

using namespace std;

namespace po = boost::program_options;

const static unsigned long magic = 0xadbfabc0;
const static unsigned long mmask = 0xf;
const static unsigned long nmask = ~mmask;


void read_master_header(ifstream* pspfile) {

    streampos size;
  char * memblock;
  
    size = pspfile->tellg();
    //memblock = new char [size];
    pspfile->seekg (0, ios::beg);
    //pspfile.read (memblock, size);
    //pspfile.close();

    cout << "the entire file content is in memory\n";

    //delete[] memblock;
}


int main(int argc, char** argv)
{

  std::string pspfilename;

  po::options_description desc(
"This routine reads in PSP files for manipulation.\n\nCommand line options");

  desc.add_options()
    ("pspfile,p",      po::value<std::string>(&pspfilename)->default_value("/Users/mpetersen/exptest/OUT.run0.00000"), "input file")
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

  ifstream pspfile (pspfilename.c_str(), ios::in|ios::binary|ios::ate);
  if (pspfile.is_open())
  {
    read_master_header(&pspfile);
  }
  
  else cout << "Unable to open file";
  return 0;
}

 
