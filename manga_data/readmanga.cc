/*
read in processed MaNGA pixel data

compile string: 
clang++ readmanga.cc -o readmanga

MSP 13 Jan 2021 prototype version

python clone for visualising:

import numpy as np
import matplotlib.pyplot as plt

def read_manga_binary(infile):
    """read in manga binary file"""
    npix = np.memmap(infile,dtype='int32',shape=(1))[0]
    tmp = np.memmap(infile,dtype='float64',shape=(5,npix),offset=4)
    
    xcoord = tmp[0]
    ycoord = tmp[1]
    flux   = tmp[2]
    gvel   = tmp[3]
    svel   = tmp[4]
    
    return xcoord,ycoord,flux,gvel,svel

infile = 'manga-8979-6102_pixels.dat'
xcoord,ycoord,flux,gvel,svel = read_manga_binary(infile)

plt.scatter(xcoord,ycoord,color=cm.RdBu_r((svel+125.)/250.,1.))


*/


#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>


using namespace std;

struct MangaData
{
  int NPIX;            // the number of recorded pixels
  
  vector<double> xcoord;    // the x coordinates  
  vector<double> ycoord;    // the x coordinates  
  vector<double> flux;      // the measured flux
  vector<double> gvel;      // the gas velocity  
  vector<double> svel;      // the stellar velocity  

};


void read_manga_file (string& manga_file, MangaData& MD) {

  /*
    read in the processed manga pixel data

  */
  ifstream in(manga_file.c_str());
  if (!in) {
        cout << "readmanga::read_manga_file: Unable to open file!\n";
        exit(1);
  }
  
  // first thing in is NPIX
  in.read((char *)&MD.NPIX, sizeof(int));

  double tmp;

  cout << "read_manga_file: reading pixels from file . . . ";
  
  // order in file is xcoord,ycoord,flux,gvel,svel
  for (int x=0;x<MD.NPIX;x++) {
    in.read((char *)&tmp, sizeof(double));
    MD.xcoord.push_back(tmp);
  }

  for (int x=0;x<MD.NPIX;x++) {
    in.read((char *)&tmp, sizeof(double));
    MD.ycoord.push_back(tmp);
  }

  for (int x=0;x<MD.NPIX;x++) {
    in.read((char *)&tmp, sizeof(double));
    MD.flux.push_back(tmp);
  }

  for (int x=0;x<MD.NPIX;x++) {
    in.read((char *)&tmp, sizeof(double));
    MD.gvel.push_back(tmp);
  }

  for (int x=0;x<MD.NPIX;x++) {
    in.read((char *)&tmp, sizeof(double));
    MD.svel.push_back(tmp);
  }
  
  cout << "success!!" << endl;

}


int main () {

  string manga_file = "manga-8979-6102_pixels.dat";

  MangaData MD;

  read_manga_file(manga_file, MD);

  // how many pixels were read in?
  cout << "n_pixels read: " << MD.NPIX << endl;

  // validate the stellar velocities
  cout << "stellar velocity spot checks: " << endl;
  for (int x=0;x<10;x++) cout << MD.svel[x] << endl;
    

}
