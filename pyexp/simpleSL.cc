/*****************************************************************************
*  Description:
*  -----------
*


***************************************************************************/

  // C++/STL headers
#include <filesystem>
#include <iostream>
#include <cstdlib>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cmath>
#include <string>

using namespace std;

  // System libs
#include <sys/time.h>
#include <sys/resource.h>

  // MDW classes
#include <numerical.H>
#include <ParticleReader.H>
#include <interp.H>
#include <massmodel.H>
#include <SphereSL.H>
#include <DataGrid.H>

// Library support
#include <localmpi.H>
#include <foarray.H>
#include <EXPini.H>
//#include <global.H>
#include <KDtree.H>
//#include <libvars.H>

//using namespace __EXP__;

// generic header includes: requires both pybind11 and eigen
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>
namespace py = pybind11;

// eigen includes
#include <Eigen/StdVector>
#include <Eigen/Dense>
using Eigen::MatrixXd;

// Globals
//
std::string OUTFILE;
double RMIN, RMAX;
int OUTR, LMAX, NMAX, L1, L2, N1, N2;
bool VOLUME, SURFACE, PROBE;


std::vector<Eigen::MatrixXd> slfunctions(string MODFILE,int LMAX=4,int NMAX=10, double logxmin=-3.0, double logxmax=0.5, int numr=2000)
{

  auto halo = std::make_shared<SphericalModelTable>(MODFILE);

  SphereSL::mpi  = false;
  SphereSL::NUMR = numr;
  double rscale  = 1.;

  SphereSL ortho(halo, LMAX, NMAX, 0, rscale, true);

  double dx = (logxmax-logxmin)/numr;

  std::vector<Eigen::MatrixXd> slfunctions;
  MatrixXd tab;

  slfunctions.resize((LMAX+1)*(LMAX+1));
  for (int l=0;l<((LMAX+1)*(LMAX+1));l++) slfunctions[l].resize(NMAX,numr);

  double xval;
  for (int ix=0;ix<numr;ix++) {

    xval = pow(10,logxmin+ix*dx);
    tab = ortho.get_pot(xval);

    for (int l=0;l<((LMAX+1)*(LMAX+1));l++) {
      for (int n=0;n<NMAX;n++){
        slfunctions[l](n,ix) = tab(l,n);
      }
    }
  }

  return slfunctions;

}


Eigen::MatrixXd
coefsl(std::vector<double> mass,
       std::vector<double> xpos, std::vector<double> ypos, std::vector<double> zpos,
       string MODFILE,
       int LMAX=4,
       int NMAX=10)
{

  // count the number of particles
  int nparticles = xpos.size();

  auto halo = std::make_shared<SphericalModelTable>(MODFILE);

  SphereSL::mpi  = false;
  SphereSL::NUMR = 4000;
  double rscale  = 1.;

  SphereSL ortho(halo, LMAX, NMAX, 0, rscale, true);

  // set up for accumulation
  ortho.reset_coefs();

  // loop through particles
  for (int i=0;i<nparticles;i++) {
     ortho.accumulate(xpos[i],ypos[i],zpos[i],mass[i]);
  }

  std::cout << "done" << std::endl;

  //------------------------------------------------------------

  std::cout << "Making coefficients . . . " << std::flush;
  ortho.make_coefs();
  std::cout << "done" << std::endl;

  // reshape the ortho.expcoef into a MatrixXd
  // allocate workspace
  MatrixXd outcoefs;

  outcoefs = ortho.retrieve_coefs();

  return outcoefs;
}


PYBIND11_MODULE(simpleSL, m) {
    m.doc() = "basic SL capabilities";

    m.def("coefsl", &coefsl, "Accumulate particles onto SL",
      py::arg("mass"),
      py::arg("x"),
      py::arg("y"),
      py::arg("z"),
      py::arg("modelfile"),
      py::arg("LMAX")    = 6,
      py::arg("NMAX")    = 20);

    m.def("slfunctions", &slfunctions, "Basic SLfunction exposure",
      py::arg("modelfile"),
      py::arg("LMAX")    = 6,
      py::arg("NMAX")    = 20,
      py::arg("logxmin") = -3,
      py::arg("logxmax") = 0.5,
      py::arg("numr")    = 2000);
}
