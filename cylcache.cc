// #define DEBUG

/*

Make a Cylindrical basis with specified parameters.

*/
                                // System libs
#include <unistd.h>
#include <getopt.h>
#include <values.h>

                                // C++/STL headers
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

//
// BOOST stuff
//
#include <boost/program_options.hpp>

namespace po = boost::program_options;


                                // MDW classes
#include <numerical.h>
#include <gaussQ.h>
#include <isothermal.h>
#include <hernquist.h>
#include <model3d.h>
#include <biorth.h>
#include <SphericalSL.H>
#include <interp.h>
#include <EmpCylSL.h>

#include <norminv.H>

#define M_SQRT1_3 (0.5773502691896257645091487)

                                // For debugging
#ifdef DEBUG
#include <fpetrap.h>
#endif

                                // Local headers
//#include "SphericalSL.h"
//#include "DiskHalo4.h" 
#include "localmpi.h"
//#include "ProgramParam.H"
//
// Parameter definition stanza
//


// Global EXP variables

#include <Particle.H>

int VERBOSE        = 4;
int nthrds         = 1;
int this_step      = 0;
unsigned multistep = 0;
unsigned maxlev    = 100;
int mstep          = 1;
int Mstep          = 1;
char threading_on  = 0;
double tpos        = 0.0;
double tnow        = 0.0;

vector<int> stepL(1, 0), stepN(1, 1);
pthread_mutex_t mem_lock;
pthread_mutex_t coef_lock;
string outdir, runtag;


// globalize some stuff
double ASCALE;
double HSCALE;
int MMAX;
int NORDER;
int LMAX;
int NMAX;
int NUMX;
int NUMY;
double RMIN;
double RMAX;
int RNUM;
int PNUM;
int TNUM;
string CACHEFILE;
bool CMAP;
bool LOGR;
bool verbose;
int NUMR;
bool evolved;
  string       dbods;
int ndisk;


//
// Analytic disk density (assuming exponential scaleheight)
//
double DiskDens(double R, double z, double phi)
{
  double ans = 0.0;


    double f = cosh(z/HSCALE);
    ans = exp(-R/ASCALE)/(4.0*M_PI*ASCALE*ASCALE*HSCALE*f*f);
    


  return ans;
}


//
// The conditioning function for the disk
//
double dcond(double R, double z, double phi, int M)
{
  //
  // No shift for M==0
  //
  if (M==0) return DiskDens(R, z, phi);

  //
  // Fold into [-PI/M, PI/M] for M>=1
  //
  double dmult = M_PI/M, phiS;
  if (phi>M_PI)
    phiS = phi + dmult*(int)((2.0*M_PI - phi)/dmult);
  else
    phiS = phi - dmult*(int)(phi/dmult);
  
  //
  // Apply a shift along the x-axis (temporarily disabled)
  //
  double x = R*cos(phiS);// - ASHIFT*ASCALE;
  double y = R*sin(phiS);
  return DiskDens(sqrt(x*x + y*y), z, atan2(y, x));
}



int 
main(int argc, char **argv)
{
  //====================
  // Inialize MPI stuff
  //====================

  local_init_mpi(argc, argv);

  /*
  //====================
  // Parse command line 
  //====================
  */

  po::options_description desc("Allowed options");
  desc.add_options()
    ("help,h", "produce this help message")
    ("verbose,v", "verbose output")
    
    ("ascale", po::value<double>(&ASCALE)->default_value(0.01), 
     "scalelength for basis construction (1.41*desired)")
    
    ("hscale", po::value<double>(&HSCALE)->default_value(0.001), 
     "scaleheight for basis construction")
    
    ("mmax", po::value<int>(&MMAX)->default_value(6), 
     "number of harmonics")
    
    ("norder", po::value<int>(&NORDER)->default_value(12),
     "number of radial functions")

    ("lmax", po::value<int>(&LMAX)->default_value(36),
     "number of harmonics to construct initial spherical basis model")

    ("nmax", po::value<int>(&NMAX)->default_value(36),
     "number of radial harmonics to construct initial spherical basis model")

    ("numx", po::value<int>(&NUMX)->default_value(128),
     "number of radial points in eof grid")

    ("numy", po::value<int>(&NUMY)->default_value(64),
     "number of vertical points in eof grid")

    ("numr", po::value<int>(&NUMR)->default_value(200),
     "number of entries in generated spherical radial basis table")
    
    ("rmin", po::value<double>(&RMIN)->default_value(0.001), 
     "minimum r value (multiples of ascale)")
    
    ("rmax", po::value<double>(&RMAX)->default_value(20.0), 
     "maximum r value (multiples of ascale)")

    ("rnum", po::value<int>(&RNUM)->default_value(200), 
     "spherical expansion radial function number")

    ("pnum", po::value<int>(&PNUM)->default_value(80), 
     "expansion phi function number")

    ("tnum", po::value<int>(&TNUM)->default_value(80), 
     "expansion theta function number")

    ("eofcache", po::value<string>(&CACHEFILE)->default_value("eof.cache.file"),
     "cachefile for EOF")

    ("cmap", po::value<bool>(&CMAP)->default_value(false),
     "use coordinate mapping")

    ("ndisk",           po::value<int>(&ndisk)->default_value(1000),                    "Number of disk particles")


    ("logr", po::value<bool>(&LOGR)->default_value(false),
     "use logarithmic radius")

    ("evolved",
    po::value<bool>(&evolved)->default_value(false), "Use existing disk body file given by <dbods> and do not create a new disk")

     ("dbods",           po::value<string>(&dbods)->default_value("disk.bods"),          "Disk particle output file")

    
    ;
  
  po::variables_map vm;

  try {
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);    
  } catch (po::error& e) {
    std::cout << "Option error: " << e.what() << std::endl;
    exit(-1);
  }

  if (vm.count("help")) {
    cout << desc << "\n";
    return 1;
  }
 
  if (vm.count("verbose")) verbose = true;

    // Vectors to contain phase space Particle structure is defined in
  // Particle.H
  //
vector<Particle> dparticles;

MPI_Barrier(MPI_COMM_WORLD);

 int n_particlesD;

  double disk_mass1 = 0.025;
  double disk_mass2 = 0.005;

  // ------------------------------------------------------------------------
  if (evolved) {		// ---------------------------
				// Use existing halo body file
    std::ifstream hin(dbods);	// ---------------------------
    
    if (hin) {
      std::string line;
      std::getline(hin, line);
      std::istringstream sin(line);

      int niatr, ndatr;
      sin >> ndisk;
      sin >> niatr;
      sin >> ndatr;
      
      // Divvy up the particles by core.  The root node gets any
      // remainder.
      //
      n_particlesD = ndisk/numprocs;

      int ibeg = 0;
      int iend = ndisk - n_particlesD*(numprocs-myid-1);
      
      if (myid>0) {
	ibeg = ndisk - n_particlesD*(numprocs-myid);
	for (int i=0; i<ibeg; i++) std::getline(hin, line);
      }

      Particle P(niatr, ndatr);

      for (int i=ibeg; i<iend; i++) {
	std::getline(hin, line);
	std::istringstream sin(line);
	sin >> P.mass;
	for (int k=0; k<3; k++)     sin >> P.pos[k];
	for (int k=0; k<3; k++)     sin >> P.vel[k];
	for (int k=0; k<niatr; k++) sin >> P.iattrib[k];
	for (int k=0; k<ndatr; k++) sin >> P.dattrib[k];
	dparticles.push_back(P);
      }

    } else {
      // Error message
      if (myid==0)
	std::cout << "Could not read halo file <" << dbods
		  << "> . . . quitting" << std::endl;

      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Finalize();
      exit(-1);
    }
  }

  
  //===========================Cylindrical expansion 1 ===========================

    if (myid==0) cout << "Cylindrical expansion . . . " << endl;


    EmpCylSL::RMIN        = RMIN;
    EmpCylSL::RMAX        = RMAX;
    EmpCylSL::NUMX        = NUMX;
    EmpCylSL::NUMY        = NUMY;
    EmpCylSL::NUMR        = NUMR;
    EmpCylSL::CMAP        = CMAP;
    EmpCylSL::VFLAG       = 16;//VFLAG;
    EmpCylSL::logarithmic = LOGR;
    EmpCylSL::DENS        = 1;//DENS;
    //SELECT      = 0;//SELECT; true if signal-to-noise methods are on
    EmpCylSL::CACHEFILE      = CACHEFILE;


    if (myid==0) {
    cout << "Process " << myid << ": "
	 << " rmin=" << EmpCylSL::RMIN
	 << " rmax=" << EmpCylSL::RMAX
	 << " a=" << ASCALE
	 << " h=" << HSCALE
	 << " nmax=" << NMAX
	 << " lmax=" << LMAX
	 << " mmax=" << MMAX
	 << " norder=" << NORDER
	 << endl << flush;
    }


  //if (basis) EmpCylSL::DENS = true;

    /* THESE SHOULD ALL BE SPECIFIED ABOVE
    int NMAX = 36; // Can't go to 72,72 here...
int LMAX = 36;
 int MMAX = 6;
int NORDER = 18;
double ASCALE = 0.0143;
 double HSCALE = 0.001;
    */

  EmpCylSL* expandd1 = NULL;
  
    expandd1 = new EmpCylSL(NMAX, LMAX, MMAX, NORDER, ASCALE, HSCALE);

if (myid==0) cout << "SCALELENGTH1:" << setw(16) << expandd1->get_ascale() << endl;



 if (myid==0) cout << "rnum=" << RNUM << " pnum=" << PNUM << " tnum="  << TNUM << endl;

 // generate the EOF

 if (evolved) {
   if (myid==0) std::cout << "Beginning disk accumulation . . . " << std::flush;
    expandd1->setup_accumulation();

      expand1->setup_eof();
      if (nthrds>1)
	expand1->accumulate_eof_thread(dparticles, true);
      // report = true means we see the EOF computation progress
      else
	expand1->accumulate_eof(dparticles, true);
      MPI_Barrier(MPI_COMM_WORLD);

      if (myid==0) std::cout << "done" << std::endl;
  
      if (myid==0) std::cout << "Making the EOF . . . " << std::flush;
      expand1->make_eof();
      MPI_Barrier(MPI_COMM_WORLD);

      //      if (myid==0) expandd1->cache_grid(1);
    
 } else {
    if (myid==0) std::cout << "Building a conditioned disc . . . " << std::flush;

	expandd1->generate_eof(RNUM, PNUM, TNUM, dcond);
 }
	


    if (myid==0) cout << "done!" <<endl;
    
 

    //====================WRAP IT ALL UP========================

  MPI_Barrier(MPI_COMM_WORLD);

  delete expandd1;

  MPI_Finalize();

  return 0;
  }

