// #define DEBUG

/*

The most basic, hard-wired initial condition generation.
MSP 12.30.16

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

                                // MDW classes
#include <numerical.h>
#include <gaussQ.h>
#include <isothermal.h>
#include <hernquist.h>
#include <model3d.h>
#include <biorth.h>
#include <SphericalSL.h>
#include <interp.h>
#include <EmpCylSL.h>
#include <exponential3.h>

// what does this do? Why does it break against boost?
#include <norminv.H>

#define M_SQRT1_3 (0.5773502691896257645091487)

                                // For debugging
#ifdef DEBUG
#include <fpetrap.h>
#endif

                                // Local headers
#include "SphericalSL.h"
#include "DiskHalo5.h" 
#include "localmpi.h"
//#include "ProgramParam.H"

                                // Also test out boost
//#include <boost/program_options.hpp>



//
// Parameter definition stanza
//


// Global variables

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

// block these out for global availability
//EmpCylSL* expandd1 = NULL;
//SphericalSL *expandh1 = NULL;

//namespace po = boost::program_options;



// I don't like that these are global: move to disk class?
// TODO
double ASCALE = 0.0143;
double HSCALE = 0.001;
double ToomreQ = 0.9;
double disk_mass1 = 0.025;

//double disk_mass1 = 0.0125;
//double disk_mass1 = 0.00625;


//double ASCALE, HSCALE, ToomreQ, disk_mass1;

//const string PACKAGE_STRING = "basicICs for Disk and Halo generation";


//
// Analytic disk density (assuming exponential scaleheight)
//
//double DiskDens(double R, double z, double phi, double ASCALE, double HSCALE)
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

  // double ASCALE = expandd1->get_ascale();
  //double HSCALE = expandd1->get_hscale();
  
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


bool expcond = 1;
int RNUM = 200;
int PNUM = 80;
int TNUM = 80;


int 
main(int argc, char **argv)
{
  //====================
  // Inialize MPI stuff
  //====================

  local_init_mpi(argc, argv);


  //====================
  // Parse command line 
  //====================

  /*

  int for a discalo1, ndisk1;
  //double ASCALE,HSCALE,ToomreQ,disk_mass1;
   
  // Setup options.
  po::options_description desc("Allowed Options");
  desc.add_options()
    ("nhalo1", po::value<int>(&nhalo1)->default_value(100), "Number of halo particles")
    ("ndisk1", po::value<int>(&ndisk1)->default_value(100), "Number of disk particles")
    //("ascale", po::value<double>(&ASCALE)->default_value(0.0143), "Scalelength for exponential disk")
    //("hscale", po::value<double>(&HSCALE)->default_value(0.001), "Scaleheight for exponential disk")
    //("toomreQ", po::value<double>(&ToomreQ)->default_value(0.9), "Toomre Q for disk")
    //("disk_mass", po::value<double>(&disk_mass1)->default_value(0.025), "Disk Mass")
    //("help,h",       "this help message");
    ;

  po::variables_map vm;

  // nice little error handler
try {
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);    
  } catch (po::error& e) {
    cout << "Option error: " << e.what() << std::endl;
    exit(-1);
  }


//  if (vm.count("help")) {
//    std::cout << desc << "\n";
//    return 1;
//  }

*/

  /*
  double n0, tol;
  unsigned T;
  int niter;
  std::string outf;

  po::options_description desc("Allowed options");
  desc.add_options()
    ("help,h", "produce this help message")
    ("density,D", po::value<double>(&n0)->default_value(1.0e-4), 
     "Density in amu/cc. Good for n0<8.5e-2")
    ("temp,T", po::value<unsigned>(&T)->default_value(25000), 
     "Temperature in Kelvin (integer value)")
    ("tol,e",     po::value<double>(&tol)->default_value(1.0e-10), 
     "error tolerance")
    ("iter,n",    po::value<int>(&niter)->default_value(1000), 
     "maximum number of iterations")
    ("outfile,o", po::value<std::string>(&outf)->default_value("IonRecombFrac.data"),
     "data file for makeIon input")
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
    std::cout << desc << "\n";
    return 1;
  }

  */

#ifdef DEBUG                    // For gdb . . . 
  sleep(20);
  set_fpu_handler();            // Make gdb trap FPU exceptions
#endif


  
  int nhalo1 = 10000000;
  // int nhalo2 = 10000;
   int ndisk1 = 1000000;

  int n_particlesH1;
  int n_particlesD1;

                       // Divvy up the particles
  n_particlesH1 = nhalo1/numprocs;
  if (myid==0) n_particlesH1 = nhalo1 - n_particlesH1*(numprocs-1);
  
  // n_particlesH2 = nhalo2/numprocs;
  //if (myid==0) n_particlesH2 = nhalo2 - n_particlesH2*(numprocs-1);

    n_particlesD1 = ndisk1/numprocs;
  if (myid==0) n_particlesD1 = ndisk1 - n_particlesD1*(numprocs-1);


#ifdef DEBUG  
  cout << "Processor " << myid << ": n_particlesH1=" << n_particlesH1 << "\n";
#endif


                                // Vectors to contain phase space
                                // Particle structure is defined in
                                // Particle.H
  vector<Particle> h1particles, d1particles;


  // set some DiskHalo Parameters
  DiskHalo::Q           = ToomreQ;
  DiskHalo::COMPRESSION = 1.0;


  //===========================Spherical expansion 1=============================

  if (myid==0) cout << "Spherical expansion 1..." << endl;

  int LMAX = 6;
  int NMAX = 20;
  double SCSPH = 0.0667;

  SphericalSL::RMIN = 0.001; //RMIN;
  SphericalSL::RMAX = 2.0; //RSPHSL;
  SphericalSL::NUMR = 4000; //NUMR;
  int nthrds = 1; // new  compatibility version
                                // Create expansion only if needed . . .
  SphericalSL *expandh1 = NULL;
  if (n_particlesH1) {
    expandh1 = new SphericalSL(nthrds,LMAX, NMAX, SCSPH);
#ifdef DEBUG
    string dumpname("debug");
    expandh1->dump_basis(dumpname);
#endif
  }

  SphericalModelTable *halo1;
  halo1 = new SphericalModelTable("SLGridSph.model", 0, 1);

  
MPI_Barrier(MPI_COMM_WORLD);

//double disk_mass1 = 0.025;

    

  
  //===========================Cylindrical expansion 1 ===========================

    if (myid==0) cout << "Cylindrical expansion 1 . . . " << endl;


  EmpCylSL::RMIN        = 0.001; //RCYLMIN;
  EmpCylSL::RMAX        = 20.; //RCYLMAX;
  EmpCylSL::NUMX        = 128; // NUMX;
  EmpCylSL::NUMY        = 64;//NUMY;
  EmpCylSL::NUMR        = 200;//NUMR;
  EmpCylSL::CMAP        = 1;//CMAP;
  EmpCylSL::VFLAG       = 16;//VFLAG;
  EmpCylSL::logarithmic = 1;//LOGR;
  EmpCylSL::DENS        = 1;//DENS;
  // EmpCylSL::SELECT      = 0;//SELECT;

    EmpCylSL::CACHEFILE      = ".eof.cache.file1";


  //if (basis) EmpCylSL::DENS = true;

                                // Create expansion only if needed . . .

    int NMAX2 = 36; // Can't go to 72,72 here...
int LMAX2 = 36;
 int MMAX = 6;
int NORDER = 18;
double ASCALE = 0.0143;
double HSCALE = 0.001;


  EmpCylSL* expandd1 = NULL;
  if (n_particlesD1) {
    expandd1 = new EmpCylSL(NMAX2, LMAX2, MMAX, NORDER, ASCALE, HSCALE);

if (myid==0) cout << "SCALELENGTH1:" << setw(16) << expandd1->get_ascale() << endl;


if (myid==0) {
    cout << "Process " << myid << ": "
	 << " rmin=" << EmpCylSL::RMIN
	 << " rmax=" << EmpCylSL::RMAX
	 << " a=" << expandd1->get_ascale() //ASCALE
	 << " h=" << expandd1->get_hscale() //HSCALE
	 << " nmax2=" << NMAX2
	 << " lmax2=" << LMAX2
	 << " mmax=" << MMAX
	 << " nordz=" << NORDER
	 << endl << flush;
}

 if (myid==0) cout << "rnum=" << RNUM << " pnum=" << PNUM << " tnum="  << TNUM << endl;

    if (expandd1->read_cache() == 0) {
      
      if (expcond) // will only do this if not accumulating from particle distribution!
	expandd1->generate_eof(RNUM, PNUM, TNUM, dcond);
    }

    if (myid==0) cout << "done!" <<endl;
    
  }


MPI_Barrier(MPI_COMM_WORLD);

  //==============================================================================================

 
  DiskHalo *diskhalo;

  string halofile1 =  "SLGridSph.model";
    string halofile2 =  "SLGridSph.fake";

    int DIVERGE=0;
    double DIVERGE_RFAC = 1.0;
    int DIVERGE2=0;
    double DIVERGE_RFAC2 = 1.0;


diskhalo = new DiskHalo(expandh1, expandd1,  
			//HSCALE,  0.7*ASCALE,  disk_mass1,
			expandd1->get_hscale(),0.7*expandd1->get_ascale(), disk_mass1, // hard-wired version for testing
	 halofile1, 0, 1,
	halofile2, 0, 1,
	 DiskHalo::Asymmetric);

	   if (myid==0) cout << "DONE GENERATING ONE-DISK MULTIMASS HALO" << endl;

	   MPI_Barrier(MPI_COMM_WORLD);

  //==============================================================================================
	   // set zero flag...not sure if this is needed
	   bool zero = 1;
  diskhalo->zero_com(zero);

  // this is a problem currently...
  //diskhalo->zero_cov(zero);

  //==================================HOUSEKEEPING=========================================
	   
                                // Open output files (make sure it exists
                                // before realizing a large phase space)
	   ofstream out_halo, out_disk1;

	   string hbods, dbods;

  hbods = "halo.input";
  dbods = "stellar.input";
  
  
  if (myid==0) {
    out_halo.open(hbods.c_str());
    if (!out_halo) {
      cout << "Could not open <" << hbods << "> for output\n";
      MPI_Abort(MPI_COMM_WORLD, 4);
      exit(0);
    }

    out_disk1.open(dbods.c_str());
    if (!out_disk1) {
      cout << "Could not open <" << dbods << "> for output\n";
      MPI_Abort(MPI_COMM_WORLD, 4);
      exit(0);
    }


  }
  //==================================HALO 1 COORDINATES============================================

 
	   // should be unchanged

    if (n_particlesH1) {
      if (myid==0) cout << "Generating halo phase space . . . " << flush;
      
      diskhalo->set_halo(h1particles, nhalo1, n_particlesH1);
      
    MPI_Barrier(MPI_COMM_WORLD);
    if (myid==0) cout << "done\n";
  }



  //======================================DISK 1 COORDINATES===============================================

	   // set disk1 coordinates

      if (n_particlesD1) {
    if (myid==0) cout << "Generating disk coordinates . . . " << flush;
    
    diskhalo->set_disk_coordinates(diskhalo->disk1, d1particles, ndisk1, n_particlesD1);
    MPI_Barrier(MPI_COMM_WORLD);
    if (myid==0) cout << "done\n";
  }

  
  MPI_Barrier(MPI_COMM_WORLD);
  /*
   //======================================DISK 2 COORDINATES===============================================


      if (n_particlesD2) {
    if (myid==0) cout << "Generating disk coordinates . . . " << flush;
    
    diskhalo->set_disk_coordinates(diskhalo->disk2, d2particles, ndisk2, n_particlesD2);
    MPI_Barrier(MPI_COMM_WORLD);
    if (myid==0) cout << "done\n";
  }

  
  MPI_Barrier(MPI_COMM_WORLD);
  */
   //======================================ACCUMULATE HALO 1===============================================


  if (n_particlesH1) {
    if (myid==0) cout << "Beginning halo accumulation . . . " << flush;
    expandh1->accumulate(h1particles);
    MPI_Barrier(MPI_COMM_WORLD);
    if (myid==0) cout << "done\n";
  }

     //======================================ACCUMULATE DISK 1===============================================

  if (n_particlesD1) {
    if (myid==0) cout << "Beginning disk accumulation I . . . " << flush;

        if (myid==0) cout << "Making disk coefficients I . . . " << flush;
    expandd1->make_coefficients();
    MPI_Barrier(MPI_COMM_WORLD);
    if (myid==0) cout << "done\n";

    if (myid==0) cout << "Reexpand I . . . " << flush;
    expandd1->accumulate(d1particles);
    expandd1->make_coefficients();
    MPI_Barrier(MPI_COMM_WORLD);
    if (myid==0) {
      cout << "done\n"; }

  }

    
  MPI_Barrier(MPI_COMM_WORLD);

  
  //====================Make the phase space velocities========================
  // only for the two disks, the halo is already complete!
  
  if (myid==0) cout << "Generating disk velocities . . . " << flush;
  diskhalo->set_vel_disk(d1particles);
  if (myid==0) cout << "done\n";


  //====================WRITE THE PARTICLES========================
  //forcechange

    if (myid==0) cout << "Writing phase space file . . . " << flush;

  diskhalo->write_disk_file(out_halo, out_disk1, h1particles, d1particles);
  if (myid==0) cout << "done\n";

  out_halo.close();
  out_disk1.close();

    //====================SPACE FOR DIAGNOSTICS========================
  // None of these have been set up for the double disk version yet
  /*

  diskhalo->virial_ratio(hparticles, dparticles);

  ofstream outprof("profile.diag");
  diskhalo->profile(outprof, dparticles, 3.0e-3*ASCALE, 5.0*ASCALE, 100);

  */
    //====================WRAP IT ALL UP========================

  MPI_Barrier(MPI_COMM_WORLD);

    delete expandh1;
  delete expandd1;

  MPI_Finalize();

  return 0;
}

