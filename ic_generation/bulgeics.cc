// #define DEBUG

/*

Test set-up of two simultaneous spherical bases

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
#include <EmpOrth9thd.h>

#include <norminv.H>

#define M_SQRT1_3 (0.5773502691896257645091487)

                                // For debugging
#ifdef DEBUG
#include <fpetrap.h>
#endif

                                // Local headers
#include "SphericalSL.h"
#include "DiskHalo6.h" 
#include "localmpi.h"
#include "ProgramParam.H"
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

// I don't like that these are global
// TODO
double ASCALE = 0.0143;
double HSCALE = 0.001;

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

  /*
  //====================
  // Parse command line 
  //====================

  try {
    if (config.parse_args(argc, argv)) return -1;
    param_assign();
  }
  catch (const char *msg) {
    cerr << msg << endl;
    return -1;
  }

  */

#ifdef DEBUG                    // For gdb . . . 
  sleep(20);
  set_fpu_handler();            // Make gdb trap FPU exceptions
#endif

  
  int nhalo1 = 10000000;
  int nhalo2 = 100000;
  int ndisk1 = 1000000;

  int n_particlesH1;
  int n_particlesH2;
  int n_particlesD1;

                       // Divvy up the particles
  n_particlesH1 = nhalo1/numprocs;
  if (myid==0) n_particlesH1 = nhalo1 - n_particlesH1*(numprocs-1);
  
  n_particlesH2 = nhalo2/numprocs;
  if (myid==0) n_particlesH2 = nhalo2 - n_particlesH2*(numprocs-1);

    n_particlesD1 = ndisk1/numprocs;
  if (myid==0) n_particlesD1 = ndisk1 - n_particlesD1*(numprocs-1);


#ifdef DEBUG  
  cout << "Processor " << myid << ": n_particlesH1=" << n_particlesH1 << "\n";
#endif

  /*
  if (n_particlesH + n_particlesD + n_particlesG <= 0) {
    if (myid==0) cout << "You have specified zero particles!\n";
    MPI_Abort(MPI_COMM_WORLD, 3);
    exit(0);
  }
  */
                                // Vectors to contain phase space
                                // Particle structure is defined in
                                // Particle.H
  vector<Particle> h1particles, h2particles, d1particles;



  //===========================Spherical expansion 1=============================

  if (myid==0) cout << "Spherical expansion 1..." << endl;

  int LMAX = 6;
  int NMAX = 20;
  double SCSPH = 0.0667;

  SphericalSL::RMIN = 0.00005; //RMIN;
  SphericalSL::RMAX = 2.0; //RSPHSL;
  SphericalSL::NUMR = 4000; //NUMR;
  SLGridSph::sph_cache_name = ".slgrid_halo_cache";
                                // Create expansion only if needed . . .
  SphericalSL *expandh1 = NULL;
  if (n_particlesH1) {
    expandh1 = new SphericalSL(LMAX, NMAX, SCSPH);
#ifdef DEBUG
    string dumpname("debug");
    expandh1->dump_basis(dumpname);
#endif
  }

  SphericalModelTable *halo1;
  halo1 = new SphericalModelTable("SLGridSph.model", 0, 1);

  
MPI_Barrier(MPI_COMM_WORLD);


 
  //===========================Spherical expansion 2=============================

  if (myid==0) cout << "Spherical expansion 2..." << endl;

  // keep same parameters for bulge
  // int LMAX = 6;
  //int NMAX = 20;
  SCSPH = 0.5;

  SphericalSL::RMIN = 0.00005; //RMIN;
  SphericalSL::RMAX = 0.1; //RSPHSL;
  SphericalSL::NUMR = 4000; //NUMR;
  SLGridSph::sph_cache_name = ".slgrid_bulge_cache";
                                // Create expansion only if needed . . .
  SphericalSL *expandh2 = NULL;
  if (n_particlesH2) {
    expandh2 = new SphericalSL(LMAX, NMAX, SCSPH);
#ifdef DEBUG
    string dumpname("debugbulge");
    expandh2->dump_basis(dumpname);
#endif
  }

  SphericalModelTable *halob;
  halob = new SphericalModelTable("SLGridBULGE.model", 0, 1);

  
MPI_Barrier(MPI_COMM_WORLD);




 

  
  //===========================Cylindrical expansion 1 ===========================

    if (myid==0) cout << "Cylindrical expansion 1 . . . " << endl;


        double disk_mass1 = 0.025;


  EmpCylSL::RMIN        = 0.0001; //RCYLMIN;
  EmpCylSL::RMAX        = 20.; //RCYLMAX;
  EmpCylSL::NUMX        = 128; // NUMX;
  EmpCylSL::NUMY        = 64;//NUMY;
  EmpCylSL::NUMR        = 1000;//NUMR;
  EmpCylSL::CMAP        = 1;//CMAP;
  EmpCylSL::VFLAG       = 16;//VFLAG;
  EmpCylSL::logarithmic = 1;//LOGR;
  EmpCylSL::DENS        = 1;//DENS;
  EmpCylSL::SELECT      = 0;//SELECT;

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
	 << " a=" << ASCALE
	 << " h=" << HSCALE
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
    string halofile3 = "SLGridBULGE.model";

    int DIVERGE=0;
    double DIVERGE_RFAC = 1.0;
    int DIVERGE2=0;
    double DIVERGE_RFAC2 = 1.0;


    diskhalo = new DiskHalo(expandh1, expandh2, expandd1,  
	   HSCALE,  0.7*ASCALE,  disk_mass1,
	 halofile1, 0, 1,
	halofile2, 0, 1,
			halofile3, 0, 1,
	 DiskHalo::Asymmetric);

	   if (myid==0) cout << "DONE GENERATING ONE-DISK MULTIMASS HALO WITH BULGE" << endl;

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
  ofstream out_halo, out_bulge, out_disk1;

  string hbods, bbods, dbods;

  hbods = "halo.input";
  bbods = "bulge.input";
  dbods = "stellar.input";
  
  
  if (myid==0) {
    out_halo.open(hbods.c_str());
    if (!out_halo) {
      cout << "Could not open <" << hbods << "> for output\n";
      MPI_Abort(MPI_COMM_WORLD, 4);
      exit(0);
    }

    out_bulge.open(bbods.c_str());
    if (!out_bulge) {
      cout << "Could not open <" << bbods << "> for output\n";
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
  //==================================HALO COORDINATES============================================

 
	   // should be unchanged

    if (n_particlesH1) {
      if (myid==0) cout << "Generating halo phase space . . . " << flush;
      
      diskhalo->set_halo(h1particles, nhalo1, n_particlesH1);
      
    MPI_Barrier(MPI_COMM_WORLD);
    if (myid==0) cout << "done\n";
  }

  //==================================BULGE COORDINATES============================================

 

    if (n_particlesH2) {
      if (myid==0) cout << "Generating bulge phase space . . . " << flush;


      //What needs to happen: bulge is used for mass and density, then composite potential is used (halo3).
      diskhalo->set_bulge_coordinates(h2particles, nhalo2, n_particlesH2);
      
    MPI_Barrier(MPI_COMM_WORLD);
    if (myid==0) cout << "done\n";
  }


    

  //======================================DISK 1 COORDINATES===============================================

	   // set disk1 coordinates

      if (n_particlesD1) {
    if (myid==0) cout << "Generating disk coordinates . . . " << flush;
    
    diskhalo->set_disk_coordinates(diskhalo->disk, d1particles, ndisk1, n_particlesD1);
    MPI_Barrier(MPI_COMM_WORLD);
    if (myid==0) cout << "done\n";
  }

  
  MPI_Barrier(MPI_COMM_WORLD);

  
   //======================================ACCUMULATE HALO 1===============================================


  if (n_particlesH1) {
    if (myid==0) cout << "Beginning halo accumulation . . . " << flush;
    expandh1->accumulate(h1particles);
    MPI_Barrier(MPI_COMM_WORLD);
    if (myid==0) cout << "done\n";
  }

    
   //======================================ACCUMULATE BULGE===============================================


  if (n_particlesH2) {
    if (myid==0) cout << "Beginning bulge accumulation . . . " << flush;
    expandh2->accumulate(h2particles);
    MPI_Barrier(MPI_COMM_WORLD);
    if (myid==0) cout << "done\n";
  }


     //======================================ACCUMULATE DISK 1===============================================

  if (n_particlesD1) {
    if (myid==0) cout << "Beginning disk accumulation  . . . " << flush;

        if (myid==0) cout << "Making disk coefficients  . . . " << flush;
    expandd1->make_coefficients();
    MPI_Barrier(MPI_COMM_WORLD);
    if (myid==0) cout << "done\n";

    if (myid==0) cout << "Reexpand  . . . " << flush;
    expandd1->accumulate(d1particles);
    expandd1->make_coefficients();
    MPI_Barrier(MPI_COMM_WORLD);
    if (myid==0) {
      cout << "done\n"; }

  }

    
  MPI_Barrier(MPI_COMM_WORLD);

  
  //====================Make the phase space velocities========================
  // only for the bulge and disk , the halo is already complete!

   if (myid==0) cout << "Generating bulge velocities . . . " << flush;
  diskhalo->set_vel_bulge(h2particles);
  if (myid==0) cout << "done\n";

  
  if (myid==0) cout << "Generating disk velocities . . . " << flush;
  diskhalo->set_vel_disk(d1particles);
  if (myid==0) cout << "done\n";


  //====================WRITE THE PARTICLES========================
  //forcechange

    if (myid==0) cout << "Writing phase space file . . . " << flush;

    diskhalo->write_file(out_halo, out_disk1, out_bulge, h1particles, d1particles, h2particles);
  if (myid==0) cout << "done\n";

  out_halo.close();
  out_disk1.close();
  out_bulge.close();

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
    delete expandh2;
  delete expandd1;

  MPI_Finalize();

  return 0;
  }

