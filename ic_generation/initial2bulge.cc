/*
  Generates a Monte Carlo realization of a halo with an embedded
  disk using Jeans' equations.

  Assumptions:

  1) Spherical halo supplied in an input table

  2) Axisymmetric (but 3-dimensional) exponential disk with sech^2(Z/z)
     vertical profile

  3) Halo as spherical velocity ellipsoid

  4) Disk as axisymmetric velocity ellipsoid in the plane (that is,
     $\sigma_r = \sigma_\phi$ and $\sigma_z$ determined by solving
     Jeans' equations in cylindrical coordinates.

 Loosely based on Enrico Vesperini's initial.cc and diskANDhalo.cc
 (mostly rewritten)

 Added the basis expansion of the disk: 12/10/01. KHB

 Rewritten and debugged by MDW between 12/28/01-12/31/01.  

        Added command line parsing.  "gendisk -h" will list parameters.

        I removed everything but the "disk and halo" case.  Removed
        multiple inheritance.  Changed interface; expansions registered
        with DiskHalo class on construction

        Switched from biortho classes to expansion classes from the
        EXP code.

        Uses a vector of particle structures rather than a Matrix to
        store an pass phase space.

        Rewrote the particle component foliation code using a more
        general algorithm.

        Solution to Jeans' equations are now computed in parallel and
        tabulated.  

        Particles are stored on local nodes and written to disk by
        master.  

        Removed lots of other cruft.

 More debugging 03/05 by MDW

        Repaired EmpCylSL scaling

	Added additional debugging output

	Compared against expsl routines

	Removed orphaned parameters

 Updated to include gas disk using local Euler solution 04/08 by MDW

 Both constant scale height and isothermal gas disks 08/08 by MDW

 Multimass gas disk 11/08 by MDW

 Updates for constructing disk velocities from an evolved halo 01/20 by MDW

*/
                                // C++/STL headers
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include <fenv.h>

// Boost stuff
//
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

namespace po = boost::program_options;

// MDW classes
//
#include <numerical.h>
#include <gaussQ.h>
#include <isothermal.h>
#include <hernquist.h>
#include <model3d.h>
#include <biorth.h>
#include <SphericalSL.H>
#include <interp.h>
#include <EmpCylSL.h>
#include <DiskModels.H>

#include <norminv.H>

#define M_SQRT1_3 (0.5773502691896257645091487)

                                // For debugging
#ifdef DEBUG
#include <fenv.h>
#include <fpetrap.h>

//===========================================
// Handlers defined in exputil/stack.cc
//===========================================

extern void mpi_print_trace(const string& routine, const string& msg,
			    const char *file, int line);

extern void mpi_gdb_print_trace(int sig);

extern void mpi_gdb_wait_trace(int sig);

//===========================================
// A signal handler to trap invalid FP only
//===========================================

void set_fpu_invalid_handler(void)
{
  // Flag invalid FP results only, such as 0/0 or infinity - infinity
  // or sqrt(-1).
  //
  feenableexcept(FE_INVALID);
  //
  // Print enabled flags to root node
  //
  if (myid==0) {
    const std::list<std::pair<int, std::string>> flags =
      {	{FE_DIVBYZERO, "divide-by-zero"},
	{FE_INEXACT,   "inexact"},
	{FE_INVALID,   "invalid"},
	{FE_OVERFLOW,  "overflow"},
	{FE_UNDERFLOW, "underflow"} };
    
    int _flags = fegetexcept();
    std::cout << "Enabled FE flags: <";
    for (auto v : flags) {
      if (v.first & _flags) std::cout << v.second << ' ';
    }
    std::cout << "\b>" << std::endl;
  }
  signal(SIGFPE, mpi_gdb_print_trace);
}

//===========================================
// A signal handler to produce a traceback
//===========================================

void set_fpu_trace_handler(void)
{
  // Flag all FP errors except inexact
  //
  // fedisableexcept(FE_ALL_EXCEPT & ~FE_INEXACT);

  // Flag invalid FP results only, such as 0/0 or infinity - infinity
  // or sqrt(-1).
  //
  feenableexcept(FE_INVALID);
  //
  // Print enabled flags to root node
  //
  if (myid==0) {
    const std::list<std::pair<int, std::string>> flags =
      {	{FE_DIVBYZERO, "divide-by-zero"},
	{FE_INEXACT,   "inexact"},
	{FE_INVALID,   "invalid"},
	{FE_OVERFLOW,  "overflow"},
	{FE_UNDERFLOW, "underflow"} };
    
    int _flags = fegetexcept();
    std::cout << "Enabled FE flags: <";
    for (auto v : flags) {
      if (v.first & _flags) std::cout << v.second << ' ';
    }
    std::cout << "\b>" << std::endl;
  }
  signal(SIGFPE, mpi_gdb_print_trace);
}

//===========================================
// A signal handler to produce stop and wait
//===========================================

void set_fpu_gdb_handler(void)
{
  // Flag all FP errors except inexact
  //
  // fedisableexcept(FE_ALL_EXCEPT & ~FE_INEXACT);

  // Flag invalid FP results only, such as 0/0 or infinity - infinity
  // or sqrt(-1).
  //
  feenableexcept(FE_INVALID);
  //
  // Print enabled flags to root node
  //
  if (myid==0) {
    const std::list<std::pair<int, std::string>> flags =
      {	{FE_DIVBYZERO, "divide-by-zero"},
	{FE_INEXACT,   "inexact"},
	{FE_INVALID,   "invalid"},
	{FE_OVERFLOW,  "overflow"},
	{FE_UNDERFLOW, "underflow"} };
    
    int _flags = fegetexcept();
    std::cout << "Enabled FE flags: <";
    for (auto v : flags) {
      if (v.first & _flags) std::cout << v.second << ' ';
    }
    std::cout << "\b>" << std::endl;
  }
  signal(SIGFPE, mpi_gdb_wait_trace);
}

#endif

                                // Local headers
#include <SphericalSL.H>
#include <DiskHalo2.H>
#include <localmpi.h>


// Hydrogen fraction
//
const double f_H = 0.76;


// Global variables
//
enum DiskType { constant, gaussian, mn, exponential, bmn };

std::map<std::string, DiskType> dtlookup =
  { {"constant",    DiskType::constant},
    {"gaussian",    DiskType::gaussian},
    {"mn",          DiskType::mn},
    {"exponential", DiskType::exponential},
    {"bmn", DiskType::bmn}
  };

DiskType     DTYPE;
double       ASCALE;
double       ASHIFT;
double       HSCALE;
double       FBULGE;
double       BSCALE;
double       RTRUNC = 1.0;
double       RWIDTH = 0.0;

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
  
double DiskDens(double R, double z, double phi)
{
  double ans = 0.0;

  switch (DTYPE) {

  case DiskType::constant:
    if (R < ASCALE && fabs(z) < HSCALE)
      ans = 1.0/(2.0*HSCALE*M_PI*ASCALE*ASCALE);
    break;

  case DiskType::gaussian:
    if (fabs(z) < HSCALE)
      ans = 1.0/(2.0*HSCALE*2.0*M_PI*ASCALE*ASCALE)*
	exp(-R*R/(2.0*ASCALE*ASCALE));
    break;

  case DiskType::mn:
    {
      double Z2 = z*z + HSCALE*HSCALE;
      double Z  = sqrt(Z2);
      double Q2 = (ASCALE + Z)*(ASCALE + Z);
      ans = 0.25*HSCALE*HSCALE/M_PI*(ASCALE*R*R + (ASCALE + 3.0*Z)*Q2)/( pow(R*R + Q2, 2.5) * Z*Z2 );
    }
    break;

  case DiskType::bmn:
    {
      double ra = sqrt(R*R+z*z)/BSCALE;
      double rc = 0; // no core allowed!
      double rho0 = (FBULGE)/(4*M_PI*BSCALE*BSCALE*BSCALE*ra*ra/(2*(1+ra)*(1+ra))); // normalisation for bulge
      double bdens = rho0/((ra+rc)*pow((1+ra),3.)); // bulge density

      double Z2 = z*z + HSCALE*HSCALE;
      double Z  = sqrt(Z2);
      double Q2 = (ASCALE + Z)*(ASCALE + Z);
      double FDISK = 1. - FBULGE;
      double ddens = 0.25*HSCALE*HSCALE*FDISK/M_PI*(ASCALE*R*R + (ASCALE + 3.0*Z)*Q2)/( pow(R*R + Q2, 2.5) * Z*Z2 );
      

      ans = bdens+ddens;
    }
    break;

    
  case DiskType::exponential:
  default:
    {
      double f = cosh(z/HSCALE);
      ans = exp(-R/ASCALE)/(4.0*M_PI*ASCALE*ASCALE*HSCALE*f*f);
    }
    break;
  }

  if (RWIDTH>0.0) ans *= erf((RTRUNC-R)/RWIDTH);

  return ans;
}

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
  // Apply a shift along the x-axis
  //
  double x = R*cos(phiS) - ASHIFT*ASCALE;
  double y = R*sin(phiS);

  return DiskDens(sqrt(x*x + y*y), z, atan2(y, x));
}

int 
main(int ac, char **av)
{
  //====================
  // Inialize MPI stuff
  //====================
  
  local_init_mpi(ac, av);

  //====================
  // Begin opt parsing
  //====================

  int          LMAX;
  int          NMAX;
  int          NUMR;
  int          SCMAP;
  double       RMIN;
  double       RCYLMIN;
  double       RCYLMAX;
  double       SCSPH;
  double       RSPHSL;
  double       DMFAC;
  double       RFACTOR;
  double       ECUT;
  double       X0;
  double       Y0;
  double       Z0;
  double       U0;
  double       V0;
  double       W0;
  int          RNUM;
  int          PNUM;
  int          TNUM;
  int          VFLAG;
  int          DFLAG;
  bool         expcond;
  bool         LOGR;
  bool         CHEBY;
  int          CMAPR;
  int          CMAPZ;
  int          NCHEB;
  int          CMTYPE;
  int          NDR;
  int          NDZ;
  int          NHR;
  int          NHT;
  double       SHFAC;
  int          NMAX2;
  int          LMAX2;
  int          MMAX;
  int          NUMX;
  int          NUMY;
  int          NOUT;
  int          NORDER;
  int          NORDER1;
  int          NODD;
  bool         SELECT;
  bool         DUMPCOEF;
  int          DIVERGE;
  double       DIVERGE_RFAC;
  int          DIVERGE2;
  double       DIVERGE_RFAC2;
  int          DF;
  double       PPower;
  double       R_DF;
  double       DR_DF;
  double       Hratio;
  double       scale_height;
  double       scale_length;
  double       scale_lenfkN;
  double       disk_mass;
  double       gas_mass;
  double       gscal_length;
  double       ToomreQ;
  double       Temp;
  double       Tmin;
  double       gen_ecut;
  bool         const_height;
  bool         images;
  bool         multi;
  bool         SVD;
  int          SEED;
  bool         DENS;
  bool         basis;
  bool         zero;
  bool         report;
  bool         ignore;
  bool         evolved;
  int          nhalo;
  int          ndisk;
  int          ngas;
  int          ngparam;
  string       hbods;
  string       dbods;
  string       gbods;
  string       suffix;
  string       centerfile;
  string       halofile1;
  string       halofile2;
  string       cachefile;
  string       config;
  string       gentype;
  string       dtype;
  string       dmodel;
  string       mtype;
  
  po::options_description desc("Allowed options");
  desc.add_options()
    ("help,h",                                                                          "Print this help message")
    ("conf,c",          po::value<string>(&config),                                     "Write template options file with current and all default values")
    ("input,f",         po::value<string>(&config),                                     "Parameter configuration file")
    ("deproject",       po::value<string>(&dmodel),                                     "The EmpCylSL deprojection from specified disk model (EXP, MN, or BDisc)")
    ("NUMR",            po::value<int>(&NUMR)->default_value(2000),                     "Size of radial grid for Spherical SL")
    ("RMIN",            po::value<double>(&RMIN)->default_value(0.005),                 "Minimum halo radius")
    ("RCYLMIN",         po::value<double>(&RCYLMIN)->default_value(0.001),              "Minimum disk radius")
    ("RCYLMAX",         po::value<double>(&RCYLMAX)->default_value(20.0),               "Maximum disk radius")
    ("SCMAP",           po::value<int>(&SCMAP)->default_value(1),                       "Turn on Spherical SL coordinate mapping (1, 2, 0=off")
    ("SCSPH",           po::value<double>(&SCSPH)->default_value(1.0),                  "Scale for Spherical SL coordinate mapping")
    ("ECUT",            po::value<double>(&ECUT)->default_value(1.0),                   "Energy cutoff for multimass ratio grid")
    ("RSPHSL",          po::value<double>(&RSPHSL)->default_value(47.5),                "Maximum halo expansion radius")
    ("ASCALE",          po::value<double>(&ASCALE)->default_value(1.0),                 "Radial scale length for disk basis construction")
    ("ASHIFT",          po::value<double>(&ASHIFT)->default_value(0.0),                 "Fraction of scale length for shift in conditioning function")
    ("HSCALE",          po::value<double>(&HSCALE)->default_value(0.1),                 "Vertical scale length for disk basis construction")
    ("BSCALE",          po::value<double>(&BSCALE)->default_value(0.1),                 "Bulge scale length")
    ("FBULGE",          po::value<double>(&FBULGE)->default_value(0.0),                 "Fraction of disc mass in bulge")
    ("RTRUNC",          po::value<double>(&RTRUNC)->default_value(0.1),                 "Maximum disk radius for erf truncation of EOF conditioning density")
    ("RWIDTH",          po::value<double>(&RWIDTH)->default_value(0.0),                 "Width for erf truncationofr EOF conditioning density (ignored if zero)")
    ("DMFAC",           po::value<double>(&DMFAC)->default_value(1.0),                  "Disk mass scaling factor for spherical deprojection model")
    ("RFACTOR",         po::value<double>(&RFACTOR)->default_value(1.0),                "Disk radial scaling factor for spherical deprojection model")
    ("X0",              po::value<double>(&X0)->default_value(0.0),                     "Disk-Halo x center position")
    ("Y0",              po::value<double>(&Y0)->default_value(0.0),                     "Disk-Halo y center position")
    ("Z0",              po::value<double>(&Z0)->default_value(0.0),                     "Disk-Halo z center position")
    ("U0",              po::value<double>(&U0)->default_value(0.0),                     "Disk-Halo x velocity center position")
    ("V0",              po::value<double>(&V0)->default_value(0.0),                     "Disk-Halo y velocity center position")
    ("W0",              po::value<double>(&W0)->default_value(0.0),                     "Disk-Halo z velocity center position")
    ("RNUM",            po::value<int>(&RNUM)->default_value(200),                      "Number of radial knots for EmpCylSL basis construction quadrature")
    ("PNUM",            po::value<int>(&PNUM)->default_value(80),                       "Number of azimthal knots for EmpCylSL basis construction quadrature")
    ("TNUM",            po::value<int>(&TNUM)->default_value(80),                       "Number of cos(theta) knots for EmpCylSL basis construction quadrature")
    ("CMAPR",           po::value<int>(&CMAPR)->default_value(1),                     "Radial coordinate mapping type for cylindrical grid (0=none, 1=rational fct)")
    ("CMAPZ",           po::value<int>(&CMAPZ)->default_value(1),                     "Vertical coordinate mapping type for cylindrical grid (0=none, 1=sech, 2=power in z")
    ("SVD",             po::value<bool>(&SVD)->default_value(false),                    "Use svd for symmetric eigenvalue problesm")
    ("LOGR",            po::value<bool>(&LOGR)->default_value(false),                   "Make a logarithmic coordinate mapping")
    ("CHEBY",           po::value<bool>(&CHEBY)->default_value(false),                  "Use Chebyshev smoothing for epicyclic and asymmetric drift")
    ("NCHEB",           po::value<int>(&NCHEB)->default_value(16),                      "Chebyshev order for smoothing")
    ("NDR",             po::value<int>(&NDR)->default_value(1600),                      "Number of points in DiskHalo radial table for disk")
    ("NDZ",             po::value<int>(&NDZ)->default_value(400),                       "Number of points in DiskHalo vertical table for disk")
    ("NHR",             po::value<int>(&NHR)->default_value(1600),                      "Number of points in DiskHalo radial table for halo")
    ("NHT",             po::value<int>(&NHT)->default_value(200),                       "Number of points in DiskHalo cos(theta) table for halo")
    ("SHFAC",           po::value<double>(&SHFAC)->default_value(16.0),                 "Scale height factor for assigning vertical table size")
    ("LMAX",            po::value<int>(&LMAX)->default_value(6),                        "Number of harmonics for Spherical SL for halo/spheroid")
    ("NMAX",            po::value<int>(&NMAX)->default_value(12),                       "Number of radial basis functions in Spherical SL for halo/spheroid")
    ("NMAX2",           po::value<int>(&NMAX2)->default_value(36),                      "Number of radial basis functions in Spherical SL for determining disk basis")
    ("LMAX2",           po::value<int>(&LMAX2)->default_value(36),                      "Number of harmonics for Spherical SL for determining disk basis")
    ("MMAX",            po::value<int>(&MMAX)->default_value(4),                        "Number of azimuthal harmonics for disk basis")
    ("NUMX",            po::value<int>(&NUMX)->default_value(256),                      "Radial grid size for disk basis table")
    ("NUMY",            po::value<int>(&NUMY)->default_value(128),                      "Vertical grid size for disk basis table")
    ("NORDER",          po::value<int>(&NORDER)->default_value(16),                     "Number of disk basis functions per M-order")
    ("NORDER1",         po::value<int>(&NORDER1)->default_value(1000),                  "Restricts disk basis function to NORDER1<NORDER after basis construction for testing")
    ("NODD",            po::value<int>(&NODD)->default_value(-1),                       "Number of vertically antisymmetric disk basis functions per M-order")
    ("NOUT",            po::value<int>(&NOUT)->default_value(1000),                     "Number of radial basis functions to output for each harmonic order")
    ("SELECT",          po::value<bool>(&SELECT)->default_value(false),                 "Enable significance selection in coefficient computation")
    ("DUMPCOEF",        po::value<bool>(&DUMPCOEF)->default_value(false),               "Dump coefficients")
    ("DIVERGE",         po::value<int>(&DIVERGE)->default_value(0),                     "Cusp extrapolation for primary halo model")
    ("DIVERGE_RFAC",    po::value<double>(&DIVERGE_RFAC)->default_value(1.0),           "Extrapolation exponent for primary mass model")
    ("DIVERGE2",        po::value<int>(&DIVERGE2)->default_value(0),                    "Cusp extrapolation for number model")
    ("DIVERGE_RFAC2",   po::value<double>(&DIVERGE_RFAC2)->default_value(1.0),          "Extrapolation exponent for number model")
    ("DF",              po::value<int>(&DF)->default_value(0),                          "Use change-over from Jeans to Eddington")
    ("R_DF",            po::value<double>(&R_DF)->default_value(20.0),                  "Change over radius for Eddington")
    ("DR_DF",           po::value<double>(&DR_DF)->default_value(5.0),                  "Width of change for to Eddington")
    ("scale_height",    po::value<double>(&scale_height)->default_value(0.1),           "Scale height for disk realization")
    ("scale_length",    po::value<double>(&scale_length)->default_value(2.0),           "Scale length for disk realization")
    ("scale_lenfkN",    po::value<double>(&scale_lenfkN)->default_value(-1.0),          "Scale for multimass gas")
    ("disk_mass",       po::value<double>(&disk_mass)->default_value(1.0),              "Mass of stellar adisk")
    ("gas_mass",        po::value<double>(&gas_mass)->default_value(1.0),               "Mass of gaseous disk")
    ("gscal_length",    po::value<double>(&gscal_length)->default_value(4.0),           "Gas disk scale length")
    ("ToomreQ",         po::value<double>(&ToomreQ)->default_value(1.2),                "Toomre Q parameter for stellar disk generation")
    ("Temp",            po::value<double>(&Temp)->default_value(2000.0),                "Gas temperature (in K)")
    ("Tmin",            po::value<double>(&Tmin)->default_value(500.0),                 "Temperature floor (in K) for gas disk generation")
    ("PPOW",            po::value<double>(&PPower)->default_value(5.0),                 "Power exponent in spherical model for deprojection")
    ("const_height",    po::value<bool>(&const_height)->default_value(true),            "Use constant disk scale height")
    ("images",          po::value<bool>(&images)->default_value(false),                 "Print out reconstructed disk profiles")
    ("multi",           po::value<bool>(&multi)->default_value(false),                  "Use multimass halo")
    ("SEED",            po::value<int>(&SEED)->default_value(11),                       "Random number seed")
    ("DENS",            po::value<bool>(&DENS)->default_value(true),                    "Compute the density basis functions")
    ("basis",           po::value<bool>(&basis)->default_value(false),                  "Print out disk and halo basis")
    ("zero",            po::value<bool>(&zero)->default_value(false),                   "zero center of mass and velocity")
    ("nhalo",           po::value<int>(&nhalo)->default_value(1000),                    "Number of halo particles")
    ("ndisk",           po::value<int>(&ndisk)->default_value(1000),                    "Number of disk particles")
    ("ngas",            po::value<int>(&ngas)->default_value(1000),                     "Number of gas particles")
    ("ngparam",         po::value<int>(&ngparam)->default_value(3),                     "Number of gas particle parameters")
    ("hbods",           po::value<string>(&hbods)->default_value("halo.bods"),          "Halo particle output file")
    ("dbods",           po::value<string>(&dbods)->default_value("disk.bods"),          "Disk particle output file")
    ("gbods",           po::value<string>(&gbods)->default_value("gas.bods"),           "Gas particle output file")
    ("suffix",          po::value<string>(&suffix)->default_value(""),                  "Suffix appended for body files")
    ("VFLAG",           po::value<int>(&VFLAG)->default_value(0),                       "Output flags for EmpCylSL")
    ("DFLAG",           po::value<int>(&DFLAG)->default_value(0),                       "Output flags for DiskHalo")
    ("threads",         po::value<int>(&nthrds)->default_value(1),                      "Number of lightweight threads")
    ("expcond",         po::value<bool>(&expcond)->default_value(true),                 "Use analytic density function for computing EmpCylSL basis")
    ("centerfile",      po::value<string>(&centerfile)->default_value("center.dat"),    "Read position and velocity center from this file")
    ("halofile1",       po::value<string>(&halofile1)->default_value("SLGridSph.model"),        "File with input halo model")
    ("halofile2",       po::value<string>(&halofile2)->default_value("SLGridSph.model.fake"),   "File with input halo model for multimass")
    ("cachefile",       po::value<string>(&cachefile)->default_value(".eof.cache.file"),        "Name of EOF cache file")
    ("runtag",          po::value<string>(&runtag)->default_value("run000"),                    "Label prefix for diagnostic images")
    ("gentype",         po::value<string>(&gentype)->default_value("Asymmetric"),               "DiskGenType string for velocity initialization (Jeans, Asymmetric, or Epicyclic)")
    ("mtype",           po::value<string>(&mtype),                                              "Spherical deprojection model for EmpCylSL (one of: Exponential, Gaussian, Plummer, Power)")
    ("condition",       po::value<string>(&dtype)->default_value("exponential"),                "Disk type for condition (one of: constant, gaussian, mn, exponential, bmn)")
    ("report",          po::value<bool>(&report)->default_value(true),                  "Report particle progress in EOF computation")
    ("evolved",         po::value<bool>(&evolved)->default_value(false),                "Use existing halo body file given by <hbods> and do not create a new halo")
    ("ignore",          po::value<bool>(&ignore)->default_value(false),                 "Ignore any existing cache file and recompute the EOF")
    ("newcache",                                                                        "Use new YAML header version for EOF cache file")
    ("ortho",                                                                           "Perform orthogonality check for basis")
    ;
        
  po::variables_map vm;
  
  // Parse command line for control and critical parameters
  //
  try {
    po::store(po::parse_command_line(ac, av, desc), vm);
    po::notify(vm);    
  } catch (po::error& e) {
    if (myid==0) std::cout << "Option error on command line: "
			   << e.what() << std::endl;
    MPI_Finalize();
    return -1;
  }
  
  // Print help message and exit
  //
  if (vm.count("help")) {
    if (myid == 0) {
      const char *mesg = "Generates a Monte Carlo realization of a halo\nwith an embedded disk using Jeans' equations.";
      std::cout << mesg << std::endl
		<< desc << std::endl << std::endl
		<< "Examples: " << std::endl
		<< "\t" << "Use parameters read from a config file in INI style"  << std::endl
		<< "\t" << av[0] << " --input=gendisk.config"  << std::endl << std::endl
		<< "\t" << "Generate a template config file in INI style from current defaults"  << std::endl
		<< "\t" << av[0] << " --conf=template.config" << std::endl << std::endl
		<< "\t" << "Override a single parameter in a config file from the command line"  << std::endl
		<< "\t" << av[0] << "--LMAX=8 --conf=template.config" << std::endl << std::endl;
    }
    MPI_Finalize();
    return 0;
  }

  // Write template config file in INI style and exit
  //
  if (vm.count("conf")) {
    // Do not overwrite existing config file
    //
    if (boost::filesystem::exists(config)) {
      if (myid == 0)
	std::cerr << av[0] << ": config file <" << config
		  << "> exists, will not overwrite" << std::endl;
      MPI_Finalize();
      return 0;
    }

    NOUT = std::min<int>(NOUT, NORDER);

    // Write template file
    //
    if (myid==0) {
      std::ofstream out(config);

      if (out) {
	// Iterate map and print out key--value pairs and description
	//
	for (const auto& it : vm) {
				// Don't write this parameter
	  if (it.first.find("conf")==0) continue;

	  out << std::setw(20) << std::left << it.first << " = ";
	  auto& value = it.second.value();
	  if (auto v = boost::any_cast<uint32_t>(&value))
	    out << std::setw(32) << std::left << *v;
	  else if (auto v = boost::any_cast<int>(&value))
	    out << std::setw(32) << std::left << *v;
	  else if (auto v = boost::any_cast<unsigned>(&value))
	    out << std::setw(32) << std::left << *v;
	  else if (auto v = boost::any_cast<float>(&value))
	    out << std::setw(32) << std::left << *v;
	  else if (auto v = boost::any_cast<double>(&value))
	    out << std::setw(32) << std::left << *v;
	  else if (auto v = boost::any_cast<bool>(&value))
	    out << std::setw(32) << std::left << std::boolalpha << *v;
	  else if (auto v = boost::any_cast<std::string>(&value))
	    out << std::setw(32) << std::left << *v;
	  else
	    out << "error";


	  //                               NO approximations -----+
	  // Add description as a comment                         |
	  //                                                      V
	  const po::option_description& rec = desc.find(it.first, false);
	  out << " # " << rec.description() << std::endl;
	}
      } else {
	if (myid==0)
	  std::cerr << av[0] << ": error opening template config file <"
		    << config << ">" << std::endl;
      }
    }
    MPI_Finalize();
    return 0;
  }

  // Read parameters fron the config file
  //
  if (vm.count("input")) {
    try {
      std::ifstream in(config);
      po::store(po::parse_config_file(in, desc), vm);
      po::notify(vm);    
    } catch (po::error& e) {
      if (myid==0) std::cout << "Option error in configuration file: "
			     << e.what() << std::endl;
      MPI_Finalize();
      return 0;
    }
  }
  
  // Enable new YAML cache header
  //
  if (vm.count("newcache")) {
    EmpCylSL::NewCache = true;
  }

  // Set EmpCylSL mtype.  This is the spherical function used to
  // generate the EOF basis.  If "deproject" is set, this will be
  // overriden in EmpCylSL.
  //
  EmpCylSL::mtype = EmpCylSL::Exponential;
  if (vm.count("mtype")) {
    if (mtype.compare("Exponential")==0)
      EmpCylSL::mtype = EmpCylSL::Exponential;
    else if (mtype.compare("Gaussian")==0)
      EmpCylSL::mtype = EmpCylSL::Gaussian;
    else if (mtype.compare("Plummer")==0)
      EmpCylSL::mtype = EmpCylSL::Plummer;
    else if (mtype.compare("Power")==0) {
      EmpCylSL::mtype = EmpCylSL::Power;
      EmpCylSL::PPOW  = PPower;
    } else {
      if (myid==0) std::cout << "No EmpCylSL EmpModel named <"
			     << mtype << ">, valid types are: "
			     << "Exponential, Gaussian, Plummer" << std::endl;
      MPI_Finalize();
      return 0;
    }
  }

  // Set DiskType.  This is the functional form for the disk used to
  // condition the basis.
  //
				// Convert dtype string to lower case
  std::transform(dtype.begin(), dtype.end(), dtype.begin(),
		 [](unsigned char c){ return std::tolower(c); });

  try {				// Check for map entry, will through if the 
    DTYPE = dtlookup.at(dtype);	// key is not in the map.

    if (myid==0)		// Report DiskType
      std::cout << "DiskType is <" << dtype << ">" << std::endl;
  }
  catch (const std::out_of_range& err) {
    if (myid==0) {
      std::cout << "DiskType error in configuration file" << std::endl;
      std::cout << "Valid options are: ";
      for (auto v : dtlookup) std::cout << v.first << " ";
      std::cout << std::endl;
    }
    MPI_Finalize();
    return 0;
  }

  //====================
  // Okay, now begin ...
  //====================

#ifdef DEBUG                    // For gdb . . . 
  sleep(20);
  // set_fpu_handler();         // Make gdb trap FPU exceptions
  set_fpu_gdb_handler();	// Make gdb trap FPU exceptions
#endif
  
  int n_particlesH, n_particlesD, n_particlesG;
  
  if (suffix.size()>0) {
    hbods = hbods + "." + suffix;
    dbods = dbods + "." + suffix;
    gbods = gbods + "." + suffix;
  }
  
  // Divvy up the particles by core
  //
  n_particlesH = nhalo/numprocs;
  if (myid==0) n_particlesH = nhalo - n_particlesH*(numprocs-1);
  
  n_particlesD = ndisk/numprocs;
  if (myid==0) n_particlesD = ndisk - n_particlesD*(numprocs-1);
  
  n_particlesG = ngas/numprocs;
  if (myid==0) n_particlesG = ngas  - n_particlesG*(numprocs-1);
  
  
#ifdef DEBUG  
  std::cout << "Processor " << myid << ": n_particlesH=" << n_particlesH
	    << std::endl
	    << "Processor " << myid << ": n_particlesD=" << n_particlesD
	    << std::endl
	    << "Processor " << myid << ": n_particlesG=" << n_particlesG
	    << std::endl;
#endif
  
  if (n_particlesH + n_particlesD + n_particlesG <= 0) {
    if (myid==0) std::cout << "You have specified zero particles!" << std::endl;
    MPI_Abort(MPI_COMM_WORLD, 3);
    exit(0);
  }
  
  // Vectors to contain phase space Particle structure is defined in
  // Particle.H
  //
  vector<Particle> dparticles, hparticles;
  
  //
  // Disk halo grid parameters
  //
  //      Prevent model evaulation inside of either grid---------+
  //                                       |                     |
  //                                       V                     V
  DiskHalo::RDMIN       = std::max<double>(RCYLMIN*scale_length, RMIN);
  DiskHalo::RHMIN       = RMIN;
  DiskHalo::RHMAX       = RSPHSL;
  DiskHalo::RDMAX       = RCYLMAX*scale_length;
  DiskHalo::NDR         = NDR;
  DiskHalo::NDZ         = NDZ;
  DiskHalo::NHR         = NHR;
  DiskHalo::NHT         = NHT;
  DiskHalo::SHFACTOR    = SHFAC;
  DiskHalo::COMPRESSION = DMFAC;
  DiskHalo::LOGSCALE    = 1;
  DiskHalo::NUMDF       = 4000;
  DiskHalo::Q           = ToomreQ;
  DiskHalo::R_DF        = R_DF;
  DiskHalo::DR_DF       = DR_DF;
  DiskHalo::ECUT_DF     = ECUT;
  DiskHalo::SEED        = SEED;
  DiskHalo::VFLAG       = static_cast<unsigned int>(DFLAG);
  DiskHalo::CHEBY       = CHEBY;
  DiskHalo::NCHEB       = NCHEB;
  if (suffix.size())
    DiskHalo::RUNTAG    = suffix;
  
  AddDisk::use_mpi      = true;
  AddDisk::Rmin         = RMIN;
  
  //===========================Spherical expansion=============================
  
  // SLGridSph::diverge = DIVERGE;
  // SLGridSph::divergexp = DIVERGE_RFAC;
  SLGridSph::model_file_name = halofile1;
  
  SphericalSL::RMIN = RMIN;
  SphericalSL::RMAX = RSPHSL;
  SphericalSL::NUMR = NUMR;
  // Create expansion only if needed . . .
  boost::shared_ptr<SphericalSL> expandh;
  if (n_particlesH) {
    expandh = boost::make_shared<SphericalSL>(nthrds, LMAX, NMAX, SCMAP, SCSPH);
  }

  //===========================Cylindrical expansion===========================


  EmpCylSL::RMIN        = RCYLMIN;
  EmpCylSL::RMAX        = RCYLMAX;
  EmpCylSL::NUMX        = NUMX;
  EmpCylSL::NUMY        = NUMY;
  EmpCylSL::NUMR        = NUMR;
  EmpCylSL::NOUT        = NOUT;
  EmpCylSL::CMAPR       = CMAPR;
  EmpCylSL::CMAPZ       = CMAPZ;
  EmpCylSL::VFLAG       = VFLAG;
  EmpCylSL::logarithmic = LOGR;
  EmpCylSL::DENS        = DENS;
  EmpCylSL::USESVD      = SVD;
  EmpCylSL::PCAVAR      = SELECT;
  EmpCylSL::CACHEFILE   = cachefile;

  if (basis) EmpCylSL::DENS = true;

                                // Create expansion only if needed . . .
  boost::shared_ptr<EmpCylSL> expandd;
  bool save_eof = false;

  if (n_particlesD) {

    expandd = boost::make_shared<EmpCylSL>(NMAX2, LMAX2, MMAX, NORDER, ASCALE, HSCALE, NODD);

#ifdef DEBUG
   std::cout << "Process " << myid << ": "
	     << " rmin="   << EmpCylSL::RMIN
	     << " rmax="   << EmpCylSL::RMAX
	     << " a="      << ASCALE
	     << " h="      << HSCALE
	     << " nmax2="  << NMAX2
	     << " lmax2="  << LMAX2
	     << " mmax="   << MMAX
	     << " nordz="  << NORDER
	     << " noddz="  << NODD
	     << std::endl  << std::flush;
#endif

    // Try to read existing cache to get EOF
    //
    if (not ignore) {
      save_eof = expandd->read_cache();
      if (myid==0) {
	if (save_eof) {
	  std::cout << "EmpCylSL requested cache read: GOOD, continuing" << std::endl;
	} else {
	  std::cout << "EmpCylSL requested cache read FAIL: exiting ..." << std::endl;
	}
      }
      if (not save_eof) {
	MPI_Finalize();
	return 0;
      }
    }

    // Use these user models to deproject for the EOF spherical basis
    //
    if (vm.count("deproject")) {
      // The scale in EmpCylSL is assumed to be 1 so we compute the
      // height relative to the length
      //
      double H = scale_height/scale_length;
      double Ab = scale_length/6.; // hard-coded for testing

      // The model instance (you can add others in DiskModels.H).
      // It's MN or Exponential if not MN.
      //
      EmpCylSL::AxiDiskPtr model;

      if (dmodel.compare("MN")==0) // Miyamoto-Nagai
	model = boost::make_shared<MNdisk>(1.0, H);
      else if (dmodel.compare("BDisc")==0) // Bulge+Miyamoto-Nagai
	model = boost::make_shared<BDisc>(Ab, 1.0, H, FBULGE);
      else			// Default to exponential
	model = boost::make_shared<Exponential>(1.0, H);

      if (RWIDTH>0.0) {
	model = boost::make_shared<Truncated>(RTRUNC/scale_length,
					      RWIDTH/scale_length,
					      model);
	if (myid==0)
	  std::cout << "Made truncated model with R=" << RTRUNC/scale_length
		    << " and W=" << RWIDTH/scale_length << std::endl;
      }

      expandd->create_deprojection(H, RFACTOR, NUMR, RNUM, model);
    }

    // Regenerate EOF from analytic density
    //
    if (expcond and not save_eof) {
      expandd->generate_eof(RNUM, PNUM, TNUM, dcond);
      save_eof = true;
    }

    // Basis orthgonality check
    //
    if (vm.count("ortho")) {
      std::ofstream out(runtag + ".ortho_check");
      expandd->ortho_check(out);
    }
  }

  //====================Create the disk & halo model===========================

  boost::shared_ptr<DiskHalo> diskhalo;

  if (multi) {
    if (myid==0) std::cout << "Initializing for a MULTI-MASS halo . . . " << std::flush;
    diskhalo =
      boost::make_shared<DiskHalo>
      (expandh, expandd,
       scale_height, scale_length, disk_mass, 
       halofile1, DIVERGE,  DIVERGE_RFAC,
       halofile2, DIVERGE2, DIVERGE_RFAC2,
       DiskHalo::getDiskGenType[gentype]);
    if (myid==0) std::cout << "done" << std::endl;

  } else {

    if (myid==0) std::cout << "Initializing for a SINGLE-MASS halo . . . " << std::flush;
    diskhalo = boost::make_shared<DiskHalo>
      (expandh, expandd,
       scale_height, scale_length, 
       disk_mass, halofile1,
       DF, DIVERGE, DIVERGE_RFAC,
       DiskHalo::getDiskGenType[gentype]);
    if (myid==0) std::cout << "done" << std::endl;
  }
  
  std::ifstream center(centerfile.c_str());
  if (center) {

    bool ok = true;

    center >> X0;
    if (center.fail()) ok = false;

    center >> Y0;
    if (center.fail()) ok = false;

    center >> Z0;
    if (center.fail()) ok = false;

    if (ok) {
      diskhalo->set_pos_origin(X0, Y0, Z0);
      if (myid==0) std::cout << "Using position origin: " 
			     << X0 << ", " << Y0 << ", " << Z0 << std::endl;
    }

    center >> U0;
    if (center.fail()) ok = false;

    center >> V0;
    if (center.fail()) ok = false;

    center >> W0;
    if (center.fail()) ok = false;

    if (ok) {
      diskhalo->set_vel_origin(U0, V0, W0);
      if (myid==0) std::cout << "Using velocity origin: " 
			     << U0 << ", " << V0 << ", " << W0 << std::endl;
    }
  }

                                // Make zero center of mass and
                                // center of velocity
  diskhalo->zero_com(zero);
  diskhalo->zero_cov(zero);
  
  //===========================================================================

                                // Open output file (make sure it exists
                                // before realizing a large phase space)
  std::ofstream out_halo, out_disk;
  if (myid==0) {
    if (not evolved and n_particlesH) {
      out_halo.open(hbods.c_str());
      if (!out_halo) {
	cout << "Could not open <" << hbods << "> for output\n";
	MPI_Abort(MPI_COMM_WORLD, 4);
	exit(0);
      }
    }

    if (n_particlesD) {
      out_disk.open(dbods.c_str());
      if (!out_disk) {
	std::cout << "Could not open <" << dbods << "> for output" << std::endl;
	MPI_Abort(MPI_COMM_WORLD, 4);
	exit(0);
      }
    }
  }

  //=================Make the phase space coordinates==========================

  if (evolved) {		// ---------------------------
				// Use existing halo body file
    std::ifstream hin(hbods);	// ---------------------------
    
    if (hin) {
      std::string line;
      std::getline(hin, line);
      std::istringstream sin(line);

      int niatr, ndatr;
      sin >> nhalo;
      sin >> niatr;
      sin >> ndatr;
      
      // Divvy up the particles by core.  The root node gets any
      // remainder.
      //
      n_particlesH = nhalo/numprocs;

      int ibeg = 0;
      int iend = nhalo - n_particlesH*(numprocs-myid-1);
      
      if (myid>0) {
	ibeg = nhalo - n_particlesH*(numprocs-myid);
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
	hparticles.push_back(P);
      }

    } else {
      // Error message
      if (myid==0)
	std::cout << "Could not read halo file <" << hbods
		  << "> . . . quitting" << std::endl;

      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Finalize();
      exit(-1);
    }

    std::cout << "Process " << myid << " has " << hparticles.size()
	      << " halo particles" << std::endl;

    if (myid==0)
      std::cout << "Generating halo tables for input halo . . . "
		<< std::flush;

    if (multi) {
      diskhalo->set_halo_table_multi(hparticles);
    } else {
      diskhalo->set_halo_table_single(hparticles);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    if (myid==0) std::cout << "done" << std::endl;

  } else {			// ---------------------------
				// Generate new halo body file
    if (n_particlesH) {		// ---------------------------
      if (multi) {
	if (myid==0) std::cout << "Generating halo phase space . . . " << std::flush;
	diskhalo->set_halo(hparticles, nhalo, n_particlesH);
      } else {
	if (myid==0) std::cout << "Generating halo coordinates . . . " << std::flush;
	diskhalo->set_halo_coordinates(hparticles, nhalo, n_particlesH);
	MPI_Barrier(MPI_COMM_WORLD);
      }
      MPI_Barrier(MPI_COMM_WORLD);
      if (myid==0) std::cout << "done" << std::endl;
    }
  }

  if (n_particlesD) {
    if (myid==0) std::cout << "Generating disk coordinates . . . " << std::flush;
    diskhalo->set_disk_coordinates(dparticles, ndisk, n_particlesD);
    MPI_Barrier(MPI_COMM_WORLD);
    if (myid==0) std::cout << "done" << std::endl;
  }

  if (n_particlesH) {
    if (myid==0) std::cout << "Beginning halo accumulation . . . " << std::flush;
    expandh->accumulate(hparticles);
    MPI_Barrier(MPI_COMM_WORLD);

    if (myid==0) {
      std::cout << "done" << std::endl;
      if (DUMPCOEF) {
	std::cout << "Dumping coefficients halo . . . " << std::flush;
	ostringstream sout;
	sout << "halo_coefs.";
	if (suffix.size()>0) sout << suffix;
	else                 sout << "dump";
	ofstream out(sout.str().c_str());
	if (out) expandh->dump_coefs(out, false);
	std::cout << "done" << std::endl;
      }
    }

  }
  
  if (n_particlesD) {
    if (myid==0) std::cout << "Beginning disk accumulation . . . " << std::flush;
    expandd->setup_accumulation();

    if (!save_eof and !expcond) {
      expandd->setup_eof();
      if (nthrds>1)
	expandd->accumulate_eof_thread(dparticles, report);
      else
	expandd->accumulate_eof(dparticles, report);
      MPI_Barrier(MPI_COMM_WORLD);

      if (myid==0) std::cout << "done" << std::endl;
  
      if (myid==0) std::cout << "Making the EOF . . . " << std::flush;
      expandd->make_eof();
      MPI_Barrier(MPI_COMM_WORLD);
    }

    if (myid==0) {
      std::cout << "done" << std::endl;
      std::cout << "Making disk coefficients . . . " << std::flush;
    }

    expandd->make_coefficients();
    MPI_Barrier(MPI_COMM_WORLD);

    if (myid==0) {
      std::cout << "done" << std::endl;
      std::cout << "Reexpand . . . " << std::flush;
    }

    if (nthrds>1)
      expandd->accumulate_thread(dparticles, 0, report);
    else
      expandd->accumulate(dparticles, 0, report);

    expandd->make_coefficients();
    MPI_Barrier(MPI_COMM_WORLD);
    if (myid==0) {
      std::cout << "done" << std::endl;
      if (DUMPCOEF) {
	std::cout << "Dumping coefficients . . . " << std::flush;
	ostringstream sout;
	sout << "disk_coefs.";
	if (suffix.size()>0) sout << suffix;
	else                 sout << "dump";
	ofstream out(sout.str().c_str());
	if (out) expandd->dump_coefs(out);
	std::cout << "done" << std::endl;
      }
    }

    if (NORDER1<NORDER) {
      if (myid==0) std::cout << "Restricting order from " << NORDER 
			     << " to " << NORDER1 << " . . . " << std::flush;
      expandd->restrict_order(NORDER1);
      if (myid==0) std::cout << "done" << std::endl;
    }

    if (images && myid==0) {
      std::cout << "Images . . . " << std::flush;
      std::ostringstream dumpname;
      dumpname << "images.0";
      expandd->dump_images(dumpname.str(), 5.0*ASCALE, 5.0*HSCALE, 64, 64, true);
      std::cout << "done" << std::endl;
    }
  }
  

  //===========================Diagnostics=====================================

                                // For examining the coverage, etc.
                                // Images can be contoured in SM using
  if (myid==0 && basis) {	// the "ch" file type
    
    std::cout << "Dumping basis images . . . " << std::flush;
    
    if (n_particlesD) {
      int nout = 200;
      string dumpstr = runtag + ".dump";
      expandd->dump_basis(dumpstr, 0);
      expandd->dump_images(runtag, 5.0*scale_length, 5.0*scale_height,
			   nout, nout, false);
      expandd->dump_images_basis(runtag, 5.0*scale_length, 5.0*scale_height,
				 nout, nout, false, 0, MMAX, 0, NORDER-1);
    }


    if (n_particlesH) {
      string extn("test");
      expandh->dump_basis(extn);
    }
    
    if (n_particlesH) {
      
      const int nstr = 5;
      const char *names[nstr] = {".dens", ".potl", ".potr", ".pott", ".potp"};
      ofstream *out = new ofstream [nstr];
      
      int nout = 200;
      double rmax = 6.0*scale_length;
      double x, y, dr = 2.0*rmax/(nout-1);
      float f;
    
      for (int i=0; i<nstr; i++) {
        string name("halo");
        name += names[i];
        out[i].open(name.c_str());
        
        out[i].write((char *)&nout, sizeof(int));
        out[i].write((char *)&nout, sizeof(int));
        out[i].write((char *)&(f=-rmax), sizeof(float));
        out[i].write((char *)&(f= rmax), sizeof(float));
        out[i].write((char *)&(f=-rmax), sizeof(float));
        out[i].write((char *)&(f= rmax), sizeof(float));
      }
      
      double r, theta, phi;
      double dens, potl, potr, pott, potp;
    
      for (int j=0; j<nout; j++) {
        y = -rmax + dr*j;
      
        for (int i=0; i<nout; i++) {
          x = -rmax + dr*i;
        
          r = sqrt(x*x + y*y);
          theta = 0.5*M_PI;
          phi = atan2(y, x);
        
          expandh->determine_fields_at_point(r, theta, phi,
                                             &dens, &potl, 
                                             &potr, &pott, &potp);
        
          out[0].write((char *)&(f=dens), sizeof(float));
          out[1].write((char *)&(f=potl), sizeof(float));
          out[2].write((char *)&(f=potr), sizeof(float));
          out[3].write((char *)&(f=pott), sizeof(float));
          out[4].write((char *)&(f=potp), sizeof(float));
        }
        
      }
    
      for (int i=0; i<nstr; i++) out[i].close();
      delete [] out;
    }

    if (n_particlesD) {

      const int nstr = 5;
      const char *names[nstr] = {".dens", ".pot", ".fr", ".fz", ".fp"};
      ofstream *out = new ofstream [nstr];
    
      int nout = 200;
      double rmax = DiskHalo::RDMAX;
      double x, y, dr = 2.0*rmax/(nout-1);
      float f;
    
      for (int i=0; i<nstr; i++) {
        string name("disk");
        name += names[i];
        out[i].open(name.c_str());
        
        out[i].write((char *)&nout, sizeof(int));
        out[i].write((char *)&nout, sizeof(int));
        out[i].write((char *)&(f=-rmax), sizeof(float));
        out[i].write((char *)&(f= rmax), sizeof(float));
        out[i].write((char *)&(f=-rmax), sizeof(float));
        out[i].write((char *)&(f= rmax), sizeof(float));
      }
    
      double z = 0.0, d0, p0, d, p, fr, fz, fp;
    
      for (int j=0; j<nout; j++) {
        y = -rmax + dr*j;
      
        for (int i=0; i<nout; i++) {
          x = -rmax + dr*i;
        
	  if (x<0.0)
	    expandd->accumulated_eval(fabs(x), y, M_PI, p0, p, fr, fz, fp);
	  else
	    expandd->accumulated_eval(x, y,  0.0, p0, p, fr, fz, fp);


          d = expandd->accumulated_dens_eval(sqrt(x*x + y*y), z, atan2(y, x), d0);
        
          
          out[0].write((char *)&(f=d ), sizeof(float));
          out[1].write((char *)&(f=p ), sizeof(float));
          out[2].write((char *)&(f=fr), sizeof(float));
          out[3].write((char *)&(f=fz), sizeof(float));
          out[4].write((char *)&(f=fp), sizeof(float));
        }
        
      }
    
      for (int i=0; i<nstr; i++) out[i].close();
      delete [] out;
    }
    
    std::cout << "done" << std::endl;
  }
  
  MPI_Barrier(MPI_COMM_WORLD);

  //====================Make the phase space velocities========================

  if (!multi and !evolved) {
    if (myid==0) std::cout << "Generating halo velocities . . . " << std::flush;
    diskhalo->set_vel_halo(hparticles);
    if (myid==0) std::cout << "done" << std::endl;
  }
  
  if (myid==0) std::cout << "Generating disk velocities . . . " << std::flush;
  diskhalo->set_vel_disk(dparticles);
  if (myid==0) std::cout << "done" << std::endl;
  

  //====================All done: write it out=================================

  if (not evolved) {
    if (myid==0) std::cout << "Writing phase space file for halo . . . " << std::flush;
    diskhalo->write_file(out_halo, hparticles);
    if (myid==0) std::cout << "done" << std::endl;
    out_halo.close();
  }

  if (myid==0) std::cout << "Writing phase space file for disk . . . " << std::flush;
  diskhalo->write_file(out_disk, dparticles);
  if (myid==0) std::cout << "done" << std::endl;
  out_disk.close();
                                // Diagnostic . . .
  diskhalo->virial_ratio(hparticles, dparticles);

  std::ofstream outprof("profile.diag");
  diskhalo->profile(outprof, dparticles, 3.0e-3*ASCALE, 5.0*ASCALE, 100);

  //====================Compute gas particles==================================

  if (myid==0 && n_particlesG) {
    std::cout << "Computing gas particles . . . " << std::endl;

    // UNITS
    // -------------------
				// cm
    const double pc = 3.08568025e18;
				// proton mass
    const double m_p = 1.67262158e-24;
				// g
    const double msun = 1.98892e33; //
				// cgs
    const double G = 6.67300e-08;
				// cgs
    const double boltz = 1.3806503e-16;

    double T = Temp;

    
    double Lunit = 3.0e5*pc;	// Virial radius
    double Munit = 1.0e12*msun;	// Virial mass
    double Tunit = sqrt(Lunit*Lunit*Lunit/(Munit*G));
    double Vunit = Lunit/Tunit;

    // Fac = kT*R_vir/(G*m_p*M_vir)
    // where M_vir = 1e12 Msun, R_vir=300 kpc
    //
    double fac = T/(G*m_p*Munit/(Lunit*boltz));

    // Thermal velocity in system units
    //
    double mm   = f_H*m_p + (1.0-f_H)*4.0*m_p;
    double vthermal = sqrt( (boltz*T)/mm ) / Vunit;
    double vmin2 = (boltz*Tmin/mm) / (Vunit*Vunit);

    // Adjust scale for multimass gas
    //
    double Scale_Length = gscal_length;
    if (scale_lenfkN > 0.0) gscal_length = scale_lenfkN;

    // Compute using Jeans theorem
    //
    double rmin = RMIN;
    double rmax = 10.0*gscal_length;
    double zmin = 0.001*scale_height;
    int   nrint = 200;
    int   nzint = 400;
    
    vector< vector<double> > zrho, zmas, vcir;
    double r, R, dR = (rmax - rmin)/(nrint-1);
    double z, dz = (log(rmax) - log(zmin))/(nzint-1);

    double p0, p, fr, fz, fp, dens, potl, potr, pott, potp;

    std::cout << "Const_height=" << (const_height ? "True" : "False") << std::endl;

    if (const_height) {

      for (int i=0; i<nrint; i++) {
	R = rmin + dR*i;

	vector<double> lrho(nzint), trho(nzint), tcir(nzint), tmas(nzint, 0);

	for (int j=0; j<nzint; j++) {
	  z = zmin*exp(dz*j);
	  r = sqrt(R*R + z*z);
	  
	  double pot=0.0, frt0=0.0, fzt0=0.0;
	  if (expandd) {
	    expandd->accumulated_eval(R, z, 0, p0, p, fr, fz, fp);
	    frt0 += -fr;
	    fzt0 += -fz;
	    pot += p0;
	  }
	  if (expandh) {
	    expandh->determine_fields_at_point(r, acos(z/(r+1.0e-8)), 0.0,
					       &dens, &potl, 
					       &potr, &pott, &potp);
	    
	    frt0 += potr;
	    fzt0 += (potr*z + pott*R*R/(r*r))/r;
	    pot += potl;
	  }
	  
	  trho[j] = fzt0*scale_height;
	  tcir[j] = sqrt(max<double>(R*frt0-R*trho[j]/Scale_Length, 0.0));
	}
	
	for (int j=0; j<nzint; j++) 
	  tmas[j] = 1.0 - exp(-zmin*exp(dz*j)/scale_height);
	
	zrho.push_back(trho);
	zmas.push_back(tmas);
	vcir.push_back(tcir);
      }

      //
      // Vertical table
      //
      string ztable("ztable.dat");
      std::cout << "Writing " << setw(15) << right << ztable
		<< " [gas] . . . " << std::flush;
      ofstream ztest(ztable.c_str());
      for (int i=0; i<nrint; i++) {
	for (int j=0; j<nzint; j++) {
	  ztest << setw(15) << rmin + dR*i
		<< setw(15) << zmin*exp(dz*j)
		<< setw(15) << zrho[i][j]
		<< setw(15) << zrho[i][j]*Vunit*Vunit*mm/boltz
		<< setw(15) << zmas[i][j]
		<< setw(15) << vcir[i][j]
		<< endl;
	}
	ztest << endl;
      }
      ztest.close();
      std::cout << "done" << std::endl;
      
    } else {

      for (int i=0; i<nrint; i++) {
	R = rmin + dR*i;



	vector<double> lrho(nzint), trho(nzint), tcir(nzint), tmas(nzint, 0);

	for (int j=0; j<nzint; j++) {
	  z = zmin*exp(dz*j);
	  r = sqrt(R*R + z*z);
	  
	  double frt0=0.0, fzt0=0.0;
	  if (expandd) {
	    expandd->accumulated_eval(R, z, 0, p0, p, fr, fz, fp);
	    frt0 += -fr;
	    fzt0 += -fz;
	  }
	  if (expandh) {
	    expandh->determine_fields_at_point(r, acos(z/(r+1.0e-8)), 0.0,
					       &dens, &potl, 
					       &potr, &pott, &potp);
	    frt0 += potr;
	    fzt0 += (potr*z + pott*R*R/(r*r))/r;
	  }
	  
	  trho[j] = -fzt0/(vthermal*vthermal);
	  tcir[j] = sqrt(max<double>(R*frt0-R*vthermal*vthermal/Scale_Length, 0.0));
	}
	
	double mass = 0.0;
	double zfac = 1.0 - exp(-dz);
				    
	lrho[0] = 0.0;
	for (int j=1; j<nzint; j++) 
	  lrho[j] = lrho[j-1] + 0.5*(trho[j-1] + trho[j]) * zmin*exp(dz*j)*zfac;
	
	for (int j=1; j<nzint; j++) 
	  tmas[j] = tmas[j-1] + 0.5*(exp(lrho[j-1]) + exp(lrho[j])) * zmin*exp(dz*j)*zfac;
	
	for (int j=0; j<nzint; j++) {
	  if (tmas[nzint-1]>0.0 && !std::isnan(tmas[nzint-1])) {
	    trho[j]  = exp(lrho[j])/tmas[nzint-1];
	    tmas[j] /= tmas[nzint-1];
	  } else {
	    trho[j] = 0.0;
	    if (j==0) tmas[j] = 0.0;
	    else      tmas[j] = 1.0;
	  }
	}
	zrho.push_back(trho);
	zmas.push_back(tmas);
	vcir.push_back(tcir);
      }


      //
      // Vertical table
      //
      std::cout << "Writing ztable.dat [gas] . . . " << std::flush;
      std::ofstream ztest("ztable.dat");
      for (int i=0; i<nrint; i++) {
	for (int j=0; j<nzint; j++) {
	  ztest << setw(15) << rmin + dR*i
		<< setw(15) << zmin*exp(dz*j)
		<< setw(15) << zrho[i][j]
		<< setw(15) << zmas[i][j]
		<< setw(15) << vcir[i][j]
		<< endl;
	}
	ztest << endl;
      }
      ztest.close();
      std::cout << "done" << std::endl;
      
    }

    // 
    // Prepare output stream
    //
    ofstream outps("gas.bods");
    if (!outps) {
      cerr << "Couldn't open <" << "gas.bods" << "> for output\n";
      exit (-1);
    }

    const int ITMAX=1000;
    const int NREPORT=1000;
    
    //
    // Maximum enclosed disk mass given rmax
    //
    double rmx2 = 1.5*rmax;
    double mmx2 = 1.0 - (1.0 + rmx2/gscal_length)*exp(-rmx2/gscal_length);
    double mmax = 1.0 - (1.0 + rmax/gscal_length)*exp(-rmax/gscal_length);
    double mfac = 1.0 - (1.0 + rmax/Scale_Length)*exp(-rmax/Scale_Length);

    //
    // Random generators
    //
    ACG gen(10, 20);
    Uniform unit(0.0, 1.0, &gen);

    //
    // Trimmed Gaussian
    //
    double minK=0.0, maxK=1.0, sigma = 3.0;
    if (sigma>0) {
      minK = 0.5*(1.0+erf(-0.5*sigma));
      maxK = 0.5*(1.0+erf( 0.5*sigma));
    }
    Uniform unitN(minK, maxK, &gen);


    double gmass, gmass0 = gas_mass/ngas;
    double KE=0.0, VC=0.0;
    vector<double> mc2(nzint);

    gmass = gmass0;
    fr = fz = potr = 0.0;

    outps << setw(8) << ngas
	  << setw(6) << 0 << setw(6) << ngparam << endl;

    for (int n=0; n<ngas; n++) {

      double F, dF, M=mmax*unit(), Z=unit();
      double R = M*rmax, phi=2.0*M_PI*unit(), x, y, z, rr, vc;
      double ax, ay, az;

				// Narrow with bisection
      double rm = 0.0, rp = rmx2;
      double fm = -M, fp = mmx2 - M;
      for (int j=0; j<15; j++) {
	R = 0.5*(rm + rp);
	F = 1.0 - M - (1.0 + R/gscal_length)*exp(-R/gscal_length);
	if (fm*F<0.0) {
	  rp = R;
	  fp = F;
	} else {
	  rm = R;
	  fm = F;
	}
      }
				// Polish with Newton-Raphson
      for (int j=0; j<ITMAX; j++) {
	F = 1.0 - M - (1.0 + R/gscal_length)*exp(-R/gscal_length);
	dF = R/(gscal_length*gscal_length)*exp(-R/gscal_length);
	R += -F/dF;
	if (fabs(F/dF)<1.0e-12) break;
      }
    
      int indr = static_cast<int>(floor(R/dR));
      if (indr<0) indr=0;
      if (indr>nrint-2) indr=nrint-2;
      double a = (dR*(indr+1) - R)/dR;
      double b = (R - indr*dR)/dR;

      vector<double> mz(nzint), vz(nzint);
      for (int j=0; j<nzint; j++) {
	mz[j] = a*zmas[indr][j] + b*zmas[indr+1][j];
	vz[j] = a*vcir[indr][j] + b*vcir[indr+1][j];
      }
      for (int j=0; j<nzint; j++) mz[j] /= mz[nzint-1];
      
      if (const_height) {
	for (int j=0; j<nzint; j++) 
	  mc2[j] = max<double>(a*zrho[indr][j] + b*zrho[indr+1][j], vmin2);
      }

      int indz = max<int>(0, min<int>(nzint-2, Vlocate(Z, mz)));

      a = (mz[indz+1] - Z)/(mz[indz+1] - mz[indz]);
      b = (Z - mz[indz  ])/(mz[indz+1] - mz[indz]);

      vc = fabs(a*vz[indr] + b*vz[indr+1]);

      z = zmin*exp(dz*(a*indz + b*(indz+1)));
      if (unit()<0.5) z *= -1.0;
      rr = sqrt(R*R + z*z);

      if (const_height) {
	vthermal = a*mc2[indz] + b*mc2[indz+1];
	vthermal = sqrt(max<double>(vmin2, vthermal));
      }

      double sinp = sin(phi), cosp = cos(phi);
      x = R*cosp;
      y = R*sinp;

      double u = -vc*sinp + vthermal*norminv(unitN());
      double v =  vc*cosp + vthermal*norminv(unitN());
      double w =  vthermal*norminv(unitN());
      
      gmass = gmass0*exp(-R*(1.0/Scale_Length - 1.0/gscal_length)) * 
	mmax*gscal_length*gscal_length/(mfac*Scale_Length*Scale_Length);

      outps << setw(18) << gmass
	    << setw(18) << R*cos(phi)
	    << setw(18) << R*sin(phi)
	    << setw(18) << z
	    << setw(18) << u
	    << setw(18) << v
	    << setw(18) << w;
      for (int k=0; k<ngparam; k++) outps << setw(18) << 0.0;
      outps << endl;
    
      if (expandd)
	expandd->accumulated_eval(R, z, phi, p0, p, fr, fz, fp);

      if (expandh)
	expandh->determine_fields_at_point(rr, acos(z/(rr+1.0e-8)), 0.0,
					   &dens, &potl, 
					   &potr, &pott, &potp);
      KE += 0.5*gmass*(u*u + v*v + w*w);

      VC += gmass*(-rr*potr + R*fr + z*fz);

      if (!((n+1)%NREPORT)) std::cout << "\r." << n+1 << std::flush;
    }

    std::cout << endl << "Done!" << std::endl;

    std::cout << "****************************" << std::endl
	      << "  Gas disk"                   << std::endl
	      << "----------------------------" << std::endl
	      << "  KE       = " << KE << std::endl
	      << "  VC       = " << VC << std::endl;
    if (VC<0.0)
      std::cout << " -2T/W     = " << -2.0*KE/VC << std::endl;
    std::cout << "****************************"  << std::endl;
  }

  //===========================================================================
  // Shutdown MPI
  //===========================================================================

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();

  return 0;
}
