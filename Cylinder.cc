using namespace std;

#include <sstream>
#include <chrono>
#include <limits>

#include "expand.h"
#include <gaussQ.h>
#include <EmpCylSL.h>
#include <Cylinder.H>
#include <MixtureBasis.H>
#include <Timer.h>

Timer timer_debug;

double EXPSCALE=1.0, HSCALE=1.0, ASHIFT=0.25;

double DiskDens(double R, double z, double phi)
{
  double f = cosh(z/HSCALE);
  return exp(-R/EXPSCALE)/(4.0*M_PI*EXPSCALE*EXPSCALE*HSCALE*f*f);
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
  double x = R*cos(phiS) - ASHIFT*EXPSCALE;
  double y = R*sin(phiS);
  return DiskDens(sqrt(x*x + y*y), z, atan2(y, x));
}


Cylinder::Cylinder(const YAML::Node& conf, MixtureBasis *m) : Basis(conf)
{
#if HAVE_LIBCUDA==1
  if (m) {
    throw std::runtime_error("Error in Cylinder: MixtureBasis logic is not yet implemented in CUDA");
  }

  // Initialize the circular storage container 
  cuda_initialize();

#endif

  id              = "Cylinder";
  geometry        = cylinder;
  mix             = m;
  dof             = 3;
				// Default values

  rcylmin         = 0.001;	// Should only change these two in
  rcylmax         = 20.0;	// extreme circumstances

  ncylnx          = 128;	// These defaults should do fine in
  ncylny          = 128;	// most cases, as well
  ncylr           = 2000;

  acyl            = 1.0;
  nmax            = 20;
  lmax            = 36;
  mmax            = 4;
  hcyl            = 1.0;
  ncylorder       = 10;
  ncylrecomp      = -1;

  rnum            = 100;
  pnum            = 40;
  tnum            = 40;
  ashift          = 0.25;

  vflag           = 0;
  eof             = 1;
  npca            = 50;
  npca0           = 0;
  self_consistent = true;
  firstime        = true;
  expcond         = true;
  cmap            = true;
  logarithmic     = false;
  pcavar          = false;
  pcavtk          = false;
  pcadiag         = false;
  pcaeof          = false;
  nvtk            = 1;
  pcainit         = true;
  density         = false;
  coef_dump       = true;
  try_cache       = true;
  dump_basis      = false;
  compute         = false;
  firstime_coef   = true;
  eof_file        = "";

  initialize();


  EmpCylSL::RMIN        = rcylmin;
  EmpCylSL::RMAX        = rcylmax;
  EmpCylSL::NUMX        = ncylnx;
  EmpCylSL::NUMY        = ncylny;
  EmpCylSL::NUMR        = ncylr;
  EmpCylSL::CMAP        = cmap;
  EmpCylSL::logarithmic = logarithmic;
  EmpCylSL::CACHEFILE   = outdir + ".eof.cache." + runtag;
  EmpCylSL::VFLAG       = vflag;

  // EOF default file name override.  Default uses runtag suffix as
  // above.  Override file must exist if explicitly specified.
  //
  if (eof_file.size()) EmpCylSL::CACHEFILE = eof_file;

  // For debugging; no use by force algorithm
  //
  if (density) EmpCylSL::DENS = true;

  // Make the empirical orthogonal basis instance
  //
  ortho = new EmpCylSL(nmax, lmax, mmax, ncylorder, acyl, hcyl);
  
  try {
    if (conf["tk_type"]) ortho->setTK(conf["tk_type"].as<std::string>());
  }
  catch (YAML::Exception & error) {
    if (myid==0) std::cout << "Error parsing parameters in Cylinder: "
			   << error.what() << std::endl
			   << std::string(60, '-') << std::endl
			   << "Config node"        << std::endl
			   << std::string(60, '-') << std::endl
			   << conf
			   << std::string(60, '-') << std::endl;
    MPI_Finalize();
    exit(-1);
  }

  
  if (expcond) {
				// Set parameters for external dcond function
    EXPSCALE = acyl;
    HSCALE   = hcyl;
    ASHIFT   = ashift;
    eof      = 0;

    bool cache_ok = false;
  
    // Attempt to read EOF file from cache with override.  Will work
    // whether first time or restart.  Aborts if overridden cache is
    // not found.
    //
    if (eof_file.size()>0) {

      cache_ok = ortho->read_cache();
      
      if (!cache_ok) {
	if (myid==0)		// Diagnostic output . . .
	  std::cerr << "Cylinder: can not read explicitly specified EOF file <"
		    << EmpCylSL::CACHEFILE << ">" << std::endl
		    << "Cylinder: shamelessly aborting . . ." << std::endl;
	
	MPI_Abort(MPI_COMM_WORLD, 12);
      }
    }


    // Attempt to read EOF file from cache on restart
    //
    if (try_cache || restart) {

      cache_ok = ortho->read_cache();

      // Diagnostic output . . .
      //
      if (!cache_ok and myid==0)
	std::cerr << "Cylinder: can not read EOF file <"
		  << EmpCylSL::CACHEFILE << ">" << std::endl
		  << "Cylinder: will attempt to generate EOF file, "
		  << "this will take some time (e.g. hours) . . ."
		  << std::endl;
    }

    // On restart, abort if the cache is gone
    //
    if (restart && !cache_ok) {
      if (myid==0) 
	std::cerr << "Cylinder: can not read cache file on restart ... aborting"
		  << std::endl;
      MPI_Abort(MPI_COMM_WORLD, 13);
    }

    // Genererate eof if needed
    //
    if (!cache_ok) ortho->generate_eof(rnum, pnum, tnum, dcond);

    firstime = false;
  }

  // Make sure that all structures are initialized to start (e.g. for
  // multi- stepping but this should be done on 1st call to determine
  // coefs by default
  //
  ortho->setup_accumulation();

  
#ifdef DEBUG
  for (int n=0; n<numprocs; n++) {
    if (myid==n) {
      cout << endl << "Process " << myid << ": Cylinder parameters: "
	   << " nmax="        << nmax
	   << " lmax="        << lmax
	   << " mmax="        << mmax
	   << " ncylorder="   << ncylorder
	   << " rcylmin="     << rcylmin
	   << " rcylmax="     << rcylmax
	   << " acyl="        << acyl
	   << " hcyl="        << hcyl
	   << " expcond="     << expcond
	   << " pcavar="      << pcavar
	   << " pcaeof="      << pcaeof
	   << " nvtk="        << nvtk
	   << " npca="        << npca
	   << " npca0="       << npca0
	   << " pcadiag="     << pcadiag
	   << " eof_file="    << eof_file
	   << " logarithmic=" << logarithmic
	   << " vflag="       << vflag
	   << endl << endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
#else
  if (myid==0) {
    cout << endl << "Cylinder parameters: "
	 << " nmax="        << nmax
	 << " lmax="        << lmax
	 << " mmax="        << mmax
	 << " ncylorder="   << ncylorder
	 << " rcylmin="     << rcylmin
	 << " rcylmax="     << rcylmax
	 << " acyl="        << acyl
	 << " hcyl="        << hcyl
	 << " expcond="     << expcond
	 << " pcavar="      << pcavar
	 << " pcaeof="      << pcaeof
	 << " nvtk="        << nvtk
	 << " npca="        << npca
	 << " npca0="        << npca0
	 << " pcadiag="     << pcadiag
	 << " eof_file="    << eof_file
	 << " logarithmic=" << logarithmic
	 << " vflag="       << vflag
	 << endl << endl;
  }
#endif
      
  ncompcyl = 0;

  pos = new Vector [nthrds];
  frc = new Vector [nthrds];
  for (int i=0; i<nthrds; i++) {
    pos[i].setsize(1, 3);
    frc[i].setsize(1, 3);
  }

#ifdef DEBUG
  offgrid.resize(nthrds);
#endif

}

Cylinder::~Cylinder()
{
  delete ortho;
  delete [] pos;
  delete [] frc;
}

void Cylinder::initialize()
{
  try {
    // These first two should not be user settable . . . but need them for now
    //
    if (conf["rcylmin"   ])    rcylmin  = conf["rcylmin"   ].as<double>();
    if (conf["rcylmax"   ])    rcylmax  = conf["rcylmax"   ].as<double>();

    if (conf["acyl"      ])       acyl  = conf["acyl"      ].as<double>();
    if (conf["hcyl"      ])       hcyl  = conf["hcyl"      ].as<double>();
    if (conf["nmax"      ])       nmax  = conf["nmax"      ].as<int>();
    if (conf["lmax"      ])       lmax  = conf["lmax"      ].as<int>();
    if (conf["mmax"      ])       mmax  = conf["mmax"      ].as<int>();
    if (conf["ncylnx"    ])     ncylnx  = conf["ncylnx"    ].as<int>();
    if (conf["ncylny"    ])     ncylny  = conf["ncylny"    ].as<int>();
    if (conf["ncylr"     ])      ncylr  = conf["ncylr"     ].as<int>();
    if (conf["ncylorder" ])  ncylorder  = conf["ncylorder" ].as<int>();
    if (conf["ncylrecomp"]) ncylrecomp  = conf["ncylrecomp"].as<int>();
    if (conf["npca"      ])       npca  = conf["npca"      ].as<int>();
    if (conf["npca0"     ])      npca0  = conf["npca0"     ].as<int>();
    if (conf["nvtk"      ])       nvtk  = conf["nvtk"      ].as<int>();
    if (conf["eof_file"  ])   eof_file  = conf["eof_file"  ].as<std::string>();
    if (conf["vflag"     ])      vflag  = conf["vflag"     ].as<int>();
    
    if (conf["rnum"      ])       rnum  = conf["rnum"      ].as<int>();
    if (conf["pnum"      ])       pnum  = conf["pnum"      ].as<int>();
    if (conf["tnum"      ])       tnum  = conf["tnum"      ].as<int>();
    if (conf["ashift"    ])     ashift  = conf["ashift"    ].as<double>();
    if (conf["expcond"   ])    expcond  = conf["expcond"   ].as<bool>();
    if (conf["logr"      ]) logarithmic = conf["logr"      ].as<bool>();
    if (conf["pcavar"    ])     pcavar  = conf["pcavar"    ].as<bool>();
    if (conf["pcaeof"    ])     pcaeof  = conf["pcaeof"    ].as<bool>();
    if (conf["pcavtk"    ])     pcavtk  = conf["pcavtk"    ].as<bool>();
    if (conf["pcadiag"   ])    pcadiag  = conf["pcadiag"   ].as<bool>();
    if (conf["try_cache" ])  try_cache  = conf["try_cache" ].as<bool>();
    if (conf["density"   ])    density  = conf["density"   ].as<bool>();
    if (conf["cmap"      ])       cmap  = conf["cmap"      ].as<bool>();
    
    if (conf["self_consistent"])
      self_consistent = conf["self_consistent"].as<bool>();
  }
  catch (YAML::Exception & error) {
    if (myid==0) std::cout << "Error parsing Cylinder parameters: "
			   << error.what() << std::endl
			   << std::string(60, '-') << std::endl
			   << "Config node"        << std::endl
			   << std::string(60, '-') << std::endl
			   << conf                 << std::endl
			   << std::string(60, '-') << std::endl;
    MPI_Finalize();
    exit(-1);
  }
}

void Cylinder::get_acceleration_and_potential(Component* C)
{
  nvTracerPtr tPtr;
  if (cuda_prof)
    tPtr = nvTracerPtr(new nvTracer("Cylinder::get_acceleration"));

  std::chrono::high_resolution_clock::time_point start0, start1, finish0, finish1;

  start0 = std::chrono::high_resolution_clock::now();

#ifdef DEBUG
  cout << "Process " << myid 
       << ": in Cylinder::get_acceleration_and_potential" << endl;
#endif
				
  cC = C;

  //====================================================
  // Accel & pot using previously computed coefficients 
  //====================================================

  if (use_external) {
    nvTracerPtr tPtr1;
    if (cuda_prof) {
      tPtr1 = nvTracerPtr(new nvTracer("Cylinder: in external"));
    }

    MPL_start_timer();
    determine_acceleration_and_potential();
    MPL_stop_timer();

    use_external = false;

    return;
  }


  //======================================
  // Determine potential and acceleration 
  //======================================

  MPL_start_timer();

  determine_acceleration_and_potential();

  MPL_stop_timer();

  //=======================
  // Recompute PCA analysis
  //=======================

				// No recomputation ever if the
  if (!expcond) {		// basis has been precondtioned

				// Only do this check only once per
				// multistep; might as well be at 
				// the end of the multistep sequence
    if ((multistep==0 || mstep==0) && !initializing) {
      ncompcyl++;
      if (ncompcyl == ncylrecomp) {
	ncompcyl = 0;
	eof = 1;
	determine_coefficients();
      }
    }

  }


  //=================
  // Debugging output
  //=================
  if (VERBOSE>3 && myid==1 && component->EJ) {
    string toutfile = homedir + "test.orientation." + runtag;
    ofstream debugf(toutfile.c_str(), ios::app);
    Vector axis = component->orient->currentAxis();
    debugf << tnow << " "
	   << component->orient->currentAxis()[1] << " " 
	   << component->orient->currentAxis()[2] << " " 
	   << component->orient->currentAxis()[3] << " " 
	   << component->orient->currentAxisVar() << " "
	   << component->orient->currentCenter()[1] << " " 
	   << component->orient->currentCenter()[2] << " " 
	   << component->orient->currentCenter()[3] << " " 
	   << component->orient->currentCenterVar() << " "
	   << component->orient->currentCenterVarZ() << " "
	   << component->orient->currentE() << " "
	   << component->orient->currentUsed()
	   << endl;
  }

}

void * Cylinder::determine_coefficients_thread(void * arg)
{
  double r, r2, phi, R2;
  double xx, yy, zz, mas;
  double Rmax2 = rcylmax*rcylmax*acyl*acyl;

  int id = *((int*)arg);
  int indx, nbeg, nend, nbodies;

  use[id] = 0;
  cylmass0[id] = 0.0;

  thread_timing_beg(id);

  vector<double> ctr;
  if (mix) mix->getCenter(ctr);

  if (eof) {

    // Will use all of the bodies independent of level
    //
    nbodies = cC->Number();
    
    if (nbodies==0) {
      thread_timing_end(id);
      return (NULL);
    }
    nbeg = nbodies*id/nthrds;
    nend = nbodies*(id+1)/nthrds;

    unsigned indx;
    PartMapItr n = cC->Particles().begin();
    for (int i=0; i<nbeg; i++) n++; // Move to beginning iterator

    for (int i=nbeg; i<nend; i++, n++) {

      indx = n->first;

      // Frozen particles don't contribute to field
      //
      if (cC->freeze(indx)) continue;
    
      if (mix) {
	for (int j=0; j<3; j++) 
	  pos[id][j+1] = cC->Pos(indx, j, Component::Local) - ctr[j];
      } else {
	for (int j=0; j<3; j++) 
	  pos[id][j+1] = cC->Pos(indx, j, 
				 Component::Local | Component::Centered);
      }
      
      if ( (cC->EJ & Orient::AXIS) && !cC->EJdryrun) 
	pos[id] = cC->orient->transformBody() * pos[id];

      xx = pos[id][1];
      yy = pos[id][2];
      zz = pos[id][3];

      r2 = xx*xx + yy*yy;
      r = sqrt(r2);
      R2 = r2 + zz*zz;
    
      if ( R2 < Rmax2) {

	mas = cC->Mass(indx);
	phi = atan2(yy, xx);

	ortho->accumulate_eof(r, zz, phi, mas, id, cC->Part(indx)->level);
	
	use[id]++;
	cylmass0[id] += mas;

      } 
    }

  } else {

    nbodies = cC->levlist[mlevel].size();
    
    if (nbodies==0) {
      thread_timing_end(id);
      return (NULL);
    }
    nbeg = nbodies*id/nthrds;
    nend = nbodies*(id+1)/nthrds;

    double adb = component->Adiabatic();

    for (int i=nbeg; i<nend; i++) {

      indx = cC->levlist[mlevel][i];

      // Frozen particles don't contribute to field
      //
      if (cC->freeze(indx)) continue;
    
      for (int j=0; j<3; j++) 
	pos[id][j+1] = cC->Pos(indx, j, Component::Local | Component::Centered);

      if ( (cC->EJ & Orient::AXIS) && !cC->EJdryrun) 
	pos[id] = cC->orient->transformBody() * pos[id];

      xx = pos[id][1];
      yy = pos[id][2];
      zz = pos[id][3];

      r2 = xx*xx + yy*yy;
      r = sqrt(r2);
      R2 = r2 + zz*zz;
    
      if ( R2 < Rmax2) {

	mas = cC->Mass(indx) * adb;
	phi = atan2(yy, xx);

	ortho->accumulate(r, zz, phi, mas, indx, id, mlevel);

	use[id]++;
	cylmass0[id] += mas;
	
      } else {

	if (VERBOSE>6) {
	  cout << "Process " << myid 
	       << ": r^2=" << R2
	       << " max r^2=" << Rmax2 
	       << " r2=" << r2 
	       << " z2=" << zz*zz 
	       << " m=" << cylmass0[id] 
	       << " eof=" << eof
	       << endl;

	  if (std::isnan(R2)) {
	    cout << endl;
	    cC->orient->transformBody().print(cout);
	    cout << endl;
	    cC->orient->currentAxis().print(cout);
	    cout << endl;
	    MPI_Abort(MPI_COMM_WORLD, -1);
	  }
	}
      }
      
    }
  }

  thread_timing_end(id);

  return (NULL);
}


void Cylinder::determine_coefficients(void)
{
  nvTracerPtr tPtr;
  if (cuda_prof)
    tPtr = nvTracerPtr(new nvTracer("Cylinder::determine_coefficients"));

  std::chrono::high_resolution_clock::time_point start0, start1, finish0, finish1;

  start0 = std::chrono::high_resolution_clock::now();

  static char routine[] = "determine_coefficients_Cylinder";

  if (!self_consistent && !firstime_coef && !initializing) return;

  if (!expcond && firstime) {
				// Try to read cache
    bool cache_ok = false;
    if (try_cache || restart) {
      cache_ok = ortho->read_cache();
				// For a restart, cache must be read
				// otherwise, abort
      if (restart && !cache_ok) {
	if (myid==0) 
	  cerr << "Cylinder: can not read cache file on restart" << endl;
	MPI_Abort(MPI_COMM_WORLD, -1);
      }
    }

    eof = cache_ok ? 0 : 1;
    firstime = false;
				// If we can't read the cache, or the cache
				// does not match the requested parameters,
				// remake the emperical orthogonal basis 
    if (eof) {
      determine_coefficients_eof();
    }
  }

  if ( (pcavar or pcaeof) and pcainit) {
    EmpCylSL::PCAVAR = pcavar;
    EmpCylSL::PCAEOF = pcaeof;
    EmpCylSL::PCAVTK = pcavtk;
    EmpCylSL::VTKFRQ = nvtk;
    std::ostringstream sout;
    if (pcadiag) 
      sout << runtag << ".pcadiag." << cC->id << "." << cC->name;
    ortho->setHall(sout.str(), component->nbodies_tot);
    if (myid==0) {
      std::cout << "Cylinder: PCA initialized";
      if (pcadiag) 
	std::cout << ", writing diagnostic output to <"
		  << sout.str() << ">";
      std::cout << std::endl;
    }
    pcainit = false;
  }

  if (pcavar or pcaeof) {
    if (this_step >= npca0)
      compute = (mstep == 0) && !( (this_step-npca0) % npca);
    else
      compute = false;
  }

  ortho->setup_accumulation(mlevel);

  cylmass0.resize(nthrds);

#ifdef LEVCHECK
  for (int n=0; n<numprocs; n++) {
    if (n==myid) {
      if (myid==0) cout << "------------------------" << endl
			<< "Level check in Cylinder:" << endl 
			<< "------------------------" << endl;
      cout << setw(4) << myid << setw(4) << mlevel << setw(4) << eof;
      if (cC->levlist[mlevel].size())
	cout << setw(12) << cC->levlist[mlevel].size()
	     << setw(12) << cC->levlist[mlevel].front()
	     << setw(12) << cC->levlist[mlevel].back() << endl;
      else
	cout << setw(12) << cC->levlist[mlevel].size()
	     << setw(12) << (int)(-1)
	     << setw(12) << (int)(-1) << endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
  MPI_Barrier(MPI_COMM_WORLD);
  if (myid==0) cout << endl;
#endif
    
#if HAVE_LIBCUDA==1
  if (component->cudaDevice>=0) {
    start1 = std::chrono::high_resolution_clock::now();
    
    if (mstep==0) {
      std::fill(use.begin(), use.end(), 0.0);
      std::fill(cylmass0.begin(), cylmass0.end(), 0.0);
    }

    if (cC->levlist[mlevel].size()) {
      determine_coefficients_cuda(compute);
      DtoH_coefs(mlevel);
    }

    finish1 = std::chrono::high_resolution_clock::now();
  } else {    
    exp_thread_fork(true);
  }
#else
				// Threaded coefficient accumulation loop
  exp_thread_fork(true);
#endif
				// Accumulate counts and mass used to
				// determine coefficients
  int use1=0, use0=0;
  double cylmassT1=0.0, cylmassT0=0.0;
  
  for (int i=0; i<nthrds; i++) {
    use1      += use[i];
    cylmassT1 += cylmass0[i];
  }
				// Turn off timer so as not bias by 
				// communication barrier
  MPL_stop_timer();

  if (tnow==resetT) {

    MPI_Allreduce ( &use1, &use0, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce ( &cylmassT1, &cylmassT0, 1, MPI_DOUBLE, MPI_SUM, 
		    MPI_COMM_WORLD );

    used    += use0;
    cylmass += cylmassT0;
  }

  MPL_start_timer();

				// Make the coefficients for this level
  if (multistep==0 || !self_consistent) {
    ortho->make_coefficients(compute);
  } else if (mlevel==multistep) {
    ortho->make_coefficients(mfirst[mstep], compute);
    compute_multistep_coefficients();
  }

    if (myid==0)
    std::cout << "Check pca 0: pcavar=" << pcavar << " pcaeof=" <<
      pcaeof << " mlevel=" << mlevel << " multistep=" << multistep << " compute=" << compute << endl;

  if ((pcavar or pcaeof) and mlevel==multistep) {
    if (myid==0)
    std::cout << "Check pca 1: pcavar=" << pcavar << " pcaeof=" <<
      pcaeof << " mlevel=" << mlevel << " multistep=" << multistep << " compute=" << compute << endl;

    ortho->pca_hall(compute);
  }
    
  //=========================
  // Dump basis on first call
  //=========================

  if ( dump_basis and (this_step==0 || (expcond and ncompcyl==0) )
       && ortho->coefs_made_all() && !initializing) {

    if (myid == 0 and multistep==0 || mstep==0) {
      
      nvTracerPtr tPtr2;
      if (cuda_prof) {
	tPtr2 = nvTracerPtr(new nvTracer("Cylinder::dump basis"));
      }

      ortho->dump_basis(runtag.c_str(), this_step);
      
      ostringstream dumpname;
      dumpname << "images" << "." << runtag << "." << this_step;
      ortho->dump_images(dumpname.str(), 5.0*acyl, 5.0*hcyl, 64, 64, true);
      //
      // This next call is ONLY for deep debug
      //
      // dump_mzero(runtag.c_str(), this_step);
    }
  }

  print_timings("Cylinder: coefficient timings");

  finish0 = std::chrono::high_resolution_clock::now();
  
#if HAVE_LIBCUDA==1
  if (component->timers) {
    std::chrono::duration<double> duration0 = finish0 - start0;
    std::chrono::duration<double> duration1 = finish1 - start1;
    
    std::cout << std::string(60, '=') << std::endl;
    std::cout << "== Coefficient evaluation [Cylinder] level="
	      << mlevel << std::endl;
    std::cout << std::string(60, '=') << std::endl;
    std::cout << "Time in CPU: " << duration0.count()-duration1.count() << std::endl;
    if (cC->cudaDevice>=0) {
      std::cout << "Time in GPU: " << duration1.count() << std::endl;
    }
    std::cout << std::string(60, '=') << std::endl;
  }
#endif

  //================================
  // Dump coefficients for debugging
  //================================

  if (false and myid==0 and mstep==0 and mlevel==multistep) {
    std::cout << std::string(60, '-') << std::endl
	      << "-- Cylinder T=" << std::setw(16) << tnow << std::endl
	      << std::string(60, '-') << std::endl;
    ortho->dump_coefs(std::cout);
    std::cout << std::string(60, '-') << std::endl;
  }

  firstime_coef = false;
}


void Cylinder::determine_coefficients_eof(void)
{
  if (eof==0) return;

  static char routine[] = "determine_coefficients_eof_Cylinder";
  
  ortho->setup_eof();
  ortho->setup_accumulation();

  cylmass = 0.0;
  if (myid==0) cerr << "Cylinder: setup for eof\n";

  cylmass0.resize(nthrds);

				// Threaded coefficient accumulation loop
  exp_thread_fork(true);

				// Accumulate counts and mass used to
				// determine coefficients
  int use0=0, use1=0;
  double cylmassT1=0.0, cylmassT0=0.0;

  for (int i=0; i<nthrds; i++) {
    use1 += use[i];
    cylmassT1 += cylmass0[i];
  }

				// Turn off timer so as not bias by 
				// communication barrier
  MPL_stop_timer();

  MPI_Allreduce ( &use1, &use0, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce ( &cylmassT1, &cylmassT0, 1, MPI_DOUBLE, MPI_SUM, 
		  MPI_COMM_WORLD );

  MPL_start_timer();

  if (myid==0) cerr << "Cylinder: eof grid mass=" << cylmassT0 
		    << ", number=" << use0 << "\n";

  ortho->make_eof();
  if (myid==0) cerr << "Cylinder: eof computed\n";

  ortho->make_coefficients();
  if (myid==0) cerr << "Cylinder: coefs computed\n";

  eof = 0;
}


void check_force_values(double phi, double p, double fr, double fz, double fp)
{
  if (
      std::isinf(phi) || std::isnan(phi) ||
      std::isinf(p  ) || std::isnan(p  ) ||
      std::isinf(fr ) || std::isnan(fr ) ||
      std::isinf(fz ) || std::isnan(fz ) ||
      std::isinf(fp ) || std::isnan(fp ) ) 
    {
      cerr << "check_force_values: Illegal value\n";
    }
}


void * Cylinder::determine_acceleration_and_potential_thread(void * arg)
{
  double r, r2, r3, phi;
  double xx, yy, zz;
  double p, p0, fr, fz, fp;

  const double ratmin = 0.75;
  const double maxerf = 3.0;
  const double midpt = ratmin + 0.5*(1.0 - ratmin);
  const double rsmth = 0.5*(1.0 - ratmin)/maxerf;

  double R2 = rcylmax*rcylmax*acyl*acyl, ratio, frac, cfrac, mfactor = 1.0;

  vector<double> ctr;
  if (mix) mix->getCenter(ctr);

  int id = *((int*)arg);

#ifdef DEBUG
  static bool firstime = true;
  ofstream out;
  int flg;
  if (firstime && myid==0 && id==0) out.open("debug.tst");
#endif

  thread_timing_beg(id);

  // If we are multistepping, compute accel only at or below <mlevel>
  //
  for (unsigned lev=mlevel; lev<=multistep; lev++) {

    unsigned nbodies = cC->levlist[lev].size();

    if (nbodies==0) continue;

    int nbeg = nbodies*id/nthrds;
    int nend = nbodies*(id+1)/nthrds;
    
#ifdef DEBUG
    cout << "Process " << myid << " id=" << id 
	 << ": nbodies=" << nbodies
	 << " lev=" << lev
	 << " nbeg=" << nbeg
	 << " nend=" << nend << endl;
#endif

    for (int q=nbeg; q<nend; q++) {

      unsigned indx = cC->levlist[lev][q];

      if (indx<1 || indx>cC->nbodies_tot) {
	cout << "Process " << myid << " id=" << id 
	     << ": index error in Cylinder q=" << q
	     << " indx=" << indx << endl;
      }

      if (mix) {

	if (use_external) {
	  cC->Pos(&pos[id][1], indx, Component::Inertial);
	  component->ConvertPos(&pos[id][1], Component::Local);
	} else
	  cC->Pos(&pos[id][1], indx, Component::Local);

	mfactor = mix->Mixture(&pos[id][1]);
	for (int k=1; k<=3; k++) pos[id][k] -= ctr[k-1];

      } else {

	if (use_external) {
	  cC->Pos(&pos[id][1], indx, Component::Inertial);
	  component->ConvertPos(&pos[id][1], Component::Local | Component::Centered);
	} else
	  cC->Pos(&pos[id][1], indx, Component::Local | Component::Centered);

      }

      if ( (component->EJ & Orient::AXIS) && !component->EJdryrun) 
	pos[id] = component->orient->transformBody() * pos[id];

      xx    = pos[id][1];
      yy    = pos[id][2];
      zz    = pos[id][3];
      
      r2    = xx*xx + yy*yy;
      r     = sqrt(r2) + DSMALL;
      phi   = atan2(yy, xx);

      ratio = sqrt( (r2 + zz*zz)/R2 );

      if (ratio >= 1.0) {
	cfrac      = 1.0 - mfactor;
	frc[id][1] = 0.0;
	frc[id][2] = 0.0;
	frc[id][3] = 0.0;
      } else if (ratio > ratmin) {
	frac  = 0.5*(1.0 - erf( (ratio - midpt)/rsmth )) * mfactor;
	cfrac = 1.0 - frac;
      } else {
	frac  = mfactor;
      }
	
      if (ratio < 1.0) {

	ortho->accumulated_eval(r, zz, phi, p0, p, fr, fz, fp);
#ifdef DEBUG
	check_force_values(phi, p, fr, fz, fp);
#endif
	frc[id][1] = ( fr*xx/r - fp*yy/r2 ) * frac;
	frc[id][2] = ( fr*yy/r + fp*xx/r2 ) * frac;
	frc[id][3] = fz * frac;
	
#ifdef DEBUG
	flg = 1;
#endif
      }

      if (ratio > ratmin) {

	r3 = r2 + zz*zz;
	p = -cylmass/sqrt(r3);	// -M/r
	fr = p/r3;		// -M/r^3

	frc[id][1] += xx*fr * cfrac;
	frc[id][2] += yy*fr * cfrac;
	frc[id][3] += zz*fr * cfrac;

#ifdef DEBUG
	offgrid[id]++;
	flg = 2;
#endif
      }
    
      if (use_external)
	cC->AddPotExt(indx, p);
      else
	cC->AddPot(indx, p);

      if ( (component->EJ & Orient::AXIS) && !component->EJdryrun) 
	frc[id] = component->orient->transformOrig() * frc[id];

      for (int j=0; j<3; j++) cC->AddAcc(indx, j, frc[id][j+1]);

#ifdef DEBUG
      if (firstime && myid==0 && id==0 && q < 5) {
	out << setw(9)  << q          << endl
	    << setw(9)  << indx       << endl
	    << setw(9)  << flg        << endl
	    << setw(18) << xx         << endl
	    << setw(18) << yy         << endl
	    << setw(18) << zz         << endl
	    << setw(18) << frc[0][1]  << endl
	    << setw(18) << frc[0][2]  << endl
	    << setw(18) << frc[0][3]  << endl;
      }
#endif
    }
  }

#ifdef DEBUG
  firstime = false;		// DEBUG
#endif

  thread_timing_end(id);

  return (NULL);
}

static int ocf = 0;

void Cylinder::determine_acceleration_and_potential(void)
{
  nvTracerPtr tPtr;
  if (cuda_prof)
    tPtr = nvTracerPtr(new nvTracer("Cylinder::determine_acceleration"));

  std::chrono::high_resolution_clock::time_point start0, start1, finish0, finish1;

  start0 = std::chrono::high_resolution_clock::now();

  static char routine[] = "determine_acceleration_and_potential_Cyl";
  
  if (use_external == false) {

    if (multistep && (self_consistent || initializing)) {
      compute_multistep_coefficients();
    }

  }

#ifdef DEBUG
  for (int i=0; i<nthrds; i++) offgrid[i] = 0;
  cout << "Process " << myid << ": about to fork" << endl;
#endif

#if HAVE_LIBCUDA==1
  if (cC->cudaDevice>=0) {
    start1 = std::chrono::high_resolution_clock::now();
    //
    // Copy coeficients from this component to device
    //
    HtoD_coefs();
    //
    // Do the force computation
    //
    determine_acceleration_cuda();
    finish1 = std::chrono::high_resolution_clock::now();
  } else {
    exp_thread_fork(false);
  }
#else
  exp_thread_fork(false);
#endif

#ifdef DEBUG
  cout << "Cylinder: process " << myid << " returned from fork" << endl;
  int offtot=0;
  for (int i=1; i<nthrds; i++) offgrid[0] += offgrid[i];
  MPI_Reduce(&offgrid[0], &offtot, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  if (myid==0) {
    if (use_external)
      cout << endl << "T=" << tnow << "  external offgrid=" << offtot << endl;
    else
      cout << endl << "T=" << tnow << "  self offgrid=" << offtot << endl;
  }    

  unsigned long imin = std::numeric_limits<unsigned long>::max();
  unsigned long imax = 0, kmin = imin, kmax = 0;

  for (auto p : cC->Particles()) {
    imin = std::min<unsigned long>(imin, p.first);
    imax = std::max<unsigned long>(imax, p.first);
    kmin = std::min<unsigned long>(kmin, p.second->indx);
    kmax = std::max<unsigned long>(kmax, p.second->indx);
  }

  cout << "Cylinder: process " << myid << " name=<" << cC->name << "> bodies ["
       << kmin << ", " << kmax << "], ["
       << kmin << ", " << kmax << "]"
       << " #=" << cC->Particles().size() << endl;
#endif

  print_timings("Cylinder: acceleration timings");


# if HAVE_LIBCUDA
  if (component->timers) {
    auto finish0 = std::chrono::high_resolution_clock::now();
  
    std::chrono::duration<double> duration0 = finish0 - start0;
    std::chrono::duration<double> duration1 = finish1 - start1;
    std::chrono::duration<double> duration2 = start1  - start0;

    std::cout << std::string(60, '=') << std::endl;
    std::cout << "== Force evaluation [Cylinder::" << cC->name
	      << "] level=" << mlevel << std::endl;
    std::cout << std::string(60, '=') << std::endl;
    std::cout << "Time in CPU: " << duration0.count()-duration1.count() << std::endl;
    if (cC->cudaDevice>=0) {
      std::cout << "Time in GPU: " << duration1.count() << std::endl;
      std::cout << "Time before: " << duration2.count() << std::endl;
    }
    std::cout << std::string(60, '=') << std::endl;
  }
#endif
}

void Cylinder::
determine_fields_at_point_sph(double r, double theta, double phi,
			      double *tdens0, double *tpotl0, 
			      double *tdens, double *tpotl, 
			      double *tpotr, double *tpott, 
			      double *tpotp)

{
  double R = r*sin(theta);
  double z = r*cos(theta);
  double tpotR, tpotZ;

  determine_fields_at_point_cyl(R, z, phi, tdens0, tpotl0, tdens, tpotl, 
				&tpotR, &tpotZ, tpotp);
  
  *tpotr =   tpotR*sin(theta) + tpotZ*cos(theta) ;
  *tpott = (-tpotZ*sin(theta) + tpotR*cos(theta) )/(r+1.0e-10);
}



void Cylinder::determine_fields_at_point_cyl(double r, double z, double phi,
					     double *tdens0, double *tpotl0, 
					     double *tdens, double *tpotl, 
					     double *tpotr, double *tpotz, double *tpotp)
{
  ortho->accumulated_eval(r, z, phi, *tpotl0, *tpotl, *tpotr, *tpotz, *tpotp);
  // Accumulated eval returns forces not potential gradients
  *tpotr *= -1.0;
  *tpotz *= -1.0;
  *tpotp *= -1.0;
  if (density)
    *tdens = ortho->accumulated_dens_eval(r, z, phi, *tdens0);
  else
    *tdens = 0.0;
}

				// Dump coefficients to a file
void Cylinder::dump_coefs(ostream& out)
{
  ortho->dump_coefs_binary(out, tnow);
}

				// Density debug
#include <fstream>

void Cylinder::dump_mzero(const string& name, int step)
{
  const double RMAX = 5.0*acyl;
  const double ZMAX = 5.0*hcyl;
  double r, dr = RMAX/(ncylnx-1);
  double z, dz = 2.0*ZMAX/(ncylny-1);

  float zz;
  string label[] = {".dens0.", ".pot0.", ".fr0.", ".fz0."};
  ofstream** out = new ofstream* [4];

  for (int i=0; i<4; i++) {
    ostringstream ins;
    ins << name << label[i] << step;
    out[i] = new ofstream(ins.str().c_str());

    out[i]->write((char *)&ncylnx, sizeof(int));
    out[i]->write((char *)&ncylny, sizeof(int));
    out[i]->write((char *)&(zz=  0.0), sizeof(float));
    out[i]->write((char *)&(zz= RMAX), sizeof(float));
    out[i]->write((char *)&(zz=-ZMAX), sizeof(float));
    out[i]->write((char *)&(zz= ZMAX), sizeof(float));
  }


				// Ok, write data
  double p, p0, d0, fr, fz, fp;

  for (int k=0; k<ncylny; k++) {

    z = -ZMAX + dz*k;
	
    for (int j=0; j<ncylnx; j++) {
	  
      r = dr*j;

      zz = ortho->accumulated_dens_eval(r, z, 0.0, d0);
      out[0]->write((char *)&zz, sizeof(float));

      ortho->accumulated_eval(r, z, 0.0, p0, p, fr, fz, fp);
      out[1]->write((char *)&(zz=p ), sizeof(float));
      out[2]->write((char *)&(zz=fr), sizeof(float));
      out[3]->write((char *)&(zz=fz), sizeof(float));
    }
  }

				// Close and delete streams
  for (int i=0; i<4; i++) {
    out[i]->close();
    delete out[i];
  }
  delete [] out;

}

void Cylinder::multistep_update(int from, int to, Component* c, int i, int id)
{
  if (!self_consistent) return;

  if (c->freeze(i)) return;

  double mass = c->Mass(i) * component->Adiabatic();

  double xx = c->Pos(i, 0, Component::Local | Component::Centered);
  double yy = c->Pos(i, 1, Component::Local | Component::Centered);
  double zz = c->Pos(i, 2, Component::Local | Component::Centered);

  double r2 = (xx*xx + yy*yy);
  double  r = sqrt(r2);
  double phi = atan2(yy, xx);

  ortho->multistep_update(from, to, r, zz, phi, mass, id);

}


void Cylinder::multistep_reset() 
{ 
  used    = 0; 
  cylmass = 0.0;
  resetT  = tnow;
  ortho->reset_mass();
  ortho->multistep_reset();
}


static int idbg = 0;
void Cylinder::multistep_debug() 
{
  if (myid==0) {
    cout << endl;
    cout << setw(70) << setfill('-') << '-' << endl;
    ostringstream sout;
    sout << "--- multistep_debug: " << idbg << endl;
    cout << setw(70) << left << sout.str() << endl << right;
    cout << setw(70) << '-' << setfill(' ') << endl;

    ostringstream sout2;
    sout2 << "cylinder.coefs." << runtag << "." << ocf++;
    ofstream out(sout2.str().c_str());
    ortho->dump_coefs(out);
  }

  ortho->multistep_debug();

  if (myid==1) ortho->dump_basis(runtag.c_str(), idbg);

  ostringstream dumpname;
  dumpname << "images" << "." << runtag << "." << idbg;
  ortho->dump_images(dumpname.str(), 5.0*acyl, 5.0*hcyl, 64, 64, true);
  dump_mzero(runtag.c_str(), idbg);
  
  idbg++;
}
