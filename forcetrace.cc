/*****************************************************************************
 *  Description:
 *  -----------
 *
 *  Trace the force of a collection of orbits through a simulation
 *
 *
 *  Call sequence:
 *  -------------
 *
 *  Parameters:
 *  ----------
 *
 *
 *  Returns:
 *  -------
 *
 *
 *  Notes:
 *  -----
 *
 *
 *  By:
 *  --
 *
 *  MSP 12/20/15
 *
 ***************************************************************************/

				// C++/STL headers
#include <cstdlib>
#include <iostream>
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
#include <Vector.h>
#include <numerical.h>
#include <Particle.h>
#include <PSP.H>
#include <interp.h>
#include <EmpOrth9thd.h>
#include <massmodel.h>
#include <SphereSL.H>


#include <localmpi.h>
#include <ProgramParam.H>
#include <foarray.H>

program_option init[] = {
  {"FILELIST",          "string",       "filelist.dat",    "string of output EXP files to analyze"},
  {"NICE",		"int",		"0",		   "system priority"},
  {"RCYLMIN",		"double",	"0.001",	   "number of scale lengths for minimum radius in table"},
  {"RCYLMAX",		"double",	"20.0",		   "number of scale lengths for maximum radius in table"},
  {"NUMX",		"int",		"128",		   "number of radial table entries"},
  {"NUMY",		"int",		"64",		   "number of vertical table entries"},
  {"NORBMAX",		"int",		"1000",		   "maximum number of orbits to output"},
  {"RSCALE",		"double",	"0.01",		   "Radial scale length for basis expansion"},
  {"VSCALE",		"double",	"0.001",	   "Vertical Scale length for basis expansion"},
  {"LMAX",		"int",		"36",		   "Maximum harmonic order for spherical expansion"},
  {"NMAX",		"int",		"8",		   "Maximum radial order for spherical expansion"},
  {"MMAX",		"int",		"4",		   "Maximum harmonic order (disk)"},
  {"LMAX",		"int",		"6",		   "Maximum harmonic order (halo)"},
  {"NMAXH",		"int",		"20",		   "Maximum radial order (halo)"},
  {"NORDER",		"int",		"12",		   "Number of basis functions per azimuthal harmonic"},
  {"INITFLAG",		"int",		"1",	  	   "Train set on Component (1=stars)"},
  {"PARTFLAG",		"int",		"1",		   "Wakes using Component(s) [1=stars | 2=gas]"},
  {"OUTFILE",		"string",	"diskprof",	   "Filename prefix"},
  {"CACHEFILE",         "string",       ".eof.cache.file", "Cachefile name"},
  {"MODFILE",		"string",	"SLGridSph.model", "Halo modelfile"},
  {"",			"",		"",		   ""}
};


const char desc[] = "Compute disk potential, force and density profiles from PSP phase-space output files\n";


ProgramParam config(desc, init);

				// Variables not used but needed for linking
int VERBOSE = 4;
int nthrds = 1;
int this_step = 0;
unsigned multistep = 0;
unsigned maxlev = 100;
int mstep = 1;
int Mstep = 1;
vector<int> stepL(1, 0), stepN(1, 1);
char threading_on = 0;
pthread_mutex_t mem_lock;
pthread_mutex_t coef_lock;
string outdir, runtag;
double tpos = 0.0;
double tnow = 0.0;
  
enum ComponentType {Star=1, Gas=2, Halo=4};


void add_particles(ifstream* in, PSPDump* psp, int& nbods, vector<Particle>& p)
{
  if (myid==0) {

    int nbody = nbods/numprocs;
    int nbody0 = nbods - nbody*(numprocs-1);

				// Send number of bodies to be received
				// by eacn non-root node
    MPI_Bcast(&nbody, 1, MPI_INT, 0, MPI_COMM_WORLD);

    vector<Particle> t(nbody);
    vector<double> val(nbody);

    SParticle *part = psp->GetParticle(in);
    Particle bod;
    
    //
    // Root's particles
    //
    for (int i=0; i<nbody0; i++) {
      if (part==0) {
	cerr << "Error reading particle [n=" << 0 << ", i=" << i << "]" << endl;
	exit(-1);
      }
      bod.mass = part->mass();
      for (int k=0; k<3; k++) bod.pos[k] = part->pos(k);
      for (int k=0; k<3; k++) bod.vel[k] = part->vel(k);
      p.push_back(bod);

      part = psp->NextParticle(in);
    }


    //
    // Send the rest of the particles to the other nodes
    //
    for (int n=1; n<numprocs; n++) {
      
      for (int i=0; i<nbody; i++) {
	if (part==0) {
	  cerr << "Error reading particle [n=" 
	       << n << ", i=" << i << "]" << endl;
	  exit(-1);
	}
	t[i].mass = part->mass();
	for (int k=0; k<3; k++) t[i].pos[k] = part->pos(k);
	for (int k=0; k<3; k++) t[i].vel[k] = part->vel(k);
	part = psp->NextParticle(in);
      }
  
      for (int i=0; i<nbody; i++) val[i] = t[i].mass;

      
      MPI_Send(&val[0], nbody, MPI_DOUBLE, n, 11, MPI_COMM_WORLD);

      for (int i=0; i<nbody; i++) val[i] = t[i].pos[0];
      MPI_Send(&val[0], nbody, MPI_DOUBLE, n, 12, MPI_COMM_WORLD);

      for (int i=0; i<nbody; i++) val[i] = t[i].pos[1];
      MPI_Send(&val[0], nbody, MPI_DOUBLE, n, 13, MPI_COMM_WORLD);

      for (int i=0; i<nbody; i++) val[i] = t[i].pos[2];
      MPI_Send(&val[0], nbody, MPI_DOUBLE, n, 14, MPI_COMM_WORLD);

      for (int i=0; i<nbody; i++) val[i] = t[i].vel[0];
      MPI_Send(&val[0], nbody, MPI_DOUBLE, n, 15, MPI_COMM_WORLD);

      for (int i=0; i<nbody; i++) val[i] = t[i].vel[1];
      MPI_Send(&val[0], nbody, MPI_DOUBLE, n, 16, MPI_COMM_WORLD);

      for (int i=0; i<nbody; i++) val[i] = t[i].vel[2];
      MPI_Send(&val[0], nbody, MPI_DOUBLE, n, 17, MPI_COMM_WORLD);
    }

  } else {

    int nbody;
    MPI_Bcast(&nbody, 1, MPI_INT, 0, MPI_COMM_WORLD);
    vector<Particle> t(nbody);
    vector<double> val(nbody);
    //vector<float> val(nbody);
				// Get and pack

    MPI_Recv(&val[0], nbody, MPI_DOUBLE, 0, 11, 
	     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    for (int i=0; i<nbody; i++) t[i].mass = val[i];

    MPI_Recv(&val[0], nbody, MPI_DOUBLE, 0, 12, 
	     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    for (int i=0; i<nbody; i++) t[i].pos[0] = val[i];

    MPI_Recv(&val[0], nbody, MPI_DOUBLE, 0, 13, 
	     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    for (int i=0; i<nbody; i++) t[i].pos[1] = val[i];

    MPI_Recv(&val[0], nbody, MPI_DOUBLE, 0, 14, 
	     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    for (int i=0; i<nbody; i++) t[i].pos[2] = val[i];

    MPI_Recv(&val[0], nbody, MPI_DOUBLE, 0, 15, 
	     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    for (int i=0; i<nbody; i++) t[i].vel[0] = val[i];

    MPI_Recv(&val[0], nbody, MPI_DOUBLE, 0, 16, 
	     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    for (int i=0; i<nbody; i++) t[i].vel[1] = val[i];

    MPI_Recv(&val[0], nbody, MPI_DOUBLE, 0, 17, 
	     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    for (int i=0; i<nbody; i++) t[i].vel[2] = val[i];

    p.insert(p.end(), t.begin(), t.end());
  }

}

void partition(ifstream* in, PSPDump* psp, int cflag, vector<Particle>& p)
{
  p.erase(p.begin(), p.end());

  int nbods = 0;
  if (cflag & Star) {
    if (myid==0) {
      nbods = psp->CurrentDump()->nstar;
      psp->GetStar();

      add_particles(in, psp, nbods, p);
    } else {
        add_particles(in, psp, nbods, p);
    }      
  }

  if (cflag & Gas) {
    if (myid==0) {
      int nbods = psp->CurrentDump()->ngas;
      psp->GetGas();

      add_particles(in, psp, nbods, p);
    } else {
      add_particles(in, psp, nbods, p);
    }      
  }

  if (cflag & Halo) {
    if (myid==0) {
      int nbods = psp->CurrentDump()->ndark;

      PSPstanza *cur = psp->GetDark();
      in->seekg(cur->pspos, ios::beg);

      add_particles(in, psp, nbods, p);
    } else {
      add_particles(in, psp, nbods, p);
    }      
  }



  
}


void force_print(SphereSL& orthoh, EmpCylSL& ortho, vector<Particle>& part, string& outf)
{


  double r, r3, costh, phi, z, mass, lz;

    
   double p0, pp, fr, fz, fp;

   double p0h, p1h, d0h, d1h, plh, frh, fth, fph;

    // particle counter number
  int ncnt=0;



  ofstream out(outf, ios::out | ios::app);


    
      
  for (auto p=part.begin(); p!=part.end(); p++) {

    mass = p->mass;
    r = sqrt(p->pos[0]*p->pos[0] + p->pos[1]*p->pos[1]);
    r3 = sqrt(p->pos[0]*p->pos[0] + p->pos[1]*p->pos[1] + p->pos[2]*p->pos[2])+1.0e-18;
    phi = atan2(p->pos[1], p->pos[0]);
    costh = p->pos[2]/r3;
    z = p->pos[2];

    lz = p->pos[0]*p->vel[1] - p->pos[1]*p->vel[0];


    ortho.accumulated_eval(r, z, phi, p0, pp, fr, fz, fp);

    orthoh.all_eval(r3, costh, phi, d0h, d1h, p0h, p1h, frh, fth, fph);

    // should I do the cylindrical force conversion here? YES

    //out[0].write((char *)&fr, sizeof(double));
	//	out[1].write((char *)&fz, sizeof(double));
	//out[2].write((char *)&fp, sizeof(double));

    // single file version
    out.write((char *)&r, sizeof(float));
    out.write((char *)&phi, sizeof(float));
    out.write((char *)&z, sizeof(float));
    out.write((char *)&lz, sizeof(float));
    out.write((char *)&fr, sizeof(float));
    out.write((char *)&fp, sizeof(float));
    out.write((char *)&fz, sizeof(float));
    out.write((char *)&frh, sizeof(float));
    out.write((char *)&fph, sizeof(float));
    out.write((char *)&fth, sizeof(float));


	    // timing print
      if ( (ncnt % 100) == 0) cout << "\r>> " << ncnt << " <<" << flush;
      ncnt++;

      // put in floor during testing
      if (ncnt >= config.get<int>("NORBMAX")) break;
      
  } // end of particle loop

} // end of force_print



int
main(int argc, char **argv)
{
  
#ifdef DEBUG
  sleep(20);
#endif  
  
  // ==================================================
  // Parse command line or input parameter file
  // ==================================================
  
  if (config.parse_args(argc,argv)) return -1;
  
  
  // ==================================================
  // MPI preliminaries
  // ==================================================

  local_init_mpi(argc, argv);
  
  // ==================================================
  // Nice process
  // ==================================================

  if (config.get<int>("NICE")>0)
    setpriority(PRIO_PROCESS, 0, config.get<int>("NICE"));

  // ==================================================
  // Assemble the list of files
  // ==================================================

  vector<string> files;


  if (myid==0) {
  ifstream ifiles(config.get<string>("FILELIST"));
      if (!ifiles) {
        // if (myid==0) {
	//	bomb("OrbTrack: provided OLIST file cannot be opened\n");
	//    }
        }
      string filenames;
      while (1) {
          ifiles >> filenames;
          if (ifiles.eof() || ifiles.bad()) break;
           files.push_back(filenames);
         }
  }
  
  

  unsigned nfiles = files.size();


    
  MPI_Bcast(&nfiles, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
  
  for (unsigned n=0; n<nfiles; n++) {
    unsigned sz;
    
    if (myid==0) sz = files[n].size();
    
    MPI_Bcast(&sz, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
    
    if (myid==0) {
      char *c = const_cast<char*>(files[n].c_str());
      MPI_Bcast(c, sz+1, MPI_CHAR, 0, MPI_COMM_WORLD);
    }
    else {
      char l[sz+1];
      MPI_Bcast(&l[0], sz+1, MPI_CHAR, 0, MPI_COMM_WORLD);
      files.push_back(l);
    }
    
    MPI_Barrier(MPI_COMM_WORLD);
  }
  


  if (myid==0) cout  << "Found " <<  nfiles  << " file(s)." << endl;

  // ==================================================
  // PSP input stream
  // ==================================================

  int iok = 1; // check to see if files were opened properly

  
  MPI_Bcast(&iok, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (iok==0) {
    MPI_Finalize();
    exit(-1);
  }

  // open output file

  string outf  = config.get<string>("OUTFILE");

  // empty out junk before appending
  ofstream oute(outf, ios::out); // | ios::app);
  oute.close();


    // ==================================================
  // Make SL expansion
  // ==================================================

        SphericalModelTable halo(config.get<string>("MODFILE"));
        SphereSL::mpi = true;
  SphereSL::NUMR = 4000;

  // Do same partition except for halo

  
  // ==================================================
  // All processes will now compute the basis functions
  // *****Using MPI****
  // ==================================================


  // get expansion values from .eof header and report...force use of cachefile

  EmpCylSL::RMIN        = config.get<double>("RCYLMIN");
  EmpCylSL::RMAX        = config.get<double>("RCYLMAX");
  EmpCylSL::NUMX        = config.get<int>("NUMX");
  EmpCylSL::NUMY        = config.get<int>("NUMY");
  EmpCylSL::CMAP        = true;
  EmpCylSL::logarithmic = true;
  EmpCylSL::DENS        = true;
  EmpCylSL::CACHEFILE   = config.get<string>("CACHEFILE");
  EmpCylSL::SELECT      = true; // do the SNR methodology

                                // Create expansion
				//
  EmpCylSL ortho(config.get<int>("NMAX"), 
		 config.get<int>("LMAX"), 
		 config.get<int>("MMAX"), 
		 config.get<int>("NORDER"),
		 config.get<double>("RSCALE"),
		 config.get<double>("VSCALE"));

  vector<Particle> particles, particlesh;
  PSPDump *psp = 0;

  if (ortho.read_cache()==0) {
    if (myid==0) {
      cout << setw(70) << "Cached file required at this time." << endl;
    }
  } // If cached file is read, jump to here.
      


  for (int n=0; n<nfiles; n++) { 

        // make the spherical model

      SphereSL orthoh(&halo, config.get<int>("LMAX"), config.get<int>("NMAXH"));


   ifstream in;

   in.open(files[n].c_str());
    
  delete psp;
  psp = new PSPDump (&in, true);

  Dump *dump = psp->GetDump();


  int icnt = 0;

  if (dump) {





    //------------------------------------------------------------

    

    if (myid==0) {
      cout << "Beginning disk partition [time="
	   << dump->header.time << "] . . . " << flush;
    }

    if (in.rdstate() & ios::eofbit) {
      in.close();
      in.open(files[n].c_str());
    }

    partition(&in, psp, 1, particles);
    if (myid==0) cout << "done" << endl;
    MPI_Barrier(MPI_COMM_WORLD);

        //------------------------------------------------------------ 

    if (myid==0) {
      cout << "Beginning halo partition [time="
	   << dump->header.time << "] . . . " << flush;
    }

    if (in.rdstate() & ios::eofbit) {
      in.close();
      in.open(files[n].c_str());
    }

    partition(&in, psp, 4, particlesh);
    if (myid==0) cout << "done" << endl;
    MPI_Barrier(MPI_COMM_WORLD);

    //------------------------------------------------------------ 

    if (myid==0) cout << "Accumulating for disk basis . . . " << flush;
    ortho.accumulate(particles, 0, true);
    MPI_Barrier(MPI_COMM_WORLD);

    if (myid==0) cout << "done" << endl;

    //------------------------------------------------------------ 

    if (myid==0) cout << "Accumulating for halo basis . . . " << flush;
    orthoh.reset_coefs();
    for (auto &i : particlesh) {
    orthoh.accumulate(i.pos[0], i.pos[1], i.pos[2], i.mass);
      }
    MPI_Barrier(MPI_COMM_WORLD);

    if (myid==0) cout << "done" << endl;



    //------------------------------------------------------------ 

    if (myid==0) cout << "Making disk coefficients . . . " << flush;
    ortho.make_coefficients();
    MPI_Barrier(MPI_COMM_WORLD);
    if (myid==0) cout << "done" << endl;

        //------------------------------------------------------------ 

    if (myid==0) cout << "Making halo coefficients . . . " << flush;
    orthoh.make_coefs();
    MPI_Barrier(MPI_COMM_WORLD);
    if (myid==0) cout << "done" << endl;




    //------------------------------------------------------------

    

    // set this up to handle a particle list
    if (myid==0) {
      cout << "Doing Force Output . . . " << flush;
      force_print(orthoh, ortho, particles, outf);
    }
    MPI_Barrier(MPI_COMM_WORLD);


    if (myid==0) cout << "done" << endl;
    

  } // end of dump if

  // clear the halo accumulation

  // PSEUDOCODE
  // delete orthoh;

  //orthoh.~SphereSL();
    
  /*
  if (myid==0) {
    cout << "Clearing ortho . . ." << flush;
  }
    ortho.reset(config.get<int>("NMAX"), 
		 config.get<int>("LMAX"), 
		 config.get<int>("MMAX"), 
		 config.get<int>("NORDER"),
		 config.get<double>("RSCALE"),
		 config.get<double>("VSCALE"));

    if (myid==0)  cout << "done" << endl;
  */
  } // end of file list array
  
  MPI_Finalize();

  return 0;
}

