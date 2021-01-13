/*****************************************************************************
 *  Description:
 *  -----------
 *
 *  Read in PSP files for a run and compute orbital diagnostics
 *  
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
 *  MSP 03/13/14
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
#include <massmodel.h>
#include <SphereSL.H>

#include <localmpi.h>
#include <ProgramParam.H>
#include <foarray.H>

program_option init[] = {
  {"PROJ",		"int",		"1",		"projection (1=cyl, 2=sphere)"},
  {"COMP",		"int",		"1",		"component (1=star, 2=gas, 4=dark)"},
  {"OUTFILE",		"string",	"lprof",	"filename prefix"},
  {"INFILE",		"string",	"OUT",		"phase space file"},
  {"RUNTAG",		"string",	"run",		"file containing desired indices for PSP output"},
  {"IBEG",		"int",		"0",		"first PSP index"},
  {"IEND",		"int",		"60",		"last PSP index"},
  {"ISKIP",		"int",		"1",		"skip PSP interval"},
  {"NSKIP",		"int",		"1",		"skip particle interval (if constructing)"},
  {"NORB",		"int",		"1000",		"number of orbits to count"},
  {"FULL",		"int",		"1",		"1=full orbits, 0=ensemble averages"},
  {"OLIST",		"string",      	"",		"orbitlist of particles to track"},
  {"",			"",		"",		""}
};


const char desc[] = "Track Orbits For a Specific Orbit List\n";

enum ProjectionType {Cylindrical=1, Spherical=2};
enum ComponentType  {Star=1, Gas=2, Halo=4};

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
  
int
main(int argc, char **argv)
{
  
  // ==================================================
  // Parse command line or input parameter file
  // ==================================================
  
  if (config.parse_args(argc,argv)) return -1;
  

  // ==================================================
  // MPI preliminaries
  // ==================================================

    local_init_mpi(argc, argv);

    //cout << setw(18) << myid <<endl;

  ofstream indx;
  ifstream in;

  vector<string> files;
				// Root nodes looks for existence of files
  if (myid==0) {
    for (int i=config.get<int>("IBEG"); i<=config.get<int>("IEND"); i += config.get<int>("ISKIP")) {
      ostringstream lab;
      lab << config.get<string>("INFILE") << "." 
	  << config.get<string>("RUNTAG") << "." 
	  << setw(5) << right << setfill('0') << i;
      ifstream in(lab.str().c_str());
      if (in) files.push_back(lab.str());
      else break;
      cout << "." << i << flush;
    }
    cout << endl;
  }

  unsigned nfiles = files.size();

  cout << nfiles << endl;

  MPI_Bcast(&nfiles, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
  for (unsigned n=0; n<nfiles; n++) {
    unsigned sz;
    if (myid==0) sz = files[n].size();
    MPI_Bcast(&sz, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
    if (myid==0) {
      char *c = const_cast<char*>(files[n].c_str());
      MPI_Bcast(c, sz+1, MPI_CHAR, 0, MPI_COMM_WORLD);
    } else {
      char l[sz+1];
      MPI_Bcast(&l[0], sz+1, MPI_CHAR, 0, MPI_COMM_WORLD);
      files.push_back(l);
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }


  int    diffskip = config.get<int>   ("NSKIP");
  int full = config.get<int> ("FULL");
  string orbitlist = config.get<string> ("OLIST");
 int nskip = config.get<int> ("NSKIP"); //100;
 int norb = config.get<int> ("NORB");//1000;
  
 





  if (myid==0) {
    cout << "Calculating spin parameter..." <<endl;

  }

 int    comp = config.get<int>   ("COMP");

  //--------------------------------------------------


  vector<double> times(nfiles);
  vector<double> spinz(nfiles);
  vector<double> spinthree(nfiles);
  

  for (int n=0; n<nfiles; n++) { // read over all the files stuffed into the array

    if (myid==0) {
    if (n==0) cout << "FILE" ;
	 cout << n << endl;
    }

    if (n % numprocs == myid) {

      // open first file
      ifstream in(files[n].c_str());
      PSPDump psp(&in, true);

      Dump *dump = psp.GetDump();
    
      if (dump) {

	times[n] = psp.CurrentTime();
	cout << setw(18) << times[n] << endl; // check where we are...

	// Do we need to close and reopen?
	if (in.rdstate() & ios::eofbit) {
	  in.close();
	  in.open(files[n].c_str());
	}

	int icnt = 0;
	vector<Particle> particles;

	psp.GetDark(); // this is only available for DM
	
	SParticle *p  =   psp.GetParticle(&in);

	  int i=0;


	    while(p) { // loop over all the particles...

	    vector<double> L(3);
	    L[0] = (p->pos[1]*p->vel[2] - p->pos[2]*p->vel[1]);
	    L[1] = (p->pos[2]*p->vel[0] - p->pos[0]*p->vel[2]);
	    L[2] = (p->pos[0]*p->vel[1] - p->pos[1]*p->vel[0]);


	      // just for z component
	      //double rtwo = p->pos[0]*p->pos[0]+p->pos[1]*p->pos[1];
	      spinz[n] += p->mass*L[2]; // z component of angular momentum
	      spinthree[n] += p->mass*(L[0]+L[1]+L[2]); //full angular momentum

	      //double rthree = p->pos[0]*p->pos[0]+p->pos[1]*p->pos[1]+p->pos[2]*p->pos[2];

				 // also should do this in r bins to see what the slope is roughly

	  p = psp.NextParticle(&in);
	  icnt++;

	  	} // while loop

      } // dump if

    } // matches the numprocs loop

  if (myid==0) {
    MPI_Reduce(MPI_IN_PLACE, &times[0], nfiles,
	       MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(MPI_IN_PLACE, &spinz[0], nfiles,
	       MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(MPI_IN_PLACE, &spinthree[0], nfiles,
	       MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  } else {
    MPI_Reduce(&times[0], 0, nfiles,
	       MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&spinz[0], 0, nfiles,
	       MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&spinthree[0], 0, nfiles,
	       MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); 

  }

  } // matches the iteration over nfiles


  if (myid==0) {
   
    for (int n=0; n<nfiles;n++) cout << setw(18) << times[n] << setw(18) << spinz[n] << setw(18) << spinthree[n] << endl;


    //         string outh = config.get<string>("OUTFILE") + ".AVG.dat";
    //    ofstream out3(outh.c_str());
    //for (int n=0; n<nbins; n++) {   
    // for (int m=0; m<nbins; m++) {
    //	 out3 << setw(18) << rvals[n] << setw(18) << bvals[m] << setw(18) << histoA[n][m] << endl;

    }





    MPI_Finalize();
  return 0;

}


