//
// Generate a two-power model with an optional core
//
// Compile string: g++ -o twopower -O3 twopower.cc
//
// (you will need g++ with a STL implementation, e.g egcs>=2.8).
//
// MDWeinberg 2/27/00
//

#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <getopt.h>

#include <string>		// from STL

using namespace std;

				// ------------------------------------------
				// Global parameters
				// ------------------------------------------

				// Default name new output model
string outputfile = "newmodel.dat";

				// ------------------------------------------
				// Global functions
				// ------------------------------------------

int 
main (int argc, char** argv)
{

  bool logarithmic = true;
  bool heggie = false;
  bool verbose = false;
  bool truncate = false;
  int number = 1000;
  double rmin = 1.0e-3;
  double rmax = 100.0;
  double rcore = 0.2;
  double alpha = 1.0;
  double beta = 2.0;
  double W = -0.5;
  double M = 1.0;
  double R = 1.0;
  double rtrunc = 1.0;
  double wtrunc = 0.2;
				// ------------------------------------------
				// Command line parsing
				// ------------------------------------------


  while (1) {

    int c = getopt(argc, argv, "o:n:c:a:b:i:r:t:w:W:M:R:SvL");
     
    if (c == -1) break;

    switch (c) {
    case 'o': outputfile.erase(); outputfile = optarg; break;
    case 'n': number = atoi(optarg); break;
    case 'c': rcore = atof(optarg); break;
    case 'a': alpha = atof(optarg); break;
    case 'b': beta = atof(optarg); break;
    case 'i': rmin = atof(optarg); break;
    case 'r': rmax = atof(optarg); break;
    case 't': rtrunc = atof(optarg); truncate=true; break;
    case 'w': wtrunc = atof(optarg); truncate=true; break;
    case 'W': W = atof(optarg); break;
    case 'M': M = atof(optarg); break;
    case 'R': R = atof(optarg); break;
    case 'S': heggie = heggie ? false : true; break;
    case 'L': logarithmic = logarithmic ? false : true; break;
    case 'v': verbose = verbose ? false : true; break;
    case '?': cerr << "usage: " << argv[0] << " -[o<name> n<number> c<rcore> a<alpha> b<beta> i<rmin> r<rmax> W<Wtot> M<Mtot> R<Rtot> t<rtrunc> w<wtrunc> L<1/0> v\n"; exit(-1);
    }
    
  }

				// ------------------------------------------
				// Allocate work space
				// ------------------------------------------

  double *r =  new double [number] - 1;
  double *d =  new double [number] - 1;
  double *m =  new double [number] - 1;
  double *pw = new double [number] - 1;
  double *p0 = new double [number] - 1;

				// ------------------------------------------
				// Open needed files
				// ------------------------------------------

  ofstream out(outputfile.c_str());
  if (!out) {
    cerr << "Error opening <" << outputfile << ">\n";
    exit(-1);
  }
				// ------------------------------------------
				// Make radial, density and mass array
				// ------------------------------------------

  double dr;
  if (logarithmic)
    dr = (log(rmax) - log(rmin))/(number - 1);
  else
    dr = (rmax - rmin)/(number - 1);

  for (int i=1; i<=number; i++) {
    if (logarithmic)
      r[i] = rmin*exp(dr*(i-1));
    else
      r[i] = rmin + dr*(i-1);
    d[i] = pow(r[i]+rcore, -alpha) * pow(r[i]+1.0, -beta);

    if (truncate) d[i] *= 0.5*(1.0 - erf((r[i]-rtrunc)/wtrunc));
  }

  m[1] = 0.0;
  pw[1] = 0.0;
  for (int i=2; i<=number; i++) {
    m[i] = m[i-1] +
      2.0*M_PI*(r[i-1]*r[i-1]*d[i-1] + r[i]*r[i]*d[i])*(r[i] - r[i-1]);
    pw[i] = pw[i-1] +
      2.0*M_PI*(r[i-1]*d[i-1] + r[i]*d[i])*(r[i] - r[i-1]);
  }

  for (int i=1; i<=number; i++) 
    p0[i] = -m[i]/(r[i]+1.0e-10) - (pw[number] - pw[i]);

  double W0 = 0.0;
  for (int i=1; i<=number; i++) 
    W0 += M_PI*(r[i-1]*r[i-1]*d[i-1]*p0[i-1] + 
		    r[i]*r[i]*d[i]*p0[i-1])*(r[i] - r[i-1]);

  cout << "orig PE = " << W0  << endl;

  double M0 = m[number];
  double R0 = r[number];


  //
  // Compute new scaling
  //

  double Gamma, Beta;

  if (heggie) {
    Beta = (W/W0) * (M0/M);
    Gamma = pow(W/W0,1.5) * pow(M0/M,3.5);
    out << "! Scaling:  W=" << W << "  M=" << M;
    }
  else {
    Beta = (M/M0) * (R0/R);
    Gamma = sqrt((M0*R0)/(M*R)) * (R0/R);
    out << "! Scaling:  R=" << R << "  M=" << M;
  }

  out << "  alpha=" << alpha << "  beta=" << beta 
      << "  rcore=" << rcore;
  if (truncate)
    out << "  rtrunc=" << rtrunc << "  wtrunc=" << wtrunc;
  out << endl;
  out << "! 1) = r   2) = rho   3) = M(r)   4) U(r) " << endl;
  out << setw(10) << number << endl;


  double rfac = pow(Beta,-0.25) * pow(Gamma,-0.5);
  double dfac = pow(Beta,1.5) * Gamma;
  double mfac = pow(Beta,0.75) * pow(Gamma,-0.5);
  double pfac = Beta;

  out.setf(ios::scientific);
  out.precision(12);

  for (int i=1; i<=number; i++) {
    out << setw(20) << r[i]*rfac;
    out << setw(20) << d[i]*dfac;
    out << setw(20) << m[i]*mfac;
    out << setw(20) << p0[i]*pfac;
    out << endl;
  }

  M0 = 0.0;
  W0 = 0.0;
  for (int i=1; i<=number; i++) {
    M0 += 2.0*M_PI*(r[i-1]*r[i-1]*d[i-1] + 
		    r[i]*  r[i]*  d[i]   ) * (r[i] - r[i-1]);
    W0 += M_PI*(r[i-1]*r[i-1]*d[i-1]*p0[i-1] + 
		r[i]*  r[i]*  d[i]*  p0[i-1]) * (r[i] - r[i-1]);
  }

  cout << "new M0 = " << M0*mfac << endl;
  cout << "new PE = " << W0*rfac*rfac*rfac*dfac*pfac << endl;
  cout << "-------- " << endl;
  cout << "Rfac = " << rfac << endl;
  cout << "Dfac = " << dfac << endl;
  cout << "Mfac = " << mfac << endl;
  cout << "Pfac = " << pfac << endl;
}
