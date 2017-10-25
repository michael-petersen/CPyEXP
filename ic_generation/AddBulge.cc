#include <math.h>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>

#include "localmpi.h"
#include <massmodel.h>
#include <interp.h>
#include "AddBulge.h"

int AddBulge::number = 4000;
double AddBulge::Rmin = 1.0e-3;
bool AddBulge::use_mpi = false;
bool AddBulge::logarithmic = false;

AddBulge::AddBulge(AxiSymModel* halo, AxiSymModel* bulge)
{

  double rmin = max<double>(halo->get_min_radius(), Rmin);
  double rmax = halo->get_max_radius();

  r =  vector<double>(number);
  d =  vector<double>(number);
  m =  vector<double>(number);
  p =  vector<double>(number);

  vector<double> dt(number);
  vector<double> dm(number);
  vector<double> dp(number);
  vector<double> mm(number);
  vector<double> pw(number);

				// ------------------------------------------
				// Make radial, density and mass array
				// ------------------------------------------
  double dr;
  if (logarithmic)
    dr = (log(rmax) - log(rmin))/(number - 1);
  else
    dr = (rmax - rmin)/(number - 1);

  for (int i=0; i<number; i++) {
    if (logarithmic)
      r[i] = rmin*exp(dr*i);
    else
      r[i] = rmin + dr*i;

    m[i] = halo->get_mass(r[i]) + bulge->get_mass(r[i]);
    dm[i] = bulge->get_mass(r[i]);
  }

  for (int i=0; i<number; i++) {
    d[i] = halo->get_density(r[i]);
    dp[i] = drv2(r[i], r, dm)/(4.0*M_PI*r[i]*r[i]);
    dt[i] = dp[i] + d[i];
  }

  m[0] = 0.0;
  mm[0] = 0.0;
  pw[0] = 0.0;
  for (int i=1; i<number; i++) {
    m[i] = m[i-1] +
      2.0*M_PI*(r[i-1]*r[i-1]*d[i-1] + r[i]*r[i]*d[i])*(r[i] - r[i-1]);
    mm[i] = mm[i-1] +
      2.0*M_PI*(r[i-1]*r[i-1]*dt[i-1] + r[i]*r[i]*dt[i])*
      (r[i] - r[i-1]);
    pw[i] = pw[i-1] +
      2.0*M_PI*(r[i-1]*dt[i-1] + r[i]*dt[i])*(r[i] - r[i-1]);
  }

  for (int i=0; i<number; i++) 
    p[i] = -mm[i]/(r[i]+1.0e-10) - (pw[number-1] - pw[i]);

#ifdef DEBUG
 {
   ostringstream outf;
   if (use_mpi)
     outf << "test_addbulge." << myid << '\0';
   else
     outf << "test_addbulge.dat\0";
   ofstream out(outf.str().c_str());
   for (int i=0; i<number; i++) {
     out 
       << setw(15) << r[i] 
       << setw(15) << d[i] 
       << setw(15) << m[i] 
       << setw(15) << p[i] 
       << setw(15) << dp[i] 
       << setw(15) << dm[i] 
       << setw(15) << mm[i] 
       << endl;
   }
   out.close();
 }
#endif

  mod = new SphericalModelTable(number, &r[0]-1, &d[0]-1, &m[0]-1, &p[0]-1);

  {
    ostringstream outf;
    if (use_mpi)
      outf << "test_addbulge_mod." << myid << '\0';
    else
      outf << "test_addbulge_mod.dat\0";
    ofstream out(outf.str().c_str());
    for (int i=1; i<number; i++) {
      out 
	<< setw(15) << 0.5*(r[i-1]+r[i])
	<< setw(15) << mod->get_density(0.5*(r[i-1]+r[i]))
	<< setw(15) << mod->get_mass(0.5*(r[i-1]+r[i]))
	<< setw(15) << mod->get_pot(0.5*(r[i-1]+r[i]))
	<< endl;
    }
    out.close();
  }

}

AddBulge::~AddBulge()
{
  delete mod;
}

