// This may look like C code, but it is really -*- C++ -*-

//
// 3D Miyamoto-Nagai disk mode
//
// load me into exp/include/
//  this will compile IF listed as exp/include/exponential.h
// more clever work is needed to make this work generally... a typical
//  3dimensional disc?
//

#ifndef _MN_h
#define _MN_h

const char rcsid_mn[] = "$Id$";

#include <massmodel.h>

#include <numerical.h>
#include <gaussQ.h>
#include <interp.h>


//class MNDisk : public AxiSymModel
class ExponentialDisk : public AxiSymModel
{
private:
  
  double a;
  double h;
  double m;
  double rmin;
  double rmax;
  double den0;
  Linear1d massfunc, sdens;

public:

  
  
  //MNDisk(double RSCALE=1.0, double HSCALE=0.1, double RMAX=20.0, double DMASS=1.0) 
  ExponentialDisk(double RSCALE=1.0, double HSCALE=0.1, double RMAX=20.0, double DMASS=1.0) 
  { 
    a       = RSCALE;
    h       = HSCALE;
    rmin    = 1.0e-8;
    rmax    = RMAX;
    m       = DMASS;
    den0    = m*0.5/M_PI/a/a;
    dim     = 3;
    ModelID = "MN3Disk"; 

    dist_defined = false;

    
    tabulate_deprojection(HSCALE/RSCALE,1.0,200,1000,massfunc,sdens);
    
  }

  double disk_density(const double r, double z) {

    double Z2 = z*z + h*h;
    double Z  = sqrt(Z2);
    double Q2 = (a + Z)*(a + Z);
    return 0.25*h*h*m/M_PI*(a*r*r + (a + 3.0*Z)*Q2)/( pow(r*r + Q2, 2.5) * Z*Z2 );
}


  // Required member functions

  // 2d mass
  //double get_mass(const double r) { return m*(1.0 - exp(-r/a)*(1.0+r/a)); }
  double get_mass(const double r) { return m*(1.0 - massfunc.eval(r)); }

  // surface density
  //double get_density(const double r) { return den0*exp(-r/a); }
  double get_density(const double r) { return sdens.eval(r); }


  
  double get_pot(const double r) {
    return 0;
  }

  void test_threedee(void) { cout << "This is a 3D Miyamoto-Nagai Disk."; };							

  double get_dpot(const double r) {
    return 0;
  }

  double get_dpot2(const double r) {
    return 0;
  }

  void get_pot_dpot(const double r, double &ur, double &dur) {
    ur = 0;
    dur = 0;
  }
  
  // Additional member functions

  double get_min_radius(void) { return rmin; }
  double get_max_radius(void) { return rmax; }
  double get_ascale(void) { return a; }
  double get_hscale(void) { return h; }  
  double get_mass(void) { return m; }

  double distf(double E, double L) {
    bomb("Dist fct not defined!");
    return 0.0;
  }

  double dfde(double E, double L) {
    bomb("Dist fct not defined!");
    return 0.0;
  }

  double dfdl(double E, double L) {
    bomb("Dist fct not defined!");
    return 0.0;
  }
  
  double d2fde2(double E, double L) {
    bomb("Dist fct not defined!");
    return 0.0;
  }





  void tabulate_deprojection(double H, double Rf, int NUMR, int NINT,
			     Linear1d &massRg, Linear1d &rsurf)
{
  LegeQuad lq(NINT);
  Linear1d densRg;//, massRg;

  std::vector<double> rr(NUMR), rl(NUMR), sigI(NUMR), rhoI(NUMR, 0.0);

  double Rmin = log(rmin);
  double Rmax = log(rmax);

  double dr = (Rmax - Rmin)/(NUMR-1);

  // Compute surface mass density, Sigma(R)
  //
  for (int i=0; i<NUMR; i++) {
    double r = Rmin + dr*i;

    // Save for finite difference
    //
    rl[i] = r;
    r = exp(r);
    rr[i] = r;

    // Interval by Legendre
    //
    sigI[i] = 0.0;
    for (int n=1; n<=NINT; n++) {
      double y   = lq.knot(n);
      double y12 = 1.0 - y*y;
      double z   = y/sqrt(y12)*H;

      sigI[i] += lq.weight(n)*2.0*H*pow(y12, -1.5)*disk_density(r*Rf, z);
    }
  }

  Linear1d surf(rl, sigI);
  //Linear1d rsurf(rr, sigI);
  rsurf = Linear1d(rr, sigI);
 
  // Now, compute Abel inversion integral
  //
  for (int i=0; i<NUMR; i++) {
    double r = rr[i];

    // Interval by Legendre
    //
    rhoI[i] = 0.0;
    for (int n=1; n<=NINT; n++) {
      double x   = lq.knot(n);
      double x12 = 1.0 - x*x;
      double z   = x/sqrt(x12);
      double R   = sqrt(z*z + r*r);
      double lR  = log(R);

      rhoI[i]   += lq.weight(n)*2.0*pow(x12, -1.5)*surf.eval(lR);
    }
  }

  std::vector<double> rho(NUMR), mass(NUMR);

  Linear1d intgr(rl, rhoI);

  for (int i=0; i<NUMR; i++)
    rho[i] = -intgr.deriv(rl[i])/(2.0*M_PI*rr[i]*rr[i]);

  mass[0] = 0.0;
  for (int i=1; i<NUMR; i++) {
    double rlst = rr[i-1], rcur = rr[i];
    mass[i] = mass[i-1] + 2.0*M_PI*(rlst*rlst*rho[i-1] + rcur*rcur*rho[i])*(rcur - rlst);
  }

  // Debug
  //
    std::ostringstream outf; outf << "deproject_sl.1.0";
    std::ofstream out(outf.str());
    if (out) {
      for (int i=0; i<NUMR; i++)
	out << std::setw(18) << rl[i]
	    << std::setw(18) << rr[i]
	    << std::setw(18) << rho[i]
	    << std::setw(18) << mass[i]
	    << std::endl;
    }

  
  // Finalize
  //
  densRg = Linear1d(rl, rho);
  massRg = Linear1d(rl, mass);
}

};

#endif
