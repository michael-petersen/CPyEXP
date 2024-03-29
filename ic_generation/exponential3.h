// This may look like C code, but it is really -*- C++ -*-

//
// 3D Exponential disk
//

#ifndef _Expon_h
#define _Expon_h

const char rcsid_expon[] = "$Id$";

#include <massmodel.h>

#include <numerical.h>
#include <gaussQ.h>
#include <interp.h>


extern "C" double i0(double);
extern "C" double i1(double);
extern "C" double k0(double);
extern "C" double k1(double);

class ExponentialDisk : public AxiSymModel
{
private:
  
  double a;
  double h;
  double m;
  double rmin;
  double rmax;
  double den0;

public:

  ExponentialDisk(double RSCALE=1.0, double HSCALE=0.1, double RMAX=20.0, double DMASS=1.0) 
  { 
    a       = RSCALE;
    h       = HSCALE;
    rmin    = 1.0e-8;
    rmax    = RMAX;
    m       = DMASS;
    den0    = m*0.5/M_PI/a/a;
    dim     = 3;
    ModelID = "Exponential3Disk"; 

    dist_defined = false;

    tabulate_deprojection(HSCALE,RSCALE,200,1000);
    
  }

  double disk_density(const double r, double z) {

    double q = 1.0/cosh(z/h);

    return m*den0*exp(-r/a)*q*q*0.5/h;
}


  // Required member functions

  // 2d mass
  double get_mass(const double r) { return m*(1.0 - exp(-r/a)*(1.0+r/a)); }

  // surface density
  double get_density(const double r) { return den0*exp(-r/a); }


  
  double get_pot(const double r) {
    double y=0.5*r/a;
    return -M_PI*den0*r*(i0(y)*k1(y) - i1(y)*k0(y));
  }

  void test_threedee(void) { cout << "This is a 3D exponential Disk."; };							

  double get_dpot(const double r) {
    double y=0.5*r/a;
    return 2.0*M_PI*den0*y*(i0(y)*k0(y) - i1(y)*k1(y));
  }

  double get_dpot2(const double r) {
    double y=0.5*r/a;
    return M_PI*den0/a*(i0(y)*k0(y) + i1(y)*k1(y) - 
			2.0*y*(i1(y)*k0(y) - i0(y)*k1(y)));
  }

  void get_pot_dpot(const double r, double &ur, double &dur) {
    double y=0.5*r/a;
    double I0=i0(y), I1=i1(y), K0=k0(y), K1=k1(y);
    ur = -M_PI*den0*r*(I0*K1 - I1*K0);
    dur = 2.0*M_PI*den0*y*(I0*K0 - I1*K1);
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





void tabulate_deprojection(double H, double Rf, int NUMR, int NINT)
{
  LegeQuad lq(NINT);
  Linear1d densRg, massRg;

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

      sigI[i] += lq.weight(n)*2.0*H*pow(y12, -1.5)*disk_density(r, z);
    }
  }

  Linear1d surf(rl, sigI);
  Linear1d rsurf(rr, sigI);
  
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
    std::ostringstream outf; outf << "deproject_sl.exp.0";
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
