// This may look like C code, but it is really -*- C++ -*-

//
// 3D Exponential disk
//

#ifndef _Expon_h
#define _Expon_h

const char rcsid_expon[] = "$Id$";

#include <massmodel.h>

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

};



#endif
