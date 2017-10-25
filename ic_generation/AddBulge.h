#ifndef _AddBulge_H
#define _AddBulge_H

#include <vector>

#include <massmodel.h>

class AddBulge
{
 private:
  vector<double> r, d, m, p;
  SphericalModelTable *mod;

 public:
  static int number;			// 4000
  static double Rmin;			// 1.0e-3
  static bool use_mpi;			// false
  static bool logarithmic;		// false

  AddBulge(AxiSymModel* halo, AxiSymModel* bulge);
  ~AddBulge();

  SphericalModelTable *get_model() { return mod; }
};


#endif
