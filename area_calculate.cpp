/*
Working compile string:

clang++ --std=c++11 -Ofast -lboost_program_options -lCGAL -lgmp mp2.cpp -o mp2

gcc --std=c++11 -I /usr/local/boost_1_61_0/include -L/usr/local/boost_1_61_0/lib -I /opt/local/include  -lboost_program_options -lCGAL -lgmp mp2.cpp -o mp2


clang++ --std=c++11 -I /usr/local/boost_1_61_0/include -L/usr/local/boost_1_61_0/lib -I /opt/local/include  -lboost_program_options -lCGAL -lgmp mp2.cpp -o mp2

// specify all libraries and includes
clang++ -Ofast --std=c++11 -v -I/usr/local/boost_1_61_0/include -L/usr/local/boost_1_61_0/lib -L/usr/local/CGAL-4.9.1/lib -I/usr/local/CGAL-4.9.1/include -L/usr/local/lib -I/usr/local/include -lCGAL -lgmp -lboost_program_options mp2.cpp -o mp2_tri

// on cluster, simple compile string works:
g++ -Ofast --std=c++11 -v -lCGAL -lgmp -lboost_program_options mp2_small.cpp -o mp2


*/

#include <boost/program_options.hpp>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>

#include <algorithm>
#include <fstream>
#include <string>
#include <array>
#include <cmath>


typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Delaunay_triangulation_2<K> DT;
typedef DT::Point Point;

namespace po = boost::program_options;

int main(int argc, char** argv)
{

  double omega, ratio, dt;
  std::string input, output, outputt;

  po::options_description desc(
"This routine uses the CGAL Delaunay triangulation to find the projected x-y\n"
"area in each orbit from the R-V grid. The high density regions are isolated\n"
"by removing the simplices with extreme axis ratio.\n\nCommand line options");

  desc.add_options()
    ("input,i",      po::value<std::string>(&input)->default_value("sample_orbits.txt"), "input data")
    ("output,o",     po::value<std::string>(&output)->default_value("processed.txt"), "output grid")
    ("outputt,ot",     po::value<std::string>(&outputt)->default_value("processed_triangles.txt"), "output triangles grid")
    ("omega,O",      po::value<double>(&omega)->default_value(37.5), "pattern speed")
    ("ratio,r",      po::value<double>(&ratio)->default_value(15.0), "maximum axis ratio")
    ("dt,t",         po::value<double>(&dt)->default_value(0.0003), "orbit time step")
    ("help,h",       "this help message");

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);
  
  if (vm.count("help")) {
    std::cout << desc << std::endl;
    return 0;
  }

  /*
    Structure of Mike's file

    #samples, R_init, Vtan_init, dt, X series (#samples long), Y series, Z
    series, VX series, VY series, VZ series
  */
  std::ifstream in(input.c_str());
  std::ofstream out(output.c_str());
  std::ofstream outt(outputt.c_str());
  
  int sz = 524288;
  char buf[sz];

  in.getline(buf, sz);

  while (in.good()) {

    std::istringstream stIn(buf);
    
    int nsample;
    double R_init, V_tan_init;

    stIn >> nsample;
    stIn >> R_init;
    stIn >> V_tan_init;
    stIn >> dt;

    std::vector<double> X(nsample), Y(nsample);
    for (int i=0; i<nsample; i++) stIn >> X[i];
    for (int i=0; i<nsample; i++) stIn >> Y[i];

    for (int i=0; i<nsample; i++) {
      double angle = omega*dt*i;
      double cosT = cos(angle);
      double sinT = sin(angle);
      double x = X[i], y = Y[i];
      X[i] =  cosT*x + sinT*y;
      Y[i] = -sinT*x + cosT*y;
     }

    std::vector<Point> points(nsample);
    for (int i=0; i<nsample; i++) points[i] = Point(X[i], Y[i]);

    // compute the Delaunay Triangulation, in two dimensions.
    DT dt(points.begin(), points.end());

    if (dt.is_valid()) {
     
      DT::Finite_faces_iterator face;
      typedef std::array<double, 3> elem;
      std::vector<elem> A;

      for (face=dt.finite_faces_begin(); face!= dt.finite_faces_end(); face++) {
	std::array<DT::Point, 3> pts;
	for (int i=0; i<3; i++) pts[i] = face->vertex(i)->point();
	
	double x1 = pts[0].cartesian(0), y1 = pts[0].cartesian(1);
	double x2 = pts[1].cartesian(0), y2 = pts[1].cartesian(1);
	double x3 = pts[2].cartesian(0), y3 = pts[2].cartesian(1);
	
	std::array<double, 3> d {
	    sqrt((x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2)),
	    sqrt((x1 - x3)*(x1 - x3) + (y1 - y3)*(y1 - y3)),
	    sqrt((x2 - x3)*(x2 - x3) + (y2 - y3)*(y2 - y3))};

	 std::sort(d.begin(), d.end());

	 double area = 0.5*((x2 - x1)*(y3 - y1) - (x3 - x1)*(y2 - y1));
	 
	 A.push_back({d[0], d[2], area});

      }



      std::sort(A.begin(), A.end());
      int indx = std::floor(0.5*A.size());
      double med = A[indx][0];

      double area = 0.0;
      for (auto e : A) {
	if (e[1]/med < ratio) area += e[2];
      }


  
      

      out << std::setw(18) << R_init
	  << std::setw(18) << V_tan_init
	  << std::setw(18) << area
	  << std::endl;

            // now loop back through the faces and record the valid ones



      outt << std::setw(18) << R_init
	   << std::setw(18) << V_tan_init;

      
      for (face=dt.finite_faces_begin(); face!= dt.finite_faces_end(); face++) {
	std::array<DT::Point, 3> pts;
	for (int i=0; i<3; i++) pts[i] = face->vertex(i)->point();
	
	double x1 = pts[0].cartesian(0), y1 = pts[0].cartesian(1);
	double x2 = pts[1].cartesian(0), y2 = pts[1].cartesian(1);
	double x3 = pts[2].cartesian(0), y3 = pts[2].cartesian(1);
	
	std::array<double, 3> d {
	    sqrt((x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2)),
	    sqrt((x1 - x3)*(x1 - x3) + (y1 - y3)*(y1 - y3)),
	    sqrt((x2 - x3)*(x2 - x3) + (y2 - y3)*(y2 - y3))};

	 std::sort(d.begin(), d.end());

	 // needs an if here for the ratio
	 if (d[2]/med < ratio) {
	 	outt  << std::setw(18) << x1 << std::setw(18) << x2 << std::setw(18) << x3
		      << std::setw(18) << y1 << std::setw(18) << y2 << std::setw(18) << y3;
	 }
      }

      	 outt << std::endl;


      
    }
       
    in.getline(buf, sz);
  }
  
  return 0;
}
