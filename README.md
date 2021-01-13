# python analysis of exp files (cpp routines)


## Version 0.1
A collection of C++ routines that will eventually be turned into something real, hopefully. These are currently mostly intended as add-ons to EXP; some codes are available as standalone implementations.

Currently, the included toolbox has
1. analysis/
   1. area_calculate: Delaunay triangulation (standalone)
   2. forcetrace: template for slicing across EXP dumps efficiently. (requires EXP)
   3. spinparam: compute the spin paramater of a halo. (requires EXP)
2. io/ (all require EXP)
3. manga_data/ (all standalone)
   1. readmanga: read in custom packaged manga data for C++ analysis.
4. ic_generation/ (all require EXP)

   2. basicics
   
	   1. DiskHalo5.cc

	   2. DiskHalo5.h
	   
   3. bulgeics:
	   
	   1. DiskHalo6.cc

	   2. DiskHalo5.h
	   
   4. flexics:
	   
   5. table_disk:
	   
   6. exponential3.h

   7. twopower: build simple two-power density models

   8. cylcache: generate empirical orthogonal functions as per Weinberg 1999

5. include/
   1. DiskModels.H: upgraded options for disc models (requires EXP)
