
void DiskHalo::
table_double_disk(vector<Particle>& part1, vector<Particle>& part2)
{
  if (disktableP != NULL) return;

  disktableP = new Matrix [NDP];
  disktableN = new Matrix [NDP];
  for (int i=0; i<NDP; i++) {
    disktableP[i].setsize(0, NDR-1, 0, NDZ-1);
    disktableN[i].setsize(0, NDR-1, 0, NDZ-1);
  }
    
  epitable.setsize(0, NDP-1, 0, NDR-1);
  dv2table.setsize(0, NDP-1, 0, NDR-1);
  asytable.setsize(0, NDP-1, 0, NDR-1);

  dP = 2.0*M_PI/NDP;

  // Sift through particles to find the maximum extent
  double maxr1 = 0.0, maxz1 = 0.0;
  double maxr =  0.0, maxz  = 0.0;
  for (auto &p : part1) {
    maxr1 = max<double>(maxr1, sqrt(p.pos[0]*p.pos[0] + p.pos[1]*p.pos[1]));
    maxz1 = max<double>(maxz1, fabs(p.pos[2]));
  }

    for (auto &p : part2) {
    maxr1 = max<double>(maxr1, sqrt(p.pos[0]*p.pos[0] + p.pos[1]*p.pos[1]));
    maxz1 = max<double>(maxz1, fabs(p.pos[2]));
  }

    maxz1 = max<double>(scaleheight1*SHFACTOR, scaleheight2*SHFACTOR, maxz1);
  

  MPI_Allreduce(&maxr1, &maxr, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(&maxz1, &maxz, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

  dR = (log(maxr) - log(RDMIN))/(NDR-1);
  dZ = maxz/(NDZ-1);

  if (myid==0) {
    cout << endl
	 << "Table disk epi parameters:" << endl
	 << "  RDMIN=" << RDMIN           << endl
	 << "   maxr=" << maxr            << endl
	 << "   maxz=" << maxz            << endl
	 << "   dR=" << dR                << endl
	 << "   dZ=" << dZ                << endl
	 << endl;
  }

  double dZ1 = maxz/(NDZ-1);

  double R, r, x, y, z, phi;
  double fr1, fz1, fp1, fr2, fz2, fp2, theta;
  double pot1, pot2;
  double dens, potl, dpr, dpt, dpp, dpdz;

				// Add no force if no component exists
  pot1 = fr1 = fz1 = fp1 = pot2 = fr2 = fz2 = fp2;
    dens = potl = dpr = dpt = dpp = 0.0;

  Vector workP (0, NDZ-1);
  Vector workN (0, NDZ-1);
  Vector workA (0, NDZ-1);
  Vector workZ (0, NDZ-1);

  Vector workR (0, NDR-1);
  Vector workE (0, NDR-1);
  Vector workE2(0, NDR-1);
  Vector workQ (0, NDR-1);
  Vector workQ2(0, NDR-1);

  Matrix workV (0, 4, 0, NDR-1);

				// For debugging
  Matrix workD (0, 7, 0, NDR-1);
  Matrix workDZ(0, 6, 0, NDZ-1);

				// Sum mass grid and make radial mesh
  // this already includes all the components naturally, so we'll just go with it
  unsigned nzcnt=0;
  vector<double> nrD(nh+1);
  for (nzero=0; nzero<=nh; nzero++) if (nhN[nzero] >= mh) break;

  if (nzero>nh) nzero=0;	// Not enough particles . . . 
  if (myid==0) cout << "Nzero=" << nzero << "/" << nh << endl;

  nzero = floor( (hDmin + nzero*dRh - log(RDMIN))/dR ) + 1;
  if (myid==0) cout << "Nzero=" << nzero << "/" << NDR << endl;

  for (int n=0; n<=nh; n++) nrD[n] = hDmin + dRh*n;
  for (int n=1; n<=nh; n++) nhD[n] += nhD[n-1];

				// Compute this table in parallel

  vector<int> ibeg(numprocs), iend(numprocs);
  int curid = -1;
  for (int i=0; i<numprocs; i++) {
    ibeg[i] = (i  )*NDP/numprocs;
    iend[i] = (i+1)*NDP/numprocs;
    if (curid<0 && iend[i]-ibeg[i]>0) curid = i;
    if (myid==0) {
      if (i==0) cout << endl << " *** Processor phi angles *** " << endl;
      cout << "# " << setw(3) << i << ": " 
	   << setw(10) << ibeg[i]
	   << setw(10) << iend[i]
	   << endl;
    }
  }

  Cheby1d *cheb = 0, *cheb2 = 0;


  for (int i=ibeg[myid]; i<iend[myid]; i++) {

    phi = dP*i;

    for (int j=0; j<NDR; j++) {

      R = RDMIN*exp(dR*j);
      x = R*cos(phi);
      y = R*sin(phi);
      workR[j] = log(RDMIN) + dR*j;

				// For epicylic frequency
				// 
      disk_eval(disk1, dmass1, R, 0.0, phi, pot1, fr1, fz1, fp1);
      disk_eval(disk2, dmass2, R, 0.0, phi, pot2, fr2, fz2, fp2);
				// Use monopole part of expansion here, only
      if (expandh)		//
	expandh->determine_fields_at_point(R, 0.5*M_PI, phi,
					   &dens, &potl, &dpr, &dpt, &dpp);

      
      workV[0][j] = log(RDMIN) + dR*j;
				// Use monopole approximation for dPhi/dr

      // might want to go back to this with the new double disk
      // workE[j] = odd2(workV[0][j], nrD, nhD, 1)/(R*R);
				// Use basis evaluation (dPhi/dr)
      workE[j]    = max<double>(-fr - fr2 + dpr, 1.0e-20);

      workV[1][j] = double_disk_surface_density(disk1, disk2, R);
				// Sigma(R)*dPhi/dr*R
      workV[2][j] = workV[1][j]*workE[j]*R;
				// [1/r dPhi/dr]^{1/2} = Omega
      workQ[j]    = sqrt(workE[j]/R);
				// r^2*[1/r dPhi/dr]^{1/2} = r^2 * Omega
      workQ2[j]   = workQ[j] * R*R;

      if (i==0) {
	workD[4][j] = -fr1 - fr2;
	workD[5][j] = dpr;
	workD[6][j] = potl;
	workD[7][j] = dens;
      }

      for (int k=0; k<NDZ; k++) {

	z = workZ[k] = dZ1*k;

	r = sqrt(R*R + z*z) + MINDOUBLE;

				// Positive
	disk_eval(disk1, dmass1, R, z, phi, pot1, fr1, fz1, fp1);
	disk_eval(disk2, dmass2, R, z, phi, pot2, fr2, fz2, fp2);
	theta = acos(z/(r + MINDOUBLE));

	if (expandh)
	  expandh->determine_fields_at_point(R, theta, phi,
					     &dens, &potl, &dpr, &dpt, &dpp);

	dpdz = -fz1 -fz2 + dpr*z/r + dpt*R*R/(r*r*r);

				// Debugging
	//disk_density(ExponentialDisk* disk, double R, double z, double scaleheight, double dmass)
	
	workDZ[0][k] = double_disk_density(disk1, disk2, R, z);
	workDZ[1][k] = -fz1 - fz2;
	workDZ[2][k] = dpr*z/r;
	workDZ[3][k] = dpt*R*R/(r*r*r);

	workP[k] = double_disk_density(disk1, disk2, R, z) * dpdz;

				// Negative
	z *= -1.0;

	disk_eval(disk1, dmass1, R, z, phi, pot1, fr1, fz1, fp1);
	disk_eval(disk2, dmass2, R, z, phi, pot2, fr2, fz2, fp2);
	theta = acos(z/(r + MINDOUBLE));

	if (expandh)
	  expandh->determine_fields_at_point(R, theta, phi,
					     &dens, &potl, &dpr, &dpt, &dpp);

	dpdz  = -fz1 - fz2 + dpr*z/r + dpt*R*R/(r*r*r);
	dpdz *= -1.0;		// No harm in changing sign

	workDZ[4][k] = -fz1 - fz2;
	workDZ[5][k] = dpr*z/r;
	workDZ[6][k] = dpt*R*R/(r*r*r);

	workN[k] = double_disk_density(disk1, disk2, R, z) * dpdz;
      }

				// Integrate positive
      if (use_spline) Splsum (workZ, workP, workA);
      else            Trapsum(workZ, workP, workA);

      for (int k=0; k<NDZ; k++)
	disktableP[i][j][k] = max<double>(workA[NDZ-1] - workA[k], MINDOUBLE);

      if (i==ibeg[myid] && j==0 && VFLAG & 4) {
	ostringstream ofile;
	ofile << "intgr_disk_P." << RUNTAG << ".d" << myid;
	ofstream dout(ofile.str().c_str());
	for (int k=0; k<NDZ; k++) {
	  dout << setw(15) << workZ[k] 
	       << setw(15) << workP[k] 
	       << setw(15) << workA[k];
	  for (int q=0; q<7; q++) dout << setw(15) << workDZ[q][k];
	  dout << endl;
	}
	dout.close();
      }
				// Integrate negative
      if (use_spline) Splsum (workZ, workN, workA);
      else            Trapsum(workZ, workN, workA);

      for (int k=0; k<NDZ; k++)
	disktableN[i][j][k] = max<double>(workA[NDZ-1] - workA[k], MINDOUBLE);

      if (i==ibeg[myid] && j==0 && VFLAG & 4) {
	ostringstream ofile;
	ofile << "intgr_disk_N." << RUNTAG << ".d" << myid;
	ofstream dout(ofile.str().c_str());
	for (int k=0; k<NDZ; k++) 
	  dout << setw(15) << workZ[k] 
	       << setw(15) << workN[k] 
	       << setw(15) << workA[k]
	       << "\n";
	dout.close();
      }

    }

    if (CHEBY) cheb = new Cheby1d(workR, workQ2, NCHEB);

				// Compute epicylic freqs
    for (int j=0; j<NDR; j++) {
      if (i==0) {
	if (CHEBY)
	  workD[0][j] += cheb->eval(workR[j]);
	else
	  workD[0][j] += workE[j];
      }

      if (CHEBY)
	epitable[i][j] = cheb->deriv(workR[j]);
      else
	epitable[i][j] = drv2(workR[j], workV[0], workQ2);

      if (i==0) workD[1][j] = epitable[0][j];
      epitable[i][j] *= 2.0*workQ[j]/exp(2.0*workR[j]);
      if (i==0) workD[2][j] = epitable[0][j];
      if (i==0) workD[3][j] = epitable[0][j];
    }
    
				// Cylindrical Jeans' equations
    for (int j=0; j<NDR; j++) {
      double vr   = 3.36*workV[1][j]*Q/sqrt(epitable[i][j]);
      workV[4][j] = log(workV[1][j]*vr*vr);
    }

    if (CHEBY) cheb2 = new Cheby1d(workV[0], workV[4], NCHEB);
  

    Trapsum(workV[0], workV[2], workV[3]);
    for (int j=0; j<NDR; j++) {
      dv2table[i][j] = (workV[3][NDR-1] - workV[3][j])/workV[1][j];
      if (CHEBY)
	asytable[i][j] = cheb2->deriv(workV[0][j]);
      else
	asytable[i][j] = drv2(workV[0][j], workV[0], workV[4], 1);
    }
  }

				// Update tables on all nodes
  for (int k=0; k<numprocs; k++) {
    for (int i=ibeg[k]; i<iend[k]; i++) {
      MPI_Bcast(&epitable[i][0], NDR, MPI_DOUBLE, k, MPI_COMM_WORLD);
      MPI_Bcast(&dv2table[i][0], NDR, MPI_DOUBLE, k, MPI_COMM_WORLD);
      MPI_Bcast(&asytable[i][0], NDR, MPI_DOUBLE, k, MPI_COMM_WORLD);
      for (int j=0; j<NDR; j++) {
	MPI_Bcast(&disktableP[i][j][0], NDZ, MPI_DOUBLE, k, MPI_COMM_WORLD);
	MPI_Bcast(&disktableN[i][j][0], NDZ, MPI_DOUBLE, k, MPI_COMM_WORLD);
      }
    }
  }

				// For debugging the solution
  if (myid==curid && expandh && VFLAG & 4) {
    ostringstream sout;
    sout << "ep_test." << RUNTAG;
    ofstream out(sout.str().c_str());
    out.setf(ios::scientific);
    out.precision(4);

    const double DR = 0.01;

    double r, r1, r2, deriv, deriv2, vrq0, vrq1, rho, lhs, rhs;

    for (int j=0; j<NDR; j++) {
      r = RDMIN*exp(dR*j);
      r1 = r*(1.0 + dR*DR);
      r2 = r*(1.0 - dR*DR);

      rho = halo->get_density(r);

      deriv = (get_disp(0.0, r1, 0.0)*halo->get_density(r1) - 
	       get_disp(0.0, r2, 0.0)*halo->get_density(r2) ) /	(r1 - r2);
      
      lhs = halo->get_mass(r);
      rhs = -r*r*deriv/rho;

      if (CHEBY)
	deriv2 = cheb->deriv(workQ2[j]);
      else
        deriv2 = drv2(workR[j], workR, workQ2);
	
      vrq0 = 3.36*(disk1->get_density(r)+disk2->get_density(r))*Q/epi(r, 0.0, 0.0);
      vrq1 = 3.36*(disk1->get_density(r)+disk2->get_density(r))*Q/sqrt(epitable[0][j]);

      out << setw(14) << r			// #1
	  << setw(14) << epitable[0][j]		// #2
	  << setw(14) << workR[j]		// #3
	  << setw(14) << workQ[j]		// #4
	  << setw(14) << workQ2[j]		// #5
	  << setw(14) << deriv2                 // #6
	  << setw(14) << vrq0                   // #7
	  << setw(14) << vrq1                   // #8
	  << setw(14) << v_circ(r, 0.0, 0.0)    // #9
	  << setw(14) << workD[0][j]		// #10  dV(tot)/dR
	  << setw(14) << workD[1][j]		// #11  d^2V(tot)/dlnR
	  << setw(14) << workD[2][j]		// #12  d^2V(tot)/dlnR + 3V(tot)
	  << setw(14) << workD[3][j]		// #13  kappa^2
	  << setw(14) << workD[4][j]		// #14  dV(disk)/dR
	  << setw(14) << workD[5][j]		// #15  dV(halo)/dR
	  << setw(14) << rho			// #16
	  << setw(14) << deriv			// #17
	  << setw(14) << lhs			// #18
	  << setw(14) << rhs			// #19
	  << setw(14) << lhs - rhs		// #20
	  << setw(14) << odd2(log(r), nrD, nhD) // #21  Enclosed mass
	  << setw(14) << epi(r, 0.0, 0.0)	// #22  Epi routine
	  << setw(14) << workD[6][j]            // #23  halo potential (from fit)
	  << setw(14) << workD[7][j]            // #24  halo density (from fit
	  << endl;
    }

    ostringstream sout2;
    sout2 << "epitable." << RUNTAG;
    ofstream dump(sout2.str().c_str());
    for (int i=0; i<NDP; i++) {
      phi = dP*i;
      for (int j=0; j<NDR; j++)
	dump << setw(18) << phi 
	     << setw(18) << RDMIN*exp(dR*j)
	     << setw(18) << epitable[i][j]
	     << endl;
      dump << endl;
    }
    dump.close();

    dump.open("dv2table.dump");
    for (int i=0; i<NDP; i++) {
      phi = dP*i;
      double c = cos(phi), s = sin(phi);
      for (int j=0; j<NDR; j++) {
	double rr = RDMIN*exp(dR*j);
	double vc = v_circ(rr*c, rr*s, 0.0);
	double vr = vr_disp2(rr*c, rr*s, 0.0);
	dump << setw(18) << phi 			// 1 Phi
	     << setw(18) << rr				// 2 R
	     << setw(18) << dv2table[i][j]		// 3 
	     << setw(18) << asytable[i][j]		// 4 d log(rho*vr^2)/ dlog(R)
	     << setw(18) << asytable[i][j]*vr/(vc*vc)	// 5 ac/vc^2
	     << setw(18) << workV[4][j]			// 6 log(rho*vr^2)
	     << setw(18) << vr				// 7 vr^2
	     << setw(18) << vc*vc;			// 8 vc^2
	if (CHEBY) dump << setw(18) << cheb2->eval(workV[0][j]);
	dump << endl;
      }
      dump << endl;
    }
    dump.close();
  }
    
  if (myid==curid && VFLAG & 4) {
    ostringstream sout;
    sout << "ep_disk." << RUNTAG;
    ofstream out(sout.str().c_str());
    out.setf(ios::scientific);
    out.precision(8);

    for (int i=0; i<=NDP; i++) {
      phi = dP*i;
      out << "# i=" << " phi=" << phi << ", " << phi*180.0/M_PI << endl;
      for (int j=0; j<NDR; j++) {
	R = RDMIN*exp(dR*j);
	x = R*cos(phi);
	y = R*sin(phi);
	out << setw(18) << x
	    << setw(18) << y
	    << setw(18) << epitable[i%NDP][j] << endl;
      }
      out << endl;
    }
    out << flush;
    out.close();
  }

  if (myid==curid && VFLAG & 4) {
    ostringstream sout;
    sout << "table_disk." << RUNTAG;
    ofstream out(sout.str().c_str());
    out.setf(ios::scientific);
    out.precision(8);
    
    for (int j=0; j<NDR; j++) {
      for (int k=0; k<NDZ; k++) {
	out << setw(18) << RDMIN*exp(dR*j)
	    << setw(18) << dZ*k
	    << setw(18) << disktableP[0][j][k]
	    << setw(18) << disktableN[0][j][k] 
	    << setw(18) << disk_density(RDMIN*exp(dR*j), dZ*k)
	    << endl;
      }
      out << endl;
    }
  }

  if (myid==0) cout << "[table] " << flush;

  delete cheb;
  delete cheb2;
}

