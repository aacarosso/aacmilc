// ------------------------------------------------------------------
// Flows link to tf given some U(ti) with offset off.
// Saves the W's obtained in the last step.

// improved scheme for variable epsilon; begins with stepsize = eps

#include "wflow_includes.h"
//-------------------------------------------------------------------

//-------------------------------------------------------------------
void wflow_imp(double eps, field_offset off, Real ti, Real tf) {
  register int dir, i;
	register site *s;
  int last=0, step;
	Real k, l=eps/epsilon, dl;
  Real t=ti, cut = 1e-7, eps_max = 0.1;
  double E, old_value=0, new_value=0, der_value, check, dS, eta, slope_E, slope_td, slope_topo;
	double E0, td0, topo0, Ek, tdk, topok, old_valuek, new_valuek, der_valuek, checkk, tk;
  double ssplaq, stplaq, td, Ps1, Pt1, Ps2, Pt2, topo, slope_newval;
	complex tc;
  su3_matrix t_mat, *S[4];
  anti_hermitmat *A[4];

  // Allocate fields used by integrator
  for (dir = 0; dir < 4; dir++){
    S[dir] = (su3_matrix *)malloc(sites_on_node * sizeof(su3_matrix));
    A[dir] = (anti_hermitmat *)malloc(sites_on_node * sizeof(anti_hermitmat));
  }
  
  // set link = U(ti)
	for (dir = 0; dir < 4; dir++){
		FORALLSITES(i,s){
			su3mat_copy((su3_matrix *)F_PT(s, off+dir*sizeof(su3_matrix)), &(s->link[dir]));
		}
	}

	// Must calculate t^2*E at ti (for interpolation) !

	// Finds 8F_munu = sum_{clover} (U - Udag)
  // Subtracts the (lattice artifact?) trace at each lattice site
	make_field_strength(F_OFFSET(link), F_OFFSET(fieldstrength));

	// The variables necessary for scale setting: used if tmax==0 only
	E = 0;
	FORALLSITES(i, s) {
		for (dir = 0; dir < 6; dir++) {
			mult_su3_nn(&(s->fieldstrength[dir]),
	                &(s->fieldstrength[dir]), &t_mat);
      tc = trace_su3(&t_mat);
      E -= (double)tc.real;
    }
  }
  g_doublesum(&E);
  E /= (volume * 64); // Normalization factor of 1/8 for each F_munu
  new_value = t * t * E;
  der_value = fabs(t) * (new_value - old_value) / fabs(eps); // nonsensical
	
	// Might as well extract topology
	topo = 0;
  FORALLSITES(i, s) {
    mult_su3_nn(&(s->fieldstrength[0]), // XYZT
                &(s->fieldstrength[5]), &t_mat);
    tc = trace_su3(&t_mat);
    topo -= (double)tc.real;

    mult_su3_nn(&(s->fieldstrength[3]), // XTYZ
                &(s->fieldstrength[2]), &t_mat);
    tc = trace_su3(&t_mat);
    topo -= (double)tc.real;

    mult_su3_na(&(s->fieldstrength[1]), // XZTY;  TY=(YT)^dag
                &(s->fieldstrength[4]), &t_mat);
    tc = trace_su3(&t_mat);
    topo -= (double)tc.real;
  }
  g_doublesum(&topo);
  // Same normalization
  topo /= (volume * 64 * 0.02533029591058444286); // 1 / (volume / 4pi^2)

  // Check with plaquette
  d_plaquette(&ssplaq, &stplaq);
  td = (ssplaq + stplaq) / 2;
  check = 12 * t * t * (3 - td);
  
	// Don't want to print these, assuming WFLOW at ti has already been calculated elsewhere.
	//node0_printf("WFLOW %g %g %g %g %g %g %g\n", t, td, E, new_value, der_value, check, topo);

	node0_printf("BEGIN WILSON FLOW (IMP) tf = %g  ti = %g\n",tf,ti);
	
  d_plaquette(&Ps1, &Pt1);

	// Wilson flow!
  for (step = 0; fabs(t) < fabs(tf) - fabs(epsilon)/2; step++){
		// Save W?
		if ( pow((t+eps-tf),2) < cut){ last = 1;}
    
		fstout_step_rk(S, A, eps, last);
		t += eps;

		// previous data for the interpolation
		E0 = E;
		td0 = td;
		topo0 = topo;

		// Finds 8F_munu = sum_{clover} (U - Udag)
    // Subtracts the (lattice artifact?) trace at each lattice site
		make_field_strength(F_OFFSET(link), F_OFFSET(fieldstrength));

		// The variables necessary for scale setting: used if tmax==0 only
		old_value = new_value;
		E = 0;
		FORALLSITES(i, s) {
			for (dir = 0; dir < 6; dir++) {
				mult_su3_nn(&(s->fieldstrength[dir]),
		                &(s->fieldstrength[dir]), &t_mat);
        tc = trace_su3(&t_mat);
        E -= (double)tc.real;
      }
    }
    g_doublesum(&E);
    E /= (volume * 64); // Normalization factor of 1/8 for each F_munu
    new_value = t * t * E;
    der_value = fabs(t) * (new_value - old_value) / fabs(eps);
    // Any negative signs in t and epsilon should cancel out anyway...
		
		// Might as well extract topology
		topo = 0;
    FORALLSITES(i, s) {
      mult_su3_nn(&(s->fieldstrength[0]), // XYZT
                  &(s->fieldstrength[5]), &t_mat);
      tc = trace_su3(&t_mat);
      topo -= (double)tc.real;

      mult_su3_nn(&(s->fieldstrength[3]), // XTYZ
                  &(s->fieldstrength[2]), &t_mat);
      tc = trace_su3(&t_mat);
      topo -= (double)tc.real;

      mult_su3_na(&(s->fieldstrength[1]), // XZTY;  TY=(YT)^dag
                  &(s->fieldstrength[4]), &t_mat);
      tc = trace_su3(&t_mat);
      topo -= (double)tc.real;
    }
    g_doublesum(&topo);
    // Same normalization
    topo /= (volume * 64 * 0.02533029591058444286); // 1 / (volume / 4pi^2)

    // Check with plaquette
    d_plaquette(&ssplaq, &stplaq);
    td = (ssplaq + stplaq) / 2;
    check = 12 * t * t * (3 - td);

		// Interpolate between eps steps
		if ( eps > epsilon ){
			//slope_E = (E - E0)/eps;
			slope_newval = (new_value - old_value)/eps;
			slope_td = (td - td0)/eps;
			slope_topo = (topo - topo0)/eps;
			tk = t - eps;
			new_valuek = tk*tk*E0;
			for (k = 1; k < l - cut; k++){
			  old_valuek = new_valuek;
				tk = t - eps + k*epsilon;
				//Ek = slope_E*(k*epsilon) + E0;
				tdk = slope_td*(k*epsilon) + td0;
				topok = slope_topo*(k*epsilon) + topo0;
				checkk = 12*tk*tk*(3 - tdk);
				new_valuek = slope_newval*k*epsilon + old_value;
				//new_valuek = tk*tk*Ek;
				Ek = new_valuek/(tk*tk);
				der_valuek = tk*(new_valuek - old_valuek)/epsilon;
				node0_printf("WFLOW %g %g %g %g %g %g %g (interp)\n", tk, tdk, Ek, new_valuek, 
				              der_valuek, checkk, topok);
			}
		}

    node0_printf("WFLOW %g %g %g %g %g %g %g\n", t, td, E, new_value, der_value, check, topo);
		
		// Setting epsilon
		d_plaquette(&Ps2, &Pt2);
		dS = fabs((Ps2 + Pt2)/2 - (Ps1 + Pt1)/2);
		if (step == 0){
			eta = dS*eps;
			//node0_printf(" eta %g\n", eta);
		}
		eps = eta/dS;
		//node0_printf(" eta/dS %g", eps);
		dl = floor(100*eps)/(epsilon*100) - l;
		if (eps < (l+0.99)*epsilon){
			eps = l*epsilon;
		}
		else {
			//node0_printf(" l += %g\n", dl);
			l += 1;
			eps = l*epsilon;
		}
		//node0_printf(" eps_new %g l_new %g\n",eps,l);

		if (eps >= eps_max){
			eps = eps_max;
			l = eps_max/epsilon;
			//node0_printf(" eps_max %g l %g\n", eps_max, l);
		}

		Ps1 = Ps2;
		Pt1 = Pt2;
	
		// Don't shoot past tmax
		if (t + eps > tf){
			eps = tf - t;
			l = eps/epsilon;
			//node0_printf(" eps %g l %g\n", eps, l);
		}

    last = 0;
  }
  node0_printf("END WILSON FLOW\n");

  for (dir = 0; dir < 4; dir++) {
    free(S[dir]);
    free(A[dir]);
  }
	return;
}

















