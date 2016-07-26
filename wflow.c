// ------------------------------------------------------------------
// Flows link to tf given some U(ti) with offset "off".
// If argument savelink = 1, it saves nsave configurations equally spaced.
// Always saves the W's obtained in the last step, for fermion flow.
// nsave must be consistent with the number of fields link0,1,2,3,... in lattice.h

#include "wflow_includes.h"
//-------------------------------------------------------------------

//-------------------------------------------------------------------
void wflow(field_offset off, Real ti, Real tf, int savelink) {
  register int dir, i;
	register site *s;
  int last=0, k=1, step;
  double E, old_value, new_value=0, der_value, check, eps = 1*epsilon;
  Real t=0, cut = 1e-7, N = floor(tmax/(eps*nsave));
  double ssplaq, stplaq, td, topo;
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
  
  node0_printf("BEGIN WILSON FLOW (FIX) tf = %g  ti = %g\n", tf, ti);
	
	if ((savelink == 1) && (ti < cut)){
		node0_printf("saving t = 0 configuration (1 of %d)\n", nsave);
		for (dir = 0; dir < 4; dir ++){
			FORALLSITES(i,s){
	  	  su3mat_copy(&(s->link[dir]),
					(su3_matrix *)F_PT(s,F_OFFSET(link0)+dir*sizeof(su3_matrix)));	
			}
  	}
  }
  
	// Wilson flow!
	t = ti;
  for (step = 0; fabs(t) < fabs(tf) - fabs(eps)/2; step++){
	  // save W?
    if ( pow((t+eps-tf),2) < cut){ last = 1;}
    
		fstout_step_rk(S, A, eps, last);
		t += eps;

		// save link?
		if ((savelink == 1) && (k < nsave) && (pow(t/eps - N*k,2) < cut )){
		  node0_printf("saving t = %g configuration (%d of %d)\n", t, k+1,nsave);
			for (dir = 0; dir < 4; dir ++){
			  FORALLSITES(i,s){
			    su3mat_copy(&(s->link[dir]),
						(su3_matrix *)F_PT(s,F_OFFSET(link0)+(4*k+dir)*sizeof(su3_matrix)));
				}
			}
			k += 1;
		}

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
    // Any negative signs in t and eps should cancel out anyway...

    //Might as well extract topology
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
		topo *= 0.02533029591058444286/64;
    //topo /= (volume * 64 * 0.02533029591058444286); // 1 / (volume / 4pi^2)

    // Check with plaquette
    d_plaquette(&ssplaq, &stplaq);
    td = (ssplaq + stplaq) / 2;
    check = 12 * t * t * (3 - td);
    node0_printf("WFLOW %g %g %g %.10g %g %g %g\n", t, td, E, new_value, der_value, check, topo);
    last = 0;
  }
  node0_printf("END WILSON FLOW\n");

  for (dir = 0; dir < 4; dir++) {
    free(S[dir]);
    free(A[dir]);
  }
}

















