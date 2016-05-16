// Function which measures data for gauge field flow

#include "wflow_includes.h"

void meas_flowdata(Real t, Real eps) {
	register int i, dir;
	register site *s;
	double E, old_value, new_value=0, der_value, check;
	double ssplaq, stplaq, td, Ps1, Pt1, Ps2, Pt2, topo;
	complex tc;
	su3_matrix t_mat;
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
	node0_printf("WFLOW %g %g %g %g %g %g %g\n", t, td, E, new_value, der_value, check, topo);
}
