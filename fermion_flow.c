// Routine for fixed epsilon fermion adjoint flow using checkpoints.
// Begins at tmax and flows to 0.

#include "wflow_includes.h"

void fermion_flow()
{
	register int i;
	register site *s;	
	int istep, n, dir, ksi, k = nsave-1, yep;
	Real t = tmax, ti, N = floor(tmax/(epsilon*nsave)), cut = 1e-7, counter = 0;
	double flowtime = 0, tvar, tvar1, ttime = 0;
	complex dot1, dot2;
	node0_printf("floor = %g\n",N);

	node0_printf("ksi %g ",t);
  for (n=0; n<npbp; n++){
    grsourcen(EVENANDODD);
    ksi = F_OFFSET(ksi1) + n*sizeof(su3_vector);
    copy_latvec(F_OFFSET(g_rand), ksi, EVENANDODD);
    dot1 = dot_su3_latvec(ksi, ksi, EVENANDODD);
    node0_printf("%.10g  ", dot1.real);
  }
	//tvar1 = dclock();
  //yep = fmeas_link(F_OFFSET(chi), F_OFFSET(link), F_OFFSET(psi), mass);
  //node0_printf("fmeas_link time = %f\n",dclock()-ttime);
  //ttime += dclock() - tvar1;

  node0_printf("\n");

  // Adjoint flow loop 

	node0_printf("\nBEGINNING FERMION ADJOINT FLOW\n\n");
	for (istep = 0; t > cut; istep++){
		// reset lattice
		for (dir = 0; dir < 4; dir++){
      FORALLSITES(i,s){
        su3mat_copy(&(s->link0[dir]),&(s->link[dir]));
      }
    }
    // flow gauge fields to t, obtain W's
    if (pow(t/epsilon - k*N,2) < cut ) k -= 1;
		ti = k*N*epsilon;
		tvar = dclock();
 		wflow(F_OFFSET(link0)+4*k*sizeof(su3_matrix), ti, t, 0);
		flowtime += dclock() - tvar;
    counter += (t-ti)/epsilon;
    node0_printf("ksi %g ", t-epsilon);
    // loop over npbp fields
    for (n=0; n<npbp; n++){
      ksi = F_OFFSET(ksi1) + n*sizeof(su3_vector);
      copy_latvec(ksi, F_OFFSET(lambda3), EVENANDODD);
      fermion_adjointstep(F_OFFSET(lambda3), epsilon);
      copy_latvec(F_OFFSET(lambda0), ksi, EVENANDODD);
      dot1 = dot_su3_latvec(ksi, ksi, EVENANDODD);
      node0_printf("%.10g  ", dot1.real);
    }
    node0_printf("\n");
    t -= epsilon;

		// Measure pbp
    tvar1 = dclock();
    yep = fmeas_link(F_OFFSET(chi),F_OFFSET(W0), F_OFFSET(psi), mass);
    ttime += dclock() - tvar1;

/*  // sequence test
    copy_latvec(F_OFFSET(g_rand),F_OFFSET(chi0),EVENANDODD);
    fermion_forwardstep(F_OFFSET(chi0));
    dot1 = dot_su3_latvec(F_OFFSET(lambda0),F_OFFSET(chi0),EVENANDODD);
    dot2 = dot_su3_latvec(F_OFFSET(lambda3),F_OFFSET(chi3),EVENANDODD);
    node0_printf("lambda0.chi0: %g  %g\n",dot1.real,dot1.imag);
    node0_printf("lambda3.chi3: %g  %g\n",dot2.real,dot2.imag);

    clear_latvec(F_OFFSET(lambda1), EVENANDODD);
    clear_latvec(F_OFFSET(lambda2), EVENANDODD);
    clear_latvec(F_OFFSET(lambda3), EVENANDODD);

    // forward step lambda0 to check with lambda3
    fermion_forwardstep(F_OFFSET(lambda0));

    dot1 = dot_su3_latvec(F_OFFSET(lambda3),F_OFFSET(lambda3),EVENANDODD);
    node0_printf("forward lambda3^2: %g\n",dot1.real);
*/
    node0_printf("\n\n");
  }

  node0_printf("\nENDING FERMION ADJOINT FLOW \n\n");

  node0_printf("total flow time = %f\n",flowtime);
  node0_printf("meas_link total time = %f\n",ttime);
  node0_printf("counter = %g\n",counter);
	return;
}
