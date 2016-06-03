// Flows lambda0 from tf to ti knowing gauge field U(ti) with
// initial step size eps at ti (for gauge flow). 

// Uses wflow_imp.c with wflow.c for gauge flow

#include "wflow_includes.h"

void fermion_flow_chunk(field_offset off, double ti, double tf, double eps)
{
	register int i;
	register site *s;
	int istep, dir, ksi, n;
  double eps_max=0.1, flowtime=0, sum_eps, ttime=0, tvar, tvar1;
  Real t = tf, cut = 1e-7, counter = 0;
  complex dot1, dot2;

  node0_printf("\n");

  for (dir = 0; dir < 4; dir++){
    FORALLSITES(i,s){
      su3mat_copy((su3_matrix *)F_PT(s,off+dir*sizeof(su3_matrix)), &(s->link[dir]));
    }
  }

	double j = (tf-ti)/eps_max - 1;
  node0_printf("\nBEGIN FERMION CHUNK FLOW tf = %g ti = %g\n\n", tf, ti);
  for (istep = 0; t > ti+cut; istep++){
		
		if (fabs(t - ti - j*eps_max) < cut ){ j -= 1;}

    // flow gauge fields to t, obtain W's
    tvar = dclock();
    wflow_imp(eps_max, off, ti, ti + j*eps_max);
    //wflow_imp(eps, off, ti, t - epsilon);
    wflow(F_OFFSET(link), ti + j*eps_max, t, 0);
    //wflow(F_OFFSET(link), t - epsilon, t, 0);
    flowtime += dclock() - tvar;

    // loop over npbp fields
    node0_printf("ksi %g ", t-epsilon);
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
    tvar1 = dclock();
   // yep = fmeas_link(F_OFFSET(chi),F_OFFSET(W0), F_OFFSET(psi), mass);
    ttime += dclock() - tvar1;
    node0_printf("\n");
  }
  node0_printf("\nEND FERMION CHUNK FLOW\n");

  //node0_printf("total flow time = %f\n",flowtime);
  //node0_printf("meas_link total time = %f\n",ttime);
	return;
}
