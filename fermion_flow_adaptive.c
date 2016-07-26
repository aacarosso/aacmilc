#include "wflow_includes.h"

// Main procedure for adaptive stepsize fermion flow.
void fermion_flow_adaptive(double eps_max)
{
  register int i;
  register site *s;
  int yep, dir, ksi, n, j, num_eps=0, num_tepsmax = 0, istep;
  double t_epsmax=0, *eps, flowtime=0, ttime=0, tvar, tvar1, ep, eta;
  Real t = tmax, cut = 1e-7;
  complex dot1, dotp;

  node0_printf("\n\nFERMION_FLOW_ADAPTIVE\n\n");

  // flow gauge fields to tmax and determine eps values
  eps = wflow_imp_epsvals(F_OFFSET(link), 0, tmax, 1, eps_max);
  for (yep = 0; yep < (int)(tmax/epsilon)-1 && eps[yep] != 0; yep++){
    node0_printf("%g\n",eps[yep]);
    if (eps[yep] > cut) { num_eps += 1;}
    if (eps[yep] < eps_max-cut && eps[yep-1] <= eps_max - cut ){
      t_epsmax += eps[yep];
      num_tepsmax += 1;
    }
  }
  node0_printf("num_eps = %d t_epsmax = %g num_tepsmax = %d\n", num_eps, t_epsmax, num_tepsmax);

  // generate npbp random fields at tmax
  node0_printf("ksi %g ",t);
  for (n=0; n<npbp; n++){
    grsourcen(EVENANDODD);
    ksi = F_OFFSET(ksi1) + n*sizeof(su3_vector);
    copy_latvec(F_OFFSET(g_rand), ksi, EVENANDODD);
    dot1 = dot_su3_latvec(ksi, ksi, EVENANDODD);
    node0_printf("%.10g  ", dot1.real);
  }
  node0_printf("\n");

  node0_printf("\nBEGIN FERMION ADJOINT FLOW\n\n");
	ep = epsilon;
	for (istep = 0; t > cut; istep++){
		
		wflow_imp(eps_max, F_OFFSET(link0)+4*4*sizeof(su3_matrix), 7.62, t-ep, eps_max);
		wflow_imp(ep, F_OFFSET(link), t-ep, t, eps_max);
		
		dotp = dot1;
		node0_printf("ep = %g\n",ep);
    // adjoint step each field
  	node0_printf("ksi %g ",t-ep);
    for (n=0; n<npbp; n++){
      ksi = F_OFFSET(ksi1) + n*sizeof(su3_vector);
      copy_latvec(ksi, F_OFFSET(lambda3), EVENANDODD);
      fermion_adjointstep(F_OFFSET(lambda3), ep);
      copy_latvec(F_OFFSET(lambda0), ksi, EVENANDODD);
      dot1 = dot_su3_latvec(ksi, ksi, EVENANDODD);
      node0_printf("%.10g  ", dot1.real);
    }
		node0_printf("\n\n");
		
		t -= ep;
		if (istep == 0){
			eta = ep*(sqrt(dotp.real) - sqrt(dot1.real));
			node0_printf("eta = %g\n", eta);
		}
		node0_printf("eta/() = %g\n", eta/(sqrt(dotp.real)-sqrt(dot1.real)));
		if (eta/(sqrt(dotp.real)-sqrt(dot1.real)) > ep + epsilon){
			ep += epsilon;
		}
		if (t-ep < -cut){
			ep = t;
		}
	}
	return;
}

