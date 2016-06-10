// Improved fermion flow using wflow and wflow_imp

// Saves nsave configurations, equal spacing in eps_max regime

#include "wflow_includes.h"

void fermion_flow_imp() {
	int istep, yep, dir, ksi, n, j, l=nsave-2, num_eps=0;
  double t_epsmax=0, eps_max=0.1, *eps, dt, flowtime=0, sum_eps, ttime=0, tvar, tvar1;
  Real t, cut = 1e-7;
  complex dot1;

	node0_printf("\n\nFERMION_FLOW_IMP\n\n");

	// flow gauge fields to tmax and determine eps values
  eps = wflow_imp_epsvals(F_OFFSET(link0), 0, tmax, 1);
  for (yep = 0; yep < (int)(tmax/epsilon)-1 && eps[yep] != 0; yep++){
    node0_printf("%g\n",eps[yep]);
    if (eps[yep] > cut) { num_eps += 1;}
    if (eps[yep] < eps_max-1e-7 && eps[yep-1] <= eps_max - 1e-7 ){ t_epsmax += eps[yep]; }
  }
  node0_printf("num_eps = %d t_epsmax = %g\n", num_eps, t_epsmax);
  dt = eps_max*floor((tmax - t_epsmax)/(eps_max*(nsave - 1)));

	// generate npbp random fields at tmax
  t = tmax;
  node0_printf("ksi %g ",t);
  for (n=0; n<npbp; n++){
    grsourcen(EVENANDODD);
    ksi = F_OFFSET(ksi1) + n*sizeof(su3_vector);
    copy_latvec(F_OFFSET(g_rand), ksi, EVENANDODD);
    dot1 = dot_su3_latvec(ksi, ksi, EVENANDODD);
    node0_printf("%.10g  ", dot1.real);
  }
  //yep = fmeas_link(F_OFFSET(chi), F_OFFSET(link), F_OFFSET(psi), mass);

  node0_printf("\n");

  node0_printf("\nBEGIN FERMION ADJOINT FLOW\n\n");
  j = num_eps-1;
  sum_eps = eps[j];
  for (istep = 0; t > cut; istep++){
    // set wflow_imp ti
    if (pow(t - (tmax - sum_eps),2) < cut){
      j -= 1;
      sum_eps += eps[j];
    }
    node0_printf(" j = %d sum_eps = %g l = %d\n", j, sum_eps, l);

    // adjust checkpoint for gauge field
    if (pow(t - (t_epsmax +l*dt),2) < cut ) l -= 1;

    // flow gauge fields to t, obtain W's
    tvar = dclock();
    if (t - t_epsmax > cut){
      //node0_printf("IMPROV1: ti = t_epsmax+l*dt = %g tf = %g\n", t_epsmax+l*dt, tmax-sum_eps);
      wflow_imp(eps[j],F_OFFSET(link1)+4*l*sizeof(su3_matrix), t_epsmax+l*dt, tmax-sum_eps);
      //wflow_imp(eps[j],F_OFFSET(link1)+4*l*sizeof(su3_matrix), t_epsmax+l*dt, t-epsilon);
    }
    else {
      //node0_printf("IMPROV2: ti = 0 tf = %g\n",tmax-sum_eps);
      wflow_imp(epsilon,F_OFFSET(link0), 0, tmax-sum_eps);
      //wflow_imp(epsilon,F_OFFSET(link0), 0, t-epsilon);
    }
    //node0_printf("FIXED: ti = %g tf = %g\n",tmax-sum_eps,t);
    wflow(F_OFFSET(link), tmax-sum_eps, t, 0);
    //wflow(F_OFFSET(link),t-epsilon,t,0);
    flowtime += dclock() - tvar;
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
    //tvar1 = dclock();
    //yep = fmeas_link(F_OFFSET(chi),F_OFFSET(W0), F_OFFSET(psi), mass);
    //ttime += dclock() - tvar1;
    node0_printf("\n\n");
  }
	tvar1 = dclock();
  yep = fmeas_link(F_OFFSET(chi),F_OFFSET(link0), F_OFFSET(psi), mass);
	ttime += dclock() - tvar1;
  node0_printf("\nEND FERMION ADJOINT FLOW\n\n");

  node0_printf("total flow time = %f\n",flowtime);
  node0_printf("meas_link total time = %f\n",ttime);
	return;
}
