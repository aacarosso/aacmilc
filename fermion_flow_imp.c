// Improved fermion flow using wflow and wflow_imp

// Fixed checkpoint saving in t > t_epsmax regime

#include "wflow_includes.h"
/*
// sums the elements of array eps up to and including the index max_index.
double sum_eps(double *eps, int max_index)
{
  int i;
  double sum = 0;
  for ( i=0; i < max_index; i++){
    sum += eps[i];
    //node0_printf("eps[%d] = %g sum = %g\n",i,eps[i],sum);
  }
  return sum;
}
*/
// Main procedure for improved fermion flow with fixed checkpoints
void fermion_flow_imp(double eps_max) {
	register int i;
	register site *s;
	int istep, yep, dir, ksi, n, j, l=nsave-2, num_eps=0, num_tepsmax = 0;
  double t_epsmax=0, *eps, dt, flowtime=0, sumeps, ttime=0, tvar, tvar1;
  Real t, cut = 1e-7;
  complex dot1;

	node0_printf("\n\nFERMION_FLOW_IMP\n\n");

	// flow gauge fields to tmax and determine eps values
  eps = wflow_imp_epsvals(F_OFFSET(link0), 0, tmax, 1, eps_max);
  for (yep = 0; yep < (int)(tmax/epsilon)-1 && eps[yep] != 0; yep++){
    node0_printf("%g\n",eps[yep]);
    if (eps[yep] > cut) { num_eps += 1;}
    if (eps[yep] < eps_max-1e-7 && eps[yep-1] <= eps_max - cut ){
			t_epsmax += eps[yep];
			num_tepsmax += 1;
		}
  }
  node0_printf("num_eps = %d t_epsmax = %g num_tepsmax = %d\n", num_eps, t_epsmax, num_tepsmax);
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

  node0_printf("\n");

  node0_printf("\nBEGIN FERMION ADJOINT FLOW\n\n");
  j = num_eps-1;
  sumeps = eps[j];
  for (istep = 0; t > /*t_epsmax +*/ cut; istep++){
    // set wflow_imp ti
    if (pow(t - (tmax - sumeps),2) < cut){
      j -= 1;
      sumeps += eps[j];
    }
    node0_printf(" j = %d sum_eps = %g l = %d\n", j, sumeps, l);

    // adjust checkpoint for gauge field
    if (pow(t - (t_epsmax +l*dt),2) < cut ) l -= 1;

    // flow gauge fields to t, obtain W's
    tvar = dclock();

    wflow_imp(epsilon,F_OFFSET(link0), 0, t-epsilon, eps_max);
/*    if (t - t_epsmax > cut){
      wflow_imp(eps[j],F_OFFSET(link0)+4*(l+1)*sizeof(su3_matrix), t_epsmax+l*dt, t-epsilon, eps_max);
    }
    else {
      wflow_imp(epsilon,F_OFFSET(link0), 0, t-epsilon, eps_max);
    }*/
    wflow(F_OFFSET(link),t-epsilon,t,0);
    flowtime += dclock() - tvar;
    node0_printf("ksi %g ", t-epsilon);

    // adjoint step each field
    for (n=0; n<npbp; n++){
      ksi = F_OFFSET(ksi1) + n*sizeof(su3_vector);
      copy_latvec(ksi, F_OFFSET(lambda3), EVENANDODD);
      fermion_adjointstep(F_OFFSET(lambda3), epsilon);
      copy_latvec(F_OFFSET(lambda0), ksi, EVENANDODD);
      dot1 = dot_su3_latvec(ksi, ksi, EVENANDODD);
      node0_printf("%.10g  ", dot1.real);
    }

    t -= epsilon;
    node0_printf("\n\n");
	}
/*
	// t < t_epsmax regime
  j = num_tepsmax;
  int sd = (int)floor(num_tepsmax/nsave);
  node0_printf("sd = %d\n", sd);
  for ( i = 0; i <= num_eps; i++){
    node0_printf("sum_eps(%d) = %f\n",i,sum_eps(eps,i));
  }
  double *saves = (double *)malloc(nsave*sizeof(double));
  for ( i = 0; i < nsave; i++){
    saves[i] = sum_eps(eps, i*sd);
    node0_printf("saves[%d] = %g\n", i, saves[i]);
  }
  wflow_imp_saves(epsilon, F_OFFSET(link0), 0, t_epsmax, saves, nsave, eps_max);
  for (i = nsave-1; i >= 0; i--){
    node0_printf("sum_eps(%d*sd) = %g, eps[%d*sd] = %g\n",i,sum_eps(eps,i*sd),i,eps[i*sd]);
    fermion_flow_chunk(F_OFFSET(link0)+i*4*sizeof(su3_matrix), sum_eps(eps, i*sd),
      sum_eps(eps, j), eps[i*sd], eps_max);
    j = i*sd;
  }
  free(saves);
*/
  // pbp with t = 0 fermions and link0
  for (dir = 0; dir < 4; dir++){
    FORALLSITES(i,s){
      su3mat_copy(&(s->link0[dir]), &(s->link[dir]));
    }
  }
  rephase(ON);
  block_and_fatten();
  ttime = dclock();
  yep = fmeas_link(F_OFFSET(chi), F_OFFSET(link), F_OFFSET(psi), mass);
  ttime = dclock() - ttime;
  rephase(OFF);

	node0_printf("\nEND FERMION ADJOINT FLOW\n\n");

  node0_printf("total flow time = %f\n",flowtime);
  node0_printf("meas_link total time = %f\n",ttime);
	return;
}
