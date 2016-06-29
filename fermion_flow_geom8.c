// Improved fermion flow using wflow and wflow_imp_saves and fermion_flow_chunk
// Geometric series checkpoint method

#include "wflow_includes.h"

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

void fermion_flow_geom8()
{
	register int i;
	register site *s;
	int yep, dir, ksi, n, j, num_eps=0, num_tepsmax = 0;
  double t_epsmax=0, eps_max=0.15, *eps, flowtime=0, ttime=0, tvar, tvar1;
  Real t, cut = 1e-7;
  complex dot1;

	node0_printf("\n\nFERMION_FLOW_GEOM8\n\n");

	// flow gauge fields to tmax and determine eps values
	//rephase(ON);
	//block_and_fatten();
	//rephase(OFF);
  eps = wflow_imp_epsvals_geom8(F_OFFSET(link), 1);
  for (yep = 0; yep < (int)(tmax/epsilon)-1 && eps[yep] != 0; yep++){
    node0_printf("%g\n",eps[yep]);
    if (eps[yep] > cut) { num_eps += 1;}
    if (eps[yep] < eps_max-cut && eps[yep-1] <= eps_max - cut ){ 
			t_epsmax += eps[yep];
			num_tepsmax += 1;
		}
  }
  node0_printf("num_eps = %d t_epsmax = %g nume_tepsmax = %d\n", num_eps, t_epsmax, num_tepsmax);

	// generate npbp random fields at tmax
  node0_printf("ksi %g ",t);
  for (n=0; n<npbp; n++){
    grsourcen(EVENANDODD);
    ksi = F_OFFSET(ksi1) + n*sizeof(su3_vector);
    copy_latvec(F_OFFSET(g_rand), ksi, EVENANDODD);
    dot1 = dot_su3_latvec(ksi, ksi, EVENANDODD);
    node0_printf("%.10g  ", dot1.real);
  }
	
	// pbp with rand fermions and linkmax
	//rephase(ON);
	//block_and_fatten();
  //yep = fmeas_link(F_OFFSET(chi), F_OFFSET(link), F_OFFSET(psi), mass);
	//rephase(OFF);

	// pbp with rand fermions and link0
/*
  for (dir = 0; dir < 4; dir++){
    FORALLSITES(i,s){
      su3mat_copy(&(s->link0[dir]), &(s->link[dir]));
    }
  }
	rephase(ON);
	block_and_fatten();
  //yep = meas_link(F_OFFSET(chi), F_OFFSET(psi), mass);
  yep = fmeas_link(F_OFFSET(chi), F_OFFSET(link), F_OFFSET(psi), mass);
	rephase(OFF);
	double p1, p2;
	d_plaquette(&p1, &p2);
	node0_printf("p1 = %g p2 = %g\n", p1, p2);
*/

  node0_printf("\n");
	
  node0_printf("\nBEGIN FERMION ADJOINT FLOW\n\n");

	int R = 8;
	double ts = R*eps_max*floor((tmax - t_epsmax - eps[num_eps-1])/(R*eps_max));
	node0_printf("ts = %g\n",ts);

	// t > t_epsmax + 7/8*ts REGIME
	fermion_flow_chunk(F_OFFSET(link1)+3*4*sizeof(su3_matrix), t_epsmax+7.0/R*ts, tmax, eps_max);

	// t > t_epsmax REGIME
	fermion_flow_chunk(F_OFFSET(link1)+2*4*sizeof(su3_matrix), 
		t_epsmax+6.0/R*ts, t_epsmax+7.0/R*ts, eps_max);

	double times1[] = {5.0/8*ts+t_epsmax};
	wflow_imp_saves(eps_max,F_OFFSET(link1)+1*4*sizeof(su3_matrix), t_epsmax+4.0/R*ts,
		t_epsmax+5.0/R*ts+epsilon, times1, 1);
	fermion_flow_chunk(F_OFFSET(link1)+3*4*sizeof(su3_matrix), t_epsmax+5.0/R*ts,
		t_epsmax+6.0/R*ts, eps_max);
	fermion_flow_chunk(F_OFFSET(link1)+1*4*sizeof(su3_matrix), t_epsmax+4.0/R*ts,
		t_epsmax+5.0/R*ts, eps_max);

	double times2[] = {1.0/R*ts+t_epsmax,2.0/R*ts+t_epsmax,3.0/R*ts+t_epsmax};
	wflow_imp_saves(eps_max,F_OFFSET(link1), t_epsmax,
		t_epsmax+3.0/R*ts+epsilon, times2, 3);
	fermion_flow_chunk(F_OFFSET(link1)+3*4*sizeof(su3_matrix), t_epsmax+3.0/R*ts,
		t_epsmax+4.0/R*ts, eps_max);
	fermion_flow_chunk(F_OFFSET(link1)+2*4*sizeof(su3_matrix), t_epsmax+2.0/R*ts,
		t_epsmax+3.0/R*ts, eps_max);
	fermion_flow_chunk(F_OFFSET(link1)+1*4*sizeof(su3_matrix), t_epsmax+1.0/R*ts,
		t_epsmax+2.0/R*ts, eps_max);
	fermion_flow_chunk(F_OFFSET(link1), t_epsmax,
		t_epsmax+1.0/R*ts, eps_max);
	
	// t < t_epsmax REGIME: fixed checkpoints
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
  wflow_imp_saves(epsilon, F_OFFSET(link0), 0, t_epsmax, saves, nsave);
  for (i = nsave-1; i >= 0; i--){
    node0_printf("sum_eps(%d*sd) = %g, eps[%d*sd] = %g\n",i,sum_eps(eps,i*sd),i,eps[i*sd]);
    fermion_flow_chunk(F_OFFSET(link0)+i*4*sizeof(su3_matrix), sum_eps(eps, i*sd), 
      sum_eps(eps, j), eps[i*sd]); 
    j = i*sd;
  }

	//fermion_flow_chunk(F_OFFSET(link0), 0, t_epsmax, epsilon);

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

  //node0_printf("total flow time = %f\n",flowtime);
  node0_printf("meas_link total time = %f\n",ttime);
	return;
}
