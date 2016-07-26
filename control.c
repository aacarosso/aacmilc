// -----------------------------------------------------------------
// Main procedure for fermionic Wilson flow
// the control file that worked for pure gauge flow is renamed to ocontrol.c
#define CONTROL
#include "wflow_includes.h"
// -----------------------------------------------------------------

int nearest_checkpt(int i, int *E){
	int j, ret = -1;
	for (j=i-1; j>=0; j--){
		if (E[j] == 1){
			ret = j;
			break;
		}
	}
	return ret;
}

void saves(int *allsaves, int *E, int R){
	int i, j=0;
	for ( i=0; i<R; i++){
		if ( E[i] == 1 ){
			allsaves[j] = i;
			node0_printf("  checkpoint[%d] = %d\n", j, allsaves[j]);
			j += 1;
		}
		if (j > nsave){
			node0_printf("NUM CHECKPOINTS > NSAVE!\n");
			break;
		}
	}
	return;
}

int link_num(int *allsaves, int checkpt){
	int i, number;
	for ( i=0; i < nsave; i++){
		if (allsaves[i] == checkpt){
			number = i;
		}
	}
	return number;
}

// -----------------------------------------------------------------
int main(int argc, char **argv) {
  int iters, prompt, dir, istep, yep;
  register int i;
  register site *s;
  double dtime;
  Real t;
  // Setup
  setlinebuf(stdout); // DEBUG
  // Remap standard I/O
  if (remap_stdio_from_args(argc, argv) == 1)
    terminate(1);

  initialize_machine(&argc, &argv);
  g_sync();
  prompt = setup();
  // Load input and run
  // readin does a rephase(ON) on the read-in lattice
  if (readin(prompt) != 0) {
    node0_printf("ERROR in readin, aborting\n");
    terminate(1);
  }
  dtime = -dclock();

#ifndef PBOUND
  node0_printf("\nANTIPERIODIC BC in time direction\n");
#else
  node0_printf("\nPERIODIC BC in time direction\n");
#endif
	//-------------------------------------------------------------------
		
	//rephase(ON);
	//block_and_fatten();
	//rephase(OFF);
	//node0_printf("BLOCKED link0 !\n");
	// copy t = 0 link into link0
	for (dir = 0; dir < 4; dir++){
		FORALLSITES(i,s){
			su3mat_copy(&(s->link[dir]), &(s->link0[dir]));
		}
	}

double *eps, t_epsmax = 0, cut = 1e-7, eps_max = 0.1;
int num_tepsmax = 0, num_eps = 0;
  eps = wflow_imp_epsvals(F_OFFSET(link0), 0, 10, 1, eps_max);
  for (yep = 0; yep < (int)(tmax/epsilon)-1 && eps[yep] != 0; yep++){
    node0_printf("%g\n",eps[yep]);
    if (eps[yep] > cut) { num_eps += 1;}
    if (eps[yep] < eps_max-cut && eps[yep-1] <= eps_max - cut ){
      t_epsmax += eps[yep];
      num_tepsmax += 1;
    }
  }
  node0_printf("num_eps = %d t_epsmax = %g num_tepsmax = %d\n", num_eps, t_epsmax, num_tepsmax);

int k, j, numsaves, num;
int R = (int)(pow(2, nsave-1)), r, p;
int nstep = R*(int)floor((double)(num_eps/R));
int *E = (int *)malloc((R+1)*sizeof(int));
E[0] = 1;
node0_printf("R = %d nstep = %d E[0] = %d\n", R, nstep, E[0]);
r = 1;
p = 0;
for (i=1; i<=R; i++){
  if (i == p + R/pow(2,r)){
    E[i] = 1;
    p += R/pow(2,r);
    r += 1;
  }
	else{ E[i] = 0;}
  node0_printf("E[%d] = %d\n",i,E[i]);
}

for ( i = 0; i <= num_eps; i++){
	node0_printf("sum_eps(%d) = %f\n",i,sum_eps(eps,i));
}

int ns = nsave, l;
p = R-1;
k = 1;
int *allsaves = (int *)malloc(nsave*sizeof(int));
saves(allsaves, E, R);
for (p=R-1; p>=0; p--){
	node0_printf("p = %d E[%d] = %d ns = %d\n",p,p,E[p],ns,k);
	if (E[p]==0){
		node0_printf(" loop\n");
		i = p;
		numsaves = 0;
		k = 1;
		while (E[i] == 0 && ns < nsave){
			//saves[ns] = i;
			ch = nearest_checkpt(i,E);
			E[i] = 1;
			numsaves += 1;
			node0_printf("  (while) ns = %d ch = %d E[%d] = %d numsaves = %d\n",
				ns, ch, i, E[i], numsaves);
			ns += 1;
			i -= 1;
		}
		saves(allsaves, E, R);
		double *newsaves = (double *)malloc(numsaves*sizeof(double));
		l = 0;
		for ( j=0; j<nsave; j++){
			if (allsaves[j] > ch && allsaves[j] <= p){
				newsaves[l] = sum_eps(eps, allsaves[j]*nstep/R);
				node0_printf(" newsaves[%d] = sum_eps(eps, %d) = %g\n", l, allsaves[j]*nstep/R, 
					sum_eps(eps, allsaves[j]*nstep/R));
				l += 1;
			}
		}
		node0_printf(" wflow ti offset: link_num = %d\n", link_num(allsaves,ch));

/*
		for (j=1; j<=numsaves; j++){
			saves[j-1] = sum_eps(eps,nstep*(p-numsaves+j)/R);
			node0_printf("  p-numsaves+j = %d saves[%d] = %g\n",p-numsaves+j, j-1,saves[j-1]);
		}
		node0_printf(" wflow: offset nsave-1-numsaves = %d, ti = %g tf = %g\n",
			nsave-1-numsaves, sum_eps(eps,ch*nstep/R), sum_eps(eps,p*nstep/R)+epsilon);
*/
		//wflow_imp_saves(eps[ch*R], F_OFFSET(link0)+4*(nsave-1-numsaves)*sizeof(su3_matrix), 
		//	sum_eps(eps,ch*R), sum_eps(eps,p*R)+epsilon, saves, numsaves, eps_max);
		free(newsaves);
	}
	//fermion_flow_chunk(F_OFFSET(link0)+4*(nsave-k)*sizeof(su3_matrix), 
	//node0_printf(" fwflow: offset (nsave-k) = %d, ti = %g tf = %g\n",nsave-k,
	//	sum_eps(eps,p*nstep/R), sum_eps(eps,(p+1)*nstep/R));
	node0_printf("fwflow: offset = %d, ti = %g, tf = %g, E[%d] = 0\n", link_num(allsaves, p),
		sum_eps(eps,p*nstep/R), sum_eps(eps,(p+1)*nstep/R), p);
	//	sum_eps[p*R], sum_eps[(p+1)*R], eps[p*R], eps_max);
	E[p] = 0;
	ns -= 1;
}




	//double *eps;
	//eps = wflow_imp_epsvals_geom8(F_OFFSET(link0), 1);

  // Flow to time tmax and save nsave configurations along the way
  //wflow(F_OFFSET(link), 0.0, tmax, 0);

	//wflow_imp(30*epsilon, F_OFFSET(link0), 0, tmax, 0.30);

	//fermion_flow_imp(0.10);

	//fermion_flow_eqlstep(0.10);

	//fermion_flow_adaptive(0.10);

	//fermion_flow();

	//fermion_flow_geom8(0.1);

	//-------------------------------------------------------------------
  rephase(ON);  

  leanlinks();
  fflush(stdout);

  node0_printf("RUNNING COMPLETED\n");
  dtime += dclock();
  node0_printf("Time = %.4g seconds\n", dtime);
  node0_printf("total_iters = %d\n\n", iters);
  fflush(stdout);

  // Save lattice if requested -- no phases in links
  if (saveflag != FORGET)
    node0_printf("not forget\n");
    save_lattice(saveflag, savefile, stringLFN);
  normal_exit(0);
  return 0;
}
// -----------------------------------------------------------------
