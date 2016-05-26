// -----------------------------------------------------------------
// Main procedure for fermionic Wilson flow - still in progress
// the control file that worked for pure gauge flow is renamed to ocontrol.c
#define CONTROL
#include "wflow_includes.h"
// -----------------------------------------------------------------


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
  complex plp;
  plp = fploop(F_OFFSET(link),3);
  node0_printf("plp(linkstart) = %.5g  %.5g\n",plp.real,plp.imag);
	
	for (dir = 0; dir < 4; dir++){
		FORALLSITES(i,s){
			su3mat_copy(&(s->link[dir]), &(s->link0[dir]));
		}
	}
  // Flow to time tmax and save nsave configurations along the way
  //wflow(F_OFFSET(link), 0.0, tmax, 1);

	// Test variable epsilon gauge flow
	int imp_steps, num_eps = 0;
	//imp_steps = wflow_imp(F_OFFSET(link0), 0, tmax, 0);

	double *eps = wflow_imp_epsvals(F_OFFSET(link0), 0, tmax, 1), t_epsmax=0, eps_max=0.1;
	for (yep = 0; yep < (int)(tmax/epsilon)-1 && eps[yep] != 0; yep++){
		node0_printf("%g\n",eps[yep]);
		if (eps[yep] != 0) { num_eps += 1;}
		if (eps[yep] < eps_max-1e-7 && eps[yep-1] <= eps_max - 1e-7 ){ t_epsmax += eps[yep]; }
	}
	node0_printf("num_eps = %d t_epsmax = %g\n", num_eps, t_epsmax);
	double dt = eps_max*floor((tmax - t_epsmax)/(eps_max*(nsave - 1)));
  // Check ploops
  //plp = fploop(F_OFFSET(link),3);
  //node0_printf("plp(linkmax) = %.5g  %.5g\n",plp.real,plp.imag);
  
	// Adjoint Flow

	// Fixed epsilon flow
	//fermion_flow();

  // generate npbp random fields at tmax
  t = tmax;//tmax;//change to tmax for adj flow
  int ksi, n, j, l=nsave-2;
	double flowtime = 0, sum_eps, ttime=dclock(), tvar;
	Real ti, cut = 1e-7, counter = 0;
  complex dot1, dot2;
  node0_printf("ksi %g ",t);
  for (n=0; n<npbp; n++){
    grsourcen(EVENANDODD);
    ksi = F_OFFSET(ksi1) + n*sizeof(su3_vector);
    copy_latvec(F_OFFSET(g_rand), ksi, EVENANDODD);
    dot1 = dot_su3_latvec(ksi, ksi, EVENANDODD);
    node0_printf("%.10g  ", dot1.real);
  }
  //yep = fmeas_link(F_OFFSET(chi), F_OFFSET(link), F_OFFSET(psi), mass);
	//node0_printf("fmeas_link time = %f\n",dclock()-ttime);
	//ttime = dclock() - ttime;

  node0_printf("\n");

  node0_printf("\nBEGINNING FERMION ADJOINT FLOW\n\n");
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
			wflow_imp(eps_max,F_OFFSET(link1)+4*l*sizeof(su3_matrix), t_epsmax+l*dt, tmax-sum_eps);
		}
		else {
			//node0_printf("IMPROV2: ti = 0 tf = %g\n",tmax-sum_eps);
			wflow_imp(epsilon,F_OFFSET(link0), 0, tmax-sum_eps);
		}
		//node0_printf("FIXED: ti = %g tf = %g\n",tmax-sum_eps,t);
		wflow(F_OFFSET(link),tmax-sum_eps,t,0);
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
		/*tvar = dclock();
		yep = fmeas_link(F_OFFSET(chi),F_OFFSET(W0), F_OFFSET(psi), mass);
		ttime += dclock() - tvar;*/
    node0_printf("\n\n");
  }

  node0_printf("\nENDING FERMION ADJOINT FLOW \n\n");
 
  node0_printf("total flow time = %f\n",flowtime);
	node0_printf("meas_link total time = %f\n",ttime);
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
