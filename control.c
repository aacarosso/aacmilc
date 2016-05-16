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
	int imp_steps;
	imp_steps = wflow_imp(F_OFFSET(link0), 0, tmax, 0);

  // Check ploops
  //plp = fploop(F_OFFSET(link),3);
  //node0_printf("plp(linkmax) = %.5g  %.5g\n",plp.real,plp.imag);
  
/*
	plp = fploop(F_OFFSET(link0),3);
  node0_printf("plp(link0) = %.5g  %.5g\n",plp.real,plp.imag);
  plp = fploop(F_OFFSET(link1), 3);
  node0_printf("plp(link1) = %.5g  %.5g\n",plp.real,plp.imag);
  plp = fploop(F_OFFSET(link2), 3);
  node0_printf("plp(link2) = %.5g  %.5g\n",plp.real,plp.imag);
  plp = fploop(F_OFFSET(link3), 3);
  node0_printf("plp(link3) = %.5g  %.5g\n",plp.real,plp.imag);

	wflow(F_OFFSET(link0),0.0,0.12,0);
	plp = fploop(F_OFFSET(link),3);
  node0_printf("plp(link1) = %.5g  %.5g\n",plp.real,plp.imag);

	wflow(F_OFFSET(link0),0.0,0.24,0);
	plp = fploop(F_OFFSET(link),3);
  node0_printf("plp(link2) = %.5g  %.5g\n",plp.real,plp.imag);

	wflow(F_OFFSET(link0),0.00,0.36,0);
	plp = fploop(F_OFFSET(link),3);
  node0_printf("plp(link3) = %.5g  %.5g\n",plp.real,plp.imag);

	wflow(F_OFFSET(link0),0.00,0.5,0);
	plp = fploop(F_OFFSET(link),3);
  node0_printf("plp(linkmax) = %.5g  %.5g\n",plp.real,plp.imag);

  node0_printf("\n");
*/
  // generate npbp random fields at tmax
  t = 0;//change to tmax for adj flow
  int ksi, n;
  complex dot1, dot2;
  node0_printf("ksi: ");
  for (n=0; n<npbp; n++){
    grsourcen(EVENANDODD);
    ksi = F_OFFSET(ksi1) + n*sizeof(su3_vector);
    copy_latvec(F_OFFSET(g_rand), ksi, EVENANDODD);
    dot1 = dot_su3_latvec(ksi, ksi, EVENANDODD);
    node0_printf("%g  ", dot1.real);
  }
	double ttime = dclock(), tvar;
  //yep = fmeas_link(F_OFFSET(chi), F_OFFSET(link), F_OFFSET(psi), mass);
	//node0_printf("fmeas_link time = %f\n",dclock()-ttime);
	ttime = dclock() - ttime;

  node0_printf("\n");

  // Adjoint flow loop 

  node0_printf("\nBEGINNING FERMION ADJOINT FLOW \n\n");
	int k = nsave-1;
	Real ti, N = floor(tmax/(epsilon*nsave)), cut = 1e-7;
	Real counter = 0;
	double flowtime = 0;
	node0_printf("floor = %g\n",N);
	imp_steps = 0;
	
  for (istep == 0; t > 0; istep++){
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
		//counter += t/epsilon;
		//wflow(F_OFFSET(link0),0,t,0);
		//imp_steps += wflow_imp(F_OFFSET(link0), 0, t, 0);
    node0_printf("\nksi: ");
    // loop over npbp fields
    for (n=0; n<npbp; n++){
      ksi = F_OFFSET(ksi1) + n*sizeof(su3_vector);
      copy_latvec(ksi, F_OFFSET(lambda3), EVENANDODD);
      fermion_adjointstep(F_OFFSET(lambda3));
      copy_latvec(F_OFFSET(lambda0), ksi, EVENANDODD);
      dot1 = dot_su3_latvec(ksi, ksi, EVENANDODD);
      node0_printf("%g  ", dot1.real);
    }
    node0_printf("\n");
		/*tvar = dclock();
		yep = fmeas_link(F_OFFSET(chi),F_OFFSET(W0), F_OFFSET(psi), mass);
		ttime += dclock() - tvar;*/

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
    t -= epsilon;
  }

  node0_printf("\nENDING FERMION ADJOINT FLOW \n\n");
 
  node0_printf("total flow time = %f\n",flowtime);
	node0_printf("meas_link total time = %f\n",ttime);
	node0_printf("counter = %g\n",counter);
	node0_printf("imp_steps = %d\n",imp_steps);
/*
  // Do a complete forward flow of lambda0
  node0_printf("\nBEGIN FORWARD FERMION FLOW\n\n");
  t = epsilon;
  for (istep = 0; t<= tmax; istep++){
    for (dir = 0; dir < 4; dir++){
      FORALLSITES(i,s){
        su3mat_copy(&(s->link0[dir]),&(s->link[dir]));
      }
    }
    wflow(t);
    clear_latvec(F_OFFSET(lambda1), EVENANDODD);
    clear_latvec(F_OFFSET(lambda2), EVENANDODD);
    clear_latvec(F_OFFSET(lambda3), EVENANDODD);
    dot1 = dot_su3_latvec(F_OFFSET(lambda0),F_OFFSET(lambda0),EVENANDODD);
    node0_printf("lambda0^2 = %g\n",dot1.real);
    fermion_forwardstep(F_OFFSET(lambda0));
    dot1 = dot_su3_latvec(F_OFFSET(lambda3),F_OFFSET(lambda3),EVENANDODD);
    node0_printf("lambda3^2 = %g\n",dot1.real);
    copy_latvec(F_OFFSET(lambda3),F_OFFSET(lambda0),EVENANDODD);
    t += epsilon;
  }
  node0_printf("END FORWARD FERMION FLOW\n\n");
*/
  
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
