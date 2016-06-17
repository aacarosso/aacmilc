// -----------------------------------------------------------------
// Main procedure for fermionic Wilson flow
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

	//double *eps;
	//eps = wflow_imp_epsvals_geom8(F_OFFSET(link0), 1);

	//wflow_imp(0.1, F_OFFSET(link2), 5.42, tmax);

  // Flow to time tmax and save nsave configurations along the way
  //wflow(F_OFFSET(link), 0.0, tmax, 1);

	//wflow_imp(1*epsilon, F_OFFSET(link0), 0, tmax);

	//fermion_flow_imp();

	fermion_flow_geom8();

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
