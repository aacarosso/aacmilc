// -----------------------------------------------------------------
// SU(3) S4b order parameters setup
#include "wflow_includes.h"
#include <string.h>

int initial_set();
void make_fields();
#define IF_OK if (status == 0)

// Each node has a params structure for passing simulation parameters
#include "params.h"
params par_buf;
// -----------------------------------------------------------------



// -----------------------------------------------------------------
int setup() {
  int prompt;

  // Print banner, get volume, seed
  prompt = initial_set();
  // Initialize the node random number generator
  initialize_prn(&node_prn, iseed, volume + mynode());
  // Initialize the layout functions, which decide where sites live
  setup_layout();
  // Allocate space for lattice, set up coordinate fields
  make_lattice();
  // Set up neighbor pointers and comlink structures
  make_nn_gathers();
  // Allocate space for fields
  make_fields();
  // Set up staggered phase vectors, boundary conditions
  phaseset();
  // Allocate space for temporary su3_matrix field
  FIELD_ALLOC(tempmat1, su3_matrix);
  return prompt;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// On node zero, read lattice size, seed and send to others
int initial_set() {
  int prompt, status;
  if (mynode() == 0) {
    // Print banner
    printf("Fermionic SU(3) Wilson Flow with PBP measurement\n");
    printf("MIMD version 7ish\n");
    printf("Machine = %s, with %d nodes\n", machine_type(), numnodes());
    printf("nHYP links, reading alpha_smear parameters from infile\n");
    printf("  IR_STAB = %.4g\n", (Real)IR_STAB);
    printf("  EPS_SQ = %.4g\n", (Real)EPS_SQ);
#ifdef NHYP_DEBUG
    printf("NHYP_DEBUG turned on\n");
#endif
#ifdef CG_DEBUG
    printf("CG_DEBUG turned on\n");
#endif
#ifdef CGTIME
    printf("CGTIME turned on\n");
#endif
    time_stamp("start");
    status = get_prompt(stdin, &prompt);
    IF_OK status += get_i(stdin, prompt, "nx", &par_buf.nx);
    IF_OK status += get_i(stdin, prompt, "ny", &par_buf.ny);
    IF_OK status += get_i(stdin, prompt, "nz", &par_buf.nz);
    IF_OK status += get_i(stdin, prompt, "nt", &par_buf.nt);
    IF_OK status += get_i(stdin, prompt, "iseed", &par_buf.iseed);

    if (status > 0)
      par_buf.stopflag = 1;
    else
      par_buf.stopflag = 0;
  }

  // Broadcast parameter buffer from node 0 to all other nodes
  broadcast_bytes((char *)&par_buf, sizeof(par_buf));
  if (par_buf.stopflag != 0)
    normal_exit(0);

  nx = par_buf.nx;
  ny = par_buf.ny;
  nz = par_buf.nz;
  nt = par_buf.nt;
  iseed = par_buf.iseed;

  this_node = mynode();
  number_of_nodes = numnodes();
  volume = nx * ny * nz * nt;
  total_iters = 0;
  return prompt;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Read in parameters for spectrum calculation
int readin(int prompt) {
  // prompt=1 indicates prompts are to be given for input
  int status, i;
  Real x;

  // On node zero, read parameters and send to all other nodes
  if (this_node == 0) {
    printf("\n\n");
    status = 0;
    // Wilson flow parameters
		IF_OK status += get_i(stdin, prompt, "nsave", &par_buf.nsave);
    IF_OK status += get_f(stdin, prompt, "epsilon", &par_buf.epsilon);
    if (par_buf.epsilon == 0) {
      node0_printf("ERROR: epsilon=%g won't get you anywhere\n",
                   par_buf.epsilon);
      status++;
    }
    IF_OK status += get_f(stdin, prompt, "tmax", &par_buf.tmax);
    if (par_buf.epsilon * par_buf.tmax < 0)
      node0_printf("WARNING: epsilon and tmax have different signs\n");

    // Limit to only one mass for now
    IF_OK status += get_f(stdin, prompt, "mass", &par_buf.mass);

/*    IF_OK status += get_i(stdin, prompt, "number_of_PF", &par_buf.num_masses); //aac

    // Hasenbusch mass /aac
    IF_OK status += get(stdin, prompt, "Hasenbusch mass", &par_buf.MH); //aac*/

    // Smearing parameters
    IF_OK status += get_f(stdin, prompt, "alpha_hyp0", &par_buf.alpha_hyp0);
    IF_OK status += get_f(stdin, prompt, "alpha_hyp1", &par_buf.alpha_hyp1);
    IF_OK status += get_f(stdin, prompt, "alpha_hyp2", &par_buf.alpha_hyp2);

    // A maximum of 100 tvalues to perfom blocking should be eough
    IF_OK status += get_i(stdin, prompt, "num_block", &par_buf.num_block);
    if (par_buf.num_block > 100) {
      node0_printf("ERROR: Need to recompile for num_block > 100\n");
      status++;
    }
    for (i = 0; i < par_buf.num_block; i++) {
      IF_OK status += get_f(stdin, prompt, "tblock", &par_buf.tblock[i]);
      // Make sure we're going in the right direction
      if (i > 0 && fabs(par_buf.tblock[i]) <= fabs(par_buf.tblock[i - 1])) {
        node0_printf("ERROR: We require tblock be sorted\n");
        node0_printf("ERROR: tblock[%d]=%g; tblock[%d]=%g\n",
                     i, par_buf.tblock[i], i - 1, par_buf.tblock[i - 1]);
        status++;
        break;
      }
    }

    // Spectrum source time slice and increment
    IF_OK status += get_i(stdin, prompt, "source_start", &par_buf.src_start);
    IF_OK status += get_i(stdin, prompt, "source_inc", &par_buf.src_inc);
    IF_OK status += get_i(stdin, prompt, "n_sources", &par_buf.n_src);

    // Number of stochastic sources
    IF_OK status += get_i(stdin, prompt, "npbp", &par_buf.npbp);

    // Maximum conjugate gradient iterations and restarts
    IF_OK status += get_i(stdin, prompt, "max_cg_iterations", &par_buf.niter);
    IF_OK status += get_i(stdin, prompt, "max_cg_restarts", &par_buf.nrestart);

    // Error per site for conjugate gradient
    IF_OK {
      status += get_f(stdin, prompt, "error_per_site", &x);
      par_buf.rsqmin = x * x;
    }

    // Find out what kind of starting lattice to use
    IF_OK status += ask_starting_lattice(stdin, prompt, &par_buf.startflag,
                                         par_buf.startfile);
    // Find out what to do with lattice at end
    IF_OK status += ask_ending_lattice(stdin, prompt, &(par_buf.saveflag),
                                       par_buf.savefile);

    if (status > 0)
      par_buf.stopflag = 1;
    else
      par_buf.stopflag = 0;
  }

  // Broadcast parameter buffer from node0 to all other nodes
  broadcast_bytes((char *)&par_buf, sizeof(par_buf));
  if (par_buf.stopflag != 0)
    normal_exit(0);
	nsave = par_buf.nsave;
  epsilon = par_buf.epsilon;
  tmax    = par_buf.tmax;
  num_block = par_buf.num_block;
  for (i=0; i<num_block;i++)
    tblock[i]=par_buf.tblock[i];

  mass = par_buf.mass;
  alpha_smear[0] = par_buf.alpha_hyp0;
  alpha_smear[1] = par_buf.alpha_hyp1;
  alpha_smear[2] = par_buf.alpha_hyp2;

  src_start = par_buf.src_start;
  src_inc = par_buf.src_inc;
  n_src = par_buf.n_src;
  npbp = par_buf.npbp;

  niter = par_buf.niter;
  nrestart = par_buf.nrestart;
  rsqmin = par_buf.rsqmin;
/*  num_masses = par_buf.num_masses; //aac
  if (num_masses > 1) //aac
    MH = par_buf.MH; //aac*/

  startflag = par_buf.startflag;
  saveflag = par_buf.saveflag; 
  strcpy(startfile, par_buf.startfile);
  strcpy(savefile, par_buf.savefile); 
  strcpy(stringLFN, par_buf.savefile);


  // Do whatever is needed to get lattice
  if (startflag == CONTINUE)
    rephase(OFF);

  startlat_p = reload_lattice(startflag, startfile);
  // If a lattice was read in, put in staggered phase factors
  // and antiperiodic boundary condition
  //phases_in = OFF; two lines commented out; don't want phases
  //rephase(ON);
  return 0;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Allocate all space for fields
// Amount Malloced is pure guesswork/imagination
void make_fields() {
  FIELD_ALLOC_VEC(gauge_field, su3_matrix, 4);
  FIELD_ALLOC_VEC(gauge_field_thin, su3_matrix, 4);

  FIELD_ALLOC_MAT_OFFDIAG(hyplink1, su3_matrix, 4);
  FIELD_ALLOC_MAT_OFFDIAG(hyplink2, su3_matrix, 4);
  FIELD_ALLOC_MAT_OFFDIAG(Staple1, su3_matrix, 4);
  FIELD_ALLOC_MAT_OFFDIAG(Staple2, su3_matrix, 4);
  FIELD_ALLOC_VEC(Staple3, su3_matrix, 4);
  FIELD_ALLOC(tempmat1, su3_matrix);
  FIELD_ALLOC(tempmat2, su3_matrix);

  node0_printf("Mallocing %.1f MBytes per node for fields\n",
               (double)sites_on_node * 61 * sizeof(su3_matrix) / 1e6);
}
// -----------------------------------------------------------------
