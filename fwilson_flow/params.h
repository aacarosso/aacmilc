// -----------------------------------------------------------------
// Structure for passing simulation parameters to each node
#ifndef _PARAMS_H
#define _PARAMS_H
#include "../include/macros.h"  // For MAXFILENAME

typedef struct {
  int stopflag;             // 1 if it is time to stop
  int nx, ny, nz, nt;       // Lattice dimensions
  int iseed;                // For random numbers
  Real beta, mass;          // Gauge coupling, fermion mass
  int startflag;            // What to do for beginning lattice
  int saveflag;             // What to do with lattice at end // added by aac   
  char startfile[MAXFILENAME], savefile[MAXFILENAME];
  char stringLFN[MAXFILENAME];  // ILDG LFN if applicable


  // Smearing parameters
  Real alpha_hyp0, alpha_hyp1, alpha_hyp2;

  // Spectrum parameters -- source timeslice and increment
  int src_start, src_inc, n_src;

  // Inversion parameters
  int npbp;                 // Number of stochastic sources
  int niter;                // Maximum number of CG iterations
  int nrestart;             // Maximum number of CG restarts
  Real rsqmin;              // For deciding on convergence

  // wilson flow paramaters
  int num_block, nsave;
  Real epsilon, tmax, tblock[100]; // deleted alpha_hyp0,1,2
} params;
#endif
// -----------------------------------------------------------------
