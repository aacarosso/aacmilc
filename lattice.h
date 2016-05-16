// -----------------------------------------------------------------
// Define global scalars and fields in the lattice
#ifndef _LATTICE_H
#define _LATTICE_H

#include "defines.h"
#include "../include/generic_quark_types.h"
#include "../include/macros.h"    // For MAXFILENAME
#include "../include/io_lat.h"    // For gauge_file
#include "../include/su3.h"
#include "../include/random.h"    // For double_prn
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// The lattice is an array of this site struct
typedef struct {
  short x, y, z, t;   // Coordinates of this site
  char parity;        // Is it even or odd?
  int index;          // Index in the array

#ifdef SITERAND
  // The state information for a random number generator
  double_prn site_prn;
  // Align to double word boundary (kludge for Intel compiler)
  int space1;
#endif

  // Gauge field
  //su3_matrix link[4]; //wilson_flow uses link[12] not [4]
  //copied-----start
  su3_matrix link[12];
  su3_matrix fieldstrength[6];
  su3_matrix hyplink1[12], hyplink2[12];//tempmat1,2 and staple already below 
  //copied---end

  // Program-dependent fields
  // Staggered phases, which have been absorbed into the matrices
  // Also the antiperiodic boundary conditions
  Real phase[4];

  // Staggered complex vectors
  su3_vector g_rand;            // Gaussian random vector
  su3_vector chi;               // Gaussian random source vector
  su3_vector psi;               // Solution vector = Kinverse * chi
  su3_vector p;                 // CG change vector
  su3_vector mp;                // CG temporary vector
  su3_vector r;                 // CG residual vector
  su3_vector ttt1[1];           // For ../generic_ks/mat_invert.c
  su3_vector M_inv;

  // Temporary vectors and matrices for gauge field and fermion flow
  su3_matrix link0[4];         // aaa
  su3_matrix link1[4];
  su3_matrix link2[4];
  su3_matrix link3[4];
  su3_vector lambda0;           // temp fields used in adjointstep
  su3_vector lambda1;
  su3_vector lambda2;
  su3_vector lambda3;
  su3_vector ksi1;              // npbp=5 fields to flow
  su3_vector ksi2;
  su3_vector ksi3;
  su3_vector ksi4;
  su3_vector ksi5;
  su3_vector chi0;              // temp fields used in forward step
  su3_vector chi1;
  su3_vector chi2;
  su3_vector chi3;
  su3_vector tempvec0;          // for fermion steps
  su3_matrix W0[4];
  su3_matrix W1[4];
  su3_matrix W2[4];
  su3_matrix W3[4];


  // Spectrum stuff
  su3_vector quark_source;      // For ../generic_ks/fpi_2.c
  su3_vector propmat[3];        // For three source colors
  su3_vector propmat2[3];       // For three source colors

  // Temporary vectors and matrices -- for gauge fixing, spectrum
  su3_vector tempvec[4];  // One for each direction
  su3_matrix tempmat1, tempmat2, staple;
} site;

// copied---start
// Defines for index on field_strength
#define FS_XY 0
#define FS_XZ 1
#define FS_YZ 2
#define FS_XT 3
#define FS_YT 4
#define FS_ZT 5
//copied----end
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Definition of global variables
#ifdef CONTROL
#define EXTERN
#else
#define EXTERN extern
#endif

// Global variables
EXTERN int nx, ny, nz, nt;  // Lattice dimensions
EXTERN int volume;          // Volume of lattice
EXTERN int iseed;           // Random number seed
EXTERN int npbp;            // Number of stochastic sources
EXTERN int niter, nrestart;
EXTERN Real beta, mass, rsqmin;

EXTERN double g_ssplaq, g_stplaq;
EXTERN double_complex linktrsum;
EXTERN u_int32type nersc_checksum;
EXTERN char stringLFN[MAXFILENAME];  // ILDG LFN if applicable (AAC added)
EXTERN char startfile[MAXFILENAME], savefile[MAXFILENAME]; // added savefile (aac)
EXTERN int startflag;   // Beginning lattice: CONTINUE, RELOAD, FRESH
EXTERN int saveflag;    // 1 if we will save the lattice; // added by aac
EXTERN int total_iters;
EXTERN int phases_in;   // 1 if KS and BC phases absorbed into matrices

// Some of these global variables are node dependent
// They are set in "make_lattice()"
EXTERN int sites_on_node;       // Number of sites on this node
EXTERN int even_sites_on_node;  // Number of even sites on this node
EXTERN int odd_sites_on_node;   // Number of odd sites on this node
EXTERN int number_of_nodes;     // Number of nodes in use
EXTERN int this_node;           // Node number of this node

// Each node maintains a structure with the pseudorandom number
// generator state
EXTERN double_prn node_prn;

EXTERN gauge_file *startlat_p;

//copied--------------start
// Loop stuff
#define nloop 6
#define nreps 1
#define max_num 300

// Global definitions for general action (checked?)
EXTERN int loop_ind[nloop][10], loop_length[nloop];
EXTERN int loop_table[nloop][max_num][10];
EXTERN int loop_num[nloop], loop_char[max_num];
EXTERN Real loop_coeff[nloop][nreps];
EXTERN int ch, loop_ch[nloop][max_num];
EXTERN Real loop_term[max_num][nreps];
EXTERN int hyp1ind[4][4][4];
EXTERN int hyp2ind[4][4];
//copied----------end

// The lattice is a single global variable
// (actually this is the part of the lattice on this node)
EXTERN site *lattice;

// Vectors for addressing
// Generic pointers, for gather routines
#define N_POINTERS 16   // Needed by ../generic_ks/nl_spectrum.c
EXTERN char **gen_pt[N_POINTERS];

EXTERN su3_matrix *gauge_field[4];
EXTERN su3_matrix *gauge_field_thin[4];

// Spectrum stuff -- source timeslice and increment
EXTERN int src_start, src_inc, n_src;

// nHYP stuff
EXTERN double alpha_smear[3];
EXTERN su3_matrix *hyplink1[4][4];
EXTERN su3_matrix *hyplink2[4][4];

EXTERN su3_matrix *Staple1[4][4];
EXTERN su3_matrix *Staple2[4][4];
EXTERN su3_matrix *Staple3[4];

EXTERN su3_matrix *tempmat1;
EXTERN su3_matrix *tempmat2;    // Used in Polyakov loop calculation

//copied----------------start
// Wilson flow stuff
EXTERN int nsave;
EXTERN double tmax;
EXTERN double epsilon;
EXTERN su3_matrix *tempmat1;

// MCRG blocking stuff
//EXTERN Real alpha_smear[3];
EXTERN int num_block;
EXTERN Real tblock[100];
//copied-----------------end

#endif // _LATTICE_H
// -----------------------------------------------------------------
