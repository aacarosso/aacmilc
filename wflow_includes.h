// -----------------------------------------------------------------
// Include files for staggered eigenvalues
#include <stdio.h>
#include <stdlib.h>
#include <string.h>             // For strlen
#include <math.h>
#include "../include/config.h"  // Keep this first
#include "../include/complex.h"
#include "../include/su3.h"
#include "../include/macros.h"
#include "lattice.h"
#include "defines.h"
#include "../include/comdefs.h"
#include "../include/io_lat.h"
#include "../include/generic.h"
#include "../include/generic_ks.h"
#include "../include/generic_nhyp.h"
#include "../include/dirs.h"
#include "../include/field_alloc.h"
// -----------------------------------------------------------------



// Prototypes for functions in high level code
int setup();
int readin(int prompt);
void grsource_imp(field_offset dest, Real M, int parity);
void grsourcen(int parity);
void dslash(field_offset chi, field_offset psi, int parity);

// Wilson Flow stuff

void fstout_step_rk(su3_matrix *S[4], anti_hermitmat *A[4], double eps, int last);
void staple(su3_matrix *stp[4]);
void wflow(field_offset off, Real ti, Real tf, int savelink); 
void wflow_imp(double eps, field_offset off, Real ti, Real tf);
double *wflow_imp_epsvals(field_offset off, Real ti, Real tf, int savelink);
void fermion_flow();
void fermion_adjointstep(field_offset flow_vec, double eps);
void fermion_forwardstep(field_offset flow_vec, double eps);
void mcrg_block(Real t, int blmax);//aac

// nHYP stuff specific to for MCRG-blocked measurements
void clear_disp(int *disp);
void diag_su3(su3_matrix* Q, complex *f);
void block_nhyp_mcrg(int num, int block);
void block_mcrg(int num, int block);

// Wilson loop stuff
void make_loop_table2();
void blocked_gauge_loops(int block, double *result);
void path(int *dir, int *sign, int length, su3_matrix *resmat);
void blocked_path(int block, int *dir, int *sign,
                  int length, su3_matrix *resmat);

// Polyakov loop stuff
complex blocked_ploop(int block, int dir);
complex fploop(field_offset W, int dir);

// CG stuff
int CG_wrapper(field_offset src, su3_vector **props, Real M, int parity);
int fCG_wrapper(field_offset src, field_offset W, su3_vector **props, Real M, int parity);
int ks_congrad(field_offset src, field_offset dest, Real M, int parity);
int fks_congrad(field_offset src, field_offset W, field_offset dest, Real M, int parity);

// Dslash stuff and associated helper functions
// (latter were in v6/generic_ks/d_congrad5.c for some reason)
void dslash_wrapper(su3_vector *src, su3_vector *dest, int parity);
void fdslash_wrapper(su3_vector *src, field_offset W, su3_vector *dest, int parity);
void dslash(field_offset chi, field_offset psi, int parity);
void fdslash(field_offset chi, field_offset W, field_offset psi, int parity);
void fdslash_special(field_offset chi, field_offset psi, int parity,
                    msg_tag **tag, int start);
void dslash_special(field_offset chi, field_offset psi, int parity,
                    msg_tag **tag, int start);
void clear_latvec(field_offset v, int parity);
void copy_latvec(field_offset src, field_offset dest, int parity);
void scalar_mult_latvec(field_offset src, Real scalar, 
                            field_offset dest, int parity);
void scalar_mult_add_latvec(field_offset src1, field_offset src2,
                            Real scalar, field_offset dest, int parity);
complex dot_su3_latvec(field_offset off1, field_offset off2, int parity);
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Order parameter routines
void d_plaquette_a(double *ss_plaq, double *st_plaq);
void meas_plaq();
int meas_link(field_offset chi, field_offset psi, Real mass);
int fmeas_link(field_offset chi, field_offset W, field_offset psi, Real mass);
// -----------------------------------------------------------------
