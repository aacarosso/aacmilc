// Performs one adjoint (backward) epsilon step for fermions
// flow_vec0,1,2,3 must be declared in that order in lattice.h
// typically an adjoint flow field is called lambda

#include "wflow_includes.h"

//-------------------------------------------------------------------


void fermion_adjointstep(field_offset flow_vec){
  Real eps = epsilon;
  register int i;
  register site *s;
  int j, off[4];

	for (j = 0; j < 4; j++) off[3-j] = flow_vec - j*sizeof(su3_vector);

  // Step EVEN sites
  
  // Obtain lambda2
  fdslash(off[3], F_OFFSET(W2), off[2], ODD);
  fdslash(off[2], F_OFFSET(W2), off[2], EVEN);
  scalar_mult_latvec(off[2], 3*eps/4, off[2], EVEN);

  // Obtain lambda1
  fdslash(off[2], F_OFFSET(W1), off[1], ODD);
  fdslash(off[1], F_OFFSET(W1), off[1], EVEN);
  scalar_mult_add_latvec(off[3], off[1], 8*eps/9,
      off[1],EVEN);

  // Obtain lambda0
  scalar_mult_add_latvec(off[1], off[2], -8/9,
      off[0], EVEN);
  fdslash(off[0], F_OFFSET(W0), off[0], ODD);
  fdslash(off[0], F_OFFSET(W0), off[0], EVEN);
  scalar_mult_add_latvec(off[2], off[0], eps/4,
      off[0], EVEN);
  scalar_mult_add_latvec(off[1], off[0], 1,
      off[0], EVEN);

  // copy even sites of lambda0
  FOREVENSITES(i,s){
    su3vec_copy((su3_vector *)F_PT(s,off[0]), &(s->tempvec0));
  }
  clear_latvec(off[0], EVENANDODD);
  clear_latvec(off[1], EVENANDODD);
  clear_latvec(off[2], EVENANDODD);

  // Step ODD sites

  // Obtain lambda2
  fdslash(off[3], F_OFFSET(W2), off[2], EVEN);
  fdslash(off[2], F_OFFSET(W2), off[2], ODD);
  scalar_mult_latvec(off[2], 3*eps/4, off[2], ODD);

  // Obtain lambda1
  fdslash(off[2], F_OFFSET(W1), off[1], EVEN);
  fdslash(off[1], F_OFFSET(W1), off[1], ODD);
  scalar_mult_add_latvec(off[3], off[1], 8*eps/9,
      off[1], ODD);

  // Obtain lambda0
  scalar_mult_add_latvec(off[1], off[2], -8/9,
      off[0], ODD);
  fdslash(off[0], F_OFFSET(W0), off[0], EVEN);
  fdslash(off[0], F_OFFSET(W0), off[0], ODD);
  scalar_mult_add_latvec(off[2], off[0], eps/4,
      off[0], ODD);
  scalar_mult_add_latvec(off[1], off[0], 1,
      off[0], ODD);
  
  // copy even sites into lambda0
  FOREVENSITES(i,s){
    su3vec_copy(&(s->tempvec0), (su3_vector *)F_PT(s,off[0]));
  }

}
