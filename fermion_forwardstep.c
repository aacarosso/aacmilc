// Performs a forward epsilon step for fermions
// uses field offset of flow_vec to find flow_vec0,1,2,3 offsets
// so the flow_vec0,1,2,3 must be declared in sequence in lattice.h

#include "wflow_includes.h"

//-------------------------------------------------------------------

void fermion_forwardstep(field_offset flow_vec){
  Real eps = epsilon;
  register int i;
  register site *s;
  int j, off[4];
  for (j = 0; j < 4; j++) off[j] = flow_vec + j*sizeof(su3_vector);

  // Step ODD sites

  // Obtain lambda1
  fdslash(off[0], F_OFFSET(W0), off[1], EVEN);
  fdslash(off[1], F_OFFSET(W0), off[1], ODD);
  scalar_mult_add_latvec(off[0], off[1], eps/4,
      off[1], ODD);
  
  // Obtain lambda2
  fdslash(off[1], F_OFFSET(W1), off[2], EVEN);
  fdslash(off[2], F_OFFSET(W1), off[2], ODD);
  scalar_mult_latvec(off[2], 8*eps/9, off[2], ODD);
  scalar_mult_add_latvec(off[2], off[1], -8/9, 
      off[2], ODD);
  scalar_mult_add_latvec(off[2], off[0], 17/9, 
      off[2], ODD);

  // Obtain lambda3
  fdslash(off[2], F_OFFSET(W2), off[3], EVEN);
  fdslash(off[3], F_OFFSET(W2), off[3], ODD);
  scalar_mult_add_latvec(off[1], off[3], 3*eps/4, 
      off[3], ODD);

  // copy odd sites of lambda3
  FORODDSITES(i,s){
    su3vec_copy((su3_vector *)F_PT(s,off[3]), &(s->tempvec0));
  }
  clear_latvec(off[1], EVENANDODD);
  clear_latvec(off[2], EVENANDODD);
  clear_latvec(off[3], EVENANDODD);

  // Step EVEN sites

  // Obtain lambda1
  fdslash(off[0], F_OFFSET(W0), off[1], ODD);
  fdslash(off[1], F_OFFSET(W0), off[1], EVEN);
  scalar_mult_add_latvec(off[0], off[1], eps/4,
      off[1], EVEN);
  
  // Obtain lambda2
  fdslash(off[1], F_OFFSET(W1), off[2], ODD);
  fdslash(off[2], F_OFFSET(W1), off[2], EVEN);
  scalar_mult_latvec(off[2], 8*eps/9, off[2], EVEN);
  scalar_mult_add_latvec(off[2], off[1], -8/9, 
      off[2], EVEN);
  scalar_mult_add_latvec(off[2], off[0], 17/9, 
      off[2], EVEN);

  // Obtain lambda3
  fdslash(off[2], F_OFFSET(W2), off[3], ODD);
  fdslash(off[3], F_OFFSET(W2), off[3], EVEN);
  scalar_mult_add_latvec(off[1], off[3], 3*eps/4, 
      off[3], EVEN);
  
  // copy odd sites into lambda3
  FORODDSITES(i,s){
    su3vec_copy(&(s->tempvec0), (su3_vector *)F_PT(s,off[3]));
  }

}
  
  
  
 
  
  
  

  
