// Function which calculates mod squared of a lattice vector

#include "wflow_includes.h"

complex meas_vec2(field_offset vec)
{
  complex tot, cc;
  tot.real = 0, tot.imag = 0;
  FORALLSITES(i,s){
    cc = su3_dot((su3_vector *)F_PT(s, vec), (su3_vector *)F_PT(s, vec))
    CSUM(tot, cc);
  }
	g_doublesum(&tot);
  return tot;
}
