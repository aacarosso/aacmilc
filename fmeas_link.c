// -----------------------------------------------------------------
// Measure fermion observables for general "even plus odd" quark actions:
//   pbp separately on even and odd sites
//   fermion action
//   links in all directions separately on even and odd sites

// Modified to accept gauge field offset as argument, and use the
// backward-stepped fermions ksi1,...,npbp

#include "ks_dyn_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Return total number of iterations
int fmeas_link(field_offset chi_off, field_offset W, field_offset psi_off, Real mass) {
  register int i;
  register site *s;
  int jpbp_reps, tot_iters = 0, miters, dir, ksi;
  double r_pbp_even, i_pbp_even, r_pbp_odd, i_pbp_odd, faction, pbp_pbp;
  complex cc;
  double_complex pbp_e, pbp_o;  // PBP
  double_complex ave_e, ave_o;

  int W_off[4];
  for (dir=XUP ; dir <= TUP; dir++){
    W_off[dir] = W + dir*sizeof(su3_matrix);
  }
  FILE *data = fopen("flow_data.txt","a");

  ave_e = dcmplx(0, 0);
  ave_o = dcmplx(0, 0);

	// loop over ksi0,1,2,3,4
  for (jpbp_reps = 0; jpbp_reps < npbp; jpbp_reps++) {
    faction = 0;
		ksi = F_OFFSET(ksi1) + jpbp_reps*sizeof(su3_vector);
    pbp_e = dcmplx(0, 0);
    pbp_o = dcmplx(0, 0);
    // Take back-flowed source ksi and do inversion
    // chi = Mdag lambda0; psi = M^{-1} lambda0
    fdslash(ksi, W, chi_off, EVENANDODD);
    scalar_mult_latvec(chi_off, -1.0, chi_off, EVENANDODD);
    scalar_mult_add_latvec(chi_off, ksi, 2.0*mass, chi_off, EVENANDODD);
    clear_latvec(psi_off, EVENANDODD);
//    miters = mat_invert_uml(F_OFFSET(g_rand), psi_off, chi_off, mass);
    // Adapted from mat_invert.c "preconditioning": chi <- M^dag * g_rand
    fdslash(ksi, W, F_OFFSET(ttt1[0]), EVENANDODD);
    scalar_mult_add_latvec(F_OFFSET(ttt1[0]), ksi,
                           -2.0 * mass, chi_off, EVENANDODD);
    scalar_mult_latvec(chi_off, -1.0, chi_off, EVENANDODD);
    miters = fks_congrad(chi_off, W, psi_off, mass, EVENANDODD);
    tot_iters += miters;
   
    // Fermion action = chi.psi
    // pbp on even sites = g_rand.psi
    FOREVENSITES(i, s) {
      cc = su3_dot((su3_vector *)F_PT(s, chi_off),
                   (su3_vector *)F_PT(s, psi_off));
      faction += cc.real;
      cc = su3_dot((su3_vector *)F_PT(s,ksi), (su3_vector *)F_PT(s, psi_off));
      CSUM(pbp_e, cc);
    }

    // pbp on odd sites
    FORODDSITES(i, s) {
      cc = su3_dot((su3_vector *)F_PT(s,ksi), (su3_vector *)F_PT(s, psi_off));
      CSUM(pbp_o, cc);
    }
/*
    // Now calculate the link differences
    tag0 = start_gather_site(psi_off, sizeof(su3_vector), TUP,
                             EVENANDODD, gen_pt[0]);
    tag1 = start_gather_site(psi_off, sizeof(su3_vector), XUP,
                             EVENANDODD, gen_pt[1]);
    tag2 = start_gather_site(psi_off, sizeof(su3_vector), YUP,
                             EVENANDODD, gen_pt[2]);
    tag3 = start_gather_site(psi_off, sizeof(su3_vector), ZUP,
                             EVENANDODD, gen_pt[3]);

    wait_gather(tag0);
    wait_gather(tag1);
    wait_gather(tag2);
    wait_gather(tag3);

    FORALLSITES(i, s) {
      mult_su3_mat_vec((su3_matrix *)F_PT(s, W_off[TUP]),
                       (su3_vector *)gen_pt[0][i], &(s->tempvec[0]));
      mult_su3_mat_vec((su3_matrix *)F_PT(s, W_off[XUP]),
                       (su3_vector *)gen_pt[1][i], &(s->tempvec[1]));
      mult_su3_mat_vec((su3_matrix *)F_PT(s, W_off[YUP]),
                       (su3_vector *)gen_pt[2][i], &(s->tempvec[2]));
      mult_su3_mat_vec((su3_matrix *)F_PT(s, W_off[ZUP]),
                       (su3_vector *)gen_pt[3][i], &(s->tempvec[3]));
    }

    FORALLSITES(i, s) {
      cc = su3_dot((su3_vector *)F_PT(s, ksi), &(s->tempvec[0]));
      if ((s->x) % 2 == 0)
        CSUM(pbp_e[1], cc);
      if ((s->x) % 2 == 1)
        CSUM(pbp_o[1], cc);
      if ((s->t) % 2 == 0)
        CSUM(pbp_e[2], cc);
      if ((s->t) % 2 == 1)
        CSUM(pbp_o[2], cc);

      cc = su3_dot((su3_vector *)F_PT(s, ksi), &(s->tempvec[1]));
      if ((s->x) % 2 == 0)
        CSUM(pbp_e[3], cc);
      if ((s->x) % 2 == 1)
        CSUM(pbp_o[3], cc);

      cc = su3_dot((su3_vector *)F_PT(s, ksi), &(s->tempvec[2]));
      if ((s->y) % 2 == 0)
        CSUM(pbp_e[4], cc);
      if ((s->y) % 2 == 1)
        CSUM(pbp_o[4], cc);

      cc = su3_dot((su3_vector *)F_PT(s, ksi), &(s->tempvec[3]));
      if ((s->z) % 2 == 0)
        CSUM(pbp_e[5], cc);
      if ((s->z) % 2 == 1)
        CSUM(pbp_o[5], cc);
    }
*/
    // Accumulate, normalize and print
    g_dcomplexsum(&pbp_e);
    g_dcomplexsum(&pbp_o);
    ave_e.real += pbp_e.real;
    ave_e.imag += pbp_e.imag;
    ave_o.real += pbp_o.real;
    ave_o.imag += pbp_o.imag;
    node0_printf("PBP: ");
    r_pbp_even = pbp_e.real * (2 / (double)volume);
    i_pbp_even = pbp_e.imag * (2 / (double)volume);
    r_pbp_odd = pbp_o.real * (2 / (double)volume);
    i_pbp_odd = pbp_o.imag * (2 / (double)volume);
    node0_printf("mass %.4g, %.8g %.8g %.8g %.8g ( %d of %d ) %d\n",
                   mass, r_pbp_even, r_pbp_odd, i_pbp_even, i_pbp_odd,
                   jpbp_reps + 1, npbp, miters);
  //  fprintf(data,"%g %g\n",r_pbp_even,r_pbp_odd);
  
    g_doublesum(&faction);
    node0_printf(" FACTION: mass = %.4g, %.8g ( %d of %d )\n",
                 mass, faction / (double)volume, jpbp_reps + 1, npbp);

    if (npbp > 1) {
      pbp_pbp = 0;
      FORALLSITES(i, s)
        su3vec_copy((su3_vector *)F_PT(s, psi_off), &(s->M_inv));

      clear_latvec(psi_off, EVENANDODD);
//      mat_invert_uml(F_OFFSET(M_inv), psi_off, chi_off, mass);
      // Adapted from mat_invert.c "preconditioning": temp <- M^dag * M_inv
      fdslash(F_OFFSET(M_inv), W, F_OFFSET(ttt1[0]), EVENANDODD);
      scalar_mult_add_latvec(F_OFFSET(ttt1[0]), F_OFFSET(M_inv),
                             -2.0 * mass, chi_off, EVENANDODD);
      scalar_mult_latvec(chi_off, -1.0, chi_off, EVENANDODD);
      miters = fks_congrad(chi_off, W, psi_off, mass, EVENANDODD);
      FORALLSITES(i, s) {
        cc = su3_dot((su3_vector *)F_PT(s, ksi), (su3_vector *)F_PT(s, psi_off));
        pbp_pbp += cc.real;
      }
      g_doublesum(&pbp_pbp);
      pbp_pbp /= (double)volume;
      node0_printf("TR_MM_INV: mass %.4g, %.8g ( %d of %d )\n",
                   mass, pbp_pbp, jpbp_reps + 1, npbp);
    }
  }

  // Print out averages, summed over the entire (half-)volume
  node0_printf("pbp: ");
  r_pbp_even = ave_e.real * (2 / ((double)volume*npbp));
  i_pbp_even = ave_e.imag * (2 / ((double)volume*npbp));
  r_pbp_odd = ave_o.real * (2 / ((double)volume*npbp));
  i_pbp_odd = ave_o.imag * (2 / ((double)volume*npbp));
  node0_printf("mass %.4g, %.8g %.8g %.8g %.8g ( ave over %d )\n",
                 mass, r_pbp_even, r_pbp_odd, i_pbp_even, i_pbp_odd, npbp);
  fprintf(data,"%g %g %g %g", r_pbp_even, r_pbp_odd, i_pbp_even, i_pbp_odd);
  fclose(data);
  return tot_iters;
}
// -----------------------------------------------------------------
