// Performs a forward fermion flow (unstable)

void fermion_forward_flow()
{	
	register int i;
	register site *s;
	int istep, dir;
	Real t;
	complex dot1;
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
	return;
}

