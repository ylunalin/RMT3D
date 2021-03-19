#include "mg_probs.hh"

double mg_prob::run(multigrid *mg, bool profile, int c) {
	const int hard_cap=20;
	mgo[c].set();
	// do initial iters, time each one
	double start=sw.elapsed(c),res,res0;

	if (profile) {
		res0 = resid(mg);
	}

    res = resid(mg);
	for (int i = 0; i < mgp[c].iters(); i++) {

		// time each v_cycle
		sw.tic(c);
        // argument mode = 1 turn on symmetric v cycle
		mg->v_cycle(0);
		sw.toc(c);
	}

    res = resid(mg);
	int extra = 0;
    // now count all the extras
	while (res > max_norm_tol && extra < hard_cap) {
		// mark increase
		extra++;
		// increase the number we do each time
        for(int n=0;n<extra;n++) {
            sw.tic(c);
            // argument mode = 1 turn on symmetric v cycle
            mg->v_cycle(0);
            sw.toc(c);
        }
		// get new residual
		res = resid(mg);
	}

    double l2res = mg->l2_error_all() / (gm->m*gm->n*gm->o);
	if(gm->rank==0 && res>max_norm_tol) printf("\nmg_prob:: rank 0, problem ID: %s maximum norm of residual=%g (tol %g), square l2 norm of residual=%g (tol %g) not reached tolerence after %d iterations\n", *id, res, max_norm_tol, l2res, tol, extra);
	if(gm->rank==0 && std::isnan(res)) printf("\nmg_prob:: rank 0, problem ID: %s maximum norm of residual is nan wtf\n", *id);
	// do two extras
	for (int i = 0; i < 2; i++) {
        // argument mode = 1 turn on symmetric v cycle
        mg->v_cycle(0);
    }

	double ltime = (sw.elapsed(c) - start)/gm->procs;
	double time;
	MPI_Allreduce(&ltime,&time,1,MPI_DOUBLE,MPI_SUM,gm->cart);
	int iters = mgp[c].iters() + extra*(extra+1)/2;
	double log_res = log10(res);
	double digs_below = log_tol - log_res;

	if (gm->rank == 0) {
		sprintf(out[c],
		"%s: (up,down)=(%d,%d), %d iters, %.3g sec/it, %.3g digs below",
		id[c],mg->v_up,mg->v_down,iters,time/iters,digs_below);
	}

	if (profile) {

		// how many iters did we take, digits did we gain,
		// and time did we take? (not counting "insurance" v_cycle)

		double digits = log10(res0) - log_res;
		if (gm->rank == 0) {
			sprintf(out[c],
				"%s, %.3g digs/it, %.3g digs/sec    ",
				out[c],digits/iters,digits/time);
		}

		if (mgo[c].searching && (iters == 1 || digs_below < 1.5*digits/iters)) mgo[c].record(digits/time);

		mgp[c].update(extra,digs_below,digits/iters);
	} else {
		mgp[c].update(extra);
	}

	if (gm->rank == 0) sprintf(out[c],"%s\n",out[c]);
    return res;
}


double * var_den_mg_prob::external_ptr(int i, int j, int k){
        int loci = i-gm->ai, locj = j-gm->aj, lock = k-gm->ak;
		return all_stencils[loci + sm*locj + smn * lock].st;
}

void var_den_mg_prob::set_rhoinv(const double *rho){
    for(int i=0;i<cellmn4*cello4;i++){
        rho_inv[i] = 1./rho[i];
    }
}

void var_den_fem::set_stencils() {
    for(int k=0;k<so;k++){
        for(int j=0;j<sn;j++){
            for(int i=0;i<sm;i++){
                int len=0;
                double * tmp = st.get(i+gm->ai,j+gm->aj,k+gm->ak,len);
                all_stencils[i + sm*j + smn*k].set(tmp, len);

                // index of region in the domain
                int si, sj, sk;
                int reg = st.reg(i+gm->ai,j+gm->aj,k+gm->ak, si, sj, sk);

                // The bounds on rhoinv we take into account
                // if si=0, then the stencil is on the low face on x, we can't access rho_inv[i-1]
                // if si=2, then it's on the hi face on x, we can't access rho_inv[i]
                int rxlo = (si==0)?i:i-1, rxhi = (si==2)?i-1:i;
                int rylo = (sj==0)?j:j-1, ryhi = (sj==2)?j-1:j;
                int rzlo = (sk==0)?k:k-1, rzhi = (sk==2)?k-1:k;

                // The size of the stencil in x, y, z direction
                int ni = (si==1)?3:2;
                int nj = (sj==1)?3:2;

                double *rho_array = new double[len];
                for(int ri = 0;ri<len;ri++) {
                    rho_array[ri] = 0.0;
                }

                for(int rk=rzlo; rk<=rzhi; rk++){
                    for(int rj=rylo; rj<=ryhi; rj++){
                        for(int ri=rxlo; ri<=rxhi; ri++){
                            // For each cell can contribute to all of its 8 nodes
                            for(int tmpk = rk-rzlo; tmpk <= rk-rzlo+1; tmpk++){
                                for(int tmpj = rj-rylo; tmpj <= rj-rylo+1; tmpj++){
                                    for(int tmpi = ri-rxlo; tmpi <= ri-rxlo+1; tmpi++){
                                        int stencil_entry_index = tmpi + tmpj*ni + tmpk * ni * nj;
                                        rho_array [ stencil_entry_index] += rho_inv[ c0 + ri + cellm4 * rj + cellmn4 * rk];
                                    }
                                }
                            }
                        }
                    }
                }

                for(int ri=0;ri<len;ri++) {
                    all_stencils[i + + sm*j + smn*k].st[ri] *= rho_array[ri] / num_rho_contrib[reg][ri];
                }

                delete [] rho_array;
            }
        }
    }
}

int const var_den_fem::num_rho_contrib[27][27] = {
{1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{1, 2, 1, 1, 2, 1, 1, 2, 1, 1, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{1, 1, 2, 2, 1, 1, 1, 1, 2, 2, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{1, 2, 1, 2, 4, 2, 1, 2, 1, 1, 2, 1, 2, 4, 2, 1, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{1, 1, 2, 2, 1, 1, 1, 1, 2, 2, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{1, 2, 1, 1, 2, 1, 1, 2, 1, 1, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{1, 1, 1, 1, 2, 2, 2, 2, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{1, 2, 1, 1, 2, 1, 2, 4, 2, 2, 4, 2, 1, 2, 1, 1, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{1, 1, 1, 1, 2, 2, 2, 2, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{1, 1, 2, 2, 1, 1, 2, 2, 4, 4, 2, 2, 1, 1, 2, 2, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{1, 2, 1, 2, 4, 2, 1, 2, 1, 2, 4, 2, 4, 8, 4, 2, 4, 2, 1, 2, 1, 2, 4, 2, 1, 2, 1},
{1, 1, 2, 2, 1, 1, 2, 2, 4, 4, 2, 2, 1, 1, 2, 2, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{1, 1, 1, 1, 2, 2, 2, 2, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{1, 2, 1, 1, 2, 1, 2, 4, 2, 2, 4, 2, 1, 2, 1, 1, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{1, 1, 1, 1, 2, 2, 2, 2, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{1, 2, 1, 1, 2, 1, 1, 2, 1, 1, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{1, 1, 2, 2, 1, 1, 1, 1, 2, 2, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{1, 2, 1, 2, 4, 2, 1, 2, 1, 1, 2, 1, 2, 4, 2, 1, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{1, 1, 2, 2, 1, 1, 1, 1, 2, 2, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{1, 2, 1, 1, 2, 1, 1, 2, 1, 1, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}
};

void var_den_mac::set_stencils() {
    for(int k=0;k<so;k++){
        for(int j=0;j<sn;j++){
            for(int i=0;i<sm;i++){
                int len=0;
                double * tmp = st.get(i+gm->ai,j+gm->aj,k+gm->ak,len);
                all_stencils[i + sm*j + smn*k].set(tmp, len);

                // index of region in the domain
                int si, sj, sk;
                int reg = st.reg(i+gm->ai,j+gm->aj,k+gm->ak, si, sj, sk);

                // The bounds on rhoinv we take into account
                // if si=0, then the stencil is on the low face on x, we can't access rho_inv[i-1]
                // if si=2, then it's on the hi face on x, we can't access rho_inv[i]
                int rxlo = (si==0)?i:i-1;
                int rylo = (sj==0)?j:j-1;
                int rzlo = (sk==0)?k:k-1;

                // The size of the stencil in x, y, z direction
                int ni = (si==1)?3:2;
                int nj = (sj==1)?3:2;

                double *rho_array = new double[len];
                for(int ri = 0;ri<len;ri++) {
                    rho_array[ri] = 0.0;
                }

                // Since the stencil is very sparse, we use seperate for loops for each cartesian direction
                // And add rho_inv to each entry in the stencil
                {
                    int rj = j, rk = k;
                    // can be written even simpler, but this is more readable
                    // if only two cells are used in a direction, the weights are just average of the two cells
                    if(si==0 || si==2) {
                        int stencil_entry_index = ni*(rj-rylo) + ni*nj*(rk-rzlo);
                        double avg = 0.5*(rho_inv[c0 + rxlo + rj*cellm4 + rk*cellmn4] + rho_inv[c0 + rxlo+1 + rj*cellm4 + rk*cellmn4]);
                        rho_array[stencil_entry_index] += avg;

                        stencil_entry_index += 1;
                        rho_array[stencil_entry_index] += avg;
                    } else if (si==1) {
                        // if 3 cels are used, the two cells at the end use the average of itself and the center cell.
                        // The center cell sums up the two averages
                        int stencil_entry_index =  ni*(rj-rylo) + ni*nj*(rk-rzlo);
                        double avg1 = 0.5*(rho_inv[c0 + rxlo + rj*cellm4 + rk*cellmn4 ] + rho_inv[c0 + rxlo+1 + rj*cellm4 + rk*cellmn4]);
                        double avg2 = 0.5*(rho_inv[c0 + rxlo +1 + rj*cellm4 + rk*cellmn4 ] + rho_inv[c0 + rxlo+2 + rj*cellm4 + rk*cellmn4]);

                        rho_array[stencil_entry_index] += avg1;

                        stencil_entry_index += 1;
                        rho_array[stencil_entry_index] +=  avg1 + avg2;

                        stencil_entry_index += 1;
                        rho_array[stencil_entry_index] += avg2;
                    }
                }

                // Y direction
                {
                    int ri = i, rk = k;
                    // if only two cells are used in a direction, the weights are just average of the two cells
                    if(sj==0 || sj==2) {
                        int stencil_entry_index = (ri-rxlo) + ni*nj*(rk-rzlo);
                        double avg = 0.5*(rho_inv[c0 + ri + rylo*cellm4 + rk*cellmn4] + rho_inv[c0 + ri + (rylo+1)*cellm4 + rk*cellmn4]);
                        rho_array[stencil_entry_index] += avg;

                        stencil_entry_index += ni;
                        rho_array[stencil_entry_index] += avg;
                    } else if (sj==1) {
                        // if 3 cels are used, the two cells at the end use the average of itself and the center cell.
                        // The center cell sums up the two averages
                        int stencil_entry_index = (ri-rxlo) + ni*nj*(rk-rzlo);
                        double avg1 = 0.5*(rho_inv[c0 + ri + rylo*cellm4 + rk*cellmn4] + rho_inv[c0 + ri + (rylo+1)*cellm4 + rk*cellmn4]);
                        double avg2 = 0.5*(rho_inv[c0 + ri + (rylo+1)*cellm4 + rk*cellmn4] + rho_inv[c0 + ri + (rylo+2)*cellm4 + rk*cellmn4]);

                        rho_array[stencil_entry_index] += avg1;

                        stencil_entry_index += ni;
                        rho_array[stencil_entry_index] +=  avg1 + avg2;

                        stencil_entry_index += ni;
                        rho_array[stencil_entry_index] += avg2;
                    }
                }

                // Z direction
                {
                    int ri = i, rj = j;
                    // if only two cells are used in a direction, the weights are just average of the two cells
                    if(sk==0 || sk==2) {
                        int stencil_entry_index = (ri-rxlo) + (rj-rylo)* ni;
                        double avg = 0.5*(rho_inv[c0 + ri + rj*cellm4 + rzlo*cellmn4] + rho_inv[c0 + ri + rj*cellm4 + (rzlo+1)*cellmn4]);
                        rho_array[stencil_entry_index] += avg;

                        stencil_entry_index += ni*nj;
                        rho_array[stencil_entry_index] += avg;
                    } else if (sk==1) {
                        // if 3 cels are used, the two cells at the end use the average of itself and the center cell.
                        // The center cell sums up the two averages
                        int stencil_entry_index = (ri-rxlo) + (rj-rylo)* ni;
                        double avg1 = 0.5*(rho_inv[c0 + ri + rj*cellm4 + rzlo*cellmn4] + rho_inv[c0 + ri + rj*cellm4 + (rzlo+1)*cellmn4]);
                        double avg2 = 0.5*(rho_inv[c0 + ri + rj*cellm4 + (rzlo+1)*cellmn4] + rho_inv[c0 + ri + rj*cellm4 + (rzlo+2)*cellmn4]);

                        rho_array[stencil_entry_index] += avg1;

                        stencil_entry_index += ni*nj;
                        rho_array[stencil_entry_index] +=  avg1 + avg2;

                        stencil_entry_index += ni*nj;
                        rho_array[stencil_entry_index] += avg2;
                    }
                }

                for(int ri=0;ri<len;ri++) {
                    if( num_rho_contrib[reg][ri]!=0) all_stencils[i + + sm*j + smn*k].st[ri] *= rho_array[ri] / num_rho_contrib[reg][ri];
                }
                delete [] rho_array;
            }
        }
    }
}


int const var_den_mac::num_rho_contrib[27][27] = {
{3, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{1, 4, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{1, 3, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{1, 0, 4, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{0, 1, 0, 1, 5, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{0, 1, 1, 4, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{1, 0, 3, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{0, 1, 0, 1, 4, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{0, 1, 1, 3, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{1, 0, 0, 0, 4, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{0, 1, 0, 0, 0, 0, 1, 5, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{0, 1, 0, 0, 1, 4, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{0, 0, 1, 0, 0, 0, 1, 0, 5, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 6, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0},
{0, 0, 0, 1, 0, 0, 0, 1, 1, 5, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{0, 0, 1, 0, 1, 0, 4, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 5, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{0, 0, 0, 1, 0, 1, 1, 4, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{1, 0, 0, 0, 3, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{0, 1, 0, 0, 0, 0, 1, 4, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{0, 1, 0, 0, 1, 3, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{0, 0, 1, 0, 0, 0, 1, 0, 4, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 5, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{0, 0, 0, 1, 0, 0, 0, 1, 1, 4, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{0, 0, 1, 0, 1, 0, 3, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 4, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{0, 0, 0, 1, 0, 1, 1, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}
};
