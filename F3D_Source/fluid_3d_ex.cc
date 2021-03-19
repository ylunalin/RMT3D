/*
A file that contains the extrapolation quadrant of the fluid_3d galaxity.
*/
#include "fluid_3d.hh"

// ############################### Extrapolation functions ###############################//

/** Gives the neighbor rules. 0 orthogonal; 1 sqrt 2; 2 sqrt 3 */
int fluid_3d::rule(){
	return 0;
}

/** Extrapolate for the first time.
 * The old set of routines are housed in fluid_3d_ex_old.cc.
 */
void fluid_3d::init_extrapolate(bool verbose){
    extrapolate(verbose);
    watch.tic(0);
	for (int k=0;k<so4;k++) for(int j=0;j<sn4;j++) for(int i=0;i<sm4;i++){
		int ind = index(i, j, k);
		for (int l=0;l<extraps.n[ind];l++){
			extraps.f[ind][l].set_full_timestep();
		}
	}
    // We blend the rho to set an updated field, based on extrapolated field values
    watch.toc(0);
}


/** Extrapolate object by object. */
void fluid_3d::extrapolate(bool verbose){
    const int comm_xpred = 1;
	for (int i=1;i<=mgmt->n_obj; i++){
		// set the index of layers to zero for every object
		n_layers=0;
		first_layer(i);

		for (n_layers=1;n_layers<=max_layers;n_layers++) {
			extrap_a_layer(i, verbose);
			// note layerID is layer index+1
            unsigned int c = (n_layers<<ldigits)+i;

            communicate_a_layer(c, comm_xpred);
		}
	}

}

/** Sweep the grid for the first time, and add appropraite mask to the extrap helper class;
 * extrapo helper then take the mask values and create the first layer of points,
 * for a particular object. */
void fluid_3d::first_layer(unsigned int oid){
	// NOTE
	// Using index() and index_unpack() wisely
	// Since reference map is communicated before extrapolation is called
	// we don't need to communicate mask after marking primary
    watch.tic(0);
	expper.reset(oid);
    expper.watch.tic(0);
	for(int k=0;k<so4;k++) for (int j=0;j<sn4;j++) for(int i=0;i<sm4;i++){
		int ind = index(i,j,k);
		ref_map &rf = rm_mem[ind];
		if(rf.c == oid ){
			expper.mark_primary(ind);
		}
	}
    expper.watch.toc(0);
	expper.add_first_layer(rule());
    watch.toc(0);
}

/** Extrapolate current layer and gather points for next layer. */
void fluid_3d::extrap_a_layer(unsigned int oid, bool verbose){
    // NOTE we can assme that n_layers is the current layer being extrapolated
    watch.tic(0);
	gather_and_list(oid, verbose);

	// we need to update_list() last, because we'll use indices array up until transfer(oid)
	expper.update_list();
    watch.toc(0);
}

/** Gather a list of grid points to consider in extrapolation routine.
 * \params[in] oid the id of the object currently being considered.
 */
void fluid_3d::gather_and_list(unsigned int oid, bool verbose){
	int * ip = expper.indices;

	expper.increase_nmask();
	for(int j=0;j<expper.num;j++){
		// center grid point, index wrt first ghost node
		int ind = ip[j];

        if(expper.within_extrapolation_domain(ind)){

            int ii,jj,kk;
            unpack_index(ind, ii, jj, kk);

            // checking 5x5x5 cube to gather data for extrapolation
            expper.watch.tic(4);
            // Get normal vector (from previous extrapolation, if available) at ind
            // since grad phi is calculated on the faces, we average over all six faces
            int num_faces = 0;
            double grad_phi[3] = {0,0,0};
            double grad_phi_norm = 0;
            // This if condition here guards that the reference map pointers are only
            // queried if they are within the non-ghost simulation domain
            if(ii >=2 && ii <= sm4-3 && jj >=2 && jj <= sn4-3 && kk >= 2 && kk <= so4-3 ) {
                ref_map *prev[4];
                prev[0] = get_refmap_prev(ind, oid);
                prev[1] = get_refmap_prev(ind + 1, oid);
                prev[2] = get_refmap_prev(ind + index(0,1,0), oid);
                prev[3] = get_refmap_prev(ind + index(0,0,1), oid);

                for(int pp=0;pp<4;pp++) {
                    if (prev[pp]!=NULL) {
                        // we average all six faces, 3 faces belong to this reference map extrap
                        if(pp==0) {
                            for(int dd = 0; dd < 3; dd++) {
                                if(prev[pp]->grad_phi[dd][dd] != -1.e200 ) {
                                    bool add = false;
                                    for(int ff = 0; ff < 3; ff++) {
                                        if(prev[pp]->grad_phi[dd][ff] != 0.0) { add = true; }
                                        grad_phi[ff] += prev[pp]->grad_phi[dd][ff];
                                    }
                                    if(add) {num_faces ++;}
                                } // end if != -1e200
                            }
                        } else {
                            // If we found a neighboring extrapolation
                            // we see if grad phi is calcuated there
                            if(prev[pp]->grad_phi[pp-1][pp-1] != -1.e200) {
                                bool add = false;
                                for(int ff=0;ff<3;ff++) {
                                    if(prev[pp]->grad_phi[pp-1][ff] != 0.0) { add = true; }
                                    grad_phi[ff] += prev[pp]->grad_phi[pp-1][ff];
                                }
                                if(add) {num_faces ++;}
                            } // end if != -1e200
                        }
                    }
                }

                // If we have any grad_phi info, we compute the average
                // If we don't, we set grad_phi = 0
                if(num_faces==6) {
                    for(int pp=0;pp<3;pp++) {
                        grad_phi[pp]/=num_faces;
                        grad_phi_norm += (grad_phi[pp] * grad_phi[pp]);
                    }
                    if(grad_phi_norm!=0) {
                        // don't forget to square root a norm again...
                        grad_phi_norm = sqrt(grad_phi_norm);
                        for(int pp=0;pp<3;pp++) {
                            grad_phi[pp]/=grad_phi_norm;
                        }
                    } else {
                        for(int pp=0;pp<3;pp++) {
                            grad_phi[pp]=0;
                        }
                    }
                } else {
                    for(int pp=0;pp<3;pp++) {
                        grad_phi[pp]=0;
                    }
                }
            }

            int search_r = 2;
            bool success = false;
#if defined(DEBUG)
            int zero_wt_data_pt=0;
            int lowp=100, lowq=100, lowr=100, hip=-100, hiq=-100, hir=-100;
#endif
            while(!success && search_r <= max_extrap_rs) {
                expper.zero_capacity();
#if defined(DEBUG)
                zero_wt_data_pt=0;
                lowp=100, lowq=100, lowr=100, hip=-100, hiq=-100, hir=-100;
#endif

                for(int r = -search_r; r<=search_r; r++) for(int q = -search_r; q<=search_r; q++) for(int p = -search_r; p<=search_r; p++){
                    // Non periodic boundary: don't need to worry about the periodic wrap around and can indeed
                    // extrapolate all the way into the ghost region, but when we get the data points
                    // we must be careful to not go beyond the bounds
                    // Since we have checked the bounds, get_refmap() is safe to use
                    if(ii+p > sm4-1 || ii+p < 0) continue;
                    if(jj+q > sn4-1 || jj+q < 0) continue;
                    if(kk+r > so4-1 || kk+r < 0) continue;
                    int step = index(p,q,r);
                    ref_map *rf = get_refmap(ind+step, oid);
                    // only if we find a reference map at index=candidate
                    // that matches the same object id, in lower layers
                    if(rf!=NULL && rf->lid()<n_layers) {
                        double lcoord[3] = {lx[ii+p], ly[jj+q], lz[kk+r]};

                        if(grad_phi_norm != 0. && n_layers<=2) {
                            // get the physical vector here, from this reference map to the current extrapolation point
                            // calculate a weight and pack it into the extrapolation routine
                            double extrap_vec[3] = {lx[ii]-lcoord[0], ly[jj]-lcoord[1], lz[kk]-lcoord[2]};
                            double data_wt = 0., extrap_vec_norm = 0.;
                            for(int pp=0;pp<3;pp++) {
                                data_wt += extrap_vec[pp] * grad_phi[pp];
                                extrap_vec_norm += extrap_vec[pp] * extrap_vec[pp];
                            }
                            extrap_vec_norm = sqrt(extrap_vec_norm);
                            data_wt /= extrap_vec_norm;
                            data_wt *= pow(expper.wt, abs(p) + abs(q) + abs(r));

                            if(data_wt > 0) {
                                expper.pack_data(lcoord,data_wt);
                                expper.pack_rval(rf->xpred);
                                expper.add_capacity(data_wt);
#if defined(DEBUG)
                                if (p<lowp) {lowp = p;}
                                else if (p>hip){hip = p;}
                                if (q<lowq) {lowq = q;}
                                else if (q>hiq){hiq = q;}
                                if (r<lowr) {lowr = r;}
                                else if (r>hir){hir = r;}
#endif
                            }
#if defined(DEBUG)
                            else{
                                zero_wt_data_pt ++;
                            }
#endif
                        }
                        else {
                            double data_wt = pow(expper.wt, abs(p) + abs(q) + abs(r));
                            expper.pack_data(lcoord,data_wt);
                            expper.pack_rval(rf->xpred);
                            expper.add_capacity(data_wt);
#if defined(DEBUG)
                            if (p<lowp) {lowp = p;}
                            else if (p>hip){hip = p;}
                            if (q<lowq) {lowq = q;}
                            else if (q>hiq){hiq = q;}
                            if (r<lowr) {lowr = r;}
                            else if (r>hir){hir = r;}
#endif
                        }
                    }
                }// end of 5x5x5 cube

                expper.watch.toc(4);

                bool v=false;
                success = expper.least_square(v);
                search_r += 1;
            }

            if(!success) {
                printf("\n\nWarning: rank %d global index (%d %d %d) [%d] should be extrapolated in %d-th layer\n"
                       "But that didn't happen due to no eligible points or degenerate matrix!\n\n",
                        rank, ai+ii-2, aj+jj-2, ak+kk-2, ind, n_layers);
                printf("Total weight to normalize: %g, number of points: %d\n", expper.wt_capacity, expper.capacity);

                double rval_avg[3];
                expper.get_rval_avg(rval_avg);
                printf("Average reference map value for the least square problem: %g %g %g\n", rval_avg[0], rval_avg[1], rval_avg[2]);
                printf("Search radius: %d\n", search_r);

#if defined(DEBUG)
                printf("Data points within search radius but non-positive weights %d\n", zero_wt_data_pt);
                printf("Range: (%d, %d) x (%d, %d) x (%d, %d)\n", lowp, hip, lowq, hiq, lowr, hir);
#endif
                p_fatal_error("Extrapolation failed. Oops", 1);
            }

            transfer(oid, ind);

        }// end if within_domain

        // looking in 3x3x3 cube for neighboring points for next layer
        // since we are not communicating indices, we need to keep around indices that are in the ghost region.
        if(n_layers < max_layers && expper.within_inner_ghost(ind)){
            expper.watch.tic(0);
            for(int r=-1;r<2;r++) for(int q=-1;q<2;q++) for(int p=-1;p<2;p++){
                // third argument 0 means use orthogonal neighbor rule
                expper.add_to_next_list(ind, p, q, r, rule(), verbose);
            }// end of 3x3x3 cube
            expper.watch.toc(0);
        }// end if n_layers < max_layers

	}// end of all indices
	if(n_layers<max_layers){
		expper.comm_mask();
		expper.scan_ghost_region();
	}
}

void fluid_3d::transfer(unsigned int oid, int ind){
    if(expper.capacity == 0) return;

    expper.watch.tic(4);

    int i,j,k;
    unpack_index(ind, i, j, k);
    double refmap_x [3], pos_subtract_avg[3], rval_avg[3];
    expper.subtract_data_avg(lx[i], ly[j], lz[k], pos_subtract_avg);
    expper.get_rval_avg(rval_avg);
	double * ep = expper.coeff;

    for(int ii=0;ii<3;ii++){
        refmap_x[ii] = rval_avg[ii];
        for(int jj=0;jj<3;jj++) refmap_x[ii] += ep[jj]*pos_subtract_avg[jj];
        ep+=3;
    }
    // n_layers starts at 1 for first layer (not zero)
    if(!(std::isnan(refmap_x[0]) || std::isnan(refmap_x[1]) || std::isnan(refmap_x[2]))){
        extraps.add_mem(ind, oid+(n_layers<<ldigits), refmap_x, 1);
    }
#if 0
        int eidtmp = index(3, 31, 33);
        if(ind == eidtmp) {
            printf("pos_subtract_avg %g %g %g\n", pos_subtract_avg[0], pos_subtract_avg[1], pos_subtract_avg[2] );
            printf("coefficient %g %g %g\n", ep[0], ep[1], ep[2]);
            printf("rval_avg%g %g %g\n", rval_avg[0], rval_avg[1], rval_avg[2]);
            printf("Index %d number of extrapolation at missing neighbor is %d\n", ind, extraps.n[eidtmp]);
            for (int tmpi = 0 ; tmpi<extraps.n[eidtmp]; tmpi ++){
                printf("extrapolation in layer %d\n", extraps.f[eidtmp][tmpi].lid());
            }
        }
#endif

    expper.watch.toc(4);
}

/* ###################### SHARED BETWEEN NEW AND OLD ########################*/

/** Mark all the cells on the primary grid inside cells if phi<0.
 * Then transfer those values to primary grid.
 * Similarly, mark those cell with phi>0 outside, change it's flag to zero.
 */
void fluid_3d::redraw_boundary(bool verbose){
	if(verbose) printf("Rank %d. Redrawing boundary.\n",rank);
	int lflag = (max_layers+1)<<ldigits;
	for (int k = -2; k < so+2; k++) {
		for (int j = -2; j < sn+2; j++) {
			for (int i = -2; i < sm+2; i++) {
				int ind = i + sm4*j + smn4*k;
				bool occupied = false;
				// need to check primary grid first
				ref_map &prime = rm0[ind];
                // reset fresh_in and fresh_out for next timestep use

				int obj_id = prime.oid();
				// values on the primary grid can only go from being inside to being outside
				if( obj_id > 0) {
                    double pphi = prime.phi(mgmt);
                    if (pphi>0) {
                        prime.c = lflag|prime.c;
						extraps.add(ind, &prime);
                        prime.clear();
					} else occupied = true;
				}

				// now check all the extrapolated maps at that cell
				int N = extraps.n0[ind];
				for(int l=0;l<N;l++) {
					ref_map &sf = extraps.f0[ind][l];
					// it's guaranteeed that obj_id won't be zero
					obj_id = sf.oid();
					// check if this extrapolated value belongs to inside of an object
					// only allow grid points that lie within the transition zone to move inside
					//bool inside = (sf.phi(mgmt) < 0 && sf.lid()<mgmt->wt_n);
                    // HACK HACK IS IT A HACK?
                    // Only allow the first layer to move inside
					bool inside = (sf.phi(mgmt) < 0 && sf.lid()==1);

					if (inside) {
						if((occupied && sf.phi(mgmt) < prime.phi(mgmt)) || !occupied) {
							if(prime.phi(mgmt) < -dx) {
								printf("redraw_boundary(): rank %d index (%d %d %d) global index (%d %d %d) prime obj %d and extrap obj %d\n",
										rank, i, j, k, ai+i, aj+j, ak+k, prime.oid(), sf.oid());
								p_fatal_error("redraw_boundary(): Two deep maps try to be inside a cell at the same time!\n", 1);
							}
							// move it to primary grid
							prime = extraps.f0[ind][l];
							prime.c = obj_id;
							occupied = true;
							extraps.rm(ind,obj_id);
						}

					}
				} // done looping over all extrapolation maps
			}
		}
	} // done for loops

    const int comm_x=0;
	communicate_a_layer((max_layers+1)<<ldigits,comm_x);
}
