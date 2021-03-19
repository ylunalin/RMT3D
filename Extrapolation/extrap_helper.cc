#include "extrap_helper.hh"

/** Transfer the set of unique extrapolation points in the next layer
 * to the indices array. */
template <int N>
void extrap_helper<N>::update_list(bool verbose){
    watch.tic(0);
	unique_indices();
	//printf("extrap_helper::update_list(): Number of indices to consider %d\n", num);
	for(int i=0;i<num;i++) {
			indices[i] = list[i];
			if(verbose){
					int ii, jj, kk;
					unpack_index(indices[i], ii, jj ,kk);
					printf("extrap_helper:: nmask %u index %d (%d %d %d)\n", nmask, indices[i]-index(2,2,2), ii-2, jj-2, kk-2);
			}
	}
	list.clear();
    watch.toc(0);
}

/** Reduce the indices stored in list vector to a set of unique integers. */
template <int N>
void extrap_helper<N>::unique_indices(){
	std::sort(list.begin(), list.end());
	std::vector<int>::iterator last = std::unique(list.begin(), list.end());
	list.erase(last, list.end());
	num = list.size();
}

/** Using the coordinate data to calcuate a linear least square fit model,
 * and store the coefficients in coeff array. */
template <int N>
bool extrap_helper<N>::least_square(bool v){

    watch.tic(1);
	double *dp = data; //the position data, length per index controlled by capacity
	double *wdp = weight_data; // the weight for each data point at this index
	double *rp = rval; //the reference map values for
    if(capacity == 0) return false;

    //first compute the averages of each component
    double det;
    double rhs[rlen], a[nmat];
    for(int k=0;k<N;k++) rval_avg[k]=0;
    for(int k=0;k<D;k++) data_avg[k]=0;
    for(int k=0;k<rlen;k++) rhs[k]=0;
    for(int k=0;k<nmat;k++) a[k]=0;
    int wk=0;
    while(wk<capacity){
        // add_avg only modifies first argument
        add_avg<N>(rval_avg, rp, wdp);
        add_avg<D>(data_avg, dp, wdp);
        wk++;
        rp+=N;
        dp+=D;
        wdp++;
    }

    // divide out the weighted number of elements
    for(int k=0;k<N;k++) rval_avg[k] /= wt_capacity;
    for(int k=0;k<D;k++) data_avg[k] /= wt_capacity;
    // now gather the matrix entires
    rp = rval;
    dp = data;
    wdp = weight_data;
    wk=0;
    while(wk<capacity){
        // lin_sys only modifies first two argument
        lin_sys(a, rhs, rval_avg, data_avg, rp, dp, wdp);
        wk++;
        rp+=N;
        dp+=D;
        wdp++;
    }

    // now build the matrix and take its inverse
    // XXX need a templated mat class if we want to template on spatial dimension
    mat3 A(a[XX], a[XY], a[XZ], a[XY], a[YY], a[YZ], a[XZ], a[YZ], a[ZZ]);
    mat3 Ainv = A.inverse(det);

    if(fabs(det)<1e-17) {
        //printf("Extrap_helper: Oops. Determinant %g\n", det);
        //A.print();
        //puts("");
        //Ainv.print();
        return false;
    }
#if 0
    if(v) printf("rhs %g %g %g\n", rhs[0], rhs[1], rhs[2]);
#endif

    dot(Ainv, rhs); // rhs now stores the coefficients
    pack_coeff(rhs); // write coefficients to storage
    // now store the avg data and reference map values
    // for safety, we should check that cp[i] is at least 1, so there's enough room to store
    //write_avg<N>(r_avg,rval);
    //write_avg<D>(data_avg, data);

    // increment point for next index
    watch.toc(1);
    return true;
}

/** Compute averages. Used in least_square(), see that function see detail data organization. */
template <int N>
template <int n>
void extrap_helper<N>::add_avg(double (&avg)[n], const double *p, const double *w, bool verbose){
	const double Q = *w;
    //verbose = true;
	for(int i=0;i<n;i++){
		avg[i]+= p[i]*Q;
		//if(verbose) printf("p[%d] = %16.13f, Q = %g, avg[%d] = %16.13f\n", i,p[i], Q, i,avg[i]);
	}
}

/** Write averages. Used in least_square(), see that function see detail data organization.
 * The average data and reference map values are used to compute extrapolated ref values. */
template <int N>
template <int n>
void extrap_helper<N>::write_avg(double (avg)[n], double *p){
	for(int i=0;i<n;i++){
		p[i] = avg[i];
	}
}

/** Build the normal equation left and right hand sides, (A^T) A = (A^T) y. */
template <int N>
void extrap_helper<N>::lin_sys(double (&melems)[6], double (&rhs)[3*N], const double (&r_avg)[N], const double (&data_avg)[3], const double *rp, const double *dp, const double *wp, bool verbose){
	const double Q = *wp;

	for(int i=0;i<3;i++){
		melems[i] += Q * (dp[i] - data_avg[i])*(dp[i] - data_avg[i]); // xx, yy, zz
	}

	melems[3] += Q* (dp[0] - data_avg[0])*(dp[1] - data_avg[1]); // xy
	melems[4] += Q* (dp[0] - data_avg[0])*(dp[2] - data_avg[2]); // xz
	melems[5] += Q* (dp[1] - data_avg[1])*(dp[2] - data_avg[2]); // yz

	for(int i=0;i<N;i++){
		double local_rval = Q *(rp[i] - r_avg[i]);
		//if(verbose) printf("rval[%d] - rval_bar[%d] = (%16.13f -%16.13f) = %16.13f\n",i, i, rp[i], r_avg[i], local_rval);
		for(int j=0;j<3;j++) rhs[N*i+j] += (dp[j] - data_avg[j]) * local_rval; // rhs[0-2]: x*rval_x, y*rval_x, z*rval_x
	}

}

/** A very specifc dop product function for 3x3 reference map components. */
template <int N>
void extrap_helper<N>::dot(mat3 A, double (&vec)[3*N]){
	double prod[3*N];
	for(int i=0; i<N; i++){
		int a=i*3;
		int b=a+1, c=a+2;
		prod[a] = A.a11 * vec[a] + A.a12 * vec[b] + A.a13 * vec[c];
		prod[b] = A.a21 * vec[a] + A.a22 * vec[b] + A.a23 * vec[c];
		prod[c] = A.a31 * vec[a] + A.a32 * vec[b] + A.a33 * vec[c];
		/*
		printf("I=%d, %g = %g * %g + %g *%g + %g *%g \n"
		"%g = %g * %g + %g *%g + %g *%g \n"
		"%g = %g * %g + %g *%g + %g *%g \n",
		i,
		prod[a] , A.a11 , vec[a] , A.a12 , vec[b] , A.a13 , vec[c],
		prod[b] , A.a21 , vec[a] , A.a22 , vec[b] , A.a23 , vec[c],
		prod[c] , A.a31 , vec[a] , A.a32 , vec[b] , A.a33 , vec[c]);
		*/
	}
	for(int i=0;i<3*N;i++) vec[i] = prod[i];
}

/** Store the linear coefficients of the linear model to coeff array. */
template <int N>
void extrap_helper<N>::pack_coeff(const double (&tmp_coeff)[rlen]){
	for(int i=0;i<rlen;i++) coeff[i] = tmp_coeff[i];
}

/** Store the coordinate data (x,y,z) in data array.
 * Also add the weight associated with that point. */
template <int N>
void extrap_helper<N>::pack_data(const double (&lcoord)[3], const double data_wt) {
	for(int n=0;n<3;n++)  {
		data[3*capacity + n] = lcoord[n];
	}

	weight_data[capacity] = data_wt;
}

/** Store the reference map data (x,y,z) in rval array.
 * Used in gather_and_list() in fluid_3d.cc. */
template <int N>
void extrap_helper<N>::pack_rval(const double (&r)[N]){
	for(int i=0;i<N;i++) rval[N*capacity+i] = r[i];
}

/** Sum up the total weight as normalization constant.
 * Increment capacity counter.
 * Used in gather_and_list() in fluid_3d.cc. */
template <int N>
void extrap_helper<N>::add_capacity(const double data_wt) {
	wt_capacity += data_wt;
    capacity++;
}


/** Special treatment for the first layer. */
template <int N>
void extrap_helper<N>::add_first_layer(const int rule){
    watch.tic(0);
	increase_nmask();
	// when checking mask, we only need to look at inner layer of ghost node
	for(int k=1;k<so4-1;k++) for(int j=1;j<sn4-1;j++) for(int i=1;i<sm4-1;i++){
		int ind = index(i,j,k);
		if(mask[ind]==cur_mask){
			for(int kk=-1;kk<2;kk++) for(int jj=-1;jj<2;jj++) for(int ii=-1;ii<2;ii++){
				add_to_next_list(ind, ii, jj, kk, rule);
			}
		}
	}
    watch.toc(0);

	comm_mask();
	scan_ghost_region();

	update_list();
}

/** Set up tables an pointers needed for communication. */
template <int N>
void extrap_helper<N>::setup_table(){
	mpi_table t(grid);
	int len = t.setup_table(grid);
	cbuf.check_buf_uint(len);
	t.copy_table(comm_table);
}

template <int N>
void extrap_helper<N>::setup_comm(){
	rank = grid->rank;
	setup_table();
	grid->set_up_neighbors(neigh = new int[27]);
	int n_neigh = 0;
	// count the number of neighbors, i=13 corresponds to oneself
	for(int i=0;i<27;i++){
		if(i!=13 && neigh[i]>=0) n_neigh ++;
	}
	if(rank == 0 && grid->procs > n_neigh) n_neigh = grid->procs;
	reqs = new MPI_Request[2*n_neigh];
	stats = new MPI_Status[2*n_neigh];
	// if we want to be able to print masks, we need to check buf
	int m=grid->m+4, n=grid->n+4, o=grid->o+4;
	cbuf.check_buf_uint(max_of_three(sm4*sn4, sm4*so4, sn4*so4));
	if(rank==0){
		cbuf.check_buf_uint(max_of_three(m*n, m*o, n*o));
	}
}

/** A working but primitive version of communication function.
 * Is there a way to piggy bag on fluid_3d communications?
 * But that's not ideal if we want a clean separation from fluid_3d core code. */
template <int N>
void extrap_helper<N>::comm_mask(){
    watch.tic(3);
	int ops = 0;
	unsigned int *buf  = reinterpret_cast<unsigned int *> (cbuf.buf);
	unsigned int *send_buf, *recv_buf;
	MPI_Request *req = reqs;

	for(int i=0;i<27;i++){
		if(i!=13 && neigh[i] != -1){
			send_buf = buf;
			comm_info cur_neigh = comm_table[i];
			// send tag is the same as neighbor tag, i.e. i
			int s_tag = cur_neigh.s_tag;
			int len = cur_neigh.len;
			int r_tag = cur_neigh.r_tag;
			// copy mask in the sending region to buffer, note buf has been incremented by cur_neigh.len
			copy_to_buf(s_tag, buf);
			// MPI nonblocking send (void *buf, int count, MPI_DATATYPE, int dest, int tag, MPI_Comm, MPI_Request *)
			MPI_Isend(send_buf, len, MPI_UNSIGNED, neigh[i], s_tag, grid->cart, req++);
			recv_buf = buf;
			MPI_Irecv(recv_buf, len, MPI_UNSIGNED, neigh[i], r_tag, grid->cart, req++);
			buf += len;
			ops += 2;
		}
	}
	MPI_Waitall(ops,reqs,stats);

	// After all communication is over, we put mask values into ghost regions
	buf = reinterpret_cast<unsigned int *> (cbuf.buf);
	for(int i=0;i<27;i++){
		if(i!=13 && neigh[i]!=-1){
			comm_info cur_neigh = comm_table[i];
			int len = cur_neigh.len;
			buf += len;
			//use neighbor tag i here because copy_from_buf reads off comm_info.g0 instead of r0
			copy_from_buf(i, buf);
		}
	}
    watch.toc(3);
}

// After masks are communicated, we add potential extrapolation points in the ghost region to list vector.
template< int N>
void extrap_helper<N>::scan_ghost_region(){
    watch.tic(3);
	for(int i=0;i<27;i++){
		if(i!=13){
			comm_info cur_neigh = comm_table[i];
			int start = cur_neigh.g0;
			unsigned int * ghost_mask = mask+start;
			int m=cur_neigh.m;
			int n=cur_neigh.n;
			int o=cur_neigh.o;
			for(int kk=0;kk<o;kk++) for(int jj=0;jj<n;jj++) for(int ii=0;ii<m;ii++){
				int ind = index(ii, jj, kk);
				/* int p, q, r;
				unpack_index(ind+start, p, q, r);
				if(rank==2 && i ==12) printf("ghost_mask[%d] (%d %d %d) = %u, nmask=%u\n", ind+start, p, q, r, ghost_mask[ind], nmask);
				*/
				if(ghost_mask[ind] == nmask) {
					list.push_back(start+ind);
				}
			}
		}
	}
    watch.toc(3);
}

template <int N>
void extrap_helper<N>::copy_to_buf(int neigh_tag, unsigned int *(&b)){
	comm_info t = comm_table[neigh_tag];
	unsigned int * send_m = mask + t.s0;
	for(int k=0;k<t.o;k++) for(int j=0;j<t.n;j++) for(int i=0;i<t.m;i++){
		int step  = index(i,j,k);
		*(b++) = send_m[step];
		/*int ii, jj, kk;
		unpack_index(step+t.s0, ii, jj, kk);
		//if(nmask ==7 && kk == 52 && rank==1) printf("rank %d send to neighbor %d  mask[%d] (%d %d %d) = %u\n",rank, neigh_tag, step+t.s0, ii, jj, kk, send_m[step]);
		*/
	}
}

template <int N>
void extrap_helper<N>::copy_from_buf(int neigh_tag, unsigned int *(&b)){
	comm_info t = comm_table[neigh_tag];
	unsigned int * recv_m = mask + t.g0;
	for(int k=0;k<t.o;k++) for(int j=0;j<t.n;j++) for(int i=0;i<t.m;i++){
		int step  = index(i,j,k);
		recv_m[step] = *(b++);
		/*int ii, jj, kk;
		unpack_index(step+t.g0, ii, jj, kk);
		if(nmask == 7 && kk==52 && rank==2) printf("rank %d receive from %d mask[%d] (%d %d %d) = %u\n", rank, neigh_tag, step+t.g0, ii, jj, kk, recv_m[step]);
		*/
	}
}

template <int N>
void extrap_helper<N>::dump_mask(){
	char *fn = new char[256];
	sprintf(fn, "new/mask%d_%d", nmask, rank);
	FILE *fm = p_safe_fopen(fn, "w");
	for(int j=0;j<sn4;j++){
			for(int i=0; i<sm4;i++) {
				fprintf(fm, "%5u ", mask [index(i, j, 52)]);
			}
			fputc('\n', fm);
	}
	fclose(fm);
	delete [] fn;
}

// explicit instantiation
template class extrap_helper<3>;
template class extrap_helper<1>;
