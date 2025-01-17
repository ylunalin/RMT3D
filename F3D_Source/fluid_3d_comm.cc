#include "fluid_3d.hh"

// ######################## SETUP ######################

const int fluid_3d::packlens[12] = {3,3,18,1,3,3,18,1,3,3,18,3};

/**
 * Set up pieces needed for communication, i.e. buffer lookup table, neighbor table
 * MPI_Request, MPI_Status objects, and comm_buffer object. */
void fluid_3d::setup_communication() {

	// setup lookup and neighbor tables
	setup_tables();
	grid->set_up_neighbors(neigh = new int[27]);

	// get number of neighbors. depends if normal proc (just count)
	int n_neigh = 0;
	for (int tag = 0; tag < 27; tag++)
		if (tag != 13 && neigh[tag] >= 0) n_neigh++;

	// or if master, which can communicate with everybody
	if (rank == 0 && grid->procs > n_neigh) n_neigh = grid->procs;

	// need MPI objects for each possible communication (2-way, each neighbor)
	reqs = new MPI_Request[2*n_neigh];
	stats = new MPI_Status[2*n_neigh];

	// check buffer lengths for each slice portion (1-way, norm -> master)
	cbuf.check_buf((sm+1)*(sn+1));
	cbuf.check_buf((sm+1)*(so+1));
	cbuf.check_buf((sn+1)*(so+1));

	// and, if master, check for each global slice portion (1-way)
	if (rank == 0) {
		cbuf.check_buf((m+1) * (n+1));
		cbuf.check_buf((m+1) * (o+1));
		cbuf.check_buf((n+1) * (o+1));
	}
}

/**
 * Cleanup function for communication objects
 */
void fluid_3d::cleanup_communication() {
	delete[] stats;
	delete[] reqs;
	delete[] neigh;
}

/** A version of the function setup_table, but the setting-up is done in mpi_table class. */
void fluid_3d::setup_tables(){
	mpi_table t(grid);
	// setup comm_table for the main grid
	int len = t.setup_table(grid);
	cbuf.check_buf(len);
	t.copy_table(comm_table);

	// setup corner table for pressure communication
	len = t.setup_table(&fem_grid);
	cbuf.check_buf(len);
	t.copy_table(corner_table);
}

//################# MPI CALLS ###################

/**
 * Flags:
 *  1: fluid velocities, cell center (3)
 *  2: pressure and auxiliary pressure, and change in pressure (3)
 *  4: fluid velocities, cell faces (18)
 *  8: object id (1)
 * 16: ref map components, cell center (3)
 * 32: ref map pred components, center (3)
 * 64: ref map components, cell faces (18)
 */
template<int flags>
void fluid_3d::communicate(int lint) {

    watch.tic(2);
	comm_info *ct = (flags==2)?corner_table:comm_table;

	// how many doubles are communicated per node?
	int f = 1,per_n=0,per_e=0,per_dvel=0;

	// per primary grid node...
	for (int i = 0;i < FFLAG_N+SFLAG_N;i++,f<<=1) {
		if (flags&f) per_n += packlens[i];
	}

	// ...and per extrapolation
	for (int i = FFLAG_N+SFLAG_N;i<FFLAG_N+SFLAG_N+EFLAG_N;i++,f<<=1) {
		if (flags&f) {
            per_e += packlens[i];
        }
	}
    
    // for dvel
    if(flags & (1<< (FFLAG_N+SFLAG_N+EFLAG_N))) {
        per_dvel = packlens[11];
        //printf("rank %d flag %d Communication communicating dvel\n", rank, flags);
    }

	// if we're communicating extraps,
	// need additional comms for lengths
	int len[27] = {};
	int *rlen;
	MPI_Request *lreq,*lr;
	MPI_Status *lstat;
	if (flags&EFLAG) {
		rlen = new int[27];
		lr=lreq = new MPI_Request[54];
		lstat = new MPI_Status[54];
	}

	// count operations
	int ops=0,lops=0;

	// get start of buffer and requests array
	double *buf = cbuf.buf,*sbuf,*rbuf;
	MPI_Request *req = reqs;

	// DO ALL SENDING
	for(int i =0; i < 27; i++) if(i!=13&&neigh[i]!=-1) {

		// get current spot in buffer and copy in, getting length of comm
		sbuf = buf;
        // get lenght of comm buff if dvel is also communicated
		len[i] = (per_n+per_dvel)*ct[i].len;

		// get num of extraps packed
		int ne_packed = copy_to_buffer<flags>(ct[i].s_tag,buf,lint,ct);
		if (flags&EFLAG) len[i] += per_e*ne_packed + ct[i].len;

		// send away
		MPI_Isend(sbuf,len[i],MPI_DOUBLE,neigh[i],ct[i].s_tag,grid->cart,req++);

		// count comm
		ops++;

		// if we're doing extrapolations, communicate the size of the comm
		if (flags&EFLAG) {

			// send length
			MPI_Isend(&len[i],1,MPI_INT,neigh[i],ct[i].s_tag+27,
				grid->cart,lr++);
			// recv length
			MPI_Irecv(rlen+i,1,MPI_INT,neigh[i],ct[i].r_tag+27,
				grid->cart,lr++);
			// count both ops
			lops += 2;
		}
	}

	// IF EXTRAPS, WAIT FOR LENGTH COMMS
	if (flags&EFLAG) MPI_Waitall(lops,lreq,lstat);
    //puts("HERE1");

	// record start of receiving buffer
	rbuf = buf;

	// DO ALL RECEIVING
	for(int i =0; i < 27; i++) if(i!=13&&neigh[i]!=-1) {

		// get length of this comm
        len[i] = (per_n + per_dvel)*ct[i].len;
		if (flags&EFLAG) len[i] = rlen[i];

		// receive it
		MPI_Irecv(buf,len[i],MPI_DOUBLE,neigh[i],ct[i].r_tag,
				grid->cart,req++);

		// update place in buffer and count op
		buf += len[i];
		ops++;
	}

	// WAIT FOR ACTUAL COMMS
	MPI_Waitall(ops,reqs,stats);
    //puts("HERE2");

	// EMPTY BUFFER
	buf = rbuf;
	for(int i =0; i<27; i++) if (i!=13&&neigh[i]!= -1) {
		copy_from_buffer<flags>(ct[i].r_tag,buf,ct);
	}

	// SHIFT REFMAP
	//DEBUG
	//if (flags&(SFLAG|EFLAG)) shift_refmap<flags>(lint,0);

	// clean up extraps tools
	if (flags&EFLAG) {
		delete[] rlen;
		delete[] lreq;
		delete[] lstat;
	}
    watch.toc(2);
}

/** Communicate a single layer of extrapolated reference
 * map values. Extrapolation done by least square model.
 */
void fluid_3d::communicate_a_layer(unsigned int c,int type) {

    watch.tic(2);
	int lops = 0,ops = 0;
	int rlen[27] = {}, len[27]={};
	comm_info (&ct)[27] = comm_table;
	double *buf = cbuf.buf;
	MPI_Request* req = reqs;
	MPI_Request lreq[54];
	MPI_Status lstat[54];
	MPI_Request *lr = lreq;

	// variables to update below
	double *send_buf;

	// go through all neighbors
	for(int i =0; i < 27; i++) {

		// if (not self) and a neighbor exists
		if (i != 13 && neigh[i] != -1) {

			// get current spots in buf, copy, and then send
			send_buf = buf;
			len[i] = copy_layer_to_buffer(ct[i].s_tag,c,buf,type);
			MPI_Isend(send_buf,len[i],
				MPI_DOUBLE,neigh[i],ct[i].s_tag,
				grid->cart,req++);

			// also send length
			MPI_Isend(&len[i],1,MPI_INT,neigh[i],ct[i].s_tag+27,
				grid->cart,lr++);
			// recv length
			MPI_Irecv(rlen+i,1,MPI_INT,neigh[i],ct[i].r_tag+27,
				grid->cart,lr++);

			ops++;
			lops += 2;
        } else if (i != 13 && type==0 && &ct != &comm_table) {
			//TODO  what do bc's look for ref_map / extraps??
		}
    } // let all lengths get delivered
	MPI_Waitall(lops,lreq,lstat);

	double *rbuf = buf; // go back here to empty buffer

	// go through all neighbors
	for(int i =0; i < 27; i++) {

		// if (not self) and a neighbor exists
		if (i != 13 && neigh[i] != -1) {

			// now receive and copy out of receiving buffer
			MPI_Irecv(buf,rlen[i],
				MPI_DOUBLE,neigh[i],ct[i].r_tag,
				grid->cart,req++);
			buf += rlen[i];
			ops++;

		} else if (i != 13 && type==0 && &ct != &comm_table) {
			// TODO same q as above...what are bcs?
		}
	}

	// let everything go above and wait to finish at end
	MPI_Waitall(ops,reqs,stats);

	// now, empty everything received into spot out of buffer
	buf = rbuf;
	for(int i =0; i<27; i++){
		if (i != 13 && neigh[i] != -1) {
			copy_extraps_from_buffer(ct[i].r_tag,buf,type);
		}
	}

	//DEBUG
	//shift_extraps(0,layer);
    watch.toc(2);
}

/** Communicates the extrapolated reference map in the
 * extrap_mmap object.
 */
void fluid_3d::communicate_extrap_fx() {

    watch.tic(2);
	int lops = 0,ops = 0;
	int rlen[27]={}, len[27]={};
	comm_info (&ct)[27] = comm_table;
	double *buf = cbuf.buf;
	MPI_Request *req = reqs;
	MPI_Request lreq[54];
	MPI_Request *lr = lreq;
	MPI_Status lstat[54];

	// variables to update below
	double *send_buf;

	// go through all neighbors
	for(int i =0; i < 27; i++) {

		// if (not self) and a neighbor exists
		if (i != 13 && neigh[i] != -1) {

			// get current spots in buf, copy, and then send
			send_buf = buf;
			len[i] = copy_extraps_to_buffer(ct[i].s_tag,buf,2);
			MPI_Isend(send_buf,len[i],
				MPI_DOUBLE,neigh[i],ct[i].s_tag,
				grid->cart,req++);

			// also send length
			MPI_Isend(&len[i],1,MPI_INT,neigh[i],ct[i].s_tag+27,
				grid->cart,lr++);
			// recv length
			MPI_Irecv(rlen+i,1,MPI_INT,neigh[i],ct[i].r_tag+27,
				grid->cart,lr++);

			ops++;
			lops += 2;

		} else if (i != 13 && &ct != &comm_table) {
			//TODO  what do bc's look like for ref_map / extraps??
		}
	}

	// let all lengths get delivered
	MPI_Waitall(lops,lreq,lstat);

	double *rbuf = buf; // go back here to empty buffer

	// go through all neighbors
	for(int i =0; i < 27; i++) {

		// if (not self) and a neighbor exists
		if (i != 13 && neigh[i] != -1) {

			// now receive and copy out of receiving buffer
			MPI_Irecv(buf,rlen[i],
				MPI_DOUBLE,neigh[i],ct[i].r_tag,
				grid->cart,req++);
			buf += rlen[i];
			ops++;

		} else if (i != 13  && &ct != &comm_table) {
			// TODO same q as above...what are bcs?
		}
	}

	// let everything go above and wait to finish at end
	MPI_Waitall(ops,reqs,stats);

	// now, empty everything received into spot out of buffer
	buf = rbuf;
	for(int i =0; i<27; i++){
		if (i != 13 && neigh[i] != -1) {
			copy_extraps_from_buffer(ct[i].r_tag,buf,2);
		}
	}

	// DEBUG
	//shift_extraps(1);
    watch.toc(2);
}

// ############ UTILITIES #############

/** Shifting the reference map with periodic boundary condition.
*/
template<int flags>
void fluid_3d::shift_refmap(int lint,double dummy){
	// type = 0 -> x
	// type = 1 -> fx
	// type = 3 -> xpred

	// orthogonal neighbors in x, y, z direction
	int ghost_m=0, ghost_n=0, ghost_o=0;
	bool prds[3] = {grid->x_prd,grid->y_prd,grid->z_prd};
	bool to_shift[6] = {ai==0,bi==m,aj==0,bj==n,ak==0,bk==o};
	int dels[3] = {1,3,9};
	double shift_by[6] = {-mgmt->lx,mgmt->lx,-mgmt->ly,mgmt->ly,-mgmt->lz,mgmt->lz};

	for (int face = 0; face < 6; face++) {
		if (prds[face/2] && to_shift[face]) {

			// need to do each ghost region on this face,
			// on this proc
			for (int c1 = 0; c1 < 3; c1++) {
				for (int c2 = 0; c2 < 3; c2++) {
					int reg = c1*dels[(face/2 + 1)%3]
						+ c2*dels[(face/2 + 2)%3]
						+ 2*(face%2)*dels[face/2];

					// get ghost region info
					ghost_m = comm_table[reg].m;
					ghost_n = comm_table[reg].n;
					ghost_o = comm_table[reg].o;
					ref_map *rf = rm_mem + comm_table[reg].g0;

					// march through ghost region and shift
					for(int k=0;k<ghost_o;k++)
						for(int j=0;j<ghost_n;j++)
							for(int i=0;i<ghost_m;i++) {
						int eid = i+j*sm4+k*smn4;

						// in primary grid, only shift if object there
						if ((flags&SFLAG) && rf[eid].oid() > 0) {
							rf[eid].shift<flags&SFLAG>(face,shift_by[face]);
						}

						// shift all extraps matching layer code
						if (flags&EFLAG) {
							int egid = eid+comm_table[reg].g0;
							for (int ii=0; ii<extraps.n[egid]; ii++) {
								if (lint<1||extraps.f[egid][ii].lid()==lint) {
									extraps.f[egid][ii].shift<flags&EFLAG>(
											face,shift_by[face]);
								}
							}
						}
					}
				}
			}
		}
	}
}

/** Shifting the reference map according to periodic boundary condition.
 * The values being shifted are stored in the extrap_mmap object. */
void fluid_3d::shift_extraps(int type,int lid){
	// type = 0 -> x
	// type = 1 -> fx
	// type = 3 -> xpred

	// orthogonal neighbors in x, y, z direction
	int ghost_m=0, ghost_n=0, ghost_o=0;
	bool prds[3] = {grid->x_prd,grid->y_prd,grid->z_prd};
	bool to_shift[6] = {ai==0,bi==m,aj==0,bj==n,ak==0,bk==o};
	int dels[3] = {1,3,9};
	double shift_by[6] = {-mgmt->lx,mgmt->lx,-mgmt->ly,mgmt->ly,-mgmt->lz,mgmt->lz};

	for (int face = 0; face < 6; face++) {
		if (prds[face/2] && to_shift[face])
			for (int c1 = 0; c1 < 3; c1++)
				for (int c2 = 0; c2 < 3; c2++) {
					int reg = c1*dels[(face/2 + 1)%3]
						+ c2*dels[(face/2 + 2)%3]
						+ 2*(face%2)*dels[face/2];

					ghost_m = comm_table[reg].m;
					ghost_n = comm_table[reg].n;
					ghost_o = comm_table[reg].o;
					int i0 = comm_table[reg].g0;
					ref_map *rm;
					for(int k=0;k<ghost_o;k++) for(int j=0;j<ghost_n;j++) for(int i=0;i<ghost_m;i++) {
						// only do extraps!
						/*
						if (rf[i+j*sm4+k*smn4].oid() > 0 && rf[i+j*sm4+k*smn4].lid() == lid) {
							rm = rf + i + j*sm4 + k*smn4;
							if(type==0)
								rm->x[face/2] += shift_by[face];
							else if (type == 3)
								rm->xpred[face/2] += shift_by[face];
							else if (type==1){
								for(int f=0;f<6;f++) rm->fx[f][face/2] += shift_by[face];
							}
						}
						*/
						for (int ii = 0; ii < extraps.n[i0 + i + j*sm4 + k*smn4]; ii++) {
							rm = extraps.f[i0+i+j*sm4+k*smn4] + ii;
							//TODO hack for now, just write separate function to do all layers
							if (lid == -1 || rm->lid() == lid) {
								if(type==0) {
									rm->x[face/2] += shift_by[face];
									rm->xpred[face/2] += shift_by[face];
								} else if (type==1){
									for(int f=0;f<6;f++) rm->fx[f][face/2] += shift_by[face];
								} else  p_fatal_error("fluid_3d::shift_extraps(): Undefined type.",  PGMG_INTERNAL_ERROR);
							}
						}
					}
				}
	}
}

template<int flags>
int fluid_3d::copy_to_buffer(int tag,double *(&b),int lint,comm_info *ct) {

	// get pointer to start of region and buffer
	field *ur0 = u_mem + ct[tag].s0;
	ref_map *rmr0 = rm_mem + ct[tag].s0;

	// march through region, filling buffer. keep
	// same ordering as regular array (i.e. i innermost)
	int n_packed = 0;
	for (int k=0; k<ct[tag].o; k++) {
		for (int j=0; j<ct[tag].n; j++) {
			for (int i=0; i<ct[tag].m; i++){
				int eid = i + sm4*j + smn4*k;
				if (flags&FFLAG) ur0[eid].pack<flags&FFLAG>(b);
				if (flags&SFLAG) rmr0[eid].pack<flags&SFLAG>(b);
                if (flags&DVELFLAG) ur0[eid].pack<flags&DVELFLAG>(b);
				if (flags&EFLAG)
					n_packed += extraps.pack<flags&EFLAG>(eid+ct[tag].s0,
						b,lint);
			}
		}
	}

	if (flags&EFLAG) return n_packed;
	else return 0;
}

template<int flags>
void fluid_3d::copy_from_buffer(int tag,double *(&b),comm_info *ct) {

	// get pointer to start of region and buffer
	field *ur0 = u_mem + ct[tag].r0;
	ref_map *rmr0 = rm_mem + ct[tag].r0;

	// march through region, filling buffer. keep
	// same ordering as regular array (i.e. i innermost)
	for (int k=0; k<ct[tag].o; k++) {
		for (int j=0; j<ct[tag].n; j++) {
			for (int i=0; i<ct[tag].m; i++){
				int eid = i + sm4*j + smn4*k;
				if (flags&FFLAG) ur0[eid].unpack<flags&FFLAG>(b);
				if (flags&SFLAG) rmr0[eid].unpack<flags&SFLAG>(b);
                if (flags&DVELFLAG) ur0[eid].unpack<flags&DVELFLAG>(b);
				if (flags&EFLAG) extraps.unpack<flags&EFLAG>(eid+ct[tag].r0,
						b);
			}
		}
	}
}

/** Copy a single layer of extrapolated values in the extrap_mmap
 * to buffer, to be communicated.
 */
int fluid_3d::copy_layer_to_buffer(int tag,unsigned int c,double *(&b),int type) {

	comm_info (&ct)[27] = comm_table;

	int ind,i0 = ct[tag].s0;
	int en,n = 0;

	// keep same ordering as regular array (i.e. i innermost)
	for (int k=0; k<ct[tag].o; k++)
		for (int j=0; j<ct[tag].n; j++)
			for (int i=0; i<ct[tag].m; i++){

				ind = i0 + i + j*sm4 + k*smn4;
				en = extraps.pack_extrap(ind,c,b,type);
				n += 4*en + 1;
			}

	return n;
}

/** Copy all extrapolated values in the extrap_mmap
 * to buffer, to be communicated.
 */
int fluid_3d::copy_extraps_to_buffer(int tag,double *(&b),int type) {

	comm_info (&ct)[27] = comm_table;

	int ind,i0 = ct[tag].s0;
	int en,n = 0;

	// keep same ordering as regular array (i.e. i innermost)
	for (int k=0; k<ct[tag].o; k++)
		for (int j=0; j<ct[tag].n; j++)
			for (int i=0; i<ct[tag].m; i++){

				ind = i0 + i + j*sm4 + k*smn4;
				en = extraps.pack_extrap(ind,b,type);
				n += ((type==2)?19:4)*en + 1;
			}

	return n;
}

/** Copy all extrapolated values from the buffer,
 * to be stored in the extrap_mmap object.
 */
void fluid_3d::copy_extraps_from_buffer(int tag,double *(&b),int type) {

	comm_info (&ct)[27] = comm_table;

	int ind,i0 = ct[tag].r0;

	// march through region, filling buffer. keep
	// same ordering as regular array (i.e. i innermost)
	for (int k=0; k<ct[tag].o; k++)
		for (int j=0; j<ct[tag].n; j++)
			for (int i=0; i<ct[tag].m; i++){

				ind = i0 + i + j*sm4 + k*smn4;
				extraps.unpack_extrap(ind,b,type);
			}
}

// velocities
template void fluid_3d::communicate<1>(int lint);
template int fluid_3d::copy_to_buffer<1>(int tag,double *(&b),int lint,comm_info *ct);
template void fluid_3d::copy_from_buffer<1>(int tag,double *(&b),comm_info *ct);
// pressures and auxiliary pressure for MAC
template void fluid_3d::communicate<2>(int lint);
template int fluid_3d::copy_to_buffer<2>(int tag,double *(&b),int lint,comm_info *ct);
template void fluid_3d::copy_from_buffer<2>(int tag,double *(&b),comm_info *ct);
// object id, x-predicted
template void fluid_3d::communicate<8|32>(int lint);
template int fluid_3d::copy_to_buffer<8|32>(int tag,double *(&b),int lint,comm_info *ct);
template void fluid_3d::copy_from_buffer<8|32>(int tag,double *(&b),comm_info *ct);
// object id, x and x-predicted
template void fluid_3d::communicate<16|32>(int lint);
template int fluid_3d::copy_to_buffer<16|32>(int tag,double *(&b),int lint,comm_info *ct);
template void fluid_3d::copy_from_buffer<16|32>(int tag,double *(&b),comm_info *ct);
// face velocities, object id and reference maps
template void fluid_3d::communicate<4>(int lint);
template void fluid_3d::communicate<8>(int lint);
template void fluid_3d::communicate<64>(int lint);
template void fluid_3d::communicate<4|64>(int lint);
template void fluid_3d::communicate<4|8|64>(int lint);
template int fluid_3d::copy_to_buffer<4|8|64>(int tag,double *(&b),int lint,comm_info *ct);
template void fluid_3d::copy_from_buffer<4|8|64>(int tag,double *(&b),comm_info *ct);
// dvel
template void fluid_3d::communicate<2048>(int lint);
template int fluid_3d::copy_to_buffer<2048>(int tag,double *(&b),int lint,comm_info *ct);
template void fluid_3d::copy_from_buffer<2048>(int tag,double *(&b),comm_info *ct);
