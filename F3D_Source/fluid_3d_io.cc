#include "fluid_3d.hh"

void fluid_3d::write_slice(write_params wp,const char* filename) {
	const int dinfo[5]={2,2,1,0,0};
	int &dim=wp.dim;

	// get global dimensions, setting slice coord to 1
	int global_dim[3] = {m, n, o};
	// if there are four items in the field, it's fluid field, hack for now
	if(wp.corner_field()) {
		global_dim[0] = fem_grid.m;
		global_dim[1] = fem_grid.n;
		global_dim[2] = fem_grid.o;
	}
	global_dim[dim] = 1;

	// find out the other two dimension to iterate through
	int iter_dim1 = dinfo[dim+2], iter_dim2=dinfo[dim];

	int gslice_len = global_dim[0] * global_dim[1] * global_dim[2];
	// declare array, instantiate if master, and pull out slice

	double *u_global = NULL;
	if (rank == 0) {
		u_global = new double [gslice_len];
		for (int i=0; i<gslice_len; i++) *(u_global+i) = 0;
	}
	slice(wp, u_global);

	// if we're not the master, we're done. wait up and then quit.
	if (rank != 0) {
		MPI_Barrier(grid->cart);
		return;
	}

	// if we're the master, open file...
	FILE *fh = p_safe_fopen(filename, "w");

	// march through writing file (starting with largest y)
	if(wp.format==0) save_matrix(fh,u_global,global_dim[iter_dim1],global_dim[iter_dim2]);
	else {
		double *coord1 = gx;
		double *coord2 = gy;
		if(dim==0) { coord1 = gy; coord2 = gz; }
		else if(dim==1) { coord2 = gz; }
		if(wp.format==1) save_text(fh,u_global,global_dim[iter_dim1],global_dim[iter_dim2],coord1,coord2);
		else save_gnuplot(fh,u_global,global_dim[iter_dim1],global_dim[iter_dim2],coord1,coord2);
	}

	// close file and deallocate array
	fclose(fh);
	delete[] u_global;

	// join barrier, then we're done
	MPI_Barrier(grid->cart);
}

void fluid_3d::save_matrix(FILE *fh,double *u_global,int d1,int d2) {
	for (int j=0;j<d2;j++) {
		for (int i = 0; i < d1; i++) fprintf(fh,"%17.14g ", u_global[i + j*d1]);
		fputc('\n', fh);
	}
}

void fluid_3d::save_text(FILE *fh,double *u_global,int d1,int d2,double *a1,double *a2) {
	for (int j=0;j<d2;j++) {
		for (int i = 0; i < d1; i++)
			fprintf(fh,"%17.14f %17.14f %17.14f\n", a1[i], a2[j], u_global[i + j*d1]);
		fputc('\n', fh);
	}
}

void fluid_3d::save_gnuplot(FILE *fh,double *u_global,int d1,int d2,double *a1,double *a2) {
	float *buf=new float[d1+1],*bp=buf+1;*buf=static_cast<float>(m);
	int i,j;

	// Output the size and x coordinates
	for(i=0;i<d1;i++) bp[i]=static_cast<float>(a1[i]);
	fwrite(buf,sizeof(float),d1+1,fh);

	// Output the y coordinate and field values
	double *up=u_global;
	for(j=0;j<d2;j++) {
		*buf=static_cast<float>(a2[j]);
		for(i=0;i<d1;i++) bp[i]=static_cast<float>(*(up++));
		fwrite(buf,sizeof(float),d1+1,fh);
	}
	delete [] buf;
}

/**
* Grab a slice of nodes at constant (global) x-, y-, or z-index
* position. Only on master node, will write into array at address
* g_val in the form u(xi,yj) -> u_global[xi + m*yj] if slicing z.
* \param[in] wp the parameters for this output.
* \param[in] g_val the pointer to write out to.
*/
void fluid_3d::slice(write_params wp, double *g_val) {
	int &dim = wp.dim, &point = wp.point;

	bool self_involved=false;
	double * send_recv_buf = cbuf.buf;
	// Put needed local and global into array of 3 elements
	// in the order of x, y, z
	const int global_dim [3] = {m,n,o};
	const int p_global_dim [3] = {fem_grid.m,fem_grid.n,fem_grid.o};
	const int local_dim [3] = {sm,sn,so};
	const int p_local_dim [3] = {fem_grid.sm,fem_grid.sn,fem_grid.so};
	const int num_proc [3] ={grid->mp,grid->np,grid->op};
	const int proc_id[3] = {grid->ip,grid->jp,grid->kp};
	const int lowest [3] = {ai,aj,ak};
	const int highest[3] = {bi,bj,bk};

	// find out the other two dimension to iterate through
	int tmp1 = (dim + 1) % 3, tmp2 = (dim + 2) % 3;

	// iter_dim1 is the smaller of the two dimensions
	// suppose we are doing z slice, iter_dim1=x, iter_dim2 = y
	// supoose we are doing y slice, iter_dim1=x, iter_dim2 = z
	const int iter_dim1 = (tmp1<tmp2)?tmp1:tmp2;
	const int iter_dim2 = (tmp1<tmp2)?tmp2:tmp1;
	// make array of sender, length, dimension in global array
	int nprocs1 = num_proc[iter_dim1];
	int nprocs2 = num_proc[iter_dim2];
	// needed for data assembly
	int global_len = global_dim[iter_dim1];
	int local_dim1 = local_dim[iter_dim1];
	int local_dim2 = local_dim[iter_dim2];

	if(wp.corner_field()) {
		global_len = p_global_dim[iter_dim1];
		local_dim1 = p_local_dim[iter_dim1];
		local_dim2 = p_local_dim[iter_dim2];
	}

	// array and field labels
	int tag_n_len [5];
	const int STAG = 0;
	const int SLEN = 1;
	const int LDIM1 = 2;
	const int LDIM2 = 3;
	const int GDIM1 = 4;
	const int SENDR = 5;

	// find the processor index containing slice along perpendicular dimension
	int slice_ind = -1;
	if(point>=lowest[dim] && point<highest[dim]) slice_ind=proc_id[dim];

    // The following 4 lines is so that if the master node doesn't belong to the slice
    // it will still know which index the slice is at
	double tmp_ind_local = (slice_ind>=0)?((double) slice_ind / double (nprocs1*nprocs2)):0;
	double tmp_ind_master;
	MPI_Reduce(&tmp_ind_local, &tmp_ind_master, 1, MPI_DOUBLE, MPI_SUM, 0, grid->cart);
	if(rank==0 && slice_ind==-1) slice_ind=(int) tmp_ind_master;

	// just random tag to send tag
	const int TAGTAG = 2;
	int slice_len, slice_tag;

	// if this is a processor that contains the slice in question...
    if(wp.o_type == 18) {
        fill_pp_rhs();
    }
	if (proc_id[dim] == slice_ind) {
		// load the info on to the writing array
		slice_len = copy_slice_to_buf(wp,point-lowest[dim]);
		// the follow algebraic gymnastic to figure out the tags
		// tag is the starting position each local slice
		// but the position is computed with respect to the global slice
		// for a z slice, my tag is ai+aj*m
		// for a y slice, my tag is ai+ak*m
		// for a x slice, my tag is aj+ak*n
		slice_tag = lowest[iter_dim1] + lowest[iter_dim2] * global_len;
		// First send the slice tags to master node, with the tag -1
		tag_n_len[STAG] = slice_tag;
		tag_n_len[SLEN] = slice_len;
		tag_n_len[LDIM1] = local_dim1;
		tag_n_len[LDIM2] = local_dim2;
		tag_n_len[GDIM1] = global_len;

		if (rank==0) {
			// since master node doesn't send stuff to itself
			// we manually increment the buffer pointer, ready to receive
			send_recv_buf += slice_len;
			self_involved=true;
		}
		else {
			MPI_Send(tag_n_len,5,MPI_INT,0,TAGTAG,grid->cart);
			// then send the slice data to master node with slice tags
			MPI_Isend(send_recv_buf, slice_len, MPI_DOUBLE, 0, slice_tag, grid->cart, reqs);
			MPI_Waitall(1,reqs,stats);
		}
	}
	// now that all's sent, we collect on the master
	if (rank == 0) {

		// define array for cartesian processor coords
		int coords[3] = {0,0,0};
		coords[dim] = slice_ind;

		int **send_list = new int* [nprocs1*nprocs2];
		// now for each processor
		for (int i = 0; i < nprocs1*nprocs2; i++) {
			// make array of size six because we also store the order
			// in which these slice messages were received, i.e. SENDR
			send_list[i] = new int[6];
			for (int j = 0; j<6; j++) send_list[i][j] = 0;
		}
		int ** sublist = send_list;
		if(self_involved){
			// we in fact did not send from master to master
			// manually load data
			for (int i=0;i<5;i++) send_list[0][i] = tag_n_len[i];
			send_list[0][5] = 0;
			sublist++;
		}

		// iterate through all the processors that intersect
		// the slice, and recv tag and length info
		int sender, comm_co=0;
		for(int i = 0; i < nprocs1; i++) {
			for(int j = 0; j < nprocs2; j++) {

				// record the processor coords and get rank
				coords[iter_dim1] = i;
				coords[iter_dim2] = j;
				MPI_Cart_rank(grid->cart,coords,&sender);
				if (sender==0) continue;

				// get pointer to list for storing info
				*(*sublist+SENDR) = sender;
				// fill in rest of list with info from sender
				MPI_Recv(*sublist,5,MPI_INT,sender,TAGTAG,grid->cart,MPI_STATUS_IGNORE);
				// printf("from sender %d, sublist : { ", sender); for (int k=0;k<5;k++) printf("%d ", *(*sublist+k)); printf(" }\n");
				// includes slice data tag and length of slice data
				int pos = *(*sublist+STAG);
				int len = *(*sublist+SLEN);
				// now pull slice data into buffer
				MPI_Irecv(send_recv_buf,len,MPI_DOUBLE,sender,pos,grid->cart,reqs+comm_co);
				send_recv_buf += len;
				sublist++;comm_co++;
			}
		}
		// wait to finish everything, then copy out of buffer
		MPI_Waitall(comm_co,reqs,stats);
		copy_slice_from_buf(send_list,nprocs1*nprocs2,g_val);

		// delete the new array we made
		for (int i = 0; i < nprocs1*nprocs2; i++) if (send_list[i]!=NULL) delete [] send_list[i];
		if(send_list!=NULL) delete[] send_list;
	}
}

/** Copy slice at const i at local index ind into buffer to send.
* \param[in] wp the parameters for this output.
* \param[in] point the local coordinate of the slice to print. */
int fluid_3d::copy_slice_to_buf(write_params wp,int point) {
	int &dim = wp.dim;

	//printf("I'm processor %d, (%d, %d, %d), slicing at %d, point %d\n", rank, ai, aj, ak, dim, point);
	if(point<0 || (dim==0 && point>=sm) || (dim==1&&point>=sn) || (dim==2 && point>=so)){
		printf("copy_slice_to_buf(): Error! Index out of bound. Nothing has been done.\n");
		return 0;
	}

	// set limits for for loop -
	// start (change dir. in question to fixed point)
	int sp[3] = {0,0,0};

	// end position (dir. in question to fix_point + 1)
	int ep[3] = {sm, sn, so};
	// if writing out pressure we need to use fem_grid dimensions
	if(wp.corner_field()) {
		ep[0] = fem_grid.sm;
		ep[1] = fem_grid.sn;
		ep[2] = fem_grid.so;
	}

    sp[dim] = point;
	ep[dim] = point+1;
	double *b = cbuf.buf;
	for (int k = sp[2]; k<ep[2]; k++) for (int j = sp[1]; j < ep[1]; j++)
		for (int i = sp[0]; i < ep[0]; i++, b++) {
		// not outputing levelset
		int eid = i + sm4*j + smn4*k;
		int eidlow = eid-smn4;

        int obj_id=-1;
		if(wp.o_type<3) {
			switch(wp.o_type) {
				case 0: *b=u0[eid].vel[0];break;
				case 1: *b=u0[eid].vel[1];break;
				case 2: *b=u0[eid].vel[2];break;
            }
        } else if(wp.o_type==3) {
                // wx d(vz)/dy - d(vy)/dz
                double tmp_divu = 0.25*dysp*(u0[eid].vel[2] - u0[eid-sm4].vel[2] + u0[eid-1].vel[2] - u0[eid-1-sm4].vel[2]
				           + u0[eidlow].vel[2] - u0[eidlow-sm4].vel[2] + u0[eidlow-1].vel[2] - u0[eidlow-1-sm4].vel[2])
				- 0.25*dzsp*(u0[eid].vel[1] - u0[eidlow].vel[1] + u0[eid-1].vel[1] - u0[eidlow-1].vel[1]
					   + u0[eid-sm4].vel[1] - u0[eidlow-sm4].vel[1] + u0[eid-1-sm4].vel[1] - u0[eidlow-1-sm4].vel[1]);
                *b = tmp_divu;
        } else if(wp.o_type==4) {
                // wy d(vx)/dz - d(vz)/dx
                double tmp_divu = 0.25*dzsp*(u0[eid].vel[0] - u0[eidlow].vel[0] + u0[eid-1].vel[0] - u0[eidlow-1].vel[0]
					   + u0[eid-sm4].vel[0] - u0[eidlow-sm4].vel[0] + u0[eid-1-sm4].vel[0] - u0[eidlow-1-sm4].vel[0])
				- 0.25*dxsp*(u0[eid].vel[2] - u0[eid-1].vel[2] + u0[eid-sm4].vel[2] - u0[eid-sm4-1].vel[2]
					   + u0[eidlow].vel[2] - u0[eidlow-1].vel[2] + u0[eidlow-sm4].vel[2] - u0[eidlow-sm4-1].vel[2]);
                *b = tmp_divu;
        } else if(wp.o_type==5) {
                // wz d(vy)/dx - d(vx)/dy
                double tmp_divu = 0.25*dxsp*(u0[eid].vel[1] - u0[eid-1].vel[1] + u0[eid-sm4].vel[1] - u0[eid-sm4-1].vel[1]
					   + u0[eidlow].vel[1] - u0[eidlow-1].vel[1] + u0[eidlow-sm4].vel[1] - u0[eidlow-sm4-1].vel[1])
				- 0.25*dysp*(u0[eid].vel[0] - u0[eid-sm4].vel[0] + u0[eid-1].vel[0] - u0[eid-1-sm4].vel[0]
				           + u0[eidlow].vel[0] - u0[eidlow-sm4].vel[0] + u0[eidlow-1].vel[0] - u0[eidlow-1-sm4].vel[0]);

                *b = tmp_divu;
        } else if(wp.o_type==6) {
                double tmp_magu = 0;
                for(int ii=0;ii<3;ii++) {
                        tmp_magu += u0[eid].vel[ii] * u0[eid].vel[ii];
                }
                *b = sqrt(tmp_magu);
        } else if(wp.o_type==7){
                *b = u0[eid].p;
		} else if (wp.o_type == 18){
            *b = mg_pres->reg[0]->r0[i + fem_grid.sm * j + fem_grid.sm * fem_grid.sn * k];
        } else if (wp.o_type == 19) {
            double tmp_divu = dxsp * (u0[eid].fvel[1][0] - u0[eid].fvel[0][0]) + dysp * (u0[eid].fvel[3][1] - u0[eid].fvel[2][1]) + dzsp * (u0[eid].fvel[5][2] - u0[eid].fvel[4][2]);
            *b = tmp_divu;
        } else if (wp.o_type == 20) {
#if defined(VAR_DEN)
            *b = lrho[eid+G0];
#else
            *b = 0.;
#endif
        } else if (wp.o_type == 21) {
#if defined (DEBUG)
            *b = lJ[eid+G0];
#else
            *b = 100.;
#endif
        } else {
			if(wp.obj_id==-1) {
				obj_id=min_phi_id(eid);
				if(obj_id==-1) {
					if(wp.o_type ==12){
						*b=0.;
					} else if(wp.o_type<12 || wp.o_type > 12){
						*b=100.;
					}
					continue;
				}
			} else obj_id=wp.obj_id;

			if(rm0[eid].oid() == obj_id) {
				switch(wp.o_type) {
					case 8: *b = rm0[eid].phi(mgmt);break;
					case 9: *b = rm0[eid].x[0];    break;
					case 10: *b = rm0[eid].x[1];   break;
					case 11: *b = rm0[eid].x[2];   break;
					case 12: *b = static_cast<double> (rm0[eid].oid()); break;
					case 13: *b = static_cast<double> (rm0[eid].lid()); break;
					case 14: *b = rm0[eid].phi(mgmt,1);      break;
					case 15: *b = rm0[eid].xpred[0];        break;
					case 16: *b = rm0[eid].xpred[1];        break;
					case 17: *b = rm0[eid].xpred[2];        break;
				}
			} else {
				bool found=false;
				/*
				if(eid==1017723 || eid==1017725 || eid==1017926){
					printf("IO: Extrapolate maps at %d : %d\n", eid, extraps.n0[eid]);
				}
				*/
				for(int n=0;n<extraps.n0[eid];n++){
					ref_map &p = extraps.f0[eid][n];
					/*
						if(eid==1017723 || eid==1017725 || eid==1017926){
								printf("IO: map %d obj id %d layer id %d\n", n, p.oid(), p.lid());
						}
					*/
					if(p.oid() == obj_id) {
						switch(wp.o_type) {
							case 8: *b = p.phi(mgmt);break;
							case 9: *b = p.x[0];break;
							case 10: *b = p.x[1];break;
							case 11: *b = p.x[2];break;
                            // we are not printout object id for extrapolated values
                            // since object id is only used for checkpoint purposes
                            // and only primary grid refmap is needed for that
                            case 12: *b = 0.; break;
                            case 13: *b = static_cast<double> (p.lid()); break;
                            case 14: *b = p.phi(mgmt,1);break;
                            case 15: *b = p.xpred[0];  break;
                            case 16: *b = p.xpred[1];  break;
                            case 17: *b = p.xpred[2];  break;
						}
						found=true;break;
					}
				}
				if(!found) {
                    if(wp.o_type==12) *b=0.;
                    else *b=100.;
                }
			}
		}

	}

	// calculate and return length of slice copied to buffer
	int slice_len = 1;
	for (int i = 0; i < 3; i++) slice_len *= ep[i] - sp[i];
	return slice_len;
}

/** Computes the object ID for the minimum phi value at a gridpoint. If no phi
 * value is available, it returns -1.
 * \param[in] eid the gridpoint to consider.
 * return The object ID corresponding to the minimum phi value. */
int fluid_3d::min_phi_id(int eid) {
	int lo=rm0[eid].oid(),obj_id=-1;
	double mphi=1e30,tphi;
	if(lo>0) {obj_id=lo;mphi=rm0[eid].phi(mgmt);}
	for(int n=0;n<extraps.n0[eid];n++) {
		ref_map &p=extraps.f0[eid][n];
		tphi=p.phi(mgmt);
		if(tphi<mphi) {
			mphi=tphi;
			obj_id=p.oid();
		}
	}
	return obj_id;
}

/** Computes the minimum phi value a gridpoint. If no phi
 * value is available, it returns 100.
 * \param[in] eid the gridpoint to consider.
 * return The minimum phi value. */
double fluid_3d::min_phi(int eid) {
	double mphi=100,tphi;
	//printf("min_phi: eid %d\n", eid);
	if(rm0[eid].oid()>0) {
        //printf("Output: %d, x here %g\n", eid, rm0[eid].x[0]);
        mphi=rm0[eid].phi(mgmt);
    }
	for(int n=0;n<extraps.n0[eid];n++) {
		ref_map &p=extraps.f0[eid][n];
		tphi=p.phi(mgmt);
		if(tphi<mphi) mphi=tphi;
	}
	return mphi;
}

/** Computes the phi value a gridpoint corresponding to a specific object. If no
 * phi value is available, it returns 100.
 * \param[in] eid the gridpoint to consider.
 * \param[in] cnum the object to consider.
 * return The minimum phi value. */
double fluid_3d::phi(int eid,int cnum) {
	if(rm0[eid].oid()==cnum) return rm0[eid].phi(mgmt);
	for(int n=0;n<extraps.n0[eid];n++) {
		ref_map &p=extraps.f0[eid][n];
		if(p.oid()==cnum) return p.phi(mgmt);
	}
	return 100;
}

/* ##################### DATA OUTPUT #####################*/
/**
 * copy slice from buffer to the global array.
 * \param[in] send_list is an double pointer, indexed by the order
 *     messaged arrived. For each message, send_list[n] includes
 *     starting position [0], ldim1[2], ldim2[3], and gdim1[4].
 *     These are the starting position of the message in global
 *     array, dimensions of the local slice, low and high
 *     (respectively), and global low dimension.
 * \param[in] g_val pointer to write out to.
 */
void fluid_3d::copy_slice_from_buf(int **send_list, int co, double *g_val) {

	// fields of parameter array
	const int STAG = 0;
	const int LDIM1 = 2;
	const int LDIM2 = 3;
	const int GDIM1 = 4;

	// info for for loop
	int slice_tag,dim1,dim2,gdim1;

	// go to buffer
	double *b = cbuf.buf;
	for (int n = 0; n < co; n++) {
		// read in fields from array
		slice_tag = send_list[n][STAG];
		dim1 = send_list[n][LDIM1];
		dim2 = send_list[n][LDIM2];
		gdim1 = send_list[n][GDIM1];

		// loop through, copying from buffer
		for (int j = 0; j < dim2; j++)
			for(int i = 0; i < dim1; i++, b++) {
				*(g_val+slice_tag+i+j*gdim1) = *b;
		}
	}
}

void fluid_3d::set_dump(int dim,int pnt) {
	dd = dim;
	ds = pnt;
}

void fluid_3d::display_stats() {
	out = new char[1024];
	if (rank == 0) sprintf(out,"multigrid statistics\n");
}

void fluid_3d::dump(int num,char *dir) {
/*		char fid[4] = {'u','v','w','p'};
		char *file = new char[256];
		for (int f = 0; f < 4; f++) {
			sprintf(file,"%s/%c%d",dir,fid[f],num);
			write_params<field> params(dd, ds, 4, f, u0, file);
			write_slice<field>(params);
		}
		if (trace) {
			sprintf(file,"%s/t%d",dir,num);
			tr->print(file);
		}
		delete[] file;*/
}

/** Sets up the array that gives the dimensions of the other processors, needed
 * for output and also for grid transfers. */
void fluid_3d::setup_output_dimensions() {

	if(rank==0) gather_sizes();
	else {

		// On other ranks, check that there will be enough space to
		// send a local slice of output. For processors that are
		// orthogonally aligned with the (0,0,0) processor, send
		// information about the grid dimension.
		if(grid->kp==0) {
			if(grid->jp==0) MPI_Send(&sm,1,MPI_INT,0,msg_trans_dims,grid->cart);
			if(grid->ip==0) MPI_Send(&sn,1,MPI_INT,0,msg_trans_dims,grid->cart);
		} else if(grid->ip==0&&grid->jp==0) MPI_Send(&so,1,MPI_INT,0,msg_trans_dims,grid->cart);
	}
}

/** A routine run by the master processor to gather information about the
 * dimensions of the other regions. The routine also calculates the maximum
 * dimensions in each direction, which is needed for allocating memory for the
 * contour plotter. */
void fluid_3d::gather_sizes() {
	int q[4],&i=*q,&j=q[1],&k=q[2],l;j=k=0;
	int &mp=grid->mp,&np=grid->np,&op=grid->op;
	osm=new int[mp+np+op+3];osn=osm+mp;oso=osn+np;max_sizes=oso+op;
	int &msm=*max_sizes,&msn=max_sizes[1],&mso=max_sizes[2];

	// Receive dimensions in the x direction
	msm=*osm=sm;
	for(i=1;i<mp;i++) {
		MPI_Cart_rank(grid->cart,q,q+3);
		MPI_Recv(osm+i,1,MPI_INT,q[3],msg_trans_dims,grid->cart,stats);
		if(osm[i]>msm) msm=osm[i];
	}

	// Receive dimensions in the y direction
	msn=*osn=sn;i=0;
	for(j=1;j<np;j++) {
		MPI_Cart_rank(grid->cart,q,q+3);
		MPI_Recv(osn+j,1,MPI_INT,q[3],msg_trans_dims,grid->cart,stats);
		if(osn[j]>msn) msn=osn[j];
	}

	// Receive dimensions in the z direction
	mso=*oso=so;j=0;
	for(k=1;k<op;k++) {
		MPI_Cart_rank(grid->cart,q,q+3);
		MPI_Recv(oso+k,1,MPI_INT,q[3],msg_trans_dims,grid->cart,stats);
		if(oso[k]>mso) mso=oso[k];
	}

	// Set up the cumulative index tables
	cu_m=new int[mp+np+op+3];cu_n=cu_m+mp+1;cu_o=cu_n+np+1;
	*cu_m=0;for(l=0;l<mp;l++) cu_m[l+1]=cu_m[l]+osm[l];
	*cu_n=0;for(l=0;l<np;l++) cu_n[l+1]=cu_n[l]+osn[l];
	*cu_o=0;for(l=0;l<op;l++) cu_o[l+1]=cu_o[l]+oso[l];

}

/** Outputs the object boundaries in a binary version of the POV-Ray mesh2
 * format.
 * \param[in] filename the output file to write. */
void fluid_3d::output_contours(const char *filename) {

	// Initialize to suppress warnings
	FILE *outf=NULL;
	std::vector<double> pv;std::vector<int> q,q2;
	MPI_Request sreq[4];

	// Assemble the output filename and open the output file
	if(rank==0) {
		int odata[1];
		*odata=mgmt->n_obj;
		outf=p_safe_fopen(filename,"wb");
		if(fwrite(odata,sizeof(int),1,outf)!=1) {
			fputs("File write error",stderr);
			MPI_Abort(world,1);
		}
	}

	// Output each contour in the list
	for(int i=0; i<mgmt->n_obj; i++) {
		send_contour_data(pv,q,q2,sreq,i+1);
		if(rank==0) output_contour(outf,i);
		MPI_Waitall(4,sreq,MPI_STATUS_IGNORE);
	}

	// Close the output file
	if(rank==0) fclose(outf);
}

/** Outputs a single contour to an open file.
 * \param[in] outf the file handle to write to.
 * \param[in] cnum the number of the contour. */
void fluid_3d::output_contour(FILE *outf,int cnum) {
	int &mp=grid->mp,&np=grid->np,&op=grid->op;
	int l,aip=mp,ajp=np,akp=op,bip=0,bjp=0,bkp=0,cr[4],&ip=*cr,&jp=cr[1],&kp=cr[2];
	int *sz=new int[6*grid->procs], *szp=sz, *szr, *cor=sz + 3*grid->procs, *corp=cor;
	MPI_Request *reqp=reqs;

	// Receive how much data each processors is going to send
	for(kp=0;kp<op;kp++) for(jp=0;jp<np;jp++) for(ip=0;ip<mp;ip++,szp+=3,reqp++) {
		MPI_Cart_rank(grid->cart,cr,cr+3);
		MPI_Irecv(szp,3,MPI_INT,cr[3],fmsg_contour1,grid->cart,reqp);
	}
	MPI_Waitall(grid->procs,reqs,stats);

	// Compute the bounding box of processors involved in this cell, which
	// is used to allocate temporary memory
	int mvert=0,nvert=0,mq2=0,ntri=0,ijk=0;
	for(szp=szr=sz;szp<sz+3*grid->procs;szp+=3,ijk++) {

		// Count global measures of the contour
		if(*szp>mvert) mvert=*szp;
		nvert+=*szp;
		if(szp[1]>mq2) mq2=szp[1];
		ntri+=szp[2];

		// Find bounds of this contour and store a list of processors
		// that have vertices to send
		if(szp[1]>0) {
			*(corp++)=l=ijk%mp;
			if(l<aip) aip=l;
			if(l>bip) bip=l;
			*(corp++)=l=(ijk/mp)%np;
			if(l<ajp) ajp=l;
			if(l>bjp) bjp=l;
			*(corp++)=l=ijk/(mp*np);
			if(l<akp) akp=l;
			if(l>bkp) bkp=l;
			*(szr++)=*szp;
			*(szr++)=szp[1];
			*(szr++)=szp[2];
		}
	}

	int odata[2];*odata=nvert;odata[1]=ntri;
	//printf("# Object %d: nvert=%d ntri=%d\n",cnum,*odata,odata[1]);
	fwrite(odata,sizeof(int),2,outf);

	// Skip completely if there are no vertices
	if(nvert==0) {
		delete [] sz;
		return;
	}

	// Allocate memory for the normal vectors and psoitions
	float *nl=new float[3*nvert],*nlp=nl,*pts=new float[3*mvert],*ptsp;

	// Allocate memory
	double *pinfo=new double[4*mvert];
	int *vt=new int[mvert],*vbx=new int[mq2],*tri=new int[3*ntri],
	    *trip=tri,

	// Compute dimensions of mapping regions, adding one if the last
	// processor in a grid is being used
	    lrm=*max_sizes+1,rm=cu_m[bip+1]-cu_m[aip]+1,
	    lrn=max_sizes[1]+1,rn=cu_n[bjp+1]-cu_n[ajp]+1,
	    lro=max_sizes[2]+1,

	// Allocate memory for mapping regions
	    mrs=lrn*lro+rm*lro+rm*rn,*mrx=new int[4*mrs],
	    *mry=mrx+4*lrn*lro,*mrz=mry+4*rm*lro,*ml1,*ml2,*ml3,*mu1,*mu2,*mu3,
	    *fm=new int[3*lrm*lrn*lro],

	// Begin reading the details of the mesh, storing some of it, and
	// printing the vertex vectors
	    vi,vj,vk,vijk,vn=0,type,*vbxp,*vbxe;
	for(int oo=0;oo<3*lrm*lrn*lro;oo++) fm[oo]=0;
	char *pp;
	double coords[3],&x=*coords,&y=coords[1],&z=coords[2];
	bool posi,posj,posk;
	for(int *szp=sz,*corp=cor;szp<szr;szp+=3,corp+=3) {
		int &i=*corp,&j=corp[1],&k=corp[2],
		    lsm4=osm[i]+4,lsn4=osn[j]+4,lsmn4=lsm4*lsn4;
		MPI_Cart_rank(grid->cart,corp,cr+3);

		ptsp=pts;

		// Set up the boundary layers
		if(i&1) {ml1=mrx;mu1=mrx+lrn*lro;} else {ml1=mrx+lrn*lro;mu1=mrx;}
		if(j&1) {ml2=mry;mu2=mry+rm*lro;} else {ml2=mry+rm*lro;mu2=mry;}
		if(k&1) {ml3=mrz;mu3=mrz+rm*rn;} else {ml3=mrz+rm*rn;mu3=mrz;}

		// Deal with the mesh points if there are any
		if(*szp>0) {

			// Determine if this box borders any of the upper
			// regions of the grid
			posi=i!=mp-1;
			posj=j!=np-1;
			posk=k!=op-1;

			// Receive index, position, and normal information
			MPI_Recv(vt,*szp,MPI_INT,cr[3],fmsg_contour2,grid->cart,MPI_STATUS_IGNORE);
			MPI_Recv(pinfo,4*(*szp),MPI_DOUBLE,cr[3],fmsg_contour3,grid->cart,MPI_STATUS_IGNORE);

			// Loop over the available points
			for(l=0;l<*szp;l++) {

				// Decode the position
				type=vt[l]&3;
				vijk=vt[l]>>2;
				vi=vijk%lsm4;
				vj=(vijk/lsm4)%lsn4;
				vk=vijk/lsmn4;

				// Calculate and store the position
				x=gx[cu_m[i]+vi];
				y=gy[cu_n[j]+vj];
				z=gz[cu_o[k]+vk];
				switch(type) {
					case 0: x+=pinfo[l<<2];break;
					case 1: y+=pinfo[l<<2];break;
					case 2: z+=pinfo[l<<2];break;
				}
				*(ptsp++)=x;*(ptsp++)=y;*(ptsp++)=z;

				// Store the normal vector
				*(nlp++)=pinfo[(l<<2)+1];
				*(nlp++)=pinfo[(l<<2)+2];
				*(nlp++)=pinfo[(l<<2)+3];

				// Store the index for later matching
				if(posk&&vk==oso[k]) {
					if(type==2) p_die();
					mu3[((cu_m[i]-cu_m[aip]+vi+rm*(cu_n[j]-cu_n[ajp]+vj))<<1)+type]=vn;
				}
				if(posj&&vj==osn[j]) {
					if(type==1) p_die();
					mu2[((cu_m[i]-cu_m[aip]+vi+rm*vk)<<1)+(type>>1)]=vn;
				}
				if(posi&&vi==osm[i]) {
					if(type==0) p_die();
					mu1[((vj+lrn*vk)<<1)+type-1]=vn;
				}
				fm[(vi+lrm*(vj+lrn*vk))*3+type]=vn++;
			}

			// Output the positions
			fwrite(pts,sizeof(float),3*(*szp),outf);
		}

		// Reindex the face indices
		if(sz[1]>0) {

			// Determine if this box borders any of the lower
			// regions of the grid
			posi=i==0;
			posj=j==0;
			posk=k==0;

			// Receive the block information
			MPI_Recv(vbx,szp[1],MPI_INT,cr[3],fmsg_contour4,grid->cart,MPI_STATUS_IGNORE);

			vbxp=vbx;vbxe=vbx+szp[1];
			int d1=3*lrm,d2=d1*lrn,key;
			while(vbxp<vbxe) {

				// Get coordinates and then convert into format
				// for looking up edge indices
				vijk=*(vbxp++);
				vi=vijk%lsm4;
				vj=(vijk/lsm4)%lsn4;
				vk=vijk/lsmn4;
				vijk=3*(vi+lrm*(vj+lrn*vk));

				// Look up the indices of the vertices involved
				key=*(vbxp++);
				for(pp=ptri_poly[key];pp<ptri_poly[key+1];pp++,trip++) switch(*pp) {
					case 0: *trip=posk||vk!=0?
						(posj||vj!=0?fm[vijk]:ml2[(cu_m[i]-cu_m[aip]+vi+rm*vk)<<1]):
						ml3[(cu_m[i]-cu_m[aip]+vi+rm*(cu_n[j]-cu_n[ajp]+vj))<<1];break;
					case 1: *trip=posk||vk!=0?fm[vijk+4]:
						ml3[((cu_m[i]-cu_m[aip]+vi+1+rm*(cu_n[j]-cu_n[ajp]+vj))<<1)+1];break;
					case 2: *trip=posk||vk!=0?fm[vijk+d1]:
						ml3[(cu_m[i]-cu_m[aip]+vi+rm*(cu_n[j]-cu_n[ajp]+vj+1))<<1];break;
					case 3: *trip=posk||vk!=0?
						(posi||vi!=0?fm[vijk+1]:ml1[((vj+lrn*vk)<<1)]):
						ml3[((cu_m[i]-cu_m[aip]+vi+rm*(cu_n[j]-cu_n[ajp]+vj))<<1)+1];break;
					case 4: *trip=posj||vj!=0?fm[vijk+d2]:
						ml2[(cu_m[i]-cu_m[aip]+vi+rm*(vk+1))<<1];break;
					case 5: *trip=fm[vijk+d2+4];break;
					case 6: *trip=fm[vijk+d2+d1];break;
					case 7: *trip=posi||vi!=0?fm[vijk+d2+1]:
						ml1[((vj+lrn*(vk+1))<<1)];break;
					case 8: *trip=posj||vj!=0?
						(posi||vi!=0?fm[vijk+2]:ml1[((vj+lrn*vk)<<1)+1]):
						ml2[((cu_m[i]-cu_m[aip]+vi+rm*vk)<<1)+1];break;
					case 9: *trip=posj||vj!=0?fm[vijk+5]:
						ml2[((cu_m[i]-cu_m[aip]+vi+1+rm*vk)<<1)+1];break;
					case 10: *trip=fm[vijk+d1+5];break;
					case 11: *trip=posi||vi!=0?fm[vijk+d1+2]:
						 ml1[((vj+1+lrn*vk)<<1)+1];break;
				}
			}
		}
	}

	// Print the normal vectors
	fwrite(nl,sizeof(float),3*nvert,outf);

	// Print out the face indices
	fwrite(tri,sizeof(int),3*ntri,outf);

	// Deallocate the dynamically allocated memory
	delete [] fm;
	delete [] mrx;
	delete [] tri;
	delete [] vbx;
	delete [] vt;
	delete [] pinfo;
	delete [] pts;
	delete [] nl;
	delete [] sz;
}

void fluid_3d::send_contour_data(std::vector<double> &pv,std::vector<int> &q,std::vector<int> &q2,MPI_Request *sreq,int cnum) {
	pv.clear();q.clear();q2.clear();
	unsigned int ttri=0;
	//geometry &g = *grid;
	bool posi=grid->ip==0,posj=grid->jp==0,posk=grid->kp==0;
	//int ui=(grid->ip==grid->mp-1)?sm-1:sm,uj=(grid->jp==grid->np-1)?sn-1:sn,uk=(grid->kp==grid->op-1)?so-1:so;

	for(int k=0;k<=so;k++) for(int j=0;j<=sn;j++) for(int i=0;i<=sm;i++) {
		int ijk=i+sm4*j+smn4*k;
		if((k>0||posk)&&(j>0||posj)&&i<sm) edge_detect(ijk,0,pv,q,cnum);
		if((k>0||posk)&&j<sn&&(i>0||posi)) edge_detect(ijk,1,pv,q,cnum);
		if(k<so&&(j>0||posj)&&(i>0||posi)) edge_detect(ijk,2,pv,q,cnum);
		if(i<sm&&j<sn&&k<so) box_detect(ijk,q2,ttri,cnum);
		//if(i<ui&&j<uj&&k<uk) box_detect(ijk,q2,ttri,cval);
	}

	// Pack the sizes and send them to the master processor
	*sizes=q.size();
	sizes[1]=q2.size();
	sizes[2]=ttri;
	MPI_Isend(sizes,3,MPI_INT,0,fmsg_contour1,grid->cart,sreq);

	// Send information about vertices to the master processor
	if(sizes[0]>0) {
		MPI_Isend(&q[0],q.size(),MPI_INT,0,fmsg_contour2,grid->cart,sreq+1);
		MPI_Isend(&pv[0],4*q.size(),MPI_DOUBLE,0,fmsg_contour3,grid->cart,sreq+2);
	} else sreq[1]=sreq[2]=MPI_REQUEST_NULL;

	// Send information about blocks to the master processor
	if(sizes[1]>0) MPI_Isend(&q2[0],q2.size(),MPI_INT,0,fmsg_contour4,grid->cart,sreq+3);
	else sreq[3]=MPI_REQUEST_NULL;
}

/* Determines whether the zero level set intersects a given edge, and if so
 * adds it to the buffer to send to the master processor.
 * \param[in] ijk the point to consider.
 * \param[in] dir the direction of the edge (0=x, 1=y, 2=z). */
void fluid_3d::edge_detect(int ijk,int dir,std::vector<double> &pv,std::vector<int> &q,int cnum) {
	const double tol=std::numeric_limits<double>::epsilon();
	int ind=(ijk<<2)|dir;
	int ijk2=ijk+(dir==0?1:(dir==1?sm4:smn4));
	double v1=phi(ijk,cnum),v2=phi(ijk2,cnum);

	if(v1>0?v2<=0:v2>0) {
		double nx,ny,nz,nx2,ny2,nz2;
		double fa=v2-v1,fb=fabs(fa)<tol?0.5:-v1/fa;
		pv.push_back(fb*(dir==0?dx:(dir==1?dy:dz)));
		q.push_back(ind);
		grad_field(ijk,nx,ny,nz,cnum);
		grad_field(ijk2,nx2,ny2,nz2,cnum);
		nx=(1-fb)*nx+fb*nx2;
		ny=(1-fb)*ny+fb*ny2;
		nz=(1-fb)*nz+fb*nz2;
		normalize(nx,ny,nz);
		pv.push_back(nx);
		pv.push_back(ny);
		pv.push_back(nz);
	}
}

/** Computes the gradient of the field at a gridpoint.
 * \param[in] ijk the point to consider.
 * \param[in] (nx,ny,nz) the components of gradient. */
void fluid_3d::grad_field(int ijk,double &nx,double &ny,double &nz,int cnum) {
	nx=dxsp*(phi(ijk+1,cnum)-phi(ijk-1,cnum));
	ny=dysp*(phi(ijk+sm4,cnum)-phi(ijk-sm4,cnum));
	nz=dzsp*(phi(ijk+smn4,cnum)-phi(ijk-smn4,cnum));
}

/* Examines whether a given box is part of a boundary, and if so, adds the
 * index to the list to send to the master processor.
 * \param[in] ijk the index fo the box to consider.
 * \param[in,out] q2 a vector to add the index and key to.
 * \param[in,out] ttri the running total number of triangles. */
void fluid_3d::box_detect(int ijk,std::vector<int> &q2,unsigned int &ttri,int cnum) {

	// Test the eight corners of the box and assemble the key
	int key=phi(ijk,cnum)>0?1:0;
	if(phi(ijk+1,cnum)>0) key|=2;
	if(phi(ijk+sm4+1,cnum)>0) key|=4;
	if(phi(ijk+sm4,cnum)>0) key|=8;
	if(phi(ijk+smn4,cnum)>0) key|=16;
	if(phi(ijk+smn4+1,cnum)>0) key|=32;
	if(phi(ijk+smn4+sm4+1,cnum)>0) key|=64;
	if(phi(ijk+smn4+sm4,cnum)>0) key|=128;

	// If contour passes through this box, then store the key value and box
	// index
	if(key!=255&&key!=0) {
		q2.push_back(ijk);
		q2.push_back(key);
		ttri+=n_poly[key];
	}
}

/** Normalizes a vector.
 * \param[in,out] (nx,ny,nz) the components of the vector, which upon
 *			     completion will have length 1. */
void fluid_3d::normalize(double &nx,double &ny,double &nz) {
	const double tol=std::numeric_limits<double>::epsilon();
	double rsq=nx*nx+ny*ny+nz*nz;
	if(rsq>tol*tol) {
		rsq=1./sqrt(rsq);
		nx*=rsq;
		ny*=rsq;
		nz*=rsq;
	}
}

/** The number of triangles in the contour configuration of each block. */
const char fluid_3d::n_poly[256]={
	0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,2,
	1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,3,
	1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,3,
	2,3,3,2,3,4,4,3,3,4,4,3,4,5,5,2,
	1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,3,
	2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,4,
	2,3,3,4,3,4,2,3,3,4,4,5,4,5,3,2,
	3,4,4,3,4,5,3,2,4,5,5,4,5,2,4,1,
	1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,3,
	2,3,3,4,3,4,4,5,3,2,4,3,4,3,5,2,
	2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,4,
	3,4,4,3,4,5,5,4,4,3,5,2,5,4,2,1,
	2,3,3,4,3,4,4,5,3,4,4,5,2,3,3,2,
	3,4,4,5,4,5,5,2,4,3,5,4,3,2,4,1,
	3,4,4,5,4,5,3,4,4,5,5,2,3,4,2,1,
	2,3,3,2,3,4,2,1,3,2,4,1,2,1,1,0
};

/** The connectivity of the triangles in the contour configuration of each
 * block. */
const char fluid_3d::tri_poly[2460]={
	0,8,3,0,1,9,1,8,3,9,8,1,1,2,10,0,8,3,1,2,10,9,2,10,0,2,9,2,8,3,
	2,10,8,10,9,8,3,11,2,0,11,2,8,11,0,1,9,0,2,3,11,1,11,2,1,9,11,9,8,11,
	3,10,1,11,10,3,0,10,1,0,8,10,8,11,10,3,9,0,3,11,9,11,10,9,9,8,10,10,8,11,
	4,7,8,4,3,0,7,3,4,0,1,9,8,4,7,4,1,9,4,7,1,7,3,1,1,2,10,8,4,7,
	3,4,7,3,0,4,1,2,10,9,2,10,9,0,2,8,4,7,2,10,9,2,9,7,2,7,3,7,9,4,
	8,4,7,3,11,2,11,4,7,11,2,4,2,0,4,9,0,1,8,4,7,2,3,11,4,7,11,9,4,11,
	9,11,2,9,2,1,3,10,1,3,11,10,7,8,4,1,11,10,1,4,11,1,0,4,7,11,4,4,7,8,
	9,0,11,9,11,10,11,0,3,4,7,11,4,11,9,9,11,10,9,5,4,9,5,4,0,8,3,0,5,4,
	1,5,0,8,5,4,8,3,5,3,1,5,1,2,10,9,5,4,3,0,8,1,2,10,4,9,5,5,2,10,
	5,4,2,4,0,2,2,10,5,3,2,5,3,5,4,3,4,8,9,5,4,2,3,11,0,11,2,0,8,11,
	4,9,5,0,5,4,0,1,5,2,3,11,2,1,5,2,5,8,2,8,11,4,8,5,10,3,11,10,1,3,
	9,5,4,4,9,5,0,8,1,8,10,1,8,11,10,5,4,0,5,0,11,5,11,10,11,0,3,5,4,8,
	5,8,10,10,8,11,9,7,8,5,7,9,9,3,0,9,5,3,5,7,3,0,7,8,0,1,7,1,5,7,
	1,5,3,3,5,7,9,7,8,9,5,7,10,1,2,10,1,2,9,5,0,5,3,0,5,7,3,8,0,2,
	8,2,5,8,5,7,10,5,2,2,10,5,2,5,3,3,5,7,7,9,5,7,8,9,3,11,2,9,5,7,
	9,7,2,9,2,0,2,7,11,2,3,11,0,1,8,1,7,8,1,5,7,11,2,1,11,1,7,7,1,5,
	9,5,8,8,5,7,10,1,3,10,3,11,5,7,0,5,0,9,7,11,0,1,0,10,11,10,0,11,10,0,
	11,0,3,10,5,0,8,0,7,5,7,0,11,10,5,7,11,5,10,6,5,0,8,3,5,10,6,9,0,1,
	5,10,6,1,8,3,1,9,8,5,10,6,1,6,5,2,6,1,1,6,5,1,2,6,3,0,8,9,6,5,
	9,0,6,0,2,6,5,9,8,5,8,2,5,2,6,3,2,8,2,3,11,10,6,5,11,0,8,11,2,0,
	10,6,5,0,1,9,2,3,11,5,10,6,5,10,6,1,9,2,9,11,2,9,8,11,6,3,11,6,5,3,
	5,1,3,0,8,11,0,11,5,0,5,1,5,11,6,3,11,6,0,3,6,0,6,5,0,5,9,6,5,9,
	6,9,11,11,9,8,5,10,6,4,7,8,4,3,0,4,7,3,6,5,10,1,9,0,5,10,6,8,4,7,
	10,6,5,1,9,7,1,7,3,7,9,4,6,1,2,6,5,1,4,7,8,1,2,5,5,2,6,3,0,4,
	3,4,7,8,4,7,9,0,5,0,6,5,0,2,6,7,3,9,7,9,4,3,2,9,5,9,6,2,6,9,
	3,11,2,7,8,4,10,6,5,5,10,6,4,7,2,4,2,0,2,7,11,0,1,9,4,7,8,2,3,11,
	5,10,6,9,2,1,9,11,2,9,4,11,7,11,4,5,10,6,8,4,7,3,11,5,3,5,1,5,11,6,
	5,1,11,5,11,6,1,0,11,7,11,4,0,4,11,0,5,9,0,6,5,0,3,6,11,6,3,8,4,7,
	6,5,9,6,9,11,4,7,9,7,11,9,10,4,9,6,4,10,4,10,6,4,9,10,0,8,3,10,0,1,
	10,6,0,6,4,0,8,3,1,8,1,6,8,6,4,6,1,10,1,4,9,1,2,4,2,6,4,3,0,8,
	1,2,9,2,4,9,2,6,4,0,2,4,4,2,6,8,3,2,8,2,4,4,2,6,10,4,9,10,6,4,
	11,2,3,0,8,2,2,8,11,4,9,10,4,10,6,3,11,2,0,1,6,0,6,4,6,1,10,6,4,1,
	6,1,10,4,8,1,2,1,11,8,11,1,9,6,4,9,3,6,9,1,3,11,6,3,8,11,1,8,1,0,
	11,6,1,9,1,4,6,4,1,3,11,6,3,6,0,0,6,4,6,4,8,11,6,8,7,10,6,7,8,10,
	8,9,10,0,7,3,0,10,7,0,9,10,6,7,10,10,6,7,1,10,7,1,7,8,1,8,0,10,6,7,
	10,7,1,1,7,3,1,2,6,1,6,8,1,8,9,8,6,7,2,6,9,2,9,1,6,7,9,0,9,3,
	7,3,9,7,8,0,7,0,6,6,0,2,7,3,2,6,7,2,2,3,11,10,6,8,10,8,9,8,6,7,
	2,0,7,2,7,11,0,9,7,6,7,10,9,10,7,1,8,0,1,7,8,1,10,7,6,7,10,2,3,11,
	11,2,1,11,1,7,10,6,1,6,7,1,8,9,6,8,6,7,9,1,6,11,6,3,1,3,6,0,9,1,
	11,6,7,7,8,0,7,0,6,3,11,0,11,6,0,7,11,6,7,6,11,3,0,8,11,7,6,0,1,9,
	11,7,6,8,1,9,8,3,1,11,7,6,10,1,2,6,11,7,1,2,10,3,0,8,6,11,7,2,9,0,
	2,10,9,6,11,7,6,11,7,2,10,3,10,8,3,10,9,8,7,2,3,6,2,7,7,0,8,7,6,0,
	6,2,0,2,7,6,2,3,7,0,1,9,1,6,2,1,8,6,1,9,8,8,7,6,10,7,6,10,1,7,
	1,3,7,10,7,6,1,7,10,1,8,7,1,0,8,0,3,7,0,7,10,0,10,9,6,10,7,7,6,10,
	7,10,8,8,10,9,6,8,4,11,8,6,3,6,11,3,0,6,0,4,6,8,6,11,8,4,6,9,0,1,
	9,4,6,9,6,3,9,3,1,11,3,6,6,8,4,6,11,8,2,10,1,1,2,10,3,0,11,0,6,11,
	0,4,6,4,11,8,4,6,11,0,2,9,2,10,9,10,9,3,10,3,2,9,4,3,11,3,6,4,6,3,
	8,2,3,8,4,2,4,6,2,0,4,2,4,6,2,1,9,0,2,3,4,2,4,6,4,3,8,1,9,4,
	1,4,2,2,4,6,8,1,3,8,6,1,8,4,6,6,10,1,10,1,0,10,0,6,6,0,4,4,6,3,
	4,3,8,6,10,3,0,3,9,10,9,3,10,9,4,6,10,4,4,9,5,7,6,11,0,8,3,4,9,5,
	11,7,6,5,0,1,5,4,0,7,6,11,11,7,6,8,3,4,3,5,4,3,1,5,9,5,4,10,1,2,
	7,6,11,6,11,7,1,2,10,0,8,3,4,9,5,7,6,11,5,4,10,4,2,10,4,0,2,3,4,8,
	3,5,4,3,2,5,10,5,2,11,7,6,7,2,3,7,6,2,5,4,9,9,5,4,0,8,6,0,6,2,
	6,8,7,3,6,2,3,7,6,1,5,0,5,4,0,6,2,8,6,8,7,2,1,8,4,8,5,1,5,8,
	9,5,4,10,1,6,1,7,6,1,3,7,1,6,10,1,7,6,1,0,7,8,7,0,9,5,4,4,0,10,
	4,10,5,0,3,10,6,10,7,3,7,10,7,6,10,7,10,8,5,4,10,4,8,10,6,9,5,6,11,9,
	11,8,9,3,6,11,0,6,3,0,5,6,0,9,5,0,11,8,0,5,11,0,1,5,5,6,11,6,11,3,
	6,3,5,5,3,1,1,2,10,9,5,11,9,11,8,11,5,6,0,11,3,0,6,11,0,9,6,5,6,9,
	1,2,10,11,8,5,11,5,6,8,0,5,10,5,2,0,2,5,6,11,3,6,3,5,2,10,3,10,5,3,
	5,8,9,5,2,8,5,6,2,3,8,2,9,5,6,9,6,0,0,6,2,1,5,8,1,8,0,5,6,8,
	3,8,2,6,2,8,1,5,6,2,1,6,1,3,6,1,6,10,3,8,6,5,6,9,8,9,6,10,1,0,
	10,0,6,9,5,0,5,6,0,0,3,8,5,6,10,10,5,6,11,5,10,7,5,11,11,5,10,11,7,5,
	8,3,0,5,11,7,5,10,11,1,9,0,10,7,5,10,11,7,9,8,1,8,3,1,11,1,2,11,7,1,
	7,5,1,0,8,3,1,2,7,1,7,5,7,2,11,9,7,5,9,2,7,9,0,2,2,11,7,7,5,2,
	7,2,11,5,9,2,3,2,8,9,8,2,2,5,10,2,3,5,3,7,5,8,2,0,8,5,2,8,7,5,
	10,2,5,9,0,1,5,10,3,5,3,7,3,10,2,9,8,2,9,2,1,8,7,2,10,2,5,7,5,2,
	1,3,5,3,7,5,0,8,7,0,7,1,1,7,5,9,0,3,9,3,5,5,3,7,9,8,7,5,9,7,
	5,8,4,5,10,8,10,11,8,5,0,4,5,11,0,5,10,11,11,3,0,0,1,9,8,4,10,8,10,11,
	10,4,5,10,11,4,10,4,5,11,3,4,9,4,1,3,1,4,2,5,1,2,8,5,2,11,8,4,5,8,
	0,4,11,0,11,3,4,5,11,2,11,1,5,1,11,0,2,5,0,5,9,2,11,5,4,5,8,11,8,5,
	9,4,5,2,11,3,2,5,10,3,5,2,3,4,5,3,8,4,5,10,2,5,2,4,4,2,0,3,10,2,
	3,5,10,3,8,5,4,5,8,0,1,9,5,10,2,5,2,4,1,9,2,9,4,2,8,4,5,8,5,3,
	3,5,1,0,4,5,1,0,5,8,4,5,8,5,3,9,0,5,0,3,5,9,4,5,4,11,7,4,9,11,
	9,10,11,0,8,3,4,9,7,9,11,7,9,10,11,1,10,11,1,11,4,1,4,0,7,4,11,3,1,4,
	3,4,8,1,10,4,7,4,11,10,11,4,4,11,7,9,11,4,9,2,11,9,1,2,9,7,4,9,11,7,
	9,1,11,2,11,1,0,8,3,11,7,4,11,4,2,2,4,0,11,7,4,11,4,2,8,3,4,3,2,4,
	2,9,10,2,7,9,2,3,7,7,4,9,9,10,7,9,7,4,10,2,7,8,7,0,2,0,7,3,7,10,
	3,10,2,7,4,10,1,10,0,4,0,10,1,10,2,8,7,4,4,9,1,4,1,7,7,1,3,4,9,1,
	4,1,7,0,8,1,8,7,1,4,0,3,7,4,3,4,8,7,9,10,8,10,11,8,3,0,9,3,9,11,
	11,9,10,0,1,10,0,10,8,8,10,11,3,1,10,11,3,10,1,2,11,1,11,9,9,11,8,3,0,9,
	3,9,11,1,2,9,2,11,9,0,2,11,8,0,11,3,2,11,2,3,8,2,8,10,10,8,9,9,10,2,
	0,9,2,2,3,8,2,8,10,0,1,8,1,10,8,1,10,2,1,3,8,9,1,8,0,9,1,0,3,8
};
