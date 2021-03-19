#include "fluid_3d.hh"

void fluid_3d::write_chk_pt(const int snum, const int fnum, const char* chk_dirname){
    // Write tracer file
	// We need to output field flag 1, 2, 4, 128, 512, 1024, 2048
	int out_ind=1, out_f=0, obj_id=-1, data_fmt=0;
	const int numf = 8;
	int output_field[numf] = {0,1,2,7,9,10,11,12};
	int Z = 2;
	write_params wp(Z, out_ind, out_f, obj_id, data_fmt);

	// We first iterate through the fields we need to output
	// then iterate through all the slices

	// some common variables
	int gslice_len = (m+1) * (n+1);

	double *u_global = NULL;
	char * filename=NULL;
    char * compress_command=NULL;
	FILE *fh=NULL;
	if(rank==0){
		filename = new char[512];
        compress_command = new char [1024];
		u_global = new double[gslice_len];
		for(int i=0;i<gslice_len;i++) u_global[i]=0;

        sprintf(compress_command, "if [ -f %s/tr.%05d ]; then cp %s/tr.%05d %s; fi;", spars->dirname, fnum, spars->dirname, fnum, chk_dirname);
        int cp_err = system(compress_command);
        if(cp_err){
            printf("fluid_3d::write_chk_pt: failed to copy tracer file.\n");
        }
	}

	for(int j=0;j<numf;j++){
		// Change write_params output type
		out_f = output_field[j];
		wp.change_otype(out_f);
		if(rank==0) {
			sprintf(filename, "%s/%s", chk_dirname, wp.filename);
			fh = p_safe_fopen(filename, "wb");
		}
		for(int i=0;i<o;i++) {
			wp.change_point(i);
			slice(wp, u_global);
			if(rank!=0) MPI_Barrier(grid->cart);
			if(rank==0) {
				fwrite(u_global, sizeof(double), gslice_len,fh);
				MPI_Barrier(grid->cart);
			}
		}
		if(rank==0){
            fclose(fh);
            // Compress the file using xz, if file exists, delete and recompress freshly generated checkpoint files
            sprintf(compress_command, "if [ -f %s.xz ]; then rm %s.xz; fi;  xz %s;", filename, filename, filename);
            int compress_err = system(compress_command);
            if(compress_err) {
                printf("fluid_3d::write_chk_pt: compression command %s went wrong!\n", compress_command);
            }
        }
	}
	if(rank==0){
		delete [] u_global;
        delete [] compress_command;
		delete [] filename;
	}

	spars->set_current_chk_num(fnum);
	spars->set_current_chk_step(snum);
	spars->set_current_time(time);
	spars->write_params(chk_dirname);
}

int fluid_3d::initialize_from_chk_point(const char * chk_dirname){
    int err=0;
	if(rank==0) printf("fluid_3d:: initializing from files in %s/\n", chk_dirname);
	time = spars->get_current_time();
    nt = spars->get_current_step();
	int out_ind=1, out_f=0, obj_id=-1, data_fmt=0;
	// We need to output field flag 1, 2, 4, 128, 512, 1024, 2048
	// take log base2 0,1,2,7 are (u,v,w,p), 9,10,11 are refx,refy,refz
	const int numf = 8;
	int output_field[numf] = {0,1,2,7,9,10,11,12};
	int Z = 2;
	write_params wp(Z, out_ind, out_f, obj_id, data_fmt);

	// We first iterate through the fields we need to output
	// then iterate through all the slices

	// some common variables
	int gslice_len = (m+1)*(n+1);
	double *u_global = NULL;
	char * filename=NULL;
    char * decompress_command=NULL;
	FILE *fh=NULL;
	if(rank==0){
		filename = new char[512];
        decompress_command = new char [1024];
		u_global = new double[gslice_len];
		for(int i=0;i<gslice_len;i++) u_global[i]=0;
	}

	for(int j=0;j<numf;j++){
		// Change write_params output type
		out_f = output_field[j];
		wp.change_otype(out_f);
		if(rank==0) {
			sprintf(filename, "%s/%s", chk_dirname, wp.filename);
            // Decompress  the file using xz only if file exists
            sprintf(decompress_command, "if [ -f %s.xz ]; then xz -d %s.xz; fi", filename, filename);
            int decompress_err = system(decompress_command);
            if(decompress_err) {
                printf("fluid_3d::initialize_from_chk_point: Decompression command %s went wrong!\n", decompress_command);
            }
			fh = fopen(filename, "rb");
            if(fh==NULL) {
                printf("fluid_3d::initialize_from_chk_point: File %s cannot be read!\n", filename);
                err = 1;
            }
		}

        // If file wasn't read in correctly, we return error status.
        MPI_Bcast(&err, 1, MPI_INT, 0, MPI_COMM_WORLD);
        if(err>0) return err;

		for(int i=0;i<o;i++) {
			wp.change_point(i);
			if(rank==0) {
				size_t bs = fread(u_global, sizeof(double), gslice_len,fh);
                if(bs==0){
                    printf("fluid_3d::initialize_from_chk_point: Checkpoint file is empty!\n");
                    err = 2;
                }
			}
			// put back to each proc
			init_from_slice(wp, u_global);
			MPI_Barrier(grid->cart);
		}
		if(rank==0) fclose(fh);
	}
	if(rank==0){
		delete [] u_global;
        delete [] decompress_command;
		delete [] filename;
	}

    // If any of the check point files readins are empty, we also return error status
    // after cleaning up
    MPI_Bcast(&err, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if(err>0) return err;

	communicate<1>();	  // fluid vels
	communicate<2>();	  // pressure
	communicate<8>();	  // obj id
	communicate<16|32>();  // refx, refxpred

    // Use memory intense extrapolation for now
    init_extrapolate(false);
#if defined(VAR_DEN)
    update_rho();
#endif
    return 0;
}

void fluid_3d::init_from_slice(write_params wp, double * g_val){
	int point = wp.point;

	bool self_involved=false;
	double * send_recv_buf = cbuf.buf;

	// make array of sender, length, dimension in global array
	int mp = grid->mp;
	int np = grid->np;
	// needed for data assembly
	int global_len = m;
	int local_dim1 = sm;
	int local_dim2 = sn;

	if(wp.corner_field()) {
		global_len = fem_grid.m;
		local_dim1 = fem_grid.sm;
		local_dim2 = fem_grid.sn;
	}

	// array and field labels, this will keep track of
	// where in the global slice a proc needs
	int tag_n_len [5];
	const int STAG = 0;
	const int SLEN = 1;
	const int LDIM1 = 2;
	const int LDIM2 = 3;
	const int GDIM1 = 4;
	const int SENDR = 5;

	// find the processor index containing slice along perpendicular dimension
	int slice_ind = -1;
	if(point>=ak && point<bk) slice_ind=grid->kp;

	// this bit makes sure if master node doesn't own a region of the slice
	// it still knows which procs to send and recv stuff
	double tmp_ind_local = (slice_ind>=0)?((double) slice_ind / double (mp*np)):0;
	double tmp_ind_master;
	MPI_Reduce(&tmp_ind_local, &tmp_ind_master, 1, MPI_DOUBLE, MPI_SUM, 0, grid->cart);
	if(rank==0 && slice_ind==-1) slice_ind=(int) tmp_ind_master;

	// just random tag to send tag
	const int TAGTAG = 2;
	int slice_len, slice_tag;

	// if this is a processor that contains the slice, we prepare the node
	// for receiving by sending information about what to receive
	if (grid->kp == slice_ind) {
		// load the info on to the writing array
		slice_len = local_dim1*local_dim2;
		// starting position is computed with respect to the global slice
		// for a z slice, my tag is ai+aj*m
		slice_tag = ai + aj * global_len;
		// First send the slice tags to master node, with the tag -1
		tag_n_len[STAG] = slice_tag;
		tag_n_len[SLEN] = slice_len;
		tag_n_len[LDIM1] = local_dim1;
		tag_n_len[LDIM2] = local_dim2;
		tag_n_len[GDIM1] = global_len;

		if (rank==0) {self_involved=true;}
		else {
			MPI_Send(tag_n_len,5,MPI_INT,0,TAGTAG,grid->cart);
			// then get ready to receive the slice data to master node with slice tags
			MPI_Irecv(send_recv_buf, slice_len, MPI_DOUBLE, 0, slice_tag, grid->cart, reqs);
			MPI_Waitall(1,reqs,stats);
		}
	}
	// now that all's received, we send from the master
	if (rank == 0) {
		// define array for cartesian processor coords
		int coords[3] = {0,0,slice_ind};

		int **recv_list = new int* [mp*np];
		// now for each processor
		for (int i = 0; i < mp*np; i++) {
			// make array of size six because we also store the order
			// in which these slice messages were received, i.e. SENDR
			recv_list[i] = new int[6];
			for (int j = 0; j<6; j++) recv_list[i][j] = 0;
		}
		int ** sublist = recv_list;
		if(self_involved){
			// we in fact did not send from master to master
			// manually load data
			for (int i=0;i<5;i++) recv_list[0][i] = tag_n_len[i];
			recv_list[0][5] = 0;
			sublist++;
		}

		// iterate through all the processors that intersect
		// the slice, and recv tag and length info
		// Afterward, send out the part of the slice requested
		int recver, comm_co=0;
		for(int i = 0; i < mp; i++) {
			for(int j = 0; j < np; j++) {

				// record the processor coords and get rank
				coords[0] = i;
				coords[1] = j;
				MPI_Cart_rank(grid->cart,coords,&recver);
				if (recver==0){
					int *master_list = recv_list[0];
					copy_chk_slice_to_buf(g_val, send_recv_buf, master_list);
					send_recv_buf +=  master_list[SLEN];
					continue;
				}

				// get pointer to list for storing info
				*(*sublist+SENDR) = recver;
				// fill in rest of list with info from sender
				MPI_Recv(*sublist,5,MPI_INT,recver,TAGTAG,grid->cart,MPI_STATUS_IGNORE);

				// includes slice data tag and length of slice data
				int pos = *(*sublist+STAG);
				int len = *(*sublist+SLEN);
				// pack the data into buffer to send to processers
				copy_chk_slice_to_buf(g_val, send_recv_buf, *sublist);
				MPI_Isend(send_recv_buf,len,MPI_DOUBLE,recver,pos,grid->cart,reqs+comm_co);
				send_recv_buf += len;
				sublist++;comm_co++;
			}
		}
		// wait to finish everything, then copy out of buffer
		MPI_Waitall(comm_co,reqs,stats);

		// delete the new array we made
		for (int i = 0; i < mp*np; i++) if (recv_list[i]!=NULL) delete [] recv_list[i];
		if(recv_list!=NULL) delete[] recv_list;
	}

	// copy from buffer to the data structure, depending on which field was being checkpointed
	if (self_involved || grid->kp == slice_ind) {
		copy_chk_to_field(wp);
	}

}

void fluid_3d::copy_chk_slice_to_buf(double *g_val, double *comm_buf, int *info){
	//const int STAG = 0;
	//const int SLEN = 1;
	//const int LDIM1 = 2;
	//const int LDIM2 = 3;
	//const int GDIM1 = 4;
	//const int SENDR = 5;
	int start = info[0];
	int lm = info[2], ln = info[3];
	double *gp = g_val+start, *cp = comm_buf;
	for(int j=0;j<ln;j++) for (int i=0;i<lm;i++,cp++){
		*(cp)= gp[i+j*m];
	}
}

void fluid_3d::copy_chk_to_field(write_params wp){
	// take log base2 0,1,2,7 are (u,v,w,p), 9,10,11 are refx,refy,refz
	int o_type = wp.get_otype();
	int k = wp.get_point()-ak;
	double *cp = cbuf.buf;
	if(o_type >=0 && o_type <=2){

		for(int j=0;j<sn;j++) for(int i=0;i<sm;i++){
			u0[index(i,j,k)].vel[o_type] = *(cp++);
		}

	} else if(o_type == 7){

		for(int j=0;j<fem_grid.sn;j++) for(int i=0;i<fem_grid.sm;i++){
			u0[index(i,j,k)].p = *(cp++);
		}

	} else if(o_type>=9 && o_type<=11){
		o_type -= 9;
		for(int j=0;j<sn;j++) for(int i=0;i<sm;i++){
			rm0[index(i,j,k)].x[o_type] = rm0[index(i,j,k)].xpred[o_type] = *cp;
			cp++;
		}

	} else if(o_type == 12){

		for(int j=0;j<sn;j++) for(int i=0;i<sm;i++){
			rm0[index(i,j,k)].c = (unsigned int) (*(cp++));
		}

	} else {
		printf("fluid_3d::copy_chk_to_field(): Check point field %d not recognized.\n", o_type);
		p_fatal_error("check point file initialization error.", 1);
	}
}
