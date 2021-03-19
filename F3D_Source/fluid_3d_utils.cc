#include "fluid_3d.hh"

/** Record the maximum length of a solid in x direction on mid z plane*/
void fluid_3d::record_magnitude(){
	double minx = 10000, maxx = -10000;
	double gminx, gmaxx;
	for(int k=0;k<so;k++) for(int j=0;j<sn;j++) for(int i=0;i<sm;i++){
		int ind = index(i,j,k);
		if(rm0[ind].oid()!=0) {
			minx = (lx0[i]<minx)?lx0[i]:minx;
			maxx = (lx0[i]>maxx)?lx0[i]:maxx;
		}
	}
	MPI_Reduce(&minx, &gminx, 1, MPI_DOUBLE, MPI_MIN,0,grid->cart);
	MPI_Reduce(&maxx, &gmaxx, 1, MPI_DOUBLE, MPI_MAX, 0, grid->cart);
	if(rank==0) printf("%g %8.6g\n", time, gmaxx-gminx);
}

/** Check if any of the passed pointers is null. */
bool fluid_3d::check_null(ref_map *r0, ref_map *r1, ref_map *r2, ref_map *r3,
			ref_map *r4, ref_map *r5, ref_map *r6, ref_map *r7){
	bool anynull = false;
	if(r0==NULL || r1==NULL || r2 ==NULL || r3 ==NULL
		|| r4==NULL || r5==NULL || r6 == NULL || r7==NULL) anynull=true;
	return anynull;
}

/** Get the reference map, either on primary or on multimap grid,
 * that matches the object id oid.
 */
ref_map * fluid_3d::get_refmap(int ind, unsigned int oid){
    if(ind<0 || ind>=smn4*so4) return NULL;

	if (rm_mem[ind].c == oid) return rm_mem+ind;
	else {
		for(int i=0;i<extraps.n[ind];i++){
			if( (extraps.f[ind][i].c & omask) == oid)
				return extraps.f[ind]+i;
		}
		return NULL;
	}
}
/** Get the reference map in the previous extrapolation. If on the primary grid, then we can just get the current value. */
ref_map * fluid_3d::get_refmap_prev(int ind, unsigned int oid){
    if(ind<0 || ind>=smn4*so4) return NULL;

	if (rm_mem[ind].c == oid) return rm_mem+ind;
	else {
		for(int i=0;i<temp_extraps.n[ind];i++){
			if( (temp_extraps.f[ind][i].c & omask) == oid)
				return temp_extraps.f[ind]+i;
		}
		return NULL;
	}
}

/** Get the reference map, either on primary or on multimap grid,
 * that matches the object id oid.
 */
ref_map * fluid_3d::get_refmap(int ind, unsigned int oid, const int (&pos)[3]){
    if (!within_local_domain(pos[0], pos[1], pos[2]) || ind<0 || ind>=smn4*so4) return NULL;

	if (rm_mem[ind].c == oid) return rm_mem+ind;
	else {
		for(int i=0;i<extraps.n[ind];i++){
			if( (extraps.f[ind][i].c & omask) == oid)
				return extraps.f[ind]+i;
		}
		return NULL;
	}
}
/** Get the reference map in the previous extrapolation. If on the primary grid, then we can just get the current value. */
ref_map * fluid_3d::get_refmap_prev(int ind, unsigned int oid, const int (&pos)[3]){
    if (!within_local_domain(pos[0], pos[1], pos[2]) || ind<0 || ind>=smn4*so4) return NULL;

	if (rm_mem[ind].c == oid) return rm_mem+ind;
	else {
		for(int i=0;i<temp_extraps.n[ind];i++){
			if( (temp_extraps.f[ind][i].c & omask) == oid)
				return temp_extraps.f[ind]+i;
		}
		return NULL;
	}
}

/** Given index i, j, k, in x, y, z direction
 * compute index from first non-ghost node.
 */
int fluid_3d::index(int i, int j, int k){
	return i+sm4*j+smn4*k;
}

/** Given index ind, unpack it into 3 indices in x,y,z direction. */
void fluid_3d::unpack_index(int ind, int &i, int &j, int &k){
	k = ind/smn4;
	j = (ind%smn4)/ sm4;
	i = ind%sm4;
}

/** Calculate the fraction of each solid, normalized if fraction sum exceeds 1
 */
void fluid_3d::solid_fractions(int ind, std::vector<solid_frac> &s_fracs){
	ref_map *prime = rm0+ind;
	int prime_id = prime->oid();
	int sum_of_fracs = 0;
	// grab the primary object fraction here
	if(prime_id > 0) {
		double tmp =  prime->frac(mgmt);
		solid_frac sf(prime,tmp);
		s_fracs.push_back(sf);
	}
	int N = extraps.n0[ind];
	for(int i=0;i<N;i++){
		ref_map *ext = extraps.f0[ind] + i;
		double tmp = ext->frac(mgmt);
		solid_frac sf(ext,tmp);
		s_fracs.push_back(sf);
		sum_of_fracs += tmp;
	}
	if(sum_of_fracs>1) {
		for(unsigned i=0;i<s_fracs.size();i++){
			s_fracs[i].normalize(sum_of_fracs);
		}
	}
}

/* ###################### WALL UTILITIES ######################## */
/*
void fluid_3d::use_walls() {
	fwalls=true;
	fwall_v[0] = fwall_v[1] = fwall_v[2] = 0;
}

void fluid_3d::use_walls(double (&fwv)[3]) {
	fwalls = true;
	fwall_v[0] = fwv[0];
	fwall_v[1] = fwv[1];
	fwall_v[2] = fwv[2];
}

void fluid_3d::fake_walls(int type) {
	comm_info *ct = comm_table;
	field *u;
	bool to_zero[6] = {ai == 0,bi == grid->m,aj == 0,bj == grid->n,ak == 0,bk == grid->o};
	int tag,tags[6] = {12,14,10,16,4,22};

	for (int t = 0; t < 6; t++) if (to_zero[t]) {
		tag = tags[t];
		field *ur0 = u_mem + ct[tag].s0;
		for (int k=0; k<ct[tag].o; k++) {
			for (int j=0; j<ct[tag].n; j++) {
				for (int i=0; i<ct[tag].m; i++){
					u = ur0 + i + j*sm4 + k*smn4;
					for (int c = 0; c < 3; c++) {
						switch (type) {
							case 0: break;
							case 2: u->vel[c] = fwall_v[c]; break;
							case 1: for (int f = 0; f < 6; f++) {u->fvel[f][c] = fwall_v[c];} break;
							default:
								p_fatal_error("fluid_3d::fake_walls(int):"
										" undefined type",1);
						}
					}
				}
			}
		}
	}
}
*/

/** solid mechanics/extrapolation functions.
 * This is essentially get_refmap */
/*
ref_map* fluid_3d::rm_step(ref_map* home,int del,int ind) {
	ref_map *rp = home->step(del,ind,extraps.f0,extraps.n0,par,rm0);
	if (rp == NULL) {
		puts("fluid_3d::rm_step(): Failed to find a pointer at given position.");
		p_fatal_error("Failure to find a reference map pointer, likely on multimap object.",1);
	}
	return rp;
}
*/

/* ###################### DEBUGGING UTILITIES ######################## */
/** Attach debugger */
void fluid_3d::attach_debug(int proc) {
	// from open-mpi.org/faq/?category=debugging
	if (grid->rank == proc) {
		volatile int i = 0;
		char hostname[256];
		gethostname(hostname,sizeof(hostname));
		printf("PID %d on %s ready to attach\n",getpid(),hostname);
		fflush(stdout);
		while (i == 0) {
			sleep(5);
		}
	}
	MPI_Barrier(grid->cart);
}

void fluid_3d::output_phi(int obj_num,char *fname,int type) {

	// suuuuper inefficient - print in order
	for (int kp = 0; kp < grid->op; kp++) {
		for (int jp = 0; jp < grid->np; jp++) {
			for (int ip = 0; ip < grid->mp; ip++) {

	if (ip != grid->ip || jp != grid->jp || kp != grid->kp) {
		MPI_Barrier(grid->cart);
	} else {

	FILE *fh;
	if (ip == jp && jp == kp && kp == 0) fh = p_safe_fopen(fname,"w");
	else fh = p_safe_fopen(fname,"a");
	double val,vx,vy,vz,vxp,vyp,vzp;
	for (int k = 0; k < so4; k++)
		for (int j = 0; j < sn4; j++)
			for (int i = 0; i < sm4; i++) {
		bool found = false;
		int ind = index(i,j, k);
		if (rm_mem[ind].oid() == obj_num) {
			val = rm_mem[ind].phi(mgmt,type);
			vxp = rm_mem[ind].xpred[0];
			vyp = rm_mem[ind].xpred[1];
			vzp = rm_mem[ind].xpred[2];
			vx = rm_mem[ind].x[0];
			vy = rm_mem[ind].x[1];
			vz = rm_mem[ind].x[2];
			found = true;
		} else {
			for (int jj = 0; jj < extraps.n[ind]; jj++) {
				if (extraps.f[ind][jj].oid() == obj_num) {
					val = extraps.f[ind][jj].phi(mgmt,type);
					vxp = extraps.f[ind][jj].xpred[0];
					vyp = extraps.f[ind][jj].xpred[1];
					vzp = extraps.f[ind][jj].xpred[2];
					vx = extraps.f[ind][jj].x[0];
					vy = extraps.f[ind][jj].x[1];
					vz = extraps.f[ind][jj].x[2];
					found = true;
				}
			}
		}
		if(found) {
			fprintf(fh,"%d %d %d %d %.9g %.9g %.9g %.9g %.9g %.9g %.9g\n",grid->rank,ai + i-2,aj+j-2,ak+k-2,val,vxp,vyp,vzp,vx,vy,vz);
		}
		else{
			fprintf(fh, "%d %d %d %d %.9g %.9g %.9g %.9g %.9g %.9g %.9g\n",grid->rank,ai+i-2,aj+j-2,ak+k-2,100.0,100.,100.,100.,100.,100.,100.);
		}
	}

	fclose(fh);
	MPI_Barrier(grid->cart);
	}
			}
		}
	}

}

void fluid_3d::debug_dump(const char *filename){

    puts("debug_dump:\nBottom ghost nodes");
	for(int k=0;k<2;k++){
		for(int j=0;j<sn4;j++){
			for(int i=0;i<sm4;i++){
				int ind = index(i,j,k);
				ref_map &f = rm_mem[ind];
                double lx = 100.;
                if(f.oid() > 0) lx = f.x[2];
                //printf("Debug dump: %s: nt=%d index %d, primary (%d %d %d) suspicious rm x=(%g %g %g) xpred=(%g %g %g) phi %g\n",
                //filename, nt, ind, i, j, k, f.x[0], f.x[1], f.x[2], f.xpred[0], f.xpred[1], f.xpred[2], f.phi(mgmt));
				if(extraps.n[ind]>0) {
					ref_map &ef = extraps.f[ind][0];
                    //printf("Debug dump: %s: nt=%d index %d, extraps (%d %d %d) suspicious layer %d x=(%g %g %g) xpred=(%g %g %g) phi %g\n",
                    //filename, nt, ind, i, j, k, ef.lid(), ef.x[0], ef.x[1], ef.x[2], ef.xpred[0], ef.xpred[1], ef.xpred[2], ef.phi(mgmt));
                    lx = ef.x[2];
				}
                printf("%g ", lx);
			}
            printf("\n");
		}
	}

    puts("debug_dump:\nTop none ghost nodes");
	for(int k=so;k<so+2;k++){
		for(int j=0;j<sn4;j++){
			for(int i=0;i<sm4;i++){
				int ind = index(i,j,k);
				ref_map &f = rm_mem[ind];
                double lx = 100.;
                if(f.oid() > 0) lx = f.x[2];
                //printf("Debug dump: %s: nt=%d index %d, primary (%d %d %d) suspicious rm x=(%g %g %g) xpred=(%g %g %g) phi %g\n",
                //filename, nt, ind, i, j, k, f.x[0], f.x[1], f.x[2], f.xpred[0], f.xpred[1], f.xpred[2], f.phi(mgmt));
				if(extraps.n[ind]>0) {
					ref_map &ef = extraps.f[ind][0];
                    //printf("Debug dump: %s: nt=%d index %d, extraps (%d %d %d) suspicious layer %d x=(%g %g %g) xpred=(%g %g %g) phi %g\n",
                    //filename, nt, ind, i, j, k, ef.lid(), ef.x[0], ef.x[1], ef.x[2], ef.xpred[0], ef.xpred[1], ef.xpred[2], ef.phi(mgmt));
                    lx = ef.x[2];
				}
                printf("%g ", lx);
			}
            printf("\n");
		}
	}
	/*
	if(rank==0) {
		char *local_fn = new char[256];
		sprintf(local_fn, "%s_X", filename);
		FILE *fx = p_safe_fopen(local_fn,"w");
		sprintf(local_fn, "%s_Y", filename);
		FILE *fy = p_safe_fopen(local_fn,"w");
		sprintf(local_fn, "%s_Z", filename);
		FILE *fz = p_safe_fopen(local_fn,"w");
		sprintf(local_fn, "%s_P", filename);
		FILE *fp = p_safe_fopen(local_fn,"w");

		int zz = 19;
		for(int j=2;j<sn4-1;j++){
			for(int i=0;i<2;i++){
				int ind = i+sm4*j+smn4*zz;
				//field &f = u0[ind];
				ref_map *rf = get_refmap(ind, 1);
				fprintf(fp, "Index (%d %d %d) Pos (%e %e %e)\n", i, j, zz, xx, yy, z0[zz]);
				fprintf(fp, "Fluid velocity (%f %f %f) pressure %f faces:\n", f.vel[0], f.vel[1], f.vel[2], f.p);
				for(int k=0;k<6;k++) fprintf(fp, "(%f %f %f)\n", f.fvel[k][0], f.fvel[k][1], f.fvel[k][2]);
				fprintf(fp, "Primary oid %d lid %d (%e %e %e) pred (%e %e %e)\n", rf.oid(), rf.lid(), rf.x[0], rf.x[1], rf.x[2] ,rf.xpred[0], rf.xpred[1], rf.xpred[2]);
				for(int l=0;l<extraps.n0[ind];l++){
					ref_map &ext_rf = extraps.f0[ind][l];
					fprintf(fp, "Extrapolated refmap oid %d lid %d (%e %e %e) pred (%e %e %e)\n", ext_rf.oid(), ext_rf.lid(), ext_rf.x[0], ext_rf.x[1], ext_rf.x[2] ,ext_rf.xpred[0], ext_rf.xpred[1], ext_rf.xpred[2]);
				}
				fputc('\n', fp);
				if(rf!=NULL) {
					if(i==0) printf("I'm (%d, %d, %d) my refmap(%f %f %f)\n", i, j, zz, rf->x[0], rf->x[1], rf->x[2]);
					fprintf(fx, "%17.14f ", rf->x[0]);
					fprintf(fy, "%17.14f ", rf->x[1]);
					fprintf(fz, "%17.14f ", rf->x[2]);
					fprintf(fp, "%17.14f ", rf->phi(mgmt));
				}
				else{
					fprintf(fx, "%17d ", 100);
					fprintf(fy, "%17d ", 100);
					fprintf(fz, "%17d ", 100);
					fprintf(fp, "%17d ", 100);
				}

				ind = i+sm4*j+smn4*zz+32;
				rf = get_refmap(ind, 1);

				if(rf!=NULL) {
					fprintf(fx, "%17.14f ", rf->x[0]);
					fprintf(fy, "%17.14f ", rf->x[1]);
					fprintf(fz, "%17.14f ", rf->x[2]);
					fprintf(fp, "%17.14f ", rf->phi(mgmt));
				}
				else{
					fprintf(fx, "%17d ", 100);
					fprintf(fy, "%17d ", 100);
					fprintf(fz, "%17d ", 100);
					fprintf(fp, "%17d ", 100);
				}

			}
				fputc('\n', fx);
				fputc('\n', fy);
				fputc('\n', fz);
				fputc('\n', fp);
		}
		fclose(fx);fclose(fy);fclose(fz);fclose(fp);
		delete [] local_fn;
	}
	*/
}

void fluid_3d::check_nans(const int id){

    puts("===============================");
    printf("Rank %d Check nans %d\n", rank, id);
	for(int k=0;k<so4;k++) for(int j=0;j<sn4;j++) for(int i=0;i<sm4;i++){
        ref_map *rp = get_refmap(index(i, j, k), 1);
        if(rp!=NULL && (
           std::isnan (rp->x[0]) ||
           std::isnan (rp->x[1]) ||
           std::isnan (rp->x[2])) ) {
            printf("rank %d Nan ref map at (%d %d %d) oid 1\n", rank, i, j, k);
        }

        if(rp!=NULL && (
           std::isnan (rp->xpred[0]) ||
           std::isnan (rp->xpred[1]) ||
           std::isnan (rp->xpred[2])) ) {
            printf("rank %d Nan ref map xpred at (%d %d %d) oid 1\n", rank, i, j, k);
        }

        field &f = u_mem[index(i,j,k)];
        if(std::isnan (f.vel[0]) ||
           std::isnan (f.vel[1]) ||
           std::isnan (f.vel[2]) ) {
            printf("rank %d Nan fluid vel at (%d %d %d) oid 1\n", rank, i, j, k);
        }

#if defined(VAR_DEN)
        if(lrho[index(i,j,k)] == 0) {
            printf("rank %d Rho field is 0 at (%d %d %d) \n", rank, i, j, k);
        }
#endif

        for(int ff=0;ff<6;ff++) {
            if(std::isnan(f.fvel[ff][0]) ||
               std::isnan(f.fvel[ff][1]) ||
               std::isnan(f.fvel[ff][2])) printf("Rank %d (%d %d %d) fvel %d is nan (%g %g %g)\n", rank, i, j, k, ff, f.fvel[ff][0], f.fvel[ff][1], f.fvel[ff][2]);
        }
    }
    
    printf("Done rank %d Check nans %d\n", rank, id);
    puts("===============================\n");
}

void fluid_3d::check_symmetry(){
	if(grid->procs >1) return;
	puts("fluid_3d::check_symmetry() x");
	printf("Timestep %d\n", nt);
	int within[2]={0,0};
	int sym[2]={0,0};
	for(int k=0;k<so;k++) for(int j=0;j<sn;j++) for(int i=0;i<sm/2;i++){
		int eid = index(i,j,k);
		int eid_symm = index(sm-1-i, j, k);
		ref_map &r1 = rm0[eid];
		ref_map &r2 = rm0[eid_symm];
		int oid1=r1.oid(), oid2=r2.oid();
		if(oid1>0) within[oid1-1] ++;
		if(oid2>0) within[oid2-1] ++;
		double phi1 = r1.phi(mgmt, 1), phi2 = r2.phi(mgmt,1);
		if( (oid1 != oid2) || (oid1>0 && oid2>0 && oid1==oid2 && fabs(phi1-phi2)>small_number)) {
            /*
			printf("oid %d (%d,%d,%d) and oid %d (%d,%d,%d) are not symmetric about mid yz plane. Phi=%16.14g %16.14g\n", oid1, i,j,k, oid2, sm-1-i, j, k, phi1, phi2);
			printf("PT1: (%16.14g, %16.14g, %16.14g), pred(%16.14g, %16.14g, %16.14g); %16.14g\n"
				"MIR: (%16.14g, %16.14g, %16.14g), pred(%16.14g, %16.14g, %16.14g); %16.14g\n"
				,
				r1.x[0], r1.x[1], r1.x[2],
				r1.xpred[0], r1.xpred[1], r1.xpred[2], r1.phi(mgmt),

				r2.x[0], r2.x[1], r2.x[2],
				r2.xpred[0], r2.xpred[1], r2.xpred[2], r2.phi(mgmt)
				);
            */
			sym[oid1-1]++;
			sym[oid2-1]++;

		} else if (oid1>0 && oid2>0 && oid1==oid2 && fabs(phi1-phi2)<small_number) {
			sym[oid1-1]+=2;
		}
        field &f1 = u0[eid];
        field &f2 = u0[eid_symm];
        double dvel_diff = 0;
        for(int s=0;s<3;s++){
            dvel_diff += (fabs(f1.dvel[s]) - fabs(f2.dvel[s])) * (fabs(f1.dvel[s]) - fabs(f2.dvel[s]));
        }
        dvel_diff = sqrt(dvel_diff);
        if(dvel_diff> small_number) {
            printf("point (%d %d %d) dvel %16.14g %16.14g %16.14g\n", i, j, k, f1.dvel[0], f1.dvel[1], f1.dvel[2]);
            printf("mirro (%d %d %d) dvel %16.14g %16.14g %16.14g\n", sm-1-i, j, k, f2.dvel[0], f2.dvel[1], f2.dvel[2]);
        }
/*
        if(fabs(f1.p-f2.p) >small_number || fabs(f1.dp-f2.dp)>small_number){
            printf("point (%d %d %d) pressure %16.14g\n", i, j, k, f1.p);
            printf("point (%d %d %d) pressure %16.14g\n", sm-1-i, j, k, f2.p);
            printf("point (%d %d %d) incr pressure %16.14g\n", i, j, k, f1.dp);
            printf("point (%d %d %d) incr pressure %16.14g\n", sm-1-i, j, k, f2.dp);
            puts("");
        }
        if(fabs(f1.q-f2.q) >small_number){
            printf("point (%d %d %d) mac pressure %16.14g\n", i, j, k, f1.q);
            printf("point (%d %d %d) mac pressure %16.14g\n", sm-1-i, j, k, f2.q);
            puts("");
        }
        if(fabs(lrho[G0+eid] - lrho[G0+eid_symm]) > small_number){
            printf("point (%d %d %d) density %16.14g\n", i, j, k, lrho[G0+eid]);
            printf("point (%d %d %d) density %16.14g\n", sm-1-i, j, k, lrho[G0+eid_symm]);
            puts("");
        }
*/
	}
	printf("Within obj 1 %d points, %d symmetric; obj 2 %d points, %d symmetric.\n", within[0], sym[0], within[1], sym[1]);
	if(within[0] != sym[0] || within[1] != sym[1]) p_fatal_error("Symmetry is violated, what is up?!\n", 1);

	puts("fluid_3d::check_symmetry() y");
	for(int k=0;k<so;k++) for(int j=0;j<sn;j++) for(int i=0;i<sm/2;i++){
		int eid = index(i,j,k);
		int eid_symm = index(i, sn-1-j, k);
        field &f1 = u0[eid];
        field &f2 = u0[eid_symm];
        double dvel_diff = 0;
        for(int s=0;s<3;s++){
            dvel_diff += (fabs(f1.dvel[s]) - fabs(f2.dvel[s])) * (fabs(f1.dvel[s]) - fabs(f2.dvel[s]));
        }
        dvel_diff = sqrt(dvel_diff);
        if(dvel_diff> small_number) {
            printf("point (%d %d %d) dvel %16.14g %16.14g %16.14g\n", i, j, k, f1.dvel[0], f1.dvel[1], f1.dvel[2]);
            printf("mirro (%d %d %d) dvel %16.14g %16.14g %16.14g\n", i, sn-1-j, k, f2.dvel[0], f2.dvel[1], f2.dvel[2]);
        }
	}
}

/** Output in the macro.dat file
 */
double fluid_3d::div_u( double_int &gmx_ext, double_int &gmn_ext ){
	double dhsp[3] = {dxsp, dysp, dzsp};

    // Here we figure out the grid L2 norm of the divergence
    // Defined as in pressure poisson problem
    // TODO is this kosher?
    double local_l2=0, global_l2=0;
    //double var_local_l2=0, var_global_l2=0;
    double max=0, mp_global_max=0;
    mg_pres->reg[0]->clear_r_field();
    fill_pp_rhs();

    int mg_ind = 0;
	for(int k=0;k<fem_grid.so;k++){
        double norm_k = (k==0 || k==fem_grid.so-1)?0.5:1.0;
        if(fem_grid.z_prd) norm_k = 1.0;

        for(int j=0;j<fem_grid.sn;j++){
            double norm_j = (j==0 || j==fem_grid.sn-1)?0.5:1.0;
            if(fem_grid.y_prd) norm_j = 1.0;

            for(int i=0;i<fem_grid.sm;i++){
                double norm_i = (i==0 || i==fem_grid.sm-1)?0.5:1.0;
                if(fem_grid.x_prd) norm_i = 1.0;

                double wt = norm_i * norm_j * norm_k;

                // Added the weight on the integral of div u squared
                local_l2 += wt*(mg_pres->reg[0]->r0[mg_ind])*(mg_pres->reg[0]->r0[mg_ind]);
                if( fabs(mg_pres->reg[0]->r0[mg_ind]) > max) max = fabs(mg_pres->reg[0]->r0[mg_ind]);

                mg_ind ++;
                /*
                int eid = index(i, j, k);
           double tmp_divu = 0.25*dxsp*(u0[eid].vel[0] - u0[eid-1].vel[0] + u0[eid-sm4].vel[0] - u0[eid-sm4-1].vel[0] +u0[eid-smn4].vel[0] - u0[eid-smn4-1].vel[0] + u0[eid-smn4-sm4].vel[0] - u0[eid-smn4-sm4-1].vel[0]);
           tmp_divu += 0.25*dysp*(u0[eid].vel[1] - u0[eid-sm4].vel[1] + u0[eid-1].vel[1] + u0[eid-sm4-1].vel[1] +u0[eid-smn4].vel[1] - u0[eid-smn4-sm4].vel[1] + u0[eid-smn4-1].vel[1] - u0[eid-smn4-sm4-1].vel[1]);
           tmp_divu += 0.25*dzsp*(u0[eid].vel[2] - u0[eid-smn4].vel[2] + u0[eid-1].vel[2] - u0[eid-smn4-1].vel[2] +u0[eid-sm4].vel[2] - u0[eid-sm4-smn4].vel[2] + u0[eid-1-sm4].vel[2] - u0[eid-smn4-sm4-1].vel[2]);
                var_local_l2 += tmp_divu*tmp_divu;
                */
            }
        }
    }

    // we can use dx*dy*dz to normalize
    // because no matter how we integrate, the whole domain is covered
    local_l2 /= (mg_pres->reg[0]->m * mg_pres->reg[0]->n * mg_pres->reg[0]->o);
    //var_local_l2 /= (mg_pres->reg[0]->m * mg_pres->reg[0]->n * mg_pres->reg[0]->o);

    MPI_Allreduce(&local_l2, &global_l2, 1, MPI_DOUBLE, MPI_SUM, grid->cart);
    //MPI_Allreduce(&var_local_l2, &var_global_l2, 1, MPI_DOUBLE, MPI_SUM, grid->cart);
	MPI_Allreduce(&max,&mp_global_max,1,MPI_DOUBLE,MPI_MAX,grid->cart);

    //double fac = dt/(12 * dx *dx *mgmt->fm.rho);
    global_l2 = sqrt(global_l2);
    //var_global_l2 = sqrt(var_global_l2);
    //printf("Calculation of approximate div u from RHS %g\n", global_l2);
    //printf("Calculation of lacal max from RHS %g, global max %g\n", max, mp_global_max);

    // The following routine figures out the max norm of divergence defined by MAC velocities
    double_int max_divu(-1e200, 0);
    double_int min_divu( 1e200, 0);
    int maxindex[3] = {0,0,0};
    int minindex[3] = {0,0,0};
	for(int k=0;k<so;k++) for (int j=0;j<sn;j++) for (int i=0;i<sm;i++){
		int elem_id = index(i,j,k);
		field *elem = u0+elem_id;
		double div=0;

		// centered difference to get grad u vector at cell center, from edge velocities
        // HERE we used mac velocities
		for(int f=0;f<3;f++) div += dhsp[f] * (elem->fvel[2*f+1][f] - elem->fvel[2*f][f]);
        if(div>max_divu.value) {
            maxindex[0]=i; maxindex[1]=j; maxindex[2]=k;
            max_divu.value = div;
            max_divu.rank = rank;
        }
        if(div<min_divu.value) {
            minindex[0]=i; minindex[1]=j; minindex[2]=k;
            min_divu.value = div;
            min_divu.rank = rank;
        }
	}

    // Add the relative lower index to the i,j,k, so we know the global index of problematic points
    maxindex[0] = ai+maxindex[0];
    maxindex[1] = aj+maxindex[1];
    maxindex[2] = ak+maxindex[2];

    minindex[0] = ai+minindex[0];
    minindex[1] = aj+minindex[1];
    minindex[2] = ak+minindex[2];

    double_int global_max = {-1e200,-1};
    double_int global_min = { 1e200, 1};
    int gmaxindex [3];
    int gminindex [3];

    MPI_Allreduce(&max_divu, &global_max, 1, MPI_DOUBLE_INT, MPI_MAXLOC, grid->cart);
    MPI_Allreduce(&min_divu, &global_min, 1, MPI_DOUBLE_INT, MPI_MINLOC, grid->cart);

    // Figure out who had the maximum divergence
    if  (rank == global_max.rank && rank!=0) {
        MPI_Request send_request;
        MPI_Status  status;
        MPI_Isend(maxindex, 3, MPI_INT, 0, 0, grid->cart, &send_request);
        MPI_Wait(&send_request, &status);
    }

    // Figure out who had the miminmum divergence
    if  (rank == global_min.rank && rank!=0) {
        MPI_Request send_request;
        MPI_Status  status;
        MPI_Isend(minindex, 3, MPI_INT, 0, 0, grid->cart, &send_request);
        MPI_Wait(&send_request, &status);
    }

    if(rank == 0 && global_max.rank!=0) {
        MPI_Status  status;
        MPI_Recv(gmaxindex, 3, MPI_INT, global_max.rank, 0, grid->cart, &status);
    }

    if(rank == 0 && global_min.rank!=0) {
        MPI_Status  status;
        MPI_Recv(gminindex, 3, MPI_INT, global_min.rank, 0, grid->cart, &status);
    }

    if(rank == 0 && global_max.rank==0) {
        for(int i=0;i<3;i++) gmaxindex[i] = maxindex[i];
    }

    if(rank == 0 && global_min.rank==0) {
        for(int i=0;i<3;i++) gminindex[i] = minindex[i];
    }

    gmx_ext  = global_max;
    gmn_ext  = global_min;
    /*
	if(rank==0) {
        printf("Maximum divergence is %g on rank %d at (%d, %d, %d)\n", global_max.value, global_max.rank, gmaxindex[0], gmaxindex[1], gmaxindex[2]);
        //for(int i=0;i<extrap_objs;i++) printf("Primary refmap OID %d. Extrapolated regions OID %d, layer %d\n", extraps.f0[ind][i].oid(), extraps.f0[ind][i].lid());

        printf("Minimum divergence is %g on rank %d at (%d, %d, %d)\n", global_min.value, global_max.rank, gminindex[0], gminindex[1], gminindex[2]);
        //for(int i=0;i<extrap_objs;i++) printf("In extrapolated regions of object %d, layer %d\n", extraps.f0[ind][i].oid(), extraps.f0[ind][i].lid());
    }
    */
    return global_l2;
}

void fluid_3d::sanity_check(){
    const double mx_amp = 10;
    printf("fluid_3d:: checking sanity at time step %d, time %g\n", nt, time);
    double_int gmx_ext, gmn_ext;
    div_u(gmx_ext, gmn_ext);

	for (int kk =0; kk < so4; kk++) {
		for (int jj =0; jj < sn4; jj++) {
			for (int ii =0; ii < sm4; ii++) {
				int eid = index(ii,jj,kk);
				ref_map *sfp = rm_mem+eid;
                if(sfp->oid()>0){
                    if (
                        fabs(sfp->x[0])>mx_amp ||
                        fabs(sfp->x[1])>mx_amp ||
                        fabs(sfp->x[2])>mx_amp
                        ) {
                        printf("fluid_3d:: sanity_check, rank %d wtf primary point for objet %d at local point (%d %d %d) (global index (%d %d %d) has ref map val bigger than 10.\n",
                        rank, sfp->oid(), ii-2, jj-2, kk-2, ai+ii-2, aj+jj-2, ak+kk-2);
                    }

                    for (int j = 0; j < extraps.n[eid]; j++) {
                            sfp = extraps.f[eid]+j;
                            if (
                                fabs(sfp->x[0])>mx_amp ||
                                fabs(sfp->x[1])>mx_amp ||
                                fabs(sfp->x[2])>mx_amp
                                ) {
                                printf("fluid_3d:: sanity_check, wtf %d local extrapolation for object %d, in layer %d at point (%d %d %d) (global index (%d %d %d) has ref map val bigger than 10.\n",
                                rank, sfp->oid(), sfp->lid(), ii-2, jj-2, kk-2, ai+ii-2, aj+jj-2, ak+kk-2);
                            }
                    }
                }
            }
        }
    }

}

int fluid_3d::check_extrapolated_region(int &ll){

    int global_tots, local_tots=0;
    int global_least_layer, least_layer = 100;
	for (int kk =0; kk < so4; kk++) {
		for (int jj =0; jj < sn4; jj++) {
			for (int ii =0; ii < sm4; ii++) {
				int eid = index(ii,jj,kk);
                if(extraps.n[eid]>1) {
                    local_tots += 1;
                    for(int s=0;s<extraps.n[eid];s++){
                        int ld = (extraps.f[eid][s]).lid();
                        //if(ld<=2) printf("Rank %d, at (%d, %d, %d) I found extrapolated map belonging to object %d, layer %d\n", rank, ii, jj, kk, extraps.f[eid][s].oid(), extraps.f[eid][s].lid());
                        if(ld < least_layer) least_layer = ld;
                    }
                }
            }
        }
    }
    MPI_Allreduce(&local_tots, &global_tots, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&least_layer, &global_least_layer, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    ll = global_least_layer;
    return global_tots;
}

int fluid_3d::check_wall_extrapolations(){
    int global_least_layer = -1;
    if(ai == 0 || aj == 0 || ak ==0 || bi == m-1 || bj == n-1 || bk == o-1){
        int least_layer = 100;
        // Hack right now, should check all non-periodic boundaries
        // check lower boundary
        for (int kk =0; kk < 2; kk++) {
            for (int jj =0; jj < sn4; jj++) {
                for (int ii =0; ii < sm4; ii++) {
                    int eid = index(ii,jj,kk);
                    for(int s=0;s<extraps.n[eid];s++){
                        int ld = (extraps.f[eid][s]).lid();
                        if(ld < least_layer) least_layer = ld;
                    }
                }
            }
        }
        MPI_Allreduce(&least_layer, &global_least_layer, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
        if(global_least_layer == 100) global_least_layer = -1;
    }
    return global_least_layer;
}

/** Compute a norm of the global speed, we can do 2 norm and inf norm (norm=-1). */
double fluid_3d::vel_mod(int norm){
    double global_tot = 0;
    double local_tot =0;
    if(norm == 2){
        for (int kk =0; kk < so; kk++) {
            for (int jj =0; jj < sn; jj++) {
                for (int ii =0; ii < sm; ii++) {
                    int eid = index(ii,jj,kk);
                    field &f = u0[eid];
                    double v_sq = f.vel[0]*f.vel[0] + f.vel[1]*f.vel[1] + f.vel[2]*f.vel[2];
                    local_tot += v_sq;
                }
            }
        }
        // since for L2 norm we sum up all the grid cell, we normalize by volume element
        local_tot *= dx*dy*dz;
    } else if (norm == -1) {
        for (int kk =0; kk < so; kk++) {
            for (int jj =0; jj < sn; jj++) {
                for (int ii =0; ii < sm; ii++) {
                    int eid = index(ii,jj,kk);
                    field &f = u0[eid];
                    double v_sq = f.vel[0]*f.vel[0] + f.vel[1]*f.vel[1] + f.vel[2]*f.vel[2];
                    if(v_sq > local_tot) {
                        local_tot = v_sq;
                    }
                }
            }
         }
    } else {
       if(rank==0) printf("vel_mod:: norm=%d isn't implemented in this function. Use norm=2 for L2 norm, and norm=-1 for infinity norm.\n",  norm);
    }

    MPI_Allreduce(&local_tot, &global_tot, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    global_tot = sqrt(global_tot);
    return global_tot;
}

/** Compute total momentum.
    Output in the macro.dat file.*/
void fluid_3d::total_momentum(double (&mv_ext)[3]){
    double global_mv[3] = {0,0,0};
    double mv[3] = {0,0,0};
    for (int kk =0; kk < so; kk++) {
        for (int jj =0; jj < sn; jj++) {
            for (int ii =0; ii < sm; ii++) {
                int eid = index(ii,jj,kk);
                field &f = u0[eid];
#if defined(VAR_DEN)
                for(int s=0;s<3;s++) mv[s] += f.vel[s]*lrho[eid+G0];
#else
                for(int s=0;s<3;s++) mv[s] += f.vel[s];
#endif
            }
        }
    }
    MPI_Allreduce(mv, global_mv, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    for(int s=0;s<3;s++) mv_ext[s] = global_mv[s]*dx*dy*dz;
}

/** Compute total non-dissipative energy, assuming constant density everywhere.
    If we compute potential energy, assume equi-potential planes is parallel to xy-plane.
 */
double fluid_3d::total_energy(double &pot_energy, double &kin_energy, double &elas_energy){
    // TODO add dissipation, taking into account
    // fluid viscous stress
    // solid artificial stress in the body and at the interface (a multiplier is applied)
    double gpot = 0, gkin = 0, gelas =0;
    double pot = 0, kin = 0, elas =0;
    double tot = 0, gtot = 0;
    const int strides[3] = {1, sm4, smn4};
    for (int kk =0; kk < so; kk++) {
        for (int jj =0; jj < sn; jj++) {
            for (int ii =0; ii < sm; ii++) {
                int eid = index(ii,jj,kk);
                int ind = G0+eid;

                field &f = u0[eid];
                ref_map &rf = rm0[eid];
                // first we get kinetic energy here
                for(int s=0;s<3;s++) kin += 0.5*f.vel[s]*f.vel[s];
                // then we compute potential energy in the case of sediment
                // assuming potential directed downward in z directikon
                if(spars->gravity!=0 && rf.oid()>0) {
                    double delta_z = lz0[kk] - mgmt->az;
                    pot +=  spars->gravity * delta_z;
                }

                if(rf.oid() > 0){
                    double tmp_elas = 0;
                    for(int s=0;s<3;s++) {
                        // get the energies own by this cell, 3 lower faces
                        tmp_elas += rf.elastic_energy[s];
                        // now get the energies own by the immediate nearby cells
                        // the 3 upper faces
                        ref_map *rp = get_refmap(ind+strides[s], rf.oid());
                        tmp_elas += rp->elastic_energy[s];
                        /*
                        if(rf.elastic_energy[s] < 0 || rp->elastic_energy[s] < 0) {
                            printf("total_energy: At (%d %d %d) %d face has negative elastic energy, primary %g (phi=%g, layer %d) neighbor %g (phi=%g, layer %d)\n", ii,jj,kk, s, rf.elastic_energy[s], rf.phi(mgmt), rf.lid(), rp->elastic_energy[s], rp->phi(mgmt), rp->lid());
                        }
                        */
                    }
                    //  average over the faces
                    tmp_elas /= 6.;
                    elas += tmp_elas;
                }
            }
        }
    }

    double dv = dx*dy*dz;
    elas*=dv;
    kin*=dv;
    pot*=dv;
    tot = elas + kin + pot;

    MPI_Allreduce(&elas, &gelas, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&kin, &gkin, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&pot, &gpot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&tot, &gtot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    elas_energy = gelas;
    kin_energy = gkin;
    pot_energy = gpot;
    return gtot;
}

#if defined(DEBUG)
double_int fluid_3d::max_solid_stress(int type){
    double_int local_max(0,-1);
    double_int global_max(0,-1);
    for (int kk =0; kk < so; kk++) {
        for (int jj =0; jj < sn; jj++) {
            for (int ii =0; ii < sm; ii++) {
                int eid = index(ii,jj,kk);
                for (int j = -1; j < extraps.n0[eid]; j++) {
                    ref_map *rp;

                    if( j == -1) rp=rm0+eid;
                    else rp = extraps.f0[eid] + j;

                    for(int F=0;F<3;F++){
                        double stress = 0;
                        switch (type) {
                            case 0: stress = fabs(rp->solid_stress_normal[F]); break;
                            case 1: stress = fabs(rp->solid_stress_shear[F]); break;
                        }
                        if(stress != 1.e200){
                            if(stress > local_max.value) {
                                local_max.value = stress;
                                local_max.rank = rp->oid();
                            }
                        }
                    }

                }
            }
        }
    }

    MPI_Allreduce(&local_max, &global_max, 1, MPI_DOUBLE_INT, MPI_MAXLOC, grid->cart);
    return global_max;
}
#endif

#if defined(DEBUG)
double_int fluid_3d::max_coll_stress(int type){
    double_int local_max(0,-1);
    double_int global_max(0,-1);
    for (int kk =0; kk < so; kk++) {
        for (int jj =0; jj < sn; jj++) {
            for (int ii =0; ii < sm; ii++) {
                int eid = index(ii,jj,kk);
                for (int j = -1; j < extraps.n0[eid]; j++) {
                    ref_map *rp;

                    if( j == -1) rp=rm0+eid;
                    else rp = extraps.f0[eid] + j;

                    for(int F=0;F<3;F++){
                        double stress = 0;
                        switch (type) {
                            case 0: stress = fabs(rp->coll_stress_normal[F]); break;
                            case 1: stress = fabs(rp->coll_stress_shear[F]); break;
                        }
                        if(stress != 1.e200){
                            if(stress > local_max.value) {
                                local_max.value = stress;
                                local_max.rank = rp->oid();
                            }
                        }
                    }

                }
            }
        }
    }

    MPI_Allreduce(&local_max, &global_max, 1, MPI_DOUBLE_INT, MPI_MAXLOC, grid->cart);
    return global_max;
}
#endif

#if defined(DEBUG)
double_int fluid_3d::max_coll_traction(){
    double_int local_max(0,-1);
    double_int global_max(0,-1);
    for (int kk =0; kk < so; kk++) {
        for (int jj =0; jj < sn; jj++) {
            for (int ii =0; ii < sm; ii++) {
                int eid = index(ii,jj,kk);
                for (int j = -1; j < extraps.n0[eid]; j++) {
                    ref_map *rp;

                    if( j == -1) rp=rm0+eid;
                    else rp = extraps.f0[eid] + j;

                    for(int F=0;F<3;F++){
                        double traction = fabs(rp->coll_traction_tot[F]);
                        if(traction != 1.e200){
                            if(traction > local_max.value) {
                                local_max.value = traction;
                                local_max.rank = rp->oid();
                            }
                        }
                    }

                }
            }
        }
    }

    MPI_Allreduce(&local_max, &global_max, 1, MPI_DOUBLE_INT, MPI_MAXLOC, grid->cart);
    return global_max;
}
#endif

#if defined(DEBUG)
/** We compute the maximum difference between consecutive calculations
 * of grad level-set
 */
double_int fluid_3d::max_grad_phi_diff(){
    double_int local_max(0,-1);
    double_int global_max(0,-1);
    for (int kk =0; kk < so; kk++) {
        for (int jj =0; jj < sn; jj++) {
            for (int ii =0; ii < sm; ii++) {
                int eid = index(ii,jj,kk);
                for (int j = -1; j < extraps.n0[eid]; j++) {
                    ref_map *rp;

                    if( j == -1) rp=rm0+eid;
                    else rp = extraps.f0[eid] + j;

                    for(int f=0;f<3;f++){
                        double diff = fabs(rp->grad_phi_diff[f]);
                        if(diff != 1.e200){
                            if(diff > local_max.value) {
                                local_max.value = diff;
                                local_max.rank = rp->oid();
                            }
                        }
                    }

                }
            }
        }
    }

    MPI_Allreduce(&local_max, &global_max, 1, MPI_DOUBLE_INT, MPI_MAXLOC, grid->cart);
    return global_max;
}
#endif

/** We compute the deviation of deformation gradient tensor determinant from 1
 * which would be the incompressible case.
 */
double fluid_3d::total_dev_in_detF(double &max_dev){
    double tot = 0, gtot = 0;
    const int strides[3] = {1, sm4, smn4};
    double local_max=-1, global_max=0;
    for (int kk =0; kk < so; kk++) {
        for (int jj =0; jj < sn; jj++) {
            for (int ii =0; ii < sm; ii++) {
                int eid = index(ii,jj,kk);
                int ind = G0+eid;

                ref_map &rf = rm0[eid];
                if(rf.oid() > 0){
                    double tmp_detF = 0;
                    for(int s=0;s<3;s++) {
                        // we get the absolute value of the deviation on all 6 faces
                        // because we don't want det <1 and det>1 to cancel out each other
                        // in the averaging process

                        // get the determinant of F own by this cell, 3 lower faces
                        tmp_detF += fabs(rf.detF[s]-1);
                        // now get the determinant of F own by the immediate nearby cells
                        // the 3 upper faces
                        ref_map *rp = get_refmap(ind+strides[s], rf.oid());
                        tmp_detF += fabs(rp->detF[s]-1);
                    }
                    //  average over the faces
                    tmp_detF /= 6.;
#if defined(DEBUG)
                    lJ[ind] = tmp_detF;
#endif
                    double tmp_sq = tmp_detF * tmp_detF;
                    // compute the deviation from 1 of the average determinant
                    tot += tmp_sq;
                    if(tmp_sq > local_max) local_max = tmp_sq;
                }
#if defined(DEBUG)
                else {
                    lJ[ind] = 100.;
                }
#endif
            }
        }
    }
    tot*= dx*dy*dz;

    MPI_Allreduce(&tot, &gtot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&local_max, &global_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    gtot = sqrt(gtot);
    if(global_max == -1.0) { global_max = 0;}
    max_dev = sqrt(global_max);

    return gtot;
}

/** Compute average detF */
double fluid_3d::avg_detF(){
    double tot = 0, gtot = 0;
    const int strides[3] = {1, sm4, smn4};
    int local_dof = 0, glob_dof;
    for (int kk =0; kk < so; kk++) {
        for (int jj =0; jj < sn; jj++) {
            for (int ii =0; ii < sm; ii++) {
                int eid = index(ii,jj,kk);
                int ind = G0+eid;

                ref_map &rf = rm0[eid];
                if(rf.oid() > 0){
                    double tmp_detF = 0;
                    local_dof ++;
                    for(int s=0;s<3;s++) {
                        // we get the absolute value of the deviation on all 6 faces
                        // because we don't want det <1 and det>1 to cancel out each other
                        // in the averaging process

                        // get the determinant of F own by this cell, 3 lower faces
                        tmp_detF += rf.detF[s];
                        // now get the determinant of F own by the immediate nearby cells
                        // the 3 upper faces
                        ref_map *rp = get_refmap(ind+strides[s], rf.oid());
                        tmp_detF += rp->detF[s];
                    }
                    //  average over the faces
                    tmp_detF /= 6.;
                    // compute the deviation from 1 of the average determinant
                    tot += tmp_detF;
                }
            }
        }
    }

    MPI_Allreduce(&local_dof, &glob_dof, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&tot, &gtot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    gtot = gtot/glob_dof;

    return gtot;
}

#if defined(DEBUG)
/** Compute maximum Gaussian curvature */
double_int fluid_3d::max_Gaussian_curvature(int & sign){
    double_int local_max(0,-1);
    double_int global_max(0,-1);
    double_int local_sign(0,-1);
    double_int global_sign(0,-1);
    const int strides[3] = {1, sm4, smn4};
	const double facs[3] ={dxsp, dysp, dzsp};
    for (int kk =0; kk < so; kk++) {
        for (int jj =0; jj < sn; jj++) {
            for (int ii =0; ii < sm; ii++) {
                int eid = index(ii,jj,kk);
                int ind = G0+eid;

                for (int j = 0; j < extraps.n0[eid]; j++) {
                    ref_map &rf = extraps.f0[eid][j];
                    // only compute these at the interfaces
                    if(rf.lid() == 1 || rf.lid() == max_layers+1) {
                        // compute gaussian curvature

                        int oid = rf.oid();
                        double grad_grad_phi [3][3];
                        double avg_norm [3] = {0,0,0};
                        int avg_dof = 0;

                        for(int s=0;s<3;s++){
                            // we first find the average surface normal here
                            ref_map * rpr = get_refmap(ind+strides[s], oid);
                            ref_map * rpl = get_refmap(ind-strides[s], oid);
                            if(rpr->lid() == 0){
                                //if the cell to the right (or up or back) is in solid
                                //we take into account its lower face normal info
                                avg_dof ++;
                                for(int i=0;i<3;i++) avg_norm[i] += rpr->grad_phi[s][i];
                            }
                            if(rpl->lid() == 0){
                                //if the cell to the left (or down or front) is in solid
                                //we take into account our own lower face normal info
                                avg_dof ++;
                                for(int i=0;i<3;i++) avg_norm[i] += rf.grad_phi[s][i];
                            }
                            // here we compute gradient of grad_phi, using cell centered centered FD
                            for(int i=0;i<3;i++){
                                grad_grad_phi[s][i] = facs[s]*(rpr->grad_phi[s][i] - rf.grad_phi[s][i]);
                            }
                        }
                        // we compute the average
                        if(avg_dof ==0 ) {
                             //if we didn't find anything, that's because the inside grid point has moved to be freshly outside
                             continue;
                        } else {
                            for(int i=0;i<3;i++) avg_norm[i] /= avg_dof;
                        }

                        // compute the SURFACE GRADIENT
                        // proj_tensor is (I - n x n)
                        double proj_tensor[3][3];
                        double curve_tensor[3][3];
                        double gaussian_curvature = 0;

                        for(int s=0;s<3;s++){
                            for(int i=0;i<3;i++){
                                proj_tensor[s][i] = -avg_norm[s] * avg_norm[i];
                            }
                            proj_tensor[s][s] += 1;
                        }

                        for(int s=0;s<3;s++){
                            for(int i=0;i<3;i++){
                                curve_tensor[s][i] = 0.;
                                for(int j=0;j<3;j++){
                                    // since we want gaussian curvature, which is negative trace
                                    // we just compute negative curvature tensor
                                    curve_tensor[s][i] += grad_grad_phi[s][j] * proj_tensor[j][i];
                                }
                            }
                            gaussian_curvature += curve_tensor[s][s];
                        } // Here we conclude gaussian curvature calculation
                        if(fabs(gaussian_curvature) > local_max.value){
                            local_max.value = fabs(gaussian_curvature);
                            local_max.rank = oid;

                            local_sign.value = fabs(gaussian_curvature);
                            local_sign.rank = (gaussian_curvature<0)?-1:1;
                        }
                    } // Here we conclude computation for a single first layer cell
                }// Here we conclude computation for all first layer cells
            }
        }
    }

    MPI_Allreduce(&local_max, &global_max, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);
    MPI_Allreduce(&local_sign, &global_sign, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);

    sign = global_sign.rank;
    return global_max;
}
#endif

void fluid_3d::write_all_Js(char* fn){
    bool write=false;
    for (int kk =0; kk < so; kk++) {
        for (int jj =0; jj < sn; jj++) {
            for (int ii =0; ii < sm; ii++) {
                int eid = index(ii,jj,kk);
                int ind = G0+eid;

                ref_map *rpp = get_refmap(ind, 1);
                if(rpp!=NULL && rpp-> lid() < 5){
                    write = true;
                    break;
                }
            }
        }
    }

    if(write){
        FILE *fh = p_safe_fopen(fn, "w");
        fprintf(fh, "#x y z level-set detF-1 face oid xpos ypos zpos wforce_x wforce_y wforce_z\n");
        const int strides[3] = {1, sm4, smn4};
        for (int kk =0; kk < so; kk++) {
            for (int jj =0; jj < sn; jj++) {
                for (int ii =0; ii < sm; ii++) {
                    int eid = index(ii,jj,kk);
                    int ind = G0+eid;

                    ref_map *rpp = get_refmap(ind, 1);
                    if(rpp!=NULL) {
                        ref_map &rf = *rpp;
                        if(rf.lid()<5){
                            for(int s=0;s<3;s++) {
                                // get the determinant of F own by this cell, 3 lower faces
                                ref_map *rp = get_refmap(ind-strides[s], rf.oid());
                                if(rp==NULL) printf("Rank %d getting a %d neighbor for %d layers, but didn't get any?\n", rank, ind-strides[s], rf.lid());
                                double phiv = (rf.phi(mgmt)+rp->phi(mgmt))*0.5;
                                // now get the determinant of F own by the immediate nearby cells
                                // the 3 upper faces
                                double tmp_detF = rf.detF[s];
                                if(tmp_detF != -1.e200)
                                fprintf(fh, "%d %d %d %10.8g %10.8g %d %d %g %g %g\n", ii+ai, jj+aj, kk+ak, phiv, tmp_detF-1, 2*s, rf.oid(), lx0[ii], ly0[jj], lz0[kk]);

                                rp = get_refmap(ind+strides[s], rf.oid());
                                if(rp==NULL) printf("Rank %d getting a %d neighbor for %d layers, but didn't get any?\n", rank, ind+strides[s], rf.lid());
                                phiv = (rf.phi(mgmt)+rp->phi(mgmt))*0.5;
                                // now get the determinant of F own by the immediate nearby cells
                                // the 3 upper faces
                                tmp_detF = rp->detF[s];
                                if(tmp_detF != -1.e200)
                                fprintf(fh, "%d %d %d %10.8g %10.8g %d %d %g %g %g\n", ii+ai, jj+aj, kk+ak, phiv, tmp_detF-1, 2*s, rf.oid(), lx0[ii], ly0[jj], lz0[kk]);
                            }
                        }
                    }
                }
            }
        }
        fclose(fh);
    }
}

void fluid_3d::compute_v_centroid(){

    mgmt->reset_avgs();
    // We calculate fractions of a cell based on the level set values at the center
    // The centroid is a weighted average of the positions in the solid body
    for (int kk =0; kk < so; kk++) {
        for (int jj =0; jj < sn; jj++) {
            for (int ii =0; ii < sm; ii++) {
                int eid = index(ii,jj,kk);
                ref_map &rf = rm0[eid];
                field &ff = u0[eid];
                if(rf.oid() > 0){
                    double phiv = rf.phi(mgmt);
                    double sfrac = mgmt->heaviside(phiv);
                    int s = rf.oid()-1;
                    mgmt->avg_velN[s] += sfrac;
                    mgmt->avg_velx[s] += ff.vel[0]*sfrac;
                    mgmt->avg_vely[s] += ff.vel[1]*sfrac;
                    mgmt->avg_velz[s] += ff.vel[2]*sfrac;

                    mgmt->avg_N[s] += sfrac;
                    mgmt->avg_x[s] += lx0[ii]*sfrac;
                    mgmt->avg_y[s] += ly0[jj]*sfrac;
                    mgmt->avg_z[s] += lz0[kk]*sfrac;
                }
            }
        }
    }

    for(int i=0;i<mgmt->n_obj;i++) {
        double glov, centx, centy, centz;
        MPI_Allreduce(&(mgmt->avg_velN[i]), &glov, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&(mgmt->avg_velx[i]), &centx, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&(mgmt->avg_vely[i]), &centy, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&(mgmt->avg_velz[i]), &centz, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        mgmt->avg_velx[i] = centx;
        mgmt->avg_vely[i] = centy;
        mgmt->avg_velz[i] = centz;
        mgmt->avg_velN[i] = glov;
    }

    for(int i=0;i<mgmt->n_obj;i++) {
        double glov, centx, centy, centz;
        MPI_Allreduce(&(mgmt->avg_N[i]), &glov, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&(mgmt->avg_x[i]), &centx, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&(mgmt->avg_y[i]), &centy, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&(mgmt->avg_z[i]), &centz, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        mgmt->avg_x[i] = centx;
        mgmt->avg_y[i] = centy;
        mgmt->avg_z[i] = centz;
        mgmt->avg_N[i] = glov;
    }

    mgmt->compute_avgs();
}

void fluid_3d::compute_centroid(double &centx, double &centy, double &centz, double &vol){
    double locx=0, locy=0, locz=0;
    // We calculate fractions of a cell based on the level set values at the center
    // The centroid is a weighted average of the positions in the solid body
    double locv=0, glov=0;
    for (int kk =0; kk < so; kk++) {
        for (int jj =0; jj < sn; jj++) {
            for (int ii =0; ii < sm; ii++) {
                int eid = index(ii,jj,kk);
                ref_map &rf = rm0[eid];
                if(rf.oid() > 0){
                    double phiv = rf.phi(mgmt);
                    double sfrac = mgmt->heaviside(phiv);
                    locv += sfrac;
                    locx += lx0[ii]*sfrac;
                    locy += ly0[jj]*sfrac;
                    locz += lz0[kk]*sfrac;
                }
            }
        }
    }

    MPI_Allreduce(&locv, &glov, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&locx, &centx, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&locy, &centy, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&locz, &centz, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    centx/= glov;
    centy/= glov;
    centz/= glov;
    vol = glov*dx*dy*dz;
}

/*
double fluid_3d::check_steady_solution() {
    double local_max=-1, global_max;
	for (int k =0; k < so; k++) {
		for (int j =0; j < sn; j++) {
			for (int i = 0; i < sm; i++) {
				int eid = index(i, j, k);
				field &f = u0[eid];
                vel_field &vf = sim_u0[eid];
				for(int dir=0; dir<3; dir++) {
                    double tmp = fabs(vf.vel[dir] - f.vel[dir]);
                    if(tmp>local_max) local_max = tmp;
                }
            }
        }
    }

    local_max /= (dt/mgmt->fm.mu); // assuming rho=1, U=1, L=1, Renum = 1/visc
    MPI_Allreduce(&local_max, &global_max, 1, MPI_DOUBLE, MPI_MAX, grid->cart);
    return global_max;
}
*/
