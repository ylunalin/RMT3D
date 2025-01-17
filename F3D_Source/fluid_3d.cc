#include "fluid_3d.hh"

/**
 * Constructor for fluid_3d class, representing incompressible
 * Navier-Stokes simulation.
 * \param[in] gm the geometry class object reference.
 * \param[in] mgmt_ pointer to the sim_manager object
 *            manages simulation parameters and the simulation itself.
 *            e.g. Reynold's number, shear modulus of solids, initial vels, etc.
 * \param[in] spars_ pointer to the simulation parameter object that handles IO
 *            of user specified parameters. */
fluid_3d::fluid_3d(geometry *gm, sim_manager *mgmt_, sim_params *spars_) :
    mgmt(mgmt_), spars(spars_),
	impl(spars->impl),godunov(spars->godunov),
#if defined(TANG)
    tang(true),
#else
    tang(false),
#endif
    trace(false),
	rank(gm->rank),m(gm->m),n(gm->n),o(gm->o),
    sm(gm->sm),sn(gm->sn),so(gm->so),
    sm4(sm+4),sn4(sn+4),so4(so+4),
	smn4(sm4*sn4),G0(index(2,2,2)),
    ai(gm->ai),aj(gm->aj),ak(gm->ak),
    bi(gm->bi),bj(gm->bj),bk(gm->bk),
    extp_err(0), dd(0),ds(m/2),nt(0),
    max_layers(spars->nlayers),
    max_extrap_rs(spars->max_extrap_rs),
	dx(spars->dx),dy(spars->dy),dz(spars->dz),
	dxsp(1.0/dx),dysp(1.0/dy),dzsp(1.0/dz),
	dt(spars->dt),
    time(0),fem_grid(*gm),temp_extraps(*gm),extraps(*gm),
    // mg_prob definitions
#if defined(VAR_DEN)
    visc(gm),
    pres(&fem_grid, sm4, sn4, so4),
    mac(gm, sm4, sn4, so4),
#else
    visc(gm),
    pres(&fem_grid),
    mac(gm),
#endif
    tr(NULL),grid(gm),
	u_mem(new field[smn4*so4]),u0(u_mem+2*(1+sm4*(1+sn4))),
    // New density and level-set field
#if defined(VAR_DEN)
    lrho(new double[smn4*so4]),
#endif

#if defined(DEBUG)
    lJ(new double[smn4*so4]),
#endif

	rm_mem(new ref_map[smn4*so4]),rm0(rm_mem+2*(1+sm4*(1+sn4))),
	gx(new double[m]), gy(new double [n]), gz(new double[o]),
	lx(new double[sm4]),ly(new double[sn4]),lz(new double[so4]),
    lx0(lx+2),ly0(ly+2),lz0(lz+2),out(NULL),
    watch("F3d", 5, "Extrapolation", "Pressure Poisson", "Communication", "Stresses", "MAC"),
	expper(*gm, mgmt->weight_fac, spars),
	osm(NULL),osn(NULL),oso(NULL),max_sizes(NULL),
	cu_m(NULL),cu_n(NULL),cu_o(NULL) {

	// Set up the polygonalizing tables
	*ptri_poly=const_cast<char*>(tri_poly);
	for(int k=1;k<=256;k++) ptri_poly[k]=ptri_poly[k-1]+3*n_poly[k-1];

	// various setup functions
	setup_faces();             // boundary conditions
	setup_communication();     // MPI communication objects (new primary exp)
	setup_cart();
	setup_output_dimensions(); // gather dimensions on master processor needed for contour output
	if(spars->display) display_stats();
	set_omp(spars->omp_num_thr);

    printf("#### I'm rank %d, my own box (%d %d %d) to (%d %d %d)\n", rank, ai, aj, ak, bi, bj, bk);
}

/**
 * Destructor
 */
fluid_3d::~fluid_3d() {

	if (trace) delete tr;
	if (impl) cleanup_impl_visc();
	if(godunov) cleanup_mac();
	cleanup_fem();
	cleanup_communication();
	if(osm!=NULL) {
		delete [] cu_m;
		delete [] osm;
	}
	if (out != NULL) delete[] out;
	delete[] lz;
	delete[] ly;
	delete[] lx;
	delete[] gz;
	delete[] gy;
	delete[] gx;
	delete[] rm_mem;
	delete[] u_mem;
#if defined(VAR_DEN)
    delete[] lrho;
#endif

#if defined(DEBUG)
    delete[] lJ;
#endif
}

const int fluid_3d::domain_face[6]= {12,14,10,16,4,22};

const bool fluid_3d::pres_update = true;

void fluid_3d::alloc_tracers(){
	tr = new tracers(*this);
	trace = true;

	/*
	puts("g\n");
	for (tracer t:tr->tr) {
		if (t.p_crd[0] > 200) {
			printf("here (predict start)!\n");
			MPI_Abort(MPI_COMM_WORLD,0);
		}
	}
	*/
}

/** Set up physical coordinates in local domain and global domain. */
void fluid_3d::setup_cart() {
	for (int i = 0; i < m; i++) gx[i] = mgmt->ax + (i + 0.5)*dx;
	for (int j = 0; j < n; j++) gy[j] = mgmt->ay + (j + 0.5)*dy;
	for (int k = 0; k < o; k++) gz[k] = mgmt->az + (k + 0.5)*dz;

	for (int i = 0; i < sm4; i++) lx[i] = mgmt->ax + (ai + i - 2 + 0.5)*dx;
	for (int j = 0; j < sn4; j++) ly[j] = mgmt->ay + (aj + j - 2 + 0.5)*dy;
	for (int k = 0; k < so4; k++) lz[k] = mgmt->az + (ak + k - 2 + 0.5)*dz;
}

/**
 * Set the ghost nodes on one of the faces of the local grid.
 * \param[in] tag the integer label of a neighbor
 */
void fluid_3d::set_bc_cc (bool verbose) {
    if(verbose) printf("Rank %d. Set cell-centered boundary conditions.\n", rank);
    bool v = false;
	for (int s = 0; s < 6; s++){
        if (neigh[domain_face[s]] == -1) faces[s].set_ghost_cc(v);
    }
	communicate<1>();
}

/** Set the cell face velocities at the boundary */
void fluid_3d::set_bc_cf (bool verbose){
    if(verbose) printf("Rank %d. Set cell face boundary conditions.\n", rank);
    bool v = false;
	for (int s = 0; s < 6; s++){
		if (neigh[domain_face[s]] == -1)
        {
			faces[s].set_ghost_cf(v);
        }
    }
	communicate<4>();
}

/**
 * A helper function to generate x, y, z tags (a number
 * between 0 and 2 inclusive) for 27 cubes
 * and determine if a certain cube is neighboring a
 * face of the center cube
 * \param[in] tag the integer number of neighbor, 0 is the neighbor sharing ai, aj, ak
 * \param[in, out] (ti, tj, tk) tags in x, y, z direction
 */
int fluid_3d::nface(int tag, int &ti, int &tj, int &tk){
	tk = tag / 9;       // z index just int divide
	ti = tag % 3;       // x index just mod
	tj = (tag / 3) % 3; // y index do both
	// return face number, if it's not a face return -1
	// Compute if the neighbor is an orthogonal neighbor, if it is, \sum_l |t_l-1| = 1
	int num_ones = 0;
	if(ti==1) num_ones++;
	if(tj==1) num_ones++;
	if(tk==1) num_ones++;
	if (num_ones == 2)
		return (ti==0)?0:((ti==2)?1:((tj==0)?2:((tj==2)?3:((tk==0)?4:5))));
	else return -1;
}

/** ##################### PROBLEM SETUP ####################*/
/** Initializes all fields and reset simulation time. */
int fluid_3d::initialize() {
    int init_err=1;
    if(spars->chk_num>0) {
	    init_err=initialize_from_chk_point(spars->chk_dirname);
    }

    if(init_err>0){
        if(rank==0) printf("# Initializing from scratch, setting chk_num to 0.\n");
        spars->chk_num = 0;
        time=0;
        init_fluid();
        init_pres();
        init_refmap();
        init_extrapolate(false);
#if defined(VAR_DEN)
        update_rho();
#endif

        // NOTE: No need to communicate here
        // initialization of field values are communicated in individual function
        // initial extrapolation gets communicated after each layer is done
    }

    fill_boundary_cc(0);

    if((spars->dump_code & 1) && spars->ntracers>0){
        alloc_tracers();
    }

    // dt related constant tables
	setup_stencils();          // stencil class
	compute_grad_basis();

    // since now we are using density in linear solves
    // we set up these solvers here:
	setup_fem();               // primary multigrid solver for projection method
	if(godunov) setup_mac();   // mac projection if we are using Godunov scheme
	if(impl) {
        setup_impl_visc();// implicit viscosity solver if applicable
    }
    return init_err;
}

void fluid_3d::init_iter(int init_err){

    // We do some iteration to get a good estimate of pressure
    if(init_err>0) {
        double tol = 1e-4;
        int num_iter = spars->num_iters;
        if(rank==0) printf("# Doing %d initial iterations, or until div(u) < %e\n", num_iter, tol);
		/*
	puts("i\n");
	for (tracer t:tr->tr) {
		if (t.p_crd[0] > 200) {
			printf("here (i)!\n");
			MPI_Abort(MPI_COMM_WORLD,0);
		}
	}
		*/
        double divu_l2;
        for (int i=0;i<num_iter;i++) {
            time=0;
            // We need this fill boundary because recalling all the initializing function clears very thing
            fill_boundary_cc(0);
            step_forward(spars->debug_flag);

            double_int di1, di2;
            divu_l2 = div_u(di1,di2);
            if(rank==0) printf("#-->%d %g\n", i, divu_l2);

            // After one iteration, we only keep pressure, but reset everything else
            init_fluid();
            init_refmap();
            extraps.reset();
            init_extrapolate(false);
#if defined (VAR_DEN)
            update_rho();
#endif
            if(divu_l2<tol) break;
        }
        nt=0;
	time=0.0;
    }
    // Need this to get F, stress, elastic energy, for initial macro or diagnostic output
    compute_stress(0);

	/*
	puts("j\n");
	for (tracer t:tr->tr) {
		if (t.p_crd[0] > 200) {
			printf("here (j)!\n");
			MPI_Abort(MPI_COMM_WORLD,0);
		}
	}
	*/
}

/**
 * Use functions specified in the sim_manager object to initialize fluid fields.
 */
void fluid_3d::init_fluid() {
	for (int k = 0; k < so4; k++) {
		double zz = lz[k];
		for (int j = 0; j < sn4; j++) {
			double yy = ly[j];
			for (int i = 0; i < sm4; i++) {
				double xx = lx[i];

                int ind = index(i,j,k);
#if defined (VAR_DEN)
                lrho[ind] = mgmt->fm.rho;
#endif
				field *f = u_mem + ind;
                f->reset();
				double p, tmp[3];
				// if we are initializing using an exact solution
				if(mgmt->has_exact){
					mgmt->exact(xx,yy,zz,time,dx,tmp[0],tmp[1],tmp[2],p);
				}
				else{
					mgmt->velocity(xx,yy,zz,tmp[0],tmp[1],tmp[2]);
				}
				f->set_vels(tmp);
			}
		}
	}
	communicate<1>();
    mgmt->obligatory_cfl_recheck();
    dt = mgmt->dt_reg;
}

/**
 * Use functions specified in the sim_manager object initialize the pressure field.
 */
void fluid_3d::init_pres() {

	// at each node...
	double hdx=dx*0.5, hdy=dy*0.5, hdz=dz*0.5;
	for (int k = 2; k < so4-1; k++) {
		double zz = lz[k] - hdz;
		for (int j = 2; j < sn4-1; j++) {
			double yy = ly[j] - hdy;
			for (int i = 2; i < sm4-1; i++) {
				double xx = lx[i] - hdx;
				field *f = u_mem + index(i, j, k);
                f->q = 0;
				mgmt->pressure(xx,yy,zz,f->p);
			}
		}
	}
	communicate<2>();
}

/**
 * Initializes the reference map.
*/
void fluid_3d::init_refmap(){

	ref_map *rm;
    // Set every reference map to be a void point to begin with
    for (int i=0;i<smn4*so4;i++) {
        rm_mem[i].reset();
    }

	for (int k = 2; k < so4-2; k++) {
		double zz=lz[k];
		for (int j = 2; j < sn4-2; j++) {
			double yy = ly[j];
			for (int i = 2; i < sm4-2; i++) {
				double xx = lx[i];
                int ind = index(i,j,k);
#if defined(DEBUG)
                lJ[ind] = 1.0;
#endif
				rm = rm_mem + ind;
				double tmpx[3] = {xx,yy,zz};

				// check each object to see if this point is in it
				bool found = false;
				for (int o = 0; (o < mgmt->n_obj) ; o++) {

					mgmt->objs[o]->rm0(tmpx,rm->x); // get initial refmap here
					mgmt->objs[o]->rm0(tmpx,rm->xpred); // get initial refmap here
					// if that gives us a negative phi...
					if (mgmt->phi(o,rm->x) < 0) {

						// mark that it belongs to this object
                        if(!found){
                            rm->c = o+1;
                            found = true;
#if defined(VAR_DEN)
                            sl_mat & tmp_sm = mgmt->sm_array[o];
                            lrho[ind] = tmp_sm.rho;
#endif
                        } else {
                            printf("fluid_3d::init_refmap: overlapping primary points!\n");
                            p_fatal_error("Initialization error!", 1);
                        }
                    }
				}
			}
		}
	}

	// communicate object id
	communicate<8>();
	// communicate reference map x/xpred
	communicate<16|32>();
}

/**
 * Set the pressure field at current time as specified by the sim_manager class..
 */
void fluid_3d::set_pres() {

	if(mgmt->has_exact) {
		for (int k = 0; k < so4; k++) {
			double zz = lz[k];
			double u_,v_,w_, yy, xx;
			for (int j = 0; j < sn4; j++) {
				yy = ly[j];
				for (int i = 0; i < sm4; i++) {
					xx = lx[i];
					field *f = u_mem + index(i,j,k);
					mgmt->exact(xx,yy,zz,time,dx,u_,v_,w_,f->p);
				}
			}
		}
	}
}

/** ##################### TIMESTEPPING ######################*/

/** An alternative way to step forward. Note: custom timestep is not fully implemented.
 *  dt is still integrated in stencil class and finite element projection solve.
 *  \param[in] debug is the timestep we want to turn debug functions on
 */
int fluid_3d::step_forward(int debug){

	// 1) compute stress at t_n
	// 2) Godunov scheme to handle advective term
	//    - calculate derivatives
	//    a) normal derivs
	//    b) tang stability
	//    - extrapolate v and ref_map to edges
	//
	// 3) compute ref map at t_{n+1}
	// 4) average to get ref map at t_{n+1/2}
	// 5) compute stress at t_{n+1/2}
	//    - note this is solid stress at t_{n+1/2} and
	//    - fluid stress at t_n. could in principle
	//    - predict velocities at t_{n+1/2} to get
	//    - representation of global stress at t_{n+1/2}
	// 6) calculate velocities at t_star
	// 7) projection step to get p at t_{n+1/2}
	// 8) use p to get velocities at t_{n+1}
	// 9) update ref_map to already-computed t_{n+1}

	nt++;
    bool basic_v  = (debug<=-1 && rank==0) ;
    bool detail_v = (debug<=-2 && rank==0) ;
    //bool verbose_max = (debug==-3) ;
    int ERR=0;

	if (out != NULL && rank == 0) out[0] = '\0';

	if (trace) {
//	if (rank==0) printf("FIRST TRACER\n");
		tr->predict();
	}
	if(basic_v) {
        int success=0;
        if(rank==0) success = system("echo \"Clock: `date`\"");
        printf("Rank %d. Print time %d. Starting timestep %d. Current sim time %g\n", rank, success, nt, time);
    }

    // Need cell-centered values at boundaries (interior or physical)
    // to compute values on cell faces
	compute_stress(detail_v);

    // Need cell-centered values at boundaries (interior or physical)
	normal_derivatives_velocity(dt, detail_v);

    // Need cell-centered values at boundaries (interior or physical)
	normal_derivatives_refmap(dt, detail_v);

#if defined(TANG)
	// tangential stability calculation
    // Upwind to determine the face velocities (ref map)
    // in tang stability calculation.
    // Since upwinding start from the lower, inner ghost layer,
    // no need to set ghost layers after upwinding if periodic.
    godunov_upwinding_set(tang, detail_v);

    // Compute the tangential derivatives terms in adv term for
    // both fluid and reference map in the same routine.
    compute_tang_derivatives(detail_v);
#endif

	compute_half_time_edge_velocities(dt, detail_v);

	compute_half_time_edge_refmap(dt, detail_v);

    // Need cell face at boundaries for upwinding (interior or physical)
    // but it handles its own communication, but when should we set physical boundary?
    // NOTE: we decide to upwind then set bc, but what values are used in upwinding at the boundary???
    time += 0.5*dt;

    fill_boundary_cf(detail_v);
	godunov_upwinding_set(false, detail_v);

    // Set boundary values at faces
    fill_boundary_cf(detail_v);

    // MAC project the edge velocities
    // mac_solve() takes care of bc on q field, the ad hoc pressure
    //check_nans(2);
    mac_solve(dt, basic_v);
    //check_nans(3);

    // In subtracting ad hoc pressure
    // we might have subtracted some bullshit values at the physical boundary
    fill_boundary_cf(detail_v);

    if(basic_v) printf("Rank %d. Godunov foo is done\n", rank);

	ERR = update_reference_map(dt, detail_v);
	if (ERR) {
		revert_reference_map_update();
		return REFMAP_UPDATE_ERR;
	}

    // at this point primary reference map and extrapolated ref map are up to date in ghost regions
    // if they are interior boundaries
    // but we still need to set physical boundaries

    // this is unnecessary, we think?
    //fill_boundary_cc(detail_v);

    compute_stress(detail_v);

	compute_ustar(dt, detail_v);

    // since ustar is computed, we need to reset interior and physical boundaries
    fill_boundary_cc(detail_v);
    // Also communicate dvel
    communicate<2048>();

	time += 0.5*dt;

    // need to assume pressure bc are handled by the fill source
	ERR = pp_solve(basic_v);
	if (ERR) {
        double_int gmx_ext, gmn_ext;
        div_u(gmx_ext, gmn_ext);
        if(rank==0) {
            printf("Maximum divergence is %g on rank %d \n", gmx_ext.value, gmx_ext.rank);
            printf("Minimum divergence is %g on rank %d \n", gmn_ext.value, gmn_ext.rank);
        }
		return PROJECTION_ERR;
	}

	update_reference_map_full(detail_v);

	redraw_boundary(detail_v);

    fill_boundary_cc(detail_v);
    // Also communicate dvel
    communicate<2048>();

	// Mutigrid info pass into output string
	if (out != NULL && rank == 0) {
		if (godunov) sprintf(out,"%s%s",out,*mac.out);
		if (impl) sprintf(out,"%s%s%s%s",out,visc.out[0],visc.out[1],visc.out[2]);
		sprintf(out,"%s%s",out,*pres.out);
	}

	if (trace) {
//	if (rank==0) printf("SECOND TRACER\n");
		tr->advect();
	}

	if(basic_v) {puts("");}

    return NORMAL_EXIT;
}

/** Pre-compute monotonocity-limited normal derivatives and extrapolated face vels.
 * \param[in] cdt the custom time step to use (full step).
 */
void fluid_3d::normal_derivatives_velocity(const double cdt, bool verbose){
	if(verbose) printf("Rank %d. Normal velocity derivatives.\n", rank);

	const int strides[3]={1,sm4,smn4};
	const double dhsp[3] = {dxsp, dysp, dzsp};
#if defined(TANG)
	const double dh2[3] ={dx/2, dy/2, dz/2};
	const double dt2 = cdt/2.0;
#endif
	// pre compute the monotonocity limited derivatives
	// and extrapolate intermediate tangential velocites
	for (int kk =0; kk < so; kk++) {
		for (int jj =0; jj < sn; jj++) {
			for (int ii = 0; ii < sm; ii++) {
				field *fp = u0 + index(ii,jj,kk);
				field &f = *fp;
				// for dx, dy, dz
				for (int i=0;i<3;i++){
					int strd=strides[i];
					// for u, v, w
					for(int j=0;j<3;j++){
						// compute monotonocity limited derivatives
						double tmp_del_r = dhsp[i]* mono_diff(fp[-2*strd].get(j), fp[-strd].get(j), f.get(j), fp[strd].get(j), fp[2*strd].get(j));
						// store the computed derivatives
						f.set_mono_derivs(i, j, tmp_del_r);
					}
				}
#if defined(TANG)
				// do tangential stability calculation if required
				// extrapolate tangetial velocities to all 6 faces
				f.inter_extrap(dh2, dt2);
#endif
			}
		}
	}
}

/** Pre-compute monotonocity-limited normal derivatives and extrapolated face reference maps.
 */
void fluid_3d::normal_derivatives_refmap(const double cdt, bool verbose){
	if(verbose) printf("Rank %d. Normal refmap derivatives.\n",rank);
	const int strides[3]={1,sm4,smn4};
	const double dhsp[3] = {dxsp, dysp, dzsp};
#if defined(TANG)
	const double dh2[3] ={dx/2, dy/2, dz/2};
	const double dt2 = cdt/2.0;
#endif
	for (int k =0; k < so; k++) {
		for (int j =0; j < sn; j++) {
			for (int i = 0; i < sm; i++) {
				int eid = index(i,j,k);
				int ind = eid + G0;

				ref_map *sfp = rm0 + eid;
#if defined(TANG)
				field &f = u0[eid];
#endif
				/////////////////////////////
				//
				// |
				// | here we calculate the the normal derivs and do
				// | the extrapolation at every ref_map inside the object
				// v
				//
				////////////////////////////
				int obj_id = sfp->oid();
				if (obj_id > 0) {

					// for dx, dy, dz
					for (int ii=0;ii<3;ii++){
                        int m1pos [3]= {i+2, j+2, k+2};
                        int m2pos [3]= {i+2, j+2, k+2};
                        int p1pos [3]= {i+2, j+2, k+2};
                        int p2pos [3]= {i+2, j+2, k+2};

                        m1pos[ii] -= 1;
                        m2pos[ii] -= 2;
                        p1pos[ii] += 1;
                        m2pos[ii] += 2;
						int strd=strides[ii];

						// for u, v, w
						for(int jj=0;jj<3;jj++){

							// compute monotonocity limited derivatives
							ref_map *m2,*m1,*p1,*p2;
							m2 = get_refmap(ind-2*strd,obj_id, m2pos);
							m1 = get_refmap(ind-strd,obj_id, m1pos);
							p1 = get_refmap(ind+strd,obj_id, p1pos);
							p2 = get_refmap(ind+2*strd,obj_id, p2pos);
							if (m2 == NULL || m1 == NULL || p1 == NULL || p2 == NULL) {
								printf("fluid_3d:: (%d %d %d) normal_derivs_refmap(): on rank %d in direction %d at mono diff\n",i, j, k, grid->rank, ii);
								printf("fluid_3d:: my level-set value here is %.10g in object %d\n", sfp->phi(mgmt), sfp->oid());
                                if(m2==NULL) printf("Missing %d neighbor\n", -2*strd);
                                if(m1==NULL) printf("Missing %d neighbor\n", -strd);
                                if(p1==NULL) printf("Missing %d neighbor\n", strd);
                                if(p2==NULL) printf("Missing %d neighbor\n", 2*strd);
                                p_fatal_error("normal_derivatives_refmap: missing necessary cells to carry out finite differences.\n", 1);
							}
							double tmp_del_r = dhsp[ii]* mono_diff(
								m2->x[jj],
								m1->x[jj],
								sfp->x[jj],
								p1->x[jj],
								p2->x[jj]);

							// store the computed derivatives
							sfp->set_mono_derivs(ii, jj, tmp_del_r);
						}
					}

#if defined(TANG)
					// here tangential stability calculation is required
					// extrapolate tangetial velocities to all 6 faces
					sfp->inter_extrap(f, dh2, dt2);
#endif
				}

				/////////////////////////////
				//
				// |
				// | and here we look through all the extraps at this spot, and
				// | if any of them have a neighbor that's inside of an object,
				// | we calculate the normal derivs and do the interp_extract there too
				// v
				//
				////////////////////////////
				for (int u = 0; u < extraps.n0[eid]; u++) {

					sfp = extraps.f0[eid]+u;
					obj_id = sfp->oid();
					if(sfp->abut_primary(eid,rm0,strides)){
						// for dx, dy, dz
						for (int ii=0;ii<3;ii++){
                            int m1pos [3]= {i+2, j+2, k+2};
                            int m2pos [3]= {i+2, j+2, k+2};
                            int p1pos [3]= {i+2, j+2, k+2};
                            int p2pos [3]= {i+2, j+2, k+2};

                            m1pos[ii] -= 1;
                            m2pos[ii] -= 2;
                            p1pos[ii] += 1;
                            m2pos[ii] += 2;
							int strd=strides[ii];

							// for u, v, w
							for(int jj=0;jj<3;jj++){

								// compute monotonocity limited derivatives
								ref_map *m2,*m1,*p1,*p2;
								m2 = get_refmap(ind-2*strd,obj_id, m2pos);
								m1 = get_refmap(ind-strd,obj_id, m1pos);
								p1 = get_refmap(ind+strd,obj_id, p1pos);
								p2 = get_refmap(ind+2*strd,obj_id, p2pos);
								if (m2 == NULL || m1 == NULL || p1 == NULL || p2 == NULL) {
									printf("fluid_3d:: normal_derivs_refmap(): on rank %d at mono diff in extrap\n"
											"index %d (%d %d %d) global index (%d %d %d) layer %d\n",
											grid->rank, ind, i, j, k, ai+i, aj+j, ak+k, sfp->lid());
									sfp->abut_primary(eid,rm0,mgmt,strides);
									if(m2==NULL) {
                                        printf("Missing %d neighbor\n", -2*strd);
                                        int eidtmp = eid-2*strd;
                                        printf("Number of extrapolation at missing neighbor is %d\n", extraps.n0[eidtmp]);
                                        for (int tmpi = 0 ; tmpi<extraps.n0[eidtmp]; tmpi ++){
                                            printf("extrapolation in layer %d\n", extraps.f0[eidtmp][tmpi].lid());
                                        }
                                    }

                                    if(m1==NULL) {
                                        printf("Missing %d neighbor\n", -strd);

                                        int eidtmp = eid-strd;
                                        printf("Number of extrapolation at missing neighbor is %d\n", extraps.n0[eidtmp]);
                                        for (int tmpi = 0 ; tmpi<extraps.n0[eidtmp]; tmpi ++){
                                            printf("extrapolation in layer %d\n", extraps.f0[eidtmp][tmpi].lid());
                                        }

                                    }

									if(p1==NULL) {
                                        printf("Missing %d neighbor\n", strd);

                                        int eidtmp = eid+strd;
                                        printf("Number of extrapolation at missing neighbor is %d\n", extraps.n0[eidtmp]);
                                        for (int tmpi = 0 ; tmpi<extraps.n0[eidtmp]; tmpi ++){
                                            printf("extrapolation in layer %d\n", extraps.f0[eidtmp][tmpi].lid());
                                        }
                                    }

									if(p2==NULL) {
                                        printf("Missing %d neighbor\n", 2*strd);

                                        int eidtmp = eid+2*strd;
                                        printf("Number of extrapolation at missing neighbor is %d\n", extraps.n0[eidtmp]);
                                        for (int tmpi = 0 ; tmpi<extraps.n0[eidtmp]; tmpi ++){
                                            printf("extrapolation in layer %d\n", extraps.f0[eidtmp][tmpi].lid());
                                        }
                                    }

                                    p_fatal_error("normal_derivatives_refmap: extrapolated rm missing necessary cells to carry out finite differences.\n", 1);
								}

								double tmp_del_r = dhsp[ii]* mono_diff(
									m2->x[jj],
									m1->x[jj],
									sfp->x[jj],
									p1->x[jj],
									p2->x[jj]);

								// store the computed derivatives
								sfp->set_mono_derivs(ii, jj, tmp_del_r);
							}
						}

#if defined (TANG)
						// here tangential stability calculation is required
						// extrapolate tangetial velocities to all 6 faces
						sfp->inter_extrap(f, dh2, dt2);
#endif
					}
				}
			}
		}
	}
}

/** Upwinding to determine face velocities.
 *\param[in] tang_stab a switch to indicate if tangential velocities are considered.
 */
void fluid_3d::godunov_upwinding_set(bool tang_stab, bool verbose){
    if(verbose) printf("Rank %d. Godunov type upwinding. Tangential stability = %d.\n",rank, tang_stab);

    // fill ghost nodes cell face values, and reference map values
	communicate<4|64>();
	communicate_extrap_fx();

	const int strides[3]={1,sm4,smn4};

	// starting from the inner ghost layer
	for (int kk =-1; kk < so; kk++) {
		for (int jj =-1; jj < sn; jj++) {
			for (int ii =-1; ii < sm; ii++) {
				int eid = index(ii,jj,kk);
				int ind = eid+G0;
				field *fp = u0+eid;
				field &f=*fp;
				ref_map *sfp = rm0+eid;

				/////////////
				//
				// |
				// |
				// | so, first we upwind any objects in the primary grid
				// |
				// v
				//
				////////////////////////
				if (sfp->oid() > 0){
                    for(int dir=0;dir<3;dir++){

                        // compare my right, top, back faces
                        // with neightbors' left, bottom, front faces
                        field &ngbr = fp[strides[dir]];
                        // grab the right neighbor
                        ref_map *tmpn = get_refmap(ind+strides[dir],sfp->oid());
                        //if (tmpn == NULL) printf("on rank %d in upwinding set\n",grid->rank);
                        ref_map &map_ngbr = *tmpn;

                        int ff = 2*dir+1, nf = ff^1;
                        int alt1 = (dir+1)%3, alt2 = (dir+2)%3;
#if defined(TANG)
                        if(tang_stab) {
                            godunov_set_tang_refmap(f.fvel[ff][dir],f.fvel[ff][alt1],f.fvel[ff][alt2],
                                ngbr.fvel[nf][dir],ngbr.fvel[nf][alt1],ngbr.fvel[nf][alt2],
                                sfp->fx[ff][dir],sfp->fx[ff][alt1],sfp->fx[ff][alt2],map_ngbr.fx[nf][dir],
                                map_ngbr.fx[nf][alt1], map_ngbr.fx[nf][alt2]);
                        } else {
                            godunov_set_refmap(f.fvel[ff][dir],f.fvel[ff][alt1],f.fvel[ff][alt2],
                                ngbr.fvel[nf][dir],ngbr.fvel[nf][alt1], ngbr.fvel[nf][alt2],
                                sfp->fx[ff][dir],sfp->fx[ff][alt1],sfp->fx[ff][alt2],map_ngbr.fx[nf][dir],
                                map_ngbr.fx[nf][alt1],map_ngbr.fx[nf][alt2]);
                        }
#else
                        godunov_set_refmap(f.fvel[ff][dir],f.fvel[ff][alt1],f.fvel[ff][alt2],
                            ngbr.fvel[nf][dir],ngbr.fvel[nf][alt1], ngbr.fvel[nf][alt2],
                            sfp->fx[ff][dir],sfp->fx[ff][alt1],sfp->fx[ff][alt2],map_ngbr.fx[nf][dir],
                            map_ngbr.fx[nf][alt1],map_ngbr.fx[nf][alt2]);
#endif
                    }
                }

				// now, also do this for anything that's been extrapolated here and is
				// also in the first layer outside of the solid, if it has a neighbor just
				// to the right

				for (int j = 0; j < extraps.n0[eid]; j++) {
					bool found = false;
					for (int d = 0; d < 3 && !found; d++) {
						if (extraps.f0[eid][j].neigh_inside(strides[d]+eid,rm0)) {

							sfp = extraps.f0[eid]+j;
							for(int dir=0;dir<3;dir++){

								// compare my right, top, back faces
								// with neightbors' left, bottom, front faces
								field &ngbr = fp[strides[dir]];
								// grab the right neighbor
								ref_map &map_ngbr = *get_refmap(ind+strides[dir], sfp->oid());

								int ff = 2*dir+1, nf = ff^1;
								int alt1 = (dir+1)%3, alt2 = (dir+2)%3;
#if defined(TANG)
								if(tang_stab) {
								godunov_set_tang_refmap(f.fvel[ff][dir],f.fvel[ff][alt1],f.fvel[ff][alt2],
									ngbr.fvel[nf][dir],ngbr.fvel[nf][alt1],ngbr.fvel[nf][alt2],
									sfp->fx[ff][dir],sfp->fx[ff][alt1],sfp->fx[ff][alt2],
									map_ngbr.fx[nf][dir],map_ngbr.fx[nf][alt1], map_ngbr.fx[nf][alt2]);
								} else {
								godunov_set_refmap(f.fvel[ff][dir],f.fvel[ff][alt1],f.fvel[ff][alt2],
									ngbr.fvel[nf][dir],ngbr.fvel[nf][alt1], ngbr.fvel[nf][alt2],
									sfp->fx[ff][dir],sfp->fx[ff][alt1],sfp->fx[ff][alt2],
									map_ngbr.fx[nf][dir],map_ngbr.fx[nf][alt1],map_ngbr.fx[nf][alt2]);
								}
#else
                            godunov_set_refmap(f.fvel[ff][dir],f.fvel[ff][alt1],f.fvel[ff][alt2],
                            ngbr.fvel[nf][dir],ngbr.fvel[nf][alt1], ngbr.fvel[nf][alt2],
                            sfp->fx[ff][dir],sfp->fx[ff][alt1],sfp->fx[ff][alt2],
                            map_ngbr.fx[nf][dir],map_ngbr.fx[nf][alt1],map_ngbr.fx[nf][alt2]);

#endif
							}
							found = true;
						}
					}
				}

				for (int dir = 0; dir < 3; dir++) {
					field &ngbr = fp[strides[dir]];
					int ff = 2*dir+1, nf = ff^1;
					int alt1 = (dir+1)%3, alt2 = (dir+2)%3;
					// after handling solids, do fluid!
#if defined(TANG)
					if(tang_stab) {
						godunov_set_tang_velocity(f.fvel[ff][dir],f.fvel[ff][alt1],
							f.fvel[ff][alt2],ngbr.fvel[nf][dir],
							ngbr.fvel[nf][alt1],ngbr.fvel[nf][alt2]);
					} else {
						godunov_set_velocity(f.fvel[ff][dir],f.fvel[ff][alt1],
							f.fvel[ff][alt2],ngbr.fvel[nf][dir],
							ngbr.fvel[nf][alt1], ngbr.fvel[nf][alt2]);
					}
#else
                    godunov_set_velocity(f.fvel[ff][dir],f.fvel[ff][alt1],
                        f.fvel[ff][alt2],ngbr.fvel[nf][dir],
                        ngbr.fvel[nf][alt1], ngbr.fvel[nf][alt2]);

#endif
				}

			}
		}
	}
}

#if defined (TANG)
/** Compute tangential stability information and store in tang[3][3] in fluid and solid fields.
 */
void fluid_3d::compute_tang_derivatives(bool verbose){
    if(verbose) printf("Rank %d. Compute tangential derivative information.\n",rank);
	const int strides[3]={1,sm4,smn4};
	const double dhsp[3] = {dxsp, dysp, dzsp};
	// after vels are determined by upwinding, compute the tangential terms in advection
	for (int kk =0; kk < so; kk++) {
		for (int jj =0; jj < sn; jj++) {
			for (int ii =0; ii < sm; ii++) {
				int eid = index(ii,jj,kk);
				double adv_vel[3] = {0,0,0};
				field &f = u0[eid];
				ref_map &sf = rm0[eid];
				// adv_vel[3] is being written by f.tang_terms()
				// and used in sf.tang_terms()
				f.tang_terms(dhsp,strides,adv_vel);
				sf.tang_terms(dhsp,strides,adv_vel);
			}
		}
	}
}
#endif

/** Compute all the terms contribute to accelearation at t=time.
 * 	Result doesn't have dt multiplied.
 * \param[in] ijk field node position relative to u0.
 * \param[in] (myx,myy,myz) the physical position of field node, for computing accelecrations.
 * \param[in,out] acc the acceleration due to stress, pressure, bc, anchors.
 * \param[in] pres if pressure gradient needs to be computed.
 * \param[in] gdn_ex if the acceleration term is meant for godunov extrapolation.
 * \param[in] verbose whether to output debug info.
 */
void fluid_3d::acceleration(int ijk,double myx,double myy, double myz,double (&acc)[3],bool pres,bool gdn_ex, bool verbose){
	field *fp = u0+ijk;
	int eid = ijk + G0;
	// compute acceleration due to stress imbalance.
	acc[0] += dxsp*(fp[1].sigma[0][0] - fp->sigma[0][0])
				+ dysp*(fp[sm4].sigma[1][0] -fp->sigma[1][0])
				+ dzsp*(fp[smn4].sigma[2][0] - fp->sigma[2][0]);

	acc[1] += dxsp*(fp[1].sigma[0][1] - fp->sigma[0][1])
				+ dysp*(fp[sm4].sigma[1][1] -fp->sigma[1][1])
				+ dzsp*(fp[smn4].sigma[2][1] - fp->sigma[2][1]);

	acc[2] += dxsp*(fp[1].sigma[0][2] - fp->sigma[0][2])
				+ dysp*(fp[sm4].sigma[1][2] -fp->sigma[1][2])
				+ dzsp*(fp[smn4].sigma[2][2] - fp->sigma[2][2]);

	//calculate pressure gradient
	double tmp[3] = {0,0,0};
	if(pres) neg_pres_grad(fp,tmp);
	for(int i=0;i<3;i++) {
        acc[i] += tmp[i];
#if defined(VAR_DEN)
        acc[i] /= lrho[eid];
#endif
    }

    for(int i=0;i<3;i++) {
        tmp[i] = 0.;
    }
	double anchor_acc[3]={0,0,0};
	double sfrac=0, lfrac, phiv;
	// solid acceleration
	for(int s=0;s<mgmt->n_obj;s++){
		ref_map *rp = get_refmap(eid,s+1);
		if(rp!=NULL){
			phiv=rp->phi(mgmt);
			if(phiv<mgmt->eps){
				lfrac = mgmt-> heaviside(phiv);
				mgmt->solid_acceleration(s, rp->x, myx, myy, myz, time, anchor_acc[0], anchor_acc[1], anchor_acc[2]);
				for(int nn=0;nn<3;nn++) tmp[nn]+=anchor_acc[nn];
				sfrac += lfrac;
			}
		}
	}
    if(sfrac>1.) for(int nn=0;nn<3;nn++) tmp[nn] /= sfrac;

    for(int nn=0;nn<3;nn++) acc[nn] += tmp[nn];

    if(sfrac<1.) {
        // get the fluid acceleration
        mgmt->fluid_acceleration(myx,myy,myz,time,anchor_acc[0], anchor_acc[1], anchor_acc[2]);
        for(int nn=0;nn<3;nn++) acc[nn]+=anchor_acc[nn]*(1-sfrac);
    }
}

/** Compute half time step, half spatial step extrapolation of velocities, and store in fvel[6][3].
 * \param[in] cdt time step
 */
void fluid_3d::compute_half_time_edge_velocities(const double cdt, bool verbose){
    if(verbose) printf("Rank %d. Compute half time edge vels.\n",rank);
	const double dh2[3] ={dx/2, dy/2, dz/2};
	const double dt2 = cdt/2.0;
    const bool pres=true, gdn_ex=true;

	for (int kk =0; kk < so; kk++) {
		double zz = lz0[kk];
		for (int jj =0; jj < sn; jj++) {
			double yy = ly0[jj];
			for (int ii = 0; ii < sm; ii++) {
				double xx = lx0[ii];

				int eid = index(ii,jj,kk);
				// get the node we're at and zero the derivs
				field *fp = u0+eid;

				double acc[3] = {0,0,0};
				// face velocities extrapolation takes pressure gradient into account
                // gdn_ex is not doing anything right now
				acceleration(eid,xx,yy,zz,acc,pres,gdn_ex);
				// add all contributions to the extrapolation term F. See Yu (2003) Eqn 3.13
                // The macro TANG controls if tangential stability is used
                // If it's turned to TRUE, it will be used.
				fp->godunov_extrap(dh2, dt2, acc);
			}
		}
	}
}

/** Compute half time step, half spatial step extrapolation of velocities, and store in fvel[6]p[3].
 * \param[in] cdt time step
 */
void fluid_3d::compute_half_time_edge_refmap(const double cdt, bool verbose){
	if(verbose) printf("Rank %d. Compute half time edge reference map.\n",rank);
	const int strides[3] = {1,sm4,smn4};
	const double dh2[3] ={dx/2, dy/2, dz/2};
	const double dt2 = cdt/2.0;

	for (int kk =0; kk < so; kk++) {
		for (int jj =0; jj < sn; jj++) {
			for (int ii = 0; ii < sm; ii++) {

				int eid = index(ii,jj,kk);
				field &f = u0[eid];
				ref_map &sf = rm0[eid];
				// if there's an object at primary grid point
				if(sf.oid()>0){
					sf.godunov_extrap(dh2, dt2, f);
				}
				int N = extraps.n0[eid];
				for(int ll=0;ll<N;ll++){
					ref_map &ext_sf = extraps.f0[eid][ll];
					if(ext_sf.abut_primary(eid, rm0, strides)){
                        // The macro TANG controls if tangential stability is used
                        // If it's turned to TRUE, it will be used.
						ext_sf.godunov_extrap(dh2,dt2,f);
					}
				}
			}
		}
	}
}

#if defined (VAR_DEN)
/** Update the density field at each grid cell. */
void fluid_3d::update_rho(){
    // TODO do we need to communicate? I don't think so, because we compute into the ghost node
    for(int k=0;k<so4;k++) {
        for(int j=0;j<sn4;j++) {
            for(int i=0;i<sm4;i++) {
                double sfrac = 0;
                double srho = 0;
                int ind = index(i,j,k);
                {
                    ref_map &rf = rm_mem[ind];

                    // Get the solid fraction for all of the solids
                    if(rf.oid()!=0) {
                        double tmp_sfrac = mgmt->heaviside( rf.phi(mgmt) );
                        sl_mat &tmp_sm = mgmt->sm_array[ rf.oid()-1 ];
                        srho += tmp_sfrac * tmp_sm.rho;
                        sfrac+= tmp_sfrac;
                    }
                }

                for(int s=0;s<extraps.n[ind];s++) {
                    ref_map &rf = extraps.f[ind][s];
                    double tmp_sfrac = mgmt->heaviside( rf.phi(mgmt) );
                    sl_mat &tmp_sm = mgmt->sm_array[ rf.oid()-1 ];
                    srho += tmp_sfrac * tmp_sm.rho;
                    sfrac+= tmp_sfrac;
                }
                // If sfrac exceeds 1, we have to normalize
                if(sfrac>1.) {
                    srho /= sfrac;
                } else {
                    srho += (1-sfrac) * mgmt->fm.rho;
                }
                // Update the rho field
                lrho[ind] =  srho;
            }
        }
    }

}
#endif

/** Compute the half timestep reference map.
 * \param[in] cdt custom time step.
 */
int fluid_3d::update_reference_map(const double cdt, bool verbose){
	if(verbose) printf("Rank %d. Update reference to half step.\n",rank);

	//compute half time step adv term. See Yu (2003) Eqn 3.11
	for (int kk =0; kk < so; kk++) {
		for (int jj =0; jj < sn; jj++) {
			for (int ii = 0; ii < sm; ii++) {
				int ind = index(ii, jj, kk);
				field &f = u0[ind];
				ref_map &sf = rm0[ind];

				// if we're in the primary grid
				if (sf.oid() > 0 ) {
					// compute prefactor
					double avgs [3] = {0,0,0};
					avgs[0] =-0.5*cdt*(f.fvel[0][0] + f.fvel[1][0])*dxsp;
					avgs[1] =-0.5*cdt*(f.fvel[2][1] + f.fvel[3][1])*dysp;
					avgs[2] =-0.5*cdt*(f.fvel[4][2] + f.fvel[5][2])*dzsp;

					double dxtmp[3] = {0,0,0};

					// apply advective operator for each reference map component
					for(int dir=0; dir<3; dir++) {
						dxtmp[dir] =
							avgs[0]*(sf.fx[1][dir]-sf.fx[0][dir])
							+ avgs[1] *(sf.fx[3][dir] - sf.fx[2][dir])
							+ avgs[2] *(sf.fx[5][dir] - sf.fx[4][dir]);
					}
#if defined(DEBUG)
                    if( std::isnan(dxtmp[0]) ||
                        std::isnan(dxtmp[1]) ||
                        std::isnan(dxtmp[2]) ){
                        printf("rank %d at (%d %d %d) dx is nan (%g %g %g)\n", rank, ii, jj, kk, dxtmp[0], dxtmp[1], dxtmp[2]);
                        printf("avgs %g %g %g\n", avgs[0], avgs[1], avgs[2]);
                        printf("face vel (%g %g %g), dhsp %g %g %g\n", f.fvel[0][0], f.fvel[0][1], f.fvel[0][2], dxsp, dysp, dzsp);
                    }
#endif

                    // sets xpred = x + dx
					sf.pre_update_ref_map(dxtmp);
				}
			}
		}
	}

	// communicate 8-obj id, 16-x and 32-xpred
    communicate<8>();
	communicate<16|32>();

	// copy old extrapolations to a temporary structure, clear the primary container
	extraps.copy(temp_extraps);
	extraps.reset();

	// carry out extrapolation, stored in xpred values
    extrapolate(false);

    int prob_co=0;
	// apply corrector to all, including ghost regions (to get communicated extraps)
	for (int kk=-2; kk<so+2; kk++) {
		for (int jj =-2; jj < sn+2; jj++) {
			for (int ii = -2; ii < sm+2; ii++) {
				int ind = index(ii,jj,kk);
				ref_map &sf = rm0[ind];
				if(sf.oid()) sf.set_half_timestep();

				// apply corrector to extrapolated fields as well
				for(int ll=0;ll<extraps.n0[ind];ll++){
					ref_map &ext_sf = extraps.f0[ind][ll];
					int obj_id = ext_sf.oid();
					bool found = false;
					for( int mm =0;(mm<temp_extraps.n0[ind])&&!found;mm++){
						ref_map &ext_sf_old = temp_extraps.f0[ind][mm];
						if(ext_sf_old.oid() == obj_id) {
							ext_sf.set_half_timestep(ext_sf_old.x);
							found = true;
						}
					}
					if(!found){
						double dist = ext_sf.phi(mgmt,1);

                        if(dist< mgmt->eps && ext_sf.lid()<=3){
							printf("fluid_3d::Rank %d update_reference_map() fatal error: "
									"New extrap can't find old extrap. At local coordinate (%d,%d,%d),"
									"global coordinate (%d,%d,%d).\n"
									"For object %d: phi=%8g, (TZW=%8g), xpred(%g %g %g), layer %d\n",
                                    rank,
									ii,jj,kk,ai+ii,aj+jj,ak+kk,
									obj_id, ext_sf.phi(mgmt,1), mgmt->eps,
									ext_sf.xpred[0], ext_sf.xpred[1], ext_sf.xpred[2],
									ext_sf.lid());
                            prob_co++;
						} else {
							ext_sf.set_x(ext_sf.xpred);
						}

					}
				}
			}
		}
	}

#if defined(VAR_DEN)
    // We blend the rho to set an updated field, based on extrapolated field values
    update_rho();
#endif
	return prob_co;
}

void fluid_3d::revert_reference_map_update(){
    // revert the mid-timestep averaging, to be output into debugging files.
	for (int kk=-2; kk<so+2; kk++) {
		for (int jj =-2; jj < sn+2; jj++) {
			for (int ii = -2; ii < sm+2; ii++) {
				int ind = index(ii,jj,kk);
				ref_map &sf = rm0[ind];
				if(sf.oid()) sf.unset_half_timestep();

				// apply corrector to extrapolated fields as well
				for(int ll=0;ll<extraps.n0[ind];ll++){
					ref_map &ext_sf = extraps.f0[ind][ll];
					int obj_id = ext_sf.oid();
                    bool found=false;
					for( int mm =0;(mm<temp_extraps.n0[ind] && !found);mm++){
						ref_map &ext_sf_old = temp_extraps.f0[ind][mm];
						if(ext_sf_old.oid() == obj_id) {
                            found=true;
							ext_sf.unset_half_timestep(ext_sf_old.x);
						}
					}
				}
			}
		}
	}
}

/** Transfer the full time step reference map update from xpred to x. */
void fluid_3d::update_reference_map_full(bool verbose){
	if(verbose) printf("Rank %d. Update reference map full.\n",rank);
	for (int kk =-2; kk < so+2; kk++) {
		for (int jj =-2; jj < sn+2; jj++) {
			for (int ii = -2; ii < sm+2; ii++) {
				int rid = index(ii, jj, kk);
				ref_map &sf = rm0[rid];
				if (sf.oid() > 0) sf.set_full_timestep();
				// apply full step to extrapolated fields as well
				for(int ll=0;ll<extraps.n0[rid];ll++){
					ref_map &ext_sf = extraps.f0[rid][ll];
					ext_sf.set_full_timestep();
				}
			}
		}
	}

#if defined(VAR_DEN)
    // We blend the rho to set an updated field, based on extrapolated field values
    update_rho();
#endif
}

/** Compute intermediate velocity, which incoporates advection term, viscous term, and source term. */
void fluid_3d::compute_ustar(const double cdt, bool verbose){

	if(verbose) printf("Rank %d. Compute u star.\n",rank);
	for (int kk =0; kk < so; kk++) {
		double zz = lz0[kk];
		for (int jj =0; jj < sn; jj++) {
			double yy = ly0[jj];
			for (int ii = 0; ii < sm; ii++) {
				double xx = lx0[ii];
				int eid = index(ii, jj, kk);
				field *fp = u0+eid;
				field &f = *fp;

				double acc[3] = {0,0,0};
				// ustar calculation doesn't take pressure gradient into account
				// also if implicit stress stencil is used, acc doesn't get doubled
                bool gdn = false;
				acceleration(eid,xx,yy,zz,acc,pres_update,gdn);
				// add all contributions to the extrapolation term F. See Yu (2003) Eqn 3.13

				// compute prefactor
				double avgs [3] = {0,0,0};
				avgs[0] =-0.5*cdt*(f.fvel[0][0] + f.fvel[1][0])*dxsp;
				avgs[1] =-0.5*cdt*(f.fvel[2][1] + f.fvel[3][1])*dysp;
				avgs[2] =-0.5*cdt*(f.fvel[4][2] + f.fvel[5][2])*dzsp;

				double dvel[3] = {0,0,0};
				// apply advective operator for each vel component
				for(int dir=0; dir<3; dir++) {
					dvel[dir] = avgs[0]*(f.fvel[1][dir]-f.fvel[0][dir]) + avgs[1] *(f.fvel[3][dir] - f.fvel[2][dir]) + avgs[2] *(f.fvel[5][dir] - f.fvel[4][dir]);
					acc[dir] *= cdt;
				}
				for(int dir=0; dir<3; dir++) {
					dvel[dir] += acc[dir];
                }
				f.reset_vel_derivs();
				f.add_vel_derivs(dvel);
				f.update_vels();
			}
		}
	}
	if(impl) impl_timestep();

}

/** New routine to compute both solid and fluid stresses.
 * Fluid stress comes from Newtonian fluid model, boils down to nabla_i u_j;
 * while solid stress comes incompressible neo-Hookian model, which stems from
 * Deformation gradient F nabla_i X_j.
 * F is the inverse of grad Xi.
 */
void fluid_3d::compute_stress(bool verbose){

    watch.tic(3);
    if(verbose) printf("Rank %d. Compute stress.\n", rank);

    // In order to extrapolate into the solid, we need grad_phi calculations in the ghost regions
	for(int k=-1; k<=so+1; k++){
		// Note possibility to be multithreaded here
		for(int j=-1;j<=sn+1;j++){
			for(int i=-1;i<=sm+1;i++){
				const int eid = index(i,j,k);
				field *fp = u0+eid;
				int s=0;
				double ss[3]={0,0,0};
				double fs[3]={0,0,0};
				double sfrac=0;
                int err = 0;
				// set left edge stress
                solid_stress<LEFT>(eid, ss, sfrac);
                if(sfrac>=1){
                    // set the stresses corresponding sigma_ij
                    // different from the f2d convention
                    // we arrange the stresses row wise
                    for(s=0;s<3;s++){
                        fp->sigma[LEFT][s] = ss[s]/sfrac;
                    }
                }
                else {
                    fluid_stress<LEFT>(eid,fs);
                    for(s=0;s<3;s++){
                        fp->sigma[LEFT][s] = (1-sfrac)*fs[s] + ss[s];
                    }
                }

                err = collision_stress<LEFT>(eid, fp->sigma[LEFT]);
                if(err!=0) {
                    printf("Something went wrong in collision_stress<LEFT> at position (%d %d %d)\n", i,j,k);
                    exit(1);
                }

                sfrac=0;
                for(s=0;s<3;s++){
                    ss[s]=0; fs[s]=0;
                }

                // set front edge stress
                solid_stress<FRONT>(eid, ss, sfrac);
                if(sfrac>=1){
                    // set the stresses corresponding sigma_ij
                    // different from the f2d convention
                    // we arrange the stresses row wise
                    for(s=0;s<3;s++){
                        fp->sigma[FRONT][s] = ss[s]/sfrac;
                    }
                }
                else {
                    fluid_stress<FRONT>(eid,fs);
                    for(s=0;s<3;s++){
                        fp->sigma[FRONT][s] = (1-sfrac)*fs[s] + ss[s];
                    }
                }

                err = collision_stress<FRONT>(eid, fp->sigma[FRONT]);
                if(err!=0){
                    printf("Something went wrong in collision_stress<FRONT> at position (%d %d %d)\n", i,j,k);
                    exit(1);
                }

                sfrac=0;
                for(s=0;s<3;s++){
                    ss[s]=0; fs[s]=0;
                }

                // set down edge stress
                solid_stress<DOWN>(eid, ss, sfrac);
                if(sfrac>=1){
                    // set the stresses corresponding sigma_ij
                    // different from the f2d convention
                    // we arrange the stresses row wise
                    for(s=0;s<3;s++){
                        fp->sigma[DOWN][s] = ss[s]/sfrac;
                    }
                }
                else {
                    fluid_stress<DOWN>(eid,fs);
                    for(s=0;s<3;s++){
                        fp->sigma[DOWN][s] = (1-sfrac)*fs[s] + ss[s];
                    }
                }
                err = collision_stress<DOWN>(eid, fp->sigma[DOWN]);

                if(err!=0){
                    printf("Something went wrong in collision_stress<DOWN> at position (%d %d %d)\n", i,j,k);
                    exit(1);
                }

			}// end of i
		}// end of j
	}// end of k

    watch.toc(3);
}

/** Routine to compute collision stress.
 * only on the lower faces. */
template<lower_faces F>
int fluid_3d::collision_stress(int eid, double (&fluid_s)[3]){
	if(F>2) p_fatal_error("collision_stress: Unknown lower_faces, try again.\n", 1);
	const int strides[3] ={1, sm4, smn4};
	const double facs[3] ={dxsp, dysp, dzsp};

	// compute strides in this direction
	const int strd = strides[F];
	const int ad1 = (F+1)%3, ad2= (F+2)%3, as1 = strides[ad1], as2 = strides[ad2];

	// pad it with the index to first non_ghost node
	int ind = G0 + eid;
    int locpos[3]={0,0,0};
    unpack_index(eid, locpos[0], locpos[1], locpos[2]);
	// Since the parent function is called on all cells, we need to be careful
	// in uing get_refmap() because if a grid point is on the edge (in outer ghost layer)
	// and get_refmap(ind-strd) is call, it has a chance to return something non NULL
	// That might just mess up the calculation here.
	// Real objects cannot overlap (but they can in fluid_2d)
	// So collision stresses need only be computed
	// when 1 extrapolated 1 real overlap, or two extrapolated regions overlap
    int p_objs = int(rm0[eid].oid()!=0);
	int e_objs = extraps.n0[eid];
    int nobjs=e_objs + p_objs;

	// if there's only one object in total, we move on
	if(nobjs<=1) return 0;

    // Put the pointers of the reference map into an array, for later
    ref_map ** collision_objs = new ref_map*[nobjs];
    // Total collision stress
    sym_matrix total_coll_sigma (0);

    for(int i=0;i<p_objs;i++){
        collision_objs[i] = rm0+eid;
    }
    for(int i=p_objs;i<nobjs;i++) {
        collision_objs[i] = extraps.f0[eid]+ (i-p_objs);
    }

	for(int ii=0;ii<nobjs-1;ii++) {

        int ngbrpos[3] = {0,0,0};
		ref_map &rf0 = *(collision_objs[ii]);
		int obj_id0 = rf0.oid();

		double G0 = mgmt->sm_array[obj_id0-1].G;
        for(int tmps=0;tmps<3;tmps++) ngbrpos[tmps] = locpos[tmps];
        ngbrpos[F] -= 1;
		ref_map *rn0 = get_refmap(ind-strd, obj_id0, ngbrpos);
		if(rn0 == NULL) continue;

		double phiv0 = 0.5*(rf0.phi(mgmt) + rn0->phi(mgmt));
		if(phiv0 > mgmt->eps) continue;

		double dphi0[3];
		// can't skip just yet, if this is outside of one object's transition zone
		// since it could be within the other objects' transition zone
		// if there are neighbors, we grab them
        for(int tmps=0;tmps<3;tmps++) ngbrpos[tmps] = locpos[tmps];
        ngbrpos[ad1] += 1;
		ref_map *r1p = get_refmap(ind+as1, obj_id0, ngbrpos);

        for(int tmps=0;tmps<3;tmps++) ngbrpos[tmps] = locpos[tmps];
        ngbrpos[ad1] += 1;
        ngbrpos[F] -= 1;
		ref_map *r1p1 = get_refmap(ind+as1-strd, obj_id0, ngbrpos);

        for(int tmps=0;tmps<3;tmps++) ngbrpos[tmps] = locpos[tmps];
        ngbrpos[ad1] -= 1;
		ref_map *r1m = get_refmap(ind-as1,obj_id0, ngbrpos);

        for(int tmps=0;tmps<3;tmps++) ngbrpos[tmps] = locpos[tmps];
        ngbrpos[ad1] -= 1;
        ngbrpos[F] -= 1;
		ref_map *r1m1 = get_refmap(ind-as1-strd, obj_id0, ngbrpos);

        for(int tmps=0;tmps<3;tmps++) ngbrpos[tmps] = locpos[tmps];
        ngbrpos[ad2] += 1;
		ref_map *r2p = get_refmap(ind+as2, obj_id0, ngbrpos);

        for(int tmps=0;tmps<3;tmps++) ngbrpos[tmps] = locpos[tmps];
        ngbrpos[ad2] += 1;
        ngbrpos[F] -= 1;
		ref_map *r2p1 = get_refmap(ind+as2-strd, obj_id0, ngbrpos);

        for(int tmps=0;tmps<3;tmps++) ngbrpos[tmps] = locpos[tmps];
        ngbrpos[ad2] -= 1;
		ref_map *r2m = get_refmap(ind-as2, obj_id0, ngbrpos);

        for(int tmps=0;tmps<3;tmps++) ngbrpos[tmps] = locpos[tmps];
        ngbrpos[ad2] -= 1;
        ngbrpos[F] -= 1;
		ref_map *r2m1 = get_refmap(ind-as2-strd, obj_id0, ngbrpos);


		if(check_null(r1p, r1p1, r1m, r1m1, r2p, r2p1, r2m, r2m1)) continue;

        // This is already calculated
        // for(int i=0;i<3;i++) dphi0[i] = rf0.grad_phi[F][i];

	    dphi0[F] = facs[F]*(rf0.phi(mgmt) - rn0->phi(mgmt));
	    dphi0[ad1] = 0.25*facs[ad1]*(r1p->phi(mgmt) + r1p1->phi(mgmt) - r1m->phi(mgmt) - r1m1->phi(mgmt));
	    dphi0[ad2] = 0.25*facs[ad2]*(r2p->phi(mgmt) + r2p1->phi(mgmt) - r2m->phi(mgmt) - r2m1->phi(mgmt));

        double norm = 0;
#if defined(DEBUG)
        double diff_in_norm = 0;
#endif
        // we normalize the gradient here first
        for(int i=0;i<3;i++) norm += dphi0[i] * dphi0[i];
        if(norm>0) {
            norm = sqrt(norm);
            for(int i=0;i<3;i++) dphi0[i] /= norm;
        } else {
            printf("collision_stress: you've got a serious issue at %d. Gradient has norm = %g! rank %d oid %d lid %d\n", eid, norm, rank, obj_id0, rf0.lid());
        }

#if defined(DEBUG)
        for(int i=0;i<3;i++) diff_in_norm += (dphi0[i]-rf0.grad_phi[F][i]) * (dphi0[i]-rf0.grad_phi[F][i]);
        diff_in_norm = sqrt(diff_in_norm);

        rf0.grad_phi_diff[F] = diff_in_norm;
#endif

		for(int jj=ii+1;jj<nobjs;jj++){

			ref_map &rf1 = *(collision_objs[jj]);
			int obj_id1 = rf1.oid();
			ref_map *rn1 = get_refmap(ind-strd, obj_id1);
			if(rn1 == NULL) continue;

			double phiv1 = 0.5*(rf1.phi(mgmt) + rn1->phi(mgmt));
            // if the width between the objects is wider than 1 transition zone
            // we don't compute collision stress
            // This is an attempt to eliminate the collision "gap"?
			if(phiv1 > mgmt->eps || phiv0 + phiv1 > mgmt-> eps) continue;

			double dphi1[3];
            for(int tmps=0;tmps<3;tmps++) ngbrpos[tmps] = locpos[tmps];
            ngbrpos[ad1] += 1;
            ref_map *s1p = get_refmap(ind+as1, obj_id1, ngbrpos);

            for(int tmps=0;tmps<3;tmps++) ngbrpos[tmps] = locpos[tmps];
            ngbrpos[ad1] += 1;
            ngbrpos[F] -= 1;
            ref_map *s1p1 = get_refmap(ind+as1-strd, obj_id1, ngbrpos);

            for(int tmps=0;tmps<3;tmps++) ngbrpos[tmps] = locpos[tmps];
            ngbrpos[ad1] -= 1;
            ref_map *s1m = get_refmap(ind-as1,obj_id1, ngbrpos);

            for(int tmps=0;tmps<3;tmps++) ngbrpos[tmps] = locpos[tmps];
            ngbrpos[ad1] -= 1;
            ngbrpos[F] -= 1;
            ref_map *s1m1 = get_refmap(ind-as1-strd, obj_id1, ngbrpos);

            for(int tmps=0;tmps<3;tmps++) ngbrpos[tmps] = locpos[tmps];
            ngbrpos[ad2] += 1;
            ref_map *s2p = get_refmap(ind+as2, obj_id1, ngbrpos);

            for(int tmps=0;tmps<3;tmps++) ngbrpos[tmps] = locpos[tmps];
            ngbrpos[ad2] += 1;
            ngbrpos[F] -= 1;
            ref_map *s2p1 = get_refmap(ind+as2-strd, obj_id1, ngbrpos);

            for(int tmps=0;tmps<3;tmps++) ngbrpos[tmps] = locpos[tmps];
            ngbrpos[ad2] -= 1;
            ref_map *s2m = get_refmap(ind-as2, obj_id1, ngbrpos);

            for(int tmps=0;tmps<3;tmps++) ngbrpos[tmps] = locpos[tmps];
            ngbrpos[ad2] -= 1;
            ngbrpos[F] -= 1;
            ref_map *s2m1 = get_refmap(ind-as2-strd, obj_id1, ngbrpos);

			if(check_null(s1p, s1p1, s1m, s1m1, s2p, s2p1, s2m, s2m1)) continue;

            //This is already calculated
            //for(int i=0;i<3;i++) dphi1[i] = rf1.grad_phi[F][i];

		    dphi1[F] = facs[F]*(rf1.phi(mgmt) - rn1->phi(mgmt));
		    dphi1[ad1] = 0.25*facs[ad1]*(s1p->phi(mgmt) + s1p1->phi(mgmt) - s1m->phi(mgmt) - s1m1->phi(mgmt));
		    dphi1[ad2] = 0.25*facs[ad2]*(s2p->phi(mgmt) + s2p1->phi(mgmt) - s2m->phi(mgmt) - s2m1->phi(mgmt));

            norm = 0;
            // we normalize the gradient here first
            for(int i=0;i<3;i++) norm += dphi1[i] * dphi1[i];
            if(norm>0) {
                norm = sqrt(norm);
                for(int i=0;i<3;i++) dphi1[i] /= norm;
            } else {
                printf("collision_stress: you've got a serious issue here. Gradient has 0 norm! rank %d grid point %d oid %d\n", rank, eid, obj_id1);
            }

#if defined (DEBUG)
            diff_in_norm = 0;
            for(int i=0;i<3;i++) diff_in_norm += (dphi1[i]-rf1.grad_phi[F][i]) * (dphi1[i]-rf1.grad_phi[F][i]);
            diff_in_norm = sqrt(diff_in_norm);
            rf1.grad_phi_diff[F] = diff_in_norm;
#endif

			double norm_vec [3], denom=0;
			for(int nn=0;nn<3;nn++){
				norm_vec[nn] = dphi1[nn] - dphi0[nn];
				denom += norm_vec[nn]*norm_vec[nn];
			}

            // Divide out the norm of the difference in grad phi vector
            denom = sqrt(denom);

            // if the this vector is very small, we ignore it
            // this means the outward normal is very co-linear
            if(denom<small_number) continue;

			for(int nn=0;nn<3;nn++){
				norm_vec[nn]/=denom;
			}

            double G1 = mgmt->sm_array[obj_id1-1].G;
            // A multiplier is included in coll_func, set in sim_params class
            // the new coll_fun also limits the linear response function between [0,1]
            // going from outside of transition zone to inside

	    double coll_fac = G0*mgmt->coll_func(phiv0) + G1*mgmt->coll_func(phiv1);

            sym_matrix coll_sigma(0);
            for(int i=0;i<3;i++){
                for(int j=0;j<3;j++){
                    coll_sigma(i,j) = -coll_fac*norm_vec[i] * norm_vec[j];
                }
            }

            double sigma_avg = coll_sigma.trace()/3.;
            for(int j=0;j<3;j++){
                coll_sigma(j,j) -= sigma_avg;
            }

            for(int i=0;i<3;i++){
                for(int j=0;j<3;j++){
                    total_coll_sigma(i,j) += coll_sigma(i,j);
                }
            }

            for(int i=0;i<3;i++) {
                fluid_s[i] += coll_sigma(F,i);
            }

		}
	}

#if defined(DEBUG)
    // After the total collision stress has been completely computed
    // we can start calculating the normal and modulus of the shear components
    // for each object
    double sum_solid_traction[3] = {0,0,0};
    for(int ii=0;ii<nobjs-1;ii++) {

        ref_map &rf = *(collision_objs[ii]);
        int obj_id = rf.oid();

        ref_map *rn = get_refmap(ind-strd, obj_id);
        if(rn == NULL) continue;

        double phiv = 0.5*(rf.phi(mgmt) + rn->phi(mgmt));
		if(phiv > mgmt->eps) continue;
        // We compute the following
        rf.coll_stress_normal[F] = 0;
        rf.coll_stress_shear[F] = 0;

        double solid_traction[3] = {0,0,0};
        for (int i=0;i<3;i++){
            for(int j=0;j<3;j++){
                rf.coll_stress_normal[F] += total_coll_sigma(i,j)*rf.grad_phi[F][j] * rf.grad_phi[F][i];
                solid_traction[i] += total_coll_sigma(i,j) * rf.grad_phi[F][j];
            }
        }

        for (int i=0;i<3;i++){
            sum_solid_traction[i] += solid_traction[i];
            solid_traction[i] -= rf.coll_stress_normal[F] * rf.grad_phi[F][i];
            rf.coll_stress_shear[F] += solid_traction[i] * solid_traction[i];
        }

        rf.coll_stress_shear[F] = sqrt(rf.coll_stress_shear[F]);
    }
    double solid_traction_norm = 0;
    for (int i=0;i<3;i++){
        solid_traction_norm = sum_solid_traction[i] * sum_solid_traction[i];
    }
    solid_traction_norm = sqrt(solid_traction_norm);

    for(int ii=0;ii<nobjs-1;ii++) {

        ref_map &rf = *(collision_objs[ii]);
        int obj_id = rf.oid();

        ref_map *rn = get_refmap(ind-strd, obj_id);
        if(rn == NULL) continue;

        double phiv = 0.5*(rf.phi(mgmt) + rn->phi(mgmt));
		if(phiv > mgmt->eps) continue;
        rf.coll_traction_tot[F] = solid_traction_norm;

    }
#endif

    delete [] collision_objs;
	return 0;
}

/** New routine to compute solid stress;
 * only on the lower faces. */
template<lower_faces F>
void fluid_3d::solid_stress(int eid, double (&solid_s)[3], double &sfrac){
	if(F>2) p_fatal_error("solid_stress: Unknown lower_faces, try again.\n", 1);

	const int strides[3] ={1, sm4, smn4};
	const double facs[3] ={dxsp, dysp, dzsp};

	// compute strides in this direction
	const int strd = strides[F];
	const int ad1 = (F+1)%3, ad2= (F+2)%3, as1 = strides[ad1], as2 = strides[ad2];

	// pad it with the index to first non_ghost node
	int ind = G0 + eid;
    int locpos[3]={0,0,0};
    unpack_index(ind, locpos[0], locpos[1], locpos[2]);
	// Since the parent function is called on all cells, we need to be careful
	// in uing get_refmap() because if a grid point is on the edge (in outer ghost layer)
	// and get_refmap(ind-strd) is call, it has a chance to return something non NULL
	// That might just mess up the calculation here.

	const field *fp = u0+ eid;
	// ii=-1 get primary ref map, then ref map from the multimap object
	for(int ii=-1;ii<extraps.n0[eid];ii++){

        int ngbrpos[3] = {0,0,0};
        ref_map *rp;

		if(ii==-1) rp=rm0+eid;
		else rp = extraps.f0[eid]+ii;

		ref_map &rf = *rp;
		int obj_id = rf.oid();

#if defined(DEBUG)
        // Initialize the debugging fields to something
        // These are set to -1.e200 for each reference map field variable
        // that exist at this grid point
        rf.solid_stress_normal[F] = -1.e200;
        rf.solid_stress_shear[F] = -1.e200;
        rf.coll_stress_normal[F] = -1.e200;
        rf.coll_stress_shear[F] = -1.e200;
        rf.grad_phi_diff[F] = -1.e200;
        rf.coll_traction_tot[F] = -1.e200;
#endif

        for(int i=0;i<3;i++) rf.grad_phi[F][i] = -1.e200;
        // Needed for macro.dat output
        rf.elastic_energy[F] = -1.e200;
        rf.detF[F] = -1.e200;
		// get pointers to all the neighbors needed
		// that match the object id
		if(obj_id==0) continue;

		// left neighbor to compute edge phi value
        for(int tmps=0;tmps<3;tmps++) ngbrpos[tmps] = locpos[tmps];
        ngbrpos[F] -= 1;
		ref_map *rn = get_refmap(ind-strd, obj_id, ngbrpos);

		if(rn==NULL) continue;
		double phiv = 0.5*(rf.phi(mgmt) + rn->phi(mgmt));
        if(phiv > mgmt->eps) continue;


        // compute the solid fraction here
		double tf = mgmt->heaviside(phiv);
		sfrac += tf;
        sl_mat & tmp_sm = mgmt->sm_array[obj_id-1];
		double G = tmp_sm.G * tf;

		// if the edge is within transition zone, need to grab some neighbors
        for(int tmps=0;tmps<3;tmps++) ngbrpos[tmps] = locpos[tmps];
        ngbrpos[ad1] += 1;
		ref_map *r1p = get_refmap(ind+as1, obj_id, ngbrpos);

        for(int tmps=0;tmps<3;tmps++) ngbrpos[tmps] = locpos[tmps];
        ngbrpos[ad1] += 1;
        ngbrpos[F] -= 1;
		ref_map *r1p1 = get_refmap(ind+as1-strd, obj_id, ngbrpos);

        for(int tmps=0;tmps<3;tmps++) ngbrpos[tmps] = locpos[tmps];
        ngbrpos[ad1] -= 1;
		ref_map *r1m = get_refmap(ind-as1,obj_id, ngbrpos);

        for(int tmps=0;tmps<3;tmps++) ngbrpos[tmps] = locpos[tmps];
        ngbrpos[ad1] -= 1;
        ngbrpos[F] -= 1;
		ref_map *r1m1 = get_refmap(ind-as1-strd, obj_id, ngbrpos);

        for(int tmps=0;tmps<3;tmps++) ngbrpos[tmps] = locpos[tmps];
        ngbrpos[ad2] += 1;
		ref_map *r2p = get_refmap(ind+as2, obj_id, ngbrpos);

        for(int tmps=0;tmps<3;tmps++) ngbrpos[tmps] = locpos[tmps];
        ngbrpos[ad2] += 1;
        ngbrpos[F] -= 1;
		ref_map *r2p1 = get_refmap(ind+as2-strd, obj_id, ngbrpos);

        for(int tmps=0;tmps<3;tmps++) ngbrpos[tmps] = locpos[tmps];
        ngbrpos[ad2] -= 1;
		ref_map *r2m = get_refmap(ind-as2, obj_id, ngbrpos);

        for(int tmps=0;tmps<3;tmps++) ngbrpos[tmps] = locpos[tmps];
        ngbrpos[ad2] -= 1;
        ngbrpos[F] -= 1;
		ref_map *r2m1 = get_refmap(ind-as2-strd, obj_id, ngbrpos);

		if(check_null(r1p, r1p1, r1m, r1m1, r2p, r2p1, r2m, r2m1)) continue;

		rf.grad_phi[F][F]   = facs[F]*(rf.phi(mgmt) - rn->phi(mgmt));
		rf.grad_phi[F][ad1] = 0.25*facs[ad1]*(r1p->phi(mgmt) + r1p1->phi(mgmt) - r1m->phi(mgmt) - r1m1->phi(mgmt));
		rf.grad_phi[F][ad2] = 0.25*facs[ad2]*(r2p->phi(mgmt) + r2p1->phi(mgmt) - r2m->phi(mgmt) - r2m1->phi(mgmt));
		rf.normalize_grad_phi(F);

		// if get to this point, all neighbors exist
		double dR[3][3];
		for(int nn=0;nn<3;nn++){
			dR[F][nn] = facs[F]*(rf.x[nn] - rn->x[nn]);
			dR[ad1][nn] = 0.25*facs[ad1]*(r1p->x[nn] + r1p1->x[nn] - r1m->x[nn] - r1m1->x[nn]);
			dR[ad2][nn] = 0.25*facs[ad2]*(r2p->x[nn] + r2p1->x[nn] - r2m->x[nn] - r2m1->x[nn]);
		}


        matrix grad_xi(0);
        grad_xi.initialize(dR[0][0], dR[1][0], dR[2][0], dR[0][1], dR[1][1], dR[2][1], dR[0][2], dR[1][2], dR[2][2]);

		matrix def_grad = grad_xi.inverse();

		sym_matrix astress(0);
		bool add_astress=false;
		// if active stress, calculate on this face
		if (mgmt->objs[obj_id-1]->act_stress_method) {

			ref_map &me = rf, &you = *rn;


			double Xx = 0.5*(me.x[0]+you.x[0]);
			double Xy = 0.5*(me.x[1]+you.x[1]);
			double Xz = 0.5*(me.x[2]+you.x[2]);

			double rmX[3] = {Xx,Xy,Xz};

			mgmt->objs[obj_id-1]->active_stress(time,rmX,F,def_grad,astress);

			add_astress=true;
		}

        // UPDATE
        // HERE WE INSERT ACTIVE STRESS
        mgmt->actuate(obj_id-1, rf.x, rn->x, time, def_grad);

        rf.detF[F] = def_grad.det();

        // We are on our way to use left Cauchy-Green tensor
        // NOT THE RIGHT one
        sym_matrix sigma(0);
		sigma = def_grad.aat();

        // At this point, we can calculate the strain energy at this face
        // Using the right Cauchy-Green tensor
        sym_matrix rCGT(0);
        rCGT = def_grad.ata();
        rf.elastic_energy[F] = 0.5*G*(rCGT.trace() - 3 -2*log(rf.detF[F]));

		// store active stress work rate
		rf.act_power[F] = 0;
		if (add_astress) {

			// grab velocity gradient tensor
			matrix gv;
			velocity_grad(F,eid,gv);

			// double-dot product
			for(int r=0;r<3;r++) for(int c=0;c<3;c++)
				rf.act_power[F] -= astress(r,c) * gv(r,c);
			
		}


		// compute the hydrostatic part of the stress
		double sigma_avg = sigma.trace()/3.0;

        for(int nn=0;nn<3;nn++) sigma(nn,nn) -= sigma_avg;
		double ex_mu = tmp_sm.ex_mu * tf *(1 + tmp_sm.ev_trans_mult * mgmt->tderiv_func_in(phiv));

        // Multiply G to the strain tensor to get the stress tensor
        sigma.scale(G);

		for(int nn=0;nn<3;nn++) {
			solid_s[nn] += sigma(F,nn);
			solid_s[nn] += ex_mu * facs[F] * (
					(fp->vel[nn] - fp[-strd].vel[nn])
			+ ( (nn==F) ? (fp->vel[nn] - fp[-strd].vel[nn]) :
			  (0.25*(fp[ strides[nn]].vel[F] + fp[ strides[nn]-strd].vel[F]
				    -fp[-strides[nn]].vel[F] - fp[-strides[nn]-strd].vel[F])))
			);

			if (add_astress) solid_s[nn] += astress(F,nn);
		}

#if defined(DEBUG)
        // We then calculate the normal and shear component of the stress in the plane defined by the grad phi here
        rf.solid_stress_normal[F] = 0;
        rf.solid_stress_shear[F] = 0;

        double solid_traction[3] = {0,0,0};
        for (int i=0;i<3;i++){
            for(int j=0;j<3;j++){
                rf.solid_stress_normal[F] += sigma(i,j)*rf.grad_phi[F][j] * rf.grad_phi[F][i];
                solid_traction[i] += sigma(i,j) * rf.grad_phi[F][j];
            }
        }

        for (int i=0;i<3;i++){
            solid_traction[i] -= rf.solid_stress_normal[F] * rf.grad_phi[F][i];
            rf.solid_stress_shear[F] += solid_traction[i] * solid_traction[i];
        }
        rf.solid_stress_shear[F] = sqrt(rf.solid_stress_shear[F]);
#endif

	}
}

/** New routine to compute fluid stress;
 * only on the lower faces. */
template<lower_faces F>
void fluid_3d::fluid_stress(int eid, double (&fluid_s)[3]){
	if(F>2) p_fatal_error("fluid_stress: Unknown lower_faces, try again.\n", 1);
	const int strides[3] ={1, sm4, smn4};
	const double facs[3] ={dxsp, dysp, dzsp};
	int strd = strides[F];
	int ad1 = (F+1)%3, ad2= (F+2)%3, as1 = strides[ad1], as2 = strides[ad2];

	field *fp = u0+eid;
    int ind = G0+eid;
    int locpos[3] = {0,0,0};
    unpack_index(ind, locpos[0], locpos[1], locpos[2]);
    locpos[F] -= 1;

	double visc = mgmt->fm.mu;
	if(true){
        if(within_local_domain(locpos[0], locpos[1], locpos[2])) {
            for(int nn=0;nn<3;nn++) fluid_s[nn] = visc* facs[F] * ( fp->vel[nn] - fp[-strd].vel[nn]);
        }
	}
	else{
        // Cannot be use unless bounds are checked!
		double du[3],tu1,tu2;
		tu1 = 0.25*facs[ad1]*(fp[as1].vel[0] + fp[as1-strd].vel[0] - fp[-as1].vel[0] -  fp[-as1-strd].vel[0]);
		tu2 = 0.25*facs[ad2]*(fp[as2].vel[0] + fp[as2-strd].vel[0] - fp[-as2].vel[0] - fp[-as2-strd].vel[0]);
		du[F] *=2;
		du[ad1] += tu1;
		du[ad2] += tu2;
		for(int nn=0;nn<3;nn++) fluid_s[nn] = du[nn] * visc;
	}
}

/** routine to compute entire velocity gradient;
 *  on one lower face */
void fluid_3d::velocity_grad(lower_faces F,int eid, matrix &grad_v){
	if(F>2) p_fatal_error("fluid_stress: Unknown lower_faces, try again.\n", 1);
	const int strides[3] ={1, sm4, smn4};
	const double facs[3] ={dxsp, dysp, dzsp};
	int strd = strides[F];
	int ad1 = (F+1)%3, ad2= (F+2)%3, as1 = strides[ad1], as2 = strides[ad2];

	field *fp = u0+eid;

	// column is vel. component
	for(int c=0;c<3;c++) {

        // Cannot be use unless bounds are checked!
		grad_v(ad1,c) = 0.25*facs[ad1]*(fp[as1].vel[c] + fp[as1-strd].vel[c] - fp[-as1].vel[c] - fp[-as1-strd].vel[c]);
		grad_v(ad2,c) = 0.25*facs[ad2]*(fp[as2].vel[c] + fp[as2-strd].vel[c] - fp[-as2].vel[c] - fp[-as2-strd].vel[c]);
		grad_v(F,c)  =  facs[F]*(fp->vel[c] - fp[-strd].vel[c]);
//		du[F] *=2;
//		du[ad1] += tu1;
//		du[ad2] += tu2;
//		for(int nn=0;nn<3;nn++) fluid_s[nn] = du[nn] * visc;
	}
}

#if defined(TANG)
/**
 * Upwind to set tangential face velocities in Godunov scheme.
 * See Yu (2003) section 3.2.2.
 * \param[in,out] (norml, vl, wl) the velocites on the "lower" face,
 *                norml is the norm direction velocity, while vl, wl are tangential
 * \param[in,out] (normr, vr, wr) the velocites on the "upper" face,
 *                normr is the norm direction velocity, while vr, wr are tangential
 */
void fluid_3d::godunov_set_tang_velocity(double &norml, double &vl, double &wl, double &normr, double &vr, double &wr) {
	// in this tangential velocities upwinding scheme, store the advection velocities on the upper face
	// and the hat velocities at the lower face, arrangement see below
	/*
			 |         v_adv  |
			 |          ⬇︎     |
	  ---------------------------
			 |		          |
	 u_hat ➡︎ |⬅︎	u_adv  		  |⬅︎ u_adv
			 |⬅︎ v_hat	      |
			 |⬅︎ w_hat         |
			 |         v_adv  |
			 |         ⬇︎      |
	  ---------------------------
			 |                |
	*/
	// flow right
	if(norml + normr > 0 && norml > 0){
		normr = norml;
		vr = vl;
		wr = wl;
	// flow left
	} else if (norml + normr < 0 && normr < 0) {
		norml = normr;
        vl = vr;
        wl = wr;
	} else {
		norml = 0.5*(norml+normr);
		vr = 0.5*(vl+vr);
		wr = 0.5*(wl+wr);
		normr = 0;
	}
}
#endif

#if defined(TANG)
/**
 * Upwind to set tangential face velocities in Godunov scheme.
 * See Yu (2003) section 3.2.2.
 * \param[in,out] (norml, vl, wl) the velocites on the "lower" face,
 *                norml is the norm direction velocity, while vl, wl are tangential
 * \param[in,out] (normr, vr, wr) the velocites on the "upper" face,
 *                normr is the norm direction velocity, while vr, wr are tangential
 */
void fluid_3d::godunov_set_tang_refmap(double &norml,double &vl,double &wl,
	double &normr,double &vr,double &wr,double &mnorml,double &mvl,
	double &mwl,double &mnormr,double &mvr,double &mwr) {

	// in this tangential velocities upwinding scheme, store the advection velocities on the upper face
	// and the hat velocities at the lower face, arrangement see below
	/*
			 |         v_adv  |
			 |          ⬇︎     |
	  ---------------------------
			 |		          |
	 u_hat ➡︎ |⬅︎	u_adv  		  |⬅︎ u_adv
			 |⬅︎ v_hat	      |
			 |⬅︎ w_hat         |
			 |         v_adv  |
			 |         ⬇︎      |
	  ---------------------------
			 |                |
	*/

	// flow right
	if (norml + normr > 0 && norml > 0) {
		mnormr = mnorml;
		mvr = mvl;
		mwr = mwl;

	// flow left
	} else if (norml + normr < 0 && normr < 0){
		mnorml = mnormr;
	} else {
		mnorml = 0.5*(mnorml+mnormr);
		mvr = 0.5*(mvl+mvr);
		mwr = 0.5*(mwl+mwr);
	}
}
#endif

/**
 * Upwind to set face velocities in Godunov scheme.
 * See Yu (2003) section 3.2.2.
 * \param[in,out] (norml, vl, wl) the velocites on the "lower" face,
 *                norml is the norm direction velocity, while vl, wl are tangential
 * \param[in,out] (normr, vr, wr) the velocites on the "upper" face,
 *                normr is the norm direction velocity, while vr, wr are tangential
 */
void fluid_3d::godunov_set_velocity(double &norml,double &vl,double &wl,
	double &normr,double &vr,double &wr) {

	if(norml + normr > 0 && norml >0){
		normr = norml;
		vr = vl;
		wr = wl;
	} else if (norml + normr < 0 && normr<0) {
		norml = normr;
		vl = vr;
		wl = wr;
	} else {
		norml = normr =0;
		double tmp = 0.5*(vl+vr);
		vl = vr = tmp;
		tmp = 0.5*(wl+wr);
		wl = wr = tmp;
	}
}

/**
 * Upwind to set face velocities in Godunov scheme.
 * See Yu (2003) section 3.2.2.
 * \param[in,out] (norml, vl, wl) the velocites on the "lower" face,
 *                norml is the norm direction velocity, while vl, wl are tangential
 * \param[in,out] (normr, vr, wr) the velocites on the "upper" face,
 *                normr is the norm direction velocity, while vr, wr are tangential
 */
void fluid_3d::godunov_set_refmap(double &norml,double &vl,double &wl,
	double &normr,double &vr,double &wr,double &mnorml,double &mvl,
	double &mwl,double &mnormr,double &mvr,double &mwr) {

	if(norml + normr > 0 && norml >0) {
		mnormr = mnorml;
		mvr = mvl;
		mwr = mwl;
	} else if (norml + normr < 0 && normr<0) {
		mnorml = mnormr;
		mvl = mvr;
		mwl = mwr;
	} else {
		double tmp = 0.5*(mnorml+mnormr);
		mnorml = mnormr = tmp;
		tmp = 0.5*(mvl+mvr);
		mvl = mvr = tmp;
		tmp = 0.5*(mwl+mwr);
		mwl = mwr = tmp;
	}
}

/**
 * Solves the exact projection problem with multigrid.
 */
void fluid_3d::mac_solve(double cdt, bool verbose){
	if(verbose) printf("Rank %d. MAC projection step, multigrid solve.\n",rank);
	if(dx!=dy || dx!=dz){
        printf("%g %g %g\n", dx, dy, dz);
		p_fatal_error("fluid_3d::mac_solve(): Spatial discretization must be homogenous.", PGMG_SETUP_ERROR);
	}

    watch.tic(4);
#if defined (VAR_DEN)
    mac.set_rhoinv(lrho);
#endif
    mac.set_stencils();
    //mg_mac->reg[0]->diagnostic_A();

    fill_mac_rhs(cdt);
    double res;

	if (out==NULL || !mac.mgo->searching) {
		res = mac.run(mg_mac, false);
	} else {
		res = mac.run(mg_mac, true);
	}
    if(std::isnan(res)) p_fatal_error("fluid_3d::mac_solve(): NaN detected in multigrid solve",PGMG_MATH_ERROR);

    extract_mac_soln();
    watch.toc(4);

    communicate<2>();

    watch.tic(4);
    subtract_q_grad(cdt);
    watch.toc(4);
}

/** Fill rhs for MAC projection solve */
void fluid_3d::fill_mac_rhs(double cdt){
    // NOTE fac contain dx, instead of 1/dx because
    // We have transfer the constant factor of 1/(dx*dx) on the Laplacian operator
    // to the rhs.
	const double fac = dx/(2*cdt);

	mg_mac->reg[0]->clear_r_field();
	double *rp = mg_mac->reg[0]->r0;

	// fill the rhs of multigrid
	for(int k=0;k<so;k++) for (int j=0;j<sn;j++) for (int i=0;i<sm;i++){
		field *elem = u0+index(i,j,k);
		double *srp = rp+i+j*sm+k*sm*sn;
		//centered difference to get grad u vector at cell center, from edge velocities
		for(int f=0;f<3;f++) *srp += fac * (elem->fvel[2*f+1][f] - elem->fvel[2*f][f]);
#if defined(DEBUG)
        if(std::isnan( *srp )) {
            printf("rank %d mac rhs is nan at (%d %d %d)\n", rank, i, j, k);
        }
#endif
	}
}

/** Extract solution from Multigrid sover. */
void fluid_3d::extract_mac_soln(){

    int mg = mg_mac->reg[0]->mg, mng = mg_mac->reg[0]->mng;
	double *x0 = mg_mac->reg[0]->x0;
	for(int k=0;k<so;k++) for (int j=0;j<sn;j++) for (int i=0;i<sm;i++){
		field *elem = u0+index(i,j,k);
        elem->q=x0[i+j*mg+k*mng];
#if defined(DEBUG)
        if(std::isnan(elem->q)) {
            printf("mac solution is nan at (%d %d %d)\n", i, j, k);
        }
#endif
	}
}

/** Subtract the gradient of auxiliary pressure field. */
void fluid_3d::subtract_q_grad(double cdt){
	const double fac = 2*cdt*dxsp;
    const int strides[3] = {1, sm4, smn4};
	for(int k=0;k<so;k++) for (int j=0;j<sn;j++) for (int i=0;i<sm;i++){
		field *elem = u0+index(i,j,k);
		//centered difference to get grad q to substract from edge velocities
        //but we don't do the boundaries, because mg_mac doesn't have ghost nodes
		for(int f=0;f<3;f++) {
            int strd = strides[f];
#if defined (VAR_DEN)
            double down = fac*( elem->q - (elem-strd)->q)/( 0.5*(lrho[G0 + index(i,j,k)] + lrho[G0 + index(i,j,k) - strd]));
            double up = fac*( (elem+strd)->q - elem->q)/(0.5*(lrho[G0+index(i,j,k) + strd] + lrho[G0 + index(i,j,k)]));
#else
            double down = fac*( elem->q - (elem-strd)->q);
            double up = fac*( (elem+strd)->q - elem->q);
#endif
            elem->fvel[2*f][f] -= down;
            elem->fvel[2*f+1][f] -= up;
		}
	}
}

/**
 * Computes monotonocity limited first derivative. See Yu(2003) section 3.2.1.
 * \param[in] (u0,u1,u2,u3,u4) the five values needed to compute finite difference derivs for comparison.
 */
double fluid_3d::mono_diff(double u0, double u1, double u2, double u3, double u4){
	// compute delta_lim and centered difference at node of interest, u2
	double dlim = del_lim(u1,u2,u3);
	double cen_diff = centered_diff(u1,u3);
	// compute delta_f at neighboring node, u1 and u3
	double neigh_delf = (del_f(u0,u1,u2) + del_f(u2,u3,u4))/6.0;
	// give 4th order limited slope
	double dr = sign(cen_diff) * min(dlim, fabs(4.0/3.0*cen_diff - neigh_delf));
	return dr;
}

/**
 * Computes the delta_lim. See Yu (2003) Eqn 3.25.
 * \param[in] (u1,u2,u3).
 */
double fluid_3d::del_lim(double u0, double u1, double u2){
	double tmp1= 2*(u1-u0);
	double tmp2= 2*(u2-u1);
	if(tmp1*tmp2>0) return min(fabs(tmp1), fabs(tmp2));
	else return 0;
}

/**
 * Compiutes delta_f. See Yu(2003) Eqn. 3.26.
 */
double fluid_3d::del_f(double u0, double u1, double u2){
	double dlim = del_lim(u0,u1,u2);
	double cen_diff = centered_diff(u0,u2);
	double sg = sign(cen_diff);
	return sg*min(fabs(cen_diff), dlim);
}

/**
 * Computes negative pressure gradient, not including density division.
 * \param[in] node a ptr to the cell center to consider.
 * \param[out] acc the array containing acceleration, gets written out to
 */
void fluid_3d::neg_pres_grad(field *node, double (&acc) [3]){
		double pressure[8] = {0,0,0,0,0,0,0,0};;
		double xavg=0, yavg=0, zavg=0;
		// grab all the pressures for 8 corners
		for(int n=0; n<8; n++) pressure[n]=node[n%2 + sm4*((n/2)%2) + (n/4)*smn4].p;
		// for 4 pairs of pressures to make gradient, we add those to the averages
		for(int n=0; n<4; n++) {
			xavg+=pressure[2*n+1]-pressure[2*n];
			int yind = (n/2)*4+n%2;
			yavg+=pressure[yind+2]-pressure[yind];
			zavg+=pressure[n+4]-pressure[n];
		}
		xavg*=-0.25*dxsp; yavg*=-0.25*dysp; zavg*=-0.25*dzsp;
		acc[0]=xavg; acc[1]=yavg; acc[2]=zavg;
}
void fluid_3d::neg_pres_dgrad(field *node, double (&acc) [3]){
		double pressure[8] = {0,0,0,0,0,0,0,0};;
		double xavg=0, yavg=0, zavg=0;
		// grab all the pressures for 8 corners
		for(int n=0; n<8; n++) pressure[n]=node[n%2 + sm4*((n/2)%2) + (n/4)*smn4].dp;
		// for 4 pairs of pressures to make gradient, we add those to the averages
		for(int n=0; n<4; n++) {
			xavg+=pressure[2*n+1]-pressure[2*n];
			int yind = (n/2)*4+n%2;
			yavg+=pressure[yind+2]-pressure[yind];
			zavg+=pressure[n+4]-pressure[n];
		}
		xavg*=-0.25*dxsp; yavg*=-0.25*dysp; zavg*=-0.25*dzsp;
		acc[0]=xavg; acc[1]=yavg; acc[2]=zavg;
}


/**
 * Update the velocities (u,v,w) after derivatives are computed.
 */
void fluid_3d::expl_timestep() {
	for (int kk = 0; kk < so; kk++)
		for (int jj = 0; jj < sn; jj++)
			for (int ii =0; ii < sm; ii++) {

		// get what node we're at and forward-step using derivative vals
		u0[index(ii,jj, kk)].update_vels();
	}
}

/**
 * If implicit viscosity is enable, use multigrid to solve for velocities at the next time step.
 */
void fluid_3d::impl_timestep() {

	double *r0 = mg_visc->reg[0]->r0;
	double *x0 = mg_visc->reg[0]->x0;
	int gm4 = mg_visc->reg[0]->mg;
	int gmn4 = gm4*mg_visc->reg[0]->ng;
    double res;

	// u step
#if defined(VAR_DEN)
    visc.set_rhoinv(lrho);
#endif
	for(int comp =0;comp<3;comp++){
		for (int kk = 0; kk < so; kk++)
			for (int jj = 0; jj < sn; jj++)
				for (int ii = 0; ii < sm; ii++) {
					x0[ii + jj*gm4 + kk*gmn4] =
					r0[ii + jj*sn + kk*sm*sn] =
					u0[index(ii,jj, kk)].vel[comp];
		}
		if (out==NULL || !visc.mgo[comp].searching) {
			res=visc.run(mg_visc,false,comp);
		} else {
			res=visc.run(mg_visc,true,comp);
		}
        if(std::isnan(res)) p_fatal_error("NaN detected in multigrid solve",PGMG_MATH_ERROR);
		for (int kk = 0; kk < so; kk++)
			for (int jj = 0; jj < sn; jj++)
				for (int ii = 0; ii < sm; ii++) {
					u0[index(ii,jj, kk)].vel[comp] = x0[ii + jj*gm4 + kk*gmn4];
		}
	}
}

/**
 * Updates velocity derivatives at provided field to reflect
 * body force at current time.
 * \param[in] node a ptr to the cell center to consider. */
inline void fluid_3d::source_step(field *node,double fx_,double fy_,
	double fz_) {
	double tmp[3] = {fx_, fy_, fz_};
	// add each function's value to derivative
	node->add_vel_derivs(tmp);
}

/** ##################### Multigrid+FE Projection #######################*/

void fluid_3d::schedule(int s) {
	pres.mgo->schedule(s);
	if (impl) for (int i = 0; i < 3; i++) visc.mgo[i].schedule(s);
	if (godunov) mac.mgo->schedule(s);
}

/** Set up grid for finite element, increment 1 in non-periodic dimension; create multigrid object */
void fluid_3d::setup_fem() {
    if(rank==0) printf("# Creating MG solver for FEM formulation of projection method.\n");
    pres.set_grid_spacings(dx, dy, dz);
	pres.setup();
#if defined(VAR_DEN)
    pres.set_rhoinv(lrho);
#endif
    pres.set_stencils();
	mg_pres = new multigrid(pres,fem_grid,cbuf);
	pres.connect(mg_pres);
	mg_pres->compute_rats();
	mg_pres->reg[0]->clear_x_field();
}

void fluid_3d::setup_impl_visc() {
    if(rank==0) printf("# Creating MG solver for implicit viscosity. BUT WAIT, constants not set in mg_visc!\n");
    visc.set_grid_spacings(dx, dy, dz);
	mg_visc = new multigrid(visc,*grid,cbuf);
	visc.connect(mg_visc);
	mg_visc->compute_rats();
	mg_visc->reg[0]->clear_x_field();
}

void fluid_3d::setup_mac() {
    if(rank==0) printf("# Creating MG solver for MAC.\n");
    mac.set_grid_spacings(dx, dy, dz);
    mac.setup();
#if defined(VAR_DEN)
    mac.set_rhoinv(lrho);
#endif
    mac.set_stencils();
	mg_mac = new multigrid(mac,*grid,cbuf);
	mac.connect(mg_mac);
	mg_mac->compute_rats();
	mg_mac->reg[0]->clear_x_field();
}

void fluid_3d::cleanup_fem() {
	delete mg_pres;
}

void fluid_3d::cleanup_impl_visc() {
	delete mg_visc;
}

void fluid_3d::cleanup_mac() {
	delete mg_mac;
}

void fluid_3d::free() {
	fem_grid.free();
	mg_pres->free();
	if (godunov) mg_mac->free();
	if (impl) mg_visc->free();
	if (trace) tr->free();
}

void fluid_3d::compute_grad_basis(){
	if(dx!=dy || dx!=dz) {
        printf("%g %g %g\n", dx, dy, dz);
		p_fatal_error("fluid_3d::compute_grad_basis(): Spatial discretization not homogenous!", PGMG_MATH_ERROR);
	}
	double fac = 3*dx/dt;
	for (int i=0; i<8; i++){
		gradphi[i][0]=(2*(i%2)-1)*fac;
		gradphi[i][1]=(2*((i/2)%2)-1)*fac;
		gradphi[i][2]=(2*(i/4)-1)*fac;
	}
}

/**
 * Maps reference nodes to local nodes in pressure grid
 * \param[in] (ei, ej, ek) the indeces of the elements, from -1-sm, -1-sn, -1-so
 * \param[in] (ei, ej, ek) goes to sm+1, sn+1, so+1 for the last processors mp, np, op
 * \param[out] local_nodes returns the indeces of nodes in flattened pressure grid. */
int fluid_3d::map_ref_to_local(int ei, int ej, int ek, int (&local_nodes) [8]){
	// there's no ghost nodes in the multigrid geometry set up
	int psm = fem_grid.sm, psmn = psm*fem_grid.sn;
	int elem_id = ei+ej*psm+ek*psmn;
	int xedge = psm-1, yedge = fem_grid.sn-1, zedge = fem_grid.so-1;

	if (within_domain(ai+ei,aj+ej,ak+ek)) {

		for(int i=0; i<8; i++){
			local_nodes[i]= elem_id + i%2 + ((i/2)%2)*psm + (i/4)*psmn;
			if(ei+i%2 >xedge) local_nodes[i]=-1;
			if(ei<0 && i%2==0) local_nodes[i]=-1;

			if(ej+((i/2)%2)>yedge) local_nodes[i]=-1;
			if (ej<0 && ((i/2)%2==0)) local_nodes[i]=-1;

			if(ek+(i/4)>zedge) local_nodes[i]=-1;
			if (ek<0 && (i/4)==0) local_nodes[i]=-1;

		}

	}
	else for(int i=0; i<8; i++) local_nodes[i]=-1;
	return index(ei, ej, ek);
}

/**
 * Fills multigrid source term by integrating u dot gra phi,
 * where phi is the nodal basis function.
 */
void fluid_3d::fill_pp_rhs(){

	mg_pres->reg[0]->clear_r_field();

	int elem_id; // element id relative to u0, i.e. it can be negative for ghost nodes
	int xhi = fem_grid.sm;
	int yhi = fem_grid.sn;
	int zhi = fem_grid.so;
	int local_nodes[8] = {0,0,0,0, 0,0,0,0};
	double *rp = mg_pres->reg[0]->r0;
	double total = 0;

	for(int k=-1;k<zhi;k++) for (int j=-1;j<yhi;j++) for (int i=-1;i<xhi;i++){
		elem_id = map_ref_to_local(i, j, k, local_nodes);

		field *elem = u0+elem_id;
		for (int n=0;n<8;n++){
			if(local_nodes[n]>=0) {
				double val=0;
				for(int l=0;l<3;l++)
                    if(pres_update) {
                        val+=elem->dvel[l]*gradphi[n][l];
                    } else {
                        val+=elem->vel[l]*gradphi[n][l];
                    }
				total += val;
				*(rp+local_nodes[n]) += val;
			}
		}

	}

	// now also go through each side if that side isn't periodic
	field *u;
	int inds[3],g_ind,off[3],g_inds[8],ind,fn[6] = {12,14,10,16,4,22};
	double atotal[6] = {0,0,0,0,0,0};

	for (int s = 0; s < 6; s++) {

		if (neigh[fn[s]] == -1) {

			inds[faces[s].n] = s%2;
			off[0] = faces[s].u0_off % sm4;
			off[1] = (faces[s].u0_off / sm4) % sn4;
			off[2] = faces[s].u0_off / smn4;

			int start[3];
			int end[3];
			for (int i = 0; i < 3; i++) {
				start[i] = -1;
				end[i] = faces[s].len[i]+1;
			}
			start[faces[s].n] = 0;
			end[faces[s].n] = 1;

			for (int k = start[2]; k < end[2]; k++)
				for (int j = start[1]; j < end[1]; j++)
					for (int i = start[0]; i < end[0]; i++) {

				map_ref_to_local(i+off[0],j+off[1],k+off[2],g_inds);

				u = faces[s].u0 + index(i,j,k);
				for (inds[faces[s].t[0]] = 0;
					inds[faces[s].t[0]] < 2; inds[faces[s].t[0]]++) {

					for (inds[faces[s].t[1]] = 0;
						inds[faces[s].t[1]] < 2; inds[faces[s].t[1]]++) {

						ind = inds[0] + 2*inds[1] + 4*inds[2];
						g_ind = g_inds[ind];

						if (g_ind >= 0) {

							// gradphi should have correct scaling,
							// and the sign of the normal
							// vector should naturally be
							// included by choosing the local
							// indices 0-7 that correspond to
							// the side in question
							double val = -u->fvel[s][s/2]*gradphi[ind][s/2];
							atotal[s] += val;
							rp[g_ind] += val;
						} else {

							map_ref_to_local(i+off[0],j+off[1],k+off[2],g_inds);
						}
					}
				}
			}
		}

	}

	MPI_Barrier(world);

}

/**
 * Solves the projection step using Multigrid. */
int fluid_3d::pp_solve(bool verbose) {
	if(verbose) printf("Rank %d. Projection step, multigrid solve.\n",rank);

    watch.tic(1);
#if defined(VAR_DEN)
    pres.set_rhoinv(lrho);
#endif
    pres.set_stencils();
	fill_pp_rhs();
    double res;

	if (out==NULL || !pres.mgo->searching) {
		res = pres.run(mg_pres, false);
	} else {
		res = pres.run(mg_pres, true);
	}
    if(std::isnan(res)) return 1;

	extract_pp_soln();
    watch.toc(1);

	communicate<2>();

    watch.tic(1);
	subtract_pres_grad();
    watch.toc(1);

    return 0;
}

/**
 * Extract solution from Multigrid solver. */
void fluid_3d::extract_pp_soln(){

	// extract the solution from multgrid solver
	int tmp_m = mg_pres->reg[0]->mg;
	int tmp_mn = mg_pres->reg[0]->mng;
	int xhi = fem_grid.sm;
	int yhi = fem_grid.sn;
	int zhi = fem_grid.so;
	double mnoxsp = 1/(double) ((fem_grid.m)*(fem_grid.n)*(fem_grid.o));
	double l_mean=0;
	double p_mean;
	double *soln = mg_pres->reg[0]->x0;

	for (int k=0;k<zhi;k++) for(int j=0;j<yhi;j++) for (int i=0;i<xhi;i++)
		l_mean+= (*(soln+i+j*tmp_m+k*tmp_mn))*mnoxsp;

	MPI_Allreduce(&l_mean,&p_mean,1,MPI_DOUBLE,MPI_SUM,fem_grid.cart);

	for (int k=0;k<zhi;k++) for(int j=0;j<yhi;j++) for (int i=0;i<xhi;i++){
        // if pres_update == true, only pressure update is computed
		if(pres_update) {
            (u0 + index(i,j,k))->p += ( *(soln+i + j*tmp_m +k*tmp_mn) -p_mean);
            (u0 + index(i,j,k))->dp = ( *(soln+i + j*tmp_m +k*tmp_mn) -p_mean);
        } else {
            (u0 + index(i,j,k))->p = ( *(soln+i + j*tmp_m +k*tmp_mn) -p_mean);
        }
    }
}

/**
 * Subtracts the pressure gradient from intermdediate velocity.
 */
void fluid_3d::subtract_pres_grad(){

    if(!pres_update) {
        for(int k=0;k<so;k++) for(int j=0;j<sn;j++) for(int i=0;i<sm;i++){
            int elem_id = index(i,j,k);
            field *f = u0+elem_id;
            double npg[3] = {0,0,0};
            neg_pres_grad(f,npg);
            for(int p=0;p<3;p++){
#if defined(VAR_DEN)
                npg[p]*=dt/lrho[elem_id+G0];
#else
                npg[p]*=dt;
#endif
                f->vel[p] += npg[p];
            }
        }
    } else {
        for(int k=0;k<so;k++) for(int j=0;j<sn;j++) for(int i=0;i<sm;i++){
            int elem_id = index(i,j,k);
            field *f = u0+elem_id;
            double npg[3] = {0,0,0};
            neg_pres_dgrad(f,npg);
            for(int p=0;p<3;p++){
#if defined(VAR_DEN)
                npg[p]*=dt/lrho[elem_id+G0];
#else
                npg[p]*=dt;
#endif
                f->vel[p] += npg[p];
                f->dvel[p] += npg[p];
            }

        }

    }
}

/**
 * A function that takes in an exact solution and computes
 * the error locally, then sends it to the master node.
 * \param[in] ref_soln pointer to the reference solution.
 * \param[out] e the array to write error in (u,v,w,p).
 */
double fluid_3d::error(field *ref_soln, double (&e)[4]) {

	// start with no errors, allocate for global sum
	for (int i = 0; i < 4; i++) e[i] = 0;
	double global_e[4] = {};

	for(int k = 0; k < fem_grid.so; k++) {
		for(int j = 0; j < fem_grid.sn; j++) {
			for(int i = 0; i < fem_grid.sm; i++) {

				// get current node
				int elem_id = index(i,j,k);
				field &node = u0[elem_id];
				field &soln = ref_soln[elem_id];
				// get vels error at this node
				if(i<sm && j<sn && i < so)
					node.add_error(soln.vel,e);
				else
					node.add_p_error(soln.p,e);
			}
		}
	}

	// sum errors and send to master processor
	MPI_Reduce(&e,&global_e,4,MPI_DOUBLE,MPI_SUM,0,grid->cart);
	if (rank == 0) for (int i = 0; i < 4; i++) e[i] = global_e[i];

	// normalize and add for return value
	double tot_err = 0;
	for (int i = 0; i < 4; i++) {

		// if master, normalize by all
		if (rank == 0) {
			e[i] /= (i<3)?(m*n*o):(fem_grid.m*fem_grid.n*fem_grid.o);
		} else {
			e[i] /= (i<3)?(sm*sn*so):(fem_grid.sm*fem_grid.sn*fem_grid.so);
		}
		// add to return scalar
		tot_err += e[i];
	}
	return tot_err;
}

/**
 * A function that takes in an exact solution and computes
 * the error locally, then sends it to the master node.
 * \param[in] (euf, evf, ewf, epf) exact analytical solutions in (u,v,w,p).
 * \param[in] e the array to write errors in each component to. */
double fluid_3d::error(double (&e)[4]) {

	// march through real nodes

	// start with no errors, allocate for global sum
	for (int i = 0; i < 4; i++) e[i] = 0;
	double global_e[4] = {0,0,0,0};

	for(int k = 0; k < fem_grid.so; k++) {
		double xx,yy,zz;
		zz = lz0[k];
		for(int j = 0; j < fem_grid.sn; j++) {
			yy = ly0[j];
			for(int i = 0; i < fem_grid.sm; i++) {
				xx = lx0[i];

				// get current node
				field &node = u0[index(i, j, k)];
				double uu,vv,ww,pp;
				// assumes isotropic
				mgmt->exact(xx,yy,zz,time,dx,uu,vv,ww,pp);
				double ex_soln[4] = {uu,vv,ww,pp};

				// get error at this node
				if (i < sm && j < sn && k < so) {
					node.add_error(ex_soln,e);
				} else {
					node.add_p_error(ex_soln[3],e);
				}
			}
		}
	}

	// sum errors and send to master processor
	MPI_Reduce(&e,&global_e,4,MPI_DOUBLE,MPI_SUM,0,grid->cart);
	if (rank == 0) for (int i = 0; i < 4; i++) e[i] = global_e[i];

	// normalize and add for return value
	double tot_err = 0;
	for (int i = 0; i < 4; i++) {

		// if master, normalize by all
		if (rank == 0) {
			e[i] /= (i<3)?(m*n*o):(fem_grid.m*fem_grid.n*fem_grid.o);
		} else {
			e[i] /= (i<3)?(sm*sn*so):(fem_grid.sm*fem_grid.sn*fem_grid.so);
		}

		// add to return scalar
		tot_err += e[i];
	}
	return tot_err;
}

/**
 * Instantiates finite difference stencils
 */
void fluid_3d::setup_stencils() {

	// 1D stencils
	// first-derivative prefix is just char.time/char.length
	double p1x = dt/dx;
	double p1y = dt/dy;
	double p1z = dt/dz;

	// second-derivative is dt/(char.length)^2
	double p2x = dt/(dx*dx);
	double p2y = dt/(dy*dy);
	double p2z = dt/(dz*dz);

	// Laplacian, because of viscous term, is 2nd deriv
	// prefix with additional factor of Reynolds num and
	// rho in denominator (as in inkjet paper)
	double slx = mgmt->fm.mu*dt/(dx*dx);
	double sly = mgmt->fm.mu*dt/(dy*dy);
	double slz = mgmt->fm.mu*dt/(dz*dz);

	// if viscosity is time-stepped implicitly, introduce
	// factor of 1/2 to stencil
	if (impl) {
		slx *= 0.5; sly *= 0.5; slz *= 0.5;
	}

	// instantiate each 1D set of stencils

	// second-order, first-deriv
	int order = 2;
	cart kind[3] = {X,Y,Z};
	type dir[3] = {forward,centered,backward};

	// outer loop is cartesian direction
	for (int i = 0; i < 3; i++) {

		// inner loop is forward/centered/backward diff
		for (int j = 0; j < 3; j++) {

			// do all of those for 1st derivs (for ENO2)
			sten1[j][i] = stencil(order,1,dir[j],kind[i],
				sm4,sn4,p1x,p1y,p1z);

			// for Godunov Scheme
			if(j==1) godu_sten1[j][i] = stencil(2,1,dir[j], kind[i],
                sm4,sn4,dxsp,dysp,dzsp);
			else godu_sten1[j][i] = stencil(1,1,dir[j], kind[i],
				sm4,sn4,dxsp,dysp,dzsp);

			if (j != 1) {
				double fac;
				if (j == 0) {
					fac = impl?0.5:1;
				} else {
					fac = 1/dt;
				}
				sten_stress[j/2][i] = stencil(1,1,dir[j],
					kind[i],sm4,sn4,fac*p1x,fac*p1y,fac*p1z);
			}

			sten_grad[j][i] = stencil((j==i)?1:2,1,(j==i)?backward:centered,
				kind[i],sm4,sn4,p1x/dt,p1y/dt,p1z/dt);
			if (i != j) sten_grad[j][i].avg_shift(-(((j==0)?1:((j==1)?sm4:smn4))));
		}

		// do all cartesian dir's for 2nd deriv, but
		// just centered diff (for ENO2)
		sten2[i] = stencil(order,2,centered,kind[i],
			sm4,sn4,p2x,p2y,p2z);
	}

	// second-order Laplacian (for viscous term)
	order = 2;
	lapl = stencil(order,2,centered,Laplacian,sm4,sn4,slx,sly,slz);

	// if doing viscosity implicitly, set that stencil too
    visc.setup_const(slx,sly,slz);
	visc.setup();
}

/**
 * Function to instantiate bc structures using double/int (value/type)
 * description of boundary conditions from constructor
 * \param[in] bct the array that describes the type of boundary condition at each face.
 * \param[in] bcv the array that prescribes the value of velocities (derivatives) at each face. */
void fluid_3d::setup_faces() {

	int bct[3];
	double bcv[3];
    double low_bound [3] = {mgmt->ax, mgmt->ay, mgmt->az};
    double hi_bound [3] = {mgmt->bx, mgmt->by, mgmt->bz};
	double h;
	for (int i = 0; i < 6; i++) {

		for (int j = 0; j < 3; j++) {
			bct[j] = spars->bct[j][i];
			bcv[j] = spars->bcv[j][i];
            if(bct[j] == BC::DIRICHLET && bcv[j] > mgmt->vmax) mgmt->vmax = bcv[j];
		}

		switch (i/2) {
		case 0: h = dx; break;
		case 1: h = dy; break;
		case 2: h = dz; break;
		}

		faces[i].setup(i,grid,bct,bcv,h,u0,lx0,ly0,lz0,low_bound,hi_bound);
	}
}
