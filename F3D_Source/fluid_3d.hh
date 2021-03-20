#ifndef FLUID_3D_HH
#define FLUID_3D_HH
#include <cstdio>
#include <vector>
#include <cmath>
#include <algorithm>
#include <unistd.h>
#ifdef _OPENMP
#include "omp.h"
#endif

#include "defs.hh"
#include "mat3.hh"
#include "custom_structs.hh"
#include "geometry.hh"
#include "common.hh"
#include "buffer.hh"
#include "stencil.hh"
#include "mg_probs.hh"
#include "multigrid.hh"
#include "fields.hh"
#include "comm_info.hh"
#include "timer.hh"
#include "tracers.hh"
#include "extrap_mmap.hh"
#include "extrap_helper.hh"

#include "write_params.hh"
#include "sim_params.hh"
#include "sim_manager.hh"

using namespace F3D_CONST;
using namespace FIELD_CONST;

class tracers;

class fluid_3d {

	public:

	static const int packlens[12];
	static const int domain_face[6];

	/** pointer to simulation info */
	sim_manager *mgmt;
    sim_params *spars;

	// ##### ALL INTEGERS #####
	// geometry info
	/** determines implicit/explicit timestepping of viscosity */
	const bool impl;
	/** to use Godunov scheme for advection update */
	const bool godunov;
    /** to compute pressure or pressure update */
    static const bool pres_update;
	/** to do tangential stability calculation in upwinding or not */
	const bool tang;
	/** whether to use tracer particles */
	bool trace;
	/** The rank of the processor */
	const int rank;
	/** global number of grid points in x,y,z direction */
	const int m, n, o;
	/** The local number of grid points in x,y,z direction*/
	const int sm, sn, so;
	/** The local number of grid points in x,y,z direction with ghost nodes*/
	const int sm4, sn4, so4;
	/** The number of grid points in x-y plane (including ghosts) */
	const int smn4;
	/** The index of the first non-ghost node */
	const int G0;
	/** The indices (x, y, z) of local lower point (inclusive) */
	const int ai, aj, ak;
	/** The indices (x, y, z) of local upper point (exclusive) */
	const int bi, bj, bk;
    /** The number of extrapolation errors */
    int extp_err;
	/** The dimension to export when dump() is called */
	int dd;
	/** Which slice to export when dump() is called */
	int ds;
	int nt;
	/** current index of layer, dynamically used in extrapolations */
	int n_layers;
	/** maximum number of layers */
	int max_layers;
	/** maximum radius of search for extrapolation */
	int max_extrap_rs;

	// ##### ALL ENUMS #####

	// ##### ALL DOUBLE #####
	/** Spatial step size in x,y,z direction */
	const double dx,dy,dz;
	/** The inverse of spatial step size in x,y,z direction */
	const double dxsp, dysp, dzsp;
	/** time discretization */
	double dt;
	/** current time */
	double time;
	/** basis of integrated basis functions for FEM solve */
	double gradphi[8][3];

	// ##### HIGH-LEVEL DATATYPES ####
	/** first-derivative stencils */
	stencil sten1[3][3]; // dim 1 f/c/b, dim 2 x/y/z
	/** second-derivative stencils */
	stencil sten2[3];	// dim x/y/z
	/** laplacian stencil */
	stencil lapl;
	/** Godunov frist-derivative stencils */
	stencil godu_sten1[3][3]; // dim1 f/c/b, dim2 x/y/z
	/** forward/backward 1st order 1st deriv stencils for stress */
	stencil sten_stress[2][3];
	/** stencils for reference map gradient */
	stencil sten_grad[3][3];
	/** communication buffer for multigrid */
	comm_buffer cbuf;
	/** pointer to geometry for fem solve in projection step */
	geometry fem_grid;
	/** multimap structure for reference map extrapolations */
	extrap_mmap temp_extraps;
	extrap_mmap extraps;
    /** multigrid solver for implicit viscosity */
    multigrid *mg_visc;
	/** multigrid solver for solenoidal projection */
	multigrid *mg_pres;
	/** multigrid solver for mac projection */
	multigrid *mg_mac;
	/** problem for mg solves */
#if defined(VAR_DEN)
	simp_visc visc;
	var_den_fem pres;
	var_den_mac mac;
#else
	simp_visc visc;
	simp_fem pres;
	simp_mac mac;
#endif

	/** array of boundary condition objects */
	face_gm faces[6];
	/** lookup table for neighbor communication information */
	comm_info comm_table[27];
	comm_info corner_table[27];
	/** tracers */
	tracers *tr;

	// ##### ALL POINTERS #####

	/** The pointer to the parent geometry class*/
	geometry *grid;
	/** array of fluid data points */
	field *u_mem;
    /** pointer to the first non-ghost fluid node */
    field *u0;

#if defined(VAR_DEN)
    /** density */
    double *lrho;
#endif

#if defined(DEBUG)
    /** determinant of F, J*/
    double *lJ;
#endif

	/** array of solid data points */
	ref_map *rm_mem;
	/** pointer to the first non-ghost solid node */
	ref_map *rm0;
	/** list of neighbors */
	int *neigh;
	/** request objects for MPI calls */
	MPI_Request *reqs;
	/** status objects for MPI calls */
	MPI_Status *stats;
	/** The arrays of global x, y, z components of grid points physical coordinates */
	double *gx, *gy, *gz;
	/** The arrays of local x, y, z components of grid points physical coordinates */
	double *lx, *ly, *lz;
	/** The pointer to first non-ghost node x,y,z coordinate */
	double *lx0, *ly0, *lz0;
	/** output string for sharing status with user */
	char *out;
	/** A timer object to keep track of old and new extrapolation routine time usage */
	timer watch;
	extrap_helper<3> expper;

	// ##### class functions #####

	/** Constructors and destructor*/
	fluid_3d(geometry *gm, sim_manager *mgmt_, sim_params *spars_);
	~fluid_3d();

	// setup functions
	void setup_communication();
	void setup_stencils();
	void setup_cart();
	void setup_faces();
	void setup_tables();
	void setup_fem();
	void setup_impl_visc();
	void setup_mac();

	// cleanup functions
	void cleanup_communication();
	void cleanup_impl_visc();
	void cleanup_fem();
	void cleanup_mac();
	void free();

	// boundary conditions
	void set_bc_cc(bool verbose);
	void set_bc_cf(bool verbose);

	// communication on f3d grid
	template<int flags>
	void shift_refmap(int type,double dummy);
	template<int flags>
	int copy_to_buffer(int tag,double *(&b),int lint,comm_info *ct);
	template<int flags>
	void copy_from_buffer(int tag,double *(&b),comm_info *ct);
	template<int flags>
	void communicate(int lint=0);

	// communication on extrap_mmap object
	void shift_extraps(int type,int lid=-1);
	int copy_layer_to_buffer(int tag,unsigned int c,double *(&b),int type);
	void copy_extraps_from_buffer(int tag,double *(&b),int type);
	int copy_extraps_to_buffer(int tag,double *(&b),int type);
	void communicate_a_layer(unsigned int c,int type);
	void communicate_extrap_fx();

    // fill boundaries
    inline void fill_boundary_cc(bool verbose){
        /** Fill boundary velocities */
        if(verbose && rank==0) printf("Rank %d. Communicate cell-centered velocities.\n", rank);

        // non-periodic boundary conditions
        // then periodic is taken care of
        set_bc_cc(verbose);
    }

    // fill boundaries
    inline void fill_boundary_cf(bool verbose){
        /** Fill boundary velocities */
        if(verbose && rank==0) printf("Rank %d. Communicate face velocities.\n", rank);

        // non-periodic boundary conditions
        set_bc_cf(verbose);
    }

	/* Utilities */
    void check_nans(const int id);
    void sanity_check();
    double check_steady_solution();
	bool check_null(ref_map *r0, ref_map *r1, ref_map *r2, ref_map *r3,
            ref_map *r4, ref_map *r5, ref_map *r6, ref_map *r7);
	int index(int i, int j, int k);
	void unpack_index(int ind, int &i, int &j, int &k);
    ref_map * get_refmap(int ind, unsigned int oid);
    ref_map * get_refmap_prev(int ind, unsigned int oid);
	ref_map * get_refmap(int ind, unsigned int oid, const int (&pos)[3]);
	ref_map * get_refmap_prev(int ind, unsigned int oid, const int (&pos)[3]);
	void solid_fractions(int ind, std::vector<solid_frac> &s_fracs);
	int nface(int tag, int& ti, int& tj, int& tk);
	void record_magnitude();

	/** Problem set up */
	int initialize();
	void init_iter(int init_err);
    int initialize_from_chk_point(const char * chk_dirname);
    void init_from_slice(write_params wp, double * g_val);
	void init_fluid();
	void init_pres();
	void set_pres();
	void init_refmap();
	void init_extrapolate(bool verbose);

	inline void set_omp(int nthreads) {
#ifdef _OPENMP
		omp_set_num_threads(nthreads);
#endif
	}

	/** Multigrid + FE projection */

	/** returns true if element is within domain */
	inline bool within_domain(int i,int j,int k) {
		return  (grid->x_prd || (i>=0 && i<=m-1)) &&
				(grid->y_prd || (j>=0 && j<=n-1)) &&
				(grid->z_prd || (k>=0 && k<=o-1));
	}
	/** returns true if element is within domain */
	inline bool within_local_domain(int i,int j,int k) {
        //printf("within local domain %d %d %d = %d\n",
        //      i,j,k, (i>=0 && i<=sm4-1) &&
        //      (j>=0 && j<=sn4-1) &&
        //      (k>=0 && k<=so4-1) );
		return  (i>=0 && i<=sm4-1) &&
				(j>=0 && j<=sn4-1) &&
				(k>=0 && k<=so4-1);
	}

	void compute_grad_basis();
	int map_ref_to_local(int ei, int ej, int ek, int (& local_nodes) [8]);
	void fill_pp_rhs();
	int pp_solve(bool verbose);
	void extract_pp_soln();
	void subtract_pres_grad();
#if defined(VAR_DEN)
    void update_rho();
#endif

	// time stepping
	void schedule(int s);
	void calc_derivs();

	// new stress functions
	template<lower_faces F>
	int collision_stress(int eid, double (&fluid_s)[3]);
	template<lower_faces F>
	void solid_stress(int eid, double (&solid_s)[3], double & sfrac);
	template<lower_faces F>
	void fluid_stress(int eid, double (&fluid_s)[3]);
	void velocity_grad(lower_faces F,int eid,matrix &grad_v);

	void compute_stress(bool verbose);

	// godunov scheme functions
	void godunov_set_refmap(double &ul, double &vl, double &wl, double &ur, double &vr, double &wr, double &mnorml, double &mvl, double &mwl, double &mnormr, double &mvr, double &mwr);
	void godunov_set_velocity(double &ul, double &vl, double &wl, double &ur, double &vr, double &wr);

#if defined(TANG)
	void godunov_set_tang_refmap(double &ul, double &vl, double &wl, double &ur, double &vr, double &wr, double &mnorml, double &mvl, double &mwl, double &mnormr, double &mvr, double &mwr);
	void godunov_set_tang_velocity(double &ul, double &vl, double &wl, double &ur, double &vr, double &wr);
#endif

	void mac_solve(double cdt, bool verbose);
    void fill_mac_rhs(double cdt);
    void extract_mac_soln();
    void subtract_q_grad(double cdt);
	void neg_pres_grad(field *node, double (&pres_acc)[3]);
	void neg_pres_dgrad(field *node, double (&pres_acc)[3]);
	double mono_diff(double u0, double u1, double u2, double u3, double u4);
	double del_lim(double u0, double u1, double u2);
	double del_f(double u0, double u1, double u2);
	inline double centered_diff(double u0, double u2){return 0.5*(u2-u0);}
	inline double min(double x, double y){ return (x<y)?x:y; }
	inline double sign(double x){return (x>=0)?1:-1; }

	// new set of extrapolation functions, under construction
	void extrapolate(bool verbose);
	void first_layer(unsigned int oid);
	void extrap_a_layer(unsigned int oid, bool verbose);
	void gather_and_list(unsigned int oid, bool verbose);
	void transfer(unsigned int oid, int ind);
	int rule();
	void redraw_boundary(bool verbose);

	// Functions for stepping forward
    // debug flag -1 is just verbose
    // but non zero values would invoke other functions, e.g. sanity check
	int step_forward(int debug);
	void normal_derivatives_velocity(const double cdt, bool verbose);
	void normal_derivatives_refmap(const double cdt, bool verbose);
	void godunov_upwinding_set(bool tang_stab, bool verbose);
#if defined(TANG)
	void compute_tang_derivatives(bool verbose);
#endif
	void compute_half_time_edge_velocities(const double cdt, bool verbose);
	void compute_half_time_edge_refmap(const double cdt, bool verbose);
	int update_reference_map(const double cdt, bool verbose);
	void revert_reference_map_update();
	void update_reference_map_full(bool verbose);
	void compute_ustar(const double cdt, bool verbose);
	void acceleration(int ijk, double myx, double myy, double myz, double (&acc) [3], bool pres, bool gdn_ex, bool verbose=false);

	// older time stepping functions
	void no_adv(bool visc=true, bool force=true);
	void expl_timestep();
	void impl_timestep();
	void step();
	void step(int n);
	void viscous_step(field *node);
	void eno_step(field *node);
	void source_step(field *node,double fx_,double fy_,double fz_);

	// data output (slices, contour plotting and errors)
	void alloc_tracers();
	void set_dump(int dim,int pnt);
	void dump(int num, char *dir);
	void display_stats();
	void copy_slice_from_buf(int **sender_list,int co,double *g_val);
	void write_slice(write_params params,const char* filename);
	void write_chk_pt(const int step_num, const int chk_num, const char* filename);
	void copy_chk_slice_to_buf(double *g_val, double *comm_buf, int *info);
	void copy_chk_to_field(write_params wp);
	void output_phi(int obj_num,char * fname,int type);
	void slice(write_params params,double *g_val);
	int copy_slice_to_buf(write_params params,int point);
	void save_matrix(FILE *fh,double *u_global,int d1,int d2);
	void save_text(FILE *fh,double *u_global,int d1,int d2,double *a1,double *a2);
	void save_gnuplot(FILE *fh,double *u_global,int d1,int d2,double *a1,double *a2);
	int min_phi_id(int eid);
	double min_phi(int eid);
	double phi(int eid,int cnum);

	// Contour plotting
	void output_contours(const char* filename);
	void output_contour(FILE *outf,int cnum);
	void assemble_contour();
	void gather_sizes();
	void setup_output_dimensions();
	void send_contour_data(std::vector<double> &pv,std::vector<int> &q,std::vector<int> &q2,MPI_Request *sreq,int cnum);
	void edge_detect(int ijk,int dir,std::vector<double> &pv,std::vector<int> &q,int cnum);
	void box_detect(int ijk,std::vector<int> &q2,unsigned int &ttri,int cnum);
	void normalize(double &nx,double &ny,double &nz);
	void grad_field(int ijk,double &nx,double &ny,double &nz,int cnum);

	// debugging tools
    int check_extrapolated_region(int &ll);
    int check_wall_extrapolations();
	void debug_dump(const char *filename=NULL);
	// TODO luxury debug item, save the current state of the sim
	// right before it fails, then rerun simulation with this
	// as initial condition
	//void save_image(const char *filename, int snap);
	void check_symmetry();
	double div_u(double_int &gmx_ext, double_int &gmn_ext);
	void attach_debug(int proc);
    double vel_mod(int norm);
    void write_all_Js(char* fn);
    void total_momentum(double (&mv_ext)[3]);
    double total_energy(double &pot_energy, double &kin_energy, double &elas_energy,double &power);
    double total_dev_in_detF(double & max_dev);
    double avg_detF();
    void compute_centroid(double &centx, double &centy, double &centz, double &vol);
    void compute_v_centroid();

#if defined(DEBUG)
    double_int max_Gaussian_curvature(int &sign);
    double_int max_solid_stress(int type);
    double_int max_coll_stress(int type);
    double_int max_coll_traction();
    double_int max_grad_phi_diff();
#endif

	// error analysis
	double error(field * ref_soln, double (&e)[4]);
	double error(double (&e)[4]);

	void set_dt(double dt_) {
		dt = dt_;
		spars->set_dt(dt_);
	}

	// Tables for contour plotting
	private:
	int sizes[3];
	int *osm,*osn,*oso,*max_sizes;
	int *cu_m,*cu_n,*cu_o;
	char* ptri_poly[257];
	static const char n_poly[256];
	static const char tri_poly[2460];
};

#endif
