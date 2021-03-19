#ifndef TRACERS_HH
#define TRACERS_HH

#include <cstdio>
#include <cmath>
#include <vector>
#include <mpi.h>
#include <gsl/gsl_rng.h>

#include "common.hh"
#include "geometry.hh"
#include "fields.hh"
#include "sim_params.hh"
#include "fluid_3d.hh"

class fluid_3d;

/**
 * struct representing a single tracer particle
 *
 * the tracer's position is represented by a dimensionless vector with
 * components vc \in [-0.5,0.5] corresponding to its position in the
 * interpolated volume, i.e. the central cube of space within a 4x4x4
 * cube of nodes. (0,0,0) corresponds to the center of the cube, and
 * the distance is normalized by the grid node spacing
 *
 * the tracer's velocity at a position is calculated by constructing a
 * tricubic interpolant from the velocities at the nodes - a precomputed
 * Vandermonde matrix is defined in tracers.cc
 */
struct tracer {

	/** id number (currently not used) */
	int id;
	/** bottom 4x4x4 left corner coordinate (i,j,k) in local processor inds */
	int p_crd[3];
	/** position w/in interpolated volume (norm'd by spacing -0.5 < rp < 0.5) */
	double rp0[3];
	/** estimated position update for improved Euler method */
	double rpf[3];
	/** velocity at start of time step */
	double v0[3];
	/** velocity at end of time step */
	double vf[3];

	// constructors
	tracer(int i) : id(i) {}
	tracer() {}

	/** print helper for debugging */
	void print_self(){
		printf("Tracer %d in box (%d %d %d), rp0(%g %g %g)\n",
			id,p_crd[0],p_crd[1],p_crd[2],rp0[0],rp0[1],rp0[2]);
	}
};

/**
 * class representing a set of tracer particles for use
 * with the fluid_3d class
 *
 * each processer owns an STL vector of tracers which are communicated
 * as needed using the custom defined MPI_TRACER datatype
 *
 *
 *
 */
class tracers {

public:
	int rank;
	/** global number of tracers FIXME currently an estimate */
	int num;
	/** memory stride lengths */
	int del[3];
	/** number of grid points (with ghosts) in each dimension */
	int lim[3];
	/** number of tracers to send to each of 27 neighbors */
	int slen[27];
	/** number of tracers to receive from each of 27 neighbors */
	int rlen[27];
	/** length of time step */
	double &dt;
	/** distance between grid points in each dimension */
	double h[3];

	/** inverted vandermonde matrix for fitting 4^3 cube to tricubic */
	static double Ivmonde[4096];
	/** cartesian coordinate of lowest grid point
	 *  (incl. ghosts) on this processor */
	double low[3];
	/** cartesian coordinate of lowest global grid point
	 *  i.e. NOT including ghost nodes */
	double glow[3];
	/** vector of tracers on this processor */
	std::vector<tracer> tr;
	/** defined tracer datatype for MPI communication */
	MPI_Datatype MPI_TRACER;
	/** for MPI communication */
	MPI_Request *reqs;
	/** for MPI communication */
	MPI_Status *stat;
	/** random number generator */
	gsl_rng *rng;
	/** pointer to array of field values */
	field *u_mem;
	/** pointer to parent geometry */
	geometry *grid;
	/** pointer to simulation parameters */
	sim_params *sim;
	/** communication buffer */
	tracer *cbuf;
	/** neighbor array for communication */
	int *neigh;

	// constructors / destructor
	tracers(fluid_3d& f3d);
	~tracers();

	// setup and freeing
	/** set up datatypes for MPI communications */
	void setup_comms();
	/** frees the memory used for the MPI tracers datatype and GSL RNG */
	void free() {MPI_Type_free(&MPI_TRACER); gsl_rng_free(rng);}

	// functions in .cc file
	void print(const char *fname,bool binary=true);
	void init_from_chkpt_files(const int ntracers, const char *fname, bool binary=true);
	void gather();
	void communicate();
	int proc_adjust(tracer &t,int st);
	void fill_buffer();
	void fill_buffer(std::vector<int> &targets);
	void empty_buffer(int st);
	void interp_vels(const tracer &t,const double (&r)[3],double (&v)[3]);

	/** returns the full ind-th cartesian coordinate of tracer t,
	* based on the initial location r0 */
	inline double loc0(const tracer &t,int ind) {
		return low[ind] + (1.5 + t.rp0[ind] + t.p_crd[ind])*h[ind];
	};
	/** returns the full ind-th cartesian coordinate of tracer t,
	* based on the half-way position rf */
	inline double locf(const tracer &t,int ind) {
		return low[ind] + (1.5 + t.rpf[ind] + t.p_crd[ind])*h[ind];
	};
	/** same as loc0, but identifies the tracer by its index
	* on this processor tr_ind */
	inline double r0(int tr_ind,int ind) {
		return loc0(tr[tr_ind],ind);
	};
	/** same as locf, but identifies the tracer by its index
	* on this processor tr_ind */
	inline double rf(int tr_ind,int ind) {
		return locf(tr[tr_ind],ind);
	};
	/** sets the "local" position info a tracer at position (x,y,z).
	* returns false if that position is not on this processor and
	* true if it is */
	inline bool set(tracer &t,double x,double y,double z, int id_) {
		double r[3] = {x,y,z};
		for (int f = 0; f < 3; f++) {
			double t_ind = (r[f] - low[f])/h[f] - 1;
			t.p_crd[f] = (int) t_ind;
			t.rp0[f] = t_ind - t.p_crd[f] - 0.5;
			if (out_of_range(f,t,1)) return false;
		}
		t.id = id_;
		return true;
	};
	/** adds a tracer to the tracers class at position (x,y,z) */
	inline void add(double x,double y,double z) {
		tracer t;
		if (set(t,x,y,z, num++)) {
			tr.push_back(t);
			if (t.p_crd[0] > 200) printf("!!!!!!\n");
//			t.print_self();
		}
	};
	/** returns true if the tracer t is out of the current processor's
	* range by being too low in the fth cartesian direction */
	inline bool too_low(int f,tracer &t,int st) {
		return (t.p_crd[f] < 0) || (t.p_crd[f]==0 && ((st==0)?(t.rpf[f]):(t.rp0[f])) < 0);
	};
	/** returns true if the tracer t is out of the current processor's
	* range by being too high in the fth cartesian direction */
	inline bool too_high(int f,tracer &t,int st) {
		return (t.p_crd[f] > (lim[f]-4)) ||
			(t.p_crd[f]==(lim[f]-4) && ((st==0)?t.rpf[f]:t.rp0[f]) >= 0);
	};
	/** returns true if the tracer t is out of the current processor's
	* range in the fth cartesian direction */
	inline bool out_of_range(int f,tracer &t,int st) {
		return too_low(f,t,st) || too_high(f,t,st);
	};
	/** returns a corner to the field node which is the front corner of
	* this tracer's local 4x4x4 cube of grid points */
	inline field* corner(tracer &t) {
		return u_mem + t.p_crd[0]*del[0]
			+ t.p_crd[1]*del[1] + t.p_crd[2]*del[2];
	};
	/** predicts the tracers position after a time step of dt using
	* forward Euler, and stores the result in t.rpf */
	inline void predict(tracer &t) {

		interp_vels(t,t.rp0,t.v0);



		for (int f = 0; f < 3; f++) {

			double dr = dx(t,f,0)/h[f];
			t.rpf[f] = t.rp0[f] + dr;
		}


//		if (rank==0) printf("pred: (%g,%g,%g) [v=(%g,%g,%g)] -> (%g,%g,%g)\n",t.rp0[0],t.rp0[1],t.rp0[2],
//				t.v0[0],t.v0[1],t.v0[2],
//				t.rpf[0],t.rpf[1],t.rpf[2]);

	};
	/** runs the predictors step on all tracers */
	inline void predict() {
		/*
		puts("bb\n");
		for (tracer t:tr) {
			if (t.p_crd[0] > 200) {
				printf("%d here (predict start)!\n",rank);
				MPI_Abort(MPI_COMM_WORLD,0);
			}
		}
		puts("b\n");
		*/
		std::vector<int> targets;
		for (unsigned i = 0; i < tr.size(); i++) {
			predict(tr[i]);
			/*
		for (tracer t:tr) {
			if (t.p_crd[0] > 200) {
				printf("%d here (before pa)!\n",rank);
				MPI_Abort(MPI_COMM_WORLD,0);
			}
		}
			*/
			int tmp = proc_adjust(tr[i],0);
			/*
		for (tracer t:tr) {
			if (t.p_crd[0] > 200) {
				printf("%d here (after pa)!\n",rank);
				MPI_Abort(MPI_COMM_WORLD,0);
			}
		}
			*/
			targets.push_back(tmp);
		}
		fill_buffer(targets);
		communicate();
		empty_buffer(0);

		/*
		puts("ccc\n");
		for (tracer t:tr) {
			if (t.p_crd[0] > 200) {
				printf("%d here (predict end)!\n",rank);
				MPI_Abort(MPI_COMM_WORLD,0);
			}
		}
		*/
	};
	/** runs the improved euler routine on tracer t, using the
	* tracer's initail position and the predicted position from predict() */
	inline void imp_euler(tracer &t) {

		/*
		for (tracer t:tr) {
			if (t.p_crd[0] > 200) {
				printf("here (imp_euler start)!\n");
				MPI_Abort(MPI_COMM_WORLD,0);
			}
		}
		*/


		interp_vels(t ,t.rpf,t.vf);
		for (int f = 0; f < 3; f++) {

			// get first update
			double dr0 = t.rpf[f]-t.rp0[f];
			double drf = dx(t,f,1)/h[f];

			t.rp0[f] += 0.5*(dr0+drf);
		}
	};
	/** runs an improved euler update on all tracers */
	inline void advect() {
//		puts("c\n");
		std::vector<int> targets;
		for (unsigned i = 0; i < tr.size(); i++) {
			imp_euler(tr[i]);
			targets.push_back(proc_adjust(tr[i],1));
		}
		fill_buffer(targets);
		communicate();
		empty_buffer(1);

		/*
		for (tracer t:tr) {
			if (t.p_crd[0] > 200) {
				printf("here (advect end)!\n");
				MPI_Abort(MPI_COMM_WORLD,0);
			}
		}
		*/
	};

	private:

	// debugging / output
	void print_vels(const tracer&,int,int);

	// helper functions:

	/** returns whether on first proc in dimension d */
	bool first_proc(int d) {
		switch (d) {
			case 0: return grid->ip==0;
			case 1: return grid->jp==0;
			case 2: return grid->kp==0;
		}
		return false;
	}
	/** returns whether on last proc in dimension d */
	bool last_proc(int d) {
		switch (d) {
			case 0: return grid->ip==grid->mp-1;
			case 1: return grid->jp==grid->np-1;
			case 2: return grid->kp==grid->op-1;
		}
		return false;
	}
	/** returns whether tracer belongs to first node in dim d */
	bool first_node(tracer &t,int d) {return first_proc(d) && t.p_crd[d]==0;}
	/** returns whether tracer belongs to last node in dim d */
	bool last_node(tracer &t,int d) {return last_proc(d) && t.p_crd[d]==(lim[d]-4);}


	double dx(tracer &t,int f,int st) {

		double v = st==0?t.v0[f]:t.vf[f];
		double update = dt*v;

		// check if this direction has a homog.
		// Dirichlet boundary condition on

		// negative vel: look at "left" wall
		// positive vel: look at "right" wall
		int lr = (v<0)?0:1;
		int    bct = sim->bct[f][2*f+lr]; // type
		double bcv = sim->bcv[f][2*f+lr]; // value

		bool homog_diri = bct==1 && bcv==0;

		// no need to check if not homogeneous dirichlet
		if (!homog_diri) return update;

		// if it is homogeneous dirichlet, see how close
		// we are to global domain boundary
		//
		// this is currently doing unneccessary checking
		// since "landlocked" procs have no chance of hosting
		// tracers near global boundaries, but that's
		// presumably not a huge waste of resources

		// lower/upper boundary cartesian coord
		double lo = glow[f]-0.5*h[f];
		int n_pts[3] = {grid->m,grid->n,grid->o};
		double hi = lo + n_pts[f]*h[f];


		// cartesian coordinate
		double xx = st==0?loc0(t,f):locf(t,f);

		const double min_sep(0.01*(hi-lo));

//		printf("%g %g %g\n",lo,xx,hi);

		// sep represents the distance from the wall
		//   (which wall? the one you're heading towards)
		// vn represents magnitude of normal velocity
		double sep,vn;
		if (v==0) {

			double v1 = fabs(xx-lo);
			double v2 = fabs(hi-xx);

			sep = v2<v1?(hi-xx):(xx-lo);

		} else sep = (v<0)?(xx-lo):(hi-xx);
		vn = fabs(v);

		// if separation is zero, set update to move to 10^-5
		if (sep <= min_sep) {

//			if(rank==0)puts("a\n");

			if (v==0) {

				double v1 = fabs(xx-lo);
				double v2 = fabs(hi-xx);

				return (v2<v1)?(hi-min_sep-xx):(lo+min_sep-xx);

			} else return (v<0)?(lo+min_sep - xx):(hi-min_sep - xx);
		}

		// do correction if update distance
		// is more than 10% distance to wall
		if (fabs(update) > 0.1*sep) {

			// if doing correction, modify update to use more knowledge of shape of
			// velocity gradient, i.e. that normal vel. is 0 on wall
//			update *= (1-0.5*dt*vn/sep);
//			if(rank==0)puts("b\n");

			// last bit is to get negative[positive] update if v<0[v>0]
			update = (exp(-dt*vn/sep)-1) * ((v<0)?sep:-sep);


		}

		// if the update will dump us out, move to 10^-5
		if ((v<0 && (xx+update) < (lo+min_sep)) || (v>0 && (xx+update) > (hi-min_sep))) {
//			if(rank==0)printf("c [%g, | %g from xx=%g, v=%g, dt=%g |, %g] \n",lo+min_sep,xx+update,
//					xx,t.v0[f],dt,hi-min_sep);


			return (v<0)?(lo+min_sep - xx):(hi-min_sep - xx);
		}
//			if(rank==0)puts("d\n");

		return update;
	}
};
#endif
