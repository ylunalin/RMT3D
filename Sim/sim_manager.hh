#ifndef SIM_MANAGER_
#define SIM_MANAGER_

#include <cmath>
#include "defs.hh"
#include "geometry.hh"
#include "object.hh"
#include "gsl/gsl_rng.h"
#include "sim_params.hh"
#include "matlib.hh"

using namespace OBJECT_DICT;

/** solid material constants */
struct sl_mat {
    /** Index */
    int index;
	/** density */
	double rho;
	/** inverse density */
	double invrho;
	/** shear modulus */
	double G;
	/** elastic wave speed */
	double c;
	/** artificial viscosity parameter */
	double ex_visc_mult;
	/** transition zone multiplier */
	double ev_trans_mult;
	/** shear wave and fluid time multiplier*/
	double dt_pad;
	/** transition zone extra viscosity time multiplier */
	double dt_ex_pad;
	/** the actual extra viscosity */
	double ex_mu;
	// constructors
	sl_mat():
		rho(1.), invrho(1./rho), G(1.),
		c(sqrt(G/rho)), ex_visc_mult(1.), ev_trans_mult(1.),
		dt_pad(0.5), dt_ex_pad(0.5) {}

	sl_mat(const sim_params * spars, const int index_):
		index(index_), rho(spars->srho[index]),invrho(1/rho),
		G(spars->shear_mod[index]),c(sqrt(G/rho)),
        ex_visc_mult(spars->ex_visc_mult),
		ev_trans_mult(spars->ev_trans_mult),
        dt_pad(spars->sdt_pad), dt_ex_pad(spars->dt_ex_pad)
        {
            mu(spars->min_dh());
        }

    inline void mu(double dx) {
        // dynamic viscosity
        // Unit: M/(LT)
        ex_mu = ex_visc_mult*dx*c*rho;
        //ex_mu = 0.01;
    }
	// assume isotropic spacing
	inline double cfl(double dx) {
		double xxsp3 = 3/(dx*dx);
		double dt1 = 0.5*rho/(ex_mu *(ev_trans_mult+1) *xxsp3) *dt_ex_pad;
		double dt0 = dx/c*dt_pad;
		return (dt0<dt1)?dt0:dt1;
	}
	void print_self(){
		printf("# Solid material property %d\n"
               "#	[Density                    %f M/L^3 ]\n"
               "#	[Shear modulus              %f M/(LT)]\n"
               "#	[Elastic wave speed         %f L/T   ]\n"
               "#	[Extra dynamic viscosity    %f M/(LT)]\n"
               "#	[Constants: ex_visc_mult %g ev_trans_mult %g dt_pad %g dt_ex_pad %g]\n"
			   , index, rho, G, c, ex_mu, ex_visc_mult, ev_trans_mult, dt_pad, dt_ex_pad);
	}
};

/** fluid material constants */
struct fl_mat {
	/** density */
	const double rho;
	/** inverse density */
	const double rhoinv;
	/** dynamic viscosity */
	const double mu;
	/** dt padding factor */
	const double dt_pad;

	// constructors
	fl_mat(const sim_params *spars)
        : rho(spars->frho),rhoinv(1/rho),
		mu(spars->fmu),dt_pad(spars->fdt_pad) {};

	void print_self() {
		printf("# Fluid material property\n"
               "#	[Density            %f M/L^3 ]\n"
               "#	[Viscosity          %f M/(LT)]\n"
               "#	[Constnat: dt_pat   %g       ]\n"
				, rho, mu, dt_pad);
	}
};

class sim_manager {

	public:
    /** Pointer to the config class */
    const sim_params* spars;
	/** whether fluid stress is computed using implicit methods */
	const bool impl;
	/** whether Godunov extrapolation is used for the advective term */
	const bool godunov;
	/** periodicity */
	const bool x_prd,y_prd,z_prd;
	/** whether an exact solution exists */
	const bool has_exact;
	/** number of immersed objects */
	int n_obj;
	/** number of layers to extrapolate */
	int nlayers;
	/** width of transition zone in units of grid spacing*/
	const double wt_n;
	/** lower x,y,z coords */
	const double ax,ay,az;
	/** upper x,y,z coords */
	const double bx,by,bz;
    /** spatial interval */
	const double dx,dy,dz;
    /** spatial interval */
	const double dxsp,dysp,dzsp;
	/** lengths */
	const double lx,ly,lz;
	/** width of transition zone */
	double eps;
	/** inverse of the width of transition zone */
	double eps_inv;
	/** Some constants to use in smooth heaviside calculations */
	double tf1, tf2, tf3;
	/** maximum dt allowed for fluid_3d simulation */
	double dt_reg;
    /** A stiffness constant for anchors */
    double K_stiff;
	/** the weight for 5x5x5 cube in extrapolation */
	double weight_fac;
    /** maximum velocity in the system.*/
    double vmax;
    /** Wall repulsion acceleration. */
    double wall_acc_default;
    double wall_acc_mult;


	/** pointers to immersed objects */
	object **objs;
    double *avg_velx;
    double *avg_vely;
    double *avg_velz;
    double *avg_velN;

    double *avg_x;
    double *avg_y;
    double *avg_z;
    double *avg_N;

	/** fluid material constants */
	fl_mat fm;
	sl_mat *sm_array;


	/** constructor */
	sim_manager(const sim_params * spars_) :
    spars(spars_),
	impl(spars->impl), godunov(spars->godunov),
	x_prd(spars->x_prd), y_prd(spars->y_prd), z_prd(spars->z_prd),
	has_exact(spars->has_exact), n_obj(spars->n_obj),
	nlayers(spars->nlayers), wt_n(spars->wt_n),
	ax(spars->ax), ay(spars->ay), az(spars->az),
	bx(spars->bx), by(spars->by), bz(spars->bz),
	dx(spars->dx), dy(spars->dy), dz(spars->dz),
	dxsp(1./dx), dysp(1./dy), dzsp(1./dz),
	lx(spars->lx), ly(spars->ly), lz(spars->lz),
    weight_fac(spars->weight_fac), vmax(1),  wall_acc_default(1), wall_acc_mult(1),
    objs(n_obj==0?NULL:new object*[n_obj]),
    avg_velx(n_obj==0?NULL:new double[n_obj]),
    avg_vely(n_obj==0?NULL:new double[n_obj]),
    avg_velz(n_obj==0?NULL:new double[n_obj]),
    avg_velN(n_obj==0?NULL:new double[n_obj]),
    avg_x(n_obj==0?NULL:new double[n_obj]),
    avg_y(n_obj==0?NULL:new double[n_obj]),
    avg_z(n_obj==0?NULL:new double[n_obj]),
    avg_N(n_obj==0?NULL:new double[n_obj]),
    fm(spars), sm_array(n_obj==0?NULL:new sl_mat[n_obj])
	{
		setup();
	};

    virtual void print_self(){
        printf("# Interface properties:\n");
        printf("#    Transition zone width %6.2g grid spacings\n"
               "#    i.e. %.6e\n"
               "#    maximum number of layers %d\n#\n",
                wt_n, eps, nlayers);
        printf("# Extrapolation weight factor %6.4f\n"
               "# Maximum velocity %g L/T\n",
            weight_fac,vmax);

    }

	virtual ~sim_manager();
	/** Set up both the maximum time step and the transition zone width */
	void setup_const();
    /** Create material constant and objects */
    void create_objs();
    inline void setup(){
        // Must be called in order, otherwise setup_const() doesn't see any objects
        reset_avgs();
        create_objs();
        setup_const();
    }
    virtual void reset_avgs(){
        for(int i=0;i<n_obj;i++){
            avg_velx[i] = 0.;
            avg_vely[i] = 0.;
            avg_velz[i] = 0.;
            avg_velN[i] = 0.;
            avg_x[i] = 0.;
            avg_y[i] = 0.;
            avg_z[i] = 0.;
            avg_N[i] = 0.;
        }
    }
    virtual void compute_avgs(){
        for(int i=0;i<n_obj;i++){
            if(avg_velN[i]!=0.) {
                avg_velx[i] /= avg_velN[i];
                avg_vely[i] /= avg_velN[i];
                avg_velz[i] /= avg_velN[i];
            } else {
                avg_velx[i] = 0;
                avg_vely[i] = 0;
                avg_velz[i] = 0;

            }

            if(avg_N[i]!=0.) {
                avg_x[i] /= avg_N[i];
                avg_y[i] /= avg_N[i];
                avg_z[i] /= avg_N[i];
            } else {
                avg_x[i] = 0;
                avg_y[i] = 0;
                avg_z[i] = 0;

            }
        }
    }
    void obligatory_cfl_recheck();

	/** utility to get the mininum of 3 numbers */
	inline double min3(double n1, double n2, double n3){
		double smst = (n1<n2)?n1:n2;
		smst = (smst<n3)?smst:n3;
		return smst;
	}
	/** returns the step size for a given number of
	 * grid points in each direction */

	inline double heaviside(double phi) const {
		if (phi >= eps) return 0;
		if (phi <= -eps) return 1;

        double tmp = 0.5  - phi*tf1 - sin(tf3*phi)*tf2;
        if(tmp<0) tmp = 0.;
        else if(tmp>1) tmp = 1;

		return tmp;
	};
	inline double tderiv_func_in(double x){
        double tmp = 0;
        if(x>-eps && x<eps) {tmp = 0.5+0.5*cos(x*tf3);}
        return tmp;
	}
	inline double delta_func_wall(double x){
		if(x<=0) return 100;
		else return (1-x)/x;
	}
	inline double coll_func(double x){
		return x>eps?0:(0.5-x*tf1);
	}
    inline double wall_trans_func(double phiv, double wtz_width, double wtzw_inv){
        if( phiv < - wtz_width) return 1;
        else if (phiv < wtz_width) return 0.5*(1-phiv*wtzw_inv);
        else return 0;
    }
	inline void set_extrap_layers(int n){ nlayers=n; }

	// Implementation see sim_manager.cc
	void add_obj(object *o);
	/** returns the max dt which satisfies a CFL-type criterion for
	 *  given spatial discretization dx */
	double cfl(double dx);

	/** specifies the intial velocity field */
	virtual void velocity(double x,double y,double z,
		double &u,double &v,double &w);

	/** specifies the intial pressure guess for multigrid */
	virtual void pressure(double x,double y,double z,double &p) {
        p=0;
    }

	/** specifies the exact velocity field (if known) at time t */
	virtual void exact(double x,double y,double z,double t,
		double h,double &u,double &v,double &w,double &p) {u=v=w=p=0;};

	/** specifies the body force at time t */
	virtual void fluid_acceleration(double x,double y,double z,double t,
		double &fx,double &fy,double &fz) { fx = fy = fz = 0; }

	virtual void solid_acceleration(const int obj_id, const double (&rval)[3], double x, double y, double z, double t, double &fx,double &fy,double &fz);
	//virtual void solid_acceleration(const int obj_id, const double (&rval)[3], double x, double y, double z, double t,	double &fx,double &fy,double &fz, const double vx=std::nan(""), const double vy=std::nan(""), const double vz=std::nan(""));

    // NOT A VIRTUAL FUNCTION
    void actuate(const int obj_id, const double (&rval)[3], const double (&rnval)[3], const double time, matrix & F){
        objs[obj_id]->actuate(rval, rnval, time, F);
    };

	virtual double phi(int obj_index, const double (&xi)[3]) const;
};

class sim_ftest : public sim_manager {
	public:
	sim_ftest(const sim_params *spars) : sim_manager(spars) {};
	/** specifies the initial velocity field */
	virtual void velocity(double x,double y,double z,
		double &u,double &v,double &w);
	/** specifies the intial pressure guess for multigrid */
	virtual void pressure(double x,double y,double z,double &p);
	/** specifies the exact velocity field (if known) at time t */
	virtual void exact(double x,double y,double z,double t,
		double h,double &u,double &v,double &w,double &p);
	/** specifies the body force at time t */
	virtual void fluid_acceleration(double x,double y,double z,double t,
		double &fx,double &fy,double &fz);
};

class sim_stest : public sim_manager {
	protected:
	double a,k,w;
	public:
	sim_stest(const sim_params *spars);
	/** specifies the initial velocity field */
	virtual void velocity(double x,double y,double z,
		double &u,double &v,double &w);
	/** specifies the intial pressure guess for multigrid */
	virtual void pressure(double x,double y,double z,double &p);
	/** specifies the exact velocity field (if known) at time t */
	virtual void exact(double x,double y,double z,double t,
		double h,double &u,double &v,double &w,double &p);
};

class sim_objects: public sim_manager {
	public:
    /** are there walls? */
    bool walls[6];
    /** wall positions */
    double wall_pos[6];
    const double wall_dist;
    const double winv;
    double gravity;

	sim_objects(const sim_params *spars) : sim_manager(spars),
    // In periodic cases, we need to not use the extrapolation from the wrong part of the object
    // In nonperiodic cases, we can safely extrapolate into the ghost
    wall_dist(spars->wall_dist), winv(1/wall_dist), gravity(spars->gravity){
        for(int i=0;i<6;i++) {
            walls[i] = spars->walls[i];
            wall_pos[i] = spars->wall_pos[i];
        }
        // NOTE this constant 20 is an arbitrary one
        wall_acc_default = 20*winv;
        set_wall_acc();
    };

    virtual void print_self(){
        printf("# Interface properties:\n");
        printf("#    Transition zone width %6.2g grid spacings\n"
               "#    i.e. %.6e\n"
               "#    maximum number of layers %d\n#\n",
                wt_n, eps, nlayers);
        printf("# Extrapolation weight factor %6.4f\n"
               "# Wall acceleration factor %6.4f (compared against %g)\n"
               "# Wall acceleration multiplier %6.4f\n"
               "# Minimum acceleration %6.4f\n"
               "# Maximum velocity %g L/T\n"
               "# Stiffness constant for anchoring %g T^{-2}\n"
               "# Wall in x [%d %d] position [%2.1f, %2.1f]\n"
               "# Wall in y [%d %d] position [%2.1f, %2.1f]\n"
               "# Wall in z [%d %d] position [%2.1f, %2.1f]\n"
               "# Wall activation distance %6.4f (dh=%6.4f)\n"
               "# Gravitational acceleration constant %6.4f (L/T^2)\n",
            weight_fac, wall_acc_default, 20*winv,
            wall_acc_mult,
            min_wall_acc,
            vmax, K_stiff,
            walls[0], walls[1], wall_pos[0], wall_pos[1],
            walls[2], walls[3], wall_pos[2], wall_pos[3],
            walls[4], walls[5], wall_pos[4], wall_pos[5],
            wall_dist, spars->min_dh(), gravity);
    }

    /** specifies the body force at time t */
    virtual void fluid_acceleration(double x,double y,double z,double t,
        double &fx,double &fy,double &fz) {
        fx = fy = fz = 0;
        fz += gravity;
    }
	virtual void solid_acceleration(const int obj_id, const double (&rval)[3], double x, double y, double z, double t,	double &fx,double &fy,double &fz);
    private:
    double min_wall_acc;
    void set_wall_acc(){
        // We compute the maximum settling speed (squared)
        // using the drag equation, assume drag coefficient is order 1
        // We are getting max( r^3/ h^3 * v^2/2h)
        double mvsq=-1;
        min_wall_acc = 20*winv;
        double mvol =-1;
        double rho_f = fm.rho;

        for (int i = 0; i < n_obj; i++) {
            double rho_s = sm_array[i].rho;
            double vsq = 0.5 * fabs((rho_s - rho_f) * gravity) / rho_f;
            vsq *= objs[i]->volume/objs[i]->primary_dim/objs[i]->primary_dim;
            if(vsq>mvsq) mvsq= vsq;

            if(objs[i]->volume > mvol) mvol = objs[i]->volume/objs[i]->primary_dim/objs[i]->primary_dim;

            // reuse vsq varialbe
            // calculate the minumum acceleration needed to support an object
            vsq = fabs((rho_s - rho_f) * gravity * objs[i]->volume / rho_s/ objs[i]->primary_dim / objs[i]->primary_dim *winv);
            if(vsq>min_wall_acc) min_wall_acc = vsq;
        }
        // The last multiplicative constant is to add some insurance
        // The accelerations in the wall contact region needs to stop the momentum of the entire object
        wall_acc_default = 0.5*mvsq*pow(winv, 4);
        wall_acc_mult = 0.5*mvol*pow(winv, 2);

        // settling velocity should also be taken into accout in cfl condition
        mvsq = sqrt(mvsq);
        if(mvsq>vmax) vmax = mvsq;
    }
};

class object_mills: public sim_objects {
	public:
	object_mills(const sim_params *spars) : sim_objects(spars) {};
	/** specifies the initial velocity field */
	virtual void velocity(double x,double y,double z,
		double &u,double &v,double &w);
	/** specifies the intial pressure guess for multigrid */
	virtual void pressure(double x,double y,double z,double &p);
	/** specifies the exact velocity field (if known) at time t */
	virtual void exact(double x,double y,double z,double t,
		double h,double &u,double &v,double &w,double &p);
	/** specifies the body force at time t */
	virtual void fluid_acceleration(double x,double y,double z,double t,
		double &fx,double &fy,double &fz);
};

#endif
