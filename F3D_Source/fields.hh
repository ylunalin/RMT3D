#ifndef F3D_FIELDS_HH
#define F3D_FIELDS_HH
#include <limits>
#include "defs.hh"
#include "common.hh"
#include "geometry.hh"
#include "stencil.hh"
#include "matlib.hh"
#include "sim_manager.hh"

using namespace FIELD_CONST;
using namespace BC;

/**
 * A simplied structure to store the velocities in the last timestep.
 * This is to compute a termination condition once a steady flow has developed.
 */

struct vel_field {
    double vel[3];
    vel_field () {
        for(int i=0;i<3;i++) vel[i] = 0.;
    }
};

/**
 * Structure storing field information at each grid point. Contains
 * velocity components u,v,w and pressure in vel[4];
 * velocitie derivatives du,dv,dw in devel[3];
 * monotonocity limited normal derivatives in nv_der[3][3];
 * extrapolated velocities at the faces in fvel[6][3];
 * tangential derivatives in the advective term in tang[3][3].
 */
struct field {
	/** velocities */
	double vel[3];
	/** pressure */
	double p,dp;
	/** auxiliary pressure in MAC solve */
    double q;
	/** time derivative of velocity */
	double dvel[3];
	/** normal spatial derivative of velocity */
	double nv_der[3][3]; // dim 1 dx/dy/dz dim 2 u, v, w
	/** extrapolated face velocities */
	double fvel[6][3];

#if defined(TANG)
	/** temporary quantities needed for tangential velocities */
	double tang[3][3]; // dim 1 dx/dy/dz dim 2 u, v, w
#endif

	/** stress (sans pressure). sigma[i][j] is partial_i u_j, located
	 * at center of face i. */
	double sigma[3][3];

    field () {
        p = 0; dp = 0;
        q = 0;
        for (int i=0;i<3;i++){
            vel[i] = 0;
            dvel[i] = 0;
            for(int j=0;j<3;j++){
                sigma[i][j] = 0;
                fvel[i][j] = 0;
                fvel[i+3][j] = 0;
                nv_der[i][j] = 0;
#if defined(TANG)
                tang[i][j] = 0;
#endif
            }
        }
    }

   void reset () { 
        // we don't reset p
        // because initial iterations may require p
        q = 0; dp = 0;
        for (int i=0;i<3;i++){
            vel[i] = 0;
            dvel[i] = 0;
            for(int j=0;j<3;j++){
                sigma[i][j] = 0;
                fvel[i][j] = 0;
                fvel[i+3][j] = 0;
                nv_der[i][j] = 0;
#if defined(TANG)
                tang[i][j] = 0;
#endif
            }
        }
    }

    /** General house keeping functions */
	inline int oid() const {
        return 0;
    }
	double get(const int f) const {
		if (f<3 && f>=0) return vel[f];
		else if( f==3) return p;
		else {
			p_fatal_error("error, field::get\n",PGMG_SETUP_ERROR);
			return std::numeric_limits<double>::quiet_NaN();
		}
	}
	inline double phi(const sim_manager *mgmt) const {
		return 0.0;
	}

    /** Functions for error calculation, when there is exact solution */
	void add_error(const double ex_soln[4], double (&e)[4]) const {
		for(int i=0;i<3;i++){
			e[i]  += (ex_soln[i] - vel[i])*(ex_soln[i] - vel[i]);
		}
		e[3] += (ex_soln[3]-p)*(ex_soln[3]-p);
	}
	inline void add_p_error(const double ep, double (&e)[4]) const {
		e[3] += (ep-p)*(ep-p);
	}
	double error(const double ex_soln[4]) const {
		double e[4] = {0,0,0,0};
		add_error(ex_soln,e);
		return e[0]+e[1]+e[2]+e[3];
	}

    /** Functions for updating velocity and velocity derivatives */
	inline void add_vel_derivs(const double dvel_[3]){
		for(int i=0;i<3;i++) dvel[i] += dvel_[i];
	}
	inline void set_vels(const double vel_[3]) {
		for(int i=0;i<3;i++) vel[i] = vel_[i];
	}
	inline void reset_vel_derivs() {
        for(int i=0;i<3;i++) dvel[i]=0.;
	}
	inline void update_vels() {
		for(int i=0;i<3;i++) vel[i] += dvel[i];
	}
	inline void set_mono_derivs(const int dim, const int vel, const double value){
		nv_der[dim][vel] = value;
	}
	/* Extrapolates velocities to the face at t=n+1/2, Yu(2003) Eqn 3.13.
	 * \param[in] (dh2[3], dt)  half spatial and time steps, respetively,
	 * \param[in] additive_f[3] addition terms from viscosity, body forces, and pressure at t=n,
	 * \param[in] tang_stab a flag indicating if tangential stability terms are computed seperately.
	 */
	void godunov_extrap(const double dh2[3], const double dt2, const double (&additive_f)[3]){
#if defined(TANG)
		double d;
		int i, j, dim, alt1, alt2;
		for(i=0;i<6;i++){
			dim = i/2;
			alt1 = (dim+1)%3;
			alt2 = (dim+2)%3;
			double norm_vel = this-> get(dim);
			d = dh2[dim] * (2*(i%2)-1);
			for(j=0;j<3;j++)
				fvel[i][j] = vel[j] + (d - dt2*norm_vel)*nv_der[dim][j] + dt2*(-tang[alt1][j] - tang[alt2][j] + additive_f[j]);
		}
#else
        double adv_vel[3] = {0,0,0};
        double d;
        int i, j, dim;
        for(i=0;i<3;i++) {
            for (j=0;j<3;j++){
                adv_vel[i]-=vel[j]*nv_der[j][i];
            }
            adv_vel[i]*=dt2;
        }
        for(i=0;i<6;i++){
            dim = i/2;
            d = dh2[dim] * (2*(i%2)-1);
            for(j=0;j<3;j++)
                fvel[i][j] = vel[j] + d*nv_der[dim][j] + adv_vel[j] + dt2*additive_f[j];
        }
#endif
	}


#if defined(TANG)
	/* Extrapolates tangential velocities to the face at t=n+1/2, Yu(2003) Eqn 3.30.
	 * \param[in] (dh2[3], dt)  half spatial and time steps, respetively,
	 */
	void inter_extrap(const double (&dh2)[3], const double dt2){
		double d;
		int i, j, dim;
		for(i=0;i<6;i++){
			dim = i/2;
			double norm_vel = vel[dim];
			d = dh2[dim] * (2*(i%2)-1);
			for(j=0;j<3;j++)
				fvel[i][j] = vel[j] + (d - dt2*norm_vel) *nv_der[dim][j];
		}
	}
#endif

#if defined(TANG)
	/* Computes tangential derivatives terms. Yu(2003) Eqn 3.29.
	 * \param[in] dhsp[3] inverse of spatial intervals,
	 * \param[in] strides[3] strides in 3 dimensions.
	 */
	void tang_terms(const double (&dhsp)[3], const int (&strides)[3]){
		field *fl = this-1, *fr = this+1, *ff = this-strides[1], *fb = this+strides[1], *fd = this-strides[2], *fu = this+strides[2];
		double adv_u = 0.5*(fr->fvel[0][0]+fvel[0][0])*dhsp[0], adv_v = 0.5*(fb->fvel[2][1]+fvel[2][1]) *dhsp[1], adv_w = 0.5*(fu->fvel[4][2]+fvel[4][2])*dhsp[2];
		tang[0][0] = adv_u*(fvel[1][0] - fl->fvel[1][0]); // u ux
		tang[0][1] = adv_u*(fr->fvel[0][1] - fvel[0][1]); // u vx
		tang[0][2] = adv_u*(fr->fvel[0][2] - fvel[0][2]); // u wx

		tang[1][0] = adv_v*(fb->fvel[2][0] - fvel[2][0]); // v uy
		tang[1][1] = adv_v*(fvel[3][1] - ff->fvel[3][1]); // v vy
		tang[1][2] = adv_v*(fb->fvel[2][2] - fvel[2][2]); // v wy

		tang[2][0] = adv_w*(fu->fvel[4][0] - fvel[4][0]); // w uz
		tang[2][1] = adv_w*(fu->fvel[4][1] - fvel[4][1]); // w vz
		tang[2][2] = adv_w*(fvel[5][2] - fd->fvel[5][2]); // w wz
	}
#endif

#if defined(TANG)
	/* Computes tangential derivatives terms. Yu(2003) Eqn 3.29.
	 * \param[in] dhsp[3] inverse of spatial intervals,
     * \param[in] strides[3] strides in 3 dimensions,
	 * \param[out] vel_adv[3] avergae advective velocities.
     */
	void tang_terms(const double (&dhsp)[3], const int (&strides)[3], double (&adv_vel) [3]){
		field *fl = this-1, *fr = this+1, *ff = this-strides[1], *fb = this+strides[1], *fd = this-strides[2], *fu = this+strides[2];
		double adv_u = 0.5*(fr->fvel[0][0]+fvel[0][0])*dhsp[0], adv_v = 0.5*(fb->fvel[2][1]+fvel[2][1]) *dhsp[1], adv_w = 0.5*(fu->fvel[4][2]+fvel[4][2])*dhsp[2];
		tang[0][0] = adv_u*(fvel[1][0] - fl->fvel[1][0]); // u ux
		tang[0][1] = adv_u*(fr->fvel[0][1] - fvel[0][1]); // u vx
		tang[0][2] = adv_u*(fr->fvel[0][2] - fvel[0][2]); // u wx

		tang[1][0] = adv_v*(fb->fvel[2][0] - fvel[2][0]); // v uy
		tang[1][1] = adv_v*(fvel[3][1] - ff->fvel[3][1]); // v vy
		tang[1][2] = adv_v*(fb->fvel[2][2] - fvel[2][2]); // v wy

		tang[2][0] = adv_w*(fu->fvel[4][0] - fvel[4][0]); // w uz
		tang[2][1] = adv_w*(fu->fvel[4][1] - fvel[4][1]); // w vz
		tang[2][2] = adv_w*(fvel[5][2] - fd->fvel[5][2]); // w wz

		adv_vel[0] = adv_u;
		adv_vel[1] = adv_v;
		adv_vel[2] = adv_w;
	}
#endif

    /** Functions for communication - putting data into buffer*/
	template<int flags>
	inline void pack(double *(&b)) const {
		if(flags&1) {
            for(int c=0;c<3;c++,b++) *b=vel[c];
        }
		if(flags&2) {
            *(b++)=p;
            *(b++)=dp;
            *(b++)=q;
        }
		if(flags&4) {
            for(int f=0;f<6;f++) for(int c=0;c<3;c++,b++) *b=fvel[f][c];
        }
		if(flags&2048) {
            for(int c=0;c<3;c++,b++) *b=dvel[c];
        }
	}
	template<int flags>
	inline void unpack(double *(&b)){
		if(flags&1) {
            for(int c=0;c<3;c++,b++) vel[c]=*b;
        }
		if(flags&2) {
            p=*(b++);
            dp=*(b++);
            q=*(b++);
        }
		if(flags&4) {
            for(int f=0;f<6;f++) for(int c=0;c<3;c++,b++) fvel[f][c]=*b;
        }
		if(flags&2048) {
            for(int c=0;c<3;c++,b++) dvel[c]=*b;
        }
	}

}; //End of field struct

/** Reference map field. */
struct ref_map {
	/** id number corresponding to which object the reference map refers to */
	unsigned int c;
	/** The reference map. */
	double x[3];
	/** The new reference map in predictive step. */
	double xpred[3];
	/** normal spatial derivative of pos vectors*/
	double nv_der[3][3]; // dim 1 dx/dy/dz dim 2 X,Y,Z
	/** extrapolated face pos vectors */
	double fx[6][3];

#if defined(TANG)
	/** temporary quantities needed for tangential pos vectors */
	double tang[3][3]; // dim 1 dx/dy/dz dim 2 X,Y,Z
#endif

	/** stress. sigma[i][j] is partial_i u_j, located
	 * at center of face i. */
	double sigma[3][3];

    // Needed for elastic energy output
    double elastic_energy[3];
    double detF[3];
    // Needed for extrapolation
    double grad_phi[3][3];
	double act_power[3];

#if defined(DEBUG)
    // Debug fields
	/** Can we not compile some of these fields to save space? */
    double grad_phi_diff[3];
    double solid_stress_normal[3];
    double solid_stress_shear[3];
    double coll_stress_normal[3];
    double coll_stress_shear[3];
    double coll_traction_tot[3];
#endif

    ref_map() {
        c = 0;
        for (int i=0;i<3;i++){
            x[i] = 0;
            xpred[i] = 0;

#if defined(DEBUG)
            grad_phi_diff[i] = 0;
            solid_stress_normal[i] = 0;
            solid_stress_shear[i] = 0;
            coll_stress_normal[i] = 0;
            coll_stress_shear[i] = 0;
            coll_traction_tot[i] = 0;
#endif

            elastic_energy[i] = 0;
            detF[i] = 0;
            for(int j=0;j<3;j++){
                grad_phi[i][j] = -1.e200;
                sigma[i][j] = 0;
                fx[i][j] = 0;
                fx[i+3][j] = 0;
                nv_der[i][j] = 0;
#if defined(TANG)
                tang[i][j] = 0;
#endif
            }
        }
    }

    void reset () {
        c = 0;
        for (int i=0;i<3;i++){
            x[i] = 0;
            xpred[i] = 0;
#if defined(DEBUG)
            grad_phi_diff[i] = 0;
            solid_stress_normal[i] = 0;
            solid_stress_shear[i] = 0;
            coll_stress_normal[i] = 0;
            coll_stress_shear[i] = 0;
            coll_traction_tot[i] = 0;
#endif
            elastic_energy[i] = 0;
            detF[i] = 0;
            for(int j=0;j<3;j++){
                grad_phi[i][j] = -1.e200;
                sigma[i][j] = 0;
                fx[i][j] = 0;
                fx[i+3][j] = 0;
                nv_der[i][j] = 0;
#if defined(TANG)
                tang[i][j] = 0;
#endif
            }
        }
    }

	/** General house keeping functions */
    inline void normalize_grad_phi(int F){
        double norm=0;
        for (int i=0;i<3;i++){
            norm += grad_phi[F][i] * grad_phi[F][i];
        }
        norm = sqrt(norm);
        if(norm!=0.){
            double ninv = 1./norm;
            for (int i=0;i<3;i++){
                grad_phi[F][i] = ninv*grad_phi[F][i];
            }
        }
    }

    inline void clear(){
        c = 0;
        for (int i=0;i<3;i++){
            x[i] = 0;
            xpred[i] = 0;
#if defined(DEBUG)
            grad_phi_diff[i] = 0;
            solid_stress_normal[i] = 0;
            solid_stress_shear[i] = 0;
            coll_stress_normal[i] = 0;
            coll_stress_shear[i] = 0;
            coll_traction_tot[i] = 0;
#endif
            elastic_energy[i] = 0;
            detF[i] = 0;
            for(int j=0;j<3;j++){
                grad_phi[i][j] = 0;
                sigma[i][j] = 0;
                fx[i][j] = 0;
                fx[i+3][j] = 0;
                nv_der[i][j] = 0;
#if defined(TANG)
                tang[i][j] = 0;
#endif
            }
        }
    }
	inline int oid() const {
		return c&omask;
	}
	inline int lid() const {
		return (c&lmask) >> ldigits;
	};
	double phi(const sim_manager *mgmt,const int type=0) const {
        int o = oid();
		if(o == 0) return nan("");
		else return (type==0)?mgmt->phi(o-1, x):mgmt->phi(o-1, xpred);
	}
	inline double frac(const sim_manager *mgmt) const {
        return mgmt->heaviside(phi(mgmt));
    }

	double get(const int f) const {
		if (f<=2 && f>=0) return x[f];
		else {
			p_fatal_error("error, field::get\n",PGMG_SETUP_ERROR);
			return std::numeric_limits<double>::quiet_NaN();
		}
	}
	inline void set_x(const double (&x_)[3]){
		for(int i =0;i<3;i++)x[i] = x_[i];
	}
	inline void set_xpred(const double (&x_)[3]){
		for(int i =0;i<3;i++) xpred[i] = x_[i];
	}
	inline void pre_update_ref_map(const double (&cx)[3]){
		for(int i=0;i<3;i++){
				xpred[i] = x[i] + cx[i];
		}
	}

    /** Functions used in computing "Godunov" extrapolations of reference map */
	inline void set_half_timestep(double oldx[3]){
		for(int i=0;i<3;i++) x[i] = 0.5*(oldx[i] + xpred[i]);
	}
    inline void unset_half_timestep(double oldx[3]){
        for(int i=0;i<3;i++) x[i] = oldx[i];
    }
	inline void set_half_timestep(){
		for(int i=0;i<3;i++) x[i] = 0.5*(x[i] + xpred[i]);
	}
    inline void unset_half_timestep(){
        for(int i=0;i<3;i++) x[i] = 2*x[i] - xpred[i];
    }
	inline void set_full_timestep(){
		for(int i=0;i<3;i++) x[i] =  xpred[i];
	}
	inline void set_mono_derivs(int dim, int vel, double value){
		nv_der[dim][vel] = value;
	};
	inline bool neigh_inside(int ind, const ref_map *rm0) const {
		return (rm0[ind].oid() == oid());
	}
	// check if a point is abutting a primary point, orthogonally
	inline bool abut_primary(int ind, ref_map *rm0, sim_manager *mgmt, const int strides[3]) const {
		// get this reference map's object id
		int o = oid();
		for(int i=0;i<3;i++){
			// if the any of orthogonal neighbor has the same object
			// id same as this reference map, return true
			printf("obj=%d, dim[%d]: stride -%d -> c=%d, phi %g \n",o,i,strides[i], rm0[ind-strides[i]].c, rm0[ind-strides[i]].phi(mgmt));
			printf("obj=%d, dim[%d]: stride %d -> c=%d, phi %g \n",o,i,strides[i], rm0[ind+strides[i]].c, rm0[ind+strides[i]].phi(mgmt));
			if(rm0[ind+strides[i]].oid() == o) return true;
			if(rm0[ind-1*strides[i]].oid() == o) return true;
		}
		// otherwise, return there's no neighbor being inside
		return false;
	}
	inline bool abut_primary(int ind, ref_map *rm0, const int strides[3]) const {
		// get this reference map's object id
		int o = oid();
		for(int i=0;i<3;i++){
			// if the any of orthogonal neighbor has the same object
			// id same as this reference map, return true
			if(rm0[ind+strides[i]].oid() == o) return true;
			if(rm0[ind-1*strides[i]].oid() == o) return true;
		}
		// otherwise, return there's no neighbor being inside
		return false;
	}

	inline ref_map* step(int del,int ind,ref_map **f0,int *n0,sim_manager *mgmt,ref_map *rm0) const {
		int o = oid();

		// need to look in mmap and object for neighbor

		// first object...
		if (rm0[ind+del].oid() == o) return rm0 + ind + del;

		// ...and mmap
		for (int j = 0; j < n0[ind+del]; j++)
			if (f0[ind+del][j].oid() == o) return f0[ind+del] + j;

		return NULL;
	};

#if defined(TANG)
	/* Extrapolates tangential velocities to the face at t=n+1/2, Yu(2003) Eqn 3.30.
	 * \param[in] f the fluid velocities at the node,
	 * \param[in] (dh2[3], dt)  half spatial and time steps, respetively,
	 */
	inline void inter_extrap(field &f, const double dh2[3], const double dt2){
		double d;
		int i, j, dim;
		for(i=0;i<6;i++){
			dim = i/2;
			double norm_vel = f.vel[dim];
			d = dh2[dim] * (2*(i%2)-1);
			for(j=0;j<3;j++){
				fx[i][j] = x[j] + (d - dt2*norm_vel) *nv_der[dim][j];
			}
		}
	};
#endif

#if defined(TANG)
	/* Computes tangential derivatives terms. Yu(2003) Eqn 3.29.
	 * \param[in] dhsp[3] inverse of spatial intervals,
     * \param[in] strides[3] strides in 3 dimensions,
	 * \param[in] vel_adv[3] avergae advective velocities.
     */
	inline void tang_terms(const double dhsp[3], const int strides[3], double adv_vel [3]){
		printf("ref_map::tang_terms(): Must implement stepping with ref_map.step().\n");
		ref_map *fl = this-1, *fr = this+1, *ff = this-strides[1], *fb = this+strides[1], *fd = this-strides[2], *fu = this+strides[2];
		double adv_u = adv_vel[0], adv_v = adv_vel[1], adv_w = adv_vel[2];

		tang[0][0] = adv_u*(fx[1][0] - fl->fx[1][0]); // u Xx
		tang[0][1] = adv_u*(fr->fx[0][1] - fx[0][1]); // u Yx
		tang[0][2] = adv_u*(fr->fx[0][2] - fx[0][2]); // u Zx

		tang[1][0] = adv_v*(fb->fx[2][0] - fx[2][0]); // v Xy
		tang[1][1] = adv_v*(fx[3][1] - ff->fx[3][1]); // v Yy
		tang[1][2] = adv_v*(fb->fx[2][2] - fx[2][2]); // v Zy

		tang[2][0] = adv_w*(fu->fx[4][0] - fx[4][0]); // w Xz
		tang[2][1] = adv_w*(fu->fx[4][1] - fx[4][1]); // w Yz
		tang[2][2] = adv_w*(fx[5][2] - fd->fx[5][2]); // w Zz
	};
#endif

	/* Extrapolates velocities to the face at t=n+1/2, Yu(2003) Eqn 3.13.
	 * \param[in] (dh2[3], dt)  half spatial and time steps, respetively,
	 * \param[in] additive_f[3] addition terms from viscosity, body forces, and pressure at t=n,
	 * \param[in] tang_stab a flag indicating if tangential stability terms are computed seperately.
	 */
	inline void godunov_extrap(const double dh2[3], const double dt2, field &f){
		double d;
		int i, j, dim, alt1, alt2;
		for(i=0;i<6;i++){
			dim = i/2;
			alt1 = (dim+1)%3;
			alt2 = (dim+2)%3;
			double norm_vel=f.get(dim);
			d = dh2[dim] * (2*(i%2)-1);
#if defined(TANG)
            for(j=0;j<3;j++)
                fx[i][j] = x[j] + (d - dt2*norm_vel)*nv_der[dim][j] - dt2*(tang[alt1][j] + tang[alt2][j]);
#else
            double alt1v = f.get(alt1);
            double alt2v = f.get(alt2);
            for(j=0;j<3;j++){
                fx[i][j] = x[j] + (d - dt2*norm_vel)*nv_der[dim][j] - dt2*(alt1v*nv_der[alt1][j] + alt2v*nv_der[alt2][j]);
            }
#endif
		}
	}

	/** Functions used in communication, either in fluid_3d class or extra_mmap class */
	template<int flags>
	inline void unpack(double *(&b)) {
		if (flags&( 8| 128)) c=(int)*(b++);
		if (flags&(16| 256)) for(int i=0;i<3;i++,b++) x[i]=*b;
		if (flags&(32| 512)) for(int i=0;i<3;i++,b++) xpred[i]=*b;
		if (flags&(64|1024)) for(int f=0;f<6;f++) for(int i=0;i<3;i++,b++) fx[f][i]=*b;
	}
	template<int flags>
	inline void pack(double *(&b)) {
		if (flags&( 8| 128)) *(b++)=c+0.5;
		if (flags&(16| 256)) for(int i=0;i<3;i++,b++) *b=x[i];
		if (flags&(32| 512)) for(int i=0;i<3;i++,b++) *b=xpred[i];
		if (flags&(64|1024)) for(int f=0;f<6;f++) for(int i=0;i<3;i++,b++) *b=fx[f][i];
	}
	template<int flags>
	inline void shift(int face,double amt) {
		int dir = face/2;
		if (flags&(16| 256)) x[dir] += amt;
		if (flags&(32| 512)) xpred[dir] += amt;
		if (flags&(64|1024)) for(int f=0;f<6;f++)fx[f][dir] += amt;
	}

	//TODO could only pack/unpack x/fx if c > 0
	inline void pack_x(double *(&b)) {
		*(b++) = c+0.5;
		for(int i=0;i<3;i++) *(b++) = x[i];
	};
	inline void unpack_x(double *(&b)) {
		c = (int) *(b++);
		for(int i=0;i<3;i++) x[i] = *(b++);
	};
	inline void pack_fxs(double *(&b)) {
		*(b++) = c+0.5;
		for(int i =0;i<6;i++) for (int j=0;j<3;j++)
			*(b++) = fx[i][j];
	};
	inline void unpack_fxs(double *(&b)) {
		c = (int) *(b++);
		for(int i =0;i<6;i++)
			for (int j=0;j<3;j++){
				fx[i][j] = *(b++);
		}
	};
	inline void pack_xpred(double *(&b)) {
		*(b++) = c+0.5;
		for(int i=0;i<3;i++) *(b++) = xpred[i];
	};
	inline void unpack_xpred(double *(&b)) {
		c = (int) *(b++);
		for(int i=0;i<3;i++) xpred[i] = *(b++);
	};
}; // End of ref_map struct

/** A tiny vector to store the id and the fraction of solid in force calculation*/
struct solid_frac {
       ref_map *rp;
       double frac;
       solid_frac(ref_map *rp_, double frac_): rp(rp_), frac(frac_) {}
       inline void normalize(double norm_factor){
               frac /= norm_factor;
       }
};

/**
 * A structure containing boundary condition information
 */
struct face_gm {

    /** normal face aligned with cartesian direction +, else - */
    int sgn;
    /** which cartesian direction is the normal */
    int n;
	/** velcoity boundary types */
	int bct[3];
	/** tangential directions */
	int t[2];
	/** strides, including ghost nodes*/
	int len4[3];
	/** bounds for the three directions */
	int len[3];
    /** the cell face to set by b.c. in the outer most real, and inner most ghost layer */
    int ghost_cf[2];
	/** offset to get to the first non-ghost node on this face */
	int u0_off;
	/** velocity boundary condition values */
	double bcv[3];
	/** spatial discretizations of cartesian direction */
	double h;
	/** pointer to first real element in face */
	field *u0;
    /** coordinates */
    double *coords[3];
    double low_bound[3];
    double hi_bound[3];

    /** Setting up the boundary condition on the face */
	inline void setup(int fnum,geometry *gm,int (&bct_)[3], double (&bcv_)[3],double h_,field *u0_, double *lx0_, double *ly0_, double *lz0_, const double (&low)[3], const double (&hi)[3]) {

		for (int i = 0; i < 3; i++) {
			bct[i] = bct_[i];
			bcv[i] = bcv_[i];
            low_bound[i] = low[i];
            hi_bound[i] = hi[i];
		}

		n = fnum/2;
		t[0] = 0;
		t[1] = 2;
		switch (n) {
		case 0: t[0] = 1; break;
		case 2: t[1] = 1; break;
		}

		len[0] = gm->sm;
		len[1] = gm->sn;
		len[2] = gm->so;
		len4[0] = 1;
		len4[1] = (len[0]+4);
		len4[2] = len4[1]*(len[1]+4);

		h = h_;

        coords[0] = lx0_;
        coords[1] = ly0_;
        coords[2] = lz0_;

		if (fnum % 2 == 0) {
			u0_off = 0;
		} else {
			u0_off = (len[fnum/2]-1)*len4[fnum/2];
            coords[fnum/2] += len[fnum/2] - 1;
		}
		u0 = u0_ + u0_off;


		len[n] = 1;
        sgn=1;
		if (fnum%2==0){
            len4[n] *= -1;
            sgn=-1;
        }

        const int inner_ghost_cf = 2*n + (sgn<0);
        ghost_cf[0]=fnum;
        ghost_cf[1]=inner_ghost_cf;
	};

    // Getting the fluid node relative to u0, in the normal direction
    // positive ind[3] is stepping in normal direction
	inline field* node(const int (&ind)[3]) const {
        int tmp_ind =  0;
        for(int i=0;i<3;i++)  tmp_ind += ind[i]*len4[i];
        return u0+tmp_ind;
	};

    // Getting the mirored node in the real domain of a given ghost node
	inline field* mirror(const int (&inds)[3]) const {
        int tmp [3] = {inds[0], inds[1], inds[2]};
        // if mirroring ghost layer one, 1-1=0
        // if mirroring ghost layer two, 1-2=-1
        tmp[n] = 1-inds[n];
		return node(tmp);
	};

    // abosolute distance to the boundary face in normal direction
	inline double abs_dist(const int (&inds)[3]) {
		return fabs(0.5*h*(2*inds[n] - 1));
	};

    // set cell center values in the ghost region
	void set_ghost_cc(bool verbose=false) {
		int start[3];
		int end[3];

		for (int i = 0; i < 3; i++) {
			start[i] = -2;
			end[i] = len[i]+2;
		}

		start[n] = 1;
		end[n] = 3;

		int inds[3];
        // The index in the normal direction can be either 1 or 2
		for (inds[2] = start[2]; inds[2] < end[2]; inds[2]++)
			for (inds[1] = start[1]; inds[1] < end[1]; inds[1]++)
				for (inds[0] = start[0]; inds[0] < end[0]; inds[0]++) {

            // get the ghost node, and the mirroring real node
			field *u = node(inds);
			field *u_mir = mirror(inds);

            for(int i=0; i<3; i++){
                if (bct[i] == DIRICHLET) {
                        u->vel[i] = 2*bcv[i] - u_mir->vel[i];
                } else if (bct[i] == NEUMANN) {
                    // get the distance between ghost node and the boundary face
                    // along normal direction
                    double dn = 2*abs_dist(inds);
                    u->vel[i] = sgn*dn*bcv[i] + u_mir->vel[i];
                } else {
                    p_fatal_error("bad boundary condition type",1);
                }
            }
		}
	}

    // set cell face values in the inner ghost layer
    // and the outer most real layer
	void set_ghost_cf(bool verbose=false) {
		int start[3];
		int end[3];

		for (int i = 0; i < 3; i++) {
			start[i] = -2;
			end[i] = len[i]+2;
		}

		start[n] = 1;
		end[n] = 2;

		int inds[3];
		for (inds[2] = start[2]; inds[2] < end[2]; inds[2]++)
			for (inds[1] = start[1]; inds[1] < end[1]; inds[1]++)
				for (inds[0] = start[0]; inds[0] < end[0]; inds[0]++) {

            double f_bc_vals[3];
            compute_ghost_cf(inds, f_bc_vals);

			field * u = node(inds);
            field * u_mir = mirror(inds);

			for (int i = 0; i < 3; i++){
                // set cell face of the real node
                u_mir->fvel[ghost_cf[0]][i] = f_bc_vals[i];
                // set cell face of the ghost node, should be setting the same face, but index is opposite of that of the real node
				u->fvel[ghost_cf[1]][i] = f_bc_vals[i];
            }
		}
	}

    void compute_ghost_cf(const int (&inds)[3], double (&f_bc_vals)[3]){
        field *y0 = mirror(inds);
        int tmp [3] = {inds[0], inds[1], inds[2]};
        tmp[n] =  2;
        field *y1 = mirror(tmp);

        for(int i=0;i<3;i++){
            f_bc_vals[i] = bcv[i];
            // If NEUMANN, use two outer most real value, and the slope at the wall, to fit a quadratic
            if(bct[i] == NEUMANN){
                f_bc_vals[i] = (3*h*bcv[i]*sgn + 9*y0->vel[i] - y1->vel[i])*0.125;
            }
        }
    }
};
#endif
