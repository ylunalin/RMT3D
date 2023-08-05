#ifndef F3D_OBJ_HH
#define F3D_OBJ_HH

#include <cmath>
#include <cstdio>
#include "defs.hh"
#include "sim_params.hh"
#include "matlib.hh"

using namespace OBJECT_DICT;
class object {
  public:

	/** whether this subclass has an active stress method */
    bool act_stress_method=false;

    static const double fuzzy_eps;
	/** Basic specificiations, center and primary dimension*/
	double c[3];
    double primary_dim;
    double omega0[3];
    double volume;
    /** Extra specifications should be extended by children class*/
	object(): primary_dim(1), omega0{0.,0.,0.}, volume(pow(primary_dim, 3)) {
		c[0] = c[1] = c[2] = 0;
	};
	object(const double *basic_specs, bool v=false):
        omega0{0.,0.,0.} {
        for(int i=0;i<3;i++) c[i] = basic_specs[i];
        primary_dim = basic_specs[3];
        if(v) print_self();
	};
	virtual ~object() {}
    virtual void print_self(){
		printf( "# Object:\n"
                "#	[Center (%f %f %f), primary dimension %f]\n",
				c[0], c[1], c[2], primary_dim);

    };
	virtual double phi(const double (&xi)[3]) {return -1;};
	virtual void rm0(const double (&x)[3],double (&xi)[3]) {
		for (int cc = 0; cc < 3; cc++) xi[cc] = x[cc];
	}
    static object* alloc(const int type_name, const double *bp, const double *ep);
    double fuzzy_max(const double d1, const double d2);
    double fuzzy_min(const double d1, const double d2);

    // used to set initial conditions
    virtual void velocity (const double (&rval) [3], double &u,double &v,double &w) {
        double L[3];
        for(int i=0;i<3;i++) L[i] = rval[i]-c[i];
        u += omega0[1]*L[2] - omega0[2]*L[1];
        v -= omega0[0]*L[2] - omega0[2]*L[0];
        w += omega0[0]*L[1] - omega0[1]*L[0];
    };
    virtual void actuate(const double (&rval)[3], const double (&rnval)[3], const double time, matrix &F){};
    virtual void passive_acc_deformation(const double  (&rval)[3], const double eps, double x, double y, double z, double t, double &dx, double &dy, double &dz, double & dtoe1, double & dtoe2){
        dx=0; dy=0; dz=0; dtoe1=2*eps; dtoe2=2*eps;
    };

	virtual void active_stress(double tt,double (&x)[3],int row,
			const matrix &F,sym_matrix &as) {

		char mesg[] = "this object doesn't have an active stress method\n";
		p_fatal_error(mesg,-1);
	}
};

class shear_cube : public object {
	public:
	double a,k;
	shear_cube(double a_,double k_) :
		object(),a(a_),k(k_) {};
	~shear_cube() {}
	virtual double phi(const double (&xi)[3]);
	virtual void rm0(const double (&x)[3],double (&xi)[3]);
};

class basic_sphere : public object {
	public:
    double R;
	basic_sphere(const double * basic_specs, bool v=false) : object(basic_specs), R(primary_dim){
        volume = M_PI * 4/3 *pow(primary_dim, 3);
        if(v) print_self();
    };
	virtual ~basic_sphere() {}
	virtual double phi(const double (&xi)[3]);
	virtual void print_self(){
		printf( "# Basic sphere:\n"
                "#	[Center (%f %f %f), radius %f]\n",
				c[0], c[1], c[2], R);
	}
};

class sphere : public basic_sphere {
	public:
    int stretch_dir;
	double sqrt_lam;
    double eigvals[3];
	sphere(const double *basic_specs, const double * extra_specs, bool v=false) :
		basic_sphere(basic_specs), stretch_dir(int(extra_specs[0])),
        sqrt_lam(extra_specs[1]),
		eigvals{sqrt_lam, sqrt_lam, sqrt_lam} {
            for(int i=0;i<3;i++) omega0[i] = extra_specs[2+i];
            eigvals[stretch_dir] = 1./(sqrt_lam*sqrt_lam);
            if(v) print_self();
    }

	~sphere() {}
	virtual void rm0(const double (&x)[3],double (&xi)[3]) {
		for(int i=0;i<3;i++){
			xi[i] = (x[i]-c[i])*eigvals[i] + c[i];
		}
	}
	virtual void print_self(){
		printf( "# General sphere:\n"
                "#	[Center (%f %f %f), radius %f]\n"
		        "#	[Initial stretch in dir %d: %f]\n"
                "#	[Initial spin omega (%f, %f, %f)]\n",
				c[0], c[1], c[2], R, stretch_dir, sqrt_lam*sqrt_lam, omega0[0], omega0[1], omega0[2]);
	}
};

class rod : public object {
	public:
	// R - radius of caps
	// L - length of rod (not including caps)
	double R,L;
    double orientation[3];
    double clow[3], chi[3];
	rod(const double *basic_specs, const double * extra_specs, bool v=false):
		object(basic_specs), R(primary_dim), L(extra_specs[0]),
        orientation{extra_specs[1], extra_specs[2], extra_specs[3]} {
        volume = M_PI * 4/3 *pow(R, 3) + M_PI * R * R * L;
        // normalize
        double norm = 0.;
        for(int i=0;i<3;i++) norm += orientation[i]*orientation[i];
        norm = sqrt(norm);
        // if the norm is zero accidental, we set it to align with y-axis
        if (norm==0.) {
            orientation[0] = 0.;
            orientation[1] = 1.;
            orientation[2] = 0.;
        } else {
            for(int i=0;i<3;i++) orientation[i] /= norm;
        }

        for(int i=0;i<3;i++) {
            clow [i] = c[i] - orientation[i]*L*0.5;
            chi [i] = c[i] + orientation[i]*L*0.5;
        }
        for(int i=0;i<3;i++) {
            omega0[i] = extra_specs[4+i];
        }
        if(v) print_self();
    };
	virtual ~rod() {}
	virtual double phi(const double (&xi)[3]);
	virtual void print_self(){
		printf(	"# Spherocynlinder:\n"
                "#	[Center (%f %f %f)]\n"
				"#	[Cap radius %g, body length %g\n"
				"#	[Orientation (%g %g %g)\n"
                "#	[Initial spin omega (%f, %f, %f)]\n"
                "#	[Cap pos (%f %f %f) (%f, %f, %f)]\n",
				c[0], c[1], c[2], R, L,
                orientation[0], orientation[1], orientation[2],
                omega0[0], omega0[1], omega0[2],
                clow[0], clow[1], clow[2], chi[0], chi[1], chi[2]);
	}

    inline double cap_plane_low(const double (&xi)[3]){
        double tmp =0;
        for(int i=0;i<3;i++) tmp += orientation[i] * (xi[i] - clow[i]);
        return tmp;
    }
    inline double cap_plane_hi(const double (&xi)[3]){
        double tmp =0;
        for(int i=0;i<3;i++) tmp += orientation[i] * (xi[i] - chi[i]);
        return tmp;
   }
};

class active_rod_as : public rod {
	public:

	double bending_dir[3];
	double om,k;
	double b_mag;
	double Rhead;

	// bending moment direction - extra_specs[4,5,6] (overwrite/read omega0 slots)
	// frequency - extra_specs[7]
	// wavenumber - extra_specs[8]
	// frequency - extra_specs[9]
	active_rod_as(const double *basic_specs,const double *extra_specs,bool v=false):
		rod(basic_specs,extra_specs,false),
		bending_dir{omega0[0],omega0[1],omega0[2]},
		om(2*M_PI*extra_specs[7]),k(2*M_PI*extra_specs[8]),
		b_mag(extra_specs[9]),Rhead(0.25)
		{if (v) print_self(); act_stress_method=true;
		for(int i=0;i<3;i++)omega0[i]=0;}

	double phi(const double (&xi)[3]) {
		// this repeates the rod phi function
		// FIXME repeated code
		
		// cap centers
		double pt_vec[3] {0,0,0}, dot_product=0, norm=0;
		for(int i=0;i<3;i++) {
			pt_vec [i]  = xi[i] - c[i];
			norm += pt_vec[i]*pt_vec[i];
			dot_product += pt_vec[i] * orientation[i];
		}

		// get dist from each
		double d_cyl=0, d_low=0, d_hi=0, d_head=0;
		for(int i =0;i<3;i++) {
			d_low += (xi[i] - clow[i])*(xi[i] - clow[i]);
			d_hi  += (xi[i] - chi[i]) *(xi[i] - chi[i]);
		}
		d_cyl = sqrt(norm - dot_product*dot_product) - R;
		d_low = sqrt(d_low)-R;
		d_head = sqrt(d_hi)- Rhead;
		d_hi  = sqrt(d_hi) -R;

		double tmp_phi = 0;

		if (fabs(dot_product)>=0.5*L) {
			tmp_phi = std::min(d_low,d_hi);
		} else {
			tmp_phi = d_cyl;
		}

		return fuzzy_min(tmp_phi,d_head);

		// Test which side of the rod we are on
//		double tmp_phi = 0;
//		if(fabs(dot_product)>L*0.5) {
//			tmp_phi = std::min(d_low, d_hi);
//		} else {
//			tmp_phi = d_cyl;
//		}
//
//
//		return tmp_phi;
	}

	void print_self() {
		printf(	"# Active rod:\n"
                "#	[Center (%f %f %f)]\n"
				"#	[Cap radius %g, body length %g\n"
				"#	[Orientation (%g %g %g)\n"
                "#	[Initial spin omega (%f, %f, %f)]\n"
                "#	[Cap pos (%f %f %f) (%f, %f, %f)]\n"
                "#	[Bending dir (%f %f %f)]\n"
                "#	[Frequency (%f)]\n"
                "#	[Wavenumber (%f)]\n"
				"#	[Bending moment (%f)]\n",
				c[0], c[1], c[2], R, L,
                orientation[0], orientation[1], orientation[2],
                omega0[0], omega0[1], omega0[2],
                clow[0], clow[1], clow[2], chi[0], chi[1], chi[2],
				bending_dir[0],bending_dir[1],bending_dir[2],om,k,b_mag);
	}

	void active_stress(double tt,double (&X)[3],int row,
			const matrix &F,sym_matrix &as) {

		// inner cylinder radius
		const double Ri=     R;
		const double Li=     L;

		// area moment of inertia
		//const double I = M_PI * Ri*Ri*Ri*Ri / 4;

		// get body-centered coord
		double x[3]; for(int i=0;i<3;i++) x[i]=X[i]-c[i];

		// body-coord sys unit vectors
		double *x_hat = orientation;
		double *y_hat = bending_dir;
		double  z_hat[3] = {x_hat[1]*y_hat[2] - y_hat[1]*x_hat[2],
                            x_hat[2]*y_hat[0] - y_hat[2]*x_hat[0],
                            x_hat[0]*y_hat[1] - y_hat[0]*x_hat[1]};

		// body-axis in deformed coords
		// ref derivative index is column in deformation gradient
		double xd_hat[3] = {0,0,0};
		for (int r=0;r<3;r++)
			for (int c=0;c<3;c++)
				xd_hat[r] += F(r,c) * x_hat[c];

		// normalize
		double mag_sq=0; for(int i=0;i<3;i++) mag_sq += xd_hat[i]*xd_hat[i];
		for(int i=0;i<3;i++) xd_hat[i] /= sqrt(mag_sq);

		// split into components
		double xx=0; for(int i=0;i<3;i++) xx += x[i]*x_hat[i];
		double yy=0; for(int i=0;i<3;i++) yy += x[i]*y_hat[i];
		double zz=0; for(int i=0;i<3;i++) zz += x[i]*z_hat[i];

		// axial stress magnitude
		// FIXME pull in shear modulus
		const double GG=1;        // shear modulus
		// b_mag is vert defl wrt wavelength
		// 2pi b_mag is defl wrt wavenumber
		double s = 9*M_PI*b_mag * GG*k * zz * cos(k*xx + om*tt);

		// squared radial distance from worm centerline
		double rr = sqrt(zz*zz+yy*yy);

		static const bool use_phi=true;

		if (use_phi) {

			if (phi(X) >= 0) s=0;

			if (rr > Ri) s=0;


			double wid=0.2;
			double zlo = (xx - (-0.5*L)) / wid;
			double zhi = (0.5*L - (Rhead-0.5*wid) - xx)    / wid;

			if (zlo < 0) s = 0;
			else if (zlo >=0 && zlo <= 1) s *= zlo;

			if (zhi < 0) s = 0;
			else if (zhi >=0 && zhi <= 1) s *= zhi;


		} else {

			if (rr > Ri) {

				/*
				double eps = 1 - (rr-Ri)/(R-Ri);
				if (eps < 0) eps=0;
				if (eps > 1) eps=1;

				s *= eps;
				*/
				s = 0;
			}

			if (fabs(xx) > 0.5*Li) {
				/*
				double eps = 1 - (fabs(xx) - 0.5*Li)/(0.5*(L-Li));
				if (eps < 0) eps=0;
				if (eps > 1) eps=1;

				s *= eps;
				*/
				s = 0;
			}
		}
		
		// "basis" stress tensor is o(x)o (outer product of orientation with itself)
		// actual is that multipled by b * zz
//		for(int c=0;c<3;c++)
//			sr[c] = s*(xd_hat[row]*xd_hat[c] - (row==c?1.:0.)/3.);
		for(int r=0;r<3;r++)
			for(int c=r;c<3;c++)
				as(r,c) = s*(xd_hat[r]*xd_hat[c] - (r==c?1.:0.)/3);
	}
};

class sheet: public object{
	public:
	double t, L; // thickness and side length (sqaure cross section)
	cube(const double *basic_specs, const double *extra_specs, bool v=false) :
		object(basic_specs), L(primary_dim),
        t(extra_specs[0]){
            volume = t*L*L;
            half_side_lengths[0] = t*0.5;
            half_side_lengths[1] = L*0.5;
            half_side_lengths[2] = L*0.5;
            if(v) print_self();
    }
	~sheet() {}
    // level-set function
	virtual double phi(const double (&xi)[3]);
    // reference map values in the initial configuration
	virtual void rm0(const double (&x)[3],double (&xi)[3]) {
		// NOTE
		// F and F inv maps vector, while xi, x
		// are coordinate. Subtract off the center
		// to get the vector
        for(int i=0;i<3;i++) xi[i] = x[i];
	}
	virtual void print_self(){
		printf("# Sheet:\n"
               "#	[Center (%f %f %f), side length %f, thickness %f]\n"
			   c[0], c[1], c[2], L, t);
	}

    protected:
    double half_side_lengths [3];
};

class cube: public object{
	public:
	double L;
    int stretch_dir;
	double sqrt_lam;
    double eigvals[3];
	cube(const double *basic_specs, const double *extra_specs, bool v=false) :
		object(basic_specs), L(primary_dim),
        stretch_dir(int(extra_specs[0])),
        sqrt_lam(extra_specs[1]),
        eigvals{sqrt_lam, sqrt_lam, sqrt_lam},
        hl(L*0.5){
            volume = pow(L, 3);
            eigvals[stretch_dir] = 1./(sqrt_lam*sqrt_lam);
            if(v) print_self();
    }
	~cube() {}
	virtual double phi(const double (&xi)[3]);
	virtual void rm0(const double (&x)[3],double (&xi)[3]) {
		// NOTE
		// F and F inv maps vector, while xi, x
		// are coordinate. Subtract off the center
		// to get the vector
        for(int i=0;i<3;i++) xi[i] = (x[i]-c[i]) * eigvals[i] + c[i];
	}
	virtual void print_self(){
		printf("# Cube:\n"
               "#	[Center (%f %f %f), side length %f]\n"
		       "#	[Initial stretch in dir %d: %f]\n"
               "#	[Initial spin omega (%f, %f, %f)]\n",
			   c[0], c[1], c[2], L, stretch_dir, sqrt_lam*sqrt_lam,
               omega0[0], omega0[1], omega0[2]);
	}

    protected:
    double hl;
};

class beam: public object{
	public:
	double W;
    double L;
    int direction;
	double th_max;
    double omega;
	beam(const double *basic_specs, const double *extra_specs, bool v=false) :
		object(basic_specs), W(primary_dim),
        L(extra_specs[0]), direction(int(extra_specs[1])),
		th_max(0.25*M_PI), omega(extra_specs[9]),
        R(W*0.5), piv_r(0.75*R), hl(L*0.5){
            volume = L*W*W;
            for (int i=0;i<3;i++) {
                chi[i] = c[i];
                clow[i] = c[i];
            }
            chi[direction] = c[direction] + hl - R;
            clow[direction] = c[direction] - hl + R;
            if(v) print_self();
    }
	~beam() {}
	virtual double phi(const double (&xi)[3]);
	virtual void print_self(){
		printf("# Beam:\n"
               "#	[Center (%f %f %f), width %f, length %f]\n"
               "#	[LOW (%f %f %f) HI (%f %f %f)]\n"
               "#	[Omega %f, theta max %f]\n"
			   , c[0], c[1], c[2], W, L, clow[0], clow[1], clow[2], chi[0], chi[1], chi[2],
			   omega, th_max);
	}
    virtual void passive_acc_deformation(const double  (&rval)[3], const double eps, double x, double y, double z, double t, double &dx, double &dy, double &dz, double & dtoe1, double & dtoe2);

    protected:
    double R;
    double piv_r;
    double hl;
    double chi[3];
    double clow[3];
};


class twist_rod: public rod{
    public:
    double omega;
    const double th_max;
    const double piv_r;
    rod* cross_rod1;
    rod* cross_rod2;
    twist_rod (const double *basic_specs, const double * extra_specs, bool v=false) :
        rod(basic_specs, extra_specs), omega(1.), th_max(0.08*M_PI), piv_r(primary_dim*0.75) {
        if(v){
            print_self();
        }

        if(extra_specs[9] != 0.0) {
            omega=extra_specs[9];
        }

        // Making sure that if nothing is specified, we get an orthogonal bar in xy plane.
        double * bp_tmp = new double [4];
        double * ep_tmp = new double [10];
        bp_tmp[3] = R;
        for(int i=0;i<10;i++) ep_tmp[i] = 0.0;

        // Let's place the first cross bar on top, length is twice the width
        ep_tmp[0] = 4*R;
        ep_tmp[1] = 1;
        ep_tmp[2] = 0;
        ep_tmp[3] = 0;

        for(int i=0;i<3;i++) {
            // places the center
            bp_tmp[i] = chi[i];
        }
        cross_rod1 = new rod(bp_tmp, ep_tmp);

        // Let's place the second cross bar at the bottom
        for(int i=0;i<3;i++) {
            // places the center
            bp_tmp[i] = clow[i];
        }
        cross_rod2 = new rod(bp_tmp, ep_tmp);
    }

    virtual ~twist_rod();
    virtual double phi(const double (&xi)[3]);
    virtual void print_self(){
           printf("# twist_rod: \n"
          "#	[Center (%f %f %f), radius %f]\n"
          "#	[anchor radius %f, theta_max %f, omega %g]\n"
          "#	[Rod 1 orientation (%g %g %g) L %f]\n"
           , c[0], c[1], c[2], R,
           piv_r, th_max, omega,
           orientation[0], orientation[1], orientation[2], L);
    }

    virtual void passive_acc_deformation(const double  (&rval)[3], const double eps, double x, double y, double z, double t, double &dx, double &dy, double &dz, double & dtoe1, double & dtoe2);

    private:
    double rod1(const double (&xi)[3]);
    double rod2(const double (&xi)[3]);
};

class rotor: public rod{
    public:
    double omega;
    const double th_max;
    const double piv_r;
    rod* cross_rod;
    rotor(const double *basic_specs, const double * extra_specs, bool v=false) :
        rod(basic_specs, extra_specs), omega(1.), th_max(M_PI), piv_r(primary_dim*0.75) {
        // Making sure that if nothing is specified, we get an orthogonal bar in xy plane.
        double norm =0.;
        double * ep_tmp = new double [10];
        for(int i=0;i<10;i++) ep_tmp[i] = extra_specs[i];

        for(int i=1;i<=3;i++){
            norm += extra_specs[i]*extra_specs[i];
        }
        // if the norm is zero accidental, we set it to align with x-axis
        if(norm==0.) {
            ep_tmp[1] = 1.;
            ep_tmp[2]=ep_tmp[3]=0.;
        } else {
            // we set the second rod perpendicular to the first rod
            ep_tmp[1] = -extra_specs[2]/norm;
            ep_tmp[2] = extra_specs[1]/norm;
        }

        cross_rod = new rod(basic_specs, ep_tmp);
        // This becomes an approximate volume
        // because of the doubly counted volume
        volume += cross_rod->volume - (M_PI * R * R) * R;

        if(v){
            print_self();
        }
        delete [] ep_tmp;

        if(extra_specs[9] != 0.0) {
            omega=extra_specs[9];
        }
    }
    virtual ~rotor();
    virtual double phi(const double (&xi)[3]);
	virtual void print_self(){
		printf("# Rotor: \n"
               "#	[Center (%f %f %f), radius %f]\n"
               "#	[anchor radius %f, theta_max %f, omega %g]\n"
               "#	[Rod 1 orientation (%g %g %g) L %f]\n"
               "#	[Rod 2 orientation (%g %g %g) L %f]\n"
			   , c[0], c[1], c[2], R,
                piv_r, th_max, omega,
                orientation[0], orientation[1], orientation[2], L,
                cross_rod->orientation[0],
                cross_rod->orientation[1],
                cross_rod->orientation[2], cross_rod->L);
	}
    virtual void passive_acc_deformation(const double  (&rval)[3], const double eps, double x, double y, double z, double t, double &dx, double &dy, double &dz, double & dtoe1, double & dtoe2);

    private:
    double h_rod(const double (&xi)[3]);
    double v_rod(const double (&xi)[3]);
};

// Red blood cell
class red_bc: public basic_sphere {
	public:
	// Coefficient for ellipsoid caps
	double ecc[3];
	red_bc(const double *basic_specs, const double * extra_specs, bool v=false) :
		basic_sphere(basic_specs),
		ecc{extra_specs[0], extra_specs[1], extra_specs[2]}
		{
			if(v) print_self();
		}

	~red_bc() {}
	virtual double phi(const double (&xi)[3]);
	virtual void print_self(){
		printf( "# Red blood cell:\n"
                "VOLUME NOT CORRECTLY COMPUTED!\n"
                "#	[Center (%f %f %f), radius %f]\n"
				"#	[Ellipsoid coeff: %f %f %f]\n",
				c[0], c[1], c[2], R, ecc[0], ecc[1], ecc[2]);
	}
};

class tryp: public object {
	public:
	double L;
	double lift;
	tryp (const double *basic_specs, const double *extra_specs, bool v=false) :
	// 0.06 is 1.5/25 -> posterior end to length ratio
	object(basic_specs), L(primary_dim), lift(extra_specs[0]),
	bl_scale(24.9206), bl_cutoff(19), R(L*0.06),
	anterior(c[0]+ bl_cutoff/bl_scale * L)
	{
		anchor_r = shape(0.)*0.8;
		if(v) print_self();
	}
	virtual void rm0(const double  (&x)[3],double (&xi)[3]);
	virtual double phi(const double (&xi)[3]);
    virtual void passive_acc_deformation(const double  (&rval)[3], const double eps, double x, double y, double z, double t, double &dx, double &dy, double &dz, double & dtoe1, double & dtoe2);
	virtual void print_self(){
		printf(	"# Tryp:\n"
                "VOLUME NOT CORRECTLY COMPUTED!\n"
                "#	[Posterior center (%f %f %f)]\n"
                "#	[Anterior center (%f %f %f)]\n"
				"#	[Length scale %f]\n"
				"#	[Sim length %f, width %f]\n"
				"#	[Length cut off %f]\n"
				"#	[Deflection %f]\n",
				c[0], c[1], c[2],
				anterior, c[1], c[2],
				bl_scale, L, R, bl_cutoff,lift);
	}

	private:
	const double bl_scale;
	const double bl_cutoff;
	const double R;
	const double anterior;
	double anchor_r;
	double shape(double bx);
};

class active_rod: public rod {
    public:
    active_rod (const double *basic_specs, const double * extra_specs, bool v=false) :
        rod(basic_specs, extra_specs) , omega(1./M_PI), lambda(2.2) {}
	~active_rod() {}
    virtual void actuate(const double (&rval)[3], const double (&rnval)[3], const double time, matrix &F);
    private:
    double omega;
    double lambda;
    double exponent(const double dist_h, const double dist_v, const double time, const double phiv);
};

class bent_beam: public object{
	public:
	double W;
    double L;
    int direction;
    double T;
	double lamb;
	bent_beam(const double *basic_specs, const double *extra_specs, bool v=false) :
		object(basic_specs), W(primary_dim),
        L(extra_specs[0]), direction(int(extra_specs[1])),
		T(extra_specs[2]), lamb(extra_specs[3]),
        R(W*0.5), edge(0.05), hl(L*0.5),
        init_w((1.0/lamb/lamb- lamb*lamb)/4/R),
        init_d(1/init_w/init_w * sqrt(1+ 4*init_w*init_w*R*R)),
        cr_lamb(0.565),
        cr_w((1.0/cr_lamb/cr_lamb- cr_lamb*cr_lamb)/4/R),
        cr_d(1/cr_w/cr_w * sqrt(1+ 4*cr_w*cr_w*R*R))
    {
            volume = L*W*W;

            for (int i=0;i<3;i++) {
                chi[i] = c[i];
                clow[i] = c[i];
            }
            chi[direction] = c[direction] + hl;
            clow[direction] = c[direction] - hl;

            bending_c[0] = c[0];
            bending_c[1] = c[1] - sqrt(init_d);
            bending_c[2] = c[2];
            if(v) print_self();
    }
	~bent_beam() {}
	virtual double phi(const double (&xi)[3]);
    virtual void rm0(const double (&x)[3],double (&xi)[3]);
	virtual void print_self(){
		printf("# Bending Instability Beam:\n"
               "#	[Center (%f %f %f), width %f, length %f]\n"
               "#	[LOW (%f %f %f) HI (%f %f %f)]\n"
               "#	[Bending center (%f %f %f) constant d = %.6f constant w = %.6f]\n"
               "#	[Initial stretch (compression) %f]\n"
               "#	[Initial bending angle %.3f rad %.3f degree]\n"
			   , c[0], c[1], c[2], W, L,
                clow[0], clow[1], clow[2], chi[0], chi[1], chi[2],
                bending_c[0], bending_c[1], bending_c[2], init_d, init_w,
			    lamb, L*init_w, L*init_w / M_PI*180);
	}
    virtual void passive_acc_deformation(const double  (&rval)[3], const double eps, double x, double y, double z, double t, double &dx, double &dy, double &dz, double & dtoe1, double & dtoe2);
    void compute_constants(const double X1, const double X2, const double X3, const double tau, double& r, double & theta, double & curr_d);
    double compute_d(const double curr_w);

    protected:
    double R;
    double edge;
    double hl;
    double chi[3];
    double clow[3];
    double init_w;
    double init_d;
    double cr_lamb;
    double cr_w;
    double cr_d;
    double bending_c[3];
};

#endif
