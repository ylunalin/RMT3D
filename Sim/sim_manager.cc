#include "sim_manager.hh"
#include <cstring>
#include <algorithm>

sim_manager::~sim_manager() {
	for (int i = 0; i < n_obj; i++) delete objs[i];
    if (avg_velx != NULL) delete[] avg_velx;
    if (avg_vely != NULL) delete[] avg_vely;
    if (avg_velz != NULL) delete[] avg_velz;
    if (avg_velN != NULL) delete[] avg_velN;
    if (avg_x != NULL) delete[] avg_x;
    if (avg_y != NULL) delete[] avg_y;
    if (avg_z != NULL) delete[] avg_z;
    if (avg_N != NULL) delete[] avg_N;
	if (objs != NULL) delete[] objs;
    if(sm_array!=NULL) delete [] sm_array;
}

double sim_manager::cfl(double dx) {

	// the minimum time scale allowed by shear wave of the solid
	double xxsp3 = 3.0/(dx*dx);
	double dt_s = 100000;
	for (int i = 0; i < n_obj; i++) {
		double odt = sm_array[i].cfl(dx);
		if (odt < dt_s) dt_s = odt;
	}
	// the time scale set by fluid viscosity
	double dt_f = 0.5* fm.rho/(fm.mu*xxsp3)*fm.dt_pad;

	return (dt_s<dt_f)?dt_s:dt_f;
}

/** Set up
 * CFL max dt and transition zone width.
 */
void sim_manager::setup_const(){
	double dh=spars->min_dh();
	// max dt allowed by cfl condition
	dt_reg = cfl(dh);
	// transition zone widht
	eps = wt_n*dh;
	eps_inv = 1./eps;
	tf1 =0.5/eps;
	tf2 =0.5/M_PI;
	tf3 =M_PI/eps;

    // User provided time step could be smaller
    if(spars->dt < dt_reg) dt_reg = spars->dt;

    // Set the anchoring acceleration
    K_stiff = 0.0025/(dt_reg*dt_reg);
}

/** This is to check for advective CFL condition after fluid velocity has been initialized */
void sim_manager::obligatory_cfl_recheck(){
        double dh=spars->min_dh();
        double dt_min = fm.dt_pad*dh/vmax;
        if(dt_reg>dt_min) {
            dt_reg = dt_min;
            // User provided time step could be smaller
            if(spars->dt < dt_reg) dt_reg = spars->dt;
        }
        // Set the anchoring acceleration
        K_stiff = 0.0025/(dt_reg*dt_reg);
}

void sim_manager::create_objs(){
    for(int i=0;i<n_obj;i++){
        sm_array[i] = sl_mat(spars, i);
        // HACK HACK HACK:
        sm_array[i].ex_mu = sm_array[i].ex_visc_mult*fm.mu;

        double *basics = spars->basic_specs + n_basic_specs*i;
        double *extras = spars->extra_specs + n_extra_specs*i;

        int obj_type_index = spars->object_list[i];

        objs [i] = object::alloc(obj_type_index, basics, extras);
    }
}

void sim_manager::add_obj(object *o) {
	// new array, copy old stuff over
	object **oarr = new object*[n_obj+1];
	for (int i = 0; i < n_obj; i++) oarr[i] = objs[i];

	// delete old array (if it exists)
	if (objs != NULL) delete[] objs;

	// assign new array to our pointer, set new val
	// and increment number of objects
	objs = oarr;
	objs[n_obj++] = o;
}

double sim_manager::phi(int obj_index, const double (&xi)[3]) const{
    if(obj_index<0 || obj_index >= n_obj) return nan("");
	return objs[obj_index]->phi(xi);
}

void sim_manager::velocity(double x,double y,double z,
	double &u,double &v,double &w) {
    u=0; v=0; w=0;
    const double coords[3] = {x,y,z};
    double rx[3] = {0,0,0};

    double * vpp = spars->vel_profile;
    for(int i=0;i<spars->vel_prof_num;i++){
        int vtype = static_cast<int> (vpp[2*i]);

        double amp = vpp[2*i+1];
        if(vtype == 0) continue;
        else if (vtype == 1) u+= amp;
        else if (vtype == 2) v+= amp;
        else if (vtype == 3) w+= amp;
        else {

        }
    }

    for(int i=0;i<n_obj;i++) {
        objs[i]->rm0(coords, rx);
        if(objs[i]->phi(rx) < 0.) {
            objs[i]->velocity(rx, u, v, w);
        }
    }

    double tmp = sqrt(u*u+v*v+w*w);
    if(tmp>vmax) vmax = tmp;
}

void sim_manager::solid_acceleration(const int obj_id, const double (&rval)[3], double x, double y, double z, double t,	double &fx,double &fy,double &fz){
//void sim_manager::solid_acceleration(const int obj_id, const double (&rval)[3], double x, double y, double z, double t,	double &fx,double &fy,double &fz, const double vx, const double vy, const double vz){
    fx=fy=fz=0;
}

// FLUID CONVERGENCE
void sim_ftest::velocity(double x,double y,double z,
	double &u,double &v,double &w) {
	u =     - sin(2*M_PI*x) * cos(2*M_PI*y) * cos(2*M_PI*z);
	v = 0.5 * cos(2*M_PI*x) * sin(2*M_PI*y) * cos(2*M_PI*z);
	w = 0.5 * cos(2*M_PI*x) * cos(2*M_PI*y) * sin(2*M_PI*z);
    double tmp = sqrt(u*u+v*v+w*w);
    if(tmp>vmax) vmax = tmp;
}

void sim_ftest::pressure(double x,double y,double z,double &p) {
	p = 3*M_PI*fm.mu*cos(2*M_PI*x)*cos(2*M_PI*y)*cos(2*M_PI*z);
}

void sim_ftest::exact(double x,double y,double z,double t,
	double h,double &u,double &v,double &w,double &p) {

	// get exact vels at cell center
	velocity(x,y,z,u,v,w);
	u *= cos(2*M_PI*t);
	v *= cos(2*M_PI*t);
	w *= cos(2*M_PI*t);

	// shift to corner for pressure
	x -= 0.5*h;
	y -= 0.5*h;
	z -= 0.5*h;
	pressure(x,y,z,p);
	p *= (cos(2*M_PI*t) - sin(2*M_PI*t)/(6*M_PI*fm.mu));
}

void sim_ftest::fluid_acceleration(double x,double y,double z,double t,
	double &fx,double &fy,double &fz) {
	fx = (12*M_PI*cos(2*M_PI*y)*cos(2*M_PI*z)*(-6*M_PI*fm.mu*cos(2*M_PI*t) +
            sin(2*M_PI*t))*sin(2*M_PI*x) + M_PI*cos(2*M_PI*t)*cos(2*M_PI*t) *
            (2+cos(4*M_PI*y) + cos(4*M_PI*z))*sin(4*M_PI*x))/4;
    fy = -(M_PI*cos(2*M_PI*t)*cos(2*M_PI*t) *
            (-1+cos(4*M_PI*x) -2*cos(4*M_PI*z))*sin(4*M_PI*y))/8;
    fz = -(M_PI*cos(2*M_PI*t)*cos(2*M_PI*t) *
            (-1+cos(4*M_PI*x) -2*cos(4*M_PI*y))*sin(4*M_PI*z))/8;
}

// SOLID CONVERGENCE / SHEAR WAVE
sim_stest::sim_stest(const sim_params *spars) : sim_manager(spars),
	a(0.01),k(2*M_PI),w(2*M_PI) {
	shear_cube *cube = new shear_cube(a,k);
	add_obj(cube);
}

void sim_stest::exact(double x,double y,double z,double t,
	double h,double &u,double &v,double &w,double &p) {
	u = this->w*a*sin(k*z - w*t);
	v = this->w*a*sin(k*z - w*t);
	w = 0;

	x -= 0.5*h;
	y -= 0.5*h;
	z -= 0.5*h;
	p = -sm_array[0].G * a*a*k * (-k*cos(2*(k*z-this->w*t)) + cos(2*this->w*t - k)*sin(k))/3;
}

void sim_stest::pressure(double x,double y,double z,double &p) {
	p = -sm_array[0].G * a*a*k * (-k*cos(2*k*z) + cos(-k)*sin(k))/3;
}

void sim_stest::velocity(double x,double y,double z,
	double &u,double &v,double &w) {
	u = this->w*a*sin(k*z);
	v = this->w*a*sin(k*z);
	w = 0;
    double tmp = sqrt(u*u+v*v+w*w);
    if(tmp>vmax) vmax = tmp;
}

void sim_objects::solid_acceleration(const int obj_id, const double (&rval)[3], double x, double y, double z, double t,	double &fx,double &fy,double &fz){
        //double k_rep_default = wall_acc_default ;
        double krepx, krepy, krepz;
        krepx = wall_acc_mult*avg_velx[obj_id]*avg_velx[obj_id];
        krepy = wall_acc_mult*avg_vely[obj_id]*avg_vely[obj_id];
        krepz = wall_acc_mult*avg_velz[obj_id]*avg_velz[obj_id];

        if(krepx<min_wall_acc) krepx = min_wall_acc;
        if(krepy<min_wall_acc) krepy = min_wall_acc;
        if(krepz<min_wall_acc) krepz = min_wall_acc;

        // We use fx,fy,fz, but they are really accelerations
        fx=fy=fz=0;
        double phiv = phi(obj_id, rval);
        // wall forces in x direction
        if(walls[0]){
            double dist = x-wall_pos[0];
            if(dist<=wall_dist){
                fx += delta_func_wall(dist*winv) * krepx *wall_trans_func(phiv,dx,dxsp);
            }
        }
        if(walls[1]){
            double dist = wall_pos[1]-x;
            if(dist<=wall_dist){
                fx -= delta_func_wall(dist*winv) * krepx *wall_trans_func(phiv,dx,dxsp);
            }
        }

        // wall forces in y direction
        if(walls[2]){
            double dist = y-wall_pos[2];
            if(dist<=wall_dist){
                fy += delta_func_wall(dist*winv) * krepy *wall_trans_func(phiv,dy,dysp);
            }
        }
        if(walls[3]){
            double dist = wall_pos[3]-y;
            if(dist<=wall_dist){
                fy -= delta_func_wall(dist*winv) * krepy *wall_trans_func(phiv,dy,dysp);
            }
        }

        // wall forces in z direction
        if(walls[4]){
            double dist = z-wall_pos[4];
            if(dist<=wall_dist) {
                fz += delta_func_wall(dist*winv) * krepz *wall_trans_func(phiv,dz,dzsp);
            }
#if 0
            // HERE TO HARD CODE A WALL AT 1/2 slope, i.e. z = 0.5*x
            double dist = z-0.5*x;
            // normal distance from a point to the plane that has a 1/2 slope
            double norm_dist = dist * 0.894;
            if(norm_dist<=wall_dist){
                double tmp_force = delta_func_wall(dist*winv) * krepz *wall_trans_func(phiv,dz,dzsp);
                fz += tmp_force * 0.8;
                fx -= tmp_force * 0.4;
            }
#endif
        }
        if(walls[5]){
            double dist = wall_pos[5]-z;
            if(dist<=wall_dist){
                fz -= delta_func_wall(dist*winv) * krepz *wall_trans_func(phiv,dz,dzsp);
            }
        }
        double f_norm = sqrt(fx*fx + fy*fy + fz*fz);
        if(f_norm >= K_stiff ){
            double normfac = K_stiff/f_norm;
            fx *= normfac;
            fy *= normfac;
            fz *= normfac;
        }
        // Sometimes the anchors edge need to be specified by 2 distances
        // e.g. a cylindrical anchor region, R, L
        double ldx=0, ldy=0, ldz=0, dtoe1=2*eps, dtoe2=2*eps;
        objs[obj_id]->passive_acc_deformation(rval, eps, x, y, z, t, ldx, ldy, ldz, dtoe1, dtoe2);
        double K = K_stiff*heaviside(dtoe1)*heaviside(dtoe2);
        fx-=K*ldx;
        fy-=K*ldy;
        fz-=K*ldz;

        // Add gravitational force (buoyancy force)
        //fz += gravity*heaviside(phiv) * (sm_array[obj_id].rho - fm.rho) / sm_array[obj_id].rho;
        fz += gravity*heaviside(phiv);
}

// FLUID CONVERGENCE
void object_mills::velocity(double x,double y,double z,
	double &u,double &v,double &w) {
	u =     - sin(2*M_PI*x) * cos(2*M_PI*y) * cos(2*M_PI*z);
	v = 0.5 * cos(2*M_PI*x) * sin(2*M_PI*y) * cos(2*M_PI*z);
	w = 0.5 * cos(2*M_PI*x) * cos(2*M_PI*y) * sin(2*M_PI*z);
    double tmp = sqrt(u*u+v*v+w*w);
    if(tmp>vmax) vmax = tmp;
}

void object_mills::pressure(double x,double y,double z,double &p) {
	p = 3*M_PI*fm.mu*cos(2*M_PI*x)*cos(2*M_PI*y)*cos(2*M_PI*z);
}

void object_mills::exact(double x,double y,double z,double t,
	double h,double &u,double &v,double &w,double &p) {

	// get exact vels at cell center
	velocity(x,y,z,u,v,w);
	u *= cos(2*M_PI*t);
	v *= cos(2*M_PI*t);
	w *= cos(2*M_PI*t);

	// shift to corner for pressure
	x -= 0.5*h;
	y -= 0.5*h;
	z -= 0.5*h;
	pressure(x,y,z,p);
	p *= (cos(2*M_PI*t) - sin(2*M_PI*t)/(6*M_PI*fm.mu));
}

void object_mills::fluid_acceleration(double x,double y,double z,double t,
	double &fx,double &fy,double &fz) {
	fx = (12*M_PI*cos(2*M_PI*y)*cos(2*M_PI*z)*(-6*M_PI*fm.mu*cos(2*M_PI*t) +
            sin(2*M_PI*t))*sin(2*M_PI*x) + M_PI*cos(2*M_PI*t)*cos(2*M_PI*t) *
            (2+cos(4*M_PI*y) + cos(4*M_PI*z))*sin(4*M_PI*x))/4;
    fy = -(M_PI*cos(2*M_PI*t)*cos(2*M_PI*t) *
            (-1+cos(4*M_PI*x) -2*cos(4*M_PI*z))*sin(4*M_PI*y))/8;
    fz = -(M_PI*cos(2*M_PI*t)*cos(2*M_PI*t) *
            (-1+cos(4*M_PI*x) -2*cos(4*M_PI*y))*sin(4*M_PI*z))/8;
}

