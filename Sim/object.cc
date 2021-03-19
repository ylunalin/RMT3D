#include "object.hh"

object* object::alloc(const int type_name, const double *bp, const double *ep){
    if(type_name == BASIC_SPHERE){
        return new basic_sphere(bp, false);

    } else if(type_name == SPHERE){
        return new sphere(bp, ep, false);

    } else if( type_name == ROD){
        return new rod(bp, ep, false);

    } else if( type_name == CUBE){
        return new cube(bp, ep, false);

    } else if( type_name == ROTOR){
        return new rotor(bp, ep, false);

    } else if( type_name == RED_BLOOD_CELL){
        return new red_bc(bp, ep, false);

    } else if( type_name == TRYP){
        return new tryp(bp, ep, false);

    } else if( type_name == ACTIVE_ROD){
        return new active_rod(bp, ep, false);

    } else if( type_name == TWIST_ROD){
        return new twist_rod(bp, ep, false);

    } else if( type_name == BEAM){
        return new beam(bp, ep, false);
    } else if( type_name == ACTIVE_ROD_AS){
        return new active_rod_as(bp, ep, false);
    } else if( type_name == BENTBEAM){
        return new bent_beam(bp, ep, false);
    } else {
        p_fatal_error("object:: Unknown object type!\n", 1);
        return NULL;
    }
}

const double object::fuzzy_eps = 1./128.;

double object::fuzzy_max(const double d1, const double d2){
    double diff = fabs(d1-d2);
    if(diff > fuzzy_eps) return std::max(d1,d2);
    else{
        double addon = (fuzzy_eps-diff)*(fuzzy_eps-diff);
        double ld1 = d1+addon, ld2 = d2+addon;
        return std::max(ld1,ld2);
    }
}

double object::fuzzy_min(const double d1, const double d2){
    double diff = fabs(d1-d2);
    if(diff > fuzzy_eps) return std::min(d1,d2);
    else{
        double takeoff = fuzzy_eps-diff;
        double ld1 = d1-takeoff, ld2 = d2-takeoff;
        return std::min(ld1,ld2);
    }
}

double shear_cube::phi(const double (&xi)[3]) {
	// for periodic cube, everything is always inside
	return -1;
}

void shear_cube::rm0(const double (&x)[3],double (&xi)[3]) {
	xi[0] = x[0] - a*cos(k*x[2]);
	xi[1] = x[1] - a*cos(k*x[2]);
	xi[2] = x[2];
}

double basic_sphere::phi(const double (&xi)[3]) {
	double sq_dist = 0;
	for (int i = 0; i < 3; i++)
		sq_dist += (xi[i]-c[i])*(xi[i]-c[i]);
    //printf("From object (%g %g %g), c (%g %g %g) sq_dist %g, phi %g\n", xi[0], xi[1], xi[2], c[0], c[1], c[2], sq_dist, sqrt(sq_dist)-R);
	return sqrt(sq_dist) - R;
}

double rod::phi(const double (&xi)[3]) {
	// cap centers
    double pt_vec[3] {0,0,0}, dot_product=0, norm=0;
    for(int i=0;i<3;i++) {
        pt_vec [i]  = xi[i] - c[i];
        norm += pt_vec[i]*pt_vec[i];
        dot_product += pt_vec[i] * orientation[i];
    }

	// get dist from each
	double d_cyl=0, d_low=0, d_hi=0;
    for(int i =0;i<3;i++) {
        d_low += (xi[i] - clow[i])*(xi[i] - clow[i]);
        d_hi  += (xi[i] - chi[i])*(xi[i] - chi[i]);
    }
    d_cyl = sqrt(norm - dot_product*dot_product) - R;
    d_low = sqrt(d_low)-R;
    d_hi = sqrt(d_hi)-R;

    // Test which side of the rod we are on
    double tmp_phi = 0;
	if(fabs(dot_product)>L*0.5) {
        tmp_phi = std::min(d_low, d_hi);
    } else {
        tmp_phi = d_cyl;
    }
    return tmp_phi;

}

double cube::phi(const double (&xi)[3]) {
    double d[3];
	for(int i=0;i<3;i++){
		d[i] = fabs(xi[i]-c[i])-hl;
	}
    double fuzzy_xy = fuzzy_max(d[0], d[1]);
    double fuzzy_xz = fuzzy_max(d[0], d[2]);
    double fuzzy_yz = fuzzy_max(d[1], d[2]);
    double max_phi = fuzzy_xy>fuzzy_xz?fuzzy_xy:fuzzy_xz;
    max_phi = fuzzy_yz>max_phi?fuzzy_yz:max_phi;

    return max_phi;
}

rotor::~rotor(){
    delete cross_rod;
}

double rotor::phi(const double (&xi)[3]) {
    double hphi = h_rod(xi);
    double vphi = v_rod(xi);

    double min_phi = fuzzy_min(hphi, vphi);
    return min_phi;
}

double rotor::h_rod(const double (&xi)[3]) {
    return cross_rod->phi(xi);
}

double rotor::v_rod(const double (&xi)[3]) {
    return rod::phi(xi);
}

// return strain, since elastic constants are set in the sim_manager class
void rotor::passive_acc_deformation(const double  (&rval)[3], const double eps, double x, double y, double z, double t, double &dx, double &dy, double &dz, double & dtoe1, double & dtoe2){
    dtoe2 = -2*eps;
    const double delx = x-c[0],  dely = y-c[1], delz = z-c[2];

    const double vecX = rval[0]-c[0],  vecY = rval[1]-c[1];
    const double r=sqrt(delx*delx + dely*dely + delz*delz);
    if (r < piv_r + eps) {
        double theta=th_max*(1-cos(t*omega))*(omega>0?1:-1);

        // rotate at constant angular velocity
        //double theta = fmod(t*omega, 2*M_PI)*(omega>0?1:-1);

        double cth=cos(theta),sth=sin(theta);
        // rotate the reference vector to the position it would have been
        // if it were rotating about the center rigidly
        double rx=c[0]+vecX*cth-vecY*sth, ry=c[1]+vecX*sth+vecY*cth;
        dtoe1 =  r-piv_r;
        dx = x-rx;
        dy = y-ry;
        // This is to fix it in place
        dz = z-rval[2];
    }
}

double red_bc::phi(const double (&xi)[3]){
	double phi_val=0;
	for (int i=0;i<3;i++) phi_val += ecc[i]*(xi[i]-c[i])*(xi[i]-c[i]);
	phi_val = sqrt(phi_val)-R;
	return phi_val;
}

void tryp::rm0(const double  (&x)[3],double (&xi)[3]) {
	xi[0] = x[0];
	xi[1] = x[1];
	xi[2] = x[2];
	if(lift!=0.) xi[2] = -(x[0] - c[0])*lift + x[2];
}

double tryp::phi(const double (&xi)[3]){
	double bx = (xi[0]-c[0])/L * bl_scale;
	double dist;
	double bw = shape(bx);
	if(bx<=0) {
		dist=0;
		for(int i=0;i<3;i++) dist+= (xi[i]-c[i])*(xi[i]-c[i]);
	} else if(bx<=bl_cutoff) {
		dist=0;
		for(int i=1;i<3;i++) dist+= (xi[i]-c[i])*(xi[i]-c[i]);
	} else {
		dist = (xi[0]-anterior)*(xi[0]-anterior);
		for(int i=1;i<3;i++) dist+= (xi[i]-c[i])*(xi[i]-c[i]);
	}
	dist = sqrt(dist)-bw;
	return dist;
}

double tryp::shape(double bx){
	if(bx<0) bx=0;
	else if(bx>bl_cutoff) bx=bl_cutoff;
	double bx2 = bx*bx;
	double bx3 = bx2*bx, bx4 = bx2*bx2, bx5=bx2*bx3;
	double bw = 0.460185*bx-0.0700506*bx2+0.00425402*bx3-0.000125024*bx4+1.47545e-06*bx5+0.29601;
	return bw*R;
}

void tryp::passive_acc_deformation(const double  (&rval)[3], const double eps, double x, double y, double z, double t, double &dx, double &dy, double &dz, double & dtoe1, double & dtoe2){
    const double delx = x-c[0],  dely = y-c[1], delz = z-c[2];
    const double r=sqrt(delx*delx+dely*dely+delz*delz);
	if(r<anchor_r+eps){
		dtoe1 = r - anchor_r;
		dx = x-rval[0]; dy = y-rval[1]; dz=z-rval[2];
	}
}

void active_rod::actuate(const double (&rval)[3], const double (&rnval)[3], const double time, matrix &F){
    // average reference map
    double avg_rval[3];
    for(int i=0;i<3;i++){
        avg_rval[i] = (rval[i] + rnval[i])/2;
    }

    double phiv = phi(avg_rval);
    // get the reference map vector to the head
    double rel_rval[3];
    for(int i=0;i<3;i++){
        rel_rval[i] = avg_rval[i] - clow[i];
    }

    /*
    double dist_h = fabs(avg_rval[0]*orientation[0] +
                    avg_rval[1]*orientation[1] +
                    avg_rval[2]*orientation[2])
    */
    // HACK HACK HACK
    // We assume the rod is horizontal and in an XZ plane.

    if(!(orientation[0]==1 && orientation[1]==0 && orientation[2]==0)) {
        //printf("orientation %g %g %g\n", orientation[0], orientation[1],orientation[2]);
        return;
    }
    double dist_h = rel_rval[0]/L;
    double dist_v = rel_rval[2]/R;
    double alpha = exponent(dist_h, dist_v, time, phiv);
#if 0
    printf("ref map (%g %g %g) \n"
            "relative to clow (%g %g %g)\n"
            "dist_h %g dist_v %g alpha %g\n",
            avg_rval[0], avg_rval[1], avg_rval[2],
            rel_rval[0], rel_rval[1], rel_rval[2],
            dist_h, dist_v, alpha);
#endif

    // Multiplied the inverse of Fa into F
    F(0,0) *= exp(alpha);
    F(2,2) *= exp(-alpha);
}

double active_rod::exponent(const double dist_h, const double dist_v, const double time, const double phiv) {
    if(dist_h < 0 || dist_h>1 || fabs(dist_v) > 0.8 || phiv>-0.001){return 0;}
    else {
        //return -lambda*dist_v*pow(sin(omega*time), 7)*sin( 2*M_PI*dist_h);
        //return -lambda*dist_v*sin(omega*time)*sin( 2*M_PI*dist_h);
        return -lambda*dist_v*sin(omega*time);
    }

}

void twist_rod::passive_acc_deformation(const double  (&rval)[3], const double eps, double x, double y, double z, double t, double &dx, double &dy, double &dz, double & dtoe1, double & dtoe2) {
    dtoe2 = -2*eps;
    const double delx_low = x-clow[0],  dely_low = y-clow[1], delz_low = z-clow[2];
    const double delx_hi = x-chi[0],  dely_hi = y-chi[1], delz_hi = z-chi[2];

    const double rlow = sqrt(delx_low*delx_low + dely_low*dely_low + delz_low*delz_low);
    const double rhi = sqrt(delx_hi*delx_hi + dely_hi*dely_hi + delz_hi*delz_hi);
    const int sign_omg = omega>0?1:-1;
    /*
    if (rhi < piv_r + eps || rlow < piv_r + eps) {
	if(rlow < piv_r + eps) {
	    double theta=th_max*(1-cos(t*omega))*sign_omg;
	    // rotate at constant angular velocity
	    double cth=cos(theta),sth=sin(theta);

            const double vecX = rval[0]-clow[0],  vecY = rval[1]-clow[1];
            // rotate the reference vector to the position it would have been
            // if it were rotating about the center rigidly
            double rx=clow[0]+vecX*cth-vecY*sth, ry=clow[1]+vecX*sth+vecY*cth;
            dtoe1 =  rlow-piv_r;
            dx = x-rx;
            dy = y-ry;
            // This is to fix it in place
            dz = z-rval[2];
        }
	if(rhi < piv_r +eps) {
	    // rotate 360 degrees then rotate back
	    double theta=-th_max*(1-cos(t*omega))*sign_omg;
	    // rotate at constant angular velocity
	    double cth=cos(theta),sth=sin(theta);

            const double vecX = rval[0]-chi[0],  vecY = rval[1]-chi[1];
            // rotate the reference vector to the position it would have been
            // if it were rotating about the center rigidly
            double rx=chi[0]+vecX*cth-vecY*sth, ry=chi[1]+vecX*sth+vecY*cth;
            dtoe1 =  rhi-piv_r;
            dx = x-rx;
            dy = y-ry;
            // This is to fix it in place
            dz = z-rval[2];
        }
    }
    */
    if (rhi < piv_r + eps || rlow < piv_r + eps) {
        if (rlow < piv_r + eps) {
            dtoe1 =  rlow-piv_r;
            // This is to fix it in place
            dx = x-rval[0];
            dy = y-rval[1];
            dz = z-rval[2];
        }
        if (rhi < piv_r +eps) {
            double theta=-th_max*(1-cos(t*omega))*sign_omg;
            double cth=cos(theta),sth=sin(theta);

            // rotate about the low anchor
            const double vecZ = rval[2]-clow[2],  vecY = rval[1]-clow[1];
            // rotate the reference vector to the position it would have been
            // if it were rotating about the center rigidly
            double ry=clow[1]+vecY*cth-vecZ*sth, rz=clow[2]+vecY*sth+vecZ*cth;
            dtoe1 =  rhi-piv_r;
            dx = x-rval[0];
            dy = y-ry;
            // This is to fix it in place
            dz = z-rz;
        }
    }
}

double twist_rod::phi(const double (&xi)[3]) {
    double phi0 = rod::phi(xi);
    double phi1 = rod1(xi);
    double phi2 = rod2(xi);

    double min_phi = fuzzy_min(phi0, phi1);
    min_phi = fuzzy_min(min_phi, phi2);
    return min_phi;
}

double twist_rod::rod1(const double (&xi)[3]) {
    return cross_rod1->phi(xi);
}

double twist_rod::rod2(const double (&xi)[3]) {
    return cross_rod2->phi(xi);
}

twist_rod::~twist_rod() {
    delete cross_rod1;
    delete cross_rod2;
}

double beam::phi(const double (&xi)[3]) {
    double d[3];
	for(int i=0;i<3;i++){
		d[i] = fabs(xi[i]-c[i])-R;
	}
    d[direction] = fabs(xi[direction] - c[direction]) - hl;
    double fuzzy_xy = fuzzy_max(d[0], d[1]);
    double fuzzy_xz = fuzzy_max(d[0], d[2]);
    double fuzzy_yz = fuzzy_max(d[1], d[2]);
    double max_phi = fuzzy_xy>fuzzy_xz?fuzzy_xy:fuzzy_xz;
    max_phi = fuzzy_yz>max_phi?fuzzy_yz:max_phi;

    return max_phi;
}

void beam::passive_acc_deformation(const double  (&rval)[3], const double eps, double x, double y, double z, double t, double &dx, double &dy, double &dz, double & dtoe1, double & dtoe2) {
    dtoe2 = -2*eps;
    const double delx_low = x-clow[0],  dely_low = y-clow[1], delz_low = z-clow[2];
    const double delx_hi = x-chi[0],  dely_hi = y-chi[1], delz_hi = z-chi[2];

    const double rlow = sqrt(delx_low*delx_low + dely_low*dely_low + delz_low*delz_low);
    const double rhi = sqrt(delx_hi*delx_hi + dely_hi*dely_hi + delz_hi*delz_hi);
    const int sign_omg = omega>0?1:-1;
    if (rhi < piv_r + eps || rlow < piv_r + eps) {
	if(rlow < piv_r + eps) {
	    double theta=th_max*(1-cos(t*omega))*sign_omg;
	    // rotate at constant angular velocity
	    double cth=cos(theta),sth=sin(theta);

            const double vecX = rval[0]-clow[0],  vecY = rval[1]-clow[1];
            // rotate the reference vector to the position it would have been
            // if it were rotating about the center rigidly
            double rx=clow[0]+vecX*cth-vecY*sth, ry=clow[1]+vecX*sth+vecY*cth;
            dtoe1 =  rlow-piv_r;
            dx = x-rx;
            dy = y-ry;
            // This is to fix it in place
            dz = z-rval[2];
        }
	if(rhi < piv_r +eps) {
	    double theta=-th_max*(1-cos(t*omega))*sign_omg;
	    // rotate at constant angular velocity
	    double cth=cos(theta),sth=sin(theta);

            const double vecX = rval[0]-chi[0],  vecY = rval[1]-chi[1];
            // rotate the reference vector to the position it would have been
            // if it were rotating about the center rigidly
            double rx=chi[0]+vecX*cth-vecY*sth, ry=chi[1]+vecX*sth+vecY*cth;
            dtoe1 =  rhi-piv_r;
            dx = x-rx;
            dy = y-ry;
            // This is to fix it in place
            dz = z-rval[2];
        }
    }
}

void bent_beam::passive_acc_deformation(const double  (&rval)[3], const double eps, double x, double y, double z, double t, double &dx, double &dy, double &dz, double & dtoe1, double & dtoe2) {
    dtoe1 = -2*eps;
    dtoe2 = -2*eps;

    double X = rval[0] - c[0];
    double theta = X * init_w;
    double xfinal, yfinal;

    //printf("X %g theta %g, slope %g\n", X, theta * 180 / M_PI, cos(theta)/sin(theta));
    if(theta!=0) {
        double m = cos(theta) / sin(theta);
        if(m!=0) {
            double invm = 1./m;
            xfinal = m * bending_c[0] + x * invm + y - bending_c[1];
            xfinal /= (m+invm);
            yfinal = m * (x-bending_c[0]) + bending_c[1];

        } else {
            xfinal = x; yfinal = bending_c[1];
        }
    } else {
            xfinal = bending_c[0]; yfinal = y;
    }

    dx = x-xfinal;
    dy = y-yfinal;
    dz = z-rval[2];
}

void bent_beam::rm0(const double (&x)[3], double (&xi)[3]) {
    double currx = x[0]-bending_c[0], curry = x[1]-bending_c[1];
    if(curry==0) return;
    double theta = atan2(currx,curry);
    xi[0] = theta/init_w + c[0];

    double Y = ((currx*currx + curry*curry) - init_d) * init_w * 0.5;
    xi[1] = Y + c[1];
    xi[2] = x[2];
}

void bent_beam::compute_constants(const double X, const double Y, const double Z, const double tau, double& r, double& theta, double& curr_d){
    // we should have checked that X1 and X2 are within bounds
    double curr_w = init_w + tau * (cr_w - init_w);
    curr_d = compute_d(curr_w);
    r = sqrt(curr_d + 2*Y/curr_w);
    theta = curr_w * X;
}

double bent_beam::compute_d(const double curr_w){
    return 1/curr_w/curr_w * sqrt(1+ curr_w*curr_w*W*W);
}

double bent_beam::phi(const double (&xi)[3]){
    double d[3];
	for(int i=0;i<3;i++){
		d[i] = fabs(xi[i]-c[i])-R;
	}

    d[direction] = fabs(xi[direction] - c[direction]) - hl;
    double fuzzy_xy = fuzzy_max(d[0], d[1]);
    double fuzzy_xz = fuzzy_max(d[0], d[2]);
    double fuzzy_yz = fuzzy_max(d[1], d[2]);
    double max_phi = fuzzy_xy>fuzzy_xz?fuzzy_xy:fuzzy_xz;
    max_phi = fuzzy_yz>max_phi?fuzzy_yz:max_phi;

    return max_phi;
}
