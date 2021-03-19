#ifndef F3D_MG_PROB
#define F3D_MG_PROB

#include <limits>
#include <cmath>
#include "timer.hh"
#include "multigrid.hh"
#include "common.hh"
#include <vector>

struct stencil27 {
    double st[27];
    stencil27() : st{0,0,0,0,0,0,0,0,0,
                     0,0,0,0,0,0,0,0,0,
                     0,0,0,0,0,0,0,0,0}

    {}
    void set(const double  *stencil, const int len_){
        for(int i=0;i<len_;i++) st[i] = stencil[i];
    }
    void print(){
        printf("[");
        for(int i=0;i<27;i++) {
            printf("%g ", st[i]);
        }
        printf("]\n");
    }
    void print(int len){
        printf("[");
        for(int i=0;i<len;i++) {
            printf("%g ", st[i]);
        }
        printf("]\n");
    }
};

struct mg_pred {

	/** iteration decrease period */
	int iP;
	/** iteration value */
	int iV;
	/** moving average */
	double V_avg;
	/** moving variance */
	double V_var;
	/** memory parameter i.e. 1/(number of interations to remember) */
	double alpha;

	mg_pred() : iP(16),iV(2*iP),V_avg(iV),V_var(1e6),alpha(1/(3*iP)) {}
	inline void update(int extras) {

		// do an extra reduction for each standard deviation we're
		// exceeding the average by
		double dV=iV-V_avg;
		int exreduce = static_cast<int>(dV/sqrt(V_var));
		exreduce=exreduce>0?exreduce:0;

		// update value (correct if too low)
		iV += iP*extras*(extras+1)/2 - 1 - exreduce;
		if (iV < iP) iV = 2*iP - 1;

		// update averages
		V_avg = (1-alpha)*V_avg + alpha*iV;
		V_var = (1-alpha)*(V_var + alpha*dV*dV);
	};
	inline void update(double extra_d,double d_per_iter) {

		int extra_iters = (int) ((extra_d / d_per_iter) - 0.5);
		iV -= iP*((int) extra_iters);
	};
	inline void update(int extras,double extra_d,double d_per_iter) {
		update(extras);
		update(extra_d,d_per_iter);
	};
	inline int iters() {return iV/iP;};

};

struct mg_opt {

	struct mg_config {
		int u,d;
		double dps;
	};

	bool searching;
	int *u,*d,n,c,r,ub,db,T,t;
	mg_config *states;
	const char *id;

	mg_opt() : searching(false),u(NULL),d(NULL),n(0),c(0),
		r(0),ub(2),db(2),T(-1),t(T),states(NULL),id(NULL) {};
	~mg_opt() {if (states != NULL) delete[] states;};

	inline void set() {
		if (--t == 0) local_search(3);
		if (searching) {
			*u = states[c].u;
			*d = states[c].d;
		} else {
			*u = ub;
			*d = db;
		}
	};

	inline void connect(multigrid *mg,const char* id_) {
		u = &mg->v_up;
		d = &mg->v_down;
		id = id_;
	};

	inline void local_search(int r_) {
		r = r_;
		states = new mg_config[9*r];
		for (int rr = 0; rr < r; rr++) {
		for (int uc = *u-1; uc <= *u+1; uc++) {
			for (int dc = *d-1; dc <= *d+1; dc++) {
				if (uc >= 0 && dc >= 0 && uc+dc > 0) {
					states[n].u = uc;
					states[n].d = dc;
					n++;
				}
			}
		}
		}

		searching = true;
	}

	inline void global_search(int r_) {
		r = r_;
		states = new mg_config[48*r];
		for (int rr = 0; rr < r; rr++) {
		for (int uc = 0; uc < 7; uc++) {
			for (int dc = 0; dc < 7; dc++) {
				if (uc+dc > 0) {
					states[n].u = uc;
					states[n].d = dc;
					n++;
				}
			}
		}
		}

		searching = true;
	}

	inline void finish() {
		FILE *fh = p_safe_fopen(id,"w");
		for (int i = 0; i < n; i++) fprintf(fh,"%2d %2d %8.5g\n",
			states[i].u,states[i].d,states[i].dps);
		fclose(fh);
		delete[] states;
		states = NULL;
		c = 0;
		n = 0;
		r = 0;
		searching = false;
		t = T;
	}

	inline void schedule(int T_) {
		if (T_ < 0) {
			global_search(3);
		} else if (T_ > 0) {
			t = T = T_;
			local_search(3);
		}
	};

	inline void record(double dps_) {
		states[c++].dps = dps_;
		if (c==n) find_best();
	};

	inline void find_best() {

		double *mins = new double[n/r];
		int *inds = new int[n/r];
		for (int i = 0; i < n/r; i++) {
			mins[i] = states[i].dps;
			inds[i] = i;
			for (int rr = 0; rr < r; rr++) {
				if (mins[i] > states[i + (n/r)*rr].dps) {
					mins[i] = states[i + (n/r)*rr].dps;
					inds[i] = i + (n/r)*rr;
				}
			}
		}

		double max = mins[0];
		ub = states[inds[0]].u;
		db = states[inds[0]].d;

		for (int i = 0; i < n/r; i++) {
			if (max < mins[i]) {
				max = mins[i];
				ub = states[inds[i]].u;
				db = states[inds[i]].d;
			}
		}

		delete[] mins;
		delete[] inds;
		finish();
	};
};

struct stencil_set {

    // TODO variable spatial dicretization
	double main[27];
	double **partial;
    // For convenience, we also store the length of each partial stencil
    int partial_len[27];
	bool use_partials,xprd,yprd,zprd;
	int m,n,o;

	stencil_set(geometry *gm) : main{
                     0,0,0,0,0,0,0,0,0,
                     0,0,0,0,0,0,0,0,0,
                     0,0,0,0,0,0,0,0,0},
        partial(NULL),
		use_partials(false),xprd(gm->x_prd),yprd(gm->y_prd),
		zprd(gm->z_prd),m(gm->m),n(gm->n),o(gm->o) { };
	~stencil_set() {
		if (partial!=NULL) {
			for (int i = 0; i < 27; i++) delete[] partial[i];
			delete[] partial;
		}
	};

	inline double *get(int i,int j,int k, int & len) {
        len = use_partials?partial_len[reg(i,j,k)]:27;
		return use_partials?partial[reg(i,j,k)]:main;
	};

	inline void print(int i,int j,int k, int len) {
        len = use_partials?partial_len[reg(i,j,k)]:27;
		if(use_partials) {
            printf("[");
            for(int s=0;s<len;s++) printf("%g ", partial[reg(i,j,k)][s]);
            printf("]\n");
        } else {
            printf("[");
            for(int s=0;s<len;s++) printf("%g ", main[s]);
            printf("]\n");
        }
	};


	inline int reg(int i,int j,int k) {
		int si,sj,sk;
		si = xprd?1:((i==0)?0:((i==m-1)?2:1));
		sj = yprd?1:((j==0)?0:((j==n-1)?2:1));
		sk = zprd?1:((k==0)?0:((k==o-1)?2:1));
		return si + 3*sj + 9*sk;
	};

    // Get the individual index in the region of the domain
	inline int reg(int i, int j, int k, int &si,int &sj,int &sk) {
		si = xprd?1:((i==0)?0:((i==m-1)?2:1));
		sj = yprd?1:((j==0)?0:((j==n-1)?2:1));
		sk = zprd?1:((k==0)?0:((k==o-1)?2:1));
		return si + 3*sj + 9*sk;
	};

	inline double fac(int i,int low,int num) {
		if (num == 3) {
			return 1;
		} else { // num ==2
			return (i!=low)?2:1;
		}
	};

	// FIXME hacky additional boolean argument to tell partial
	// stencil generator to subtract missing "face" entries from
	// central element, representing homogenous Neumann BC on
	// MAC projection solve
	void set_partials(bool centers,bool macfix=false) {
//		prnt = true;
		double **p = partial = new double*[27];
		double *pi;
		double tfac; // total factor
		int ni,nj,nk;
		int li,lj,lk;
		for (int k = 0; k < 3; k++)
			for (int j = 0; j < 3; j++)
				for (int i = 0; i < 3; i++) {
			ni = (i==1)?3:2;
			nj = (j==1)?3:2;
			nk = (k==1)?3:2;

			li = (i==0)?1:0;
			lj = (j==0)?1:0;
			lk = (k==0)?1:0;

            partial_len[i+j*3+k*9] = ni*nj*nk;
            //printf("set partial length (%d %d %d) = %d\n", i, j, k, partial_len[i+j*3+k*9]);
            //printf("set partial [");
            //printf("{");
			pi = *(p++) = new double[ni*nj*nk];
			for (int kk = 0; kk < nk; kk++)
				for (int jj = 0; jj < nj; jj++)
					for (int ii = 0; ii < ni; ii++) {

				tfac = centers?1:(fac(ii,li,ni)*fac(jj,lj,nj)*fac(kk,lk,nk));
				*(pi) = main[(li+ii) + 3*(lj+jj) + 9*(lk+kk)]/tfac;

				// if doing MAC tweak and at central element
				if (macfix && (li+ii)==1 && (lj+jj)==1 && (lk+kk)==1) {

					// add "face" element to center element for
					// each direction corresponding to a ghost node
					int ns[3] = {ni,nj,nk}; // n nodes
					int ls[3] = {li,lj,lk}; // first index into main array
					int st[3] = {1,3,9};    // strides
					int ind = (li+ii) + 3*(lj+jj) + 9*(lk+kk); // current main index

					// for each direction...
					for (int d=0;d<3;d++) {

						// if on a wall...
						if (ns[d] < 3) {

							// add left or right stencil as appropriate
							// if ls[d]==1, we're "below" (e.g. on left wall)
							// if ls[d]==0, we're "above" (e.g. on right wall)
							(*pi) += main[ind + (ls[d]==1?-1:1)*st[d]];
						}
					}

				}

                //printf("%g, ", fabs(*pi));
                pi++;
			}
            //for(int ii=ni*nj*nk; ii<27;ii++) printf("0, ");
            //puts("}");
		}
		use_partials = true;
	};
};

/** A base class with basic function to define a simple problem for multigrid*/
class mg_prob {
public:
	int channels;
	stencil_set st;
	double tol;
	double log_tol;
    double max_norm_tol;
	mg_opt *mgo;
	mg_pred *mgp;
	geometry *gm;
	char (*out)[256];
	const char **id;
    double dx, dy, dz;
	timer sw;

	mg_prob(geometry *gm_, int ch_) :
        channels(ch_),st(gm_),
		mgo(new mg_opt[channels]),
		mgp(new mg_pred[channels]),gm(gm_),
		out(new char[channels][256]),
		id(new const char*[channels]),
        sw("MG_PROBS", channels){
            //puts("Base mg_prob");
        };
	virtual ~mg_prob() {
        //puts("Base destructor mg_prob");
		delete[] mgo;
		delete[] mgp;
		delete[] id;
		delete[] out;
	};

    // Definition of stencils, can be changed by children classes
    virtual void set_grid_spacings(const double dx_, const double dy_, const double dz_){
        dx = dx_;
        dy = dy_;
        dz = dz_;
    }

	virtual bool internal(int i,int j,int k) {
        return false;
    }

	virtual double *external_ptr(int i,int j,int k) {
        int len;
		return st.get(i,j,k,len);
	}

	virtual int mem_size(int i,int j,int k) {
		p_fatal_error("mg_prob should never call mem_size",1);
		return 0;
	}

	virtual void fill_entries(int i,int j,int k) {
		p_fatal_error("mg_prob should never call fill_entries",1);
	}

    // TODO variable spacial discretization
    // Pure
    virtual void setup() = 0;

    // virtual functions that have to do with variable densities
    // A pure virtual one for now
    // TODO make Pure
    virtual void set_rhoinv(const double *rho) = 0;

    // Pure
    virtual void set_stencils() = 0;

    inline bool interior(int i,int j,int k) {
        return (gm->x_prd||(i>0&&i<gm->m-1)) &&
            (gm->y_prd||(j>0&&j<gm->n-1)) &&
            (gm->z_prd||(k>0&&k<gm->o-1));
    }
	inline int range_xd(int i,int j,int k) {return (gm->x_prd||i>0)?-1:0;}
    inline int range_xu(int i,int j,int k) {return (gm->x_prd||i<(gm->m-1))?2:1;}
    inline int range_yd(int i,int j,int k) {return (gm->y_prd||j>0)?-1:0;}
    inline int range_yu(int i,int j,int k) {return (gm->y_prd||j<(gm->n-1))?2:1;}
    inline int range_zd(int i,int j,int k) {return (gm->z_prd||k>0)?-1:0;}
    inline int range_zu(int i,int j,int k) {return (gm->z_prd||k<(gm->o-1))?2:1;}

	inline double x_field(int i,int j,int k) {
		p_fatal_error("mg_prob should never call x_field",1);
		return 0;
	}

	inline double r_field(int i,int j,int k) {
		p_fatal_error("mg_prob should never call r_field",1);
		return 0;
	}

	inline void connect(multigrid *mg) {
		for (int c = 0; c < channels; c++) {
			mgo[c].connect(mg,id[c]);
		}
	}

	inline double resid(multigrid *mg) {
		//double resid = mg->l2_error_all() / (gm->m*gm->n*gm->o);
		//MPI_Bcast(&resid,1,MPI_DOUBLE,0,gm->cart);
		double resid = mg->l0_error_all();
        /*
		if (std::isnan(resid)) p_fatal_error(
			"NaN detected in multigrid solve",MG3D_MATH_ERROR);
        */
		return resid;
	}

	double run(multigrid *mg, bool profile, int c=0);
};

class simp_fem: public mg_prob {
    public:
        simp_fem(geometry *gm_) : mg_prob(gm_, 1) {
            if(gm->rank==0) puts("# Creating simple FEM projection MG problem.");
        }
        ~simp_fem() {
            //puts("Destructor simple fem");
        }
        virtual void setup() {

            // even are -1, odd are -2
            for (int i = 0; i < 27; i++) st.main[i] = -1-(i%2);

            // "faces" are 0
            st.main[4] = st.main[22] = st.main[10] =
                st.main[16] = st.main[12] = st.main[14] = 0;

            // central node is 32
            st.main[13] = 32;

            // set up boundary stencils if not periodic in all dirs
            if (!(gm->x_prd && gm->y_prd && gm->z_prd))
                st.set_partials(false);

            // tol scales with central node (32)
            double eps = std::numeric_limits<double>::epsilon();
            tol = 32*32*eps*eps*1e4;
            max_norm_tol = 32*eps*1e2;
            //max_norm_tol = 1e-8;
            log_tol = log10(tol);

            // id strings
            *id = "proj";

        }
        virtual void set_rhoinv(const double *rho) {};
        virtual void set_stencils() {};
};

class simp_mac: public mg_prob {
    public:
        simp_mac(geometry *gm_) : mg_prob(gm_, 1) {
            if(gm->rank==0) puts("# Creating simple MAC projection MG problem.");
        }
        ~simp_mac() {
            //puts("Destructor simp mac");
        }
        virtual void setup() {
            for (int i = 0; i < 27; i++) st.main[i] = 0;

            st.main[4]  = st.main[22] = 1;
            st.main[10] = st.main[16] = 1;
            st.main[12] = st.main[14] = 1;

            // central node is 1 plus sum
            st.main[13] = -6;

            // tol scales with central node (6)
            double eps = std::numeric_limits<double>::epsilon();
            tol =6*6*eps*eps*1e4;
            max_norm_tol = 6*eps*1e2;
            //max_norm_tol = 1e-8;
            log_tol = log10(tol);

            *id = "mac";
            if (!(gm->x_prd && gm->y_prd && gm->z_prd)){
                st.set_partials(true,true);
            }

        }
        virtual void set_rhoinv(const double *rho) {};
        virtual void set_stencils() {};
};

class simp_visc: public mg_prob {
    public:
        simp_visc(geometry *gm_) : mg_prob(gm_, 3) {
            if(gm->rank==0) puts("# Creating simple implicit viscosity MG problem.");
        }
        ~simp_visc() {
            //puts("Desctructor simple visc");
        }
        virtual void setup(){

            // start with everything zero
            for (int i = 0; i < 27; i++) st.main[i] = 0;

            // "faces" take provided prefactors
            st.main[4] = st.main[22] = -lz;
            st.main[10] = st.main[16] = -ly;
            st.main[12] = st.main[14] = -lx;

            // central node is 1 plus sum
            st.main[13] = 1 + 2*(lx+ly+lz);

            // tol scales with central node (1)
            double eps = std::numeric_limits<double>::epsilon();
            tol = eps*eps*1e4;
            max_norm_tol = eps*1e2;
            //max_norm_tol = 1e-8;
            log_tol = log10(tol);

            *id = "u_visc";
            id[1] = "v_visc";
            id[2] = "w_visc";
            if (!(gm->x_prd && gm->y_prd && gm->z_prd)){
                st.set_partials(true);
            }
        }
        virtual void set_rhoinv(const double *rho) {};
        virtual void set_stencils() {};
        void setup_const(const double lx_, const double ly_, const double lz_){
            lx=lx_;
            ly=ly_;
            lz=lz_;
        }
    private:
        double lx, ly, lz;
};

class var_den_mg_prob: public mg_prob {
    public:
        int sm, sn, so;
        int cellm4, celln4, cello4;
        int smn;
        int cellmn4;
        int c0;
        double * rho_inv;
        stencil27 * all_stencils;
        int *   stencil_len;

        var_den_mg_prob(geometry *gm_, const int cellm4_, const int celln4_, const int cello4_):
        mg_prob(gm_,1),
        sm(gm->sm), sn(gm->sn), so(gm->so),
        cellm4(cellm4_), celln4(celln4_), cello4(cello4_),
        smn(sm*sn),
        cellmn4 (cellm4*celln4),
        c0 (2+ 2*cellm4 +2*cellmn4),
        rho_inv(new double[cellm4 * celln4 * cello4]),
        all_stencils(new stencil27[sm*sn*so])
        {
             for(int i=0;i<cellm4 * celln4 * cello4;i++) rho_inv[i]=1;
            //puts("var den mg prob base");
        }
        virtual ~var_den_mg_prob(){
            //puts("Destructor var den mg prob base");
            delete [] all_stencils;
            delete [] rho_inv;
        }
        virtual double * external_ptr(int i, int j, int k);
        virtual void set_rhoinv(const double *rho);
        // Pure
        virtual void set_stencils() = 0;

};

/** Inherit from abstract class var_den_mg_prob */
class var_den_fem: public var_den_mg_prob{
    public:
        static int const num_rho_contrib [27][27];
        var_den_fem(geometry *gm_, const int cellm4_, const int celln4_, const int cello4_) :
        var_den_mg_prob(gm_, cellm4_, celln4_, cello4_)
        {
            if(gm->rank==0) puts("# Creating variable density FEM projection MG problem.");
        }
        ~ var_den_fem() {
            //puts("Destructor var den mg prob fem");
        }
        virtual void setup(){
            // even are -1, odd are -2
            for (int i = 0; i < 27; i++) st.main[i] = -1-(i%2);

            // "faces" are 0
            st.main[4] = st.main[22] = st.main[10] =
                st.main[16] = st.main[12] = st.main[14] = 0;

            // central node is 32
            st.main[13] = 32;

            // set up boundary stencils if not periodic in all dirs
            if (!(gm->x_prd && gm->y_prd && gm->z_prd))
                st.set_partials(false);

            // tol scales with central node (32)
            double eps = std::numeric_limits<double>::epsilon();
            tol = 32*32*eps*eps*1e4;
            max_norm_tol = 32*eps*1e2;
            //max_norm_tol = 1e-8;
            log_tol = log10(tol);

            // id strings
            *id = "proj";
        }
        virtual void set_stencils();
};

/** Inherit from abstract class var_den_mg_prob */
class var_den_mac: public var_den_mg_prob{
    public:
        static int const num_rho_contrib [27][27];
        var_den_mac(geometry *gm_, const int cellm4_, const int celln4_, const int cello4_) :
        var_den_mg_prob(gm_, cellm4_, celln4_, cello4_)
        {
            if(gm->rank==0) puts("# Creating variable density MAC projection MG problem.");
        }
        ~ var_den_mac() {
            //puts("Destructor var den mg prob mac");
        }
        virtual void setup(){
            for (int i = 0; i < 27; i++) st.main[i] = 0;

            st.main[4]  = st.main[22] = 1;
            st.main[10] = st.main[16] = 1;
            st.main[12] = st.main[14] = 1;

            // central node is 1 plus sum
            st.main[13] = -6;

            // tol scales with central node (6)
            double eps = std::numeric_limits<double>::epsilon();
            tol =6*6*eps*eps*1e4;
            max_norm_tol = 6*eps*1e2;
            //max_norm_tol = 1e-8;
            log_tol = log10(tol);

            *id = "mac";
            if (!(gm->x_prd && gm->y_prd && gm->z_prd)){
                st.set_partials(true,true);
            }

        }
        virtual void set_stencils();
};

#endif
