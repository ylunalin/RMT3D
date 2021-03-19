#ifndef F3D_EXTRAP_HELPER_
#define F3D_EXTRAP_HELPER_
#include <cstdlib>
#include <cstdio>
#include <vector>
#include <algorithm>
#include <cmath>
#include <climits>
#include <mpi.h>

#include "geometry.hh"
#include "comm_info.hh"
#include "buffer.hh"
#include "mat3.hh"
#include "timer.hh"
#include "sim_params.hh"

// A helper class for extrapolation, builds a linear model for a given object.
// This is used for each object, one by one, controlled by fluid_3d.cc.
// This class has its own communication functions to handle the boundaries.
// This class is templated on extrapolation dimesnion N.
template<int N>
class extrap_helper{
	public:
    const static int cube_num = 125;
	const static int D=3;
    const static int data_num = cube_num*D;
	const static int nmat = (D+1)*D/2;
	const static int XX = 0;
	const static int YY = 1;
	const static int ZZ = 2;
	const static int XY = 3;
	const static int XZ = 4;
	const static int YZ = 5;
	const static int rlen = D*N;

	// geometry information
	const int sm, sn, so;
	const int sm4, sn4, so4;
	const int smn4, smno4, smno;
	// defining top 4 corners of the 3x3x3 cube, used in neighbor rules
	const int corner0;
	const int corner1;
	const int corner2;
	const int corner3;
	// the current highest integer value in mask array
	unsigned int nmask;
	// the starting index for the current object
	unsigned int cur_mask;
	// current object being extrapolated
	int cur_obj;
	// the number of points in the layer for an object in this domain
	// guarantee < smno (total num of grid points)
	int num;
	int primary_num;
	// MPI rank
	int rank;
    // number of surrounding points for each index
    int capacity;
	// the total weights at a given extrapolation point
	double wt_capacity;
    // the averages of data and rhs
    double data_avg [D];
    double rval_avg [N];

	// the indices of grid points to extrapolate
	int * indices;
	// neighbor id
	int * neigh;
	// Mask grid to indicate the interior of an object
	unsigned int * mask;
	// extrapolation weights
	const double wt;
	// the coefficients for each spatial coordinate in linear model
	double * coeff;
	// contains x,y,z weighted components of ref map, 3xcapacity elements
	double * data;
	// the weight for each data points
	double * weight_data;
	// the reference map values (x, [y], [z]) correponding to indices array
	double * rval;
	// A list of grid indices for the next layer, temporary storage
	std::vector<int> list;
	// pointer to parent geometry object
	geometry *grid;
	// mpi communication request and status objects
	MPI_Request *reqs;
	MPI_Status *stats;
	// communication table
	comm_info comm_table[27];
	// communication buffer
	comm_buffer cbuf;
    // A timer
    timer watch;

	extrap_helper ():
		sm(0), sn(0), so(0), sm4(0), sn4(0), so4(0),
		smn4(0), smno4(0), smno(0),
		corner0(0), corner1(0), corner2(0), corner3(0),
		nmask(1), cur_mask(1), cur_obj(0), num(0),
		primary_num(0), rank(0), capacity(0), wt_capacity(0.),
		indices(NULL), neigh(NULL), mask(NULL),
		wt(0.5), coeff(NULL), data(NULL), weight_data(NULL),
		rval(NULL), grid(NULL), reqs(NULL), stats(NULL), watch("Extrap_Helper_Lite", 0)
		{}

	extrap_helper (geometry &gm, const double wt_, const sim_params *sp):
		sm(gm.sm), sn(gm.sn), so(gm.so), sm4(sm+4), sn4(sn+4), so4(so+4),
		smn4(sm4*sn4), smno4(smn4*so4), smno(sm*sn*so),
		corner0(smn4+sm4-1), corner1(smn4+sm4+1),
		corner2(smn4-sm4-1), corner3(smn4-sm4+1),
		nmask(1), cur_mask(1), cur_obj(0), num(0), primary_num(0), rank(gm.rank),
        capacity(0), wt_capacity(0.),
        // TODO capacity - need 1
        // coeff - 3*N
        // data maxi  125*3
        // weight_data max 125
        // rval max  125*N
        // wt_capacity = capacity?
        indices(NULL), neigh(NULL), mask(NULL),
        wt(wt_), coeff(NULL), data(NULL), weight_data(NULL),
        rval(NULL), grid(&gm), reqs(NULL), stats(NULL),
        watch("Extrap_Helper_Lite", 5, "Add_layer", "Least_square", "Reset", "Communication", "Data_transfer")
		{

            indices=new int [smno4];
            mask=new unsigned int[smno4];
            coeff=new double [3*N];
            data=new double [375];
            weight_data=new double [125];
            rval=new double [125*N];

            for(int i=0;i<smno4;i++){
                indices[i]=0;
                mask[i]=0;
            }
            for(int i=0;i<3*N;i++) coeff[i]=0;
            for(long i=0;i<data_num;i++) data[i]=0;
            for(long i=0;i<cube_num; i++) weight_data[i]=0;
            for(long i=0;i<cube_num*N;i++) rval[i]=0;

            setup_comm();
		}

	~extrap_helper (){
		if(stats!=NULL) delete [] stats;
		if(reqs!=NULL) delete [] reqs;
		if(neigh!=NULL) delete [] neigh;
		if(rval!=NULL) delete [] rval;
		if(weight_data!=NULL) delete [] weight_data;
		if(data!=NULL) delete [] data;
		if(coeff!=NULL) delete [] coeff;
		if(mask!=NULL) delete [] mask;
		if(indices!=NULL) delete [] indices;
	}

	/** Use to add spatial coordiante data into data arrray. */
	void pack_data(const double (&lcoord)[3], const double data_wt);
	/** Use to add right hand side value into rval array. */
	void pack_rval(const double (&r)[N]);
	void add_capacity(const double data_wt);
	void add_first_layer(const int rule);
	bool least_square(bool v=false);
	void update_list(bool verbose=false);
	/** Add an index into the list vector, according to a rule. 0-orthogonal neighbor, etc. */
	inline bool add_to_next_list(int ind, int p, int q, int r, const int rule, bool verbose=false){
		int step = index(p, q, r);
		int abs_step = index(abs(p), abs(q), abs(r));
		int candidate = ind+step;
		// if the candidate position we are considering is in the ghost region,
		// some other processor will take care of it.
		// if(!within_domain(candidate)) return;

		// here we choose the rule of having orthogonal neighbors for now;
		// later we can change it so that it can between different rule
		if(rule == 0){
			if(orthogonal(step, candidate)){
				list.push_back(candidate);
				mask[candidate]=nmask;
                int ii,jj,kk;
                unpack_index(ind, ii, jj, kk);
				return true;
			}
			else return false;
		}
		else if(rule==1){
			if(sqrt2nbs(step, abs_step, candidate, verbose)){
				list.push_back(candidate);
				mask[candidate]=nmask;
				return true;
			}
			else return false;
		}
		else if(rule==2){
			if(sqrt3nbs(abs_step, candidate, verbose)){
				list.push_back(candidate);
				mask[candidate]=nmask;
				return true;
			}
			else return false;
		}
		else{
			printf("Extrap_helper::add_to_next_list(): Not a rule.\n");
			exit(1);
			return false;
		}
	}
	/** Mark the inside of the current object, set starting index */
	inline void mark_primary(int ind){
		mask[ind]= nmask;
		//if(within_domain(ind)) primary_num++;
	}
	/** Reset nmask when it's close to overflow unsigned integer. */
	void reset(int oid){
        watch.tic(2);
		cur_obj = oid;
		// if the starting mask integer is within a certain range of the unsigned int max
		// we reset the mask array to all zeros, and start from nmask=1
		if(nmask > UINT_MAX-50){
			nmask=1;
			for(int i=0;i<smno4;i++) mask[i]=0;
		}
		// let cur_mask (mask in the primary object) be set by current nmask
		cur_mask = nmask;
        primary_num=0;
        watch.toc(2);
	}
	inline void increase_nmask(){
		nmask++;
	}
	/** Communication functions. */
	void comm_mask();
	void scan_ghost_region();

	//.private:

	/** Communication functions. */
	void setup_table();
	void setup_comm();
	void copy_to_buf(int s_tag, unsigned int *(&b));
	void copy_from_buf(int r_tag, unsigned int *(&b));

	// functions that defines the neighbor cells

	// a rule to find orthogonal extrapolation neighbors
	inline bool orthogonal(int step, int candidate){
		//int ii, jj, kk;
		//unpack_index(candidate-step, ii, jj, kk);
		//if(step == 1 && ii==32 && jj==16 && kk==15) printf("\n\n\n\nexaming in (%d %d %d)'s neighbor (%d) mask %d cur_mask %d\n", ii, jj, kk, candidate, mask[candidate-step], cur_mask);
		int diff = abs(step);
		if(( diff==1 || diff==sm4 || diff==smn4) && (mask[candidate] < cur_mask) ){
			return true;
		}
		else return false;
	}

	// a rule to find all the neighbors within sqrt(2)
	inline bool sqrt2nbs(int step, int abs_step, int candidate, bool verbose =false){
		int diff = abs(step);
		if ( (mask[candidate] < cur_mask)
			&& sqrt3nbs(abs_step, candidate)
			&& diff!=corner0
			&& diff!=corner1
			&& diff!=corner2
			&& diff!=corner3) return true;
		else return false;
	}

	// a rule to find all the neighbors within sqrt(3)
	inline bool sqrt3nbs(int abs_step, int candidate, bool verbose =false){
		int i, j, k;
		unpack_index(abs_step, i,j, k);
		//printf("step %d (%d %d %d), mask here %u, cur_mask %u\n", abs_step, i, j, k, mask[candidate], cur_mask);
		i = abs(i); j=abs(j), k=abs(k);
		if(i>1) i = sm4-i;
		if(j>1) j = sn4-j;
		if(k>1) k = so4-k;
		if( (mask[candidate] < cur_mask) && i<2 && j<2 && k<2) {
			return true;
		}
		else return false;
	}
	inline int index(int i, int j, int k){
		return i+sm4*j+smn4*k;
	}
	inline void unpack_index(int ind, int &i, int &j, int &k){
		k = ind/smn4;
		j = (ind % smn4) / sm4;
		i = ind%sm4;
	}
	// We extrapolate everywhere except
    // 1. ghost regions for processor boundaries
    // 2. ghost regions for domain boundaries in periodic boundary conditions
	// The index used used is pos w.r.t. to first ghost node
	inline bool within_extrapolation_domain(int ind){
		int i, j, k;
        bool within = true;
		unpack_index(ind, i, j, k);
        int right_edge = sm4-3, back_edge = sn4-3, top_edge = so4-3;
        // On the left face
        if(i<2){
            if(j<2) {
                // if there's neighbor procs at the vertical thin long ghost strip
                // then it's either processor boundary or periodic boundary
                // don't extrapolate
                if(k>=2 && k<= top_edge) {
                    if (neigh[9] >=0) within=false;
                } else if (k<2) {
                    if(neigh[0] >=0) within=false;
                } else {
                    if(neigh[18]>=0) within=false;
                }
            } else if(j>back_edge){
                if(k>=2 && k<= top_edge) {
                    if(neigh[15]>=0) within=false;
                } else if (k<2) {
                    if(neigh[6] >=0) within=false;
                } else {
                    if(neigh[24]>=0) within=false;
                }
            } else {
                if(k<2) {
                    if(neigh[3]>=0) within=false;
                } else if(k>top_edge) {
                    if(neigh[21] >=0) within=false;
                } else {
                    if(neigh[12]>=0) within=false;
                }
            }
        } else if(i>right_edge) {
            if(j<2) {
                // if there's neighbor procs at the vertical thin long ghost strip
                // then it's either processor boundary or periodic boundary
                // don't extrapolate
                if(k>=2 && k<= top_edge) {
                    if(neigh[11] >=0) within=false;
                } else if (k<2) {
                    if(neigh[2] >=0) within=false;
                } else {
                    if(neigh[20]>=0) within=false;
                }
            } else if(j>back_edge){
                if(k>=2 && k<= top_edge){
                    if(neigh[17]>=0) within=false;
                } else if(k<2) {
                    if(neigh[8]>=0) within=false;
                } else {
                    if(neigh[26]>0) within=false;
                }
            } else {
                if(k<2) {
                    if(neigh[5]>=0) within=false;
                } else if(k>top_edge) {
                    if(neigh[23] >=0) within=false;
                } else {
                    if(neigh[14]>=0) within=false;
                }
            }

        }

        if(j<2){
            if(k<2) {
                if(i>=2 && i<= right_edge && neigh[1]>=0) {
                    within=false;
                }
            } else if(k>top_edge){
                if(i>=2 && i<= right_edge && neigh[19]>=0) {
                    within=false;
                }
            } else{
                if(i>=2 && i<= right_edge && neigh[10]>=0) {
                    within=false;
                }
            }
        } else if(j>back_edge){
            if(k<2) {
                if(i>=2 && i<= right_edge && neigh[7]>=0) {
                    within=false;
                }
            } else if(k>top_edge){
                if(i>=2 && i<= right_edge && neigh[25]>=0) {
                    within=false;
                }
            } else{
                if(i>=2 && i<= right_edge && neigh[16]>=0) {
                    within=false;
                }
            }

        }

        if(k<2) {
            if(i>=2 && i<= right_edge && j>=2 && j<= back_edge){
                if(neigh[4]>=0) within=false;
            }
        } else if (k>top_edge) {
            if(i>=2 && i<= right_edge && j>=2 && j<= back_edge){
                if(neigh[22]>=0) within=false;
            }
        }

		return within;
	}
	// check if a grid point is within the inner ghost region
	// this is used for getting next layers spawn from points might be in this layer
	inline bool within_inner_ghost(int ind){
		int i, j, k;
		unpack_index(ind, i, j, k);
		if(i>=1 && i<= sm4-2 && j>=1 && j<= sn4-2 && k>=1 && k<=so4-2) return true;
		else return false;
	}
    // zero out capacity and wt_capacity
    inline void zero_capacity(){
        capacity=0;
        wt_capacity=0.;
    }

	void unique_indices();
	template<int n>
	void add_avg(double (&avg)[n], const double *p, const double *w, bool verbose=false);
	template<int n>
	void write_avg(double (avg)[n], double *p);
	void lin_sys(double (&melems)[6], double (&rhs)[3*N], const double (&r_avg)[N], const double (&arb_data_avg)[3], const double *rp, const double *dp, const double *wp, bool verbose=false);
	void dot(mat3 A, double (&vec)[3*N]);
	void pack_coeff(const double (&tmp_coeff)[rlen]);
	inline static int max_of_three(int x1, int x2, int x3){
		int max =(x1>x2)?x1:x2;
		max = (max>x3)?max:x3;
		return max;
	}

    inline void subtract_data_avg(const double x, const double y, const double z, double (&out)[D]){
        out[0] = x-data_avg[0];
        out[1] = y-data_avg[1];
        out[2] = z-data_avg[2];
    }

    inline void get_rval_avg(double (&out)[N]){
        for(int i=0;i<N;i++) out[i] = rval_avg[i];
    }
	void dump_mask()	;
};

#endif
