/** A header file for the sim_params class, which contains
 * simulation parameters for the bac_sim class. Part of it
 * is adapted from fileinfo class in Chris's code. The * general structure of the code is similar to modularized dld.
 *
 * Author	: Y Luna Lin
 * Email	 : y.lin2783@gmail.com
 * Date	  : Nov 11, 2018
 */

#ifndef SIM_PARAMS_HH
#define SIM_PARAMS_HH

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "omp.h"
#include "defs.hh"
#include "common.hh"
using namespace OBJECT_DICT;

class sim_params{
	public:
	/** Parameters for the overall system. */
	bool suppress;
	bool display;
	bool impl, godunov, x_prd, y_prd, z_prd, has_exact;
    bool walls[6];
	static const int buf_size = 1024;
    int chk_step;
	int chk_num;
	int sys_size[3];
	int frames;
    int stop_short;
    int chkpt_freq;
	int omp_num_thr;
	int debug_flag;
	int debug_obj_id;
	int obj_body;
	int data_fmt;
	// Add any slice outputs
	// 1 velocity u_x (0)
	// 2 velocity u_y (1)
	// 4 velocity u_z (2)
	// 8 vorticity w_x (3)
	// 16 vorticity w_y (4)
	// 32 vorticity w_z (5)
	// 64 vorticity magnitude (6)
	// 128 pressure (7)
	// 256 level set phi (8)
	// 512 reference map xi_x (9)
	// 1024 reference map xi_y (10)
	// 2048 reference map xi_z (11)
	int out_flag;
	// output dimension, 0-x, 1-y, 2-z
	int output_dim;
	// the index of the output slice
	int output_ind;
	// The flag to indicate what type of output to write
	// Every bit is used
	// first bit  => write tracers
	// second bit => write contours
	// third bit  => write slices of various fields
	// fourth bit => write checkpoint files for every frame
    // fifth bit => write diagnostic file, turn on DEBUG macro to work
    // sixth bit => write macro.dat
	unsigned int dump_code;
	// number of tracers to use
	int ntracers;
    // number of initial iterations
    int num_iters;
	double dt, T;
	double cur_time;
	double ax,ay,az,bx,by,bz;
	double lx,ly,lz;
    double wall_pos[6];
    // Distance to the wall, used to calculate the repulsive force
    // from the wall. sim_manager class calculates 5*dx as the minimal
    // so anything lower than that will not be taken by the simulation
    double wall_dist;
    double dx, dy, dz;
	// Output direcotry
	char * dirname;
	char * chk_dirname;
	char * recover_dirname;

	/** Params for fluid. */
    int vel_prof_num;
	const double frho;
    double fmu, fdt_pad;
    //VELOCITY PROFILE:
    //NO INITIAL VELOCITIES=0
    //CONSTANT X 1 vx
    //CONSTANT Y 2 vy
    //CONSTANT Z 3 vz
    double * vel_profile;
	/** Params for solids, shared or global settings */
	int n_obj, nlayers, max_extrap_rs;
	double wt_n;
	double ex_visc_mult, ev_trans_mult, sdt_pad, dt_ex_pad;
    double weight_fac;
	double gravity;
	/** boundary condition types */
    int bct[3][6];
    double bcv[3][6];


    /** Arrays of inidividual object settings, when sim_type == "objects" */
    //TYPE OF OBJECT:
    //SPHERE=0
    //DEFORMED_SPHERE=1
    //ROD=2
    //CUBE=3
    //THREE_POINT_ROTOR=4
    int * object_list;
    double * srho;
    double * shear_mod;
	double * basic_specs; // user specified center and dimensions of the objects
    // merge into extra_specs
    // double sqrt_lambda;
	double * extra_specs; // user specified extra dimensions, such as rod length
	char * sim_type;


	/** Constructor that scans a config file. */
	sim_params(const char *config);

	/** Destructor */
	~sim_params(){
		if (chk_dirname!=NULL) delete [] chk_dirname;
		if (recover_dirname!=NULL) delete [] recover_dirname;
		if (sim_type!=NULL) delete [] sim_type;
		if (extra_specs!=NULL) delete [] extra_specs;
		if (basic_specs!=NULL) delete [] basic_specs;
		if (srho!=NULL) delete [] srho;
		if (shear_mod!=NULL) delete [] shear_mod;
		if (object_list!=NULL) delete [] object_list;
		if (vel_profile!=NULL) delete [] vel_profile;
		if (dirname!=NULL) delete [] dirname;
	}

	//private:
	/** Tests to see if two strings are equal.
	 * \param[in] p1 a pointer to the first string.
	 * \param[in] p2 a pointer to the second string.
	 * \return True if they are equal, false otherwise. */
	inline bool se(const char *p1,const char *p2) {
		return strcmp(p1,p2)==0;
	}
	/** Finds the next token in a string and interprets it as a
	 * double precision floating point number. If none is availble,
	 * it gives an error message.
	 * \param[in] ln the current line number. */
	inline double next_double(int ln) {
		return atof(next_token(ln));
	}
	/** Finds the next token in a string, interprets it as a double
	 * precision floating point number, and checks that there are
	 * no subsequent values.
	 * \param[in] ln the current line number. */
	inline double final_double(int ln) {
		double temp=next_double(ln);
		check_no_more(ln);
		return temp;
	}
	/** Finds the next token in a string and interprets it as an
	 * integer. If none is availble, it gives an error message.
	 * \param[in] ln the current line number. */
	inline int next_int(int ln) {
		return atoi(next_token(ln));
	}
	/** Finds the next token in a string, interprets it as an
	 * integer, and checks that there are no subsequent values.
	 * \param[in] ln the current line number. */
	inline int final_int(int ln) {
		int temp=next_int(ln);
		check_no_more(ln);
		return temp;
	}
	/** Finds the next token in a string and interprets it as an
	 * unsigned int. If none is availble, it gives an error message.
	 * \param[in] ln the current line number. */
	inline unsigned int next_uint(int ln) {
		long lnum = atol(next_token(ln));
		return (unsigned int) (lnum);
	}
	/** Finds the next token in a string, interprets it as an
	 * integer, and checks that there are no subsequent values.
	 * \param[in] ln the current line number. */
	inline unsigned int final_uint(int ln) {
		unsigned int temp=next_uint(ln);
		check_no_more(ln);
		return temp;
	}
	char* next_token(int ln);
	void invalid_error(const char *cs,int ln);
	void check_no_more(int ln);
	void print_params();
	void write_params(const char* chk_dirname);
	void set_current_time(double t){cur_time=t;}
	double get_current_time(){return cur_time;}

	/**
	 * Method by which other classes can set dt (i.e. for
	 * synchronization across multiple info accessors)
	 */
	void set_dt(double dt_) {dt=dt_;}


    void set_current_chk_step(const int cn){chk_step=cn;}
	int get_current_step(){return chk_step;}

	void set_current_chk_num(const int cn){chk_num=cn;}
    inline double min_dh() const{
        double dh = (dx<dy)?dx:dy;
        dh = (dh<dz)?dh:dz;
        return dh;
    }

	private:
	int rank;
};

#endif
