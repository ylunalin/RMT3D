#ifndef STEPPER_HH
#define STEPPER_HH

#include <cstdio>
#include <cstring>
#include <vector>
#include <mpi.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <float.h>

#include "fields.hh"
#include "common.hh"
#include "custom_structs.hh"
#include "write_params.hh"
#include "sim_params.hh"

using namespace F3D_CONST;
/** class for managing progress through simulation, such that
 * some number of frames are outputted consisting of equal time
 * steps */
class stepper {
	public:
	/** length of loading bar */
	static const int bar_length = 40;
	/** whether to show bar */
	bool show_bar;
    /** whether there are tracers in the system */
    bool trace;
	/** To debug or not */
	int debug;
	/** Debugging solid id */
	int debug_obj_id;
    /** inital step */
    int init_step;
	/** current step within frame */
	int curr_step;
    /** initial frame */
    int init_frame;
	/** current frame */
	int curr_frame;
	/** current total step */
	int curr_tot_step;
	/** number of frames */
	int frames;
    /** stop-short frames */
    int stop_short;
	/** number of steps per frame (not counting remainders) */
	int steps;
	/** number of extra steps in the last frame */
	int last_steps;
	/** total number of steps (with remainders) */
	int total;
	/** integers for displaying time */
	int times[6];
	/** time to reach */
	int lines;

    /** frequency of checkpoint files. */
    int chkpt_freq;
	/** output operator code */
	// 1 == dump tracers
	// 2 == dump contours
	// 4 == write slices
	// 8 == check point files
	// 16 == diagnostics
	// 32 == macro info
	const unsigned int dump_code;
	double T;
	/** step time */
	double dt;
	/** frame time */
	double df;
	/** start time */
	double start;
	/** elapsed time */
	double elapse;
	/** estimated time to completion */
	double eta;
	/** format string for bar */
	char *bar_string;
	/** output directory */
	char *dir;
	/** A vector of write_params to dump each time. */
	std::vector<write_params> wparams;

	// constructors / destructors
	stepper(sim_params &spars);
	template<class p_class>
	stepper(sim_params &spars,p_class &p);
	~stepper() {delete[] bar_string;};

	// functions in cc file
	template<class p_class, typename F>
	void run(p_class &pr);

	template<class p_class>
	void run(p_class &pr);

	template<class p_class>
	int step(p_class &pr);

	template<class p_class>
	void write_files(p_class &pr,const int snum, const int fnum);

#if defined (DEBUG)
    template<class p_class>
    void run_diagnostics(p_class &pr);
#endif

    template<class p_class>
    void collect_macro_data(p_class &pr);

	void set_str();
	void add_output(int dim, int point, int otype, int obj_id, int format) {
		wparams.push_back(write_params(dim,point,otype,obj_id,format));
	}
	void add_output(int dim, int point, int otype, int obj_id, int format,const char* filename) {
		wparams.push_back(write_params(dim,point,otype,obj_id,format,filename));
	}
	template<class p_class>
	void dump_all(p_class &pr);

	// inline functions
	/** update the elapsed and estimated remaining times */
	inline void clock() {
		elapse = MPI_Wtime() - start;
		eta = elapse*((double) total/curr_tot_step - 1);
	}
	/** convert times to integer HH:MM:SS representations */
	inline void time_ints() {
		times[0] = elapse/3600; //hour
		times[1] = (int) (elapse/60 - 60*times[0]);
		times[2] = (int) (elapse - 3600*times[0] - 60*times[1]);
		times[3] = eta/3600;
		times[4] = (int) (eta/60 - 60*times[3]);
		times[5] = (int) (eta - 3600*times[3] - 60*times[4]);
	}
	/** push terminal cursors back to start of output */
	inline void reset_cursor() {
		for (int i = 0; i < lines; i++) printf("\033[F\33[2K");
		lines = 0;
	}
	/** print out the loading bar with the given number of blips */
	inline void print_blips(int blips) {
		printf("[");
		for (int i = 0; i < blips; i++) printf("*");
		for (int i = 0; i < bar_length-blips; i++) printf(" ");
		printf("]\n");
		lines++;
	}
	/** print out the problem classes's status string, if it exits */
	template<class p_class>
	inline void print_str(p_class &pr) {
		if (pr.out == NULL) return;
		printf("%s",pr.out);
		for (int i = 0; (unsigned)i < strlen(pr.out); i++) {
			if (pr.out[i] == '\n') lines++;
		}
	}
	/** print info on progress through the simulation */
	inline void print_info() {
		set_str();
		printf("\33[2K%s",bar_string);
		lines++;
	}
	/** print an initial message at simulation start */
	inline void print_empty_info() {
		if (frames > 2) {
			sprintf(bar_string,"initializing %d frames of %d time steps\n",
				frames,steps);
		} else {
			sprintf(bar_string,"initializing %d time steps\n",steps);
		}
		printf("%s",bar_string);
		lines++;
	}
	/** print the loading bar and progress info */
	inline void print_bar() {
		print_blips(curr_tot_step * bar_length / total);
		print_info();
		fflush(stdout);
	}
	/** print an loading bar and progress info */
	inline void print_empty_bar() {
		print_blips(0);
		print_empty_info();
		fflush(stdout);
	}
	/** print a complete progress update, with bar, progress info
	 * 	and problem status string */
	template <class p_class>
	inline void show_progress(p_class &pr,bool empty) {
		if (pr.rank != 0) return;
		if (!empty) reset_cursor();
		print_str(pr);
		if (empty) {
			print_empty_bar();
		} else {
			print_bar();
		}
	}
	private:
	template<class p_class> void link(p_class &pr);
	char fname_buf[256];
};

/** constructor for a loading bar stepper with the provided number
 * of frames, total simulation time T_, maximum time step dt, and
 * outputted data director dir_ */
template<class p_class>
stepper::stepper(sim_params &spars,p_class &pr) :
	show_bar(!spars.suppress), trace(spars.ntracers>0), debug(spars.debug_flag),
	debug_obj_id(spars.debug_obj_id),
	init_step(spars.chk_step),curr_step(0),
    init_frame(spars.chk_num),curr_frame(0),
	curr_tot_step(init_step),frames(spars.frames),
    stop_short(spars.stop_short),lines(0),
    chkpt_freq(spars.chkpt_freq),
	dump_code(spars.dump_code), T(spars.T),dt(spars.dt),df(T/(frames-1)),
	dir(spars.dirname) {

	steps = (int) (df/dt);
	double rem = df - steps*dt;
	if (rem>0) dt -= (dt-rem)/(++steps);

	total = (frames-1)*steps;

    // If user doesn't provide any step to stop short at,
    // set it to be the last step
    if(stop_short <=0) stop_short = total+1;
    // otherwise, add the number of steps before stopping
    // to the initial step
    else stop_short += init_step;

	if(frames-1>total) frames = total+1;

		/*

	// total is the number of steps
	total = (int) (T/dt);
	double rem = T - total*dt;
	if(rem>0) {
		total++;
		dt -= (dt-rem)/total;
	}

    // If user doesn't provide any step to stop short at,
    // set it to be the last step
    if(stop_short <=0) stop_short = total+1;
    // otherwise, add the number of steps before stopping
    // to the initial step
    else stop_short += init_step;

	if(frames-1>total) frames = total+1;

	steps = int ( total/(frames-1));
	last_steps = total%(frames-1);
    if(last_steps>steps){
        int extra_frames = (int) (last_steps/steps);
        last_steps =  last_steps%steps;
        frames += extra_frames;
        spars.frames= frames;
    }

	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);

	if (rank==0) {

		printf("-----------------------\n");

		printf("frames: %d, duration: %g\n",frames,T);

		printf("-----------------------\n");
	}

	MPI_Finalize();
	exit(0);
		*/

	bar_string = new char[512];

	wparams.clear();
	if(dump_code&4) {
        for (int i = 0; i < write_params::numf; i++) {
            if ((1<<i) & spars.out_flag) {
                add_output(spars.output_dim,spars.output_ind,i,spars.obj_body,spars.data_fmt);
            }
        }
    }

	link(pr);
}

/** synchronizes dt between stepper and p_class */
template<class p_class> void stepper::link(p_class &pr) {
	pr.set_dt(dt);
}

/** steps the problem class forward one time step */
template<class p_class>
int stepper::step(p_class &pr) {
	curr_tot_step++;
	int err = pr.step_forward(debug);
#if defined(DEBUG)
    run_diagnostics(pr);
#endif

    // debug flag was designed to indicate which time step you would like to dump everything out
    // negative values just print out a statement when each function in the step_forward routine
    // is called
    if(curr_tot_step == debug) dump_all(pr);

    // If step_forward has exited abnormally, we dump all and quit
    if(err != NORMAL_EXIT && err != STEADY_SOLN_EXIT) {
        dump_all(pr);
        clock();
        if(pr.rank==0) printf("Stepper: ! DID DUMP ALL IN THE END. Error code %d.\nTime spent on this failed simulation %10.5g seconds.\n", err, elapse);
        if(pr.rank==0) printf("Normal exit code %d steady soln exit code %d\n", NORMAL_EXIT, STEADY_SOLN_EXIT);
        p_fatal_error("Stepper: Something went wrong, sorry.", err);
    }


	clock();
	if (pr.rank == 0) time_ints();
	if (show_bar) show_progress(pr,false);

    return err;
}

template<class p_class>
void stepper::run(p_class &pr) {

	// update pr class dt
	pr.set_dt(dt);

	// Output the initial fields
	int lsteps;
	write_files(pr,init_step, init_frame);

#if defined(DEBUG)
    // While diagnostics are run per step
    // Macro data, like L2 norm of divergence
    // is collected every frame
    run_diagnostics(pr);
#endif
    collect_macro_data(pr);

    MPI_Barrier(pr.grid->cart);
	// start the races and show we've started
	start = MPI_Wtime();
	clock();
	if (show_bar) show_progress(pr,true);

	// go through all the frames we're doing
    int errcode = 0;
	for (curr_frame = init_frame+1; curr_frame <= frames-1 && curr_tot_step<stop_short; curr_frame++) {
		lsteps = (curr_frame == frames-1)?(steps+last_steps):steps;
		// do all the frame's steps
		for (curr_step = 0; curr_step < lsteps && curr_tot_step < stop_short; curr_step++) {
            errcode = step(pr);
            if(errcode == STEADY_SOLN_EXIT) break;
        }
		// Output the fields
		write_files(pr, curr_tot_step, curr_frame);
        collect_macro_data(pr);
        if(errcode == STEADY_SOLN_EXIT) break;
	}
	// drop a newline at the end
    clock();
	if (pr.rank == 0) printf("# Stepper: Done with everything, took a total %10.5g seconds.\n"
                            "Simulated to T=%g, a total of %d steps and %d frames\n"
                            "Error code %d\n"
                            , elapse, pr.time, curr_tot_step, curr_frame, errcode);
}

template<class p_class>
void stepper::collect_macro_data(p_class &pr){
    // dump_code 32 turns on the macro data output
    if(dump_code & 32) {
        // only one rank writes it out
        sprintf(fname_buf, "%s/macro.dat", dir);
        FILE *fh;
        // Check if the file already exists
        // if not, we ask rank 0 to write a new one
        if(pr.rank==0) {
            bool file_exist = false;
            fh = fopen(fname_buf, "r");
            if(fh !=NULL ){
                fclose(fh);
                file_exist = true;
            }
            // Write out diagnostic file header
            if(!file_exist || curr_tot_step == init_step) {
                fh = p_safe_fopen(fname_buf, "w");
                fprintf(fh, "%10s %6s %6s "
                            "%10s %6s %10s %6s %10s "
                            "%10s %10s %10s %10s "
                            "%10s %10s %10s "
							"%10s "
                            "%10s %10s %16s %16s %10s\n",
                            "#Time", "frame", "steps",
                            "div_u_max", "rank", "div_u_min", "rank", "divu_L2",
                            "tot_mmt_x","tot_mmt_y","tot_mmt_z","tot_energy",
                            "kinetic", "elastic", "potential",
							"power",
                            "max_dev_detF", "dev_detF_L2", "sing. obj. cent.", "volume", "change in u");
                fclose(fh);
            }
        }

        double_int gmx_ext, gmn_ext;
        double divu_l2 = pr.div_u(gmx_ext, gmn_ext);

        // Total momentum
        double mmt[3];
        pr.total_momentum(mmt);

        // Total energy
        double pot, kin, elas, pow;
        double energy = pr.total_energy(pot, kin, elas,pow);

        // Deviation of det(F) from 1
        double max_dev;
        double tot_dev_in_detF = pr.total_dev_in_detF(max_dev);

        int digs = DBL_DIG;
        if(pr.rank==0) {
            fh = p_safe_fopen(fname_buf, "a");
            fprintf(fh, "%.*e %06d %06d %.*e %06d %.*e %06d ", digs, pr.time, curr_frame, curr_tot_step, digs, gmx_ext.value, gmx_ext.rank, digs, gmn_ext.value, gmn_ext.rank);
            fprintf(fh, "%.*e ", digs, divu_l2);
            fprintf(fh, "%.*e %.*e %.*e %.*e ", digs, mmt[0], digs, mmt[1], digs, mmt[2], digs, energy);
            fprintf(fh, "%.*e %.*e %.*e %.*e ", digs, kin, digs, elas, digs, pot,digs,pow);
            fprintf(fh, "%.*e %.*e ", digs, max_dev, digs, tot_dev_in_detF);
            fclose(fh);
        }

        if(pr.mgmt->n_obj>0) {
            pr.compute_v_centroid();
            // only one rank writes it out
            sprintf(fname_buf, "%s/centroid.dat", dir);
            FILE *fh;
            // Check if the file already exists
            // if not, we ask rank 0 to write a new one
            if(pr.rank==0) {
                bool file_exist = false;
                fh = fopen(fname_buf, "r");
                if(fh !=NULL ){
                    fclose(fh);
                    file_exist = true;
                }
                // Write out diagnostic file header
                if(!file_exist || curr_tot_step == init_step) {
                    fh = p_safe_fopen(fname_buf, "w");
                    fprintf(fh, "%10s %6s %6s %6s %6s %6s %6s %6s %6s->\n",
                                "#Time", "frame", "steps", "x", "y", "z", "vx", "vy", "vz");
                    fclose(fh);
                }

                fh = p_safe_fopen(fname_buf, "a");
                fprintf(fh, "%.*e %d %d ", digs, pr.time, curr_frame, curr_tot_step);
                for (int i=0;i<pr.mgmt->n_obj;i++) {
                    fprintf(fh, "%.*e %.*e %.*e ", digs, pr.mgmt->avg_x[i], digs, pr.mgmt->avg_y[i], digs, pr.mgmt->avg_z[i]);
                    fprintf(fh, "%.*e %.*e %.*e ", digs, pr.mgmt->avg_velx[i], digs, pr.mgmt->avg_vely[i], digs, pr.mgmt->avg_velz[i]);
                }
                fprintf(fh, "\n");
                fclose(fh);
             }
         }
    } else {
        return;
    }
}


#if defined(DEBUG)
template<class p_class>
void stepper::run_diagnostics(p_class &pr){
    // dump_code 16 turns on the diagnostic output
    if(dump_code & 16) {
        // only one rank writes it out
        sprintf(fname_buf, "%s/diagnostics.dat", dir);
        FILE *fh;
        // Check if the file already exists
        // if not, we ask rank 0 to write a new one
        if(pr.rank==0) {
            bool file_exist = false;
            fh = fopen(fname_buf, "r");
            if(fh !=NULL ){
                fclose(fh);
                file_exist = true;
            }
            // Write out diagnostic file header
            if(!file_exist || curr_tot_step == init_step) {
                fh = p_safe_fopen(fname_buf, "w");
                fprintf(fh, "%10s %6s %6s %10s %6s %10s %6s "
                            "%10s %10s %10s %10s %10s %10s "
                            "%10s %10s %10s %10s "
                            "%10s %6s %10s %6s "
                            "%10s %6s %10s %6s "
                            "%10s %6s %10s %6s "
                            "%10s %10s %10s %10s %6s %6s %16s %16s\n",
                            "#Time", "frame", "steps", "div_u_max", "rank", "div_u_min", "rank",
                            "u_2norm","u_inf_norm","tot_mmt_x","tot_mmt_y","tot_mmt_z","tot_energy",
                            "kinetic", "elastic", "potential","power",
                            "mx_norm_str", "oid", "mx_shear_str", "oid",
                            "mx_coll_norm", "oid", "mx_coll_shear", "oid",
                            "mx_coll_traction", "oid", "mx_grad_phi_diff", "oid",
                            "avg_detF", "tot_dev_detF", "max_dev_detF", "max_gauss_curv", "oid", "sign", "sing. obj. cent.", "volume");
                fclose(fh);
            }
        }

        double_int gmx_ext, gmn_ext;
        pr.div_u(gmx_ext, gmn_ext);

        double vel_infty = pr.vel_mod(-1);
        double vel_l2 = pr.vel_mod(2);

        double mmt[3];
        pr.total_momentum(mmt);

        double_int gmx_solid_norm  = pr.max_solid_stress(0);
        double_int gmx_solid_shear = pr.max_solid_stress(1);
        double_int gmx_coll_norm   = pr.max_coll_stress(0);
        double_int gmx_coll_shear  = pr.max_coll_stress(1);
        double_int gmx_coll_tract  = pr.max_coll_traction();
        double_int gmx_grad_phi_diff = pr.max_grad_phi_diff();
        int sign = 0;
        double_int gmx_gaussian_curvature = pr.max_Gaussian_curvature(sign);

        double pot, kin, elas,pow;
        double energy = pr.total_energy(pot, kin, elas, pow);
        double max_dev;
        double tot_dev_in_detF = pr.total_dev_in_detF(max_dev);
        double avg_detF = pr.avg_detF();

        double centx, centy, centz, vol;
        pr.compute_centroid(centx, centy, centz, vol);
        if(pr.rank==0) {
            fh = p_safe_fopen(fname_buf, "a");
            fprintf(fh, "%+6.4g %06d %06d %+6.4e %06d %+6.4e %06d ", pr.time, curr_frame, curr_tot_step, gmx_ext.value, gmx_ext.rank, gmn_ext.value, gmn_ext.rank);
            fprintf(fh, "%+10.4e %+10.4e %+10.4e %+10.4e %+10.4e %+10.4e ", vel_l2, vel_infty, mmt[0], mmt[1], mmt[2], energy);
            fprintf(fh, "%+10.4e %+10.4e %+10.4e %+10.4e ", kin, elas, pot,pow);
            fprintf(fh, "%+10.4g %06d %+10.4e %06d ", gmx_solid_norm.value, gmx_solid_norm.rank, gmx_solid_shear.value, gmx_solid_shear.rank);
            fprintf(fh, "%+10.4g %06d %+10.4e %06d ", gmx_coll_norm.value, gmx_coll_norm.rank, gmx_coll_shear.value, gmx_coll_shear.rank);
            fprintf(fh, "%+10.4g %06d %+10.4e %06d %+10.4e %+10.4e %+10.4e ", gmx_coll_tract.value, gmx_coll_tract.rank, gmx_grad_phi_diff.value, gmx_grad_phi_diff.rank, avg_detF, tot_dev_in_detF, max_dev);
            fprintf(fh, "%+10.4g %06d %06d ", gmx_gaussian_curvature.value, gmx_gaussian_curvature.rank, sign);
            fprintf(fh, "%+10.4g %+10.4g %+10.4g %+10.8g\n",  centx, centy, centz, vol);
            fclose(fh);
        }
	/*

        MPI_Barrier(pr.grid->cart);
        if(gmx_solid_norm.value>0.5){
            sprintf(fname_buf, "%s/chk.%d", dir, -1);
            mkdir(fname_buf,S_IRWXU|S_IRWXG|S_IROTH|S_IXOTH);

            pr.write_chk_pt(curr_tot_step, -1, fname_buf);
            p_fatal_error("Stress value too hight, check the chkpoint files\n", 1);
        }

        if(elas<0){
            sprintf(fname_buf, "%s/determinant.%d", dir, pr.rank);
            pr.write_all_Js(fname_buf);

            MPI_Barrier(pr.grid->cart);
            p_fatal_error("Negative elastic energy, check determinant files\n", 1);
        }
        */

    } else {
        return;
    }
}
#endif

template<class p_class>
void stepper::write_files(p_class &pr, const int  snum, const int fnum) {

	if(trace && dump_code&1) {
		sprintf(fname_buf,"%s/tr.%05d",dir,fnum);
		pr.tr->print(fname_buf);
		sprintf(fname_buf,"%s/tr_st.%05d",dir,fnum);
		pr.tr->print(fname_buf,false);
        //printf("Rank %d. Write tracer %d\n", pr.rank, fnum);
	}

	if(dump_code&2) {
		sprintf(fname_buf,"%s/ctr.%05d",dir,fnum);
		pr.output_contours(fname_buf);
        //printf("Rank %d. Write contour %d\n", pr.rank, fnum);
	}

	// Write some slices if specified
	for(unsigned int i=0;i<wparams.size();i++) {
		sprintf(fname_buf,"%s/%s.%05d",dir,wparams[i].filename,fnum);
		pr.write_slice(wparams[i],fname_buf);
	}

	// Write checkpoint files if specified
	if( (dump_code & 8) && (fnum%chkpt_freq == 0)) {
        sprintf(fname_buf, "%s/chk.%05d", dir, fnum);
        mkdir(fname_buf,S_IRWXU|S_IRWXG|S_IROTH|S_IXOTH);

        sprintf(fname_buf, "%s/chk.%05d", dir, fnum);
		pr.write_chk_pt(snum, fnum, fname_buf);
	}
}

template<class p_class>
void stepper::dump_all(p_class &pr) {
	write_params dps(2,0,8,debug_obj_id,0);
	for(int i=0;i<=2;i++) {
		dps.change_otype(i);
		for(int o=0;o<pr.o;o++){
			dps.point = o;
			sprintf(fname_buf,"%s/debug_%s_%d.%d",dir,dps.filename,curr_tot_step, o);
			pr.write_slice(dps,fname_buf);
		}
	}
	for(int i=8;i<=17;i++) {
        if (i==12) continue;
		dps.change_otype(i);
		for(int o=0;o<pr.o;o++){
			dps.point = o;
			sprintf(fname_buf,"%s/debug_%s_%d.%d",dir,dps.filename,curr_tot_step, o);
			pr.write_slice(dps,fname_buf);
		}
	}
}

#endif
