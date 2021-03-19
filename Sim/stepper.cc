#include "stepper.hh"

/** constructor for a loading bar stepper with the provided number
 * of frames, total simulation time T_, maximum time step dt, and
 * outputted data director dir_ */
stepper::stepper(sim_params &spars) :
	show_bar(!spars.suppress), trace(spars.ntracers>0), debug(spars.debug_flag),
	debug_obj_id(spars.debug_obj_id),
	init_step(spars.chk_step),curr_step(0),
    init_frame(spars.chk_num),curr_frame(0),
	curr_tot_step(init_step),frames(spars.frames),
    stop_short(spars.stop_short),lines(0),
    chkpt_freq(spars.chkpt_freq),
	dump_code(spars.dump_code),T(spars.T),dt(spars.dt),df(T/(frames-1)),
	dir(spars.dirname) {


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
	steps = int ( total/(frames));
	last_steps = total%(frames);
    if(last_steps>steps){
        int extra_frames = (int) (last_steps/steps);
        last_steps =  last_steps%steps;
        frames += extra_frames;
        spars.frames= frames;
    }

	bar_string = new char[512];

	wparams.clear();
	if(dump_code&4) {
        for (int i = 0; i < write_params::numf; i++) {
            if ((1<<i) & spars.out_flag) {
                add_output(spars.output_dim,spars.output_ind,i,spars.obj_body,spars.data_fmt);
            }
        }
    }
}

/** sets the progress string for display under the bar */
void stepper::set_str() {
	if (frames > 2) {
		sprintf(bar_string,"Frame %d/%d, Step %d/%d,"
			" Elapsed: %d:%02d:%02d, ETA: %d:%02d:%02d\n",
			curr_frame+1,frames,curr_step+1,steps,
		times[0],times[1],times[2],times[3],times[4],times[5]);
	} else {
		sprintf(bar_string,"Step %d/%d, Elapsed: %d:%02d:%02d, ETA: %d:%02d:%02d\n",
			curr_step+1,steps,times[0],times[1],times[2],times[3],times[4],times[5]);
	}
}
