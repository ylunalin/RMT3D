/** Implementation of sim_params class.
 *
 * Author    : Y Luna Lin
 * Email     : y.lin2783@gmail.com
 * Date      : Nov 11, 2018
 */

#include "sim_params.hh"

sim_params::sim_params(const char *fn):
    // Default values of necessary parameters
    suppress(false), display(false), impl(false), godunov(true),
    x_prd(true), y_prd(true), z_prd(true), has_exact(false),
    walls{false, false, false, false, false, false},
    chk_step(0), chk_num(0), frames(100), stop_short(-1),
    chkpt_freq(1), omp_num_thr(1),
    debug_flag(0), debug_obj_id(-1),
    obj_body(-1), data_fmt(0), out_flag(1|2|4|128),
    output_dim(2), output_ind(-1), dump_code(1|2),
    ntracers(1000), num_iters(5),
    dt(1.), T(1.), cur_time(0),
    ax(0), ay(0), az(0), bx(1),
    by(1), bz(1), lx(bx-ax), ly(by-ay), lz(bz-az),
    wall_pos{ax, bx, ay, by, az, bz}, wall_dist (0),
    // fluid velocity profile and material constants
    vel_prof_num(0),
    frho(1), fmu (0.002), fdt_pad(0.5),
	vel_profile(NULL),
    // solid params, assume no solid, only add as user provides
    n_obj(0), nlayers(8),
    max_extrap_rs(5),
    wt_n(2.5),
    ex_visc_mult(1.), ev_trans_mult(1.), sdt_pad(0.25),
    dt_ex_pad(0.8), weight_fac(0.25),
    gravity(0),
    bct{{0,0,0,0,0,0},{0,0,0,0,0,0},{0,0,0,0,0,0}},bcv{{0,0,0,0,0,0},{0,0,0,0,0,0},{0,0,0,0,0,0}},
    object_list(NULL),
    srho(NULL), shear_mod(NULL),
    basic_specs(NULL), extra_specs(NULL)
    {
    for(int i=0;i<3;i++) sys_size[i] = 96;

    MPI_Comm_rank(world, &rank);
    sim_type = new char[256];
    recover_dirname = new char[256];
    sprintf(sim_type, "%s", "ftest");
    int i, ln=1, len = strlen(fn);
    char buf[buf_size], *bp;
    if(len<4 || fn[len-4]!='.' || fn[len-3]!='c' || fn[len-2]!='f' || fn[len-1]!='g'){
            fatal_error("sim_params::sim_params(const char*): Filename must end in '.cfg'\n", 1);
    }

    // Use the root of config filename as output directory name
    dirname = new char[len+1];
    for(i=0;i<len-3;i++) dirname[i]=fn[i];
    dirname[len-3]='o';
    dirname[len-2]='u';
    dirname[len-1]='t';
    dirname[len]='\0';
    // Initially set recovery name to the current directory name
    for(i=0;i<len+1;i++) recover_dirname[i] = dirname[i];

    // We also write down the recovery directory name
    FILE *f=safe_fopen(fn,"r");
    if(f==NULL) {
        exit(1);
    }
    while(!feof(f)){
        // Reads in config file one line at a time
        if(fgets(buf, buf_size, f)==NULL) break;

        // Replace comment character with a null character, skip commented lines
        bp=buf;
        while((*bp)!=0){
            if(*bp=='#') {*bp=0; break;}
            bp++;
        }

        // Search for a keyword, skip if none found
        bp = strtok(buf, " \t\n");
        if(bp==NULL){
            ln++;
            continue;
        }
        // Compare keyword with known keywords
        if(se(bp, "suppress")){
            suppress = final_int(ln);
        } else if(se(bp, "display")){
            display = final_int(ln);
        } else if(se(bp, "implicit")){
            impl = final_int(ln);
        } else if(se(bp, "godunov")){
            godunov = final_int(ln);
        } else if(se(bp, "x_prd")){
            x_prd = final_int(ln);
        } else if(se(bp, "y_prd")){
            y_prd = final_int(ln);
        } else if(se(bp, "z_prd")){
            z_prd = final_int(ln);
        } else if(se(bp, "has_exact")){
            has_exact = final_int(ln);
        } else if(se(bp, "walls")){
            for(int i=0;i<5;i++) walls[i] = next_int(ln);
            walls[5] = final_int(ln);
        } else if(se(bp, "wall_pos")){
            for(int i=0;i<5;i++) wall_pos[i] = next_double(ln);
            wall_pos[5] = final_double(ln);
        } else if(se(bp, "wall_dist")){
            wall_dist = next_double(ln);
        } else if(se(bp, "recover_dirname")){
            sprintf(recover_dirname, "%s", next_token(ln));
        } else if(se(bp, "chk_num")){
            chk_num = final_int(ln);
        } else if(se(bp, "chk_step")){
            chk_step = final_int(ln);
        } else if(se(bp, "sys_size")){
            int tmp = next_int(ln);
            if(tmp>0) sys_size[0] = tmp;
            tmp = next_int(ln);
            if(tmp>0) sys_size[1] = tmp;
            tmp = final_int(ln);
            if(tmp>0) sys_size[2] = tmp;
        } else if(se(bp, "frames")){
            frames = final_int(ln);
        } else if(se(bp, "stop_short")){
            stop_short = final_int(ln);
        } else if(se(bp, "chkpt_freq")){
            chkpt_freq = final_int(ln);
        } else if(se(bp, "omp_num_thr")){
            int tmp = final_int(ln);
            if(tmp>0) omp_num_thr = tmp;
        } else if(se(bp, "debug_flag")){
            debug_flag = final_int(ln);
        } else if(se(bp, "debug_obj_id")){
            debug_obj_id = final_int(ln);
        } else if(se(bp, "out_flag")){
            out_flag = final_int(ln);
        } else if(se(bp, "obj_body")){
            obj_body = final_int(ln);
        } else if(se(bp, "data_fmt")){
            data_fmt = final_int(ln);
        } else if(se(bp, "output_dim")){
            output_dim = final_int(ln);
            if(output_dim < 0 || output_dim > 2 ) output_dim = -1;
        } else if(se(bp, "output_ind")){
            output_ind = final_int(ln);
            if(output_ind < 0 || output_ind >=sys_size[output_dim]) output_ind = -1;
        } else if(se(bp, "dump_code")){
            dump_code =  final_uint(ln);
        } else if(se(bp, "ntracers")){
            ntracers = final_int(ln);
        } else if(se(bp, "num_iters")){
            num_iters = final_int(ln);
        } else if(se(bp, "dt")){
            dt = final_double(ln);
        } else if(se(bp, "T")){
            double tmp = final_double(ln);
            if(tmp>0) T = tmp;
        } else if(se(bp, "current_time")){
            double tmp = final_double(ln);
            if(tmp>0) cur_time = tmp;
        } else if(se(bp, "xyz_bounds")){
            ax = next_double(ln);
            bx = next_double(ln);
            ay = next_double(ln);
            by = next_double(ln);
            az = next_double(ln);
            bz = final_double(ln);
            lx = bx-ax; ly = by-ay; lz = bz-az;
        } else if(se(bp, "fmu")){
            fmu = final_double(ln);
        } else if(se(bp, "fdt_pad")){
            fdt_pad = final_double(ln);
        } else if(se(bp, "vel_prof_num")){
            vel_prof_num = final_int(ln);
            if(vel_prof_num>0){
                vel_profile = new double[2*vel_prof_num];
            }
        } else if(se(bp, "n_obj")){
            n_obj = final_int(ln);
            if(n_obj>0 && object_list==NULL){
                object_list = new int[n_obj];
                basic_specs = new double[n_obj*n_basic_specs];
				extra_specs = new double[n_obj*n_extra_specs];
                srho = new double[n_obj];
                shear_mod = new double[n_obj];
                for(i=0;i<n_obj;i++) {
                    object_list[i] = 0;
                    srho[i] = 1.; shear_mod[i] = 1.;
                    extra_specs[n_extra_specs*i] = 1.;
                    for(int j=1;j<n_extra_specs;j++) extra_specs[n_extra_specs*i+j] = 0.;
                }
            }
        } else if(se(bp, "nlayers")){
            nlayers = final_int(ln);
        } else if(se(bp, "max_extrap_rs")){
            max_extrap_rs = final_int(ln);
            if(max_extrap_rs < 2) {
                max_extrap_rs = 2;
            }
        } else if(se(bp, "wt_n")){
            wt_n = final_double(ln);
        } else if(se(bp, "ex_visc_mult")){
            ex_visc_mult = final_double(ln);
        } else if(se(bp, "ev_trans_mult")){
            ev_trans_mult = final_double(ln);
        } else if(se(bp, "sdt_pad")){
            sdt_pad = final_double(ln);
        } else if(se(bp, "dt_ex_pad")){
            dt_ex_pad = final_double(ln);
        } else if(se(bp, "weight_fac")){
            weight_fac = final_double(ln);
        } else if(se(bp, "gravity")){
            gravity = final_double(ln);
        } else if(se(bp, "bc_xl")){
            int tmp_type[3]={0,0,0};
            double tmp_value[3]={0,0,0};
            for(i=0;i<3;i++) tmp_type[i] = next_int(ln);
            for(i=0;i<2;i++) tmp_value[i] = next_double(ln);
            tmp_value[2] = final_double(ln);
            for(i=0;i<3;i++) {
                bct[i][0] = tmp_type[i];
                bcv[i][0] = tmp_value[i];
            }
        } else if(se(bp, "bc_xh")){
            int tmp_type[3]={0,0,0};
            double tmp_value[3]={0,0,0};
            for(i=0;i<3;i++) tmp_type[i] = next_int(ln);
            for(i=0;i<2;i++) tmp_value[i] = next_double(ln);
            tmp_value[2] = final_double(ln);
            for(i=0;i<3;i++) {
                bct[i][1] = tmp_type[i];
                bcv[i][1] = tmp_value[i];
            }
        } else if(se(bp, "bc_yl")){
            int tmp_type[3]={0,0,0};
            double tmp_value[3]={0,0,0};
            for(i=0;i<3;i++) tmp_type[i] = next_int(ln);
            for(i=0;i<2;i++) tmp_value[i] = next_double(ln);
            tmp_value[2] = final_double(ln);
            for(i=0;i<3;i++) {
                bct[i][2] = tmp_type[i];
                bcv[i][2] = tmp_value[i];
            }
        } else if(se(bp, "bc_yh")){
            int tmp_type[3]={0,0,0};
            double tmp_value[3]={0,0,0};
            for(i=0;i<3;i++) tmp_type[i] = next_int(ln);
            for(i=0;i<2;i++) tmp_value[i] = next_double(ln);
            tmp_value[2] = final_double(ln);
            for(i=0;i<3;i++) {
                bct[i][3] = tmp_type[i];
                bcv[i][3] = tmp_value[i];
            }
        } else if(se(bp, "bc_zl")){
            int tmp_type[3]={0,0,0};
            double tmp_value[3]={0,0,0};
            for(i=0;i<3;i++) tmp_type[i] = next_int(ln);
            for(i=0;i<2;i++) tmp_value[i] = next_double(ln);
            tmp_value[2] = final_double(ln);
            for(i=0;i<3;i++) {
                bct[i][4] = tmp_type[i];
                bcv[i][4] = tmp_value[i];
            }
        } else if(se(bp, "bc_zh")){
            int tmp_type[3]={0,0,0};
            double tmp_value[3]={0,0,0};
            for(i=0;i<3;i++) tmp_type[i] = next_int(ln);
            for(i=0;i<2;i++) tmp_value[i] = next_double(ln);
            tmp_value[2] = final_double(ln);
            for(i=0;i<3;i++) {
                bct[i][5] = tmp_type[i];
                bcv[i][5] = tmp_value[i];
            }
        } else if(se(bp, "basic_specs")){
            if(basic_specs==NULL) printf("sim_params:: warning: basic_specs is not allocated but you attempted to access it!\n");
            int index = next_int(ln);
            if(index < n_obj){
				for(i=0;i<n_basic_specs-1;i++) basic_specs[index*n_basic_specs+i] = next_double(ln);
                basic_specs[(index+1)*n_basic_specs-1] = final_double(ln);
            }
        } else if(se(bp, "vel_profile")){
            if(vel_profile==NULL) printf("sim_params:: warning: vel_profile is not allocated but you attempted to access it!\n");
            int index = next_int(ln);
            if(index < vel_prof_num){
                int vtype = next_int(ln);
                vel_profile[index*2] = static_cast<double> (vtype);
                if(vtype !=0) vel_profile[index*2+1] = final_double(ln);
            }
        } else if(se(bp, "srho")){
            if(srho==NULL) printf("sim_params:: warning: srho is not allocated but you attempted to access it!\n");
            int index = next_int(ln);
            if(index < n_obj){
                srho[index] = final_double(ln);
            }
        } else if(se(bp, "shear_mod")){
            if(shear_mod==NULL) printf("sim_params:: warning: shear_mod is not allocated but you attempted to access it!\n");
            int index = next_int(ln);
            if(index < n_obj){
                shear_mod[index] = final_double(ln);
            }
        } else if(se(bp, "object_list")){
            if(object_list==NULL) printf("sim_params:: warning: object_list is not allocated but you attempted to access it!\n");
            int index = next_int(ln);
            if(index < n_obj){
                object_list[index] = final_int(ln);
            }
        } else if(se(bp, "extra_specs")){
			if(extra_specs==NULL) printf("sim_params:: warning: extra_specs is not allocated but you attempted to access it!\n");
			int index = next_int(ln);
            if(index < n_obj){
				for(i=0;i<n_extra_specs-1;i++) extra_specs[index*n_extra_specs+i] = next_double(ln);
                extra_specs[(index+1)*n_extra_specs-1] = final_double(ln);
            }
        } else if(se(bp, "sim_type")){
                sprintf(sim_type, "%s", next_token(ln));
        } else {
                printf("sim_params::sim_params(const char*): Unrecognized keyword '%s' at line %d of file %s\n", bp, ln, fn);
    print_params();
                exit(1);
        }
        ln++;

    }

    // some paramters depend on others, make sure to update them
    if(output_ind==-1) output_ind = sys_size[output_dim]/2;
    lx=bx-ax; ly=by-ay; lz=bz-az;

    int bc_sum[3]={0,0,0};
    for(i=0;i<3;i++){
        bc_sum[i] = bct[0][2*i] + bct[1][2*i] + bct[2][2*i];
    }
    for(i=0;i<3;i++){
        // the least non-zero code for bc is 3, for Dirichlet BC all around
        // if this condition is not met, we zero everything out and consider things periodic
        if(bc_sum[i]>=0 && bc_sum[i]<3) {
            bct[0][2*i] = bct[1][2*i] = bct[2][2*i] = 0;
            bct[0][2*i+1] = bct[1][2*i+1] = bct[2][2*i+1] = 0;
            switch (i) {
                case 0: x_prd=1; break;
                case 1: y_prd=1; break;
                case 2: z_prd=1; break;
            }
        } else{
        // otherwise, we make sure to turn periodic flags off accordingly
        // user might have set it already, but we double make sure
            switch (i) {
                case 0: x_prd=0; break;
                case 1: y_prd=0; break;
                case 2: z_prd=0; break;
            }
        }
    }

    if(!x_prd){
        walls[0] = walls[1] = true;
    }
    if(!y_prd){
        walls[2] = walls[3] = true;
    }
    if(!z_prd){
        walls[4] = walls[5] = true;
    }

    // system size is the number of cells
    // dx, dy, dz then is unaffected by prd boundary conditions
    dx = lx/sys_size[0];
    dy = ly/sys_size[1];
    dz = lz/sys_size[2];
    wall_dist = std::max(min_dh(), wall_dist);

    // after recovery directory is read in, we check if checkpoint files exist
    chk_dirname = NULL;
    if(chk_num>0 && recover_dirname!=NULL) {
        chk_dirname = new char[256];
        sprintf(chk_dirname, "%s/chk.%05d", recover_dirname, chk_num);
        int exist = 1;
        if(rank ==0) {
            char tmp_filename[256];
            // let's test a file to see if the directory is accessible
            sprintf(tmp_filename, "%s/recover.cfg", chk_dirname);
            FILE *fh = fopen(tmp_filename, "rb");
            if(fh==NULL) {
                printf("sim_params::Cannot find checkpt directory %s. Reset to none.\n", chk_dirname);
                exist=0;
            }
            if(fh!=NULL) fclose(fh);
        }
        MPI_Bcast(&exist, 1, MPI_INT, 0, MPI_COMM_WORLD);
        if(!exist){
            chk_num = 0;
            chk_step = 0;
            delete [] chk_dirname;
            chk_dirname=NULL;
        }
    }

    if(strcmp(sim_type, "ftest") == 0) n_obj=0;
    if(strcmp(sim_type, "stest") == 0) n_obj=0;
    if(strcmp(sim_type, "objects") == 0 && rank == 0 && n_obj==0) printf("sim_params:: warning, configuring for a simulation of objects, but there is NO settings for objects!\n");
#if !defined(VAR_DEN)
    bool var_den=false;
    for(i=0;i<n_obj;i++) {
        if(srho[i]!=1.0) var_den=true;
    }
    if(var_den && rank==0) printf("sim_params:: warning, solid density is different from fluid density, but multigrid solves do not use variable density formulation.\n");
#endif
}

/** Finds the next token in a string, and if none is availble, gives an error
 * message.
 * \param[in] ln the current line number. */
char* sim_params::next_token(int ln) {
    char *temp=strtok(NULL," \t\n");
    if(temp==NULL) {
        fprintf(stderr,"Not enough arguments at input line %d\n",ln);
        exit(1);
    }
    return temp;
}

/** Prints a message about an invalid quantity and exits.
 * \param[in] cs the invalid quantity.
 * \param[in] ln the current line number. */
void sim_params::invalid_error(const char* cs,int ln) {
    fprintf(stderr,"Invalid %s at line %d",cs,ln);
    exit(1);
}

/** Checks that there are no subsequent values.
 * \param[in] ln the current line number. */
void sim_params::check_no_more(int ln) {
    if(strtok(NULL," \t\n")!=NULL) {
        fprintf(stderr,"Too many arguments at input line %d\n",ln);
        exit(1);
    }
}

void sim_params::print_params(){
    printf("            Keyword                   Input(default)           Type\n");
    printf("            suppress                  0/1(0)                   bool\n"
           "            display                   0/1(0)                   bool\n"
           "            implicit                  0/1(0)                   bool\n"
           "            godunov                   0/1(0)                   bool\n"
           "            x_prd                     0/1(1)                   bool\n"
           "            y_prd                     0/1(1)                   bool\n"
           "            z_prd                     0/1(1)                   bool\n"
           "            has_exact                 0/1(0)                   bool\n"
           "            chk_step                  (0)                      int\n"
           "            chk_num                   (0)                      int\n"
           "            sys_size                  ({96,96,96})             int\n"
           "            frames                    (96)                     int\n"
           "            stop_short                (96)                     int\n"
           "            chkpt_freq                (1)                      int\n"
           "            omp_num_thr               (1)                      int\n"
           "            debug_flag                (0)                      int\n"
           "            debug_obj_id              (-1)                     int\n"
           "            obj_body                  (-1)                     int\n"
           "            data_fmt                  (0)                      int\n"
           "            out_flag                  (135)                    int\n"
           "            output_dim                0/1/2(2)                 int\n"
           "            output_ind                (48)                     int\n"
           "            dump_code                 (3)                      unsigned int\n"
           "            ntracers                  (1000)                   int\n"
           "            dt                        (1.)                     double\n"
           "            T                         (1.)                     double\n"
           "            current_time              (0.)                     double\n"
           "            xyz_bounds                ax bx ay by az bz        doubles\n"
           "                                      (0 1 0 1 0 1)\n"
           "            fmu                       (0.002)                  double\n"
           "            fdt_pad                   (0.5)                    double\n"
           "            n_obj                     (0)                      int\n"
           "            nlayers                   (8)                      int\n"
           "            wt_n                      (2.5)                    double\n"
           "            srho                      (ind, sr)                double\n"
           "            shear_mod                 (ind, sm)                double\n"
           "            ex_visc_mult              (1.)                     double\n"
           "            ev_trans_mult             (1.)                     double\n"
           "            sdt_pad                   (0.25)                   double\n"
           "            dt_ex_pad                 (0.8)                    double\n"
           "            gravity                   (10.)                    double\n"
           "            basic_specs               (ind, x0,y0,z0,r0)       double\n"
           "            extra_specs                                        double\n"
           "            sim_type                  \"ftest\"                string\n");
}

void sim_params::write_params(const char * chk_dirname){
    if(rank==0){
        char buf[buf_size];
        sprintf(buf, "%s/recover.cfg", chk_dirname);
        FILE *fh = p_safe_fopen(buf, "w");
        fprintf(fh, "#Keyword                  Input(default)           Type\n");
        fprintf(fh, "#System parameters\n");
        fprintf(fh, "suppress                  %d#(0)                   bool\n"
                    "display                   %d#(0)                   bool\n"
                    "implicit                  %d#(0)                   bool\n"
                    "godunov                   %d#(0)                   bool\n"
                    "x_prd                     %d#(0)                   bool\n"
                    "y_prd                     %d#(0)                   bool\n"
                    "z_prd                     %d#(0)                   bool\n"
                    "has_exact                 %d#(0)                   bool\n"
                    "recover_dirname           %s#                      string\n"
                    "chk_num                   %d#(0)                   int\n"
                    "chk_step                  %d#(0)                   int\n"
                    "sys_size                  %d %d %d#({96,96,96})    int\n"
                    "frames                    %d#(96)                  int\n"
                    "stop_short                %d#(96)                  int\n"
                    "chkpt_freq                %d#(1)                   int\n"
                    "omp_num_thr               %d#(1)                   int\n"
                    "debug_flag                %d#(0)                   int\n"
                    "debug_obj_id              %d#(-1)                   int\n"
                    "obj_body                  %d#(-1)                  int\n"
                    "data_fmt                  %d#(0)                   int\n"
                    "out_flag                  %d#(135)                 int\n"
                    "output_dim                %d#0/1/2(2)              int\n"
                    "output_ind                %d#(48)                  int\n"
                    "dump_code                 %u#(3)                   unsigned int\n"
                    "ntracers                  %d#(1000)                int\n"
                    "num_iters                 %d#(5)                   int\n"
                    "dt                        %g#(1.)                  double\n"
                    "T                         %g#(1.)                  double\n"
                    "current_time              %g#(0.)                  double\n"
                    "xyz_bounds                %g %g %g %g %g %g\n"
                    "                          #(0 1 0 1 0 1)         doubles\n"
                    "# SIMULATION TYPE CONFIGURATIONS\n"
                    "sim_type                  %s#(\"ftest\")           string\n"
                    "# FLUID PROPERTIES\n"
                    "fmu                       %g#(0.002)               double\n"
                    "fdt_pad                   %g#(0.5)                 double\n"
                    "# SOLID PROPERTIES\n"
                    "n_obj                     %d#(0)                   int\n"
                    "nlayers                   %d#(8)                   int\n"
                    "wt_n                      %g#(2.5)                 double\n"
                    "ex_visc_mult              %g#(1.)                  double\n"
                    "ev_trans_mult             %g#(1.)                  double\n"
                    "sdt_pad                   %g#(0.25)                double\n"
                    "dt_ex_pad                 %g#(0.8)                 double\n"
                    "gravity                   %g#(10.)                 double\n"
                    ,suppress, display, impl, godunov, x_prd, y_prd, z_prd, has_exact,
                    dirname, chk_num, chk_step, sys_size[0], sys_size[1], sys_size[2],
                    frames, stop_short,
                    chkpt_freq, omp_num_thr, debug_flag, debug_obj_id,
                    obj_body, data_fmt, out_flag, output_dim,
                    output_ind, dump_code, ntracers, num_iters,
                    dt, T, cur_time,
                    ax,bx,ay,by,az,bz,
                    sim_type, fmu, fdt_pad,
                    n_obj, nlayers, wt_n, ex_visc_mult,
                    ev_trans_mult, sdt_pad, dt_ex_pad, gravity);
                    fprintf(fh, "# SOLIDS CONFIGURATIONS\n");
                    if(object_list!=NULL){
                        for(int i=0;i<n_obj;i++){
                            fprintf(fh, "object_list               %d %d\n", i, object_list[i]);
                        }
                    }
                    if(basic_specs!=NULL){
                        for(int i=0;i<n_obj;i++){
                            fprintf(fh, "basic_specs               %d ", i);
							for(int j=0;j<n_basic_specs;j++) fprintf(fh, "%g ", basic_specs[n_basic_specs*i+j]);
							fprintf(fh, "\n");
                        }
                    }
                    if(extra_specs!=NULL){
                        for(int i=0;i<n_obj;i++){
                            fprintf(fh, "extra_specs               %d ", i);
							for(int j=0;j<n_extra_specs;j++) fprintf(fh, "%g ", extra_specs[n_extra_specs*i+j]);
							fprintf(fh, "\n");
                        }
                    }
                    if(srho!=NULL){
                        for(int i=0;i<n_obj;i++){
                            fprintf(fh, "srho                      %d %g\n", i, srho[i]);
                        }
                    }
                    if(shear_mod!=NULL){
                        for(int i=0;i<n_obj;i++){
                            fprintf(fh, "shear_mod                 %d %g\n", i, shear_mod[i]);
                        }
                    }

                    fprintf(fh, "vel_prof_num              %d#(1)                   int\n", vel_prof_num);
                    if(vel_profile!=NULL){
                        for(int i=0;i<vel_prof_num;i++){
                            fprintf(fh, "vel_profile               %d %d %g\n", i, static_cast<int>(vel_profile[i*2]), vel_profile[i*2+1]);
                        }
                    }
                    fprintf(fh, "# DOMAIN BOUNDARY CONDITIONS\n");
                    fprintf(fh, "bc_xl                     %d %d %d %g %g %g\n", bct[0][0], bct[1][0], bct[2][0], bcv[0][0], bcv[1][0], bcv[2][0]);
                    fprintf(fh, "bc_xh                     %d %d %d %g %g %g\n", bct[0][1], bct[1][1], bct[2][1], bcv[0][1], bcv[1][1], bcv[2][1]);

                    fprintf(fh, "bc_yl                     %d %d %d %g %g %g\n", bct[0][2], bct[1][2], bct[2][2], bcv[0][2], bcv[1][2], bcv[2][2]);
                    fprintf(fh, "bc_yh                     %d %d %d %g %g %g\n", bct[0][3], bct[1][3], bct[2][3], bcv[0][3], bcv[1][3], bcv[2][3]);

                    fprintf(fh, "bc_zl                     %d %d %d %g %g %g\n", bct[0][4], bct[1][4], bct[2][4], bcv[0][4], bcv[1][4], bcv[2][4]);
                    fprintf(fh, "bc_zh                     %d %d %d %g %g %g\n", bct[0][5], bct[1][5], bct[2][5], bcv[0][5], bcv[1][5], bcv[2][5]);
                    fprintf(fh, "# WALL BOUNDARIES\n");
                    fprintf(fh, "walls                     %d %d %d %d %d %d\n", walls[0], walls[1], walls[2], walls[3], walls[4], walls[5]);
                    fprintf(fh, "wall_pos                  %g %g %g %g %g %g\n", wall_pos[0], wall_pos[1], wall_pos[2], wall_pos[3], wall_pos[4], wall_pos[5]);


                    fprintf(fh, "wall_dist                 %g\n", wall_dist);
        fclose(fh);
    }
}
