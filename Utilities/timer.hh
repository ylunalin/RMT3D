#ifndef TIMER
#define TIMER

#include <mpi.h>
#include <vector>
#include <cstdarg>

enum unit{
	second = 's',
	minute = 'm',
	hour = 'h',
	day = 'd'
};

#define TIME_FUNC(TIMER_NAME, TIMER_NUM, FUNC_NAME){\
    TIMER_NAME.tic(TIMER_NUM); \
    FUNC_NAME(); \
    TIMER_NAME.toc(TIMER_NUM);\
}

#define TIME_FUNC_ARGS(TIMER_NAME, TIMER_NUM, FUNC_NAME, ...) {\
    TIMER_NAME.tic(TIMER_NUM); \
    FUNC_NAME(__VA_ARGS__); \
    TIMER_NAME.toc(TIMER_NUM);\
}

class timer {

	public:
	int ntypes;
    int name_size;
    double total_time;
	double *times;
	double *t0;
    const char* watch_name;
    std::vector<const char*> names;

	/** constructs timer of n different types */
	timer(const char* wname, int n, ...):
	ntypes(n), total_time(0.),
    times(new double[ntypes]),
	t0(new double[ntypes]),
    watch_name(wname)
	{
        va_list temp_names;
        va_start(temp_names, n);
		// The last entries is the total of all processes
		for (int i = 0; i < ntypes; i++) {
            times[i] = 0;
            t0[i] = 0;
            names.push_back( va_arg(temp_names, const char*));
        }
        va_end(temp_names);
        name_size = static_cast<int> (names.size());
	}
	/** destructor */
	~timer() {
		delete[] times;
		delete[] t0;
	}
	/** starts timer for given type */
	void tic(int type) {if (type<ntypes) t0[type] = MPI_Wtime();};
	/** returns time since type's last tic(), adds to running total */
	void toc(int type) {
		double t =0;
        if (type<ntypes){
            t = MPI_Wtime() - t0[type];
            times[type] += t;
            total_time += t;
        }
	}
	/** returns current value of running total for given type */
	double elapsed(int type) {
        if (type<ntypes) return times[type];
        else return 0.;
    }
	/** sets total elapsed time to 0 */
	void reset(int type) {if (type<ntypes) times[type] = 0;};
    /** set all timer to zero */
    void reset_all() {
        for(int i=0;i<ntypes;i++){ times[i] = 0.; }
    }
    const char* get_name(int type){
        if (type<ntypes) return names[type];
        else return NULL;
    }
	/** converts time into appropriate unit */
	unit convert(double &t_out, double t){
        unit u = second;
        if(t<100) u = second;
        else if(t>=100 && t<6000){
            t /= 60.0;
            u = minute;
        }
        else if(t>=6000 && t<360000){
            t/=3600.;
            u = hour;
        }
        else{
            t/=86400.;
            u = day;
        }

        t_out = t>0?t:-1;
		return u;
	}
    void report(){
        printf("# Timer  %s report:\n", watch_name);
        double tout;
        unit u = convert(tout, total_time);
        printf("#     Total tracked compute time %6.4g %c\n", tout, u);
        for(int i=0;i<ntypes;i++){
            u = convert(tout, times[i]);
            printf("#      %s time %6.4g %c\n", names[i], tout, u);
        }
    }
};

#endif
