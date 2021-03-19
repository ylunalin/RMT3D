#ifndef WRITE_PARAMS_HH
#define WRITE_PARAMS_HH
#include <cstdlib>
#include <cstdio>

// Output type key
// 0: Fluid velocity u_x
// 1: Fluid velocity u_y
// 2: Fluid velocity u_z
// 3: Vorticity w_x
// 4: Vorticity w_y
// 5: Vorticity w_z
// 6: Velocity magnitude
// 7: Pressure
// 8: Level set phi field
// 9: Reference map xi_x
// 10: Reference map xi_y
// 11: Reference map xi_z
// 12: Obj ID
// 13: Obj layers
// 14: Level set phi field evaluted on xpred
// 15: Predictive reference map xi_x
// 16: Predictive reference map xi_y
// 17: Predictive reference map xi_z
// 18: Divergence of u
// 19: Divergence of u edge
// 20: Density -> only output nonzero values when variable density
// 21: ||J-1|| -> only output sensible values when debug macro is on

/** An object that holds the paramters for writing out slices */
class write_params{
	public:
	// the dimension to slice
	int dim;
	// the coordinate along dimension to slice, in unit of grid points
	int point;
	// the output type
	int o_type;
	// the object id to consider (only relevant for solid fields). If this is
	// -1, then combined output is done based on the minimum phi value
	int obj_id;
	// the output format (0: matrix, 1: text, 2: gnuplot)
	int format;
	static const int numf=22;
	// the filename to write out
	const char* filename;
	static const char* default_names[numf];

	write_params(int dim_, int point_, int o_type_, int obj_id_, int format_, const char* filename_)
			: dim(dim_), point(point_), o_type(o_type_), obj_id(obj_id_),
			format(format_), filename(filename_) {}
	write_params(int dim_, int point_, int o_type_, int obj_id_, int format_)
			: dim(dim_), point(point_), o_type(o_type_), obj_id(obj_id_),
			format(format_), filename(default_names[o_type]) {}
	~write_params() {}

	// WRITE FUNCTIONS
	inline void change_otype(int t){
		o_type = t;
		filename = default_names[t];
	}
	inline void change_point(int p){
		point = p;
	}
	// READ ONLY FUNCTIONS
	inline int get_otype(){ return o_type; }
	inline int get_point(){ return point; }
	inline bool corner_field() {
		bool tmp = (o_type>=3&&o_type<6 ) || o_type == 7 || o_type == 18;
		return tmp;
	}
};

#endif
