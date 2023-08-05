/* A header file that defines some commonly used
 * constants.
 */
#ifndef DEFS_HH
#define DEFS_HH
#include <limits>

namespace FIELD_CONST {
    const int U=0;
    const int V=1;
    const int W=2;
    const int ldigits = 10;
    const unsigned int omask = 1023;
    const unsigned int lmask = 63<<ldigits;
}

namespace BC{
    const int PERIODIC=0;
    const int DIRICHLET=1;
    const int NEUMANN=2;
}

namespace F3D_CONST {
    // F3D MPI COMM
    const int FFLAG_N=3;
    const int FFLAG=(1<<FFLAG_N)-1;
    const int SFLAG_N=4;
    const int SFLAG=((1<<SFLAG_N)-1)<<FFLAG_N;
    // EFLAG and EFLAG_N are obsolete
    const int EFLAG_N=4;
    const int EFLAG=((1<<EFLAG_N)-1)<<(FFLAG_N+SFLAG_N);
    // We have to communicate dvel too, for approximate projection
    const int DVELFLAG=1<<(FFLAG_N + SFLAG_N + EFLAG_N);

    // F3D IO
    const int fmsg_trans_dims=1<<6;
    const int fmsg_contour1=2<<6;
    const int fmsg_contour2=3<<6;
    const int fmsg_contour3=4<<6;
    const int fmsg_contour4=5<<6;

    // F3D Misc
    enum lower_faces {LEFT=0, FRONT=1, DOWN=2};
    const double small_number = 100.*std::numeric_limits<double>::epsilon();

    // ERROR CODES
    const int STEADY_SOLN_EXIT  = 100;
    const int NORMAL_EXIT       = 0;
    const int REFMAP_UPDATE_ERR = 1;
    const int PROJECTION_ERR    = 2;
}

namespace EXTRAP_MMAP_CONST{
    const int max_maps=16;
    const int init_map_mem=256;
    const int max_map_mem=256;
}

namespace OBJECT_DICT {
	const int n_basic_specs=4;
	const int n_extra_specs=10;

    const int object_type_total = 7;
    enum object_enum {
        BASIC_SPHERE=0,
        SPHERE=1,
        ROD=2,
        CUBE=3,
        ROTOR=4,
		RED_BLOOD_CELL=5,
		TRYP=6,
        ACTIVE_ROD=7,
        TWIST_ROD=8,
        BEAM=9,
        BENTBEAM=10,
		ACTIVE_ROD_AS=11,
        SHEET=12
    };
}

#endif
