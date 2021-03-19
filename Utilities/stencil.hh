#include <cstdlib>
#ifndef STENCIL_HH
#define STENCIL_HH

/**
 * Data type representing a point in the stencil, containing the relative
 * step from the differentiation point to here and the prefactor by which
 * it should be multiplied.
 */
struct grid_pt{

    // relative step to the derivative point, in flattened array
	int rel_step;

    // prefactor
	double pre;

	// struct constructors
	grid_pt(int rel_step_, double pre_):rel_step(rel_step_), pre(pre_){}
	grid_pt(){}
};

/**
 * enum for stencil direction types
 */
enum type {forward, backward, centered};

/**
 * enum for stencil cartesian direction types
 */
enum cart {X,Y,Z,Laplacian};

/** A class representing a finite difference stencil. */
class stencil {

	public:

	// coefficients for finite difference stencils.
	// rows are order of accuracy (half thereof for
	// central and laplacian), and columns are the
	// actual coefficient values
	static const int n_derivs = 2;
	static const int n_orders = 4;
	static const int n_coeffs = 6;

	double forward_coeffs[n_derivs][n_orders][n_coeffs];
	double backward_coeffs[n_derivs][n_orders][n_coeffs];
	double centered_coeffs[n_derivs][n_orders][n_coeffs];

	// stencil's determined by order of convergence, direction and derivative
	int order;
	int deriv;
	type direction;
	cart kind;

	// number of points in stencil
	int npts;

	// prefactor in each cartesian direction
	double pres[3];

	// the stencil itself
	grid_pt* points;

	// constructor/destructor

    // constructors with three prefix vals (default type is centered diff)
	stencil(int ord_,int deriv_,type dir_,cart kind_,int sm4_,int sn4_,
			double pre_x, double pre_y, double pre_z);
	stencil();

	~stencil() {if (points != NULL) delete[] points;};

	// setup and print functions
	void construct_coeff_arrays();
	void set_1D_stencil(int sm4, int sn4);
	void set_laplace_stencil(int sm4, int sn4);
	void print_stencil();
	void avg_shift(int step);
	stencil& operator=(const stencil &that);
};
#endif
