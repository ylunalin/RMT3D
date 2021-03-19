#include <cstdlib>
#include <cstdio>
#include "stencil.hh"

/**
 * Construct stencil coefficients (hard-coded here)
 */
void stencil::construct_coeff_arrays() {

	// coefficients for finite difference stencils.
	// rows are order of accuracy (half thereof for
	// central and laplacian), and columns are the
	// actual coefficient values
	double f[n_derivs][n_orders][n_coeffs] = {
		{
			{	-1.,	1.},
			{	-1.5,	2.,	-0.5},
			{	-11./6.,	3.,	-1.5,	1./3.},
			{	-25./12.,	4.,	-3.,	4./3.,	-1./4.},
		}, {
			{	1.,	-2.,	1.},
			{	2.,	-5.,	4.,	-1.},
			{	35./12.,	-26./3.,	19./2.,	-14./3.,	11./12.},
			{	15./4.,	-77./6.,	107./6.,	-13.,	61./12.,	-5./6.}
		}
	};

	// backward stencils are determined by forward
	// (just a sign flip for odd derivatives plus
	// left-right reversal)
	double b[n_derivs][n_orders][n_coeffs];
	for (int i = 0; i < n_derivs; i++)
		for (int j = 0; j < n_orders; j++)
			for (int k = 0; k < j + 2 + i; k++)
				b[i][j][k] = (-1 + 2*i) * f[i][j][j+1+i-k];

	// central stencils are their own thing (could actually
	// get by summing one-way stencils, but ehh.)
	double c[n_derivs][n_orders][n_coeffs] = {
		{
			{},
			{	-0.5,	0,	0.5},
			{},
			{	1./12.,	-2./3.,	0,	2./3.,	-1./12.}
		}, {
			{},
			{	1.,	-2.,	1.},
			{},
			{	-1./12.,	4./3.,	-5./2.,	4./3.,	-1./12.}
		}
	};

	for (int i = 0; i < n_derivs; i++)
		for (int j = 0; j < n_orders; j++)
			for (int k = 0; k < n_coeffs; k++) {
				forward_coeffs[i][j][k] = f[i][j][k];
				backward_coeffs[i][j][k] = b[i][j][k];
				centered_coeffs[i][j][k] = c[i][j][k];
	}
}

/**
 * General stencil constructor with various prefactors
 */
stencil::stencil(int ord_,int deriv_,type dir_,cart kind_,
		int sm4_,int sn4_,double pre_x,double pre_y,double pre_z) :
		order(ord_),deriv(deriv_),direction(dir_),kind(kind_) {

	// load in coefficients
	construct_coeff_arrays();

	// load prefactor array
	pres[0] = pre_x;
	pres[1] = pre_y;
	pres[2] = pre_z;

	// fill in stencil grid points
	if (kind == Laplacian) {
		set_laplace_stencil(sm4_,sn4_);
	} else {
		set_1D_stencil(sm4_,sn4_);
	}
}

stencil::stencil() {
	construct_coeff_arrays();
	points = NULL;
}

/**
 * Helper method for constructing and arranging grid points to
 * obtain 1D stencil of desired type / order / cartesian direction
 */
void stencil::set_1D_stencil(int sm4, int sn4) {

	// check that we have correct stencil here
	if (order > 4) {
		printf("error: only stencils up to order 4 are implemented.\n");
		exit(-1);
	} else if (direction == centered && order % 2 == 1) {
		printf("error: centered difference requires even order.\n");
		exit(-1);
	}

	// relative step in flattened array for one step in cart. directions
	int dels[3] = {1,sm4,sm4*sn4};

	// variables for making new grid points
	int rel_step;
	double pref;

	// if we're not doing a Laplacian, figure out which
	// cartesian index we're using, our starting ste dist.,
	// and which set of stencil coefficients we're using
	int cart_ind = -1;
	int start_ind = -1;
	double (*coeffs)[n_coeffs] = NULL;

	// which cartesian direction?
	switch (kind) {
		case X: cart_ind = 0; break;
		case Y: cart_ind = 1; break;
		case Z: cart_ind = 2; break;
		default:
			puts("bad cartesian enum for stencil");
			exit(-1);
			break;
	}

	// which stencil type?
	switch (direction) {

		case forward:
			start_ind = 0;
			coeffs = &forward_coeffs[deriv-1][order-1];
			npts = order + deriv;
			break;

		case backward:
			coeffs = &backward_coeffs[deriv-1][order-1];
			npts = order + deriv;
			start_ind = 1 - npts;
			break;

		case centered:
			start_ind = -order/2;
			coeffs = &centered_coeffs[deriv-1][order-1];
			npts = order + 1;
			break;

		default:
			puts("bad direction enum for stencil");
			exit(-1);
			break;
	}

	// make array of grid points
	points = new grid_pt[npts];

	// now, for each point in the stencil...
	for (int i = 0; i < npts; i++) {

		// get relative step and prefactor...
		rel_step = (start_ind + i)*dels[cart_ind],
		pref = pres[cart_ind] * (*coeffs)[i];

		// ...and make a grid point
		points[i] = grid_pt(rel_step,pref);
	}
}

/**
 * Helper method for constructing and arranging grid points
 * to obtain Laplace stencil of desired order
 */
void stencil::set_laplace_stencil(int sm4, int sn4) {

	// check that we have correct stencil here
	if (order > 4) {
		printf("error: only stencils up to order 4 are implemented.\n");
		exit(-1);
	} else if (deriv != 2) {
		printf("error: Laplace operator is second derivative.\n");
		exit(-1);
	} else if (order % 2 == 1) {
		printf("error: Laplace operator requires even order of accuracy.\n");
		exit(-1);
	}

	// set up number of points, get pointer to coeffs
	npts = 1 + 3*order;
	double (*coeffs)[n_coeffs] = &centered_coeffs[deriv-1][order-1];

	// relative step in flattened array for one step in cart. directions
	int dels[3] = {1,sm4,sm4*sn4};

	// make array of grid points
	points = new grid_pt[npts];

	// variables for making new grid points
	int rel_step;
	double pref;

	// instantiate relative step dist. and points array index
	int p = 0;
	int start_ind = -order/2;

	// total up prefactor of center point as we go
	double center_pref = 0;

	// for each cartesian direction...
	for (int cart_ind = 0; cart_ind < 3; cart_ind++) {

		// and each point in the 1D stencil...
		for (int i = 0; i < order + 1; i++) {

			// if we're not in the center...
			if (i != -start_ind) {

				// calculate relative step and prefactor
				rel_step = (start_ind + i)*dels[cart_ind];
				pref = pres[cart_ind] * (*coeffs)[i];

				// make the gridpoint
				points[p] = grid_pt(rel_step,pref);
				p++;

			} else {

				// otherwise, add to our center prefactor
				center_pref += pres[cart_ind] * (*coeffs)[i];
			}
		}
	}

	// after going through, last gridpoint is center
	points[npts-1] = grid_pt(0,center_pref);
}

/**
 * Helper function to print stencil for debugging.
 */
void stencil::print_stencil(){

	// just print each point
	for (int i=0; i<npts; i++){
		printf("Point relative to differentiation pt %d,"
			" prefactor %f\n", points[i].rel_step, points[i].pre);
	}
}

stencil& stencil::operator=(const stencil &that) {

	// check for self-assignment
	if (this == &that) return *this;

	// delete points if allocated
	if (points != NULL) delete[] points;

	//  reallocate and copy values
	points = new grid_pt[that.npts];
	for (int i = 0; i < that.npts; i++)
		points[i] = that.points[i];

	// copy all other values
	// (coeff arrays already constructed in default constructor)
	order = that.order;
	deriv = that.deriv;
	direction = that.direction;
	kind = that.kind;
	npts = that.npts;
	for (int i = 0; i < 3; i++) pres[i] = that.pres[i];

	return *this;
}

void stencil::avg_shift(int step) {

	grid_pt *newpts = new grid_pt[2*npts];
	for (int i = 0; i < 2*npts; i++) {
		newpts[i] = points[i%npts];
		newpts[i].pre *= 0.5;
		if (i >= npts) newpts[i].rel_step += step;
	}
	npts *= 2;

	delete[] points;
	points = newpts;
}

/**
 * Driver method
 */
/*
int main() {

	stencil* st = NULL;
	int order, deriv;

	deriv = 2;
	order = 2;
	printf("\norder %d, deriv %d, Laplace, centered stencil:\n\n",
		order,deriv);
	st = new stencil(order, Laplacian,104,104,1.);
	st->print_stencil();
	delete st;

	order = 4;
	printf("\norder %d, deriv %d, Laplace, centered stencil:\n\n",
		order,deriv);
	st = new stencil(order, Laplacian,104,104,1.);
	st->print_stencil();
	delete st;

	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 2; j++) {
			order = i + 1;
			deriv = j + 1;

			printf("\norder %d, deriv %d, X-dir, forward stencil:\n\n",
				i+1,j+1);
			st = new stencil(order, deriv,
				forward,X,104,104,1.);
			st->print_stencil();
			delete st;

			printf("\norder %d, deriv %d, X-dir, backward stencil:\n\n",
				i+1,j+1);
			st = new stencil(i + 1, j+1,
				backward,X,104,104,1.);
			st->print_stencil();
			delete st;

			if (order % 2 == 0) {
				printf("\norder %d, deriv %d, X-dir, centered stencil:\n\n",
				i+1,j+1);
				st = new stencil(i + 1, j+1,
				centered,X,104,104,1.);
				st->print_stencil();
				delete st;
			}

			printf("\norder %d, deriv %d, Y-dir, forward stencil:\n\n",
				i+1,j+1);
			st = new stencil(order, deriv,
				forward,Y,104,104,1.);
			st->print_stencil();
			delete st;

			printf("\norder %d, deriv %d, Y-dir, backward stencil:\n\n",
				i+1,j+1);
			st = new stencil(i + 1, j+1,
				backward,Y,104,104,1.);
			st->print_stencil();
			delete st;

			if (order % 2 == 0) {
				printf("\norder %d, deriv %d, Y-dir, centered stencil:\n\n",
				i+1,j+1);
				st = new stencil(i + 1, j+1,
				centered,Y,104,104,1.);
				st->print_stencil();
				delete st;
			}

			printf("\norder %d, deriv %d, Z-dir, forward stencil:\n\n",
				i+1,j+1);
			st = new stencil(order, deriv,
				forward,Z,104,104,1.);
			st->print_stencil();
			delete st;

			printf("\norder %d, deriv %d, Z-dir, backward stencil:\n\n",
				i+1,j+1);
			st = new stencil(i + 1, j+1,
				backward,Z,104,104,1.);
			st->print_stencil();
			delete st;

			if (order % 2 == 0) {
				printf("\norder %d, deriv %d, Z-dir, centered stencil:\n\n",
				i+1,j+1);
				st = new stencil(i + 1, j+1,
				centered,Z,104,104,1.);
				st->print_stencil();
				delete st;
			}

		}
	}
}
*/
