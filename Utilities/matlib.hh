#ifndef MATLIB
#define MATLIB

#include <cstdio>
#include <cmath>

struct sym_matrix {

	public:
	/** values */
	double v[6];

	sym_matrix() {
//		for (int i = 0; i < 6; i++) v[i] = 0;
	};
	sym_matrix(double d) {
		for (int i = 0; i < 6; i++) v[i] = 0;
		v[0] = v[3] = v[5] = d;
	};
	sym_matrix(double d1,double d2,double d3) {
		for (int i = 0; i < 6; i++) v[i] = 0;
		v[0] = d1;
		v[3] = d2;
		v[5] = d3;
	};
	sym_matrix(double a11,double a12,double a13,
			double a22,double a23,double a33) {
		v[0]=a11;
		v[1]=a12;
		v[2]=a13;
		v[3]=a22;
		v[4]=a23;
		v[5]=a33;
	};
	inline int ind(int i,int j) const {
		if (i > j) {
			int tmp = j;
			j = i;
			i = tmp;
		}
		//XXX is this inefficient?
		return j + i*(5-i)/2;
	};
	inline double get(int i,int j) {
		return v[ind(i,j)];
	};
	inline void set(int i,int j,double v_) {
		v[ind(i,j)] = v_;
	};
	inline sym_matrix trans() {
		return *this;
	};
	inline void scale(double s) {
		for (int i = 0; i < 6; i++) v[i] *= s;
	};
	inline double frobenius(){
		double tmp=0;
		for(int j=0;j<3;j++) for(int i=0;i<3;i++) tmp+=get(i,j)*get(i,j);
		return sqrt(tmp);
	}
	inline sym_matrix times(double s) {
		sym_matrix m2;
		for (int i = 0; i < 6; i++) m2.v[i] = v[i]*s;
		return m2;
	};
	inline sym_matrix plus(sym_matrix m1) {
		sym_matrix m2;
		for (int i = 0; i < 6; i++) m2.v[i] = v[i] + m1.v[i];
		return m2;
	};
	inline sym_matrix minus(sym_matrix m1) {
		sym_matrix m2;
		for (int i = 0; i < 6; i++) m2.v[i] = v[i] - m1.v[i];
		return m2;
	};
	inline sym_matrix operator- (const sym_matrix& m1) {
		return minus(m1);
	};
	inline sym_matrix operator* (const double& s) {
		return times(s);
	};
	inline sym_matrix times(sym_matrix m1) {
		sym_matrix m2;
		for (int i = 0; i < 3; i++) {
			for (int j = i; j < 3; j++) {
				double val = 0;
				for (int k = 0; k < 3; k++) val += get(i,k)*m1.get(k,j);
				m2.set(i,j,val);
			}
		}
		return m2;

	};
	inline void get_row(int i,double (&d)[3]) {
//		print();
//		printf("getting row %d\n",i);
//		printf("start: [%f,%f,%f]\n",d[0],d[1],d[2]);
		for (int j = 0; j < 3; j++) d[j] = get(i,j);
//		printf("got: [%f,%f,%f]\n",d[0],d[1],d[2]);
	};
	inline void set_row(int i,double (&d)[3]) {
		for (int j = 0; j < 3; j++) set(i,j,d[j]);
	};
	inline double trace() {
		return v[0]+v[3]+v[5];
	};
	inline void print() {
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				printf(" %12.6g ",get(i,j));
			}
			printf("\n");
		}
	};

	// subscripting
	inline double& operator() (int row,int col) {
		return v[ind(row,col)];
	};
	inline double operator() (int row,int col) const {
		return v[ind(row,col)];
	};
	inline double det() {
		return v[0]*(v[3]*v[5] - v[4]*v[4]) +
			v[1]*(2*v[2]*v[4]-v[1]*v[5]) - v[2]*v[2]*v[3];
	};
};

struct matrix {

	public:
	/** values */
	double v[9];

	matrix() {
	//	for (int i = 0; i < 9; i++) v[i] = 0;
	};
	matrix(double d) {
		for (int i = 0; i < 9; i++) v[i] = 0;
		v[0] = v[4] = v[8] = d;
	};
	matrix(double d1,double d2,double d3) {
		for (int i = 0; i < 9; i++) v[i] = 0;
		v[0] = d1;
		v[4] = d2;
		v[8] = d3;
	};
	matrix(double a11,double a12,double a13,double a21,double a22,
		double a23,double a31,double a32,double a33) {
		v[0] = a11;
		v[1] = a12;
		v[2] = a13;
		v[3] = a21;
		v[4] = a22;
		v[5] = a23;
		v[6] = a31;
		v[7] = a32;
		v[8] = a33;
	};
	void initialize(double a11,double a12,double a13,double a21,double a22,
		double a23,double a31,double a32,double a33) {
		v[0] = a11;
		v[1] = a12;
		v[2] = a13;
		v[3] = a21;
		v[4] = a22;
		v[5] = a23;
		v[6] = a31;
		v[7] = a32;
		v[8] = a33;
	};
	inline int ind(int i,int j) const {
		return 3*i + j;
	};
	inline double get(int i,int j) {
		return v[ind(i,j)];
	};
	inline void set(int i,int j,double v_) {
		v[ind(i,j)] = v_;
	};
	inline double frobenius(){
		double tmp=0;
		for(int j=0;j<3;j++) for(int i=0;i<3;i++) tmp+=get(i,j)*get(i,j);
		return sqrt(tmp);
	}
	inline matrix trans() {
		matrix m;
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
				m.set(i,j,get(j,i));
		return m;
	};
	inline matrix plus(const matrix& m1) {
		matrix m2;
		for (int i = 0; i < 9; i++) m2.v[i] = v[i] + m1.v[i];
		return m2;
	};
	inline matrix minus(const matrix& m1) {
		matrix m2;
		for (int i = 0; i < 9; i++) m2.v[i] = v[i] - m1.v[i];
		return m2;
	};
	inline matrix times(matrix m1) {
		matrix m2;
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				double val = 0;
				for (int k = 0; k < 3; k++) val += get(i,k)*m1.get(k,j);
				m2.set(i,j,val);
			}
		}
		return m2;
	};
	inline sym_matrix sym_times(matrix m1) {
		sym_matrix m2;
		for (int i = 0; i < 3; i++) {
			for (int j = i; j < 3; j++) {
				double val = 0;
				for (int k = 0; k < 3; k++) val += get(i,k)*m1.get(k,j);
				m2.set(i,j,val);
			}
		}
		return m2;
	};
	inline matrix operator+(const matrix& m1) {return plus(m1);};
	inline matrix operator-(const matrix& m1) {return minus(m1);};
	inline sym_matrix aat() {
		return (*this).sym_times(this->trans());
	};
	inline sym_matrix ata() {
		return this->trans().sym_times(*this);
	};
	inline double trace() {
		return v[0]+v[4]+v[8];
	};
	inline void copy(matrix &m) {
		for (int i = 0; i < 9; i++) v[i] = m.v[i];
	};
	inline void print() {
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				printf(" %12.6g ",get(i,j));
			}
			printf("\n");
		}
	};
	inline void lu_decomp(matrix &l,matrix &u,matrix &p) {
		matrix eye(1);
		u.copy(*this);
		l.copy(eye);
		p.copy(eye);

		double tmp;
		double max;
		int ind;
		for (int j = 0; j < 3; j++) {
			max = 0;
			ind = -1;
			for (int i = j; i < 3; i++) {
				if (fabs(u.get(i,j)) >= max) {
					ind = i;
					max = fabs(u.get(i,j));
				}
			}

			for (int k = 0; k < 3; k++) {
				if (k < j) {
					tmp = l.get(j,k);
					l.set(j,k,l.get(ind,k));
					l.set(ind,k,tmp);
				} else {
					tmp = u.get(j,k);
					u.set(j,k,u.get(ind,k));
					u.set(ind,k,tmp);
				}

				tmp = p.get(j,k);
				p.set(j,k,p.get(ind,k));
				p.set(ind,k,tmp);
			}

			for (int i = j+1; i < 3; i++) {
				l.set(i,j,u.get(i,j)/u.get(j,j));
				for (int k = j; k < 3; k++) {
					u.set(i,k,u.get(i,k)-l.get(i,j)*u.get(j,k));
				}
			}
		}
	};
	inline void get_col(int j,double (&d)[3]) {
		for (int i = 0; i < 3; i++) d[i] = get(i,j);
	};
	inline void set_col(int j,double (&d)[3]) {
		for (int i = 0; i < 3; i++) set(i,j,d[i]);
	};
	inline void get_row(int i,double (&d)[3]) {
		for (int j = 0; j < 3; j++) d[j] = get(i,j);
	};
	inline void set_row(int i,double (&d)[3]) {
		for (int j = 0; j < 3; j++) set(i,j,d[j]);
	};

	inline matrix inverse() {
/*
		matrix l,u,p,inv;
		lu_decomp(l,u,p);

		double pcol[3];
		double icol[3];
		for (int i = 0; i < 3; i++) {
			p.get_col(i,pcol);
			solve(l,u,pcol,icol);
			inv.set_col(i,icol);
		}
		return inv;
*/
		matrix m2;
		double d = det();
		for (int j = 0; j < 3; j++) for (int i = 0; i < 3; i++) {

			// THIS IS SETTING THE TRANSPOSE OF THE COFACTOR MATRIX
			m2.v[ind(j,i)] =
				v[ind((i+1)%3,(j+1)%3)]*v[ind((i+2)%3,(j+2)%3)] -
				v[ind((i+1)%3,(j+2)%3)]*v[ind((i+2)%3,(j+1)%3)];

			m2.v[ind(j,i)] /= d;
		}
		return m2;
	};
	inline static void bsub(matrix &u, double (&b)[3], double (&x)[3]) {
		for (int i = 2; i >= 0; i--) {
			x[i] = b[i];
			for (int j = i+1; j < 3; j++) x[i] -= u.get(i,j)*x[j];
			x[i] /= u.get(i,i);
		}
	};
	inline static void fsub(matrix &l, double (&b)[3], double (&x)[3]) {
		for (int i = 0; i < 3; i++) {
			x[i] = b[i];
			for (int j = 0; j < i; j++) x[i] -= l.get(i,j)*x[j];
			x[i] /= l.get(i,i);
		}
	};
	inline static void solve(matrix &l,matrix &u,double (&b)[3],
			double (&x)[3]) {

		double y[3];

		// LUx = b, Ux = y, Ly = b
		fsub(l,b,y);
		bsub(u,y,x);
	};
	inline sym_matrix dbl_sym_part() {
		sym_matrix m(2*v[0],v[1]+v[3],v[2]+v[6],
				2*v[4],v[5]+v[7],2*v[8]);
		return m;
	};

	// subscripting
	inline double &operator() (int row,int col) {
		return v[ind(row,col)];
	};
	inline double operator() (int row,int col) const {
		return v[ind(row,col)];
	};
	inline double det() {
		return v[0]*(v[4]*v[8] - v[5]*v[7])
			+ v[3]*(v[2]*v[7] - v[1]*v[8])
			+ v[6]*(v[1]*v[5] - v[2]*v[4]);
	};
};

#endif
