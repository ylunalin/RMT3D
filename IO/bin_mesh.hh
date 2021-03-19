#ifndef BIN_MESH_HH
#define BIN_MESH_HH

#include <cstdio>
#include <cmath>

class bin_mesh {
	public:
		/** number of vertices */
		int nvert;
		/** number of triangles */
		int ntri;
		/** number of floats of vertex storage */
		int mem_vert;
		/** number of ints of triangle storage */
		int mem_tri;
		/** vertex values and normal vectors */
		float *finfo;
		/** triangle values */
		int *tinfo;

		bin_mesh() : finfo(NULL), tinfo(NULL) {}
		~bin_mesh();

		/** print out solid angle subtended by triangle
		 *  if vertices are assumed to lie on a sphere
		 *  with the provided curvature */
		float solid_angle(float K);
		/** print out average mean curvature value */
		float mean_curvature();
		/** get sphere fit params for a single triangle */
		void sphere_fit(int t,float &K,float (&nf)[3][3]);
		/** print out average mean curvature value */
		static void sphere_fit(const float (&v)[3][3],const float (&n)[3][3],
				float &K,float (&nf)[3][3]);
		/** print out the volume enclosed by the mesh */
		float volume(bool cap=true);
		/** print out the surface area of the mesh */
		float area(bool cap=true);


		void load(FILE *inf);
		void draw_mesh_pov(int i,FILE *outf,bool normals=true);
		void draw_cylinders_pov(int i,FILE *outf);
		void draw_nodes_pov(int i,FILE *outf);
		void draw_nodes(FILE *outf);
		void safe_fread(void *ptr,size_t size,size_t count,FILE *fp,const char* p);
		FILE* safe_fopen(const char* filename,const char* mode);
	private:
		/** calculate volume of one triangle */
		float volume(int,bool);
		/** calculate surface area of one triangle */
		float area(int,bool);
		/** calculate curvature * area over one triangle */
		float mean_curvature(int);
		/** calculate sides "quadratic magnitude" */
		float side_mag(int t) {
			float s[3][3];
			sides(t,s);

			return (dot(s[0],s[0]) + dot(s[1],s[1]) + dot(s[2],s[2]))/8;
		}

		void cyl_print(FILE *outf,float *p,float *q);

		/**
		 * grab vertices of t'th triangle, with
		 * v[i][c] being c'th component of i'th vertex
		 */
		void vertices(int t,float (&v)[3][3]) {

			int *tp = tinfo+3*t; // point to triangle

			// read in
			for(int i=0;i<3;i++) {
				for(int c=0;c<3;c++) {
					v[i][c] = finfo[3*tp[i]+c];
				}
			}
		}

		/**
		 * grab all three side vectors end-to-end
		 *
		 * i.e. end point of one is start point of next
		 *
		 * s[i][c] is c'th component of i'th side
		 *
		 * start point of s[i] is i'th vertex
		 */
		void sides(int t,float (&s)[3][3]) {

			// grab vertices
			float v[3][3];
			vertices(t,v);

			for(int i=0; i<3; i++) {

				// look at "next" vertex (cyclic perm)
				int j = (i+1)%3;
				for (int c=0;c<3;c++) s[i][c] = v[j][c]-v[i][c];
			}
		}

		/**
		 * grab the normal vector at each endpoint of triangle t
		 *
		 * n[i][c] is the c'th component of the normal at the i'th vertex
		 */
		void normals(int t,float (&n)[3][3]) {

			int *tp = tinfo+3*t; // point to triangle

			// read in
			for(int i=0;i<3;i++) {
				for(int c=0;c<3;c++) {

					// normals are second half of vertex array
					n[i][c] = finfo[3*tp[i]+c+mem_vert];
				}

				normalize(n[i]);
			}
		}

		static float dot(const float (&v1)[3],const float (&v2)[3]) {

			// do dot product
			float sum=0;
			for(int i=0;i<3;i++) sum += v1[i]*v2[i];
			return sum;
		}
		static void cross(const float (&v1)[3],const float(&v2)[3],
				float (&out)[3]) {

			// do cross product
			for(int i=0;i<3;i++) {

				// just a cyclic permutation
				int j = (i+1)%3;
				int k = (i+2)%3;

				out[i] = v1[j]*v2[k] - v2[j]*v1[k];
			}
		}
		static float mag(const float (&v)[3]) {return sqrt(dot(v,v));}
		static void normalize (float (&v)[3]) {
			float vmag = mag(v);
			for(int i=0;i<3;i++) v[i] /= vmag;
		}
		static void mean_vec(const float (&n)[3][3],float (&nm)[3]);
};

#endif
