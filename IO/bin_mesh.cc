#include <cstdlib>
#include <cstring>

#include "bin_mesh.hh"
#include <cmath>

bin_mesh::~bin_mesh() {
	if(tinfo!=NULL) {
		delete [] tinfo;
		delete [] finfo;
	}
}

void bin_mesh::sphere_fit(int t,float &K,float (&nf)[3][3]) {

	// get vertices
	float v[3][3];
	vertices(t,v);

	// normals
	float n[3][3];
	normals(t,n);

	// do fit
	sphere_fit(v,n,K,nf);
}

void bin_mesh::mean_vec(const float (&v)[3][3],float (&mv)[3]) {

	for(int i=0;i<3;i++) {
		mv[i]=0;
		for (int j=0;j<3;j++) mv[i] += v[j][i]/3;
	}
}

void bin_mesh::sphere_fit(const float (&v)[3][3],const float (&n)[3][3],
		float &K,float (&nf)[3][3]) {

	// mean normal
	float mn[3];
	mean_vec(n,mn);

	// mean vertex
	float mv[3];
	mean_vec(v,mv);

	// avg dot with normals
	float vn_dot=0;
	for(int i=0;i<3;i++) vn_dot += dot(v[i],n[i])/3;

	// curvature
	float num = 1-dot(mn,mn);
	float den = vn_dot-dot(mn,mv);
	K = den==0 ? 0 : (num/den);

	if (std::isnan(K)) {
		printf("%g / %g\n",num,den);

		printf("v:\n");
		for(int i=0;i<3;i++)
			printf("%g %g %g\n",v[i][0],v[i][1],v[i][2]);
		printf("n:\n");
		for(int i=0;i<3;i++)
			printf("%g %g %g\n",n[i][0],n[i][1],n[i][2]);

	}

	// fitted normal vecs
	for (int i=0;i<3;i++) for (int j=0;j<3;j++)
		nf[i][j] = mn[j] + K * (v[i][j]-mv[j]);
}

/*
float bin_mesh::solid_angle(int t,float K) {
	float num = 6*volume(t,false);
}
*/

/**
 * Print out surface area of t'th triangle
 */
float bin_mesh::area(int t,bool cap) {

	// check index
	if (t >= ntri || t < 0) {
		printf("triangle index %d out of range (%d triangles)\n",t,ntri);
		exit(-1);
	}

	// get sides
	float s[3][3];
	sides(t,s);

	// do cross product
	float cp[3];
	cross(s[0],s[1],cp);

	// area is half of cross product magnitude
	float A = 0.5*mag(cp);

	if (cap) {
		float K = mean_curvature(t);
//		A += A*K*K*side_mag(t)* 2./3; // AD HOC

		A *= (1 + K*K*side_mag(t) * 2./3);
	}

	return A;
}

/**
 * Print out the surface area of the mesh
 */
float bin_mesh::area(bool cap) {

	// total area
	float sum=0;

	// add from each triangle
	for(int t=0;t<ntri;t++) sum += area(t,cap);

	return sum;
}

/**
 * Return the volume of a single tetrahedron corresponding
 * to a triangle from the mesh
 *
 * wlog all tetrahedrons share a vertex at the origin
 */
float bin_mesh::volume(int t,bool cap) {

	// check index
	if (t >= ntri || t < 0) {
		printf("triangle index %d out of range (%d triangles)\n",t,ntri);
		exit(-1);
	}

	// grab vertices
	float v[3][3];
	vertices(t,v);

	// cross-product of second two vertices
	float cp[3];
	cross(v[1],v[2],cp);

	// vol is 1/6 of dot product of that with first vertex
	float V = dot(v[0],cp)/6;

	if (cap) {
		// get flat area and curvature
		float A = area(t,false), K = mean_curvature(t);

		V += A*K*side_mag(t)/3 * 6./5; // AD HOC factor 6/5
//		V *= (1 + K*K*side_mag(t));
	}

	return V;
}


/**
 * Return the volume of the mesh
 */
float bin_mesh::volume(bool cap) {

	// total volume
	float sum=0;
	for(int i=0;i<ntri;i++) sum += volume(i,cap);
	return sum;
}

float bin_mesh::mean_curvature() {

	double AK = 0;

	float At =0;
	float err1=0;
	float err2 = 0;
	for (int i=0;i<ntri;i++) {

		float K1 = mean_curvature(i);
		float A = area(i,false);

		AK += K1*A;

		float K2;
		float nf[3][3];
		sphere_fit(i,K2,nf);

		if (std::isnan(K2)) printf("%d\n",i);


		err1 += (K1-5)*A;
		err2 += (K2-5)*A;
		At    += A;

	}

	printf("%g %g\n",err1/At,err2/At);

	return AK/area(false);
}

/**
 * Return the mean curvature over the t'th triangle
 */
float bin_mesh::mean_curvature(int t) {

	// check index
	if (t >= ntri || t < 0) {
		printf("triangle index %d out of range (%d triangles)\n",t,ntri);
		exit(-1);
	}

	// grab three sides and normals
	float s[3][3],n[3][3];
	sides(t,s);
	normals(t,n);

	// calculate curvature along each side
	//
	// weight it by length of side
	//
	float w_kappa=0, S_tot=0;
	for (int i=0;i<3;i++) {

		float S = mag(s[i]);

		S_tot += S;

		// sometimes floats aren't large enough to
		// resolve side lengths
		if (S == 0) continue;

		// difference in normal vectors on this side
		float n_diff[3];
		for(int c=0;c<3;c++) n_diff[c] = n[(i+1)%3][c] - n[i][c];

		// curvature along this edge (dot with side itself)
		float kappa = (dot(n_diff,s[i]) / (S*S));

		// add into weighted total
		w_kappa += S*kappa;
	}

	// use geometric mean
	// float K = pow(kappa[0]*kappa[1]*kappa[2],1./3);
	
	// use arithmetic mean
	return w_kappa/S_tot;
}


void bin_mesh::load(FILE *inf) {
	int ibuf[2];

	// Do some sanity checks on the size of the mesh
	safe_fread(ibuf,sizeof(int),2,inf,"header information");
	if(*ibuf>16777216||ibuf[1]>16777216) {
		fputs("Mesh size out of range\n",stderr);
		exit(1);
	}

	// Set mesh constants
	nvert=*ibuf;ntri=ibuf[1];
	mem_vert=3*nvert;mem_tri=3*ntri;

	// Remove any memory that was previously allocated
	if(tinfo!=NULL) {
		delete [] tinfo;
		delete [] finfo;
	}

	// Read in the floating point information
	finfo=new float[2*mem_vert];
	safe_fread(finfo,sizeof(float),2*mem_vert,inf,"vertex vectors and normal vectors");

	// Read in the face information
	tinfo=new int[mem_tri];
	safe_fread(tinfo,sizeof(int),mem_tri,inf,"face indices");
}

void bin_mesh::draw_mesh_pov(int i,FILE *outf,bool normals) {

	// Print header information
	fprintf(outf,"#mesh2 {\nvertex_vectors {\n%d,\n",nvert);

	// Print vertex information
	float *fp=finfo;
	for(;fp<finfo+mem_vert;fp+=3) fprintf(outf,"<%g,%g,%g>\n",*fp,fp[1],fp[2]);

	// Print normal information (if needed)
	if(normals) {
		fprintf(outf,"}\nnormal_vectors {\n%d,\n",nvert);
		for(;fp<finfo+2*mem_vert;fp+=3) fprintf(outf,"<%g,%g,%g>\n",*fp,fp[1],fp[2]);
	}

	// Print triangle information
	fprintf(outf,"}\nface_indices {\n%d,\n",ntri);
	for(int *tp=tinfo;tp<tinfo+mem_tri;tp+=3) fprintf(outf,"<%d,%d,%d>\n",*tp,tp[1],tp[2]);
	fprintf(outf,"}\ntexture{t_msh%d}\n}\n",i);
}

void bin_mesh::draw_cylinders_pov(int i,FILE *outf) {

	fprintf(outf,"#declare cyls%d=union {\n",i);
	for(float *fp=finfo;fp<finfo+mem_vert;fp+=3) fprintf(outf,"sphere{<%g,%g,%g>,s}\n",*fp,fp[1],fp[2]);

	for(int *tp=tinfo;tp<tinfo+mem_tri;tp+=3) {
		cyl_print(outf,finfo+*tp*3,finfo+tp[1]*3);
		cyl_print(outf,finfo+*tp*3,finfo+tp[2]*3);
		cyl_print(outf,finfo+tp[1]*3,finfo+tp[2]*3);
	}
	fputs("}\n",outf);
}

void bin_mesh::draw_nodes_pov(int i,FILE *outf) {
	fprintf(outf,"#declare sphs%d=union {\n",i);
	for(float *fp=finfo;fp<finfo+mem_vert;fp+=3) fprintf(outf,"sphere{<%g,%g,%g>,r}\n",*fp,fp[1],fp[2]);
	fputs("}\n",outf);
}

void bin_mesh::draw_nodes(FILE *outf) {
	for(float *fp=finfo;fp<finfo+mem_vert;fp+=3) fprintf(outf,"%g %g %g\n",*fp,fp[1],fp[2]);
}

/** Prints a cylinder connecting two vertices to a file, skipping any with
 * duplicate start and end points, since they cause an error in POV-Ray.
 * \param[in] outf the file handle to write to.
 * \param[in] (p,q) pointers to the vertex pair to consider. */
void bin_mesh::cyl_print(FILE *outf,float *p,float *q) {
	char vbuf1[80],vbuf2[80];
	sprintf(vbuf1,"%g,%g,%g",*p,p[1],p[2]);
	sprintf(vbuf2,"%g,%g,%g",*q,q[1],q[2]);
	if(strcmp(vbuf1,vbuf2)!=0) fprintf(outf,"cylinder{<%s>,<%s>,s}\n",vbuf1,vbuf2);
}

/** \brief Opens a file and checks the operation was successful.
 *
 * Opens a file, and checks the return value to ensure that the operation
 * was successful.
 * \param[in] filename the file to open.
 * \param[in] mode the cstdio fopen mode to use.
 * \return The file handle. */
FILE* bin_mesh::safe_fopen(const char* filename,const char* mode) {
	FILE *temp=fopen(filename,mode);
	if(temp==NULL) {
		fprintf(stderr,"Error opening file \"%s\"",filename);
		exit(1);
	}
	return temp;
}

/** \brief Reads data from a file and checks that the operation was successful.
 *
 * Reads data from a file and checks the return value to ensure that the
 * operation was successful. If not successful, it prints an error message and
 * exits.
 * \param[in] ptr the memory to write to.
 * \param[in] size the size of each element to be read.
 * \param[in] count the number of elements to be read.
 * \param[in] fp the file handle to read from.
 * \param[in] p a description of what is being read, used in the error message
 *              if the operation isn't successful. */
void bin_mesh::safe_fread(void *ptr,size_t size,size_t count,FILE *fp,const char* p) {
	if(fread(ptr,size,count,fp)!=count) {
	    fprintf(stderr,"Can't read %s from file\n",p);
	    exit(1);
	}
}
