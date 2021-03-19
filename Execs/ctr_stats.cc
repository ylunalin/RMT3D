#include "bin_mesh.hh"

#include <cstdio>
#include <cstdlib>
#include <cstring>

int main(int argc,char** argv) {

	if(argc < 3 || argc > 4) {
		fputs("Usage: ./ctr_stats [-c] index infile\n",stderr);
		fputs("  -c: suppress contribution from surface curvature\n",stderr);
		exit(0);
	}

	// filename from args
	bool no_cap = strcmp(argv[1],"-c")==0;
	int index = atoi(argv[no_cap?2:1]);
	const char *fname = argv[no_cap?3:2];

	// open file, load bin mesh
	bin_mesh mesh;
	FILE *fh = mesh.safe_fopen(fname,"rb");

	// mesh form:
	// 1 integer w/number of meshes, then each mesh continuous

	// get number of meshes
	int n_meshes;
	mesh.safe_fread(&n_meshes,sizeof(int),1,fh,"number of meshes");
	double tot_area=0;
	double tot_vol=0;
	for(int i=0;i<n_meshes;i++) {

		// load mesh
		mesh.load(fh);

		// print area/volume to stdout
//		printf("%d %g %g %g\n", index, mesh.area(!no_cap),mesh.volume(!no_cap),mesh.mean_curvature());
		tot_area += mesh.area(!no_cap);
		tot_vol  += mesh.volume(!no_cap);
	}
	printf("%d %g %g\n", index, tot_area, tot_vol);

	fclose(fh);
}
