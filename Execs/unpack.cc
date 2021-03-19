#include <cstdio>
#include <cstring>

#include "bin_mesh.hh"

int main(int argc,char** argv) {

	// Check for the correct number of command-line arguments
	if(argc!=4) {
		fputs("Syntax: ./unpack <mode> <binary_contour_file> <output_filename>"
		      "\n\nMode: \"msh\" for POV-Ray mesh (with normals)\n"
		      "      \"mtr\" for POV-Ray mesh (with flat triangles, no normals)\n"
		      "      \"cyl\" for POV-Ray cylinders\n"
		      "      \"sph\" for POV-Ray spheres\n"
		      "      \"txt\" for plain text vertices\n\n"
		"If the output filename is \"-\" then the data will written to standard\n"
		"output\n",stderr);
		return 1;
	}

	// Find output mode
	int mode;
	if(strcmp(argv[1],"msh")==0) mode=0;
	else if(strcmp(argv[1],"mtr")==0) mode=1;
	else if(strcmp(argv[1],"cyl")==0) mode=3;
	else if(strcmp(argv[1],"sph")==0) mode=4;
	else if(strcmp(argv[1],"txt")==0) mode=5;
	else {
		fprintf(stderr,"unpack: Mode type \"%s\" not known\n",argv[1]);
		return 1;
	}

	// Open input and output files
	int ibuf[1],i;
	bin_mesh bm;
	FILE *inf=bm.safe_fopen(argv[2],"rb");
	bool std=strcmp(argv[3],"-")==0;
	FILE *outf=std?stdout:bm.safe_fopen(argv[3],"w");

	// Loop over the meshes in the input file, and do the requested output
	bm.safe_fread(ibuf,sizeof(int),1,inf,"number of meshes");
	for(i=0;i<*ibuf;i++) {
		bm.load(inf);
		if(bm.nvert==0) continue;
		switch(mode) {
			case 0: bm.draw_mesh_pov(i,outf);break;
			case 1: bm.draw_mesh_pov(i,outf,false);break;
			case 3: bm.draw_cylinders_pov(i,outf);break;
			case 4: bm.draw_nodes_pov(i,outf);break;
			case 5: bm.draw_nodes(outf);
		}
	}

	// Close and/or flush the input and output files
	if(!std) fclose(outf);
	else fflush(stdout);
	fclose(inf);
}
