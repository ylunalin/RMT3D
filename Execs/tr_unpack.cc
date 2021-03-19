#include <cstdio>
#include <cstdlib>
#include <cstring>

FILE* safe_fopen(const char* filename,const char* mode) {
	FILE *temp=fopen(filename,mode);
	if(temp==NULL) {
		fprintf(stderr,"tr_unpack: Error opening file \"%s\"",filename);
		exit(1);
	}
	return temp;
}

void safe_fread(void *ptr,size_t size,size_t count,FILE *fp,const char* p) {
	if(fread(ptr,size,count,fp)!=count) {
	    fprintf(stderr,"tr_unpack: Can't read %s from file\n",p);
	    exit(1);
	}
}

int main(int argc,char** argv) {

	// Check for the correct number of command-line arguments
	if(argc!=4) {
		fputs("Syntax: ./tr_unpack <mode> <binary_tracer_file> <output_filename>"
		      "\n\nMode: \"sph\" for POV-Ray spheres\n"
		      "      \"txt\" for plain text\n\n"
		"If the output filename is \"-\" then the data will written to standard\n"
		"output\n",stderr);
		return 1;
	}

	// Find output mode
	int mode;
	if(strcmp(argv[1],"sph")==0) mode=0;
	else if(strcmp(argv[1],"txt")==0) mode=1;
	else {
		fprintf(stderr,"tr_unpack: Mode type \"%s\" not known\n",argv[1]);
		return 1;
	}

	// Open input and output files
	int ibuf[1],&n=*ibuf;
	FILE *inf=safe_fopen(argv[2],"rb");
	bool std=strcmp(argv[3],"-")==0;
	FILE *outf=std?stdout:safe_fopen(argv[3],"w");

	// Loop over the meshes in the input file, and do the requested output
	safe_fread(ibuf,sizeof(int),1,inf,"number of tracers");
	if(n<0||n>16777216) {
		fputs("tr_unpack: Number of tracers out of range\n",stderr);
		return 1;
	}

	float *tout=new float[3*n],*tp=tout,*te=tout+3*n;
	safe_fread(tout,sizeof(float),3*n,inf,"tracer positions");
	if(mode==0) {
		for(;tp<te;tp+=3) fprintf(outf,"sphere{<%g,%g,%g>,r}\n",*tp,tp[1],tp[2]);
	} else for(;tp<te;tp+=3) fprintf(outf,"%g %g %g\n",*tp,tp[1],tp[2]);

	// Close and/or flush the input and output files
	if(!std) fclose(outf);
	else fflush(stdout);
	fclose(inf);
    delete [] tout;
}
