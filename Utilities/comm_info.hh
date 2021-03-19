#ifndef F3D_COMM_INFO
#define F3D_COMM_INFO
#include "geometry.hh"
/**
 * Struct containing useful communication info
 */
struct comm_info {
	/** x-dimension of communicated region */
	int m;
	/** y-dimension of communicated region */
	int n;
	/** z-dimension of communicated region */
	int o;
	/** total number of nodes in communicated region */
	int len;
	/** index of first node of sent region */
	int s0;
	/** index of first node of received region */
	int r0;
	/** index of first node of ghost region adjacent to sent region */
	int g0;
	/** tag of send message */
	int s_tag;
	/** tag of recv message */
	int r_tag;

	void dump(){
		printf("(%dx%dx%d)=%d;\ts0=%d, r0=%d, g0=%d, s_tag=%d, r_tag=%d\n", m, n, o, len, s0, r0, g0, s_tag, r_tag);
	}
	void reset(){
		m = n = o = len =s0= r0 = g0 =s_tag=r_tag =0;
	}
};

/** The comm table that uses comm_info struct */
class mpi_table{
	public:
	const int sm,sn,so;
	const int sm4,sn4,smn4;
	comm_info ct[27];

	mpi_table (geometry *gm):
	sm(gm->sm), sn(gm->sn), so(gm->so),
	sm4(sm+4), sn4(sn+4), smn4(sm4*sn4)
	{
        for(int i=0; i<27; i++) ct[i].reset();
    }
	~mpi_table() {}

	/**
	 * A helper function to generate x, y, z tags (a number
	 * between 0 and 2 inclusive) for 27 cubes
	 * and determine if a certain cube is neighboring a
	 * face of the center cube
	 * \param[in] tag the integer number of neighbor, 0 is the neighbor sharing ai, aj, ak
	 * \param[in, out] (ti, tj, tk) tags in x, y, z direction
	 */
	int nface(int tag, int &ti, int &tj, int &tk){
		tk = tag / 9;       // z index just int divide
		ti = tag % 3;       // x index just mod
		tj = (tag / 3) % 3; // y index do both
		// return face number, if it's not a face return -1
		// Compute if the neighbor is an orthogonal neighbor, if it is, \sum_l |t_l-1| = 1
		int num_ones = 0;
		if(ti==1) num_ones++;
		if(tj==1) num_ones++;
		if(tk==1) num_ones++;
		if (num_ones == 2)
			return (ti==0)?0:((ti==2)?1:((tj==0)?2:((tj==2)?3:((tk==0)?4:5))));
		else return -1;
	}

	/**
	 * Generates the lookup table for communications between this
	 * processor and its neighbors. Table is 2D, with a row
	 * corresponding to communication tag and column
	 * corresponding to info field. Columns are (left to right):
	 * - length of comm buffer array (LEN)
	 * - x-dim length of region (M)
	 * - y-dim length of region (N)
	 * - z-dim length of region (O)
	 * - integer which must be added to umem pointer to reach
	 *   starting index of region which is be sent (SEND)
	 * - integer which must be added to umem pointer to reach
	 *   starting index of region which is be filled (RECV)
	 *
	 * As table is generated, buffer lengths are compared and
	 * the class variable max_buf is set to be equal to the largest
	 * region size (except for tag=13)
	 */
	int setup_table(geometry *gm){
		// will add up combined length and check that fits in buffer
		int len = 0;

		// go through and allocate each row by tag
		for (int tag = 0; tag < 27; tag++) {

			// allocate set of columns
			ct[tag].s_tag = tag;
			ct[tag].r_tag = 26 - tag;

			// get 0,1,2 indices for tag
			int ti,tj,tk;
			nface(tag,ti,tj,tk);

			// store dimensions
			ct[tag].m = (ti==1)?gm->sm:2;
			ct[tag].n = (tj==1)?gm->sn:2;
			ct[tag].o = (tk==1)?gm->so:2;

			// store length
			ct[tag].len = ct[tag].m * ct[tag].n * ct[tag].o;

			// check fits in buffer if not middle region
			// multiple of 2 -> (3 vels + 1 pressure or 18 face velocities,
			// + one object id int cast as a double) x 2 way comms
			if (tag != 13)
				len += 2*38*ct[tag].len;

			int sm4 = sm+4;
			int sn4 = sn+4;
			int smn4 = sm4 * (sn+4);
			// store start
			ct[tag].s0 =
				((ti<2)?0:(sm-2)) +
				sm4*((tj<2)?0:(sn-2)) +
				smn4*((tk<2)?0:(so-2)) + 2*(1+sm4*(1+sn4));

			// store receive
			ct[tag].g0 =
				((ti<2)?(2*ti):(sm+2)) +
				sm4*((tj<2)?(2*tj):(sn+2)) +
				smn4*((tk<2)?(2*tk):(so+2));

			// the start of ghost region of the receiving neighbor
			ct[tag].r0 =
				((ti>0)?(2*(2-ti)):(sm+2)) +
				sm4*((tj>0)?(2*(2-tj)):(sn+2)) +
				smn4*((tk>0)?(2*(2-tk)):(so+2));
		}
		return len;
	}
	void copy_table( comm_info (&out)[27]){
		for(int i=0;i<27;i++){
			out[i] = ct[i];
		}
	}
};
#endif
