#ifndef F3D_EXTRAP_MMAP
#define F3D_EXTRAP_MMAP
#include <cstdio>
#include <cstdlib>
#include <vector>
#include "defs.hh"
#include "fields.hh"
#include "geometry.hh"
#include "common.hh"


using namespace EXTRAP_MMAP_CONST;

class extrap_mmap {
	public:
		/** number of spots in array (sm4*sn4*so4) */
		int m;
		/** difference between processor corner and ghost corner */
		int i0;
		/** The number of maps at a position. */
		int *n;
		/** analogous to u0 */
		int *n0;
		/** A array of pointers to the maps at a position. */
		ref_map **f;
		/** analogous to u0 */
		ref_map **f0;
		/** An array of arrays of multimaps. */
		ref_map **mec;
		/** An array of backpointers. */
		int **bptr;
		/** The number of entries in each multimap array. */
		int *nmap;
		/** The sizes of the multimap arrays. */
		int *msz;

		extrap_mmap(const geometry &gm);
		~extrap_mmap();

        // Add a reference map variable in extrapolation process, creaing a new value
		void add_mem(int ind,unsigned int oid,double (&xi)[3],int type);
        // Add a reference map any other time based on another variable
		void add_mem(int ind, ref_map * rmp_ptr);

        template<int flags>
		void add_mem(int ind,double *(&b));

		inline void add(int ind, ref_map *rmp_ptr){
			add_mem(i0+ind,rmp_ptr);
		}
		bool rm_mem(int ind,int oid);
		inline bool rm(int ind,int oid) {return rm_mem(i0+ind,oid);};
		void expand(int i);
		void get_obj_flag(int ind, std::vector<int> &existing_obj_id);
		void copy(extrap_mmap &another){
			another.reset();
			for(int i=1;i<max_maps;i++){
				for(int j=0;j<nmap[i]*i;j++){
					int bptr_ind = int(j/i);
					//another.add_mem(bptr[i][bptr_ind], mec[i][j].c, mec[i][j].x, 0);
					another.add_mem(bptr[i][bptr_ind], mec[i]+j);
				}
			}
		}
		template<int flags>
		int pack(int ind,double *(&b),int lint);
		template<int flags>
		void unpack(int ind,double *(&b));


		/** reset multimap structure for next round of extrapolations */
		inline void reset() {
			for (int i = 1; i < max_maps; i++) nmap[i] = 0;
			for (int i = 0; i < m; i++) n[i]=0;
		}

		// TODO make this consistent with unpack ... where do faces fit in?
		/** pack a set of extrapolations into a communication buffer */
		inline int pack_extrap(int ind, unsigned int c, double *(&cb), int type) {

			// pack number of extrapolations at this point
			int n_here = 0;
			for (int j = 0; j < n[ind]; j++)
				if (f[ind][j].c == c) n_here++;

			// only packing those in this layer
			*(cb++) = 0.5 + n_here;

			// now, for each one...
			for (int j = 0; j < n[ind]; j++) {
				if (f[ind][j].c == c) {

					// store ref map
					if(type==0) f[ind][j].pack_x(cb);
					else if(type==1) f[ind][j].pack_xpred(cb);
					else if(type==2) f[ind][j].pack_fxs(cb);
					else p_fatal_error("extrap_mmap::pack_extrap(): Undefined type.",MG3D_INTERNAL_ERROR);
				}
			}
			return n_here;
		}

		/** pack a set of extrapolations into a communication buffer */
		inline int pack_extrap(int ind,double *(&cb),int type) {

			// pack number of extrapolations at this point
			*(cb++) = 0.5 + n[ind];

			// now, for each one...
			for (int j = 0; j < n[ind]; j++) {

				// store ref map
					if(type==0) f[ind][j].pack_x(cb);
					else if(type==1) f[ind][j].pack_xpred(cb);
					else if(type==2) f[ind][j].pack_fxs(cb);
					else  p_fatal_error("extrap_mmap::pack_extrap(): Undefined type.",MG3D_INTERNAL_ERROR);
				/*switch (type) {
					case 0: f[ind][j].pack_x(cb); break;
					case 1: f[ind][j].pack_xpred(cb); break;
					case 2: f[ind][j].pack_fxs(cb); break;
				}
				*/
			}

			return n[ind];
		}

		/** unpack a set of extrapolations from a communication buffer */
		inline void unpack_extrap(int ind,double *(&cb),int type) {

			// unpack number of extraps at this point
			int n_here = (int) *(cb++);

			// ref map components
			ref_map rm;

			// now, for each extrap here,
			for (int i = 0; i < n_here; i++) {
				if(type==0){
						rm.unpack_x(cb);
						add_mem(ind,rm.c,rm.x,type);
				}
				else if(type==1){
						rm.unpack_xpred(cb);
						add_mem(ind,rm.c,rm.xpred,type);
				}
				else if(type==2){
						rm.unpack_fxs(cb);
						for (int j = 0; j < n[ind]; j++) {
							if (f[ind][j].c == rm.c) {
								for (int g = 0; g < 6; g++)
									for (int v = 0; v < 3; v++)
										f[ind][j].fx[g][v] = rm.fx[g][v];
							}
						}
				}
				else p_fatal_error("extrap_mmap::unpack_extrap(): Undefined type.", MG3D_INTERNAL_ERROR);
			}
		}

		/** Prints a list of maps at a real gridpoint.
		 * \param[in] i the gridpoint to consider. */
		inline void print_maps(int i) {
			for(int j=0;j<n[i];j++)
				printf(" (layer %d, object %d: %g,%g,%g)",f0[i][j].lid(),
						f0[i][j].oid(),f0[i][j].x[0],f0[i][j].x[1],f0[i][j].x[2]);
			puts("");
		}
};

#endif
