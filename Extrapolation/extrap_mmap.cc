#include "extrap_mmap.hh"

extrap_mmap::extrap_mmap(const geometry &gm) :
	m((gm.sm+4)*(gm.sn+4)*(gm.so+4)),i0(2*(1 + (gm.sm+4)*(1 + (gm.sn+4)))),
	n(new int[m]),n0(n+i0),f(new ref_map*[m]),f0(f+i0),
	mec(new ref_map*[max_maps]),bptr(new int*[max_maps]),
	nmap(new int[max_maps]),msz(new int[max_maps]) {
	// Allocate memory for multi-map structure
	for(int i=1;i<max_maps;i++) {
		msz[i]=init_map_mem;
		mec[i]=new ref_map[i*msz[i]];
		bptr[i]=new int[msz[i]];
	}

	reset();
}

extrap_mmap::~extrap_mmap() {
	for(int i=max_maps-1;i>0;i--) {
		delete [] bptr[i];
		delete [] mec[i];
	}
	delete [] msz;
	delete [] nmap;
    delete [] bptr;
	delete [] mec;
	delete [] f;
	delete [] n;
}

template<int flags>
int extrap_mmap::pack(int ind,double *(&b),int lint) {
	double &n_here = *(b++);
	n_here = 0.5;
	for (int j = 0; j < n[ind]; j++) {
		if (lint<1||f[ind][j].lid()==lint) {
			n_here++;
			f[ind][j].pack<flags>(b);
		}
	}
	return n_here;
}

template<int flags>
void extrap_mmap::unpack(int ind,double *(&b)) {
	int n_here = *(b++);
	for (int j = 0; j < n_here; j++) {
		bool found = false;
		for (int i = 0; i < n[ind] && (!found); i++) {
			if (f[ind][i].c == *b) {
				f[ind][i].unpack<flags>(b);
				found=true;
			}
		}
		if (!found) add_mem<flags>(ind,b);
	}
}


template<int flags>
void extrap_mmap::add_mem(int ind,double *(&b)) {

	int q=++n[ind];

	// XXX this hasn't been tested
	if(nmap[q]==msz[q]) {
		expand(q);
	}

	// We are going to store into slot mec[n[i]][n[i]*nmap[n[i]]]
	ref_map *fp=mec[q]+q*nmap[q],*fp2=fp;

	// If there were more than zero maps at i, then we need to copy
	// them into the first n[i]-1 entries of the slot. Set the
	// backpointer of the new entry to be i.
	for(int j=0;j<q-1;j++) {
		*(fp++)=f[ind][j];
	}

	// read in the information, using communication flags
	fp->unpack<flags>(b);

	// also set up backpointer
	bptr[q][nmap[q]++]=ind;

	// Rejigger mec[q-1] to be continuous. Copy the last entry of
	// mec[q-1] over the one that was moved into mec[q]
	if(q>1) {
		nmap[q-1]--;
		for(int j=0;j<q-1;j++) f[ind][j]=mec[q-1][(q-1)*nmap[q-1]+j];
		int l=static_cast<int>(f[ind]-mec[q-1])/(q-1);
		bptr[q-1][l]=bptr[q-1][nmap[q-1]];

		// Update the pointer on the moved slot
		int k=bptr[q-1][nmap[q-1]];
		if(k!=ind) f[k]=f[ind];
	}

	// Finally, set the pointer at the current gridpoint to the new
	// slot in mec[q] that was created
	f[ind]=fp2;
}

/** Adds a map to the structure.
 */
void extrap_mmap::add_mem(int ind, ref_map * rmp_ptr){

	int q=++n[ind];
	// XXX this hasn't been tested
	if(nmap[q]==msz[q]) {
		expand(q);
	}

	// We are going to store into slot mec[n[i]][n[i]*nmap[n[i]]]
	ref_map *fp=mec[q]+q*nmap[q],*fp2=fp;

	// If there were more than zero maps at i, then we need to copy
	// them into the first n[i]-1 entries of the slot. Set the
	// backpointer of the new entry to be i.
	for(int j=0;j<q-1;j++) {
		*(fp++)=f[ind][j];
	}

	// introduce ref map components
    *fp = *(rmp_ptr);

	// also set up backpointer
	bptr[q][nmap[q]++]=ind;

	// Rejigger mec[q-1] to be continuous. Copy the last entry of
	// mec[q-1] over the one that was moved into mec[q]
	if(q>1) {
		nmap[q-1]--;
		for(int j=0;j<q-1;j++) f[ind][j]=mec[q-1][(q-1)*nmap[q-1]+j];
		int l=static_cast<int>(f[ind]-mec[q-1])/(q-1);
		bptr[q-1][l]=bptr[q-1][nmap[q-1]];

		// Update the pointer on the moved slot
		int k=bptr[q-1][nmap[q-1]];
		if(k!=ind) f[k]=f[ind];
	}

	// Finally, set the pointer at the current gridpoint to the new
	// slot in mec[q] that was created
	f[ind]=fp2;
}

/** Adds a map to the structure.
 * \param[in] ind the gridpoint to add to.
 * \param[in] o_id the object ID of the map.
 * \param[in] xi_x the x-coordinate of the reference map to store
 * \param[in] xi_y the y-coordinate of the reference map to store
 * \param[in] xi_z the z-coordinate of the reference map to store
 */
void extrap_mmap::add_mem(int ind,unsigned int c,double (&xi)[3],int type) {

	int q=++n[ind];
	// XXX this hasn't been tested
	if(nmap[q]==msz[q]) {
		expand(q);
	}

	// We are going to store into slot mec[n[i]][n[i]*nmap[n[i]]]
	ref_map *fp=mec[q]+q*nmap[q],*fp2=fp;

	// If there were more than zero maps at i, then we need to copy
	// them into the first n[i]-1 entries of the slot. Set the
	// backpointer of the new entry to be i.
	for(int j=0;j<q-1;j++) {
		*(fp++)=f[ind][j];
	}

	// introduce ref map components
	if (type == 0 || type==1) {
		for (int k = 0; k < 3; k++) fp->xpred[k] = fp->x[k] = xi[k];
	} else {
		p_fatal_error("extrap_mmap::add_mem(): Undefined type!\n", PGMG_INTERNAL_ERROR);
	}
	fp->c = c;

	// also set up backpointer
	bptr[q][nmap[q]++]=ind;

	// Rejigger mec[q-1] to be continuous. Copy the last entry of
	// mec[q-1] over the one that was moved into mec[q]
	if(q>1) {
		nmap[q-1]--;
		for(int j=0;j<q-1;j++) f[ind][j]=mec[q-1][(q-1)*nmap[q-1]+j];
		int l=static_cast<int>(f[ind]-mec[q-1])/(q-1);
		bptr[q-1][l]=bptr[q-1][nmap[q-1]];

		// Update the pointer on the moved slot
		int k=bptr[q-1][nmap[q-1]];
		if(k!=ind) f[k]=f[ind];
	}

	// Finally, set the pointer at the current gridpoint to the new
	// slot in mec[q] that was created
	f[ind]=fp2;
}

/** Removes a map from the structure.
 * \param[in] ind the gridpoint to remove from.
 * \param[in] oid the object ID of the map (just object, no layer!)
 * \param[out] whether such an object was there to be removed
 */
bool extrap_mmap::rm_mem(int ind, int oid) {

	// can't remove if there was nothing there.
	if (n[ind] == 0) return false;

	int j0 = -1;
	for (int j = 0; j < n[ind]; j++)
		if (f[ind][j].oid() == oid && j0 == -1) j0 = j;

	// can't remove if this object wasn't there
	if (j0 == -1) return false;

	// now we know it's there, so take it out
	int q = --n[ind];

	// need to copy over if there were other extraps there
	ref_map *fp = NULL,*fp2 = NULL;
	if (q > 0) {

		//TODO is this working?
		if(nmap[q]==msz[q]) expand(q);

		// We are going to store into slot mec[n[i]][n[i]*nmap[n[i]]]
		fp2 = fp =mec[q]+q*nmap[q];

		// If there were more than one maps at i, then we need to copy
		// the remaining maps into the n[i]-1 entries on the lower level. Set the
		// backpointer of the new entry to be i.
		for(int j=0;j<q+1;j++) if (j != j0) *(fp++)=f[ind][j];

		// also set up backpointer
		bptr[q][nmap[q]++]=ind;
	}

	// Rejigger mec[q+1] to be continuous. Copy the last entry of
	// mec[q+1] over the one that was moved into mec[q]
	nmap[q+1]--;
	for(int j=0;j<q+1;j++) f[ind][j]=mec[q+1][(q+1)*nmap[q+1]+j];
	int l=static_cast<int>(f[ind]-mec[q+1])/(q+1);
	bptr[q+1][l]=bptr[q+1][nmap[q+1]];

	// Update the pointer on the moved slot
	int k=bptr[q+1][nmap[q+1]];
	if(k!=ind) f[k]=f[ind];

	// Finally, set the pointer at the current gridpoint to the new
	// slot in mec[q] that was created
	if (q > 0) f[ind]=fp2;

	return true;
}
void extrap_mmap::expand(int q) {

	// allocate new array of ref maps and backpointers
	ref_map *nr = new ref_map[2*q*msz[q]];
	int *nb = new int[2*msz[q]];

	// copy elements over
	for (int i = 0; i < msz[q]; i++) {
		for (int j = 0; j < q; j++) nr[i*q+j] = mec[q][i*q+j];
		nb[i] = bptr[q][i];
	}

	// delete old arrays, assign new arrays
	delete[] mec[q];
	delete[] bptr[q];
	mec[q] = nr;
	bptr[q] = nb;

	// update where f is pointing
	for (int i = 0; i < msz[q]; i++)
		f[bptr[q][i]] = mec[q]+q*i;

	// note size has doubled
	msz[q] *= 2;
}

void extrap_mmap::get_obj_flag(int ind, std::vector<int> &existing_obj_id){
	for (int i =0; i< n0[ind]; i++) {
		int cur_flag = f0[ind][i].c;
		existing_obj_id.push_back(cur_flag);
	}
}

// object id and x
template int extrap_mmap::pack<128|512>(int ind,double *(&b),int lint);
template void extrap_mmap::unpack<128|512>(int ind,double *(&b));
// object id, x, x-pred
template int extrap_mmap::pack<128|256|512>(int ind,double *(&b),int lint);
template void extrap_mmap::unpack<128|256|512>(int ind,double *(&b));
// object id and face x
template int extrap_mmap::pack<128|1024>(int ind,double *(&b),int lint);
template void extrap_mmap::unpack<128|1024>(int ind,double *(&b));

template void extrap_mmap::add_mem<128|512>(int ind,double *(&b));
template void extrap_mmap::add_mem<128|256|512>(int ind,double *(&b));
template void extrap_mmap::add_mem<128|1024>(int ind,double *(&b));
