/*
 * decomp_lu.h
 *
 *  Created on: Jul 5, 2009
 *      Author: matheus
 */

#ifndef DECOMP_LU_H_
#define DECOMP_LU_H_

void decomp_lu(float* L[], float* U[], unsigned n);
void solve_lu(float solution[], float* L[], float* U[], float C[], unsigned n);

#endif /* DECOMP_LU_H_ */
