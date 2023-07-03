
//@HEADER
// ***************************************************
//
// HPCG: High Performance Conjugate Gradient Benchmark
//
// Contact:
// Michael A. Heroux ( maherou@sandia.gov)
// Jack Dongarra     (dongarra@eecs.utk.edu)
// Piotr Luszczek    (luszczek@eecs.utk.edu)
//
// ***************************************************
//@HEADER

/*!
 @file GenerateGeometry.cpp

 HPCG routine
 */

#include <cmath>
#include <cstdlib>
#include <cassert>

#include "ComputeOptimalShapeXYZ.hpp"
#include "GenerateGeometry.hpp"

#ifdef HPCG_DEBUG
#include <fstream>
#include "hpcg.hpp"
using std::endl;

#endif

Geometry::Geometry(const int size_a, const int rank_a, const int numThreads_a,
  const int pz_a, const local_int_t zl, const local_int_t zu,
  const local_int_t nx_a, const local_int_t ny_a, const local_int_t nz_a,
  const int npx_a, const int npy_a, const int npz_a)
    : size{size_a}, rank{rank_a}, numThreads{numThreads_a}, nx{nx_a}, ny{ny_a}, nz{nz_a},
    npx{npx_a}, npy{npy_a}, npz{npz_a}, pz{pz_a}, npartz{0}
{

  if (npx * npy * npz <= 0 || npx * npy * npz > size)
    ComputeOptimalShapeXYZ( size, npx, npy, npz );

  if (pz==0) { // No variation in nz sizes
    npartz = 1;
    partz_ids.resize(1);
    partz_nz.resize(1);
    partz_ids[0] = npz;
    partz_nz[0] = nz;
  }
  else {
    npartz = 2;
    partz_ids.resize(2);
    partz_ids[0] = pz;
    partz_ids[1] = npz;
    partz_nz.resize(2);
    partz_nz[0] = zl;
    partz_nz[1] = zu;
  }
//  partz_ids[npartz-1] = npz; // The last element of this array is always npz
  int ipartz_ids = 0;
  for (int i=0; i< npartz; ++i) {
    assert(ipartz_ids<partz_ids[i]);  // Make sure that z partitioning is consistent with computed npz value
    ipartz_ids = partz_ids[i];
  }

  // Now compute this process's indices in the 3D cube
  ipz = rank/(npx*npy);
  ipy = (rank-ipz*npx*npy)/npx;
  ipx = rank%npx;

#ifdef HPCG_DEBUG
  if (rank==0)
    HPCG_fout   << "size = "<< size << endl
        << "nx  = " << nx << endl
        << "ny  = " << ny << endl
        << "nz  = " << nz << endl
        << "npx = " << npx << endl
        << "npy = " << npy << endl
        << "npz = " << npz << endl;

  HPCG_fout    << "For rank = " << rank << endl
      << "ipx = " << ipx << endl
      << "ipy = " << ipy << endl
      << "ipz = " << ipz << endl;

  assert(size>=npx*npy*npz);
#endif

  // These values should be defined to take into account changes in nx, ny, nz values
  // due to variable local grid sizes
  gnx = npx*nx;
  gny = npy*ny;
  //global_int_t gnz = npz*nz;
  // We now permit varying values for nz for any nx-by-ny plane of MPI processes.
  // npartz is the number of different groups of nx-by-ny groups of processes.
  // partz_ids is an array of length npartz where each value indicates the z process of the last process in the ith nx-by-ny group.
  // partz_nz is an array of length npartz containing the value of nz for the ith group.

  //        With no variation, npartz = 1, partz_ids[0] = npz, partz_nz[0] = nz

  gnz = 0;
  ipartz_ids = 0;

  for (int i=0; i< npartz; ++i) {
    ipartz_ids = partz_ids[i] - ipartz_ids;
    gnz += partz_nz[i]*ipartz_ids;
  }
  //global_int_t giz0 = ipz*nz;
  giz0 = 0;
  ipartz_ids = 0;
  for (int i=0; i< npartz; ++i) {
    int ipart_nz = partz_nz[i];
    if (ipz < partz_ids[i]) {
      giz0 += (ipz-ipartz_ids)*ipart_nz;
      break;
    } else {
      ipartz_ids = partz_ids[i];
      giz0 += ipartz_ids*ipart_nz;
    }

  }
  gix0 = ipx*nx;
  giy0 = ipy*ny;
}
