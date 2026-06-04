#ifndef MEZA_MATRIX_POOL_OPS_HPP
#define MEZA_MATRIX_POOL_OPS_HPP

#include <algorithm>
#include <execution>
#include <vector>

// #include "log/log.h"
// #include "meza/ops/ops.hpp"
// #include <chrono>
// #include <cstddef>
// #include <future>
// #include <numeric>

#include <quilt/types.hpp>

#include "meza/pool/hap_comp.hpp"
#include "meza/pool/split.hpp"
#include "meza/pool/split_pool_types.hpp"

#if MEZA_USE_CUDA
#include "meza/ops/ops.cuh"
#include "meza/pool/hap_comp.cuh"
#include "meza/pool/split_cuda.cuh"
#endif

namespace meza::pool_ops
{
using namespace meza::pool::hap_comp; // for haps_comp_set
using qt::u32, qt::u8;

void haps_xor_cpu(const meza::pool::matrix_pool<u8> &p,
		  meza::pool::hap_comp::hap_comp_matrix<u32> &cmp_mat, u32 len,
		  u32 col_shift, u32 res_shift);

void haps_sum_cpu(const meza::pool::matrix_pool<u8> &p,
		  meza::pool::hap_comp::hap_comp_matrix<u32> &cmp_mat, u32 len,
		  u32 col_shift, u32 res_shift);

#if MEZA_USE_CUDA
void haps_xor_cuda(const meza::pool::matrix_pool_cuda<u8> &p,
		   meza::pool::hap_comp::hap_comp_matrix_cuda<u32> &cmp_mat,
		   u32 len, u32 col_shift, u32 res_shift,
		   cudaStream_t stream = 0);

void haps_sum_cuda(const meza::pool::matrix_pool_cuda<u8> &p,
		   meza::pool::hap_comp::hap_comp_matrix_cuda<u32> &cmp_mat,
		   u32 len, u32 col_shift, u32 res_shift,
		   cudaStream_t stream = 0);
#endif

#if MEZA_USE_CUDA
template <typename T>
meza::pool::hap_comp::hap_comp_matrix<T> &
base_mat(meza::pool::hap_comp::hap_comp_matrix_cuda<T> &cmp_mat)
{
	return cmp_mat.base_mut();
}
#endif

template <typename T>
meza::pool::hap_comp::hap_comp_matrix<T> &
base_mat(meza::pool::hap_comp::hap_comp_matrix<T> &cmp_mat)
{
	return cmp_mat;
}

/**
 *  T = meza::pool::matrix_pool_cuda<u8>, or
 *  T = meza::pool::matrix_pool<u8>
 *
 *  U = meza::pool::hap_comp::hap_comp_matrix_cuda<u8>, or
 *  U = meza::pool::hap_comp::hap_comp_matrix<u8>
 */
template <typename T, typename U>
void run_in_haps_par(const T &p, U &f, meza::pool::comparison_op op)
{
	meza::pool::hap_comp::hap_comp_matrix<qt::u32> &cmp_mat = base_mat(f);
	const auto &b = f.get_filter().base();

	const u32 J = cmp_mat.cols();

	// 1. Create a vector of your 'k' indices to iterate over
	const std::vector<u32> &ks = cmp_mat.get_valid_ks();

	// 2. Run the loop in parallel
	std::for_each(std::execution::par, ks.begin(), ks.end(),
		      [&](u32 k)
		      {
			      u32 col_shift = k * J;
			      u32 xor_shift = cmp_mat.k_offset(k) * J;
			      u32 len = cmp_mat.k_len(k) * J;

			      if (op == comparison_op::bitwise_xor)
				      haps_xor_cpu(p, cmp_mat, len, col_shift,
						   xor_shift);
			      else if (op == comparison_op::sum)
				      haps_sum_cpu(p, cmp_mat, len, col_shift,
						   xor_shift);
		      });
}

/**
 *  T = meza::pool::matrix_pool_cuda<u8>, or
 *  T = meza::pool::matrix_pool<u8>
 *
 *  U = meza::pool::hap_comp::hap_comp_matrix_cuda<u8>, or
 *  U = meza::pool::hap_comp::hap_comp_matrix<u8>
 *
 * p is the pool containing the reference and filter matrices
 * f is the hap_comp_matrix containing the results of the comparisons
 * (xor and sum)
 *
 * Returns a haps_comp_set containing the results of the comparisons, including
 * the set of reversals and the matches/mismatches for each pair of haplotypes.
 */
template <typename T, typename U>
void handle_set(T &p, U &f)
{
	meza::pool::hap_comp::hap_comp_matrix<u32> &cmp_mat = base_mat(f);

	cmp_mat.comp_comparable_coords();

	run_in_haps_par(p, cmp_mat, comparison_op::bitwise_xor);
	run_in_haps_par(p, cmp_mat, comparison_op::sum);

	base_mat(f).find_mx();
	base_mat(f).find_rx();

	// base_mat(f).clear();
}

} // namespace meza::pool_ops

#endif // MEZA_MATRIX_POOL_OPS_HPP
