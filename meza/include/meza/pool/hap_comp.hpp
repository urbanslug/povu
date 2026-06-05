#ifndef MZ_MATRIX_POOL_HAP_COMP_HPP
#define MZ_MATRIX_POOL_HAP_COMP_HPP

#include <cassert>
#include <cstddef>
#include <cstdio>
#include <cstring>
#include <stdexcept>
#include <string>
#include <vector>

// #include <algorithm>
// #include <chrono>
// #include <future>
// #include <optional>
// #include <set>
// #include <string_view>
// #include <unordered_map>
// #include <unordered_set>

#include "meza/ops/ops.hpp" // for prefix_sum_cpu
#include <log/log.h>
#include <quilt/shim.hpp>  // for qt::shim::contains
#include <quilt/types.hpp> // for qt::u32, qt::u8, qt::op_t

#include "meza/pool/split.hpp"		  // for matrix_pool
#include "meza/pool/split_pool_types.hpp" // for ov_mat_t

namespace meza::pool::hap_comp
{

// -------
// aliases
// -------
using meza::pool::comparison_op;
using meza::pool::matrix_pool;
using qt::u32, qt::u8;

/**
 * haplotype comparison matrix
 *
 *
 * a square matrix
 *
 * stores from diagonal 1 upwards
 */
template <typename T>
struct hap_comp_matrix {
public:
	// ----------------
	// helpers (public)
	// ----------------

	/*  ----- accessors ----- */

	/**
	 * @brief Given k & k-offset, calculate the row & col pair in
	 * the matrix for the haplotype comparison.
	 */
	[[nodiscard]] inline static qt::up_t<u32> comp_hap_pair(u32 k,
								u32 k_offset)
	{
		return {k_offset, k + k_offset}; // {h1, h2}
	}

	[[nodiscard]] inline u32 comp_rm_idx(u32 k, u32 k_off) const
	{
		return k_offset(k) + k_off;
	}

	/**
	 * calculates the number of comparisons in the k-th diagonal of the
	 * haplotype comparison matrix. This is based on the fact that the
	 * matrix is upper triangular and has H_ rows, so the length of the
	 * k-th diagonal is H_ - k.
	 */
	[[nodiscard]] u32 k_len(u32 k) const
	{
		return this->H_ - k;
	}

	/**
	 * Number of elements from diagonal k upwards in the upper triangle of
	 * the matrix.
	 */
	[[nodiscard]] u32 elements(u32 k) const
	{
		u32 H = this->rows();
		if (k >= H)
			throw std::out_of_range(
				"k value out of range for haplotype comparison "
				"matrix. k: " +
				std::to_string(k) + " H: " + std::to_string(H));

		return ((H - k) * (H - k + 1)) / 2;
	}

	/**
	 * calculates the offset in the haplotype comparison matrix for
	 * a given k value. This is based on the number of comparisons
	 * that come before the k-th diagonal in the upper triangle of
	 * the matrix.
	 */
	[[nodiscard]] u32 k_offset(u32 k) const
	{
		return elements(1) - elements(k);
	}

	[[nodiscard]] u32 pool_offset() const
	{
		return pool_offset_;
	}

	/**
	 * number of elements in the upper triangle of the matrix (excluding
	 * diagonal) multiplied by the number of elements per row in the
	 * filter matrix. This is used to determine how much space to
	 * allocate for the xor and sum results in the pool.
	 */
	[[nodiscard]] u32 cols() const
	{
		return J_;
	}

	/**
	 * also the no. of haplotypes (H)
	 */
	[[nodiscard]] u32 rows() const
	{
		return H_;
	}

	[[nodiscard]] const T *get_xor_data() const
	{
		return xor_data_;
	}

	[[nodiscard]] const T *get_sum_data() const
	{
		return sum_data_;
	}

	[[nodiscard]] u32 capacity() const
	{
		return this->capacity_;
	}

	[[nodiscard]] std::size_t hap_cmp_count() const
	{
		return this->exp_comparisons_;
	}

	[[nodiscard]] std::size_t allele_cmp_count() const
	{
		return this->hap_cmp_count() * this->cols();
	}

	[[nodiscard]] const std::vector<qt::u8> &get_matches() const
	{
		return this->matches_;
	}

	[[nodiscard]] const std::vector<qt::u8> &get_mismatches() const
	{
		return this->mismatches_;
	}

	[[nodiscard]] const std::vector<qt::u8> &get_reversals() const
	{
		return this->reversals_;
	}

	/*  ----- modifiers ----- */

	[[nodiscard]] T *get_sum_data_mut()
	{
		return sum_data_;
	}

	[[nodiscard]] T *get_xor_data_mut()
	{
		return xor_data_;
	}

	[[nodiscard]] const std::vector<u32> &get_valid_ks() const
	{
		return valid_ks_;
	}

	[[nodiscard]]
	const std::vector<qt::op_t<u32>> &get_hap_comparable() const
	{
		return this->hap_comparable_;
	}

	void clear()
	{
		u32 N = this->exp_comparisons_;
		std::memset(xor_data_, T{}, N * sizeof(T));
		std::memset(sum_data_, T{}, N * sizeof(T));
	}

	void clear_hap_cmp_data()
	{
		std::size_t N_comparisons = sizeof(u8) * exp_comparisons_;
		std::memset(reversals_.data(), 0, N_comparisons);
		std::memset(matches_.data(), 0, N_comparisons);
		std::memset(mismatches_.data(), 0, N_comparisons);

		this->valid_ks_.clear();
		this->hap_comparable_.clear();
	}

	void comp_comparable_coords()
	{
		const auto &b = filter_->base();

		u32 K = this->rows();

		for (u32 k = 1; k < K; k++) {

			if (b.is_row_blank(k))
				continue;

			this->valid_ks_.emplace_back(k);

			u32 k_len = this->k_len(k);
			for (u32 k_off = 0; k_off < k_len; k_off++) {
				auto [h1, h2] = comp_hap_pair(k, k_off);
				if (b.is_row_blank(h1) || b.is_row_blank(h2))
					continue;

				this->hap_comparable_.emplace_back(k, k_off);
			}
		}
	}

	/**
	 * finds reversals based on the special prefix sum calculated from the
	 * sum results of the comparisons.
	 */
	void find_rx()
	{
		const u32 J = this->cols();

		for (auto [k, k_off] : this->hap_comparable_) {
			u32 rm_idx = comp_rm_idx(k, k_off);

			if (mismatches_[rm_idx] == 0)
				continue;

			u32 running = 0;
			std::size_t start = rm_idx * J;
			std::size_t end = (start + J) - 1;

			for (std::size_t i = start; i <= end; i++) {
				if (sum_data_[i] > 2)
					running += sum_data_[i];

				sum_data_[i] = running;
			}

			if (running > 0)
				reversals_[rm_idx] = 1;
		}
	}

	/**
	 * finds matches and mismatches based on the xor results.
	 */
	void find_mx()
	{
		const u32 J = this->cols();
		auto no_inc = [&](u32 start, u32 end) -> bool
		{
			if (start == 0)
				return xor_data_[end] == 0;

			// start > 0
			return (xor_data_[end] - xor_data_[start - 1]) == 0;
		};

		for (auto [k, k_off] : this->hap_comparable_) {
			// row major idx for the pair of haplotypes
			// (h1, h2) in the upper triangle of the matrix
			u32 rm_idx = comp_rm_idx(k, k_off);

			u32 start = rm_idx * J;
			u32 end = (start + J) - 1;

			T *start_ptr = &xor_data_[start];
			meza::cpu_ops::prefix_sum_cpu(start_ptr, J);

			if (no_inc(start, end))
				matches_[rm_idx] = 1;
			else
				mismatches_[rm_idx] = 1;
		}
	}

	// ------------
	// constructors
	// ------------

	void set_filter(const ov_mat_t *f, std::size_t pool_offset)
	{
		this->filter_ = f;
		this->pool_offset_ = pool_offset;
		this->J_ = filter_->base().cols();
	}

	[[nodiscard]] const ov_mat_t &get_filter() const
	{
		return *this->filter_;
	}

	// delete copy, default, move, assignment constructors
	hap_comp_matrix() = delete;
	hap_comp_matrix(const hap_comp_matrix &) = delete;
	hap_comp_matrix &operator=(const hap_comp_matrix &) = delete;
	hap_comp_matrix(hap_comp_matrix &&) = delete;
	hap_comp_matrix &operator=(hap_comp_matrix &&) = delete;

	static hap_comp_matrix<T> create(std::size_t capacity, u32 H)
	{
		return hap_comp_matrix<T>{capacity, H};
	}

	~hap_comp_matrix()
	{
		delete[] xor_data_;
		delete[] sum_data_;
	}

private:
	/* ================ private data members ====================== */

	// at most ~ 10,000 * 10,000 haps (10M elements)
	static constexpr std::size_t max_comparisons_ = 100 * 1024 * 1024;

	// a reference to the filter matrix
	const ov_mat_t *filter_ = nullptr;

	// The offset in the pool where the haplotype comparison matrix starts.
	// This is used to calculate the correct offsets for
	// accessing the data in the pool when performing comparisons.
	std::size_t pool_offset_ = 0;

	const u32 H_ = 0; // set once in constructor
	u32 J_ = 0;

	std::size_t capacity_; // in bytes

	// number of comparisons in the upper triangle of the matrix
	std::size_t exp_comparisons_ = 0;

	// pointers to the data in the pool for the xor and sum results
	// these are used to store the results of the comparisons for
	// each pair of haplotypes. The data is stored in a flattened
	// format, where the comparisons for each pair of haplotypes are
	// stored contiguously in memory. The offsets for accessing the
	// correct data for each pair of haplotypes are calculated based
	// on the k value and the number of comparisons that come before
	// it in the upper triangle of the matrix.
	T *xor_data_;
	T *sum_data_;

	std::vector<qt::u8> matches_;
	std::vector<qt::u8> mismatches_;
	std::vector<qt::u8> reversals_;

	// valid k values for comparison (i.e., those that are not blank in M_f)
	std::vector<u32> valid_ks_;
	std::vector<qt::op_t<u32>> hap_comparable_;

	/* ================= private helper functions ================== */

	// ---------------------
	// constructor (private)
	// ---------------------

	hap_comp_matrix(std::size_t capacity, u32 H)
	    : H_(H), capacity_(capacity), xor_data_(new T[capacity]),
	      sum_data_(new T[capacity]), matches_(max_comparisons_, 0),
	      mismatches_(max_comparisons_, 0), reversals_(max_comparisons_, 0)
	{
		std::size_t N = this->max_comparisons_;
		this->exp_comparisons_ = this->elements(1);

		if (max_comparisons_ < exp_comparisons_) {
			std::string err = qs::format(
				"H is too large. Expected {}, Actual {}",
				exp_comparisons_, max_comparisons_);
			throw std::out_of_range(err);
		}

		valid_ks_.reserve(H_);
		hap_comparable_.reserve(N);
	}
};

} // namespace meza::pool::hap_comp
#endif // MZ_MATRIX_POOL_HAP_COMP_HPP
