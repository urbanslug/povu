#ifndef MEZA_MATRIX_POOL_HPP
#define MEZA_MATRIX_POOL_HPP

#include <cstddef>

#include "meza/pool/hap_comp.hpp"
#include "meza/pool/joint.hpp"
#include "meza/pool/mem.hpp"
#include "meza/pool/pool_ops.hpp"
#include "meza/pool/split.hpp"

#if MEZA_USE_CUDA
#include "meza/ops/ops.cuh"
#include "meza/pool/hap_comp.cuh"
#include "meza/pool/split_cuda.cuh"
#endif

namespace meza::pool
{
using namespace meza::pool::hap_comp; // for haps_comp_set
using namespace meza::pool::joint;    // for joint_pool

struct pool_mem_split {
	/*
	 * u8
	 * --------
	 * 1 byte per u8 value
	 * 1024*1024 = 1,048,576 u8 values are
	 *
	 * u32
	 * --------
	 * 4 bytes per u32 value
	 * (1024*1024) / 4 = 262,144
	 * 262144 u32 values are ~1M
	 * 1024 values of u32 are 1M
	 * 2,621,440 u32 values are ~10M
	 *
	 *
	 * hap comp elements
	 * 20ull * 1024 * 1024 * 1024 = 20G elements
	 * 10ull * 1024 * 1024 = 10M elements
	 *
	 */

	// std::size_t split_pool = 10ull * 1024 * 1024 * 1024;	  // 10G
	// std::size_t hap_comp = 50ull * 1024 * 1024 * 1024;	  // 50 G
	// std::size_t joint_pool_size = 20ull * 1024 * 1024 * 1024; // 20 G

	std::size_t split_pool = 0;	 // 5 G
	std::size_t hap_comp = 0;	 // 5 G
	std::size_t joint_pool_size = 0; // 1 G

	pool_mem_split() = delete;

	[[nodiscard]] static inline std::double_t to_g(unsigned long long bytes)
	{
		return static_cast<std::double_t>(bytes) / (1024 * 1024 * 1024);
	}

	void print_split()
	{
		log_info("Memory pool split (in GB): Split pool: %.2f G, "
			 "Hap comp: %.2f G, Joint pool: %.2f G",
			 to_g(split_pool), to_g(hap_comp),
			 to_g(joint_pool_size));
	}

	static pool_mem_split init()
	{
		unsigned long long total = meza::mem::get_total_ram();
		unsigned long long available = meza::mem::get_available_ram();

		unsigned long long min_g = 2ull;
		unsigned long long min = min_g * 1024 * 1024 * 1024; // 2 G

		if (total < min) {
			log_fatal("Total RAM %.2fG is less than the minimum "
				  "required: %uG",
				  to_g(total), min_g);
			std::exit(EXIT_FAILURE);
		}

		// pick the lesser of
		//  - 10% of total RAM
		//  - 90% of available RAM
		unsigned long long target =
			std::min(total / 10, available * 9 / 10);

		// split target into 3 pools with a 4:4:2 ratio
		unsigned long long split_pool = target * 4 / 10;
		unsigned long long hap_comp = target * 4 / 10;
		unsigned long long joint_pool_size = target * 2 / 10;

		// log_info("Total %2f, Available %.2f", to_g(total),
		//	 to_g(available));
		// log_info("Allocating. %.2f  %.2f  %.2f ", to_g(split_pool),
		//	 to_g(hap_comp), to_g(joint_pool_size));

		// std::cerr << split_pool << " " << hap_comp << " "
		//	  << joint_pool_size << "\n";

		return pool_mem_split{split_pool, hap_comp, joint_pool_size};
	}

private:
	pool_mem_split(unsigned long long sp, unsigned long long hp,
		       unsigned long long jp)
	    : split_pool(sp), hap_comp(hp), joint_pool_size(jp)
	{
		if (split_pool + hap_comp + joint_pool_size >
		    meza::mem::get_total_ram()) {
			log_fatal("Total pool size exceeds total RAM. "
				  "Split pool: %.2f G, Hap comp: %.2f G, "
				  "Joint pool: %.2f G, Total RAM: %.2f G",
				  to_g(split_pool), to_g(hap_comp),
				  to_g(joint_pool_size),
				  to_g(meza::mem::get_total_ram()));
			std::exit(EXIT_FAILURE);
		}
	}
};

/**
 * T for the matrix pool (reference, filter, xor),
 * S for the joint pool (depth matrix)
 */
template <typename T, typename S>
struct pool {
public:
	void hap_compare(const meza::pool::ov_mat_t &filter_mat,
			 qt::u32 pool_offset)
	{
		// meza::pool::hap_comp::haps_comp_set cmp_set;
#if MEZA_USE_CUDA
		cmp_mat_cuda.base_mut().set_filter(&filter_mat, pool_offset);
		meza::pool_ops::handle_set(mat_pool_cuda, cmp_mat_cuda);
#else
		cmp_mat_cpu.set_filter(&filter_mat, pool_offset);
		meza::pool_ops::handle_set(mat_pool_cpu, cmp_mat_cpu);
#endif
		// return cmp_set;
	}

	const hap_comp_matrix<S> &get_hap_comp_matrix() const
	{
		return cmp_mat_cpu;
	}

	// [[nodiscard]] std::size_t hap_cmp_count() const
	// {
	//	return cmp_mat_cpu.hap_cmp_count();
	// }

	// [[nodiscard]] inline qt::up_t<qt::u32> comp_rm_idx(std::size_t i)
	// const
	// {
	//	return cmp_mat_cpu.comp_rm_idx(i);
	// }

	void clear_hap_cmp_data()
	{
		cmp_mat_cpu.clear_hap_cmp_data();
	}

	void run_convolutions(qt::u32 pool_j_offset)
	{
#if MEZA_USE_CUDA
		mat_pool_cuda.copy_to_device();
		meza::pool::matrix_pool_cuda<T> &p = this->mat_pool_cuda;

		const u8 *d_ref = p.ref_ptr_mut();
		const u8 *d_filter = p.filter_ptr_mut();
		u8 *d_xor = p.xor_ptr_mut();

		meza::cuda_ops::cuda_mat_xor(d_ref, d_filter, d_xor,
					     pool_j_offset);

		mat_pool_cuda.copy_region_to_host(meza::pool::pool_region::Xor);
#else
		meza::pool::matrix_pool<T> &p = this->mat_pool_cpu;

		const u8 *h_ref = p.ref_start_ptr();
		const u8 *h_filter = p.filter_start_ptr();
		u8 *h_xor = p.xor_start_ptr();

		meza::cpu_ops::cpu_mat_xor(h_ref, h_filter, h_xor,
					   pool_j_offset);
#endif
	}

	[[nodiscard]] bool is_full() const
	{
#if MEZA_USE_CUDA
		return mat_pool_cuda.base().is_full();
#else
		return mat_pool_cpu.is_full();
#endif
	}

	void clear_split_pool()
	{
#if MEZA_USE_CUDA
		mat_pool_cuda.clear();
#else
		mat_pool_cpu.clear();
#endif
	}

	void reset_depth_matrix()
	{
		joint_pool_cpu.clear();
	}

	template <typename U, typename W>
	[[nodiscard]] meza::view::ov_matrix<T, U, W>
	alloc_ov_matrix(qt::u32 I, qt::u32 J, pool_region region)
	{
		return mat_pool_cpu.template alloc_ov_matrix<U, W>(I, J,
								   region);
	}

	meza::pool::joint::full_view<S> alloc_depth_matrix(qt::u32 I, qt::u32 J)
	{
		return joint_pool_cpu.alloc_full(I, J);
	}

	static pool init(u32 H)
	{
		auto pool_split = pool_mem_split::init();
		return pool(pool_split.split_pool, pool_split.hap_comp,
			    pool_split.joint_pool_size, H);
	}

private:
	/* ================= constructor ============== */

	pool(std::size_t sp_sz, std::size_t hc_sz, std::size_t jp_sz, u32 H)
	    : mat_pool_cpu(matrix_pool<T>::create(sp_sz)),
	      cmp_mat_cpu(hap_comp_matrix<S>::create(hc_sz, H)),
	      joint_pool_cpu(joint_pool<S>::init(jp_sz))
#if MEZA_USE_CUDA
	      ,
	      mat_pool_cuda(mat_pool_cpu),
	      cmp_mat_cuda(hap_comp_matrix_cuda<T>{cmp_mat_cpu})
#endif
	{}

	meza::pool::matrix_pool<T> mat_pool_cpu;
	meza::pool::hap_comp::hap_comp_matrix<S> cmp_mat_cpu;
	meza::pool::joint::joint_pool<S> joint_pool_cpu;

#if MEZA_USE_CUDA
	meza::pool::matrix_pool_cuda<T> mat_pool_cuda;
	meza::pool::hap_comp::hap_comp_matrix_cuda<S> cmp_mat_cuda;
#endif
};

}; // namespace meza::pool

#endif // MEZA_MATRIX_POOL_HPP
