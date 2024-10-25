#ifndef CORE_HPP
#define CORE_HPP

#include <cstdint>
#include <string>
#include <unistd.h>
#include <vector>
#include <utility>
#include <iostream>
#include <filesystem>

#include "../common/utils.hpp"

namespace core {
namespace pu = povu::utils;


/*
 * =========
 * App types
 * =========
 */

// rename to task_e
enum class task_t {
  call,        // call variants
  deconstruct, // deconstruct a graph
  info,        // print graph information
  unset        // unset
};

// rename to input_format_e
enum class input_format_t {
  file_path,
  params,  // CLI params
  unset // unset
};
std::ostream& operator<<(std::ostream& os, const task_t& t);

enum class subgraph_category {
  bubble,
  component,
  unset
};

/**
 * @brief app config
 */
struct config {
  task_t task;
  std::string input_gfa;
  std::filesystem::path forest_dir;
  std::string chrom;
  //std::optional<std::filesystem::path> pvst_path;
  std::filesystem::path output_dir; // output directory for task and deconstruct

  // general
  unsigned char v; // verbosity
  bool print_dot_ { true }; // generate dot format graphs

  unsigned int thread_count_ {1}; // number of threads to use

  // references
  std::string references_txt; // the path to the file containing the reference paths
  input_format_t ref_input_format;
  std::vector<std::string> reference_paths;
  bool undefined_vcf;

  // -------------
  // Contructor(s)
  // -------------

  // constructor(s)
  config()
      : task(task_t::unset),
        chrom(""), // default is empty string
        output_dir("."),                // default is current directory
        v(0),
        thread_count_(1),
        ref_input_format(input_format_t::unset),
        reference_paths(std::vector<std::string>{}),
        undefined_vcf(false)
    {}

  // ---------
  // getter(s)
  // ---------
  std::string get_input_gfa() const { return this->input_gfa; }
  std::filesystem::path get_forest_dir() const { return this->forest_dir; }
  std::filesystem::path get_output_dir() const { return this->output_dir; }
  const std::string& get_chrom() const { return this->chrom; }
  std::vector<std::string> const& get_reference_paths() const { return this->reference_paths; }
  std::vector<std::string>* get_reference_ptr() { return &this->reference_paths; }
  const std::string& get_references_txt() const { return this->references_txt; }
  std::size_t verbosity() const { return this->v; } // can we avoid this being a size_t?
  unsigned int thread_count() const { return this->thread_count_; }
  bool print_dot() const { return this->print_dot_; }
  bool gen_undefined_vcf() const { return this->undefined_vcf; }
  task_t get_task() const { return this->task; }

  // ---------
  // setter(s)
  // ---------
  // as it os from the user not the handlegraph stuff
  void set_chrom(std::string&& s) { this->chrom = s; }
  void set_ref_input_format(input_format_t f) { this->ref_input_format = f; }
  void add_reference_path(std::string s) { this->reference_paths.push_back(s); }
  void set_reference_paths(std::vector<std::string>&& v) { this->reference_paths = std::move(v); }
  void set_reference_txt_path(std::string&& s) { this->references_txt = std::move(s); }
  void set_references_txt(std::string s) { this->references_txt = s; }
  void set_verbosity(unsigned char v) { this->v = v; }
  void set_thread_count(uint8_t t) { this->thread_count_ = t; }
  void set_print_dot(bool b) { this->print_dot_ = b; }
  void set_input_gfa(std::string s) { this->input_gfa = s; }
  void set_forest_dir(std::string s) { this->forest_dir = s; }
  void set_output_dir(std::string s) { this->output_dir = s; }
  void set_task(task_t t) { this->task = t; }
  void set_undefined_vcf(bool b) { this->undefined_vcf = b; }

  // --------
  // other(s)
  // --------
  void dbg_print() {
    std::cerr << "CLI parameters: " << std::endl;
    std::cerr << "\t" << "verbosity: " << this->verbosity() << "\n";
    std::cerr << "\t" << "thread count: " << this->thread_count() << "\n";
    std::cerr << "\t" << "print dot: " << (this->print_dot() ? "yes" : "no") << "\n";
    std::cerr << "\t" << "task: " << this->task << std::endl;
    std::cerr << "\t" << "input gfa: " << this->input_gfa << std::endl;
    std::cerr << "\t" << "forest dir: " << this->forest_dir << std::endl;
    std::cerr << "\t" << "output dir: " << this->output_dir << std::endl;
    std::cerr << "\t" << "chrom: " << this->chrom << std::endl;
    std::cerr << "\t" << "Generate undefined vcf: " << std::boolalpha << this->undefined_vcf << std::endl;
    if (this->ref_input_format == input_format_t::file_path) {
      std::cerr << "\t" << "Reference paths file: " << this->references_txt << std::endl;
    }

    std::cerr << "\t" << "Reference paths (" << this->reference_paths.size() << "): ";
    pu::print_with_comma(std::cerr, this->reference_paths, ',');

    std::cerr << std::endl;
    }
};
} // namespace core
#endif
