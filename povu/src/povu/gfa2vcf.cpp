#include "povu/gfa2vcf.hpp"

#include <cstdlib>    // for exit, mkdtemp, EXIT_FAILURE, size_t
#include <filesystem> // for remove_all, path
#include <string>     // for basic_string, char_traits, opera...

#include <log/log.h>	      // for log_info, log_fatal
#include <povu/call.hpp>      // for do_call
#include <povu/decompose.hpp> // for do_decompose
#include <quilt/shim.hpp>     // for format

namespace povu::subcommands::gfa2vcf
{

void do_gfa2vcf(const core::config &app_config)
{
	std::size_t ll = app_config.verbosity();

	// Create a temporary directory for the forest files
	std::string temp_template = "/tmp/povu_gfa2vcf_XXXXXX";
	char *temp_dir = mkdtemp(temp_template.data());
	if (temp_dir == nullptr) {
		log_fatal("Error: Could not create temporary directory");
		exit(EXIT_FAILURE);
	}
	std::string temp_dir_str(temp_dir);

	if (ll > 0)
		log_info("Using temporary directory: %s", temp_dir_str.c_str());

	// ------------------------------------------------------
	// Step 1: Run decompose to generate forest of PVST files
	// ------------------------------------------------------
	if (ll > 0)
		log_info("Step 1: Decomposing graph");

	// Create a config for decompose with the temp directory
	core::config decompose_config = app_config;
	decompose_config.set_task(core::task_e::decompose);
	decompose_config.set_output_dir(temp_dir_str);

	decompose::do_decompose(decompose_config); // Run decompose

	// --------------------------------
	// Step 2: Run call to generate VCF
	// --------------------------------
	if (ll > 0)
		log_info("Step 2: Calling variants");

	// Create a config for call with stdout output & the temp forest dir
	core::config call_config = app_config;
	call_config.set_task(core::task_e::call);
	call_config.set_forest_dir(temp_dir_str);

	// Run call (it will handle everything including stdout output)
	call::do_call(call_config);

	fs::remove_all(temp_dir_str); // Clean up the temporary directory

	return;
}

} // namespace povu::subcommands::gfa2vcf
