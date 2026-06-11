// #include "log/log.h"
// #include <cstddef>
#include <iomanip>
#include <iostream>

// Detect the Operating System
#if defined(_WIN32) || defined(_WIN64)
#define OS_WINDOWS
#include <windows.h>
#elif defined(__APPLE__) && defined(__MACH__)
#define OS_MACOS
#include <mach/mach.h>
#include <sys/sysctl.h>
#include <sys/types.h>
#elif defined(__linux__) || defined(__linux) || defined(__unix__)
#define OS_LINUX
#include <sys/sysinfo.h>
#include <unistd.h>
#else
#define OS_UNKNOWN
#endif

namespace meza::mem
{
// Returns the total physical RAM in bytes
unsigned long long get_total_ram()
{
#if defined(OS_WINDOWS)
	MEMORYSTATUSEX status;
	status.dwLength = sizeof(status);
	if (GlobalMemoryStatusEx(&status)) {
		return status.ullTotalPhys;
	}
#elif defined(OS_MACOS)
	int mib[2] = {CTL_HW, HW_MEMSIZE};
	int64_t physical_memory = 0;
	size_t length = sizeof(physical_memory);
	if (sysctl(mib, 2, &physical_memory, &length, NULL, 0) == 0) {
		// return static_cast<std::size_t>(physical_memory);
		return (unsigned long long)physical_memory;
	}
#elif defined(OS_LINUX)
	long pages = sysconf(_SC_PHYS_PAGES);
	long page_size = sysconf(_SC_PAGE_SIZE);
	if (pages > 0 && page_size > 0) {
		return (unsigned long long)pages * page_size;
	}
#endif
	return 0; // Unsupported or failed
}

// Returns the currently available (free) physical RAM in bytes
unsigned long long get_available_ram()
{
#if defined(OS_WINDOWS)
	MEMORYSTATUSEX status;
	status.dwLength = sizeof(status);
	if (GlobalMemoryStatusEx(&status)) {
		return status.ullAvailPhys;
	}
#elif defined(OS_MACOS)
	// On macOS, we query the Host VM statistics to get free pages
	vm_size_t page_size;
	mach_port_t mach_port = mach_host_self();
	mach_msg_type_number_t count = HOST_VM_INFO_COUNT;
	vm_statistics_data_t vm_stats;

	if (host_page_size(mach_port, &page_size) == KERN_SUCCESS &&
	    host_statistics(mach_port, HOST_VM_INFO, (host_info_t)&vm_stats,
			    &count) == KERN_SUCCESS) {
		return (unsigned long long)vm_stats.free_count * page_size;
	}
#elif defined(OS_LINUX)
	// On Linux, sysinfo is the standard way to grab freeram
	struct sysinfo info;
	if (sysinfo(&info) == 0) {
		return (unsigned long long)info.freeram * info.mem_unit;
	}
	// Fallback if sysinfo fails
	long pages = sysconf(_SC_AVPHYS_PAGES);
	long page_size = sysconf(_SC_PAGE_SIZE);
	if (pages > 0 && page_size > 0) {
		return (unsigned long long)pages * page_size;
	}
#endif
	return 0; // Unsupported or failed
}

int print_mem_info()
{
	unsigned long long total = get_total_ram();
	unsigned long long available = get_available_ram();

	const double bytes_to_gb = 1024.0 * 1024.0 * 1024.0;

	std::cerr << "--- System Memory Query ---" << std::endl;

	if (total > 0) {
		std::cerr << "Total Physical RAM: " << std::fixed
			  << std::setprecision(2) << (total / bytes_to_gb)
			  << " GB (" << total << " bytes)" << std::endl;
	}
	else {
		std::cerr << "Could not determine Total RAM." << std::endl;
	}

	if (available > 0) {
		std::cerr << "Available Free RAM: " << std::fixed
			  << std::setprecision(2) << (available / bytes_to_gb)
			  << " GB (" << available << " bytes)" << std::endl;
	}
	else {
		std::cerr << "Could not determine Available RAM." << std::endl;
	}

	return 0;
}

}; // namespace meza::mem
