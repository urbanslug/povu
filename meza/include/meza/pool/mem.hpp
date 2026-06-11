#ifndef MEZA_MEM_HPP
#define MEZA_MEM_HPP

namespace meza::mem
{

// Returns the total physical RAM in bytes
unsigned long long get_total_ram();

// Returns the currently available (free) physical RAM in bytes
unsigned long long get_available_ram();

int print_mem_info();
}; // namespace meza::mem
#endif // MEZA_MEM_HPP
