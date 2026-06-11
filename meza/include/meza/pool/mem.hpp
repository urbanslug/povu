#ifndef MEZA_MEM_HPP
#define MEZA_MEM_HPP

// #include "log/log.h"
// #include <cstddef>
// #include <iomanip>
// #include <iostream>

namespace meza::mem
{

// Returns the total physical RAM in bytes
unsigned long long get_total_ram();

// Returns the currently available (free) physical RAM in bytes
unsigned long long get_available_ram();

int print_mem_info();
}; // namespace meza::mem
#endif // MEZA_MEM_HPP
