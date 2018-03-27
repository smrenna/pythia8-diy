#ifndef DIYGEN_CalcConfig_HPP
#define DIYGEN_CalcConfig_HPP 1

#include <vector>
#include "config.hpp"

std::vector<PointConfig> mkSingleRunConfigs(size_t n_ranks,  size_t n_events, size_t base_seed, std::vector<std::string> cnf, std::string f_out);
#endif
