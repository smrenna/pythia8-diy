#ifndef DIYGEN_CalcConfig_HPP
#define DIYGEN_CalcConfig_HPP 1

#include <vector>
#include "config.hpp"

std::vector<PointConfig> mkSingleRunConfigs(size_t n_ranks,  size_t n_events, size_t base_seed, std::vector<std::string> cnf, std::vector<std::string> analyses, std::string f_out);
PointConfig mkRunConfig(size_t n_ranks,  size_t n_events, size_t base_seed, std::vector<std::string> cnf, std::vector<std::string> analyses, std::string f_out, std::string f_det, std::vector<std::string> mg5_cnf, bool use_mg5);
PointConfig mkRunConfig(size_t n_ranks,  size_t n_events, size_t base_seed, std::vector<std::string> cnf, std::vector<std::string> analyses, std::string f_out, std::string f_det="");
#endif
