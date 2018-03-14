#ifndef DIYGEN_CalcConfig_HPP
#define DIYGEN_CalcConfig_HPP 1

#include <vector>
#include "config.hpp"

std::vector<PointConfig> calculate_block_configs(std::vector<PointConfig> const& pc,
						 std::vector<float> const& time_weights,
						 size_t blocks);

std::vector<PointConfig> mkSingleRunConfigs(size_t n_ranks, size_t n_events_per_job, size_t n_events);
#endif
