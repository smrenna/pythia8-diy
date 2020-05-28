
#include "CalcConfig.hpp"

#include <iterator>
#include <cmath>
#include <numeric>
#include <stdexcept>
#include <iostream>
#include <algorithm>

using namespace std;


// This function calculates how many blocks we need to accumulate
// the desired number of events given number of available ranks
// for one physics configuration.
PointConfigs mkSingleRunConfigs(size_t n_ranks, size_t n_events, size_t base_seed, std::vector<std::string> conf, std::vector<std::string> analyses, std::string f_out)
{
  PointConfigs revised;

  auto n_perrank = ceil(n_events/n_ranks);

  auto n_blocks = n_ranks;


  // Distribute the desdired number of events to the blocks
  // The last block on occasion will have to generate slightly more events
  // just so the total number of events is exactly what the user requests.
  for (size_t i=0; i<n_blocks;++i)
  {
    size_t seed = base_seed;
    if (i==n_blocks-1)
    {
      revised.push_back({ 1, static_cast<size_t>(n_events - i*n_perrank), seed, 1, conf, analyses, f_out});
    }
    else
    {
      revised.push_back({ 1, static_cast<size_t>(n_perrank), seed, 1, conf, analyses, f_out});
    }
  }

  // Set config
  return revised;
}

PointConfig mkRunConfig(size_t n_ranks, size_t n_events, size_t base_seed, std::vector<std::string> conf, std::vector<std::string> analyses, std::string f_out)
{
  PointConfig pc;

  auto n_perrank = ceil(n_events/n_ranks);

  // Set config
  return {1, static_cast<size_t>(n_perrank), base_seed, 1, conf, analyses, f_out};
}
