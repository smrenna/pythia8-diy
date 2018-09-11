
#include "CalcConfig.hpp"

#include <iterator>
#include <cmath>
#include <numeric>
#include <stdexcept>
#include <iostream>
#include <algorithm>

using namespace std;

namespace {

  typedef std::vector<float> Weights;

  size_t calc_sum(PointConfigs const& pcc, Weights const& tw)
  {
    size_t w_event_total = 0;
    for(size_t i=0;i<pcc.size();++i)
      w_event_total += ceil(tw[i] * pcc[i].num_events);
    return w_event_total;
  }
  
  std::vector<size_t> calc_portions(size_t total_w_events,
				    float weight_total,
				    size_t total_blocks,
				    Weights const& w)
  {
    std::vector<size_t> portions;
    
    for(auto& it : w)
      {
	auto tmp = round(it/weight_total * (float)total_blocks);
	portions.push_back( tmp );
      }
    return portions;
  }

}

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
  return std::move(revised);
}

PointConfig mkRunConfig(size_t n_ranks, size_t n_events, size_t base_seed, std::vector<std::string> conf,
		std::vector<std::string> analyses, std::string f_out, std::string f_det, std::vector<std::string> mg5_cnf, bool use_mg5)
{
  auto n_perrank = ceil(n_events/n_ranks);

  // Set config
  return {1, static_cast<size_t>(n_perrank), base_seed, 1, conf, analyses, f_out, f_det, mg5_cnf, use_mg5};
}

PointConfig mkRunConfig(size_t n_ranks, size_t n_events, size_t base_seed, std::vector<std::string> conf,
		std::vector<std::string> analyses, std::string f_out, std::string f_det)
{
  auto n_perrank = ceil(n_events/n_ranks);

  std::vector<std::string> mg5_cnf;
  mg5_cnf.push_back("NONE");
  // Set config
  return {1, static_cast<size_t>(n_perrank), base_seed, 1, conf, analyses, f_out, f_det, mg5_cnf, false};
}
