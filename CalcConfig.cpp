
#include "CalcConfig.hpp"

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

PointConfigs calculate_block_configs(PointConfigs const& pc, Weights const& time_weights,
				     size_t blocks)
{
  PointConfigs revised;
  
  // need weights, configs, blocks as input, PointConfigs as output
  auto w_event_total = calc_sum(pc,time_weights);
  auto tot_weights = std::accumulate(std::cbegin(time_weights),
				     std::cend(time_weights),
				     0.0f);
  auto portions = calc_portions(w_event_total,tot_weights,
				blocks,time_weights);
  auto tot_needed = std::accumulate(std::cbegin(portions),
				    std::cend(portions),
				    0UL);
  // if there are any zeros in the portions, then there are not enough blocks.
  auto min_p = std::find(std::cbegin(portions),
			 std::cend(portions),0UL);
  
  if(time_weights.size() >blocks || min_p!=std::cend(portions) || tot_needed>blocks)
    {
      cerr << "too few blocks to process data.  Need at least "
	   << time_weights.size() << "\n";
      cerr << "or not enough blocks " << tot_needed << ">" << blocks << "\n";
      throw std::runtime_error("cannot configure job for block count and threads");
    }
  
  // go through the point configs, figuring out how many blocks are
  // needed for each one, given the max/block and weight constraints.
  
  auto per_block = [](size_t count, size_t blocks)
    { return ceil(count / blocks); };
  auto excess = [](size_t count, size_t blocks, size_t per)
    { return (per * blocks) - count; };
  
  size_t pos = 0;
  for(size_t i=0;i<pc.size();++i)
    {
      auto events_per = per_block(pc[i].num_events, portions[i]);
      auto over = excess(pc[i].num_events, portions[i], events_per);
      
      for(size_t j=0; j<portions[i]; ++j)
	{
	  revised.push_back({ pc[i].psp_id, static_cast<size_t>(events_per),
		pc[i].seed, pc[i].physics_id });
	  ++pos;
	}
      revised.back().num_events -= over;
    }
  return std::move(revised);
}

// This function calculates how many blocks we need to accumulate
// the desired number of events given number of available ranks
// for one physics configuration.
PointConfigs mkSingleRunConfigs(size_t n_ranks, size_t n_events_per_job, size_t n_events)//, std::string out_file)
{
  PointConfigs revised;

  auto n_perrank = ceil(n_events/n_ranks);

  auto n_blocks = n_ranks;

 
  for (size_t i=0; i<n_blocks;++i)
  {
   int seed = i + 1234;
   revised.push_back({ 1, static_cast<size_t>(n_perrank), seed, 1});
  }
  return std::move(revised);
}
