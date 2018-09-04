
#include "config.hpp"

#include <vector>
#include <stdexcept>
#include <sstream>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <iterator>

#include <cstdio>
#include <unistd.h>

PointConfigs read_config(std::string const& fname)
{
  std::ifstream ist(fname.c_str());
  if(!ist) throw std::runtime_error("Bad file name in read_config");
  std::istream_iterator<PointConfig> iter(ist);
  std::istream_iterator<PointConfig> eof;
  PointConfigs pc(iter,eof);
  return pc;
}

std::vector<float> read_weights(std::string const& fname)
{
  std::ifstream ist(fname.c_str());
  if(!ist) throw std::runtime_error("Bad file name in read_weights");
  std::istream_iterator<double> iter(ist);
  std::istream_iterator<double> eof;
  std::vector<float> f(iter,eof);
  return f;
}

namespace {
  const PointConfigs pc = { {1,10000, 1,1}, {2,20000, 2,1},
			    {3,10000, 3,1}, {4,15000, 4,1} };
  
  const std::vector<float> time_weights = { 1.21, 1.43, 2.13, 0.54 }; 

  template <typename T, typename F>
  void test_one(std::string const& s, F f, T const& c)
  {
    char name[100] = "testconfig_reading_XXXXXXX";
    char* tname = mktemp(name);
    if(tname !=name) throw std::runtime_error("no temp file made in test_one");
    std::string c_file(tname);
    std::ofstream c_ost(c_file.c_str());
    if(!c_ost) throw std::runtime_error("cannot open the output file in test_one");
    c_ost << s << "\n";
    c_ost.close();

    T ans = f(c_file);
    std::remove(tname);

    if(ans!=c)
    {
      std::cerr << "configuration read not what was expected\n";
      for(auto const& i : c) std::cerr << i << " ";
      std::cerr << "\n not equal to \n";
      for(auto const& i : ans) std::cerr << i << " ";
      std::cerr <<"\n";
      throw std::runtime_error("Bad config reader");
    }
  } 
}

void run_config_test()
{
  std::string pcs = "1 10000 1 1\n2 20000 2 1\n3 10000 3 1\n4 15000 4 1";
  std::string ws  = "1.21\n1.43\n2.13\n0.54";

  std::istringstream ist_c(pcs.c_str());
  std::istringstream ist_w(ws.c_str());

  std::istream_iterator<PointConfig> iter(ist_c);
  std::istream_iterator<PointConfig> eof;

  PointConfigs ans(iter,eof);

  if(ans!=pc)
    {
      std::cerr << "configuration read not what was expected\n";
      for(auto const& i : pc) std::cerr << i << " ";
      std::cerr << "\n not equal to \n";
      for(auto const& i : ans) std::cerr << i << " ";
      std::cerr <<"\n";
      throw std::runtime_error("Bad config reader");
    }

  try{
    test_one(pcs, read_config, pc);
    test_one(ws, read_weights, time_weights);
  }
  catch( std::exception& e)
    {
      std::cerr << e.what();
      throw;
    }
}
