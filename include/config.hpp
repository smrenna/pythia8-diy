#ifndef DG_CONFIG_HPP
#define DG_CONFIG_HPP 1

#include <vector>
#include <ostream>
#include <istream>
#include <string>

struct PointConfig
{
  size_t psp_id;
  size_t num_events;
  size_t seed; // not correct type
  size_t physics_id; // not correct information
  std::vector<std::string> conf;
  std::vector<std::string> analyses;
  std::string f_out;
};

typedef std::vector<PointConfig> PointConfigs;

inline std::ostream& operator<<(std::ostream& ost, PointConfig const& pc)
{
  ost << pc.psp_id << " " << pc.num_events << " "
      << pc.seed << " " << pc.physics_id << " " << pc.f_out;
  return ost;
}

inline std::istream& operator>>(std::istream& ist, PointConfig& pc)
{
  ist >> pc.psp_id >> pc.num_events >> pc.seed >> pc.physics_id;
  return ist;
}

inline bool operator==(PointConfig const& a, PointConfig const& b)
{
  return a.psp_id==b.psp_id && a.seed==b.seed
    && a.num_events==b.num_events && a.physics_id==b.physics_id;
}

PointConfigs read_config(std::string const& fname);
std::vector<float> read_weights(std::string const& fname);

void run_config_test();

#endif
