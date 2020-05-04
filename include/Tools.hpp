#ifndef DIYGEN_Tools_HPP
#define DIYGEN_Tools_HPP


#include <glob.h>
// Get all files from subdirectories following a certain pattern
// https://stackoverflow.com/questions/8401777/simple-glob-in-c-on-unix-system
inline std::vector<std::string> glob(const std::string& pat)
{
  using namespace std;
  glob_t glob_result;
  glob(pat.c_str(),GLOB_TILDE,NULL,&glob_result);
  vector<string> ret;
  for(unsigned int i=0;i<glob_result.gl_pathc;++i){
      ret.push_back(string(glob_result.gl_pathv[i]));
  }
  globfree(&glob_result);
  return ret;
}

// Read a config file, ignore empty lines and # commented lines
bool readConfig(std::string fname, std::vector<std::string> & vconf, bool verbose)
{
  std::ifstream f(fname.c_str());

    if (verbose) std::cerr << "Opening file '"<<fname<<"'\n";
  // Check if object is valid
  if(!f)
  {
    if (verbose) std::cerr << "Error opening file '"<<fname<<"'\n";
    return false;
  }

  std::string temp;
  // Read the next line from File untill it reaches the end.
  while (std::getline(f, temp))
  {
    // Line contains string of length > 0 then save it in vector
    if (temp.size() > 0 && temp.find("#")!=0) vconf.push_back(temp);
  }
  f.close();
  return true;
}

constexpr bool is_powerof2(int v) {
    return v && ((v & (v - 1)) == 0);
}
#endif
