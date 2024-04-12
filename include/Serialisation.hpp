#ifndef DIYGEN_Serialisation_HPP
#define DIYGEN_Serialisation_HPP

#include "YODA/ReaderYODA.h"
#include "YODA/IO.h"
#include "YODA/WriterYODA.h"
#include "YODA/AnalysisObject.h"
#include "Rivet/AnalysisHandler.hh"
#include "Rivet/Tools/RivetYODA.hh"

//typedef std::shared_ptr<YODA::AnalysisObject> AO_ptr;
//typedef YODA::AO* AO_ptr;
typedef std::vector<YODA::AnalysisObjectPtr> AnalysisObjects;

typedef Rivet::AnalysisHandler AH_ptr;
//typedef std::vector<AH_ptr> AnalysisHandlers;

// ----This is common to all introduced types - can it be automated? ------
namespace diy {
  // Serialisation for the PointConfig (seed and nEvent info mostly))
  template <> struct Serialization<PointConfig>
  {
    static void save(diy::BinaryBuffer& bb, const PointConfig& m)
    {

       diy::save(bb, m.psp_id);
       diy::save(bb, m.num_events);
       diy::save(bb, m.seed);
       diy::save(bb, m.physics_id);
       diy::save(bb, m.conf.size());
       for (auto s : m.conf) diy::save(bb, s);
       diy::save(bb, m.analyses.size());
       for (auto s : m.analyses) diy::save(bb, s);
       diy::save(bb, m.f_out);
    }
    static void load(diy::BinaryBuffer& bb, PointConfig& m)
    {
       size_t size;
       diy::load(bb, size);
       m.psp_id=size;
       diy::load(bb, size);
       m.num_events=size;
       diy::load(bb, size);
       m.seed=size;
       diy::load(bb, size);
       m.physics_id=size;
       diy::load(bb, size);
       std::string temp;
       m.conf.resize(size);
       for (size_t i=0;i<size;++i) { // Iteration over now known number of elements
          // Elementwise loading of string from bb using diy's native load for std::string
          diy::load(bb, temp);
          m.conf[i] = temp;
       }
       diy::load(bb, size);
       m.analyses.resize(size);
       for (size_t i=0;i<size;++i) { // Iteration over now known number of elements
          // Elementwise loading of string from bb using diy's native load for std::string
          diy::load(bb, temp);
          m.analyses[i] = temp;
       }
       diy::load(bb, temp);
       m.f_out = temp;
    }
  };
  

  // This is JBKs tricking of the low level DIY serialisation to work with PointConfig
  //namespace mpi {
    //namespace detail {
      //template<> struct mpi_datatype<PointConfig>
      //{
        //static MPI_Datatype datatype() { return MPI_BYTE; }
        //static const void* address(PointConfig const& x) { return &x; }
        //static void* address(PointConfig& x) { return &x; }
        //static int count(PointConfig const&)
        //{ return sizeof(PointConfig); }
      //};
    //}
  //}

  // Serialisation for all kinds of YODA objects
  // TODO: Can we move this into a different file???
  template <> struct Serialization<AnalysisObjects>
  {
    typedef AnalysisObjects::value_type::element_type data_type;
    //typedef std::vector<data_type*> Ptrs;
    typedef std::vector<YODA::AO*> Ptrs;

    static void save(diy::BinaryBuffer& bb, AnalysisObjects const& m)
    {

      std::ostringstream stream_a;

//      YODA::WriterYODA::write(stream_a, begin(m), end(m), "yoda");
      YODA::write(stream_a, begin(m), end(m), "yoda");

    //  std::string s = stream_a.str();

      // Unclear to SM why this was here
      
//      const vector<YODA::AnalysisObjectPtr> out(m.size());
//      std::transform(m.cbegin(), m.cend(), out.begin(),
//		  [](AnalysisObjects::value_type const& x) { return x.get(); });

//      YODA::WriterYODA::write(stream_a, begin(out), end(out) );
//      YODA::write(stream_a, begin(out), end(out), "");

      std::string s = stream_a.str(); 
      

      diy::save(bb, s.size());
      diy::save(bb, s.c_str(), s.size());
    }
    static void load(diy::BinaryBuffer& bb, AnalysisObjects& m)
    {
      size_t sz;
      diy::load(bb,sz); // sz is a return argument
      std::string str;
      str.resize(sz); // TODO use Marc's way

      diy::load(bb, str.data(), sz);
      std::istringstream stream(str);
      auto tmp = YODA::ReaderYODA::read(stream);
      AnalysisObjects in(tmp.begin(),tmp.end());
      m.swap(in);
    }
  };
 
}

namespace YODA {
  void addVectors(AnalysisObjects &a, AnalysisObjects const &b) {  

    std::ostringstream stream_a;
    //YODA::WriterYODA::write(stream_a, begin(a), end(a));
    YODA::write(stream_a, begin(a), end(a), "yoda");
    std::string s_a = stream_a.str();

    std::ostringstream stream_b;
    //YODA::WriterYODA::write(stream_b, begin(b), end(b));
    YODA::write(stream_b, begin(b), end(b), "yoda");
    std::string s_b = stream_b.str();

    //MRENNA
    bool preload(false);
    string fmt("yoda");

    std::istringstream stream(s_a);

    AH_ptr ahmerge;
    ahmerge.readData(stream, fmt, preload);

    std::istringstream streamp(s_b);

    AH_ptr ahtemp;
    ahtemp.readData(streamp, fmt, preload);

    ahmerge.merge( ahtemp );
    ahmerge.finalize();

    a = ahmerge.getYodaAOs(true);
    /*
    */
    
  };
}   

#endif
