#ifndef DIYGEN_Serialisation_HPP
#define DIYGEN_Serialisation_HPP

#include "YODA/ReaderYODA.h"
#include "YODA/WriterYODA.h"
#include "YODA/AnalysisObject.h"
typedef std::shared_ptr<YODA::AnalysisObject> AO_ptr;
typedef std::vector<AO_ptr> AnalysisObjects;


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
    typedef std::vector<data_type*> Ptrs;

    static void save(diy::BinaryBuffer& bb, AnalysisObjects const& m)
    {
      std::ostringstream stream;
      Ptrs out(m.size());
      std::transform(m.cbegin(), m.cend(), out.begin(),
		     [](AnalysisObjects::value_type const& x) { return x.get(); });
      YODA::WriterYODA::write(stream, out);
      std::string s = stream.str();

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

// What follows is for reduction of YODA analysis objects
// TODO: Can we move this into a different file???

// should be able to define
// operator+(AnalysisObject_ptr, AnalysisObject_ptr)
// here and the generic reduce should work.

namespace YODA {
template <typename T>
bool addThisKind(AO_ptr& copy, AO_ptr const& other)
{
  auto const& bt = other.get();

  if(typeid(*bt).hash_code() == typeid(T).hash_code())
    {
      // TODO: this test is causing a lot of warnings at compile time, maybe rewrite?
      try {dynamic_cast<T&>(*copy) ;} catch (const std::exception& e)
      {
        fmt::print(stderr, "\n\n exception in add h1d 1: {} {} {}\n", e.what(), copy->path(), copy->title());
        copy->reset();
        //return false;
      }
      try {dynamic_cast<T&>(*other);} catch (const std::exception& e)
      {
        fmt::print(stderr, "\n\n exception in add h1d 2: {} {} {}\n", e.what(), other->path(), other->title());
        //return false;
      }
      auto& nh = dynamic_cast<T&>(*copy);
      auto const& bh = dynamic_cast<T&>(*other); // Cannot be const when calling scaleW
      nh+=bh;
      if (nh.hasAnnotation("OriginalScaledBy") && bh.hasAnnotation("OriginalScaledBy"))
      {
        double sc_n = std::stod(nh.annotation("OriginalScaledBy"));
        double sc_b = std::stod(bh.annotation("OriginalScaledBy"));
        nh.setAnnotation("OriginalScaledBy", sc_n+sc_b);
      }
      return true;
    }
  else
    return false;
}
}

template <typename T>
bool addCounter(AO_ptr& copy, AO_ptr const& other)
{
  auto const& bt = other.get();
  if(typeid(*bt).hash_code() == typeid(T).hash_code())
    {
      try {dynamic_cast<T&>(*copy) ;} catch (const std::exception& e)
      {
        fmt::print(stderr, "\n\n exception in add h1d 1: {} {} {}\n", e.what(), copy->path(), copy->title());
        copy->reset();
        //return false;
      }
      try {dynamic_cast<T&>(*other);} catch (const std::exception& e)
      {
        fmt::print(stderr, "\n\n exception in add h1d 2: {} {} {}\n", e.what(), other->path(), other->title());
        //return false;
      }
      auto& nh = dynamic_cast<T&>(*copy);
      auto& bh = dynamic_cast<T&>(*other);
      nh+=bh;
      return true;
    }
  else
    return false;
}

// Scatters do not have a += operator as the operation
// is not well defined
template <typename T>
bool addScatters(AO_ptr& copy, AO_ptr const& other)
{
  auto const& bt = other.get();

  if(typeid(*bt).hash_code() == typeid(T).hash_code())
    {
      // auto& nh = dynamic_cast<T&>(*copy);
      // auto const& bh = dynamic_cast<T&>(*other);
     //std::cerr << "Warning, no operator += defined for " << bt->type() << "\n";
      return true;
    }
  else
    return false;
}

// Definiton of a + operator to make general DIY reduction work with yoda objects
namespace YODA {
AO_ptr operator+(AO_ptr const& a, AO_ptr const& b)
{
  AO_ptr n(a->newclone());

  if(!addThisKind<YODA::Histo1D>(n,b)   &&
     !addThisKind<YODA::Histo2D>(n,b)   &&
     !addThisKind<YODA::Profile1D>(n,b) &&
     !addThisKind<YODA::Profile2D>(n,b) &&
     !addCounter<YODA::Counter>(n,b)    &&
     !addScatters<YODA::Scatter1D>(n,b) &&
     !addScatters<YODA::Scatter2D>(n,b) &&
     !addScatters<YODA::Scatter3D>(n,b))
      {
	std::cerr << "in op+ - but no match!!\n";
	throw std::runtime_error("no YODA type match in op+");
      }

  return n;
}
}

#endif
