
#include <cmath>
#include <vector>
#include <array>
#include <algorithm>
#include <functional>
#include <numeric>
#include <stdexcept>
#include <random>
#include <typeinfo>

#include <diy/master.hpp>
#include <diy/reduce.hpp>
#include <diy/partners/merge.hpp>
#include <diy/decomposition.hpp>
#include <diy/assigner.hpp>
#include <diy/mpi.hpp>
#include <diy/serialization.hpp>

#include "config.hpp"
#include "GenericBlock.hpp"
#include "Reduce.hpp"
#include "CalcConfig.hpp"

#include "YODA/ReaderYODA.h"
#include "YODA/WriterYODA.h"
#include "YODA/AnalysisObject.h"

#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/HepMC2.h"
#include "Rivet/AnalysisHandler.hh"
#undef foreach // This line prevents a clash of definitions of rivet's legacy foreach with that of DIY

#include "HepMC/IO_GenEvent.h"

using namespace std;
using namespace Pythia8;
#include "opts.h"

using namespace std;

/*
  Note: This is not done.  The reduction is across all configurations.   
  It really needs to be only over same physics / model space parameters.
  Discussed doing this by constructing a custom neighbor relationship
  that only include like physics and model space parameters. 

  Templates that handle standard boilderplate code are needed for:
  1) blocks
  2) management of master, assigner, and decomposer
  3) reduction operation
  4) perhaps defining custom data types to be sent around (see PointConfig)
  5) helpers for the "foreach" processing
 */

typedef diy::DiscreteBounds Bounds;
typedef diy::RegularGridLink RCLink;

typedef std::vector<float> Weights;

typedef std::shared_ptr<YODA::AnalysisObject> AO_ptr;
typedef std::vector<AO_ptr> AnalysisObjects;

// ----This is common to all introduced types - can it be automated? ------
namespace diy {
  template <> struct Serialization<PointConfig>
  {
    static void save(diy::BinaryBuffer& bb, const PointConfig& m)
    { diy::save(bb, &m,sizeof(PointConfig)); }
    static void load(diy::BinaryBuffer& bb, PointConfig& m)
    { diy::load(bb, &m, sizeof(PointConfig)); }
  };

  template <> struct Serialization<AnalysisObjects>
  {
    typedef AnalysisObjects::value_type::element_type data_type;
    typedef std::vector<data_type*> Ptrs;
    
    static void save(diy::BinaryBuffer& bb, AnalysisObjects const& m)
    { 
      std::ostringstream stream;
      Ptrs out(m.size());
      std::transform(std::cbegin(m),std::cend(m), std::begin(out),
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

  namespace mpi {
    namespace detail {
      template<> struct mpi_datatype<PointConfig>
      {
	static MPI_Datatype datatype() { return MPI_BYTE; }
	static const void* address(PointConfig const& x) { return &x; }
	static void* address(PointConfig& x) { return &x; }
	static int count(PointConfig const&)
	{ return sizeof(PointConfig); }
      };
    }
  }
 }

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
      auto& nh = dynamic_cast<T&>(*copy);
      auto const& bh = dynamic_cast<T&>(*other); // Cannot be const when calling scaleW
      nh+=bh;
      if (nh.hasAnnotation("OriginalScaledBy"))
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
      auto& nh = dynamic_cast<T&>(*copy);
      auto const& bh = dynamic_cast<T&>(*other);
     //std::cerr << "Warning, no operator += defined for " << bt->type() << "\n";
      return true;
    }
  else
    return false;
}

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

typedef GenericBlock<Bounds, PointConfig, AnalysisObjects> Block;
typedef ConfigBlockAdder<Bounds, RCLink, Block, PointConfigs> AddBlock;

  
void print_block(Block* b,                             // local block
                 const diy::Master::ProxyWithLink& cp, // communication proxy
                 bool verbose)                         // user-defined additional arguments
{
  if (verbose && cp.gid() == 0)
    {
      //      for (size_t i = 0; i < b->data.size(); ++i)
      //fmt::print(stderr, "({},{}) ", b->data[i], b->buffer[i]);
      //fmt::print(stderr, "\n");
    }
}

void process_block(Block* b, diy::Master::ProxyWithLink const& cp, int rank, std::vector<std::string> analyses)
{
  Rivet::Log::setLevel("Rivet", Rivet::Log::ERROR);
    b->pythia.readString("Print:quiet = on");
    b->pythia.readString("Random:setSeed = on");
    b->pythia.readString("Beams:idA = 11");
    b->pythia.readString("Beams:idB = -11");
    b->pythia.readString("Beams:eCM = 91.2");
    b->pythia.readString("WeakSingleBoson:ffbar2gmZ = on");

    // Py8 random seed for this block read from point config
    b->pythia.readString("Random:setSeed = on");
    std::string seedConf("Random:seed = " + std::to_string(b->state.seed));
    b->pythia.readString(seedConf);

    
    b->pythia.readString("Main:numberOfEvents = " + std::to_string(b->state.num_events));
    //
    int nEvents = b->pythia.mode("Main:numberOfEvents");
    b->pythia.init();

    for (auto a : analyses) 
    {
       //std::cerr << "Adding analysis " << a <<"\n";
       b->ah.addAnalysis(a);
    }

    for (int iEvent = 0; iEvent < nEvents; ++iEvent) {
      // Generate events. Quit if many failures.
      if (!b->pythia.next()) {
        //if (++iAbort < nAbort) continue;
        //cout << " Event generation aborted prematurely, owing to error!\n";
        break;
      }
      HepMC::GenEvent* hepmcevt = new HepMC::GenEvent();
      b->ToHepMC.fill_next_event( b->pythia, hepmcevt );

      try {b->ah.analyze( *hepmcevt ) ;} catch (const std::exception& e) 
      {
     //   fmt::print(stderr, "[{}] exception in analyze: {}\n", cp.gid(), e.what());
      }
      delete hepmcevt;
    }
    b->ah.setCrossSection(b->pythia.info.sigmaGen() * 1.0E9);  
    b->ah.finalize();
    b->data = b->ah.getData();
    
    // Write out so we can sanity check with yodamerge
   
    // This is a bit annoying --- we need to unscale Histo1D and Histo2D beforge the reduction
    for (auto ao : b->data) {
      if (ao->hasAnnotation("ScaledBy"))
      {
        double sc = std::stod(ao->annotation("ScaledBy"));
        if (ao->type()=="Histo1D") 
        {
           dynamic_cast<YODA::Histo1D&>(*ao).scaleW(1./sc);
           dynamic_cast<YODA::Histo1D&>(*ao).addAnnotation("OriginalScaledBy", 1./sc);
        }
        else if (ao->type()=="Histo2D") 
        {
           dynamic_cast<YODA::Histo2D&>(*ao).scaleW(1./sc);
           dynamic_cast<YODA::Histo2D&>(*ao).addAnnotation("OriginalScaledBy", 1./sc);
        }
      }
    }
    //fmt::print(stderr, "[{}] finalised\n", cp.gid());
}


void write_yoda(Block* b, diy::Master::ProxyWithLink const& cp, std::string out_file, int rank)
{
  if (rank==0 && cp.gid()==0) {
    for (auto ao : b->buffer) {
       if (ao->hasAnnotation("OriginalScaledBy"))
       {
         double sc = std::stod(ao->annotation("OriginalScaledBy"));
         if (ao->type()=="Histo1D") 
         {
            dynamic_cast<YODA::Histo1D&>(*ao).scaleW(1./sc);
         }
         else if (ao->type()=="Histo2D") 
         {
            dynamic_cast<YODA::Histo2D&>(*ao).scaleW(1./sc);
         }
       }
    }
    YODA::WriterYODA::write(out_file, b->buffer);
  }
}

// --- main program ---//
int main(int argc, char* argv[])
{
    diy::mpi::environment env(argc, argv); 
    diy::mpi::communicator world;

    int threads = 1;
    int nEventsPerBlock=1000;
    int nEvents=1000;
    size_t seed=1234;
    vector<std::string> analyses;
    std::string out_file="diy.yoda";
    // get command line arguments
    using namespace opts;
    Options ops(argc, argv);
    bool                      verbose     = ops >> Present('v',
                                                           "verbose",
                                                           "verbose output");
    ops >> Option('t', "thread",   threads,       "Number of threads");
    ops >> Option('n', "nevents",  nEvents,       "Number of events to generate in total");
    ops >> Option('a', "analysis",  analyses,       "Rivet analyses");
    ops >> Option('m', "nmin",  nEventsPerBlock,  "Minimum number of events per block.");
    ops >> Option('o', "output",  out_file,  "Output filename.");
    ops >> Option('s', "seed",  seed,  "The Base seed.");
    if (ops >> Present('h', "help", "Show help"))
    {
        std::cout << "Usage: " << argv[0] << " [OPTIONS]\n";
        std::cout << ops;
        return 1;
    }
  
    int mem_blocks  = -1;  // all blocks in memory, if value here then that is how many are in memory

    size_t blocks = world.size() * threads;

    PointConfigs revised;

    if( world.rank()==0 )
      {
        revised = mkSingleRunConfigs(world.size()*threads,nEventsPerBlock, nEvents, seed);
      }

    diy::mpi::broadcast(world, revised, 0);

    if(world.rank()==0 && verbose)
      {
        for(size_t it=0;it<blocks;++it)
          {
            cout << revised[it].psp_id << " " << revised[it].num_events << " "
                 << revised[it].seed << " " << revised[it].physics_id << "\n";
          }
      }

    // -------- above this point is all initialization custom for this application

    // probably do not need this if the domain is used (Bounds)
    // which of the blocks are mine?
    auto interval = blocks/world.size();
    auto my_ndx = interval*(world.rank());
    auto my_start = std::cbegin(revised)+my_ndx;
    auto my_end = (world.rank()==world.size()-1) ? std::end(revised) : my_start+interval;
    
    // ----- starting here is a lot of standard boilerplate code for this kind of
    //       application.
    
    // diy initialization
    diy::FileStorage storage("./DIY.XXXXXX"); // used for blocks moved out of core
    diy::Master master(world, // master is the top-level diy object
  		     threads,
  		     mem_blocks,
  		     &Block::create, &Block::destroy,
  		     &storage,
  		     &Block::save, &Block::load);

    // an object for adding new blocks to master
    AddBlock create(master, my_start, my_end);
    
    //  -------

    int dim(1); // TODO get rid of this
    // changed to discrete bounds and give range of revised
    // set some global data bounds
    // each block is given the configuration slot it processes
    Bounds domain;
    for (int i = 0; i < dim; ++i)
      {
        domain.min[i] = 0;
        domain.max[i] = revised.size()-1;
      }
    
    // choice of contiguous or round robin assigner
    diy::ContiguousAssigner   assigner(world.size(), blocks);
    
    // decompose the domain into blocks
    // This is a DIY regular way to assign neighbors. You can do this manually.
    diy::RegularDecomposer<Bounds> decomposer(dim, domain, blocks);
    decomposer.decompose(world.rank(), assigner, create);

    // ----------- below is the processing for this application
    // threads active here
    master.foreach([world, analyses](Block* b, const diy::Master::ProxyWithLink& cp)
                     {process_block(b, cp, world.rank(), analyses); });

     //this is MPI
     //merge-based reduction: create the partners that determine how groups are formed
     //in each round and then execute the reduction
    int k = 2;       // the radix of the k-ary reduction tree
    // partners for merge over regular block grid
    diy::RegularMergePartners  partners(decomposer,  // domain decomposition
                                        k,           // radix of k-ary reduction
                                        true); // contiguous = true: distance doubling
    // contiguous = false: distance halving
    // reduction
    // this assumes that all blocks are participating in the reduction
    diy::reduce(master,              // Master object
                assigner,            // Assigner object
                partners,            // RegularMergePartners object
                &reduceData<Block>);

    // This is where the write out happens
    master.foreach([world,out_file](Block* b, const diy::Master::ProxyWithLink& cp)
                   { write_yoda(b, cp, out_file, world.rank()); });

    return 0;
}
