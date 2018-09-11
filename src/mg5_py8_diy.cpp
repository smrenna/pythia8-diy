
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
#include <diy/partners/broadcast.hpp>
#include <diy/reduce-operations.hpp>



#include "config.hpp"
#include "GenericBlock.hpp"
#include "Reduce.hpp"
#include "Tools.hpp"
#include "CalcConfig.hpp"
#include "Serialisation.hpp"

#include "YODA/ReaderYODA.h"
#include "YODA/WriterYODA.h"
#include "YODA/AnalysisObject.h"

#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/HepMC2.h"
#include "Pythia8Plugins/LHAMadgraph.h" // For MadGraph5 
#include "Rivet/AnalysisHandler.hh"
#undef foreach // This line prevents a clash of definitions of rivet's legacy foreach with that of DIY

#include "HepMC/IO_GenEvent.h"




using namespace std;
using namespace Pythia8;
#include "opts.h"

using namespace std;


typedef diy::DiscreteBounds Bounds;
typedef diy::RegularGridLink RCLink;

typedef GenericBlock<Bounds, PointConfig, AnalysisObjects> Block;
typedef ConfigBlockAdder<Bounds, RCLink, Block> AddBlock;


void print_block(Block* b,                             // local block
                 const diy::Master::ProxyWithLink& cp, // communication proxy
                 bool verbose)                         // user-defined additional arguments
{
	fmt::print(stderr, "Pythia configuration\n");
   for (auto s : b->state.conf) {
      fmt::print(stderr, "\t[{}]: {}\n", cp.gid(), s);
   }
   fmt::print(stderr, "\n");
   fmt::print(stderr, "Rivet analyses:\n");
   for (auto s : b->state.analyses) {
      fmt::print(stderr, "\t[{}]: {} ", cp.gid(), s);
   }
   fmt::print(stderr, "\ndetector configuration:");
   fmt::print(stderr, "\t[{}]: {}\n", cp.gid(), b->state.detector_conf);
   fmt::print(stderr, "output:");
   fmt::print(stderr, "\t[{}]: {}\n", cp.gid(), b->state.f_out);
   fmt::print(stderr, "MadGraph5:\n");
   for (auto s : b->state.mg5_conf) {
      fmt::print(stderr, "\t[{}]: {}\n", cp.gid(), s);
   }
   fmt::print(stderr, "\n");
}

struct MatchPathSeparator
{
	bool operator()( char ch ) const
	{
		return ch == '/';
	}
};

std::string basepath(const std::string& pathname){
	return std::string(
			pathname.begin(),
			std::find_if(pathname.rbegin(), pathname.rend(), MatchPathSeparator()).base()
			);
}


void process_block(Block* b, diy::Master::ProxyWithLink const& cp, int rank,  bool verbose, int npc)
{
	fmt::print(stderr, "{} Start to process\n", cp.gid());
  // This make rivet only report ERRORs
  // TODO: can we have a global flag to steer verbosity of all moving parts?
  if (!verbose) Rivet::Log::setLevel("Rivet", Rivet::Log::ERROR);

  // Minimise pythia's output
  if (verbose) b->pythia.readString("Print:quiet = off");
  else b->pythia.readString("Print:quiet = on");

  // Add MadGraph
  if(b->state.use_mg5) {
	  fmt::print(stderr, "{} I am HERE\n", cp.gid());
	  // if(b->mg5) delete b->mg5;
	  if(!b->mg5) {
		  try{
			  fmt::print(stderr, "{} LHAupMadgraph WAAA\n", cp.gid());
			  b->mg5 = new LHAupMadgraph(&(b->pythia), true, "amcatnlorun_"+cp.gid(), "mg5_aMC");
		  } catch(const std::bad_alloc& e) {
			  fmt::print(stderr, "{} cannot iniatiate LHAupMadgraph\n", cp.gid());
		  }
	  } else {
		  fmt::print(stderr, "{} LHAupMadgraph already there!\n", cp.gid());
	  }
	  for(auto s: b->state.mg5_conf) b->mg5->readString(s);
	  try {
		  b->pythia.setLHAupPtr(b->mg5);
	  } catch (const std::bad_alloc& e){
		  fmt::print(stderr, "{} HERE?", cp.gid());
	  }
  } else {
	  // Configure pythia with a vector of strings
	  for (auto s  : b->state.conf) b->pythia.readString(s);
  }

  fmt::print(stderr, "{} HERE 3 \n", cp.gid());

  // Py8 random seed for this block read from point config
  b->pythia.readString("Random:setSeed = on");
  b->pythia.readString("Random:seed = " + std::to_string(b->state.seed+cp.gid()));
  b->pythia.readString("Main:numberOfEvents = " + std::to_string(b->state.num_events));

  // All configurations done, initialise Pythia
  // b->pythia.initPtrs(); // TODO --- is this really necessary here?
  fmt::print(stderr, "{} HERE 3.2 \n", cp.gid());
  b->pythia.init();
  fmt::print(stderr, "{} HERE 3.3 \n", cp.gid());
  b->pythia.settings.listChanged();

  // Delete the AnalysisHandlerPtr to ensure there is no memory
  if (b->ah)
  {
    delete b->ah;
  }
  b->ah = new Rivet::AnalysisHandler;

  // Add all anlyses to rivet
  // TODO: we may want to feed the "ignore beams" switch as well
  for (auto a : b->state.analyses) b->ah->addAnalysis(a);


  // Update Rivet simulation
  b->ah->setSimulationFile(b->state.detector_conf.c_str());

  // The event loop
  int nAbort = b->pythia.mode("Main:timesAllowErrors");
  int iAbort = 0;
  if (verbose) fmt::print(stderr, "[{}] generating {} events\n", cp.gid(),  b->state.num_events);
  for (int iEvent = 0; iEvent < b->state.num_events; ++iEvent) {
    if (!b->pythia.next()) {
      if (++iAbort < nAbort) continue;
      break;
    }
    HepMC::GenEvent* hepmcevt = new HepMC::GenEvent();
    b->ToHepMC.fill_next_event( b->pythia, hepmcevt );

    try {b->ah->analyze( *hepmcevt ) ;} catch (const std::exception& e)
    {
      if (verbose) fmt::print(stderr, "[{}] exception in analyze: {}\n", cp.gid(), e.what());
    }
    delete hepmcevt;
    if (iEvent%1000 == 0 && cp.gid()==0) fmt::print(stderr, "[{}]  {}/{} \n", cp.gid(),  iEvent, b->state.num_events);;
  }
  if (verbose) fmt::print(stderr, "[{}] finished generating {} events\n", cp.gid(),  b->state.num_events);

  // Event loop is done, set xsection correctly and normalise histos
  b->ah->setCrossSection(b->pythia.info.sigmaGen() * 1.0E9);
  b->ah->finalize();

  // Push histos into block
  b->data = b->ah->getData();

  // Debug write out --- uncomment to write each block's YODA file
  //b->ah->writeData(std::to_string((1+npc)*(b->state.seed+cp.gid()))+".yoda");


  // This is a bit annoying --- we need to unscale Histo1D and Histo2D beforge the reduction
  // TODO: Figure out whether this is really necessary
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
}


void write_yoda(Block* b, diy::Master::ProxyWithLink const& cp, int rank, bool verbose)
{
 if (verbose) fmt::print(stderr, "[{}] -- rank {} sees write_yoda \n", cp.gid(), rank);
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
    if (verbose) fmt::print(stderr, "[{}] -- writing to file {}  \n", cp.gid(), b->state.f_out);
    YODA::WriterYODA::write(b->state.f_out, b->buffer);
  }
}

void clear_buffer(Block* b, diy::Master::ProxyWithLink const& cp, bool verbose)
{
  if (verbose) fmt::print(stderr, "[{}] -- clear buffer  \n", cp.gid());
  b->buffer.clear();
}

void set_pc(Block* b,                             // local block
                 const diy::Master::ProxyWithLink& cp, // communication proxy
                 PointConfig pc)                         // user-defined additional arguments
{
  if (cp.gid() == 0) b->state = pc;
}

// --- main program ---//
int main(int argc, char* argv[])
{
    diy::mpi::environment env(argc, argv);
    diy::mpi::communicator world;

    size_t nBlocks = 0;
    int threads = 1;
    int nEvents=1000;
    size_t seed=1234;
    vector<std::string> analyses;
    std::string out_file="diy.yoda";
    std::string pfile="runPythia.cmd";
    std::string dfile="detector_config.yoda";
    std::string indir=".";
	std::string mfile = "mg5.cmd";
    // get command line arguments
    using namespace opts;
    Options ops(argc, argv);
    bool                      verbose     = ops >> Present('v',
                                                           "verbose",
                                                           "verbose output");
    ops >> Option('t', "thread",    threads,   "Number of threads");
    ops >> Option('n', "nevents",   nEvents,   "Number of events to generate in total");
    ops >> Option('b', "nblocks",   nBlocks,   "Number of blocks");
    ops >> Option('a', "analysis",  analyses,  "Rivet analyses --- can be given multiple times");
    ops >> Option('o', "output",    out_file,  "Output filename.");
    ops >> Option('s', "seed",      seed,      "The Base seed --- this is incremented for each block.");
    ops >> Option('p', "pfile",     pfile,     "Parameter config file for testing __or__ input file name to look for in directory.");
    ops >> Option('d', "dfile",     dfile,     "detector config file for testing __or__ input file name to look for in directory.");
    ops >> Option('m', "mfile",     mfile,     "MadGraph5 config file for testing __or__ input file name to look for in directory.");
    ops >> Option('i', "indir",     indir,     "Input directory with hierarchical configuration files");
    if (ops >> Present('h', "help", "Show help"))
    {
        std::cout << "Usage: " << argv[0] << " [OPTIONS]\n";
        std::cout << ops;
        return 1;
    }

    int nConfigs;
    
    std::vector<std::string> physConfig;
    std::vector<std::vector<std::string> > physConfigs;
    std::vector<std::string> out_files;
    std::vector<std::string> detectorConfigs;
	// MadGraph Configureations
    std::vector<std::string> mg5Config;
    std::vector<std::vector<std::string> > mg5Configs;
    bool f_ok;
	bool use_mg5 = true;
    if( world.rank()==0 ) {
       // Program logic: check whether a single parameter file has been given or
       // a directory.
       f_ok = readConfig(pfile, physConfig, verbose);
       if (!f_ok) {
          // Use glob to look for files in directory
          for (auto f : glob(indir + "/*/" + pfile)) {
             physConfig.clear();
             bool this_ok = readConfig(f, physConfig, verbose);
             if (this_ok) {
                physConfigs.push_back(physConfig);
                out_files.push_back(f+".yoda");
				if(verbose) {
					fmt::print(stderr, "detector config: {}\n", basepath(f)+dfile);
				}
				detectorConfigs.push_back(basepath(f)+dfile);
             }
          }
       }
       else {
          physConfigs.push_back(physConfig);
          out_files.push_back(out_file);
		  detectorConfigs.push_back(dfile);
       }
       nConfigs=physConfigs.size();
	   // check if MadGraph configuration files are there
	   f_ok = readConfig(mfile, mg5Config, verbose);
	   if(!f_ok){
		   for (auto f : glob(indir + "/*/" + mfile)) {
				mg5Config.clear();
				bool this_ok = readConfig(f, mg5Config, verbose);
				if(this_ok) mg5Configs.push_back(mg5Config);
		   }
	   } else {
		   mg5Configs.push_back(mg5Config);
	   }
	   if(mg5Configs.size() != physConfigs.size()){
		   fmt::print(stderr, "Configurations for MadGraph and Pythia do not match\n Not using MadGraph5");
		   mg5Configs.clear();
		   mg5Config.push_back("NONE");
		   mg5Configs.push_back(mg5Config);
		   use_mg5 = false;
	   }
    }

    MPI_Bcast(&nConfigs,   1, MPI_INT, 0, world);



    // ----- starting here is a lot of standard boilerplate code for this kind of
    //       application.
    int mem_blocks  = -1;  // all blocks in memory, if value here then that is how many are in memory
    int dim(1);
    
    size_t blocks;
    if (nBlocks==0) blocks= world.size() * threads;
    else blocks=nBlocks;
    // diy initialization
    diy::FileStorage storage("./DIY.XXXXXX"); // used for blocks moved out of core
    Bounds domain;
    for (int i = 0; i < dim; ++i) {
      domain.min[i] = 0;
      domain.max[i] = blocks-1;
    }
    ////// choice of contiguous or round robin assigner
    diy::ContiguousAssigner   assigner(world.size(), blocks);
    //// decompose the domain into blocks
    //// This is a DIY regular way to assign neighbors. You can do this manually.
    diy::RegularDecomposer<Bounds> decomposer(dim, domain, blocks);

    int k = 2;       // the radix of the k-ary reduction tree
    diy::RegularBroadcastPartners comm(decomposer, k, true);

    diy::RegularMergePartners  partners(decomposer,  // domain decomposition
                                        k,           // radix of k-ary reduction
                                        true); // contiguous = true: distance doubling


    diy::Master master(world, threads, mem_blocks, &Block::create, &Block::destroy,
                     &storage, &Block::save, &Block::load);
    AddBlock create(master);
    decomposer.decompose(world.rank(), assigner, create); // Note: only decompose once!

    if( world.rank()==0 ) {
      fmt::print(stderr, "\n*** This is diy running Pythia8 ***\n");
      fmt::print(stderr, "\n    Blocks:                  {}\n", blocks);
      fmt::print(stderr, "\n    Physics configurations:  {}\n", nConfigs);
      fmt::print(stderr, "\n    Number of events/config: {}\n", nEvents);
      fmt::print(stderr, "\n    Total number of events:  {}\n", nEvents*nConfigs);
      fmt::print(stderr, "\n    World size:  {}\n", world.size());
      fmt::print(stderr, "\n    MadGraph configurations:  {}\n", mg5Configs.size());
      fmt::print(stderr, "***********************************\n");
    }

    PointConfig pc;
    for (size_t ipc=0;ipc<nConfigs;++ipc) {
       if (world.rank()==0) {
		   if(use_mg5) {
			   pc = mkRunConfig(blocks, nEvents, seed, physConfigs[ipc], analyses, out_files[ipc], detectorConfigs[ipc], mg5Configs[ipc], true);
		   } else {
			   pc = mkRunConfig(blocks, nEvents, seed, physConfigs[ipc], analyses, out_files[ipc], detectorConfigs[ipc], mg5Configs[0], false);
		   }
       }

       // We need to tell the first block about the new configuration
       master.foreach([world, pc](Block* b, const diy::Master::ProxyWithLink& cp)
                        {set_pc(b, cp, pc); });

       // Broadcast the runconfig to the other blocks
       diy::reduce(master, assigner, comm, &bc_pointconfig<Block>);

       if (verbose) master.foreach([world](Block* b, const diy::Master::ProxyWithLink& cp)
                        {print_block(b, cp, world.rank()); });

       master.foreach([world, verbose, ipc](Block* b, const diy::Master::ProxyWithLink& cp)
                        {process_block(b, cp, world.rank(), verbose, ipc); });


       diy::reduce(master,              // Master object
                   assigner,            // Assigner object
                   partners,            // RegularMergePartners object
                   &reduceData<Block>);

       //////This is where the write out happens --- the output file name is part of the blocks config
       master.foreach([world, verbose](Block* b, const diy::Master::ProxyWithLink& cp)
                      { write_yoda(b, cp, world.rank(), verbose); });

       // Wipe the buffer to prevent double counting etc.
       master.foreach([world, verbose](Block* b, const diy::Master::ProxyWithLink& cp)
                      { clear_buffer(b, cp, verbose); });
   }


    return 0;
}
