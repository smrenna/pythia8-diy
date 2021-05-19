
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

//#include "PtjTMSdefinitionHooks.h" // Stefan's unlops hook

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
#include "Pythia8/LHEF3.h"
#include "Pythia8Plugins/HepMC2.h"
#include "Pythia8Plugins/LHEH5.h"
#include "Pythia8Plugins/LHAHDF5.h"
// Include UserHooks for Jet Matching.
#include "Pythia8Plugins/CombineMatchingInput.h"
// Include UserHooks for randomly choosing between integrated and
// non-integrated treatment for unitarised merging.
#include "Pythia8Plugins/aMCatNLOHooks.h"
//
//
#include "Rivet/AnalysisHandler.hh"
#undef foreach // This line prevents a clash of definitions of rivet's legacy foreach with that of DIY

#include "HepMC/IO_GenEvent.h"

#include "fastjet/ClusterSequence.hh" // This is to quieten fastjet

using namespace std;
using namespace Pythia8;

#include "opts.h"

typedef diy::DiscreteBounds Bounds;
typedef diy::RegularGridLink RCLink;

typedef GenericBlock<Bounds, PointConfig, AnalysisObjects> Block;
typedef ConfigBlockAdder<Bounds, RCLink, Block> AddBlock;

void print_block(Block* b,                             // local block
                 const diy::Master::ProxyWithLink& cp, // communication proxy
                 bool verbose)                         // user-defined additional arguments
{
   for (auto s : b->state.conf) {
      fmt::print(stderr, "[{}]: {}\n", cp.gid(), s);
   }
   fmt::print(stderr, "\n");
   for (auto s : b->state.analyses) {
      fmt::print(stderr, "[{}]: {} ", cp.gid(), s);
   }
   fmt::print(stderr, "\n");
}
void process_block_lhe_performance_same(Block* b, diy::Master::ProxyWithLink const& cp, int size, int rank,  bool verbose, int npc, string in_file, int nMax)
{
  // Minimise pythia's output
  //b->pythia.readString("Print:quiet = on");

  // Configure pythia with a vector of strings
  for (auto s  : b->state.conf) b->pythia.readString(s);
  // Py8 random seed for this block read from point config
  b->pythia.readString("Random:setSeed = on");
  b->pythia.readString("Random:seed = " + std::to_string(b->state.seed));
  // TODO seed mode 3, i.e. draw nrank random numbers

  HighFive::File file(in_file, HighFive::File::ReadOnly);  
  hid_t dspace = H5Dget_space(file.getDataSet("event/start").getId());
  
  size_t eventOffset = 0;

  int readSize(nMax);

  shared_ptr<LHAupH5> lhaUpPtr = make_shared<LHAupH5>(&file, eventOffset, readSize);
  
  // Create an LHAup object that can access relevant information in pythia.
 // LHAupH5* LHAup = new LHAupH5( &file , eventOffset, nMax, nMax, verbose, true, true);

  b->pythia.settings.mode("Beams:frameType", 5);
  // Give the external reader to Pythia
  b->pythia.setLHAupPtr( lhaUpPtr );

//  double sigmaTotal  = 0.;
//  double errorTotal  = 0.;
  double xs = 0.;
  for (int i=0; i < b->pythia.info.nProcessesLHEF(); ++i)
    xs += b->pythia.info.sigmaLHEF(i);
  if (verbose) fmt::print(stderr, "[{}] xs: {}\n", cp.gid(), xs);
//  double sigmaSample = 0., errorSample = 0.; 

  // All configurations done, initialise Pythia
  b->pythia.init();

  if (verbose) fmt::print(stderr, "[{}] starting event loop\n", cp.gid());
  // The event loop
  int nAbort = 5;
  int iAbort = 0;
  if (verbose)  fmt::print(stderr, "[{}] generating  {} events\n", cp.gid(),lhaUpPtr->getSize());
  for (size_t iEvent = 0; iEvent < nMax; ++iEvent) {
    if (verbose) fmt::print(stderr, "[{}] is at event {}\n", cp.gid(), iEvent);
    if (!b->pythia.next()) {
       // Gracefully ignore events with 0 weight
       if (lhaUpPtr->weight() == 0) {
          if (verbose) fmt::print(stderr, "[{}] encounters and ignores event {} as it has zero weight\n", cp.gid(), iEvent);
          continue;
       }
       else {
          if (++iAbort < nAbort) continue; // All other errors contribute to the abort counter
       }
      break;
    }
    //LHAup->listEvent();
    if (verbose) lhaUpPtr->listEvent();
  }

}
void process_block_lhe_performance(Block* b, diy::Master::ProxyWithLink const& cp, int size, int rank,  bool verbose, int npc, string in_file, int nMax)
{
  // Minimise pythia's output
  b->pythia.readString("Print:quiet = on");

  // Configure pythia with a vector of strings
  for (auto s  : b->state.conf) b->pythia.readString(s);
  // Py8 random seed for this block read from point config
  b->pythia.readString("Random:setSeed = on");
  b->pythia.readString("Random:seed = " + std::to_string(b->state.seed+cp.gid()));
  // TODO seed mode 3, i.e. draw nrank random numbers

  HighFive::File file(in_file, HighFive::File::ReadOnly);  
  hid_t dspace = H5Dget_space(file.getDataSet("event/start").getId());
  size_t nEvents  =  H5Sget_simple_extent_npoints(dspace);
  if (nMax > 0 && nMax < nEvents) {
     nEvents = nMax;
  }
  size_t ev_rank = floor(nEvents/size);
  size_t eventOffset = rank*ev_rank;
  
  // Create an LHAup object that can access relevant information in pythia.
  //  LHAupH5* LHAup = new LHAupH5( &file , eventOffset, ev_rank, nEvents, verbose);
  shared_ptr<LHAupH5> lhaUpPtr = make_shared<LHAupH5>(&file, eventOffset, ev_rank); 

  b->pythia.settings.mode("Beams:frameType", 5);
  // Give the external reader to Pythia
  b->pythia.setLHAupPtr(lhaUpPtr);

//  double sigmaTotal  = 0.;
//  double errorTotal  = 0.;
  double xs = 0.;
  for (int i=0; i < b->pythia.info.nProcessesLHEF(); ++i)
    xs += b->pythia.info.sigmaLHEF(i);
  if (verbose) fmt::print(stderr, "[{}] xs: {}\n", cp.gid(), xs);
//  double sigmaSample = 0., errorSample = 0.;

  // All configurations done, initialise Pythia
  b->pythia.init();


  if (verbose) fmt::print(stderr, "[{}] starting event loop\n", cp.gid());
  // The event loop
  int nAbort = 5;
  int iAbort = 0;
  if (verbose)  fmt::print(stderr, "[{}] generating  {} events\n", cp.gid(), lhaUpPtr->getSize());
  for (size_t iEvent = 0; iEvent < lhaUpPtr->getSize(); ++iEvent) {
    if (verbose) fmt::print(stderr, "[{}] is at event {}\n", cp.gid(), iEvent);
    if (!b->pythia.next()) {
       // Gracefully ignore events with 0 weight
       if (lhaUpPtr->weight() == 0) {
          if (verbose) fmt::print(stderr, "[{}] encounters and ignores event {} as it has zero weight\n", cp.gid(), iEvent);
          continue;
       }
       else {
          if (++iAbort < nAbort) continue; // All other errors contribute to the abort counter
       }
      break;
    }
    if (verbose ) lhaUpPtr->listEvent();
  }
}

//void process_block(Block* b, diy::Master::ProxyWithLink const& cp, int rank, std::vector<std::string> physConfig, std::vector<std::string> analyses, bool verbose)
void process_block_lhe(Block* b, diy::Master::ProxyWithLink const& cp, int size, int rank,  bool verbose, int npc, string in_file, int nMax, size_t batchSize)
{
   // TODO does that work???
   //if (cp.gid()>0)  {
      //fastjet::ClusterSequence ___;
      //___.set_fastjet_banner_stream(0);
   //}
  // This makes rivet only report ERRORs
  // TODO: can we have a global flag to steer verbosity of all moving parts?
  if (!verbose) Rivet::Log::setLevel("Rivet", Rivet::Log::WARNING);

  // Minimise pythia's output
  if (verbose) b->pythia.readString("Print:quiet = off");
  else b->pythia.readString("Print:quiet = on");

  // Configure pythia with a vector of strings
  for (auto s  : b->state.conf) b->pythia.readString(s);
  // Py8 random seed for this block read from point config
  b->pythia.readString("Random:setSeed = on");
  b->pythia.readString("Random:seed = " + std::to_string(b->state.seed+cp.gid()));
  // TODO seed mode 3, i.e. draw nrank random numbers

  HighFive::File file(in_file, HighFive::File::ReadOnly);  
  hid_t dspace = H5Dget_space(file.getDataSet("event/start").getId());
  size_t nEvents  =  H5Sget_simple_extent_npoints(dspace);
  if (nMax > 0 && nMax < nEvents) nEvents = nMax;
  size_t ev_rank = floor(nEvents/size); // Total number of events this block/rank processes
  size_t eventOffset = rank*ev_rank; //

  fmt::print(stderr, "[{}] reads {}/{} events starting at {} in batches of {}\n", cp.gid(), ev_rank, nEvents, eventOffset, batchSize);

  b->pythia.settings.mode("Beams:frameType", 5);

  shared_ptr<LHAupH5> lhaUpPtr = make_shared<LHAupH5>(&file, eventOffset, ev_rank);  
  // Create an LHAup object that can access relevant information in pythia.
  //  LHAupH5* LHAup = new LHAupH5( &file , eventOffset, ev_rank, nEvents, verbose, true, false,  batchSize);
  // Give the external reader to Pythia
  b->pythia.setLHAupPtr(lhaUpPtr);

  lhaUpPtr->setBatchSize( batchSize );
  // All configurations done, initialise Pythia
  b->pythia.init();

  double sigmaTotal  = 0.;
  double errorTotal  = 0.;
  double xs = 0.;
  for (int i=0; i < b->pythia.info.nProcessesLHEF(); ++i)
    xs += b->pythia.info.sigmaLHEF(i);
  if (verbose) fmt::print(stderr, "[{}] xs: {}\n", cp.gid(), xs);
  double sigmaSample = 0., errorSample = 0.; // NOTE in Stefan's unlops this is reset for each to be merged multi

  // Delete the AnalysisHandlerPtr to ensure there is no memory confusion
  if (b->ah) delete b->ah;
  
  b->ah = new Rivet::AnalysisHandler;
  b->ah->setIgnoreBeams();

  // Add all analyses to rivet
  // TODO: we may want to feed the "ignore beams" switch as well
  for (auto a : b->state.analyses) {
     b->ah->addAnalysis(a);
     if (verbose) fmt::print(stderr, "[{}] add  ######## analysis {}\n", cp.gid(), a);
  }

  if (verbose) fmt::print(stderr, "[{}] starting event loop\n", cp.gid());
  // The event loop
  int nAbort = 5;
  int iAbort = 0;


  for (auto w : lhaUpPtr->_work) {
     fmt::print(stderr, "[{}] next batch\n", cp.gid());
     lhaUpPtr->setBatch(w);

     fmt::print(stderr, "[{}] generating  {} events\n", cp.gid(),  lhaUpPtr->getSize());
     for (size_t iEvent = 0; iEvent < lhaUpPtr->getSize(); ++iEvent) {
       if (verbose) fmt::print(stderr, "[{}] is at event {}\n", cp.gid(), iEvent+w[0]);
       if (b->pythia.next()) {
          if (verbose) b->pythia.stat();
          if (verbose) lhaUpPtr->listEvent();
          // Gracefully ignore events with 0 weight
          if (lhaUpPtr->weight() == 0) {
             if (verbose) fmt::print(stderr, "[{}] encounters and ignores event {} as it has zero weight\n", cp.gid(), iEvent+w[0]);
             continue;
          }
       } else {
             if (++iAbort < nAbort) {
                continue; // All other errors contribute to the abort counter
             } else {
                break;
             }
       }
       //if (verbose && iEvent < 2 ) LHAup->listEvent();
       if (verbose ) lhaUpPtr->listEvent();
       if (verbose) fmt::print(stderr, "[{}] event weight {} \n", cp.gid(), b->pythia.info.weight());
       //if (verbose) fmt::print(stderr, "[{}] event weight {} {} {}\n", cp.gid(), LHAup->weight(), b->pythia.info.weight(), b->pythia.info.eventWeightLHEF);
         
       // Get event weight(s).
       double evtweight         = b->pythia.info.weight();
       // Additional PDF/alphaS weight for internal merging.
       evtweight               *= lhaUpPtr->weight();
       // Additional weight due to random choice of reclustered/non-reclustered
       // treatment. Also contains additional sign for subtractive samples.
                                   //*hardProcessBookkeepingPtr->getNormFactor(); // NOTE this would be necessary for NLO
       if (verbose) fmt::print(stderr, "[{}] after weight {} \n", cp.gid(), evtweight);

       HepMC::GenEvent* hepmcevt = new HepMC::GenEvent();

       
       // weight push_back
       // Work with weighted (LHA strategy=-4) events.
       //double normhepmc = 1.;
       //if (abs(b->pythia.info.lhaStrategy()) == 4)  // Triggering in l 128
         ////normhepmc = 1. / double(1e9*nEvent);
         //normhepmc = 1. / double(LHAup->nTrials());
         //// Work with unweighted events.
       //else
         ////normhepmc = xs / double(1e9*nEvent);
         //normhepmc = xs / double(LHAup->nTrials());

       // SM probably not getTrials, but total number of events in H5 file
       double normhepmc = 1. / double(lhaUpPtr->getTrials()); //
       normhepmc = 1.0;


// SM comment out temporarily until I understand weights
       hepmcevt->weights().push_back(evtweight*normhepmc);
       if (verbose) fmt::print(stderr, "[{}] norm weight {} \n", cp.gid(), evtweight*normhepmc);
       b->ToHepMC.set_print_inconsistency(false);
       b->ToHepMC.set_free_parton_exception(false);
       b->ToHepMC.fill_next_event( b->pythia, hepmcevt );

       sigmaTotal  += evtweight*normhepmc;
       sigmaSample += evtweight*normhepmc;
       errorTotal  += pow2(evtweight*normhepmc);
       errorSample += pow2(evtweight*normhepmc);
       // Report cross section to hepmc
       HepMC::GenCrossSection xsec;
// SM comment out temporarily until I understand weights
//       xsec.set_cross_section( sigmaTotal, b->pythia.info.sigmaErr() );
//       hepmcevt->set_cross_section( xsec );
       if (verbose) fmt::print(stderr, "[{}] xsec {} \n", cp.gid(), sigmaTotal);

       // Here more
       try {b->ah->analyze( *hepmcevt ) ;} catch (const std::exception& e)
       {
         if (verbose) fmt::print(stderr, "[{}] exception in analyze: {}\n", cp.gid(), e.what());
       }
       delete hepmcevt;
       // TODO: make 1000 a free parameter a la msg_every
       if (iEvent%100 == 0 && cp.gid()==0) {
          if ( (b->state.num_events <0) | b->state.num_events > lhaUpPtr->getSize()) {
             fmt::print(stderr, "[{}]a  {}/{} \n", cp.gid(),  iEvent+w[0], lhaUpPtr->getSize());
          }
          else {
             fmt::print(stderr, "[{}]b  {}/{} \n", cp.gid(),  iEvent+w[0], b->state.num_events);
          }
       }

       //if (iEvent >= b->state.num_events && b->state.num_events>=0) {
          //fmt::print(stderr, "[{}] exiting event loop after: {}/{}\n", cp.gid(), iEvent, b->state.num_events);
          //break;
       //}
     }
  }

  // Event loop is done, set xsection correctly and normalise histos
  // TODO: check that this setting of the xs is really correct
  b->pythia.stat();
//  b->ah->setCrossSection( sigmaTotal, b->pythia.info.sigmaGen()*1.0e9);
  b->ah->finalize();

  // Push histos into block
  // For Rivet 2.X, uncomment
  //b->data = b->ah->getData();
  // For Rivet 3.x, uncomment
  
  auto raos = b->ah->getRivetAOs();
  // This is where we store the AOs to be written.
  vector<YODA::AnalysisObjectPtr> output;

  // First get all multiwight AOs
  output.reserve(raos.size());

  for ( auto rao : raos ) {
    rao.get()->setActiveFinalWeightIdx(0);
    if ( rao->path().find("/TMP/") != string::npos ) continue;
    double sc = 1.0;
    if (rao->hasAnnotation("ScaledBy")) {
      sc = std::stod(rao->annotation("ScaledBy"));
    }
   // fmt::print(stderr, "sc [{}]   \n", sc);
    auto yptr = rao.get()->activeYODAPtr();
    output.push_back(yptr);
  }

  b-> data = output;

  // SM TEST
  // Debug write out --- uncomment to write each block's YODA file
  //b->ah->writeData(std::to_string((1+npc)*(b->state.seed+cp.gid()))+".yoda");

  // This is a bit annoying --- we need to unscale Histo1D and Histo2D beforge the reduction
  // TODO: Figure out whether this is really necessary
  
  for (auto ao : b->data) {
    if (ao->hasAnnotation("ScaledBy")) {
      double sc = std::stod(ao->annotation("ScaledBy"));
      //fmt::print(stderr, "sc data0 [{}]   \n", sc);
    } 
    if (ao->hasAnnotation("ScaledBy")) {
      double sc = std::stod(ao->annotation("ScaledBy"));
      if (ao->type()=="Histo1D") {
         if (sc>0) {
            dynamic_cast<YODA::Histo1D&>(*ao).scaleW(1./sc);
//uncommenting
     //       dynamic_cast<YODA::Histo1D&>(*ao).rmAnnotation("ScaledBy");
            dynamic_cast<YODA::Histo1D&>(*ao).addAnnotation("OriginalScaledBy",1./sc);
         }
      } else if (ao->type()=="Histo2D") {
         if (sc>0) {
            dynamic_cast<YODA::Histo2D&>(*ao).scaleW(1./sc);
//uncommenting
    //        dynamic_cast<YODA::Histo2D&>(*ao).rmAnnotation("ScaledBy");
            dynamic_cast<YODA::Histo2D&>(*ao).addAnnotation("OriginalScaledBy",1./sc);
         }
      }
    } 
    if (ao->hasAnnotation("ScaledBy")) {
      double sc = std::stod(ao->annotation("ScaledBy"));
      //fmt::print(stderr, "sc data1 [{}]   \n", sc);
    } 
  }
  for (auto ao : b->data) {
    if (ao->hasAnnotation("ScaledBy")) {
      double sc = std::stod(ao->annotation("ScaledBy"));
     // fmt::print(stderr, "sc data [{}]   \n", sc);
    }
  } 
}


void write_yoda(Block* b, diy::Master::ProxyWithLink const& cp, int rank, bool verbose)
{
 if (verbose) fmt::print(stderr, "[{}] -- rank {} sees write_yoda \n", cp.gid(), rank);
 fmt::print(stderr, "[{}] -- rank {} sees write_yoda \n", cp.gid(), rank);
  if (rank==0 && cp.gid()==0) {
    for (auto ao : b->buffer) {
      //if (ao->hasAnnotation("OriginalScaledBy")) fmt::print(stderr,"got here\n");
      if (ao->hasAnnotation("OriginalScaledBy"))
      {
        double sc = std::stod(ao->annotation("OriginalScaledBy"));
        //double sc = std::stod(ao->annotation("ScaledBy"));
 fmt::print(stderr, "[{}] -- rank {} sees write_yoda {} \n", cp.gid(), rank, sc);
        if (ao->type()=="Histo1D")
        {
         if (sc>0) {
           dynamic_cast<YODA::Histo1D&>(*ao).scaleW(1./sc);
         }
        }
        else if (ao->type()=="Histo2D")
        {
         if (sc>0) {
           dynamic_cast<YODA::Histo2D&>(*ao).scaleW(1./sc);
         }
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
    int runmode = 0;
    int nEvents=1000;
    size_t seed=1234;
    size_t batch=100000;
    vector<std::string> analyses ={"MC_XS"};
    std::string out_file="diy.yoda";
    std::string pfile="runPythia.cmd";
    std::string in_file="test.h5";
    std::string indir="";
    // get command line arguments
    using namespace opts;
    Options ops(argc, argv);
    bool                      verbose     = ops >> Present('v',
                                                           "verbose",
                                                           "verbose output");
    ops >> Option('t', "thread",    threads,   "Number of threads");
    ops >> Option('n', "nevents",   nEvents,   "Number of events to generate in total");
    ops >> Option('m', "runmode",   runmode,   "Runmode - 0 normal, 1 means performance");
    ops >> Option('b', "nblocks",   nBlocks,   "Number of blocks");
    ops >> Option("batch",   batch,   "Batch size when reading from H5");
    ops >> Option('a', "analysis",  analyses,  "Rivet analyses --- can be given multiple times");
    ops >> Option('o', "output",    out_file,  "Output filename.");
    ops >> Option('s', "seed",      seed,      "The Base seed --- this is incremented for each block.");
    ops >> Option('p', "pfile",     pfile,     "Parameter config file for testing __or__ input file name to look for in directory.");
    ops >> Option('i', "indir",     indir,     "Input directory with hierarchical pfiles");
    ops >> Option('H', "h5input",   in_file,   "HDF5 input file");
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
    bool f_ok;
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
             }
          }
       } 
       else {
          physConfigs.push_back(physConfig);
          out_files.push_back(out_file);
       }
       nConfigs=physConfigs.size();
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
      fmt::print(stderr, "\n    LHE will be read from {}\n", in_file);
      fmt::print(stderr, "\n    Blocks:                  {}\n", blocks);
      fmt::print(stderr, "\n    Physics configurations:  {}\n", nConfigs);
      fmt::print(stderr, "\n    Number of events/config: {}\n", nEvents);
      fmt::print(stderr, "\n    Total number of events:  {}\n", nEvents*nConfigs);
      if (runmode==0) {
         fmt::print(stderr, "\n  S T A N D A R D  \n");
      }
      if (runmode==1) {
         fmt::print(stderr, "\n  P E R F O R M A N C E  \n");
      }
      if (runmode==2) {
         fmt::print(stderr, "\n  P E R F O R M A N C E S A M E \n");
      }
      fmt::print(stderr, "***********************************\n");
    }

    double starttime, endtime;

    PointConfig pc;
    for (size_t ipc=0;ipc<nConfigs;++ipc) {
       if (world.rank()==0) {
          pc = mkRunConfig(blocks, nEvents, seed, physConfigs[ipc], analyses, out_files[ipc]);
       }

       // We need to tell the first block about the new configuration
       master.foreach([world, pc](Block* b, const diy::Master::ProxyWithLink& cp)
                        {set_pc(b, cp, pc); });

       // Broadcast the runconfig to the other blocks
       diy::reduce(master, assigner, comm, &bc_pointconfig<Block>);

       if (verbose) master.foreach([world](Block* b, const diy::Master::ProxyWithLink& cp)
                        {print_block(b, cp, world.rank()); });

       if (runmode==1) {
         std::streambuf *old = cout.rdbuf();
         stringstream ss;
         cout.rdbuf (ss.rdbuf());
         starttime = MPI_Wtime();
         master.foreach([world, verbose, ipc, in_file, nEvents](Block* b, const diy::Master::ProxyWithLink& cp)
                           {process_block_lhe_performance(b, cp, world.size(), world.rank(), verbose, ipc, in_file, nEvents); });
         endtime   = MPI_Wtime(); 
         cout.rdbuf (old);
         printf("[%i] took %.3f seconds\n",world.rank(), endtime-starttime);
       }
       else if (runmode==2) {
         //std::streambuf *old = cout.rdbuf();
         //stringstream ss;
         //cout.rdbuf (ss.rdbuf());
         starttime = MPI_Wtime();
         master.foreach([world, verbose, ipc, in_file, nEvents](Block* b, const diy::Master::ProxyWithLink& cp)
                           {process_block_lhe_performance_same(b, cp, world.size(), world.rank(), verbose, ipc, in_file, nEvents); });
         endtime   = MPI_Wtime(); 
         //cout.rdbuf (old);
         printf("[%i] took %.3f seconds\n",world.rank(), endtime-starttime);
       }
       else {

          master.foreach([world, verbose, ipc, in_file, nEvents, batch](Block* b, const diy::Master::ProxyWithLink& cp)
                           {process_block_lhe(b, cp, world.size(), world.rank(), verbose, ipc, in_file, nEvents, batch); });


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
   }


    return 0;
}
