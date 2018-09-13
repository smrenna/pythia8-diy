
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
#include "Pythia8/LHEF3.h"
#include "Pythia8Plugins/HepMC2.h"
#include "Rivet/AnalysisHandler.hh"
#undef foreach // This line prevents a clash of definitions of rivet's legacy foreach with that of DIY

#include "HepMC/IO_GenEvent.h"


#include "lheh5.h"


using namespace std;
using namespace Pythia8;
using namespace lheh5;


class LHAupH5 : public Pythia8::LHAup {
  public:
    LHAupH5( HighFive::File* file_, size_t firstEvent, size_t readSize, bool verbose=false) : _numberRead(0),  _sumW(0) {
      file = file_;
      _index       = file->getGroup("index");
      _particle    = file->getGroup("particle");
      _event       = file->getGroup("event");
      _init        = file->getGroup("init");
      _procInfo    = file->getGroup("procInfo");
      // This reads and holds the information of readSize events, starting from firstEvent
      //setInit();
      lheevents = lheh5::readEvents(_index, _particle, _event, firstEvent, readSize);
     // Sum of trials for this block
     DataSet _trials     =  file->getDataSet("event/trials");
     std::vector<int>    _vtrials;
     _trials    .select({firstEvent}, {readSize}).read(_vtrials);
     _nTrials = std::accumulate(_vtrials.begin(), _vtrials.end(), 0);
     if (verbose) fmt::print(stderr, "sum trials {}\n", _nTrials);
    }
  
       
    void setScalesFromLHEF(bool b) { setScalesFromLHEF_ = b; }
    
    // Read and set the info from init and procInfo
    bool setInit();// override;
    bool setEvent(int idProc=0);// override;

    int getSize() { return lheevents._vnparticles.size(); }
  
  private:

    HighFive::File*                         file;
    // Connect with groups
    HighFive::Group                         _index, _particle, _event, _init, _procInfo;
    lheh5::Events                        lheevents;
    int                                     _numberRead;
    int                                     _nTrials;
    double                                  _sumW;
  

    // Flag to set particle production scales or not.
    bool setScalesFromLHEF_;
    LHAscales scalesNow;


};

bool LHAupH5::setInit()
{
  /*   if (!runInfo) return false;
       const HEPRUP &heprup = *runInfo->getHEPRUP();*/
   
   int beamA, beamB;
   double energyA, energyB;
   int PDFgroupA, PDFgroupB;
   int PDFsetA, PDFsetB;

   _init.getDataSet("beamA")    .read(beamA);
   _init.getDataSet("energyA")  .read(energyA);
   _init.getDataSet("PDFgroupA").read(PDFgroupA);
   _init.getDataSet("PDFsetA")  .read(PDFsetA);
   
   _init.getDataSet("beamB")    .read(beamB);
   _init.getDataSet("energyB")  .read(energyB);
   _init.getDataSet("PDFgroupB").read(PDFgroupB);
   _init.getDataSet("PDFsetB")  .read(PDFsetB);

   setBeamA(beamA, energyA, PDFgroupA, PDFsetA);
   setBeamB(beamB, energyB, PDFgroupB, PDFsetB);
   
   int weightingStrategy;
   _init.getDataSet("weightingStrategy").read(weightingStrategy);
   setStrategy(weightingStrategy);
   
   int numProcesses;
   _init.getDataSet("numProcesses").read(numProcesses);
   fmt::print(stderr, "numProcesses {}\n", numProcesses);

   // NOTE this is a hack for testing only
   //numProcesses = 1;


   vector<int> procId;        // NOTE: C++17 allows int[numProcesses]
   vector<double> xSection;   // NOTE: C++17 allows double[numProcesses]
   vector<double> error;      // NOTE: C++17 allows double[numProcesses]
   vector<double> unitWeight; // NOTE: C++17 allows double[numProcesses]
   _procInfo.getDataSet("procId").read(procId);
   _procInfo.getDataSet("xSection").read(xSection);
   _procInfo.getDataSet("error").read(error);
   _procInfo.getDataSet("unitWeight").read(unitWeight);
   for (size_t np=0; np<numProcesses;++np) {
     addProcess(procId[np], xSection[np], error[np], unitWeight[np]);
     xSecSumSave += xSection[np];
     xErrSumSave += pow2(error[np]);
   }


  return true;
}

bool LHAupH5::setEvent(int idProc)
{

  lheh5::EventHeader eHeader = lheevents.mkEventHeader( _numberRead );

  setProcess(eHeader.pid,eHeader.weight*(1. / (1e9*_nTrials)),eHeader.scale,eHeader.aqed,eHeader.aqcd);

  //nupSave    = eHeader.nparticles;
  //idprupSave = eHeader.pid;
  //xwgtupSave = eHeader.weight*(1. / (1e9*_nTrials));
  //scalupSave = eHeader.scale; // TODO which scale?
  //aqedupSave = eHeader.aqed;
  //aqcdupSave = eHeader.aqcd;
  // Set directly!  what is scale?   
  //getpro >> nupSave >> idprupSave >> xwgtupSave >> scalupSave
    //>> aqedupSave >> aqcdupSave;

  double scalein = -1.;
  for (auto part : lheevents.mkEvent( _numberRead ) ) {
    fmt::print(stderr, "Adding particle {}\n", part);
    addParticle(part.id,part.status,part.mother1,part.mother2,part.color1,part.color2,
		part.px,part.py,part.pz,part.e,part.m,part.lifetime,part.spin,scalein);
  }
    
  // Scale setting
  scalesNow.clear();
  scalesNow.muf   = eHeader.fscale;
  scalesNow.mur   = eHeader.rscale;
  scalesNow.mups  = eHeader.scale;

  infoPtr->scales = &scalesNow;

  fmt::print(stderr, "numread {}\n", _numberRead);
  // Trials --- TODO ask Stefan again
  //infoPtr->setLHEF3EventInfo( &reader.hepeup.attributes, 0, 0, 0, 0, 0,
       //vector<double>(), "", 1.0);

  //setIdX(this->id1(), this->id2(), this->x1(), this->x2());
  //setPdf(this->id1pdf(), this->id2pdf(), this->x1pdf(), this->x2pdf(),
         //this->scalePDF(), this->pdf1(), this->pdf2(), this->pdfIsSet());
  _numberRead++;


  return true;
}



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


//void process_block(Block* b, diy::Master::ProxyWithLink const& cp, int rank, std::vector<std::string> physConfig, std::vector<std::string> analyses, bool verbose)
void process_block_lhe(Block* b, diy::Master::ProxyWithLink const& cp, int size, int rank,  bool verbose, int npc, string in_file)
{
  // This make rivet only report ERRORs
  // TODO: can we have a global flag to steer verbosity of all moving parts?
  if (!verbose) Rivet::Log::setLevel("Rivet", Rivet::Log::WARNING);

  // Minimise pythia's output
  if (verbose) b->pythia.readString("Print:quiet = off");
  else b->pythia.readString("Print:quiet = on");

  //b->pythia.readString("Print:quiet = on");
  // Tell Pythia that the we are handling the LHE information
  b->pythia.settings.mode("Beams:frameType", 5);

  // Configure pythia with a vector of strings
  for (auto s  : b->state.conf) b->pythia.readString(s);

  // Py8 random seed for this block read from point config
  b->pythia.readString("Random:setSeed = on");
  b->pythia.readString("Random:seed = " + std::to_string(b->state.seed+cp.gid()));


  HighFive::File file(in_file, HighFive::File::ReadOnly);  
  hid_t dspace = H5Dget_space(file.getDataSet("index/start").getId());
  size_t nEvents  =  H5Sget_simple_extent_npoints(dspace);
  size_t ev_rank = floor(nEvents/size);
  // Detect testing
  if (ev_rank > 1e6) ev_rank = 1e6;
  size_t eventOffset = rank*ev_rank;
  if (rank == size-1 && size>1) {
     ev_rank = nEvents-eventOffset;
  }

  // TODO: can't hurt to test whether this logic ^^^ is correct

  if (verbose) fmt::print(stderr, "[{}] reads {} events starting at {}\n", cp.gid(), ev_rank, eventOffset);
  // Create an LHAup object that can access relevant information in pythia.
  LHAupH5* LHAup = new LHAupH5( &file , eventOffset, ev_rank, verbose);

  if (verbose) LHAup->listInit();
  if (verbose) fmt::print(stderr, "[{}] read {} events\n", cp.gid(), LHAup->getSize());

  // Give the external reader to Pythia
  b->pythia.setLHAupPtr(LHAup);


  // All configurations done, initialise Pythia
  b->pythia.init();

  // Delete the AnalysisHandlerPtr to ensure there is no memory
  if (b->ah)
  {
    delete b->ah;
  }
  b->ah = new Rivet::AnalysisHandler;
  b->ah->setIgnoreBeams();

  // Add all anlyses to rivet
  // TODO: we may want to feed the "ignore beams" switch as well
  for (auto a : b->state.analyses) {
     b->ah->addAnalysis(a);
     if (verbose) fmt::print(stderr, "[{}] add  ######## analysis {}\n", cp.gid(), a);
  }

  if (verbose) fmt::print(stderr, "[{}] starting event loop\n", cp.gid());
  // The event loop
  int nAbort = 5;
  int iAbort = 0;
  if (verbose) fmt::print(stderr, "[{}] generating {} events\n", cp.gid(),  LHAup->getSize());
  for (int iEvent = 0; iEvent < LHAup->getSize(); ++iEvent) {
    if (!b->pythia.next()) {
      if (++iAbort < nAbort) continue; // TODO investigate influenec of apbort on sum trials
      break;
    }
    if (verbose) fmt::print(stderr, "[{}] event weight {} \n", cp.gid(), b->pythia.info.weight());
    //if (verbose) fmt::print(stderr, "[{}] event weight {} {} {}\n", cp.gid(), LHAup->weight(), b->pythia.info.weight(), b->pythia.info.eventWeightLHEF);
    if (verbose && iEvent < 5 ) LHAup->listEvent();
    HepMC::GenEvent* hepmcevt = new HepMC::GenEvent();
    b->ToHepMC.fill_next_event( b->pythia, hepmcevt );

    // Here more
    try {b->ah->analyze( *hepmcevt ) ;} catch (const std::exception& e)
    {
      if (verbose) fmt::print(stderr, "[{}] exception in analyze: {}\n", cp.gid(), e.what());
    }
    delete hepmcevt;
    // TODO: make 1000 a free parameter a la msg_every
    if (iEvent%1000 == 0 && cp.gid()==0) {
       if (b->state.num_events <0 | b->state.num_events > LHAup->getSize()) {
          fmt::print(stderr, "[{}]  {}/{} \n", cp.gid(),  iEvent, LHAup->getSize());
       }
       else {
          fmt::print(stderr, "[{}]  {}/{} \n", cp.gid(),  iEvent, b->state.num_events);
       }
    }

    if (iEvent >= b->state.num_events && b->state.num_events>=0) {
       fmt::print(stderr, "[{}] exiting event loop after: {}/{}\n", cp.gid(), iEvent, b->state.num_events);
       break;
    }
  }

  //fmt::print(stderr, "[{}] xs after: {}\n", cp.gid(), b->pythia.info.sigmaGen());
  // Event loop is done, set xsection correctly and normalise histos
  // TODO: check that this setting of the xs is really correct
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
    ops >> Option('b', "nblocks",   nBlocks,   "Number of blocks");
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
      fmt::print(stderr, "***********************************\n");
    }

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

       master.foreach([world, verbose, ipc, in_file](Block* b, const diy::Master::ProxyWithLink& cp)
                        {process_block_lhe(b, cp, world.size(), world.rank(), verbose, ipc, in_file); });


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
