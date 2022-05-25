
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
#include "Pythia8Plugins/ResonanceDecayFilterHook.h"
#include "Rivet/AnalysisHandler.hh"
//#include "Rivet/Analysis.hh"
#undef foreach // This line prevents a clash of definitions of rivet's legacy foreach with that of DIY

#include "HepMC/IO_GenEvent.h"

#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5File.hpp>

#include "H5FDmpio.h"


using namespace std;
using namespace Pythia8;
using namespace HighFive;
using namespace Rivet;
#include "opts.h"

using namespace std;

typedef diy::DiscreteBounds Bounds;
typedef diy::RegularGridLink RCLink;

typedef GenericBlock<Bounds, PointConfig, AnalysisObjects> Block;
typedef ConfigBlockAdder<Bounds, RCLink, Block> AddBlock;

void treeVertex(HepMC::GenEvent& event) {
    // Iterate over all vertices to find PS vertices
    int vtx_id = -1;

    for (HepMC::GenEvent::vertex_const_iterator vit=event.vertices_begin();
       vit!=event.vertices_end(); ++vit) {
      // Is this a PS Vertex?
      if ((*vit)->id()==4) {
        std::vector<HepMC::GenParticle*> remove;
        //// Loop over outgoing particles
        for (HepMC::GenVertex::particles_out_const_iterator pout
               =(*vit)->particles_out_const_begin();
             pout!=(*vit)->particles_out_const_end(); ++pout) {
            if ( (*pout)->end_vertex() ) {
              vtx_id = (*pout)->end_vertex()->id(); //
              // Disconnect outgoing particle from end/decay vertex of type (1,2,3)
              if (vtx_id==1 || vtx_id==2 || vtx_id==3 )
                  remove.push_back((*pout));
            }
        }
        // Loop over incoming particles
        for (HepMC::GenVertex::particles_in_const_iterator pin
               =(*vit)->particles_in_const_begin();
             pin!=(*vit)->particles_in_const_end(); ++pin) {
          vtx_id = (*pin)->production_vertex()->id();
          // Disconnect incoming particle from production vertex of type (1,2,3)
          if (vtx_id==1 || vtx_id==2 || vtx_id==3 )
                  remove.push_back((*pin));
        }
        // Iterate over Genparticle pointers to remove from current vertex (*vit)
        for (unsigned int nrem=0;nrem<remove.size();++nrem) {
            (*vit)->remove_particle(remove[nrem]);
        }
      } // Close if statement (vertex id==4)
    } // Close loop over vertices
}

void print_block(Block* b, const diy::Master::ProxyWithLink& cp)
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


void process_block(Block* b, diy::Master::ProxyWithLink const& cp, bool verbose)
{
  // This make rivet only report ERRORs
  // TODO: can we have a global flag to steer verbosity of all moving parts?
  if (!verbose) Rivet::Log::setLevel("Rivet", Rivet::Log::ERROR);

  // Explicit desctruction and recreation of pythia --- this is important when running multiple
  // physics configs!!! https://stackoverflow.com/questions/1124634/call-destructor-and-then-constructor-resetting-an-object
  (&b->pythia)->~Pythia();
  new (&b->pythia) Pythia();


  // Minimise pythia's output
  if (verbose) b->pythia.readString("Print:quiet = off");
  else b->pythia.readString("Print:quiet = on");


  // UserHooks wrapper
  auto myUserHooks = make_shared<ResonanceDecayFilterHook>(b->pythia.settings);
  b->pythia.setUserHooksPtr(myUserHooks);

  // Configure pythia with a vector of strings
  for (auto s  : b->state.conf) b->pythia.readString(s);

  // Py8 random seed for this block read from point config
  b->pythia.readString("Random:setSeed = on");
  b->pythia.readString("Random:seed = " + std::to_string(b->state.seed+cp.gid()));


  // All configurations done, initialise Pythia
//  b->pythia.initPtrs(); // TODO --- is this really necessary here?
  b->pythia.init();
  //fmt::print(stderr, "[{}] ### W {} - T {} - S {}  \n", cp.gid(), b->pythia.info.weight(), b->pythia.info.nTried(), b->pythia.info.sigmaGen() * 1.0E9);

  // Delete the AnalysisHandlerPtr to ensure there is no memory
  if (b->ah)
  {
    delete b->ah;
  }
  b->ah = new Rivet::AnalysisHandler;

  // Add all anlyses to rivet
  // TODO: we may want to feed the "ignore beams" switch as well
  for (auto a : b->state.analyses) b->ah->addAnalysis(a);

  // The event loop
  int nAbort = 5;
  int iAbort = 0;
  if (verbose) fmt::print(stderr, "[{}] generating {} events\n", cp.gid(),  b->state.num_events);
  for (unsigned int iEvent = 0; iEvent < b->state.num_events; ++iEvent) {
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

  // Event loop is done, set xsection correctly and normalise histos
  //SM b->ah->setCrossSection(b->pythia.info.sigmaGen() * 1.0E9);
  b->ah->setCrossSection(b->pythia.info.sigmaGen() * 1.0E9,b->pythia.info.sigmaErr() * 1.0E9);
  b->ah->finalize();

  b-> data = b->ah->getYodaAOs(true);

}


void write_yoda(Block* b, diy::Master::ProxyWithLink const& cp, int rank, bool verbose)
{
 if (verbose) fmt::print(stderr, "[{}] -- rank {} sees write_yoda \n", cp.gid(), rank);
  if (rank==0 && cp.gid()==0) {
  
    if (verbose) fmt::print(stderr, "[{}] -- writing to file {}  \n", cp.gid(), b->state.f_out);
    YODA::WriterYODA::write(b->state.f_out, b->buffer);
  }
}

//void process_block_hepmc(Block* b, diy::Master::ProxyWithLink const& cp, bool verbose)
//{
  //// Explicit desctruction and recreation of pythia --- this is important when running multiple
  //// physics configs!!! https://stackoverflow.com/questions/1124634/call-destructor-and-then-constructor-resetting-an-object
  //(&b->pythia)->~Pythia();
  //new (&b->pythia) Pythia();

  //// Minimise pythia's output
  //b->pythia.readString("Print:quiet = on");

  //// Configure pythia with a vector of strings
  //for (auto s  : b->state.conf) b->pythia.readString(s);

  //// Py8 random seed for this block read from point config
  //b->pythia.readString("Random:setSeed = on");
  //b->pythia.readString("Random:seed = " + std::to_string(b->state.seed+cp.gid()));

  //// All configurations done, initialise Pythia
  //b->pythia.init();
  //// Delete the AnalysisHandlerPtr to ensure there is no memory
  //if (b->ah) { delete b->ah; }


  //// The event loop
  //b->events.reserve(b->state.num_events);
  //int nAbort = 5;
  //int iAbort = 0;
  //if (verbose) fmt::print(stderr, "[{}] generating {} events\n", cp.gid(),  b->state.num_events);
  //for (unsigned int iEvent = 0; iEvent < b->state.num_events; ++iEvent) {
    //if (!b->pythia.next()) {
      //if (++iAbort < nAbort) continue;
      //break;
    //}
    //HepMC::GenEvent* hepmcevt = new HepMC::GenEvent();
    //b->ToHepMC.fill_next_event( b->pythia, hepmcevt );
    //b->events.push_back(hepmcevt);
    //if (iEvent%1000 == 0 && cp.gid()==0) fmt::print(stderr, "[{}]  {}/{} \n", cp.gid(),  iEvent, b->state.num_events);;
  //}

  //for (auto e : b->events) {
     ////treeVertex(&e);
     //e->print();
  //}
  //// Push histos into block
  ////b->data = b->ah->getData();

  //// Debug write out --- uncomment to write each block's YODA file
  ////b->ah->writeData(std::to_string((1+npc)*(b->state.seed+cp.gid()))+".yoda");


  //// This is a bit annoying --- we need to unscale Histo1D and Histo2D beforge the reduction
  //// TODO: Figure out whether this is really necessary
//}




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
    

    File *hdfoutput_file;
    size_t nBlocks = 0;
    int threads = 1;
    int runnum = -1;
    int nEvents=1000;
    size_t seed=1234;
    vector<std::string> analyses;
    std::string out_file="diy.yoda";
    std::string pfile="runPythia.cmd";
    std::string indir="";
    // get command line arguments
    using namespace opts;
    Options ops(argc, argv);
    bool                      verbose     = ops >> Present('v',
                                                           "verbose",
                                                           "verbose output");
    bool writehepmc   = ops >> Present('w', "writehepmc", "write hepmc events");
    
    ops >> Option('t', "thread",    threads,   "Number of threads");
    ops >> Option('n', "nevents",   nEvents,   "Number of events to generate in total");
    ops >> Option('b', "nblocks",   nBlocks,   "Number of blocks");
    ops >> Option('a', "analysis",  analyses,  "Rivet analyses --- can be given multiple times");
    ops >> Option('o', "output",    out_file,  "Output filename.");
    ops >> Option('s', "seed",      seed,      "The Base seed --- this is incremented for each block.");
    ops >> Option('p', "pfile",     pfile,     "Parameter config file for testing __or__ input file name to look for in directory.");
    ops >> Option('i', "indir",     indir,     "Input directory with hierarchical pfiles");
    ops >> Option('r', "runnum",    runnum,    "Use only runs ending runnum");
    if (ops >> Present('h', "help", "Show help"))
    {
        std::cout << "Usage: " << argv[0] << " [OPTIONS]\n";
        std::cout << ops;
        return 1;
    }

    bool eventOverride = (nEvents != 1000) ? true : false;
    bool  seedOverride = (seed != 1234) ? true: false;
    
    cout << "override = " << eventOverride << " " << seedOverride << endl;

    int nConfigs;

    int batchsize=10;


    //int data[1][2] = {{world.rank(), world.rank()}};
    std::vector<int> data(50, world.rank());

/*    if (writehepmc) {
             hdfoutput_file = new File(out_file,
			File::ReadWrite|File::Create|File::Truncate,
			MPIOFileDriver(MPI_COMM_WORLD,MPI_INFO_NULL));
             Group g_event = hdfoutput_file->createGroup("event");

             size_t size(world.size());
             std::vector<size_t> min(1,size);
             min.front()*=nEvents;
             std::vector<size_t> max(min);
             DataSetCreateProps props;
             max.front()=DataSpace::UNLIMITED;
	     props.add(Chunking(std::vector<hsize_t>(1,batchsize*size)));
             DataSet ds_nparticles = g_event.createDataSet<int>("nparticles",DataSpace(min,max),props);
             ds_nparticles.select({std::size_t(world.rank()), 0}).write(data);
             std::cout<<"kkk\n";
             //ds_nparticles.write(data)();
             //hdfoutput_file->close();
             //

    } */

    std::vector<int> iSeeds;
    std::vector<int> numEvents;
    std::vector<std::string> physConfig;
    std::vector<std::vector<std::string> > physConfigs;
    std::vector<std::string> out_files;
    bool f_ok;
    // copy of Pythia for reading parameters
    Pythia pyDumb;
    if( world.rank()==0 ) {
       // Program logic: check whether a single parameter file has been given or
       // a directory.
       f_ok = readConfig(pfile, physConfig, verbose);
       if( f_ok ) {
	 int iNumSave  = pyDumb.settings.mode("Main:numberOfEvents");
         int iSeedSave = pyDumb.settings.mode("Random:seed");
         bool iSetSave = pyDumb.settings.flag("Random:setSeed");
	 pyDumb.readFile(pfile);
	 int iNum  = pyDumb.settings.mode("Main:numberOfEvents");
	 int iSeed = pyDumb.settings.mode("Random:seed");
	 bool iSet = pyDumb.settings.flag("Random:setSeed");
	 cout << iNum << " " << iSeed << " " << iSet << endl;
	 // Reset to default settings
	 pyDumb.settings.mode("Main:numberOfEvents",iNumSave);
	 pyDumb.settings.mode("Random:seed",iSeedSave);
	 pyDumb.settings.flag("Random:setSeed",iSetSave);
	 // Work out the logic for using the seed/events in the params
         if( seedOverride || !iSet ) {
           iSeeds.push_back( seed );
         } else {
           iSeeds.push_back( iSet );
         }
         if( eventOverride ) {
           numEvents.push_back( nEvents );
         } else {
           numEvents.push_back( iNum );
         }
       }
       if (!f_ok) {
          // Use glob to look for files in directory
          std::string pattern = indir + "/*/" + pfile;
          if (runnum >=0) {
             pattern =  indir + "/*" + std::to_string(runnum) + "/" + pfile;
          } 
          for (auto f : glob(pattern)) {
             physConfig.clear();
             bool this_ok = readConfig(f, physConfig, verbose);
             if (this_ok) {
	       int iNumSave  = pyDumb.settings.mode("Main:numberOfEvents");
	       int iSeedSave = pyDumb.settings.mode("Random:seed");
	       bool iSetSave = pyDumb.settings.flag("Random:setSeed");
	       pyDumb.readFile(f);
	       int iNum  = pyDumb.settings.mode("Main:numberOfEvents");
	       int iSeed = pyDumb.settings.mode("Random:seed");
	       bool iSet = pyDumb.settings.flag("Random:setSeed");
	       cout << iNum << " " << iSeed << " " << iSet << endl;
	       // Reset to default settings
	       pyDumb.settings.mode("Main:numberOfEvents",iNumSave);
	       pyDumb.settings.mode("Random:seed",iSeedSave);
	       pyDumb.settings.flag("Random:setSeed",iSetSave);
	       // Work out the logic for using the seed/events in the params
               if( seedOverride || !iSet ) {
                 iSeeds.push_back( seed );
               } else {
                 iSeeds.push_back( iSet );
               }
               if( eventOverride ) {
                 numEvents.push_back( nEvents );
               } else {
                 numEvents.push_back( iNum );
               }               
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
      fmt::print(stderr, "\n    Blocks:                  {}\n", blocks);
      fmt::print(stderr, "\n    Physics configurations:  {}\n", nConfigs);
      fmt::print(stderr, "\n    Number of events/config: {}\n", nEvents);
      fmt::print(stderr, "\n    Total number of events:  {}\n", nEvents*nConfigs);
      fmt::print(stderr, "***********************************\n");
    }

    PointConfig pc;
    for (int ipc=0;ipc<nConfigs;++ipc) {
       if (world.rank()==0) {
          pc = mkRunConfig(blocks, numEvents[ipc], iSeeds[ipc], physConfigs[ipc],
			   analyses, out_files[ipc]);
       }

       // We need to tell the first block about the new configuration
       master.foreach([world, pc](Block* b, const diy::Master::ProxyWithLink& cp)
                        {set_pc(b, cp, pc); });

       // Broadcast the runconfig to the other blocks
       diy::reduce(master, assigner, comm, &bc_pointconfig<Block>);

       if (verbose) master.foreach([world](Block* b, const diy::Master::ProxyWithLink& cp)
                        {print_block(b, cp); });

       //if (writehepmc) {
          //master.foreach([world, verbose, ipc](Block* b, const diy::Master::ProxyWithLink& cp)
                           //{process_block_hepmc(b, cp, verbose); });
       //}

       //else {

          master.foreach([world, verbose, ipc](Block* b, const diy::Master::ProxyWithLink& cp)
                           {process_block(b, cp, verbose); });


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
       //}
   }


    return 0;
}
