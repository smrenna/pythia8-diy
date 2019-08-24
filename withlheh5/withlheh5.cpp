
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

#include "PtjTMSdefinitionHooks.h" // Stefan's unlops hook

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

#include "fastjet/ClusterSequence.hh" // This is to quieten fastjet

#include "lheh5.h"


using namespace std;
using namespace Pythia8;
using namespace lheh5;


class LHAupH5 : public Pythia8::LHAup {
  public:
    LHAupH5( HighFive::File* file_, size_t firstEvent, size_t readSize, size_t nTotal, bool verbose=false, bool readNTrials=true, bool noRead=false, size_t batchSize=100000) : _numberRead(0),  _sumW(0), _batchSize(batchSize) {
      file = file_;
      _particle    = file->getGroup("particle");
      _event       = file->getGroup("event");
      _init        = file->getGroup("init");
      _procInfo    = file->getGroup("procInfo");
      try {
         _weightvar   = file->getGroup("weightvariations");
         _hasWeights=true;
      }
       catch (Exception& err) {
        // catch and print any HDF5 error
        std::cerr << "NOTE: Group weightvariations does not exist." << std::endl;
         _hasWeights=false;
      }
      _readSize = readSize;
      _noRead=noRead;
      _firstEvent=firstEvent;
      _readSize=readSize;
      // This reads and holds the information of readSize events, starting from firstEvent
      setInit();
      //if (noRead){
         //lheevents = lheh5::readEvents(_index, _particle, _event, firstEvent, 1);
      DataSet _npLO   =  _procInfo.getDataSet("npLO");
      DataSet _npNLO  =  _procInfo.getDataSet("npNLO");
      _npLO.read(npLO);
      _npNLO.read(npNLO);
      //}
     
      size_t control=0;
      while (control < readSize) {
        vector<size_t> temp;
        temp.push_back(firstEvent + control);
        if (control+_batchSize <= readSize) temp.push_back(_batchSize);
        else temp.push_back(readSize - control);
        control+=_batchSize;
        _work.push_back(temp);
      }
  for (auto w : _work) {
     std::cerr << "offset: " << w[0] << " reads: " << w[1] << " bs: " << _batchSize <<"\n";
  }

      if (readNTrials) {
         // Sum of trials for ALL to be processed events!!!
         DataSet _trials     =  _event.getDataSet("trials");
         std::vector<size_t>    _vtrials;
         _trials    .select({0}, {nTotal}).read(_vtrials);
         _nTrials = std::accumulate(_vtrials.begin(), _vtrials.end(), 0.0);
      }
      // This is for measurement only
      else {
         _nTrials = 100;
      }
    }
    
    // Read and set the info from init and procInfo
    bool setInit() override;// override;
    bool setEvent(int idProc=0) override;// override;

    void setBatch(std::vector<size_t> mywork) {
       _numberRead = 0;
      lheevents = lheh5::readEvents(_particle, _event, mywork[0], mywork[1], npLO, npNLO);
    }
    size_t nTrials() { return _nTrials; }

    int getSize() { return lheevents._vnparticles.size(); }
    int npLO, npNLO;
    vector< vector<size_t> > _work; // Tuple : offset and readsize
  
  private:

    HighFive::File*                         file;
    // Connect with groups
    HighFive::Group                         _index, _particle, _event, _init, _procInfo, _weightvar;
    lheh5::Events2                        lheevents;
    size_t                                     _numberRead;
    size_t                                     _nTrials;
    size_t                                     _batchSize;
    double                                  _sumW;
    bool                                  _noRead;
    size_t                                _readSize;
    bool _hasWeights;
    vector<string>                     _weightnames;
    vector<vector<double>>            _weightvalues;
    vector<double>                    _eventweightvalues;
    size_t _firstEvent;


    // Flag to set particle production scales or not.
    LHAscales scalesNow;


};

bool LHAupH5::setInit()
{
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

   _weightnames.clear();
   _weightvalues.clear();
   if (_hasWeights) {
      for (auto wname : _weightvar.listObjectNames()) {
         _weightnames.push_back(wname);
      }

      std::vector<double> temp;
      temp.reserve(_readSize);
      for (auto w : _weightnames) {

            DataSet _wtemp     =  _weightvar.getDataSet(w);
            _wtemp    .select({_firstEvent}, {_readSize}).read(temp);
            _weightvalues.push_back(temp);
      }
   }





   setBeamA(beamA, energyA, PDFgroupA, PDFsetA);
   setBeamB(beamB, energyB, PDFgroupB, PDFsetB);
   
   int weightingStrategy;
   _init.getDataSet("weightingStrategy").read(weightingStrategy);
   setStrategy(-4);
   
   int numProcesses;
   _init.getDataSet("numProcesses").read(numProcesses);

   // NOTE this is a hack for testing only
   numProcesses = 1;


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

  //setProcess(eHeader.pid,eHeader.weight*(1. / (1e9*_nTrials)),eHeader.scale,eHeader.aqed,eHeader.aqcd);
  setProcess(eHeader.pid,eHeader.weight,eHeader.scale,eHeader.aqed,eHeader.aqcd);

  nupSave    = eHeader.nparticles;
  idprupSave = eHeader.pid;
  xwgtupSave = eHeader.weight;//*(1. / (1e9*_nTrials));
  scalupSave = eHeader.scale; // TODO which scale?
  aqedupSave = eHeader.aqed;
  aqcdupSave = eHeader.aqcd;
  std::vector<lheh5::Particle> particles;
  double scalein = -1.;

   //std::cerr << "INFO: " << infoPtr << "\n";
   //std::cerr << "INFOWCMN: " << infoPtr->weights_compressed_names << "\n";
   infoPtr->weights_compressed_names = &_weightnames;
  _eventweightvalues.clear();;
  for (size_t  i=0; i<_weightnames.size(); ++i) _eventweightvalues.push_back(_weightvalues[i][_numberRead]);

  infoPtr->weights_compressed = &_eventweightvalues;

  // TEMPorary hack for mothers not being set in Sherpa
  if (_noRead) {
     particles = lheevents.mkEvent( 0 );
  }
  else {
     particles = lheevents.mkEvent( _numberRead );
  }

  for (unsigned int ip=0;ip< particles.size(); ++ip) {
     lheh5::Particle part = particles[ip];
     if (ip < 2) {
       addParticle(part.id,part.status,0, 0,part.color1,part.color2,
                   part.px,part.py,part.pz,part.e,part.m,part.lifetime,part.spin,scalein);
     }
     else {
    addParticle(part.id,1,part.mother1,part.mother2,part.color1,part.color2,
                part.px,part.py,part.pz,part.e,part.m,part.lifetime,part.spin,scalein);
     }
  }
    
  // Scale setting
  scalesNow.clear();
  scalesNow.muf   = eHeader.fscale;
  scalesNow.mur   = eHeader.rscale;
  scalesNow.mups  = eHeader.scale;

  infoPtr->scales = &scalesNow;
  
  infoPtr->setEventAttribute("npLO",  std::to_string(eHeader.npLO));
  infoPtr->setEventAttribute("npNLO", std::to_string(eHeader.npNLO));

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
  
  // Create an LHAup object that can access relevant information in pythia.
  LHAupH5* LHAup = new LHAupH5( &file , eventOffset, nMax, nMax, verbose, true, true);

  b->pythia.settings.mode("Beams:frameType", 5);
  // Give the external reader to Pythia
  b->pythia.setLHAupPtr(LHAup);

   // hier unlops zeug
   // Allow to set the number of addtional partons dynamically. TODO here important
  int scheme = ( b->pythia.settings.flag("Merging:doUMEPSTree")
              || b->pythia.settings.flag("Merging:doUMEPSSubt")) ?
              1 :
               ( ( b->pythia.settings.flag("Merging:doUNLOPSTree")
              || b->pythia.settings.flag("Merging:doUNLOPSSubt")
              || b->pythia.settings.flag("Merging:doUNLOPSLoop")
              || b->pythia.settings.flag("Merging:doUNLOPSSubtNLO")) ?
              2 :
              0 );
  HardProcessBookkeeping* hardProcessBookkeepingPtr
    = new HardProcessBookkeeping(scheme);

  b->pythia.setUserHooksPtr(hardProcessBookkeepingPtr);

  double sigmaTotal  = 0.;
  double errorTotal  = 0.;
  double xs = 0.;
  for (int i=0; i < b->pythia.info.nProcessesLHEF(); ++i)
    xs += b->pythia.info.sigmaLHEF(i);
  if (verbose) fmt::print(stderr, "[{}] xs: {}\n", cp.gid(), xs);
  double sigmaSample = 0., errorSample = 0.; // NOTE in Stefan's unlops this is reset for each to be merged multi

  b->pythia.readString("Merging:unlopsTMSdefinition = 1");
  int unlopsType = b->pythia.settings.mode("Merging:unlopsTMSdefinition");

   MergingHooks* ptjTMSdefinitionPtr = (unlopsType<0)
    ? NULL
    : new PtjTMSdefinitionHooks(b->pythia.parm("Merging:TMS"),6.0);
  if (unlopsType >0) b->pythia.setMergingHooksPtr( ptjTMSdefinitionPtr );

  // All configurations done, initialise Pythia
  b->pythia.init();


  if (verbose) fmt::print(stderr, "[{}] starting event loop\n", cp.gid());
  // The event loop
  int nAbort = 5;
  int iAbort = 0;
  if (verbose)  fmt::print(stderr, "[{}] generating  {} events\n", cp.gid(),  LHAup->getSize());
  for (size_t iEvent = 0; iEvent < nMax; ++iEvent) {
    if (verbose) fmt::print(stderr, "[{}] is at event {}\n", cp.gid(), iEvent);
    if (!b->pythia.next()) {
       // Gracefully ignore events with 0 weight
       if (LHAup->weight() == 0) {
          if (verbose) fmt::print(stderr, "[{}] encounters and ignores event {} as it has zero weight\n", cp.gid(), iEvent);
          continue;
       }
       else {
          if (++iAbort < nAbort) continue; // All other errors contribute to the abort counter
       }
      break;
    }
    //LHAup->listEvent();
    if (verbose ) LHAup->listEvent();
  }

  delete hardProcessBookkeepingPtr;
  if (unlopsType>0) delete ptjTMSdefinitionPtr;
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
  LHAupH5* LHAup = new LHAupH5( &file , eventOffset, ev_rank, nEvents, verbose);

  b->pythia.settings.mode("Beams:frameType", 5);
  // Give the external reader to Pythia
  b->pythia.setLHAupPtr(LHAup);

   // hier unlops zeug
   // Allow to set the number of addtional partons dynamically. TODO here important
  HardProcessBookkeeping* hardProcessBookkeeping = NULL;
  int scheme = ( b->pythia.settings.flag("Merging:doUMEPSTree")
              || b->pythia.settings.flag("Merging:doUMEPSSubt")) ?
              1 :
               ( ( b->pythia.settings.flag("Merging:doUNLOPSTree")
              || b->pythia.settings.flag("Merging:doUNLOPSSubt")
              || b->pythia.settings.flag("Merging:doUNLOPSLoop")
              || b->pythia.settings.flag("Merging:doUNLOPSSubtNLO")) ?
              2 :
              0 );
  HardProcessBookkeeping* hardProcessBookkeepingPtr
    = new HardProcessBookkeeping(scheme);

  b->pythia.setUserHooksPtr(hardProcessBookkeepingPtr);

  double sigmaTotal  = 0.;
  double errorTotal  = 0.;
  double xs = 0.;
  for (int i=0; i < b->pythia.info.nProcessesLHEF(); ++i)
    xs += b->pythia.info.sigmaLHEF(i);
  if (verbose) fmt::print(stderr, "[{}] xs: {}\n", cp.gid(), xs);
  double sigmaSample = 0., errorSample = 0.; // NOTE in Stefan's unlops this is reset for each to be merged multi

  b->pythia.readString("Merging:unlopsTMSdefinition = 1");
  int unlopsType = b->pythia.settings.mode("Merging:unlopsTMSdefinition");

   MergingHooks* ptjTMSdefinitionPtr = (unlopsType<0)
    ? NULL
    : new PtjTMSdefinitionHooks(b->pythia.parm("Merging:TMS"),6.0);
  if (unlopsType >0) b->pythia.setMergingHooksPtr( ptjTMSdefinitionPtr );

  // All configurations done, initialise Pythia
  b->pythia.init();


  if (verbose) fmt::print(stderr, "[{}] starting event loop\n", cp.gid());
  // The event loop
  int nAbort = 5;
  int iAbort = 0;
  if (verbose)  fmt::print(stderr, "[{}] generating  {} events\n", cp.gid(),  LHAup->getSize());
  for (size_t iEvent = 0; iEvent < LHAup->getSize(); ++iEvent) {
    if (verbose) fmt::print(stderr, "[{}] is at event {}\n", cp.gid(), iEvent);
    if (!b->pythia.next()) {
       // Gracefully ignore events with 0 weight
       if (LHAup->weight() == 0) {
          if (verbose) fmt::print(stderr, "[{}] encounters and ignores event {} as it has zero weight\n", cp.gid(), iEvent);
          continue;
       }
       else {
          if (++iAbort < nAbort) continue; // All other errors contribute to the abort counter
       }
      break;
    }
    if (verbose ) LHAup->listEvent();
  }

  delete hardProcessBookkeepingPtr;
  if (unlopsType>0) delete ptjTMSdefinitionPtr;
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
  // Create an LHAup object that can access relevant information in pythia.
  //LHAupH5* LHAup = new LHAupH5( &file , eventOffset, ev_rank, nEvents, verbose);
  //// Give the external reader to Pythia
  //b->pythia.setLHAupPtr(LHAup);

   // hier unlops zeug
   // Allow to set the number of addtional partons dynamically. TODO here important
  HardProcessBookkeeping* hardProcessBookkeeping = NULL;
  int scheme = ( b->pythia.settings.flag("Merging:doUMEPSTree")
              || b->pythia.settings.flag("Merging:doUMEPSSubt")) ?
              1 :
               ( ( b->pythia.settings.flag("Merging:doUNLOPSTree")
              || b->pythia.settings.flag("Merging:doUNLOPSSubt")
              || b->pythia.settings.flag("Merging:doUNLOPSLoop")
              || b->pythia.settings.flag("Merging:doUNLOPSSubtNLO")) ?
              2 :
              0 );
  HardProcessBookkeeping* hardProcessBookkeepingPtr
    = new HardProcessBookkeeping(scheme);

  if (scheme!=0) b->pythia.setUserHooksPtr(hardProcessBookkeepingPtr);

  double sigmaTotal  = 0.;
  double errorTotal  = 0.;
  double xs = 0.;
  for (int i=0; i < b->pythia.info.nProcessesLHEF(); ++i)
    xs += b->pythia.info.sigmaLHEF(i);
  if (verbose) fmt::print(stderr, "[{}] xs: {}\n", cp.gid(), xs);
  double sigmaSample = 0., errorSample = 0.; // NOTE in Stefan's unlops this is reset for each to be merged multi

  b->pythia.readString("Merging:unlopsTMSdefinition = 1");
  int unlopsType = b->pythia.settings.mode("Merging:unlopsTMSdefinition");

   MergingHooks* ptjTMSdefinitionPtr = (unlopsType<0)
    ? NULL
    : new PtjTMSdefinitionHooks(b->pythia.parm("Merging:TMS"),6.0);
  if (unlopsType >0) b->pythia.setMergingHooksPtr( ptjTMSdefinitionPtr );


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

  // Create an LHAup object that can access relevant information in pythia.
  LHAupH5* LHAup = new LHAupH5( &file , eventOffset, ev_rank, nEvents, verbose, true, false,  batchSize);
  // Give the external reader to Pythia
  b->pythia.setLHAupPtr(LHAup);
  // All configurations done, initialise Pythia
  b->pythia.init();

  for (auto w : LHAup->_work) {
     fmt::print(stderr, "[{}] next batch\n", cp.gid());
     LHAup->setBatch(w);

     fmt::print(stderr, "[{}] generating  {} events\n", cp.gid(),  LHAup->getSize());
     for (size_t iEvent = 0; iEvent < LHAup->getSize(); ++iEvent) {
       if (verbose) fmt::print(stderr, "[{}] is at event {}\n", cp.gid(), iEvent+w[0]);
       if (!b->pythia.next()) {
          if (verbose) b->pythia.stat();
          if (verbose) LHAup->listEvent();

          // Gracefully ignore events with 0 weight
          if (LHAup->weight() == 0) {
             if (verbose) fmt::print(stderr, "[{}] encounters and ignores event {} as it has zero weight\n", cp.gid(), iEvent+w[0]);
             continue;
          }
          else {
             if (++iAbort < nAbort) continue; // All other errors contribute to the abort counter
          }
         break;
       }
       //if (verbose && iEvent < 2 ) LHAup->listEvent();
       if (verbose ) LHAup->listEvent();
       if (verbose) fmt::print(stderr, "[{}] event weight {} \n", cp.gid(), b->pythia.info.weight());
       //if (verbose) fmt::print(stderr, "[{}] event weight {} {} {}\n", cp.gid(), LHAup->weight(), b->pythia.info.weight(), b->pythia.info.eventWeightLHEF);
         
       // Get event weight(s).
       double evtweight         = b->pythia.info.weight();
       // Additional PDF/alphaS weight for internal merging.
       //evtweight               *= b->pythia.info.mergingWeightNLO() // commented out
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

       double normhepmc = 1. / double(LHAup->nTrials()); //


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
       xsec.set_cross_section( sigmaTotal, b->pythia.info.sigmaErr() );
       hepmcevt->set_cross_section( xsec );
       if (verbose) fmt::print(stderr, "[{}] xsec {} \n", cp.gid(), sigmaTotal);

       // Here more
       try {b->ah->analyze( *hepmcevt ) ;} catch (const std::exception& e)
       {
         if (verbose) fmt::print(stderr, "[{}] exception in analyze: {}\n", cp.gid(), e.what());
       }
       delete hepmcevt;
       // TODO: make 1000 a free parameter a la msg_every
       if (iEvent%100 == 0 && cp.gid()==0) {
          if (b->state.num_events <0 | b->state.num_events > LHAup->getSize()) {
             fmt::print(stderr, "[{}]  {}/{} \n", cp.gid(),  iEvent+w[0], LHAup->getSize());
          }
          else {
             fmt::print(stderr, "[{}]  {}/{} \n", cp.gid(),  iEvent+w[0], b->state.num_events);
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
  b->ah->setCrossSection( sigmaTotal);//b->pythia.info.sigmaGen() * 1.0E9);
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
         if (sc>0) {
            dynamic_cast<YODA::Histo1D&>(*ao).scaleW(1./sc);
            //dynamic_cast<YODA::Histo1D&>(*ao).addAnnotation("OriginalScaledBy", 1./sc);
         }
      }
      else if (ao->type()=="Histo2D")
      {
         if (sc>0) {
            dynamic_cast<YODA::Histo2D&>(*ao).scaleW(1./sc);
            //dynamic_cast<YODA::Histo2D&>(*ao).addAnnotation("OriginalScaledBy", 1./sc);
         }
      }
    }
  }
  // Clean-up TODO important
  delete hardProcessBookkeepingPtr;
  if (unlopsType>0) delete ptjTMSdefinitionPtr;
}


void write_yoda(Block* b, diy::Master::ProxyWithLink const& cp, int rank, bool verbose)
{
 if (verbose) fmt::print(stderr, "[{}] -- rank {} sees write_yoda \n", cp.gid(), rank);
  if (rank==0 && cp.gid()==0) {
    for (auto ao : b->buffer) {
      //if (ao->hasAnnotation("OriginalScaledBy"))
      if (ao->hasAnnotation("ScaledBy"))
      {
        //double sc = std::stod(ao->annotation("OriginalScaledBy"));
        double sc = std::stod(ao->annotation("ScaledBy"));
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
