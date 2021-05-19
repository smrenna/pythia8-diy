// LHAHDF5.h is a part of the PYTHIA event generator.
// Copyright (C) 2021 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.
// Author: Christian Preuss, January 2021.

#ifndef Pythia8_LHAHDF5_H
#define Pythia8_LHAHDF5_H

// Interface includes.
#include "Pythia8Plugins/LHEH5.h"

// Generator includes.
#include "Pythia8/Pythia.h"

namespace Pythia8 {

//==========================================================================

// HDF5 file reader. Converts to Pythia-internal events by acting as
// replacement Les Houches Event reader.

class LHAupH5 : public Pythia8::LHAup {

 public:

  LHAupH5(HighFive::File* h5fileIn, size_t firstEventIn, size_t readSizeIn,
    string versionIn = "") :
    nTrials(0), nRead(0), isOldFormat(h5fileIn->exist("/index")),
    particleSav(h5fileIn->getGroup("particle")),
    eventSav(h5fileIn->getGroup("event")),
    initSav(h5fileIn->getGroup("init")),
    procInfoSav(h5fileIn->getGroup("procInfo")),
    readSizeSav(readSizeIn), firstEventSav(firstEventIn), fileSav(h5fileIn),
    batchSizeSav(10000000) {

    // A dataset that should always exist
    DataSet v = eventSav.getDataSet("nparticles");

    hid_t dset = v.getId();

    hid_t dcpl = H5Dget_create_plist (dset);

    H5D_layout_t layout = H5Pget_layout (dcpl);

    is_virtual = (H5D_VIRTUAL == layout) ? true : false;

    if (is_virtual)  {
        printf(" Dataset has a virtual layout \n");
    }
    else
        printf(" No layout found \n");

    if( is_virtual ) {

    size_t num_map;
    ssize_t len;
    herr_t status = H5Pget_virtual_count (dcpl, &num_map);
    printf(" Number of mappings is %lu\n", (unsigned long)num_map);

    hid_t vspace;
    char         *filename;                  
    char         *dsetname;
    hsize_t      start_out[2],
                 stride_out[2],
                 count_out[2],
                 block_out[2];

    for (int i = 0; i < (int)num_map; i++) {   
        printf(" Mapping %d \n", i);
        printf("         Selection in the virtual dataset ");
        vspace = H5Pget_virtual_vspace (dcpl, (size_t)i);

        if (H5Sget_select_type(vspace) == H5S_SEL_HYPERSLABS) { 
            if (H5Sis_regular_hyperslab(vspace)) {
                status = H5Sget_regular_hyperslab (vspace, start_out, stride_out, count_out, block_out);
                printf("\n         start  = [%llu, %llu] \n", (unsigned long long)start_out[0], (unsigned long long)start_out[1]);
                printf("         stride = [%llu, %llu] \n", (unsigned long long)stride_out[0], (unsigned long long)stride_out[1]);
                printf("         count  = [%lli, %llu] \n", (signed long long)count_out[0], (unsigned long long)count_out[1]);
                printf("         block  = [%llu, %llu] \n", (unsigned long long)block_out[0], (unsigned long long)block_out[1]);
            }
     
        }
        len = H5Pget_virtual_filename (dcpl, (size_t)i, NULL, 0);
        filename = (char *)malloc((size_t)len*sizeof(char)+1);
	// This vector will contain ths start points for each of the datasets
	vstart.push_back( start_out[0] );

        H5Pget_virtual_filename (dcpl, (size_t)i, filename, len+1);
        printf("         Source filename %s\n", filename);
        len = H5Pget_virtual_dsetname (dcpl, (size_t)i, NULL, 0);
        dsetname = (char *)malloc((size_t)len*sizeof(char)+1);
        H5Pget_virtual_dsetname (dcpl, (size_t)i, dsetname, len+1);
        printf("         Source dataset name %s\n", dsetname);

        /* printf("         Selection in the source dataset "); */
        /* hid_t src_space = H5Pget_virtual_srcspace (dcpl, (size_t)i); */
        /* if(H5Sget_select_type(src_space) == H5S_SEL_ALL) { */
	/*   hsize_t *dims, *maxdims; */
	/*   int nDim = H5Sget_simple_extent_dims(src_space,dims,maxdims); */
	/*   std::cout << " nDim " << nDim << " " << dims << " " << maxdims << std::endl; */
        /* } */

        status = H5Sclose(vspace);
	//        status = H5Sclose(src_space);
        free(filename);
        free(dsetname);
    }

    } else {
      vstart.push_back(0L);
    }

    // Check if version number exists in event file.
    valid = false;
    if (h5fileIn->exist("/init/version")) {
      DataSet version = initSav.getDataSet("version");
      vector<int> versionNo;
      version.read(versionNo);
      versionSav = to_string(versionNo[0]) + "." + to_string(versionNo[1])
        + "." + to_string(versionNo[2]);
      isOldFormat = (versionSav=="0.1.0");
    } else if (isOldFormat) {
      // Based on the existence of the index group, we can guess this is 0.1.0.
      versionSav = "0.1.0";
    } else {
      // Have to rely on input if we could not guess the version.
      versionSav = versionIn;
      if (versionSav=="") {
        cout << " LHAupH5 error - cannot determine version number.\n";
        return;
      }
    }
    cout << " LHAupH5 version " << versionSav << ".\n";

    // Check if version supports multiweights. (Starting from 1.0.0.)
    hasMultiWts = !(versionSav=="0.1.0" || versionSav=="0.2.0");
    cout << " LHAupH5 file format "
         << (hasMultiWts?"supports":"does not support")
         << " multi-weights." << endl;

    hid_t dspace;
    if( !isOldFormat ) {
      DataSet npLO  = procInfoSav.getDataSet("npLO");
      DataSet npNLO = procInfoSav.getDataSet("npNLO");
      npLO.read(npLOSav);
      npNLO.read(npNLOSav);
      dspace = H5Dget_space(fileSav->getDataSet("event/start").getId());
    } else {
      indexSav      = fileSav->getGroup("index");
      dspace = H5Dget_space(fileSav->getDataSet("index/start").getId());
    }
    // Check if size of read is compatible.
    if (readSizeIn > H5Sget_simple_extent_npoints(dspace) ) {
      cout << "H5 size request incompatible with file.\n";
      return;
    }

    // Check for multiweights.
    if ( fileSav->exist("event/weight") && hasMultiWts ) {
      DataSet weights = eventSav.getDataSet("weight");
      auto attr_keys = weights.listAttributeNames();
      Attribute a = weights.getAttribute(attr_keys[0]);
      a.read(weightsNames);
    }

    // This reads and holds the information of readSize events,
    // starting from firstEvent.
    if (!isOldFormat) {
      lheEvts2Sav = LHEH5::readEvents2(particleSav, eventSav,
	firstEventSav, readSizeSav, npLOSav, npNLOSav, hasMultiWts,
	is_virtual);
      particlesSizeSav = lheEvts2Sav._vnparticles.size();
    } else {
      lheEvtsSav = LHEH5::readEvents(indexSav, particleSav, eventSav,
        firstEventSav, readSizeSav);
      particlesSizeSav = lheEvtsSav._vnparticles.size();
    }

    size_t control=0;
    while (control < readSizeSav) {
      vector<size_t> temp;
      temp.push_back(firstEventSav + control);
      if (control+batchSizeSav <= readSizeSav) temp.push_back(batchSizeSav);
      else temp.push_back(readSizeSav - control);
      control+=batchSizeSav;
      _work.push_back(temp);
    }
    valid = true;

  }

  // Read and set the info from init and procInfo.
  bool setInit() override;
  bool setEvent(int idProc=0) override;
  void forceStrategy(int strategyIn) {setStrategy(strategyIn);}
  void setBatchSize( size_t batchSizeIn ) { batchSizeSav = batchSizeIn; }
  size_t getTrials() {return nTrials;}
  size_t getSize() {return particlesSizeSav;}
  void setBatch(std::vector<size_t> mywork) {
    nRead = 0;
    // This reads and holds the information of readSize events,
    // starting from firstEvent.
    std::cout << " in set Batch ... " << isOldFormat << std::endl;
    if (!isOldFormat) {
      lheEvts2Sav = LHEH5::readEvents2(particleSav, eventSav,
        firstEventSav, readSizeSav, npLOSav, npNLOSav, 
        hasMultiWts, is_virtual);
      particlesSizeSav = lheEvts2Sav._vnparticles.size();

    } else {
      lheEvtsSav = LHEH5::readEvents(indexSav, particleSav, eventSav,
        firstEventSav, readSizeSav);
      particlesSizeSav = lheEvtsSav._vnparticles.size();
    }
  }
  vector< vector<size_t> > _work; // Tuple : offset and readsize

 private:

  // HDF5 file.
  HighFive::File* fileSav{};

  // HDF5 Groups.
  HighFive::Group particleSav, eventSav, initSav, procInfoSav, indexSav;

  // Events from HDF5 file.
  LHEH5::Events  lheEvtsSav;
  LHEH5::Events2 lheEvts2Sav;

  // Info for reader.
  size_t         readSizeSav, firstEventSav, particlesSizeSav,
                 batchSizeSav;
  int            nTrials, nRead;
  int            npLOSav, npNLOSav;
  bool           valid, isOldFormat, hasMultiWts, is_virtual;
  string         versionSav;

  // Vector contain start points for virtual datasets
  vector<unsigned long long> vstart;

  // Multiweight vector. Reset each event.
  vector<double> weightsSav;
  vector<string> weightsNames;

  // Particle production scales.
  LHAscales scalesNow;

};

//--------------------------------------------------------------------------

// HDF5 file reader. Converts to Pythia-internal events by acting as
// replacement Les Houches Event reader.

bool LHAupH5::setInit() {
  int beamA, beamB;
  double energyA, energyB;
  int PDFgroupA, PDFgroupB;
  int PDFsetA, PDFsetB;

  if (!valid) return false;
  initSav.getDataSet("beamA").read(beamA);
  initSav.getDataSet("energyA").read(energyA);
  initSav.getDataSet("PDFgroupA").read(PDFgroupA);
  initSav.getDataSet("PDFsetA").read(PDFsetA);

  initSav.getDataSet("beamB").read(beamB);
  initSav.getDataSet("energyB").read(energyB);
  initSav.getDataSet("PDFgroupB").read(PDFgroupB);
  initSav.getDataSet("PDFsetB").read(PDFsetB);

  setBeamA(beamA, energyA, PDFgroupA, PDFsetA);
  setBeamB(beamB, energyB, PDFgroupB, PDFsetB);

  int weightingStrategy = 3;
  initSav.getDataSet("weightingStrategy").read(weightingStrategy);
  setStrategy(weightingStrategy);

  int nProcesses = 1;
  initSav.getDataSet("numProcesses").read(nProcesses);

  vector<int> procId;
  vector<double> xSection, error, unitWeight;
  procInfoSav.getDataSet("procId").read(procId);
  procInfoSav.getDataSet("xSection").read(xSection);
  procInfoSav.getDataSet("error").read(error);
  procInfoSav.getDataSet("unitWeight").read(unitWeight);
  infoPtr->sigmaLHEFSave.resize(0);
  for (int iProc(0); iProc<nProcesses; ++iProc) {
    addProcess(procId[iProc], xSection[iProc], error[iProc],
               unitWeight[iProc]);
    infoPtr->sigmaLHEFSave.push_back(xSection[iProc]);
  }
  return true;

}

//--------------------------------------------------------------------------

// Read an event.

bool LHAupH5::setEvent(int idProc) {

  // Equivalent of end of file.
  if (!valid) return false;
  if (nRead >= readSizeSav) return false;
  LHEH5::EventHeader evtHeader = !isOldFormat ?
    lheEvts2Sav.mkEventHeader(nRead) : lheEvtsSav.mkEventHeader(nRead);
  weightsSav = evtHeader.weights;
  nTrials += evtHeader.trials;
  // Skip zero-weight events, but add trials.
  while (weightsSav[0] == 0. && nRead < readSizeSav - 1) {
    ++nRead;
    evtHeader  = !isOldFormat ?
      lheEvts2Sav.mkEventHeader(nRead) : lheEvtsSav.mkEventHeader(nRead);
    weightsSav = evtHeader.weights;
    nTrials   += evtHeader.trials;
  }
  // Communicate event weight to Info.
  infoPtr->weightContainerPtr->setWeightNominal( weightsSav[0] );
  infoPtr->weightContainerPtr->weightsLHEF.bookVectors(
    weightsSav, weightsNames );

  xwgtupSave = weightsSav[0];
  idprupSave = evtHeader.pid;
  scalupSave = evtHeader.scale;
  aqedupSave = evtHeader.aqed;
  aqcdupSave = evtHeader.aqcd;
  setProcess(idprupSave, xwgtupSave, scalupSave, aqedupSave, aqcdupSave);
  double scalein = -1.;
  vector<LHEH5::Particle> particles;
  particles = !isOldFormat ? lheEvts2Sav.mkEvent(nRead)
    : lheEvtsSav.mkEvent(nRead);

  // Set particles.
  int nPtcls = 0;
  for (int iPtcl(0); iPtcl<int(particles.size()); ++iPtcl) {
    LHEH5::Particle ptcl = particles.at(iPtcl);
    if (ptcl.id == 0) continue;
    nPtcls++;
    if (iPtcl < 2) addParticle(ptcl.id, ptcl.status, 0, 0,
        ptcl.color1, ptcl.color2, ptcl.px, ptcl.py, ptcl.pz, ptcl.e, ptcl.m,
        ptcl.lifetime, 9*ptcl.spin, scalein);
    else addParticle(ptcl.id, ptcl.status, ptcl.mother1, ptcl.mother2,
        ptcl.color1, ptcl.color2, ptcl.px, ptcl.py, ptcl.pz, ptcl.e, ptcl.m,
        ptcl.lifetime, 9*ptcl.spin, scalein);
  }
  nupSave = nPtcls;

  // Scale setting
  scalesNow.clear();
  scalesNow.muf   = evtHeader.fscale;
  scalesNow.mur   = evtHeader.rscale;
  scalesNow.mups  = evtHeader.scale;
  infoPtr->scales = &scalesNow;
  ++nRead;
  return true;

}

//==========================================================================

} // end namespace Pythia8

#endif // Pythia8_LHAHDF5_H
