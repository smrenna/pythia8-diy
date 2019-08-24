// PtjTMSdefinitionHooks.h is intended as header for main programs of
// the PYTHIA event generator.

// This header is written by Stefan Prestel.
// It includes bookkeeping of hard process states read from input files
// and intended for CKKW-L/UMEPS/UNLOPS merging in PYTHIA 8

#ifndef Pythia8_PtjTMSdefinitionHooks_H
#define Pythia8_PtjTMSdefinitionHooks_H

#include "Pythia8/Pythia.h"

//==========================================================================

// Functions for histogramming
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/CDFMidPointPlugin.hh"
#include "fastjet/CDFJetCluPlugin.hh"
#include "fastjet/D0RunIIConePlugin.hh"

namespace Pythia8 {

//==========================================================================

// Use userhooks to set the number of requested partons dynamically, as
// needed when running CKKW-L or UMEPS on a single input file that contains
// all parton multiplicities.

class HardProcessBookkeeping : public UserHooks {

public:

  // Constructor and destructor.
  HardProcessBookkeeping() : mergingScheme(0), normFactor(1.) {}
  HardProcessBookkeeping(int mergingSchemeIn)
    : mergingScheme(mergingSchemeIn), normFactor(1.) {}
 ~HardProcessBookkeeping() {}

  double getNormFactor(){return normFactor;}

  // Allow to set the number of partons.
  bool canVetoProcessLevel() { return true;}
  // Set the number of partons.
  bool doVetoProcessLevel( Event& process) {

    int nPartons = 0;
    normFactor = 1.;

    // Do not include resonance decay products in the counting.
    omitResonanceDecays(process);

    // Get the maximal quark flavour counted as "additional" parton.
    int nQuarksMerge = settingsPtr->mode("Merging:nQuarksMerge");

    // Dynamically set the process string.
    if ( settingsPtr->word("Merging:Process") == "guess" ) {
      string processString = "";
      // Set incoming particles.
      int beamAid = beamAPtr->id();
      int beamBid = beamBPtr->id();
      if (abs(beamAid) == 2212) processString += "p";
      if (beamAid == 11)        processString += "e-";
      if (beamAid ==-11)        processString += "e+";
      if (abs(beamBid) == 2212) processString += "p";
      if (beamBid == 11)        processString += "e-";
      if (beamBid ==-11)        processString += "e+";
      processString += ">";
      // Set outgoing particles.
      bool foundOutgoing = false;
      for(int i=0; i < int(workEvent.size()); ++i)
        if ( workEvent[i].isFinal()
          && ( workEvent[i].colType() == 0
            || workEvent[i].idAbs() > 21
            || ( workEvent[i].id() != 21
              && workEvent[i].idAbs() > nQuarksMerge) ) ) {
          foundOutgoing = true;
          ostringstream proc;
          proc << "{" << workEvent[i].name() << "," << workEvent[i].id()
               << "}";
          processString += proc.str();
        }
      // Set the process string.
      if (foundOutgoing) settingsPtr->word("Merging:Process", processString);
    }

    // Loop through event and count.
    for(int i=0; i < int(workEvent.size()); ++i)
      if ( workEvent[i].isFinal()
        && workEvent[i].colType()!= 0
        && ( workEvent[i].id() == 21 || workEvent[i].idAbs() <= nQuarksMerge))
        nPartons++;

    // Store merging scheme.
    bool isumeps  = (mergingScheme == 1);
    bool isunlops = (mergingScheme == 2);

    // Get number of requested partons.
    string nps_nlo = infoPtr->getEventAttribute("npNLO",true);
    int np_nlo     = (nps_nlo != "") ? atoi((char*)nps_nlo.c_str()) : -1;
    string nps_lo  = infoPtr->getEventAttribute("npLO",true);
    int np_lo      = (nps_lo != "") ? atoi((char*)nps_lo.c_str()) : -1;

    if ( (settingsPtr->flag("Merging:doUNLOPSTree")
       || settingsPtr->flag("Merging:doUNLOPSSubt")) && np_lo == 0)
       return true;

    if (settingsPtr->word("Merging:process").compare("pp>aj") == 0)
      nPartons -= 1;
    if (settingsPtr->word("Merging:process").compare("pp>jj") == 0)
      nPartons -= 2;

    // Set number of requested partons.
    if (np_nlo > -1){
      settingsPtr->mode("Merging:nRequested", np_nlo);
      np_lo = -1;
    } else if (np_lo > -1){
      settingsPtr->mode("Merging:nRequested", np_lo);
      np_nlo = -1;
    } else {
      settingsPtr->mode("Merging:nRequested", nPartons);
      np_nlo = -1;
      np_lo = nPartons;
    }

    // Reset the event weight to incorporate corrective factor.
    bool updateWgt = settingsPtr->flag("Merging:includeWeightInXsection");
    double norm    = (abs(infoPtr->lhaStrategy()) == 4) ? 1./1e9 : 1.;

    if (settingsPtr->flag("Merging:doXSectionEstimate")) {
      if (settingsPtr->flag("Merging:doUNLOPSSubt")) {
        settingsPtr->flag("Merging:doUNLOPSTree", true);
        settingsPtr->flag("Merging:doUNLOPSSubt", false);
      }
      if (settingsPtr->flag("Merging:doUNLOPSSubtNLO")) {
        settingsPtr->flag("Merging:doUNLOPSLoop", true);
        settingsPtr->flag("Merging:doUNLOPSSubtNLO", false);
      }
      return false;
    }

    // Choose randomly if this event should be treated as subtraction or
    // as regular event. Put the correct settings accordingly.
    if (isunlops && np_nlo == 0 && np_lo == -1) {
      settingsPtr->flag("Merging:doUNLOPSTree", false);
      settingsPtr->flag("Merging:doUNLOPSSubt", false);
      settingsPtr->flag("Merging:doUNLOPSLoop", true);
      settingsPtr->flag("Merging:doUNLOPSSubtNLO", false);
      settingsPtr->mode("Merging:nRecluster",0);
      normFactor *= 1.;
    } else if (isunlops && np_nlo > 0 && np_lo == -1) {
      normFactor *= 2.;
      if (rndmPtr->flat() < 0.5) {
        normFactor *= -1.;
        settingsPtr->flag("Merging:doUNLOPSTree", false);
        settingsPtr->flag("Merging:doUNLOPSSubt", false);
        settingsPtr->flag("Merging:doUNLOPSLoop", false);
        settingsPtr->flag("Merging:doUNLOPSSubtNLO", true);
        settingsPtr->mode("Merging:nRecluster",1);
      } else {
        settingsPtr->flag("Merging:doUNLOPSTree", false);
        settingsPtr->flag("Merging:doUNLOPSSubt", false);
        settingsPtr->flag("Merging:doUNLOPSLoop", true);
        settingsPtr->flag("Merging:doUNLOPSSubtNLO", false);
        settingsPtr->mode("Merging:nRecluster",0);
      }
    } else if (isunlops && np_nlo == -1 && np_lo > -1) {
      normFactor *= 2.;
      if (rndmPtr->flat() < 0.5) {
        normFactor *= -1.;
        settingsPtr->flag("Merging:doUNLOPSTree", false);
        settingsPtr->flag("Merging:doUNLOPSSubt", true);
        settingsPtr->flag("Merging:doUNLOPSLoop", false);
        settingsPtr->flag("Merging:doUNLOPSSubtNLO", false);
        settingsPtr->mode("Merging:nRecluster",1);

        // Double reclustering for exclusive NLO cross sections.
        bool isnlotilde = settingsPtr->flag("Merging:doUNLOPSTilde");
        int nmaxNLO = settingsPtr->mode("Merging:nJetMaxNLO");
        if ( isnlotilde
          && nmaxNLO > 0
          && np_lo <= nmaxNLO+1
          && np_lo > 1 ){
          normFactor *= 2.;
          if (rndmPtr->flat() < 0.5)
            settingsPtr->mode("Merging:nRecluster",2);
        }
      } else {
        settingsPtr->flag("Merging:doUNLOPSTree", true);
        settingsPtr->flag("Merging:doUNLOPSSubt", false);
        settingsPtr->flag("Merging:doUNLOPSLoop", false);
        settingsPtr->flag("Merging:doUNLOPSSubtNLO", false);
        settingsPtr->mode("Merging:nRecluster",0);
      }
    } else if (isumeps && np_lo == 0) {
      settingsPtr->flag("Merging:doUMEPSTree", true);
      settingsPtr->flag("Merging:doUMEPSSubt", false);
      settingsPtr->mode("Merging:nRecluster",0);
    } else if (isumeps && np_lo > 0) {
      normFactor *= 2.;
      if (rndmPtr->flat() < 0.5) {
        normFactor *= -1.;
        settingsPtr->flag("Merging:doUMEPSTree", false);
        settingsPtr->flag("Merging:doUMEPSSubt", true);
        settingsPtr->mode("Merging:nRecluster",1);
      } else {
        settingsPtr->flag("Merging:doUMEPSTree", true);
        settingsPtr->flag("Merging:doUMEPSSubt", false);
        settingsPtr->mode("Merging:nRecluster",0);
      }
    }

    // Reset the event weight to incorporate corrective factor.
    if ( updateWgt) {
      infoPtr->updateWeight(infoPtr->weight()*norm*normFactor);
      normFactor = 1.;
    }

    // Done
    return false;
  }

private:

  int mergingScheme;
  double normFactor;
};

//==========================================================================

// Class for user interaction with the merging

class PtjTMSdefinitionHooks : public MergingHooks {

private:

public:

  // Default constructor
  PtjTMSdefinitionHooks() {}
  PtjTMSdefinitionHooks(double ptminIn, double ymaxIn)
    : ptmin(ptminIn), ymax(ymaxIn) {}
  // Destructor
  ~PtjTMSdefinitionHooks(){}

  // Functional definition of the merging scale
  virtual double tmsDefinition( const Event& event);

//--------------------------------------------------------------------------

  double ptmin, ymax;

};

//--------------------------------------------------------------------------

// Definition of the merging scale

double PtjTMSdefinitionHooks::tmsDefinition( const Event& event){

  // Only consider first emissions.
  if (!isFirstEmission(event)) { setVeto(-99); return -infoPtr->eCM(); }

  double yPartonMax = 100.;
  double Rparam = 0.4;

  // Fastjet analysis - select algorithm and parameters
  fastjet::Strategy               strategy = fastjet::Best;
  fastjet::RecombinationScheme    recombScheme = fastjet::E_scheme;
  fastjet::JetDefinition         *jetDef
     = new fastjet::JetDefinition(fastjet::kt_algorithm, Rparam,
       recombScheme, strategy);
  // Fastjet input
  std::vector <fastjet::PseudoJet> fjInputs;
  // Reset Fastjet input
  fjInputs.resize(0);

  //SlowJet slowJet(1, 0.4, 0.0, 100., 2, 2, NULL, true);

  // No more shower vetoes if a hard MPI has been produced.
  if (partonSystemsPtr->sizeSys() > 1) { setVeto(-99); return -infoPtr->eCM();}

  //// Construct input for the jet algorithm.
  //Event jetInput;
  //jetInput.init("jet input",particleDataPtr);
  //jetInput.clear();

  // Loop over event record to decide what to pass to FastJet
  for (int i = 0; i < event.size(); ++i) {
    // (Final state && coloured+photons) only!
    if ( !event[i].isFinal()
      || event[i].isLepton()
      || event[i].id() == 23
      || abs(event[i].id()) == 24
      || abs(event[i].y()) > yPartonMax)
      continue;

    // Skip MPI.
    if (partonSystemsPtr->getSystemOf(i,true) > 0) continue;

    // Store as input to Fastjet
    fjInputs.push_back( fastjet::PseudoJet (event[i].px(),
            event[i].py(), event[i].pz(),event[i].e() ) );
    //jetInput.append(Particle(22, 23, 0, 0, 0, 0, 0, 0, event[i].p(), 0, 0));
  }

  // Do nothing for empty input
  if (int(fjInputs.size()) == 0) {
    delete jetDef;
    setVeto(-99);
    return -infoPtr->eCM();
  }

  // Run Fastjet algorithm
  fastjet::ClusterSequence clusterSeq(fjInputs, *jetDef);
  vector<fastjet::PseudoJet> jets
    = fastjet::sorted_by_pt(clusterSeq.inclusive_jets());
  //// Run jet algorithm.
  //slowJet.analyze(jetInput);
  int n(0);
  for (size_t i(0);i<jets.size();++i) {
    Vec4 pj(jets[i].px(),jets[i].py(),jets[i].pz(),jets[i].E());
    if (pj.pT()> ptmin && abs(pj.rap())< ymax) n++;
  }
  //int n1(0);
  //for (size_t i(0);i<slowJet.sizeJet();++i) {
  //  if (slowJet.pT(i)> ptmin && abs(slowJet.y(i))< ymax) n1++;
  //}

  // Get number of requested partons.
  string nps_nlo = infoPtr->getEventAttribute("npNLO",true);
  int np_nlo     = (nps_nlo != "") ? atoi((char*)nps_nlo.c_str()) : -1;
  string nps_lo  = infoPtr->getEventAttribute("npLO",true);
  int np_lo      = (nps_lo != "") ? atoi((char*)nps_lo.c_str()) : -1;
  int n_partons_in_born= (np_nlo >= 0) ? np_nlo : np_lo;

  // For shower veto, remember to keep track of reclusterings in counting.
  n_partons_in_born += -(((doUNLOPSLoop() || doUNLOPSSubtNLO())
                        && (nRecluster()!=0)) ? 1 : 0)
                       + nInProcessNow - n_partons_in_born;

  // Set veto.
  if (n >  n_partons_in_born) setVeto( 1);
  if (n == n_partons_in_born) setVeto( 0);
  if (n <  n_partons_in_born) setVeto(-1);

  return 0;

}

} // end namespace Pythia8

#endif // end Pythia8_PtjTMSdefinitionHooks_H
