// main89.cc is a part of the PYTHIA event generator.
// Copyright (C) 2018 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This program is written by Stefan Prestel.
// It illustrates how to do run PYTHIA with LHEF input, allowing a
// sample-by-sample generation of
// a) Non-matched/non-merged events
// b) MLM jet-matched events (kT-MLM, shower-kT, FxFx)
// c) CKKW-L and UMEPS-merged events
// d) UNLOPS NLO merged events
// see the respective sections in the online manual for details.

#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/HepMC2.h"
#include <unistd.h>
#include "PtjTMSdefinitionHooks.h"

using namespace Pythia8;

//==========================================================================

// Example main programm to illustrate merging.

int main( int argc, char* argv[] ){

  // Check that correct number of command-line arguments
  if (argc != 3) {
    cerr << " Unexpected number of command-line arguments ("<<argc<<"). \n"
         << " You are expected to provide the arguments" << endl
         << " 1. Input file for settings" << endl
         << " 2. Output file for HepMC events" << endl
         << " Program stopped. " << endl;
    return 1;
  }

  Pythia pythia;

  // Input parameters:
  pythia.readFile(argv[1],0);

  pythia.readString("Merging:unlopsTMSdefinition = 1");
  int unlopsType = pythia.settings.mode("Merging:unlopsTMSdefinition");

  // Histograms combined over all jet multiplicities.
  Hist ptlund("pTlund", 100, 0., 200.);

  MergingHooks* ptjTMSdefinitionPtr = (unlopsType<0)
    ? NULL
    : new PtjTMSdefinitionHooks(pythia.parm("Merging:TMS"),6.0, &ptlund);
  if (unlopsType >0) pythia.setMergingHooksPtr( ptjTMSdefinitionPtr );

  // Interface for conversion from Pythia8::Event to HepMC one.
  HepMC::Pythia8ToHepMC ToHepMC;
  // Specify file where HepMC events will be stored.
  HepMC::IO_GenEvent ascii_io(argv[2], std::ios::out);
  // Switch off warnings for parton-level events.
  ToHepMC.set_print_inconsistency(false);
  ToHepMC.set_free_parton_exception(false);
  // Do not store cross section information, as this will be done manually.
  ToHepMC.set_store_pdf(false);
  ToHepMC.set_store_proc(false);
  ToHepMC.set_store_xsec(false);

  // Get number of subruns.
  int nMerge = pythia.mode("Main:numberOfSubruns");

  // Number of events. Negative numbers mean all events in the LHEF will be
  // used.
  long nEvent = pythia.settings.mode("Main:numberOfEvents");
  if (nEvent < 1) nEvent = 1000000000000000;

  // Allow to set the number of addtional partons dynamically.
  HardProcessBookkeeping* hardProcessBookkeeping = NULL;
  // Store merging scheme.
  int scheme = ( pythia.settings.flag("Merging:doUMEPSTree")
              || pythia.settings.flag("Merging:doUMEPSSubt")) ?
              1 :
               ( ( pythia.settings.flag("Merging:doUNLOPSTree")
              || pythia.settings.flag("Merging:doUNLOPSSubt")
              || pythia.settings.flag("Merging:doUNLOPSLoop")
              || pythia.settings.flag("Merging:doUNLOPSSubtNLO")) ?
              2 :
              0 );
  HardProcessBookkeeping* hardProcessBookkeepingPtr
    = new HardProcessBookkeeping(scheme);
    pythia.setUserHooksPtr(hardProcessBookkeepingPtr);

  vector<double> xss;

  // Cross section and error.
  double sigmaTotal  = 0.;
  double errorTotal  = 0.;

  // Allow abort of run if many errors.
  int  nAbort  = pythia.mode("Main:timesAllowErrors");
  int  iAbort  = 0;
  bool doAbort = false;

  cout << endl << endl << endl;
  cout << "Start generating events" << endl;

  // Loop over subruns with varying number of jets.
  for (int iMerge = 0; iMerge < nMerge; ++iMerge) {

    double sigmaSample = 0., errorSample = 0.;

    // Read in name of LHE file for current subrun and initialize.
    pythia.readFile(argv[1], iMerge);

    // If the process string is "guess", temporarily set it to something safe
    // for initialization.
    bool doGuess = pythia.settings.word("Merging:process") == "guess";
    if (doGuess) pythia.settings.word("Merging:process","pp>e+e-");
    // Initialise.
    pythia.init();
    // Reset the process string to "guess" if necessary.
    if (doGuess) pythia.settings.word("Merging:process","guess");

    // Get the inclusive x-section by summing over all process x-sections.
    double xs = 0.;
    for (int i=0; i < pythia.info.nProcessesLHEF(); ++i)
      xs += pythia.info.sigmaLHEF(i);

    // Start generation loop
    while( pythia.info.nSelected() < nEvent ){

      // Generate next event
      if( !pythia.next() ) {
        if ( pythia.info.atEndOfFile() ) break;
        else if (++iAbort > nAbort) {doAbort = true; break;}
        else continue;
      }

      // Get event weight(s).
      double evtweight         = pythia.info.weight();
      // Additional PDF/alphaS weight for internal merging.
      evtweight               *= pythia.info.mergingWeightNLO()
      // Additional weight due to random choice of reclustered/non-reclustered
      // treatment. Also contains additional sign for subtractive samples.
                                *hardProcessBookkeepingPtr->getNormFactor();

      // Do not print zero-weight events.
      if ( evtweight == 0. ) continue;
      // Construct new empty HepMC event.
      HepMC::GenEvent* hepmcevt = new HepMC::GenEvent();

      // Work with weighted (LHA strategy=-4) events.
      double normhepmc = 1.;
      if (abs(pythia.info.lhaStrategy()) == 4)
        normhepmc = 1. / double(1e9*nEvent);
      // Work with unweighted events.
      else
        normhepmc = xs / double(1e9*nEvent);

      // Set event weight
      hepmcevt->weights().push_back(evtweight*normhepmc);
      // Fill HepMC event
      ToHepMC.fill_next_event( pythia, hepmcevt );
      // Add the weight of the current event to the cross section.
      sigmaTotal  += evtweight*normhepmc;
      sigmaSample += evtweight*normhepmc;
      errorTotal  += pow2(evtweight*normhepmc);
      errorSample += pow2(evtweight*normhepmc);
      // Report cross section to hepmc
      HepMC::GenCrossSection xsec;
      xsec.set_cross_section( sigmaTotal*1e9, pythia.info.sigmaErr()*1e9 );
      hepmcevt->set_cross_section( xsec );
      // Write the HepMC event to file. Done with it.
      ascii_io << hepmcevt;
      delete hepmcevt;

    } // end loop over events to generate.
    if (doAbort) break;

    // print cross section, errors
    pythia.stat();

    cout << endl << " Contribution of sample " << iMerge
         << " to the inclusive cross section : "
         << scientific << setprecision(8)
         << sigmaSample << "  +-  " << sqrt(errorSample)  << endl;

  }

  cout << endl << endl << endl;
  if (doAbort)
    cout << " Run was not completed owing to too many aborted events" << endl;
  else
    cout << "Inclusive cross section: " << scientific << setprecision(8)
         << sigmaTotal << "  +-  " << sqrt(errorTotal) << " mb " << endl;
  cout << endl << endl << endl;

  ofstream write;
  write.open("ptlund.dat");
  ptlund.table(write);
  write.close();

  // Clean-up
  delete hardProcessBookkeepingPtr;
  if (unlopsType>0) delete ptjTMSdefinitionPtr;

  // Done
  return 0;

}
