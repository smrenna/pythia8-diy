Tune:ee = 7
Beams:frameType             = 5

! Specify merging parameters for CKKW-L, UMEPS, UNLOPS.
Merging:TMS                 = 20.                  ! merging scale value
Merging:Process             = pp>LEPTONS,NEUTRINOS ! process definition
Merging:nJetMax             = 8        ! maximal number of additional LO jets --- globally!!!
Merging:nJetMaxNLO          = 0        ! maximal number of additional NLO jets

! Wimpy shower
TimeShower:pTmaxMatch       = 1
SpaceShower:pTmaxMatch      = 1

! Factorisation/renormalisation scales in the 2->2 process
Merging:muFac               = 91.188
Merging:muRen               = 91.188
Merging:muFacInME           = 91.188
Merging:muRenInME           = 91.188

! Use same PDFs / alpha_s value as in ME calculation (not necessary!)
SpaceShower:alphaSvalue   = 0.118
TimeShower:alphaSvalue    = 0.118
PartonLevel:MPI = Off
#PartonLevel:ISR = Off

! Do not include rapidity ordering (not necessary!)
SpaceShower:rapidityOrder = off

! Be more forgiving with momentum mismatches.
Check:epTolErr               = 2e-2

#Merging:enforceCutOnLHE = off



Merging:applyVeto = off # prevent py8 from sucking in events
Merging:includeWeightInXsection = off # lass lieber aus



! Subruns for UNLOPS NLO merging
#Merging:doUNLOPSTilde     = on
#LHEFInputs:nSubruns       = 1
#LHEFInputs:nSubruns       = 2
Main:subrun               = 0
#Merging:doUNLOPSTree      = off
#Merging:doUNLOPSSubt      = off
#Merging:doUNLOPSLoop      = on
#Merging:doUNLOPSSubtNLO   = off
#Beams:LHEF                = zProduction_UnlopsLoop_01.lhe.gz
#Main:subrun               = 1
#Merging:doUNLOPSTree      = on
#Merging:doUNLOPSSubt      = off
#Merging:doUNLOPSLoop      = off
#Merging:doUNLOPSSubtNLO   = off
Merging:doUserMerging   = on
#Beams:LHEF                = zProduction_UnlopsTree_12.lhe.gz

