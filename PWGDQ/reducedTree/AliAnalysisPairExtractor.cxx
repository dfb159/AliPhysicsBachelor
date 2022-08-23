//
// Implementation of the DataExtractor for further TMVA analysis for Pairs
//
// Creation date: 2022/05/05
// Author: Jonathan Sigrist, j.sigrist@web.de

#include <TClonesArray.h>
#include "AliAnalysisPairExtractor.h"
#include <iostream>
#include <dirent.h>
#include <string.h>

ClassImp(AliAnalysisPairExtractor)

using std::cout;
using std::endl;
using std::flush;

//________________________________________________________________________________

AliAnalysisPairExtractor::AliAnalysisPairExtractor(AliReducedPairInfo::CandidateType type) {
  switch(type) {
  case AliReducedPairInfo::CandidateType::kJpsiToEE: AliAnalysisPairExtractor(443, -11, 11); break; // positron is Leg1, electron is Leg2
  default: Error("AliAnalysisPairExtractor", "Unknown Candidate Type. Please set PDGs manually!"); AliAnalysisPairExtractor(0, 0, 0);
  }
}

AliAnalysisPairExtractor::AliAnalysisPairExtractor(int mother, int leg1, int leg2) : 
  pdgMother(mother), pdgLeg1(leg1), pdgLeg2(leg2) {}

void AliAnalysisPairExtractor::SetUp(TString outpath) {
	outfile = TFile::Open(outpath, "RECREATE");
	
	if (!outfile || !outfile->IsOpen()) {cout << "File could not be opened!" << endl; return;}
	
  dataTree = new TTree("data", "candidates from experimental data"); createBranches(dataTree);
  signalTree = new TTree("mcSignal", "generated candidates with confirmed origin"); createBranches(signalTree);
  backgroundTree = new TTree("mcBackground", "generated candidates with rejected origin"); createBranches(backgroundTree);
}

void AliAnalysisPairExtractor::Write() {
  outfile->cd();
  dataTree->Write();
  signalTree->Write();
  backgroundTree->Write();
  outfile->Close();
}

Bool_t AliAnalysisPairExtractor::checkTreeIntegrity(TTree* tree, Bool_t debug) {
  Bool_t state = kTRUE;
    
  if (!state && debug) Error("AliAnalysisPairExtractor::checkTreeIntegrity", "The given TTree has an invalid data structure");
  // check if tree is of the correct class
  // check if the given tree contains all of the needed information
  return state;
}

void AliAnalysisPairExtractor::extractDataDirectory(TString path, TString treeName) {
  DIR *dir;
  struct dirent *ent;
  if ((dir = opendir(path.Data())) != NULL) {
    while ((ent = readdir(dir)) != NULL) {
      string name(ent->d_name);
      if (name.find(".") != std::string::npos) continue;
      extractDataFile(path + name + "/JpsiCandidates_data.root");
      cout << "Data Run: " << name << endl;
    }
    closedir (dir);
  }
}

void AliAnalysisPairExtractor::extractDataFile(TString path, TString treeName) {
  TFile* fin = TFile::Open(path.Data(), "READ");
  if (fin->GetListOfKeys()->GetEntries() > 0) {
    TTree* intree = fin->Get<TTree>(treeName.Data());
    if (checkTreeIntegrity(intree)) extractData(intree);

  }
  fin->Close();
}

void AliAnalysisPairExtractor::extractData(TTree* intree) { // extraction method for data variable names
  TClonesArray* fTracks;
  TClonesArray* fPairs;
  AliReducedTrackInfo* leg1;
  AliReducedTrackInfo* leg2;
  AliReducedTrackInfo* track;
  AliReducedPairInfo* pair;
  AliReducedEventInfo* event = new AliReducedEventInfo();
  intree->SetBranchAddress("Event",&event);
  
  
  int n = intree->GetEntries();
  for(int i = 0; i < n; i++) { // for every event
  	//if (i > 1000) break;
    intree->GetEntry(i);
 
    fTracks = event->GetTracks();
    fPairs = event->GetPairs();
    TIter nextTrack(fTracks);
    TIter nextPair(fPairs);
    for (int j = 0; j < fPairs->GetEntries(); j++) { // for every candidate Pair
      pair = (AliReducedPairInfo*) nextPair();
      
      //cout << "Event " << i << " : Pair " << j << ", leg1 = " << pair->LegId(0) << ", leg2 = " << pair->LegId(1) << endl;
      
      //cout << "TrackIDs: ";
      // search for both legs: leg1.charge > 0, leg2.charge < 0
      nextTrack.Reset(); leg1 = 0x0; leg2 = 0x0;
      for (int k = 0; k < fTracks->GetEntries(); k++) {
        track = (AliReducedTrackInfo*) nextTrack(); // TODO does this overwrite data or create new object (keep leg1, leg2 data)
        if (track->IsMCTruth()) continue;
        if (track->TrackId() == pair->LegId(0)) {leg1 = (AliReducedTrackInfo*) fTracks->At(k);} // does this create a nice copy
        if (track->TrackId() == pair->LegId(1)) {leg2 = (AliReducedTrackInfo*) fTracks->At(k);}
        //cout << track->TrackId() << ", ";
      }
      //cout << endl;
      
      if (!leg1 || !leg2) {Error("AliAnalysisPairExtractor::extractData", "Could not find a corresponding track for both legs"); continue;}
      
      fillVars(event, pair, leg1, leg2);
      dataTree->Fill();
      //cout << "fill Data" << endl;
    }
  }
  intree->ResetBranchAddresses();
}

void AliAnalysisPairExtractor::extractMCDirectory(TString path, TString treeName) {
  DIR *dir;
  struct dirent *ent;
  if ((dir = opendir(path.Data())) != NULL) {
    while ((ent = readdir(dir)) != NULL) {
      string name(ent->d_name);
      if (name.find(".") != std::string::npos) continue;
      extractMCFile(path + name + "/JpsiCandidates_mc.root");
      cout << "MC Run: " << name << endl;
    }
    closedir (dir);
  }
}

void AliAnalysisPairExtractor::extractMCFile(TString path, TString treeName) {
  TFile* fin = TFile::Open(path.Data(), "READ");
  if (fin->GetListOfKeys()->GetEntries() > 0) {
    TTree* intree = fin->Get<TTree>(treeName.Data());
    extractMC(intree);
  }
  fin->Close();
}
void AliAnalysisPairExtractor::extractMC(TTree* intree) { // extraction method for MC variable names

  TClonesArray* fTracks;
  TClonesArray* fPairs;
  AliReducedTrackInfo* leg1;
  AliReducedTrackInfo* leg2;
  AliReducedTrackInfo* track;
  AliReducedPairInfo* pair;
  AliReducedEventInfo* event = new AliReducedEventInfo();
  intree->SetBranchAddress("Event",&event);
  
  
  int n = intree->GetEntries();
  for(int i = 0; i < n; i++) { // for every event
    intree->GetEntry(i);
 
    fTracks = event->GetTracks();
    fPairs = event->GetPairs();
    TIter nextTrack(fTracks);
    TIter nextPair(fPairs);
    for (int j = 0; j < fPairs->GetEntries(); j++) { // for every candidate Pair
      pair = (AliReducedPairInfo*) nextPair();
      
      //cout << "Event " << i << " : Pair " << j << ", leg1 = " << pair->LegId(0) << ", leg2 = " << pair->LegId(1) << endl;

      //cout << "TrackIDs: ";
      // search for both legs: leg1.charge > 0, leg2.charge < 0
      nextTrack.Reset(); leg1 = 0x0; leg2 = 0x0;
      for (int k = 0; k < fTracks->GetEntries(); k++) {
        track = (AliReducedTrackInfo*) nextTrack(); // TODO does this overwrite data or create new object (keep leg1, leg2 data)
        if (track->IsMCTruth()) continue; // not a reconstructed track
        if (track->TrackId() == pair->LegId(0)) {leg1 = (AliReducedTrackInfo*) fTracks->At(k);} // does this create a nice copy
        if (track->TrackId() == pair->LegId(1)) {leg2 = (AliReducedTrackInfo*) fTracks->At(k);}
        //cout << track->TrackId() << ", ";
      }
      //cout << endl;
      
      if (!leg1 || !leg2) {Error("AliAnalysisPairExtractor::extractMC", "Could not find a corresponding track for both legs"); continue;}
      //if (pdgMother != leg1->MCPdg(1) || pdgMother != leg2->MCPdg(1)) {Error("AliAnalysisPairExtractor::extractMC", "One of the legs has the wrong mother track ID"); continue;}

      /*
      cout << "leg1pdg: " << leg1->MCPdg(0) << " & " << pdgLeg1 << " = " << (leg1->MCPdg(0) == pdgLeg1) << endl;
      cout << "leg2pdg: " << leg2->MCPdg(0) << " & " << pdgLeg2 << " = " << (leg2->MCPdg(0) == pdgLeg2) << endl;
      cout << "mother1pdg: " << leg1->MCPdg(1) << " & " << pdgMother << " = " << (leg1->MCPdg(1) == pdgMother) << endl;
      cout << "mother2pdg: " << leg2->MCPdg(1) << " & " << pdgMother << " = " << (leg2->MCPdg(1) == pdgMother) << endl;
      cout << "motherlabel: " << leg1->MCLabel(1) << " & " << leg2->MCLabel(1) << " = " <<  (leg1->MCLabel(1) == leg2->MCLabel(1)) << endl;
      
      //*/

      fillVars(event, pair, leg1, leg2);
      if (pdgLeg1 == leg1->MCPdg(0) && pdgLeg2 == leg2->MCPdg(0) // both tracks are electrons
        && pdgMother == leg1->MCPdg(1) && pdgMother == leg2->MCPdg(1) // both mothers are JPsi
        && leg1->MCLabel(1) == leg2->MCLabel(1)) { // both mothers are the same particle
        signalTree->Fill();
        //cout << "fill Signal MC" << endl;
      } else { // is a wrongly identified JPsi2ee decay
        backgroundTree->Fill();
        //cout << "fill Background MC" << endl;
      }
    }
  }
  intree->ResetBranchAddresses();
}

void AliAnalysisPairExtractor::createBranches(TTree* tree) {
  
  // Event Information
  tree->Branch("event_runNo", &runNo, "Run No/I");
  tree->Branch("event_vtx_X", &vtxX, "Vertex X/F");
  tree->Branch("event_vtx_Y", &vtxY, "Vertex Y/F");
  tree->Branch("event_vtx_Z", &vtxZ, "Vertex Z/F");
  tree->Branch("event_vtx_N", &vtxN, "Vertex N contributors/I");
  tree->Branch("event_diamond_X", &diaX, "Diamond X/F");
  tree->Branch("event_diamond_Y", &diaY, "Diamond Y/F");
  tree->Branch("event_diamond_Z", &diaZ, "Diamond Z/F");
  tree->Branch("event_totalTracks", &totalTracks, "Total reconstructed Tracks/I");
  tree->Branch("event_centrality_V0", &centV0, "Centrality V0/F");
  tree->Branch("event_centrality_SPD", &centSPD, "Centrality SPD/F");
  tree->Branch("event_centrality_TPC", &centTPC, "Centrality TPC/F");
  tree->Branch("event_centrality_ZEMvsZDC", &centZEMvsZDC, "Centrality ZEMvsZDC/F");
  tree->Branch("event_centrality_V0A", &centV0A, "Centrality V0A/F");
  tree->Branch("event_centrality_V0C", &centV0C, "Centrality V0C/F");
  tree->Branch("event_centrality_ZNA", &centZNA, "Centrality ZNA/F");

  tree->Branch("event_centrality_V0M_new", &centV0Mnew, "Centrality V0M new/F");
  tree->Branch("event_centrality_V0M_new+05", &centV0MnewPlus05, "Centrality V0M new + 05/F");
  tree->Branch("event_centrality_V0M_new-05", &centV0MnewMinus05, "Centrality V0M new - 05/F");
  tree->Branch("event_centrality_V0M_new+10", &centV0MnewPlus10, "Centrality V0M new + 10/F");
  tree->Branch("event_centrality_V0M_new-05", &centV0MnewMinus10, "Centrality V0M new - 10/F");
  tree->Branch("event_centrality_V0M+05", &centV0MPlus05, "Centrality V0M + 05/F");
  tree->Branch("event_centrality_V0M-05", &centV0MMinus05, "Centrality V0M - 05/F");
  tree->Branch("event_centrality_V0M+10", &centV0MPlus10, "Centrality V0M + 10/F");
  tree->Branch("event_centrality_V0M-10", &centV0MMinus10, "Centrality V0M - 10/F");
  
  // Pair Information
  tree->Branch("pair_P", &pairP, "Pair P/F");
  tree->Branch("pair_Pt", &pairPt, "Pair Pt/F");
  tree->Branch("pair_Phi", &pairPhi, "Pair Phi/F");
  tree->Branch("pair_Theta", &pairTheta, "Pair Theta/F");
  tree->Branch("pair_Eta", &pairEta, "Pair Eta/F");
  tree->Branch("pair_Mass", &pairMass, "Pair Mass/F");
  tree->Branch("pair_Energy", &pairEnergy, "Pair Energy/F");
  tree->Branch("pair_Rapidity", &pairRapidity, "Pair Rapidity/F");
  tree->Branch("pair_DecayRadius", &pairDecayRadius, "Pair Decay Radius/F");
  tree->Branch("pair_PsProper", &pairPsProper, "Pair Proper Decay Radius/F");
  tree->Branch("pair_PointingAngle", &pairPointingAngle, "Pair Pointing Angle/F");
  tree->Branch("pair_Chi2", &pairChi2, "Pair Chi2/F");
  
  // Track Information Leg1
  tree->Branch("leg1_P", &leg1p, "Leg1 P/F");
  tree->Branch("leg1_Pt", &leg1pt,"Leg1 Pt/F");
  tree->Branch("leg1_Phi", &leg1phi, "Leg1 Phi/F");
  tree->Branch("leg1_Theta", &leg1theta, "Leg1 Theta/F");
  tree->Branch("leg1_Eta", &leg1eta, "Leg1 Eta/F");
  tree->Branch("leg1_Charge", &leg1charge, "Leg1 Charge/I");
  
  tree->Branch("leg1_massTracking", &leg1massTracking, "Leg1 Mass Tracking/F");
  tree->Branch("leg1_trackLength", &leg1trackLength, "Leg1 Track Length/F");
  tree->Branch("leg1_dcaXY", &leg1dcaXY, "Leg1 DCA XY/F");
  tree->Branch("leg1_dcaZ", &leg1dcaZ, "Leg1 DCA Z/F");
  tree->Branch("leg1_helixX", &leg1helixX, "Leg1 HelixX/F");
  tree->Branch("leg1_helixY", &leg1helixY, "Leg1 HelixY/F");
  tree->Branch("leg1_helixR", &leg1helixR, "Leg1 HelixR/F");
  
  tree->Branch("leg1_its_nSigElec", &leg1itsnSigElec, "Leg1 ITS nSig Electron/F");
  tree->Branch("leg1_its_nSigPion", &leg1itsnSigPion, "Leg1 ITS nSig Pion/F");
  tree->Branch("leg1_its_nSigKaon", &leg1itsnSigKaon, "Leg1 ITS nSig Kaon/F");
  tree->Branch("leg1_its_nSigProt", &leg1itsnSigProt, "Leg1 ITS nSig Proton/F");
  tree->Branch("leg1_its_signal", &leg1itsSignal, "Leg1 ITS Signal/F");
  tree->Branch("leg1_its_Chi2", &leg1itsChi2, "Leg1 ITS Chi2/F");
  tree->Branch("leg1_its_Ncls", &leg1itsNcls, "Leg1 ITS N clusters/i");
  tree->Branch("leg1_its_NclsShared", &leg1itsNclsShared, "Leg1 ITS N shared clusters/i");
  
  tree->Branch("leg1_tpc_nSigElec", &leg1tpcnSigElec, "Leg1 TPC nSig Electron/F");
  tree->Branch("leg1_tpc_nSigPion", &leg1tpcnSigPion, "Leg1 TPC nSig Pion/F");
  tree->Branch("leg1_tpc_nSigKaon", &leg1tpcnSigKaon, "Leg1 TPC nSig Kaon/F");
  tree->Branch("leg1_tpc_nSigProt", &leg1tpcnSigProt, "Leg1 TPC nSig Proton/F");
  tree->Branch("leg1_tpc_signal", &leg1tpcSignal, "Leg1 TPC dEdx/F");
  tree->Branch("leg1_tpc_signal_tuned", &leg1tpcSignalTuned, "Leg1 TPC dEdx Signal tuned/F");
  tree->Branch("leg1_tpc_signalN", &leg1tpcSignalN, "Leg1 TPC no clusters dEdx/i");
  tree->Branch("leg1_tpc_Chi2", &leg1tpcChi2, "Leg1 TPC Chi2/F");
  tree->Branch("leg1_tpc_Ncls", &leg1tpcNcls, "Leg1 TPC N clusters/i");
  tree->Branch("leg1_tpc_NclsShared", &leg1tpcNclsShared, "Leg1 TPC N shared clusters/i");
  tree->Branch("leg1_tpc_Rows", &leg1tpcRows, "Leg1 TPC crossed rows/i");
  tree->Branch("leg1_tpc_activeLength", &leg1tpcActiveLength, "Leg1 TPC active length/F");
  tree->Branch("leg1_tpc_geomLength", &leg1tpcGeomLength, "Leg1 TPC geometric length/F");
  tree->Branch("leg1_tpc_dEdx_Qmax_IROC", &leg1tpcdEdxQmaxIROC, "Leg1 TPC dEdx Qmax IROC/F");
  tree->Branch("leg1_tpc_dEdx_Qmax_medOROC", &leg1tpcdEdxQmaxMedOROC, "Leg1 TPC dEdx Qmax medium OROC/F");
  tree->Branch("leg1_tpc_dEdx_Qmax_longOROC", &leg1tpcdEdxQmaxLongOROC, "Leg1 TPC dEdx Qmax long OROC/F");
  tree->Branch("leg1_tpc_dEdx_Qmax_allOROC", &leg1tpcdEdxQmaxAllOROC, "Leg1 TPC dEdx Qmax all OROC/F");
  tree->Branch("leg1_tpc_dEdx_Qtot_IROC", &leg1tpcdEdxQtotIROC, "Leg1 TPC dEdx Qtot IROC/F");
  tree->Branch("leg1_tpc_dEdx_Qtot_medOROC", &leg1tpcdEdxQtotMedOROC, "Leg1 TPC dEdx Qtot medium OROC/F");
  tree->Branch("leg1_tpc_dEdx_Qtot_longOROC", &leg1tpcdEdxQtotLongOROC, "Leg1 TPC dEdx Qtot long OROC/F");
  tree->Branch("leg1_tpc_dEdx_Qtot_allOROC", &leg1tpcdEdxQtotAllOROC, "Leg1 TPC dEdx Qtot all OROC/F");
  
  tree->Branch("leg1_tof_nSigElec", &leg1tofnSigElec, "Leg1 TOF nSig Electron/F");
  tree->Branch("leg1_tof_nSigPion", &leg1tofnSigPion, "Leg1 TOF nSig Pion/F");
  tree->Branch("leg1_tof_nSigKaon", &leg1tofnSigKaon, "Leg1 TOF nSig Kaon/F");
  tree->Branch("leg1_tof_nSigProt", &leg1tofnSigProt, "Leg1 TOF nSig Proton/F");
  tree->Branch("leg1_tof_Beta", &leg1tofBeta, "Leg1 TPC Beta/F");
  tree->Branch("leg1_tof_Time", &leg1tofTime, "Leg1 TPC Time/F");
  tree->Branch("leg1_tof_MismatchProbab", &leg1tofMisProbab, "Leg1 TPC Mismatch Probability/F");
  tree->Branch("leg1_tof_Chi2", &leg1tofChi2, "Leg1 TPC Chi2/F");
  
  tree->Branch("leg1_trd_Ntracklets", &leg1trdNtracklets, "Leg1 TRD N tracklets/i");
  tree->Branch("leg1_trd_NtrackletsPID", &leg1trdNtrackletsPID, "Leg1 TRD N tracklets PID/i");
  tree->Branch("leg1_trd_Sig1D_Elec", &leg1trdSig1DElec, "Leg1 TRD Sig 1D Electron/i");
  tree->Branch("leg1_trd_Sig1D_Pion", &leg1trdSig1DPion, "Leg1 TRD Sig 1D Pion/i");
  tree->Branch("leg1_trd_Sig2D_Elec", &leg1trdSig2DElec, "Leg1 TRD Sig 2D Electron/i");
  tree->Branch("leg1_trd_Sig2D_Pion", &leg1trdSig2DPion, "Leg1 TRD Sig 2D Pion/i");
  
  // Track Information Leg2
  tree->Branch("leg2_P", &leg2p, "Leg2 P/F");
  tree->Branch("leg2_Pt", &leg2pt,"Leg2 Pt/F");
  tree->Branch("leg2_Phi", &leg2phi, "Leg2 Phi/F");
  tree->Branch("leg2_Theta", &leg2theta, "Leg2 Theta/F");
  tree->Branch("leg2_Eta", &leg2eta, "Leg2 Eta/F");
  tree->Branch("leg2_Charge", &leg2charge, "Leg2 Charge/I");
  
  tree->Branch("leg2_massTracking", &leg2massTracking, "Leg2 Mass Tracking/F");
  tree->Branch("leg2_trackLength", &leg2trackLength, "Leg2 Track Length/F");
  tree->Branch("leg2_dcaXY", &leg2dcaXY, "Leg2 DCA XY/F");
  tree->Branch("leg2_dcaZ", &leg2dcaZ, "Leg2 DCA Z/F");
  tree->Branch("leg2_helixX", &leg2helixX, "Leg2 HelixX/F");
  tree->Branch("leg2_helixY", &leg2helixY, "Leg2 HelixY/F");
  tree->Branch("leg2_helixR", &leg2helixR, "Leg2 HelixR/F");
  
  tree->Branch("leg2_its_nSigElec", &leg2itsnSigElec, "Leg2 ITS nSig Electron/F");
  tree->Branch("leg2_its_nSigPion", &leg2itsnSigPion, "Leg2 ITS nSig Pion/F");
  tree->Branch("leg2_its_nSigKaon", &leg2itsnSigKaon, "Leg2 ITS nSig Kaon/F");
  tree->Branch("leg2_its_nSigProt", &leg2itsnSigProt, "Leg2 ITS nSig Proton/F");
  tree->Branch("leg2_its_signal", &leg2itsSignal, "Leg2 ITS Signal/F");
  tree->Branch("leg2_its_Chi2", &leg2itsChi2, "Leg2 ITS Chi2/F");
  tree->Branch("leg2_its_Ncls", &leg2itsNcls, "Leg2 ITS N clusters/i");
  tree->Branch("leg2_its_NclsShared", &leg2itsNclsShared, "Leg2 ITS N shared clusters/i");
  
  tree->Branch("leg2_tpc_nSigElec", &leg2tpcnSigElec, "Leg2 TPC nSig Electron/F");
  tree->Branch("leg2_tpc_nSigPion", &leg2tpcnSigPion, "Leg2 TPC nSig Pion/F");
  tree->Branch("leg2_tpc_nSigKaon", &leg2tpcnSigKaon, "Leg2 TPC nSig Kaon/F");
  tree->Branch("leg2_tpc_nSigProt", &leg2tpcnSigProt, "Leg2 TPC nSig Proton/F");
  tree->Branch("leg2_tpc_signal", &leg2tpcSignal, "Leg2 TPC dEdx/F");
  tree->Branch("leg2_tpc_signal_tuned", &leg2tpcSignalTuned, "Leg2 TPC dEdx Signal tuned/F");
  tree->Branch("leg2_tpc_signalN", &leg2tpcSignalN, "Leg2 TPC no clusters dEdx/i");
  tree->Branch("leg2_tpc_Chi2", &leg2tpcChi2, "Leg2 TPC Chi2/F");
  tree->Branch("leg2_tpc_Ncls", &leg2tpcNcls, "Leg2 TPC N clusters/i");
  tree->Branch("leg2_tpc_NclsShared", &leg2tpcNclsShared, "Leg2 TPC N shared clusters/i");
  tree->Branch("leg2_tpc_Rows", &leg2tpcRows, "Leg2 TPC crossed rows/i");
  tree->Branch("leg2_tpc_activeLength", &leg2tpcActiveLength, "Leg2 TPC active length/F");
  tree->Branch("leg2_tpc_geomLength", &leg2tpcGeomLength, "Leg2 TPC geometric length/F");
  tree->Branch("leg2_tpc_dEdx_Qmax_IROC", &leg2tpcdEdxQmaxIROC, "Leg2 TPC dEdx Qmax IROC/F");
  tree->Branch("leg2_tpc_dEdx_Qmax_medOROC", &leg2tpcdEdxQmaxMedOROC, "Leg2 TPC dEdx Qmax medium OROC/F");
  tree->Branch("leg2_tpc_dEdx_Qmax_longOROC", &leg2tpcdEdxQmaxLongOROC, "Leg2 TPC dEdx Qmax long OROC/F");
  tree->Branch("leg2_tpc_dEdx_Qmax_allOROC", &leg2tpcdEdxQmaxAllOROC, "Leg2 TPC dEdx Qmax all OROC/F");
  tree->Branch("leg2_tpc_dEdx_Qtot_IROC", &leg2tpcdEdxQtotIROC, "Leg2 TPC dEdx Qtot IROC/F");
  tree->Branch("leg2_tpc_dEdx_Qtot_medOROC", &leg2tpcdEdxQtotMedOROC, "Leg2 TPC dEdx Qtot medium OROC/F");
  tree->Branch("leg2_tpc_dEdx_Qtot_longOROC", &leg2tpcdEdxQtotLongOROC, "Leg2 TPC dEdx Qtot long OROC/F");
  tree->Branch("leg2_tpc_dEdx_Qtot_allOROC", &leg2tpcdEdxQtotAllOROC, "Leg2 TPC dEdx Qtot all OROC/F");
  
  tree->Branch("leg2_tof_nSigElec", &leg2tofnSigElec, "Leg2 TOF nSig Electron/F");
  tree->Branch("leg2_tof_nSigPion", &leg2tofnSigPion, "Leg2 TOF nSig Pion/F");
  tree->Branch("leg2_tof_nSigKaon", &leg2tofnSigKaon, "Leg2 TOF nSig Kaon/F");
  tree->Branch("leg2_tof_nSigProt", &leg2tofnSigProt, "Leg2 TOF nSig Proton/F");
  tree->Branch("leg2_tof_Beta", &leg2tofBeta, "Leg2 TPC Beta/F");
  tree->Branch("leg2_tof_Time", &leg2tofTime, "Leg2 TPC Time/F");
  tree->Branch("leg2_tof_MismatchProbab", &leg2tofMisProbab, "Leg2 TPC Mismatch Probability/F");
  tree->Branch("leg2_tof_Chi2", &leg2tofChi2, "Leg2 TPC Chi2/F");
  
  tree->Branch("leg2_trd_Ntracklets", &leg2trdNtracklets, "Leg2 TRD N tracklets/i");
  tree->Branch("leg2_trd_NtrackletsPID", &leg2trdNtrackletsPID, "Leg2 TRD N tracklets PID/i");
  tree->Branch("leg2_trd_Sig1D_Elec", &leg2trdSig1DElec, "Leg2 TRD Sig 1D Electron/i");
  tree->Branch("leg2_trd_Sig1D_Pion", &leg2trdSig1DPion, "Leg2 TRD Sig 1D Pion/i");
  tree->Branch("leg2_trd_Sig2D_Elec", &leg2trdSig2DElec, "Leg2 TRD Sig 2D Electron/i");
  tree->Branch("leg2_trd_Sig2D_Pion", &leg2trdSig2DPion, "Leg2 TRD Sig 2D Pion/i");
}

void AliAnalysisPairExtractor::fillVars(AliReducedEventInfo* event, AliReducedPairInfo* pair, AliReducedTrackInfo* leg1, AliReducedTrackInfo* leg2) {
  
  // Event Information
  runNo = event->RunNo();
  vtxX = event->Vertex(0);
  vtxY = event->Vertex(1);
  vtxZ = event->Vertex(2);
  vtxN = event->VertexNContributors();
  diaX = event->DiamondY();
  diaY = event->DiamondX();
  diaZ = event->DiamondZ();
  totalTracks = event->NTracksTotal();
  centV0 = event->CentralityVZERO();
  centSPD = event->CentralitySPD();
  centTPC = event->CentralityTPC();
  centZEMvsZDC = event->CentralityZEMvsZDC();
  centV0A = event->CentralityVZEROA();
  centV0C = event->CentralityVZEROC();
  centZNA = event->CentralityZNA();
  centV0Mnew = event->CentralityV0MNew();
  centV0MnewPlus05 = event->CentralityV0MNewPlus05();
  centV0MnewMinus05 = event->CentralityV0MNewMinus05();
  centV0MnewPlus10 = event->CentralityV0MNewPlus10();
  centV0MnewMinus10 = event->CentralityV0MNewMinus10();
	centV0MPlus05 = event->CentralityV0MPlus05();
  centV0MMinus05 = event->CentralityV0MMinus05();
  centV0MPlus10 = event->CentralityV0MPlus10();
  centV0MMinus10 = event->CentralityV0MMinus10();

  // Pair Information
  pairP = pair->P();
  pairPt = pair->Pt();
  pairPhi = pair->Phi();
  pairTheta = pair->Theta();
  pairEta = pair->Eta();
  pairMass = pair->Mass();
  pairEnergy = pair->Energy();
  pairRapidity = pair->Rapidity();
  pairDecayRadius = pair->DecayRadius();
  pairPsProper = pair->PsProper();
  pairPointingAngle = pair->PointingAngle();
  pairChi2 = pair->Chi2();

  // Track Information Leg1
  leg1p = leg1->P();
  leg1pt = leg1->Pt();
  leg1phi = leg1->Phi();
  leg1theta = leg1->Theta();
  leg1eta = leg1->Eta();
  leg1charge = leg1->Charge();
  
  leg1massTracking = leg1->MassForTracking();
  leg1trackLength = leg1->TrackLength();
  leg1dcaXY = leg1->DCAxy();
  leg1dcaZ = leg1->DCAz();
  leg1helixX = leg1->HelixX();
  leg1helixY = leg1->HelixY();
  leg1helixR = leg1->HelixR();

  leg1itsnSigElec = leg1->ITSnSig(0);
  leg1itsnSigPion = leg1->ITSnSig(1);
  leg1itsnSigKaon = leg1->ITSnSig(2);
  leg1itsnSigProt = leg1->ITSnSig(3);
  leg1itsSignal = leg1->ITSsignal();
  leg1itsChi2 = leg1->ITSchi2();
  leg1itsNcls = leg1->ITSncls();
  leg1itsNclsShared = leg1->ITSnSharedCls();

  leg1tpcnSigElec = leg1->TPCnSig(0);
  leg1tpcnSigPion = leg1->TPCnSig(1);
  leg1tpcnSigKaon = leg1->TPCnSig(2);
  leg1tpcnSigProt = leg1->TPCnSig(3);
  leg1tpcSignal = leg1->TPCsignal();
  leg1tpcSignalTuned = leg1->TPCsignalTunedOnData();
  leg1tpcChi2 = leg1->TPCchi2();
  leg1tpcSignalN = leg1->TPCsignalN();
  leg1tpcNcls = leg1->TPCncls();
  leg1tpcNclsShared = leg1->TPCnclsShared();
  leg1tpcRows = leg1->TPCCrossedRows();

	leg1tpcdEdxQmaxIROC = leg1->TPCdEdxInfoQmax(0);
  leg1tpcdEdxQmaxMedOROC = leg1->TPCdEdxInfoQmax(1);
  leg1tpcdEdxQmaxLongOROC = leg1->TPCdEdxInfoQmax(2);
  leg1tpcdEdxQmaxAllOROC = leg1->TPCdEdxInfoQmax(3);
	leg1tpcdEdxQtotIROC = leg1->TPCdEdxInfoQtot(0);
  leg1tpcdEdxQtotMedOROC = leg1->TPCdEdxInfoQtot(1);
  leg1tpcdEdxQtotLongOROC = leg1->TPCdEdxInfoQtot(2);
  leg1tpcdEdxQtotAllOROC = leg1->TPCdEdxInfoQtot(3);
	leg1tpcActiveLength = leg1->TPCActiveLength();
  leg1tpcGeomLength = leg1->TPCGeomLength();
  
  leg1tofnSigElec = leg1->TOFnSig(0);
  leg1tofnSigPion = leg1->TOFnSig(1);
  leg1tofnSigKaon = leg1->TOFnSig(2);
  leg1tofnSigProt = leg1->TOFnSig(3);
  leg1tofBeta = leg1->TOFbeta();
  leg1tofTime = leg1->TOFtime();
  leg1tofMisProbab = leg1->TOFmismatchProbab();
  leg1tofChi2 = leg1->TOFchi2();
  
  leg1trdNtracklets = leg1->TRDntracklets(0);
  leg1trdNtrackletsPID = leg1->TRDntracklets(1);
  leg1trdSig1DElec = leg1->TRDpidLQ1D(0);
  leg1trdSig1DPion = leg1->TRDpidLQ1D(1);
  leg1trdSig2DElec = leg1->TRDpidLQ2D(0);
  leg1trdSig2DPion = leg1->TRDpidLQ2D(1);
  
  // Track Information Leg2
  leg2p = leg2->P();
  leg2pt = leg2->Pt();
  leg2phi = leg2->Phi();
  leg2theta = leg2->Theta();
  leg2eta = leg2->Eta();
  leg2charge = leg2->Charge();
  
  leg2massTracking = leg2->MassForTracking();
  leg2trackLength = leg2->TrackLength();
  leg2dcaXY = leg2->DCAxy();
  leg2dcaZ = leg2->DCAz();
  leg2helixX = leg2->HelixX();
  leg2helixY = leg2->HelixY();
  leg2helixR = leg2->HelixR();

  leg2itsnSigElec = leg2->ITSnSig(0);
  leg2itsnSigPion = leg2->ITSnSig(1);
  leg2itsnSigKaon = leg2->ITSnSig(2);
  leg2itsnSigProt = leg2->ITSnSig(3);
  leg2itsSignal = leg2->ITSsignal();
  leg2itsChi2 = leg2->ITSchi2();
  leg2itsNcls = leg2->ITSncls();
  leg2itsNclsShared = leg2->ITSnSharedCls();

  leg2tpcnSigElec = leg2->TPCnSig(0);
  leg2tpcnSigPion = leg2->TPCnSig(1);
  leg2tpcnSigKaon = leg2->TPCnSig(2);
  leg2tpcnSigProt = leg2->TPCnSig(3);
  leg2tpcSignal = leg2->TPCsignal();
  leg2tpcSignalTuned = leg2->TPCsignalTunedOnData();
  leg2tpcChi2 = leg2->TPCchi2();
  leg2tpcSignalN = leg2->TPCsignalN();
  leg2tpcNcls = leg2->TPCncls();
  leg2tpcNclsShared = leg2->TPCnclsShared();
  leg2tpcRows = leg2->TPCCrossedRows();

	leg2tpcdEdxQmaxIROC = leg2->TPCdEdxInfoQmax(0);
  leg2tpcdEdxQmaxMedOROC = leg2->TPCdEdxInfoQmax(1);
  leg2tpcdEdxQmaxLongOROC = leg2->TPCdEdxInfoQmax(2);
  leg2tpcdEdxQmaxAllOROC = leg2->TPCdEdxInfoQmax(3);
	leg2tpcdEdxQtotIROC = leg2->TPCdEdxInfoQtot(0);
  leg2tpcdEdxQtotMedOROC = leg2->TPCdEdxInfoQtot(1);
  leg2tpcdEdxQtotLongOROC = leg2->TPCdEdxInfoQtot(2);
  leg2tpcdEdxQtotAllOROC = leg2->TPCdEdxInfoQtot(3);
	leg2tpcActiveLength = leg2->TPCActiveLength();
  leg2tpcGeomLength = leg2->TPCGeomLength();
  
  leg2tofnSigElec = leg2->TOFnSig(0);
  leg2tofnSigPion = leg2->TOFnSig(1);
  leg2tofnSigKaon = leg2->TOFnSig(2);
  leg2tofnSigProt = leg2->TOFnSig(3);
  leg2tofBeta = leg2->TOFbeta();
  leg2tofTime = leg2->TOFtime();
  leg2tofMisProbab = leg2->TOFmismatchProbab();
  leg2tofChi2 = leg2->TOFchi2();
  
  leg2trdNtracklets = leg2->TRDntracklets(0);
  leg2trdNtrackletsPID = leg2->TRDntracklets(1);
  leg2trdSig1DElec = leg2->TRDpidLQ1D(0);
  leg2trdSig1DPion = leg2->TRDpidLQ1D(1);
  leg2trdSig2DElec = leg2->TRDpidLQ2D(0);
  leg2trdSig2DPion = leg2->TRDpidLQ2D(1);
}
