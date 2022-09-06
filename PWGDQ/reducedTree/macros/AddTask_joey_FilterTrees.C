/* 
 * takes a data file from gsi enriched dataset
 * or MC data from old dataset
 * and outputs the JPsi candidate pairs in it
 */

#include <iostream>
#include <fstream>
#include <string>

#include "TF1.h"

#define ISPP // if process is pp or else
//#def DO_RUNS // if runs should be counted

void Setup(AliReducedAnalysisFilterTrees* processor, TString prod, Bool_t isEnriched);
void SetupCuts(AliReducedAnalysisFilterTrees* processor, TString prod, Bool_t isEnriched);
void SetupHistogramManager(AliReducedAnalysisFilterTrees* task,  TString prod, Bool_t isEnriched);
void DefineHistograms(AliReducedAnalysisFilterTrees* task, TString prod, Bool_t isEnriched);

//__________________________________________________________________________________________
AliAnalysisTask* AddTask_joey_FilterTrees(Bool_t isAliRoot=kTRUE, Int_t runMode=1, Bool_t isMC = kFALSE, Bool_t isEnriched = kFALSE, TString prod="LHC10h") {    
   // isAliRoot={kTRUE for ESD/AOD analysis on the grid, kFALSE for local analysis on reduced trees}
   // runMode={AliAnalysisTaskReducedEventProcessor::kUseOnTheFlyReducedEvents=1, AliAnalysisTaskReducedEventProcessor::kUseEventsFromTree=2}

	TString className = "AddTask_joey_FilterTrees";

  cout << "Running " << className << ": (isAliRoot = " << isAliRoot << ", runMode = " << runMode << ", isMC = " << isMC << ", prod = " << prod << ")" << endl;

  AliReducedAnalysisFilterTrees* filterTask = new AliReducedAnalysisFilterTrees("FilterTrees","filter DST trees");
  filterTask->Init();
  filterTask->SetFilteredTreeWritingOption(AliReducedAnalysisTaskSE::kFullEventsWithFullTracks);
  filterTask->SetWriteFilteredTracks(kTRUE); // keep tracks for IsTrackSelected
  filterTask->SetWriteFilteredPairs(kTRUE); // keep candidates for IsPairSelected
  filterTask->SetBuildCandidatePairs(AliReducedPairInfo::kJpsiToEE); // set correct mass and symmetric conversion
  filterTask->SetBuildCandidateLikePairs(kFALSE);
  filterTask->SetRunCandidatePrefilter(kTRUE); // kill pion background
  filterTask->SetRunCandidatePrefilterOnSameCharge(kFALSE);
  filterTask->SetRejectEmptyPairs(kTRUE); // Do not save events with no pair candidates
  filterTask->SetRunOverMC(isMC);

  filterTask->UseTracks1(kTRUE); // from Ailec's AddTask -> Full Tracks
  filterTask->UseTracks2(kFALSE); // Base Tracks no usefull info
  filterTask->UseFullTracksOnly(kTRUE);
  filterTask->DebugEvery(1000); // -1 for no debug

  Setup(filterTask, prod, isEnriched);
  // initialize an AliAnalysisTask which will wrap the AliReducedAnalysisFilterTrees such that it can be run in an aliroot analysis train (e.g. LEGO, local analysis etc)
  AliAnalysisTaskReducedEventProcessor* task = new AliAnalysisTaskReducedEventProcessor("ReducedEventAnalysisManager", runMode, kTRUE); // (name, runMode, writeTrees)
  task->AddTask(filterTask);
  
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) { Error(className, "No analysis manager found."); return 0; }
  
  AliAnalysisDataContainer* cReducedEvent = NULL;
  switch(runMode) {
  case AliAnalysisTaskReducedEventProcessor::kUseOnTheFlyReducedEvents:
    cReducedEvent = (AliAnalysisDataContainer*) mgr->GetContainers()->FindObject("ReducedEventDQ");
    break;
  case AliAnalysisTaskReducedEventProcessor::kUseEventsFromTree:
    cReducedEvent = mgr->GetCommonInputContainer();
    break;
  default:
    Error(className, "No exchange container with ReducedEvent found.");
    return 0;
	}
	
	mgr->AddTask(task);
  mgr->ConnectInput(task, 0, cReducedEvent);

  AliAnalysisDataContainer *cOutputHist = mgr->CreateContainer("Histos", THashList::Class(), AliAnalysisManager::kOutputContainer, Form("AnalysisHistograms_Jpsi2ee_%s.root", isMC ? "MC" : "data"));      
  mgr->ConnectOutput(task, 1, cOutputHist );

  AliAnalysisDataContainer *cOutputTree = mgr->CreateContainer("filteredDstTree", TTree::Class(), AliAnalysisManager::kOutputContainer, Form("JpsiCandidates_%s.root", isMC ? "MC" : "data"));
  mgr->ConnectOutput(task, 2, cOutputTree );
  
  return task;
}


//_________________________________________________________________
void Setup(AliReducedAnalysisFilterTrees* processor, TString prod /*="LHC10h"*/, Bool_t isEnriched) {
  SetupCuts(processor, prod, isEnriched);
  SetupHistogramManager(processor, prod, isEnriched);
}
  
void SetupCuts(AliReducedAnalysisFilterTrees* processor, TString prod /*="LHC10h"*/, Bool_t isEnriched) {
  // Configure cuts
  bool isMC = processor->GetRunOverMC();
  cout << "isMC : " << isMC << endl;

  // Set event cuts
  AliReducedEventCut* evCut1 = new AliReducedEventCut("Centrality","Centrality selection"); // in one AliCut ALL added cuts have to be passed
  evCut1->AddCut(AliReducedVarManager::kVtxZ, -10.0, 10.0); // vertex centrality
  evCut1->AddCut(AliReducedVarManager::kIsPhysicsSelection, 0.1, 2.); // request physics selection
  if (isEnriched) evCut1->AddEventTagFilterBit(14); // gsi enriched data: needs event->EventTag(14) for MB
  evCut1->AddCut(AliReducedVarManager::kINT7Triggered, 0.1, 2.);
  //evCut1->AddCut(AliReducedVarManager::kHighMultV0Triggered, -0.1, 0.1);
  //evCut1->AddCut(AliReducedVarManager::kHighMultSPDTriggered, -0.1, 0.1);
  processor->AddEventCut(evCut1);

  // Set basic track cuts
  AliReducedTrackCut* standardCut = new AliReducedTrackCut("standard","");
  standardCut->AddCut(AliReducedVarManager::kP, .1,100.);
  standardCut->AddCut(AliReducedVarManager::kEta, -0.9,0.9);
  standardCut->AddCut(AliReducedVarManager::kDcaXY, -1.0,1.0);
  standardCut->AddCut(AliReducedVarManager::kDcaZ, -3.0,3.0);
  standardCut->AddCut(AliReducedVarManager::kTPCncls, 70.,160.0);
  standardCut->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kElectron, -3.0, 3.0); // accept electrons
  //standardCut->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kPion, -3.0, 3.0, kTRUE); // exclude pions
  //standardCut->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kKaon, -3.0, 3.0, kTRUE); // exclude kaons
  //standardCut->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kProton, -3.0, 3.0, kTRUE); // exclude protons
  standardCut->SetRejectKinks();
  standardCut->SetRequestITSrefit();
  standardCut->SetRequestTPCrefit();
  //standardCut->SetRequestTOFout();
  standardCut->SetRequestSPDany();
  //standardCut->SetRejectPureMC(kTRUE); // YES keep mc data...
  //standardCut->SetTrackFilterBit(0); // full track info, fQualityFlag bit 32 // disabled, because unknown in used trees
  standardCut->AddCut(AliReducedVarManager::kTPCchi2, 0., 4.0);
  processor->AddTrackCut(standardCut); 
  processor->AddCandidateLeg1Cut(standardCut); // need just Leg1, because ee-Pair is symmetric // IsCandidateLegSelected()
  
  
  // Setup Pion background rejection prefilter
  AliReducedTrackCut* prefTrackCut1 = new AliReducedTrackCut("prefCutPt07","prefilter P selection");
  prefTrackCut1->AddCut(AliReducedVarManager::kP, 0.2,100.0); // TODO why?
  prefTrackCut1->SetRequestTPCrefit();
  //prefTrackCut1->SetRejectPureMC(kTRUE);
  processor->AddCandidateLeg1PrefilterCut(prefTrackCut1);  
  
  AliReducedVarCut* prefPairCut = new AliReducedVarCut("prefCutM50MeV","prefilter pair cuts");
  prefPairCut->AddCut(AliReducedVarManager::kMass, 0.0, 0.05, kTRUE); // kill possible pion conversions -> reduce background
  processor->AddCandidateLeg1PairPrefilterCut(prefPairCut);


  
  AliReducedTrackCut* pairCut = new AliReducedTrackCut("JpsiCut","Pt pair selection");
  pairCut->AddCut(AliReducedVarManager::kPt, 0.0,100.0);
  pairCut->AddCut(AliReducedVarManager::kRap, -0.9,0.9);
  processor->AddCandidatePairCut(pairCut);


	// Select one of the following cuts

  AliReducedTrackCut* jpsiMCCut = new AliReducedTrackCut("jpsiCutMC","Rapidity JpsiSelection");
  jpsiMCCut->AddCut(AliReducedVarManager::kPt, 0.0, 100.0);
  jpsiMCCut->AddCut(AliReducedVarManager::kRap, -0.9, 0.9);
  jpsiMCCut->AddCut(AliReducedVarManager::kPdgMC, 442.5, 443.5); // select JPsi Particles only
  //jpsiMCCut->AddCut(AliReducedVarManager::kPdgMC+1, 442.5, 443.5, kTRUE);
  //jpsiMCCut->SetRejectPureMC(kFALSE);
  AliReducedTrackCut* electronMCCut = new AliReducedTrackCut("standardMC","Pt electron Selection");
  processor->AddJpsiMotherMCCut(jpsiMCCut,electronMCCut);
  
  AliReducedTrackCut* jpsiMCCutEta = (AliReducedTrackCut*) jpsiMCCut->Clone("jpsiCutMCEta");
  AliReducedTrackCut* electronMCCutEta = new AliReducedTrackCut("standardMCEta","Eta electron Selection");
  electronMCCutEta->AddCut(AliReducedVarManager::kEta, -0.9, 0.9);
  //processor->AddJpsiMotherMCCut(jpsiMCCutEta,electronMCCutEta);


  AliReducedTrackCut* jpsiMCCutEtaPt = (AliReducedTrackCut*) jpsiMCCut->Clone("jpsiCutMCEtaPt");
  AliReducedTrackCut* electronMCCutEtaPt = new AliReducedTrackCut("standardMCEtaPt","Eta and Pt electron Selection");
  electronMCCutEtaPt->AddCut(AliReducedVarManager::kEta, -0.9, 0.9);
  electronMCCutEtaPt->AddCut(AliReducedVarManager::kPt, 1., 1e5);
  //processor->AddJpsiMotherMCCut(jpsiMCCutEtaPt,electronMCCutEtaPt);
}


//_________________________________________________________________
void SetupHistogramManager(AliReducedAnalysisFilterTrees* task, TString prod /*="LHC10h"*/, Bool_t isEnriched) {
  //
  // setup the histograms manager
  //
  AliReducedVarManager::SetDefaultVarNames();
  
  DefineHistograms(task, prod, isEnriched);
  
  AliReducedVarManager::SetUseVars(task->GetHistogramManager()->GetUsedVars());
  AliReducedVarManager::SetUseVariable(AliReducedVarManager::kRunID);
}

//_________________________________________________________________
void DefineHistograms(AliReducedAnalysisFilterTrees* task, TString prod /*="LHC10h"*/, Bool_t isEnriched) {
  //
  // define histograms
  // NOTE: The DefineHistograms need to be called after the track cuts are defined because the name of the cuts are used in the histogram lists
  // NOTE: The name of the pair histogram lists for event mixing need to contain "PairMEPP", "PairMEPM", "PairMEMM", in their name, see below.
  
  // TODO: make needed changes such that this becomes less prone to mistakes
  
  
  #ifdef DO_RUNS
  
  #ifdef ISPP
  const char* fileRuns = "/home/glegras/alice/data/runLists/runsMC.txt";
  const char* fileRuns2 = "/gluster1/glegras/runsMC.txt";
  #else
  const char* fileRuns = "/home/glegras/alice/data/runLists/runspPb.txt";
  const char* fileRuns2 = "/gluster1/glegras/runspPb.txt";
  #endif // ISPP
  
  int nRuns = 0;
  TString runNumbers = "";  fileRuns; const char* fileRuns2;
  fstream newfile;
  newfile.open(fileRuns,ios::in);
  if(!(newfile.is_open())) newfile.open(fileRuns2,ios::in);
  if(newfile.is_open()) {
    string temp;
    while(getline(newfile, temp)){
      runNumbers+=temp;
      runNumbers+=";";
      nRuns++;
    }
    newfile.close(); 
  }
  AliReducedVarManager::SetRunNumbers(runNumbers);
  #endif // DO_RUNS


  AliHistogramManager* man = task->GetHistogramManager(); 
  bool isMC=task->GetRunOverMC();
   
  TString histClasses = "";
  histClasses += "Event_BeforeCuts;"; // FillEventInfo
  histClasses += "Event_AfterCuts;";

  histClasses += "EventTag_BeforeCuts;"; // FillEventTagInput
  histClasses += "EventTag_AfterCuts;"; 
    
  histClasses += "EventTriggers_BeforeCuts;"; // FillEventOnlineTrigger
  histClasses += "EventTriggers_AfterCuts;";   
  
  if(task->GetWriteFilteredTracks()) {
    histClasses += "Track_BeforeCuts;"; // FillTrackInfo + FillClusterMatchedTrackInfo
    for(Int_t i=0; i<task->GetNTrackCuts(); ++i) {
      TString cutName = task->GetTrackCutName(i);
      histClasses += Form("Track_%s;", cutName.Data());
    }
  }
  
  if(task->GetWriteFilteredPairs()) {
    histClasses += "Pair_BeforeCuts;"; // FillPairInfo
    histClasses += "PairQualityFlags_BeforeCuts;"; // FillPairQualityFlag
    for(Int_t i=0; i<task->GetNPairCuts(); ++i) {
      TString cutName = task->GetPairCutName(i);
      histClasses += Form("Pair_%s;", cutName.Data());
      histClasses += Form("PairQualityFlags_%s;", cutName.Data());
    }
  }
  
  if(task->GetBuildCandidatePairs()) { // FillTrackInfo + FillClusterMatchedTrackInfo
    for(Int_t i=0; i<task->GetNCandidateLegCuts(); ++i) {
      histClasses += Form("Track_LEG1_BeforePrefilter_%s;", task->GetCandidateLegCutName(i,1));
      histClasses += Form("Track_LEG2_BeforePrefilter_%s;", task->GetCandidateLegCutName(i,2));
      if(task->GetRunCandidatePrefilter()) {
        histClasses += Form("Track_LEG1_AfterPrefilter_%s;", task->GetCandidateLegCutName(i,1));
        histClasses += Form("Track_LEG2_AfterPrefilter_%s;", task->GetCandidateLegCutName(i,2));
      }
      histClasses += Form("Pair_Candidate12_%s%s;", 
			  task->GetCandidateLegCutName(i,1), (task->IsAsymmetricDecayChannel() ? Form("_%s", task->GetCandidateLegCutName(i,2)) : ""));
      if(task->GetBuildCandidateLikePairs()) {
        	histClasses += Form("Pair_Candidate11_%s%s;", 
			    task->GetCandidateLegCutName(i,1), (task->IsAsymmetricDecayChannel() ? Form("_%s", task->GetCandidateLegCutName(i,2)) : ""));
	        histClasses += Form("Pair_Candidate22_%s%s;", 
			    task->GetCandidateLegCutName(i,1), (task->IsAsymmetricDecayChannel() ? Form("_%s", task->GetCandidateLegCutName(i,2)) : ""));
      }
    }
    if(task->GetRunCandidatePrefilter()) {
        histClasses += "Track_LEG1_PrefilterTrack;";
        histClasses += "Track_LEG2_PrefilterTrack;";
		}
  }

  if(isMC) {
    for (int i=0; i<task->GetNJpsiMotherMCCuts();i++) 
      histClasses += Form("PureMCTRUTH_AfterSelection_%s;", task->GetJpsiMotherMCcutName(i));
    
    //histClasses += Form("PureMCTRUTH_DetectedDaughters_%s;", task->GetJpsiMotherMCcutName(task->GetNJpsiMotherMCCuts()-1));
  }

	#ifdef ISPP
  int nMultBins = 201; float minMult = -0.5; float maxMult = 200.5;
  Int_t runNBins = 42800;
  Double_t runHistRange[2] = {252200.,295000.};
  #else // maybe use elif
  int nMultBins = 301; float minMult = -0.5; float maxMult = 300.5;
  Int_t runNBins = 300;
  Double_t runHistRange[2] = {265300.,265600.};
  #endif
  
  if (!prod.CompareTo("LHC10h")) {runNBins = 2500; runHistRange[0] = 137100.; runHistRange[1] = 139600.;} // Pb-Pb of 2010 // default
  if (!prod.CompareTo("LHC11h")) {runNBins = 2700; runHistRange[0] = 167900.; runHistRange[1] = 170600.;} // Pb-Pb of 2011 
  if (!prod.CompareTo("LHC15o")) {runNBins = 2100; runHistRange[0] = 244900.; runHistRange[1] = 247000.;} // Pb-Pb of 2015
  if (!prod.CompareTo("LHC13b")) {runNBins =  400; runHistRange[0] = 195300.; runHistRange[1] = 195700.;} // p-Pb of 2013
  if (!prod.CompareTo("LHC13c")) {runNBins =  400; runHistRange[0] = 195300.; runHistRange[1] = 195700.;} // p-Pb of 2013 part 2
  if (!prod.CompareTo("LHC16l")) {runNBins = 1140; runHistRange[0] = 258880.; runHistRange[1] = 260020.;} // pp at 13 TeV
  if (!prod.CompareTo("LHC16r")) {runNBins = 1000; runHistRange[0] = 265400.; runHistRange[1] = 266400.;} // p-Pb at 8.16 TeV
  
  const Int_t kNBCBins = 3600;
  Double_t bcHistRange[2] = {-0.5,3599.5};
  
  Double_t massBinWidth = 0.01;     // *GeV/c^2
  Double_t massRange[2] = {0.0,5.0};
  Int_t nMassBins = TMath::Nint((massRange[1]-massRange[0])/massBinWidth);
  
  TString candidateNames = "#gamma#rightarrow e^{+}e^{-};K^{0}_{S}#rightarrow#pi^{+}#pi^{-};";
  candidateNames += "#Lambda#rightarrow p#pi^{-};#bar{#Lambda}#rightarrow #bar{p}#pi^{+};";
  
  TString classesStr(histClasses);
  TObjArray* arr=classesStr.Tokenize(";");
  
  cout << "\nHistogram classes included in the Histogram Manager:" << endl;
  for(Int_t iclass=0; iclass<arr->GetEntries(); ++iclass) {
    TString classStr = arr->At(iclass)->GetName();
    
    // Event wise histograms
    if (classStr.Contains("EventTag_")) {
      man->AddHistClass(classStr.Data());
      cout << "  " << classStr.Data() << endl;
      TString tagNames = "";
      tagNames += "AliAnalysisUtils 2013 selection;";
      tagNames += "AliAnalysisUtils MV pileup;";
      tagNames += "AliAnalysisUtils MV pileup, no BC check;";
      tagNames += "AliAnalysisUtils MV pileup, min wght dist 10;";
      tagNames += "AliAnalysisUtils MV pileup, min wght dist 5;";
      tagNames += "IsPileupFromSPD(3,0.6,3.,2.,5.);";
      tagNames += "IsPileupFromSPD(4,0.6,3.,2.,5.);";
      tagNames += "IsPileupFromSPD(5,0.6,3.,2.,5.);";
      tagNames += "IsPileupFromSPD(6,0.6,3.,2.,5.);";
      tagNames += "IsPileupFromSPD(3,0.8,3.,2.,5.);";
      tagNames += "IsPileupFromSPD(4,0.8,3.,2.,5.);";
      tagNames += "IsPileupFromSPD(5,0.8,3.,2.,5.);";
      tagNames += "IsPileupFromSPD(6,0.8,3.,2.,5.);";
      tagNames += "vtx distance selected;";
      man->AddHistogram(classStr.Data(), "EventTags", "Event tags", kFALSE,
			20, -0.5, 19.5, AliReducedVarManager::kEventTag, 0, 0.0, 0.0, AliReducedVarManager::kNothing, 0, 0.0, 0.0, AliReducedVarManager::kNothing, tagNames.Data());
      man->AddHistogram(classStr.Data(), "EventTags_CentVZERO", "Event tags vs VZERO centrality", kFALSE,
			20, -0.5, 19.5, AliReducedVarManager::kEventTag, 9, 0.0, 90.0, AliReducedVarManager::kCentVZERO, 0, 0.0, 0.0, AliReducedVarManager::kNothing, tagNames.Data());
      continue;
    }
    
    if (classStr.Contains("EventTriggers_")) {
      man->AddHistClass(classStr.Data());
      cout << "  " << classStr.Data() << endl;
      TString triggerNames = "";
      for(Int_t i=0; i<64; ++i) {triggerNames += AliReducedVarManager::fgkOfflineTriggerNames[i]; triggerNames+=";";}
       
      man->AddHistogram(classStr.Data(), "Triggers2", "", kFALSE,
			64, -0.5, 63.5, AliReducedVarManager::kOnlineTriggerFired2, 0, 0.0, 0.0, AliReducedVarManager::kNothing, 0, 0.0, 0.0, AliReducedVarManager::kNothing, triggerNames.Data());
      man->AddHistogram(classStr.Data(), "CentVZERO_Triggers2", "", kFALSE,
			64, -0.5, 63.5, AliReducedVarManager::kOnlineTriggerFired2, 9, 0.0, 90.0, AliReducedVarManager::kCentVZERO, 0, 0.0, 0.0, AliReducedVarManager::kNothing, triggerNames.Data());
      continue;
    }
    
    if (classStr.Contains("Event_")) {
      man->AddHistClass(classStr.Data());
      cout << "  " << classStr.Data() << endl;
      man->AddHistogram(classStr.Data(),"RunNo","Run numbers",kFALSE, runNBins, runHistRange[0], runHistRange[1], AliReducedVarManager::kRunNo);
      #ifdef DO_RUNS
      man->AddHistogram(classStr.Data(),"RunID","Run ID",kFALSE, nRuns, 0, nRuns, AliReducedVarManager::kRunID);
      #endif
      man->AddHistogram(classStr.Data(),"TimeFromSOR","Events vs time",kFALSE, 450, 0., 450., AliReducedVarManager::kTimeRelativeSOR);
      man->AddHistogram(classStr.Data(),"VtxX","Vtx X",kFALSE,300,-0.4,0.4,AliReducedVarManager::kVtxX);
      man->AddHistogram(classStr.Data(),"VtxX_TimeFromSOR_prof","<Vtx X> vs time from SOR",kTRUE, 90, 0., 450., AliReducedVarManager::kTimeRelativeSOR, 300,-0.4,0.4,AliReducedVarManager::kVtxX);
      man->AddHistogram(classStr.Data(),"VtxY","Vtx Y",kFALSE,300,-0.4,0.4,AliReducedVarManager::kVtxY);
      man->AddHistogram(classStr.Data(),"VtxY_TimeFromSOR_prof","<Vtx Y> vs time from SOR",kTRUE, 90, 0., 450., AliReducedVarManager::kTimeRelativeSOR, 300,-0.4,0.4,AliReducedVarManager::kVtxY);
      man->AddHistogram(classStr.Data(),"VtxZ","Vtx Z",kFALSE,300,-15.,15.,AliReducedVarManager::kVtxZ);
      man->AddHistogram(classStr.Data(),"VtxZ_TimeFromSOR_prof","<Vtx Z> vs time from SOR",kTRUE, 90, 0., 450., AliReducedVarManager::kTimeRelativeSOR, 300,-15.0,15.0,AliReducedVarManager::kVtxZ);
      man->AddHistogram(classStr.Data(),"CentVZERO","Centrality(VZERO)",kFALSE, 100, 0.0, 100.0, AliReducedVarManager::kCentVZERO);
      man->AddHistogram(classStr.Data(),"CentVZERO_VtxZ","Centrality(VZERO) vs vtxZ",kFALSE, 50, 0.0, 100.0, AliReducedVarManager::kCentVZERO, 50, -10., 10., AliReducedVarManager::kVtxZ);
      man->AddHistogram(classStr.Data(),"CentVZERO_TimeFromSOR","Centrality(VZERO) vs time from SOR",kFALSE, 90, 0., 450., AliReducedVarManager::kTimeRelativeSOR, 50, 0.0, 100.0, AliReducedVarManager::kCentVZERO);
      man->AddHistogram(classStr.Data(),"CentSPD","Centrality(SPD)",kFALSE, 100, 0.0, 100.0, AliReducedVarManager::kCentSPD);
      man->AddHistogram(classStr.Data(),"CentSPD_VtxZ","Centrality(SPD) vs vtxZ",kFALSE, 50, 0.0, 100.0, AliReducedVarManager::kCentSPD, 50, -10., 10., AliReducedVarManager::kVtxZ);
      man->AddHistogram(classStr.Data(),"CentTPC","Centrality(TPC)",kFALSE, 100, 0.0, 100.0, AliReducedVarManager::kCentTPC);
      man->AddHistogram(classStr.Data(),"CentTPC_VtxZ","Centrality(TPC) vs vtxZ",kFALSE, 100, 0.0, 100.0, AliReducedVarManager::kCentTPC, 50, -10., 10., AliReducedVarManager::kVtxZ);
      man->AddHistogram(classStr.Data(),"CentQuality","Centrality quality",kFALSE, 500, 0., 500., AliReducedVarManager::kCentQuality);
      man->AddHistogram(classStr.Data(),"IRIntClosestInt1Map","Closest out of bunch interactions Int1",kFALSE,7000,-3500.,3500.,AliReducedVarManager::kIRIntClosestIntMap);
      man->AddHistogram(classStr.Data(),"IRIntClosestInt2Map","Closest out of bunch interactions Int2",kFALSE,7000,-3500.,3500.,AliReducedVarManager::kIRIntClosestIntMap+1);
      man->AddHistogram(classStr.Data(),"NTracksTotal","Number of total tracks per event",kFALSE,500,0.,30000.,AliReducedVarManager::kNtracksTotal);
      man->AddHistogram(classStr.Data(),"NTracksTotal_BeamIntensity0_prof","Number of total tracks per event",kTRUE,100, 3.0e+12, 9.0e+12, AliReducedVarManager::kBeamIntensity0, 500,0.,20000.,AliReducedVarManager::kNtracksTotal);

      man->AddHistogram(classStr.Data(),"NTracksSelected_TimeFromSOR","Averaged number of selected tracks per event vs time from SOR",kTRUE, 90, 0., 450., AliReducedVarManager::kTimeRelativeSOR, 100,0.,100.,AliReducedVarManager::kNtracksSelected);
      man->AddHistogram(classStr.Data(),"NTracksSelected_CentVZERO_TimeFromSOR","Averaged number of selected tracks per event per centrality vs time from SOR",kTRUE, 20, 0.0, 100.0, AliReducedVarManager::kCentVZERO, 90, 0., 450., AliReducedVarManager::kTimeRelativeSOR, 100,0.,100.,AliReducedVarManager::kNtracksSelected);
      man->AddHistogram(classStr.Data(),"DeltaVtxZ","Z_{global}-Z_{TPC}",kFALSE,300,-1.,1.,AliReducedVarManager::kDeltaVtxZ);
      man->AddHistogram(classStr.Data(),"DeltaVtxZ_TimeFromSOR_prof","Z_{global}-Z_{TPC} vs time from SOR",kTRUE,90, 0.0, 450.,  AliReducedVarManager::kTimeRelativeSOR, 300,-1.,1.,AliReducedVarManager::kDeltaVtxZ);
      man->AddHistogram(classStr.Data(),"NVtxTPCContributors","",kFALSE,500,0.,10000.,AliReducedVarManager::kNVtxTPCContributors);
      man->AddHistogram(classStr.Data(),"NVtxTPCContributors_BeamIntensity0","",kTRUE, 100, 3.0e+12, 9.0e+12, AliReducedVarManager::kBeamIntensity0, 500,0.,10000.,AliReducedVarManager::kNVtxTPCContributors);
      man->AddHistogram(classStr.Data(),"SPDntracklets", "SPD #tracklets in |#eta|<1.0", kFALSE, 200, 0., 5000., AliReducedVarManager::kSPDntracklets);
      man->AddHistogram(classStr.Data(),"SPDntracklets_TimeFromSOR_prof", "SPD <#tracklets> in |#eta|<1.0 vs time from SOR", kTRUE, 90, 0., 450., AliReducedVarManager::kTimeRelativeSOR, 200, 0., 5000., AliReducedVarManager::kSPDntracklets);
      //Bool_t isCalibrated = kTRUE;
      man->AddHistogram(classStr.Data(), "VZEROA_NEmptyChannels_VtxCent_prof", "No. VZERO-A empty channels per event vs. centrality SPD and vertex Z", kTRUE,
			24, -12.0, +12.0, AliReducedVarManager::kVtxZ, 20, 0.0, 100., AliReducedVarManager::kCentSPD, 32, 0.0, 32.0, AliReducedVarManager::kVZEROAemptyChannels);
      man->AddHistogram(classStr.Data(), "VZEROC_NEmptyChannels_VtxCent_prof", "No. VZERO-C empty channels per event vs. centrality SPD and vertex Z", kTRUE,
			24, -12.0, +12.0, AliReducedVarManager::kVtxZ, 20, 0.0, 100., AliReducedVarManager::kCentSPD, 32, 0.0, 32.0, AliReducedVarManager::kVZEROCemptyChannels);
      continue;
    }  // end if className contains "Event" 

    // Track histograms
    if (classStr.Contains("TrackITSclusterMap_")) {
      man->AddHistClass(classStr.Data());
      cout << "  " << classStr.Data() << endl;
      man->AddHistogram(classStr.Data(), "ITSlayerHit", "Hits in the ITS layers", kFALSE, 6, 0.5, 6.5, AliReducedVarManager::kITSlayerHit);
      man->AddHistogram(classStr.Data(), "ITSlayerHit_Phi", "Hits in the ITS layers vs #varphi", kFALSE, 180, 0.0, 6.29, AliReducedVarManager::kPhi, 6, 0.5, 6.5, AliReducedVarManager::kITSlayerHit);
      man->AddHistogram(classStr.Data(), "ITSlayerHit_Eta", "Hits in the ITS layers vs #eta", kFALSE, 100, -1.0, 1.0, AliReducedVarManager::kEta, 6, 0.5, 6.5, AliReducedVarManager::kITSlayerHit);
      //man->AddHistogram(classStr.Data(), "ITSlayerHit_DCAxy", "Hits in the ITS layers vs DCA_{xy}", kFALSE, 1000, -0.5, 0.5, AliReducedVarManager::kDcaXY, 6, 0.5, 6.5, AliReducedVarManager::kITSlayerHit);
      //man->AddHistogram(classStr.Data(), "ITSlayerHit_DCAz", "Hits in the ITS layers vs DCA_{z}", kFALSE, 1800, -1.0, 1.0, AliReducedVarManager::kDcaZ, 6, 0.5, 6.5, AliReducedVarManager::kITSlayerHit);
      continue;
    }  // end of ITSclusterMap histogram definitions
    
    if (classStr.Contains("TrackTPCclusterMap_")) {
      man->AddHistClass(classStr.Data());
      cout << "  " << classStr.Data() << endl;
      man->AddHistogram(classStr.Data(), "TPCclusterMap", "TPC cluster map", kFALSE,
			8, -0.5, 7.5, AliReducedVarManager::kTPCclusBitFired);
      man->AddHistogram(classStr.Data(), "TPCclusterMap_Phi", "TPC cluster map vs #varphi", kFALSE,
			180, 0.0, 6.29, AliReducedVarManager::kPhi, 8, -0.5, 7.5, AliReducedVarManager::kTPCclusBitFired);
      man->AddHistogram(classStr.Data(), "TPCclusterMap_Eta", "TPC cluster map vs #eta", kFALSE,
			100, -1.0, 1.0, AliReducedVarManager::kEta, 8, -0.5, 7.5, AliReducedVarManager::kTPCclusBitFired);
      man->AddHistogram(classStr.Data(), "TPCclusterMap_Pt", "TPC cluster map vs p_{T}", kFALSE,
			100, 0.0, 10.0, AliReducedVarManager::kPt, 8, -0.5, 7.5, AliReducedVarManager::kTPCclusBitFired);
      continue;
    }  // end of TPCclusterMap histogram definitions
    
    TString trkStatusNames = "";
    for(Int_t iflag=0; iflag<AliReducedVarManager::kNTrackingStatus; ++iflag) {
      trkStatusNames += AliReducedVarManager::fgkTrackingStatusNames[iflag];
      trkStatusNames += ";";
    }
    if (classStr.Contains("TrackStatusFlags_")) {
      man->AddHistClass(classStr.Data());
      cout << "  " << classStr.Data() << endl;
      man->AddHistogram(classStr.Data(), "TrackingFlags", "Tracking flags;;", kFALSE,
			AliReducedVarManager::kNTrackingStatus, -0.5, AliReducedVarManager::kNTrackingStatus-0.5, AliReducedVarManager::kTrackingFlag, 
			0, 0.0, 0.0, AliReducedVarManager::kNothing, 0, 0.0, 0.0, AliReducedVarManager::kNothing, trkStatusNames.Data());
      /*man->AddHistogram(classStr.Data(),"TPCnsigElectronCorrected_Pin_TrackingFlags","Corrected TPC N_{#sigma} electron vs. inner param P vs tracking flags;;",kFALSE,
	50,0.0,5.0,AliReducedVarManager::kPin,50,-5.0,5.0,AliReducedVarManager::kTPCnSigCorrected+AliReducedVarManager::kElectron, AliReducedVarManager::kNTrackingStatus, -0.5, AliReducedVarManager::kNTrackingStatus-0.5, AliReducedVarManager::kTrackingFlag, "", "", trkStatusNames.Data());
      man->AddHistogram(classStr.Data(),"TPCncls_Eta_NTracksTPCoutFromPileup_TrackingFlags_prof","",kTRUE, 20,-1.0,1.0,AliReducedVarManager::kEta,
	43,-1500., 20000., AliReducedVarManager::kNTracksTPCoutFromPileup, AliReducedVarManager::kNTrackingStatus, -0.5, AliReducedVarManager::kNTrackingStatus-0.5, AliReducedVarManager::kTrackingFlag, "", "", trkStatusNames.Data(),  AliReducedVarManager::kTPCncls);*/
      continue;
    }
    
    if (classStr.Contains("Track_")) {
      man->AddHistClass(classStr.Data());
      cout << "  " << classStr.Data() << endl;
      man->AddHistogram(classStr.Data(), "P_Pin", "P vs Pin", kFALSE, 200, 0.0, 10.0, AliReducedVarManager::kP, 200, 0.0, 10.0, AliReducedVarManager::kPin);
      man->AddHistogram(classStr.Data(), "Pt", "p_{T} distribution", kFALSE, 1000, 0.0, 50.0, AliReducedVarManager::kPt);
      man->AddHistogram(classStr.Data(), "Pt_TimeFromSOR", "<p_{T}> vs time from SOR", kTRUE, 90, 0., 450., AliReducedVarManager::kTimeRelativeSOR, 1000, 0.0, 50.0, AliReducedVarManager::kPt);
      man->AddHistogram(classStr.Data(), "Eta", "", kFALSE, 1000, -1.5, 1.5, AliReducedVarManager::kEta);
      man->AddHistogram(classStr.Data(), "Eta_Phi", "", kFALSE, 100, -1.0, 1.0, AliReducedVarManager::kEta, 100, 0., 6.3, AliReducedVarManager::kPhi);
      man->AddHistogram(classStr.Data(), "Eta_TimeFromSOR", "<#eta> vs time from SOR", kTRUE, 90, 0., 450., AliReducedVarManager::kTimeRelativeSOR, 1000, -1.5, 1.5, AliReducedVarManager::kEta);
      man->AddHistogram(classStr.Data(), "Phi", "", kFALSE, 1000, 0.0, 6.3, AliReducedVarManager::kPhi);
      man->AddHistogram(classStr.Data(), "Phi_TimeFromSOR", "<#varphi> vs time from SOR", kTRUE, 90, 0., 450., AliReducedVarManager::kTimeRelativeSOR, 1000, 0.0, 6.3, AliReducedVarManager::kPhi);
      man->AddHistogram(classStr.Data(), "DCAxy_Pt", "DCAxy", kFALSE, 100, -2.0, 2.0, AliReducedVarManager::kDcaXY, 50, 0.0, 5.0, AliReducedVarManager::kPt);
      man->AddHistogram(classStr.Data(), "DCAxy_TPCchi2", "DCAxy vs TPC chi2", kFALSE, 100, -2.0, 2.0, AliReducedVarManager::kDcaXY, 50, 0.0, 5.0, AliReducedVarManager::kTPCchi2);
      man->AddHistogram(classStr.Data(), "DCAz_TPCchi2", "DCAz vs TPC chi2", kFALSE, 100, -5.0, 5.0, AliReducedVarManager::kDcaZ, 50, 0.0, 5.0, AliReducedVarManager::kTPCchi2);
      man->AddHistogram(classStr.Data(), "DCAxy_ITSchi2", "DCAxy vs ITS chi2", kFALSE, 100, -2.0, 2.0, AliReducedVarManager::kDcaXY, 100, 0.0, 25.0, AliReducedVarManager::kITSchi2);
      man->AddHistogram(classStr.Data(), "DCAz_ITSchi2", "DCAz vs ITS chi2", kFALSE, 100, -5.0, 5.0, AliReducedVarManager::kDcaZ, 100, 0.0, 25.0, AliReducedVarManager::kITSchi2);
      man->AddHistogram(classStr.Data(), "DCAxy_goldenChi2", "DCAxy vs golden chi2", kFALSE, 100, -2.0, 2.0, AliReducedVarManager::kDcaXY, 100, 0.0, 100.0, AliReducedVarManager::kChi2TPCConstrainedVsGlobal);
      man->AddHistogram(classStr.Data(), "DCAz_goldenChi2", "DCAz vs golden chi2", kFALSE, 100, -5.0, 5.0, AliReducedVarManager::kDcaZ, 100, 0.0, 100.0, AliReducedVarManager::kChi2TPCConstrainedVsGlobal);
      man->AddHistogram(classStr.Data(), "DCAxy", "DCAxy", kFALSE, 200, -5.0, 5.0, AliReducedVarManager::kDcaXY);
      man->AddHistogram(classStr.Data(), "DCAxy_TimeFromSOR", "<DCAxy> vs time from SOR", kTRUE, 90, 0., 450., AliReducedVarManager::kTimeRelativeSOR, 200, -5.0, 5.0, AliReducedVarManager::kDcaXY);
      man->AddHistogram(classStr.Data(), "DCAz", "DCAz", kFALSE, 200, -5.0, 5.0, AliReducedVarManager::kDcaZ);
      man->AddHistogram(classStr.Data(), "DCAz_Pt", "DCAz", kFALSE, 100, -3.0, 3.0, AliReducedVarManager::kDcaZ, 50, 0.0, 5.0, AliReducedVarManager::kPt);
      man->AddHistogram(classStr.Data(), "DCAz_TimeFromSOR", "<DCAz> vs time from SOR", kTRUE, 90, 0., 450., AliReducedVarManager::kTimeRelativeSOR, 200, -5.0, 5.0, AliReducedVarManager::kDcaZ);
      man->AddHistogram(classStr.Data(),"DCAxy_Eta_MultVZERO_NTracksTPCoutFromPileup_prof","",kTRUE, 10, 0., 40000., AliReducedVarManager::kVZEROTotalMult, 22,-2000., 20000., AliReducedVarManager::kNTracksTPCoutFromPileup, 20,-1.0,1.0, AliReducedVarManager::kEta, "","","", AliReducedVarManager::kDcaXY);
      man->AddHistogram(classStr.Data(),"DCAz_Eta_MultVZERO_NTracksTPCoutFromPileup_prof","",kTRUE, 10, 0., 40000., AliReducedVarManager::kVZEROTotalMult, 22,-2000., 20000., AliReducedVarManager::kNTracksTPCoutFromPileup, 20,-1.0,1.0, AliReducedVarManager::kEta, "","","", AliReducedVarManager::kDcaZ);
      man->AddHistogram(classStr.Data(), "NsigmaTPC", "Nsigma TPC", kFALSE, 200, -3.0, 3.0, AliReducedVarManager::kTPCnSig);
      man->AddHistogram(classStr.Data(), "NsigmaTPC_P", "Nsigma TPC vs p", kFALSE, 50, -3.0, 3.0, AliReducedVarManager::kTPCnSig,50,1.,10., AliReducedVarManager::kP);

        
      if(classStr.Contains("MCTruth")) {
          man->AddHistogram(classStr.Data(), "PtMC", "p_{T} MC", kFALSE, 150, 0., 15.0, AliReducedVarManager::kPtMC);
          man->AddHistogram(classStr.Data(), "PtRec_PtMC", "p_{T} MC vs p_{T} reconstructed", kFALSE, 150, 0., 15.0, AliReducedVarManager::kPtMC, 150, 0., 15.0, AliReducedVarManager::kPt);
          man->AddHistogram(classStr.Data(), "PhiMC", "#varphi MC", kFALSE, 180, 0., 6.3, AliReducedVarManager::kPhiMC);
          man->AddHistogram(classStr.Data(), "PhiRec_PhiMC", "#varphi MC vs #varphi reconstructed", kFALSE, 180, 0., 6.3, AliReducedVarManager::kPhiMC, 180, 0., 6.3, AliReducedVarManager::kPhi);
          man->AddHistogram(classStr.Data(), "EtaMC", "#eta MC", kFALSE, 100, -1.0, 1.0, AliReducedVarManager::kEtaMC);
          man->AddHistogram(classStr.Data(), "EtaRec_EtaMC", "#eta MC vs #eta reconstructed", kFALSE, 100, -1.0, 1.0, AliReducedVarManager::kEtaMC, 100, -1.0, 1.0, AliReducedVarManager::kEta);          
          man->AddHistogram(classStr.Data(), "PDGcode0", "PDG code of the track", kFALSE, 12000, -6000., 6000.0, AliReducedVarManager::kPdgMC);
          man->AddHistogram(classStr.Data(), "PDGcode1", "PDG code of the track's mother", kFALSE, 12000, -6000., 6000.0, AliReducedVarManager::kPdgMC+1);
          man->AddHistogram(classStr.Data(), "PDGcode2", "PDG code of the track's grand-mother", kFALSE, 12000, -6000., 6000.0, AliReducedVarManager::kPdgMC+2);
          man->AddHistogram(classStr.Data(), "PDGcode3", "PDG code of the track's grand-grand mother", kFALSE, 12000, -6000., 6000.0, AliReducedVarManager::kPdgMC+3);
      }
      continue;
    }  // end if "TrackQA"
    
    if (classStr.Contains("PairQualityFlags")) {
       man->AddHistClass(classStr.Data());
       cout << "  " << classStr.Data() << endl;
       TString pairQualityFlagNames = " ;K^{0}_{S}#rightarrow#pi^{+}#pi^{-};#Lambda#rightarrow p#pi^{-};#bar{#Lambda}#rightarrow #bar{p}#pi^{+};#gamma#rightarrow e^{+}e^{-};";
       man->AddHistogram(classStr.Data(), "PairQualityFlags", "Pair quality flags;;", kFALSE,
                         32, -0.5, 31.5, AliReducedVarManager::kPairQualityFlag, 0, 0.0, 0.0, AliReducedVarManager::kNothing, 0, 0.0, 0.0, AliReducedVarManager::kNothing, pairQualityFlagNames.Data());
       continue;
    }
    
    // Histograms for pairs
    if (classStr.Contains("Pair_")) {
       man->AddHistClass(classStr.Data());
       cout << "  " << classStr.Data() << endl;
       man->AddHistogram(classStr.Data(), "CandidateId", "Candidate id", kFALSE,
                         5, -0.5, 4.5, AliReducedVarManager::kCandidateId, 0, 0.0, 0.0, AliReducedVarManager::kNothing, 0, 0.0, 0.0, AliReducedVarManager::kNothing, candidateNames.Data());
       man->AddHistogram(classStr.Data(), "PairType", "Pair type", kFALSE, 4, -0.5, 3.5, AliReducedVarManager::kPairType);
       man->AddHistogram(classStr.Data(), "PairChi2", "Pair #chi^{2}", kFALSE, 200, 0.0, 50, AliReducedVarManager::kPairChisquare);
       man->AddHistogram(classStr.Data(), "Mass", "Invariant mass", kFALSE, nMassBins, massRange[0], massRange[1], AliReducedVarManager::kMass);
       man->AddHistogram(classStr.Data(), "PseusoproperDecayTime", "Pseudoproper decay time", kFALSE,200, -1.0, 1.0, AliReducedVarManager::kPseudoProperDecayTime);
      /* man->AddHistogram(classStr.Data(), "Mass_V0K0s", "Invariant mass, K^{0}_{s} assumption", kFALSE,
                         nMassBins, massRange[0], massRange[1], AliReducedVarManager::kMassV0);
       man->AddHistogram(classStr.Data(), "Mass_V0Lambda", "Invariant mass, #Lambda^{0} assumption", kFALSE,
                         nMassBins, massRange[0], massRange[1], AliReducedVarManager::kMassV0+1);
       man->AddHistogram(classStr.Data(), "Mass_V0ALambda", "Invariant mass, #bar{#Lambda^{0}} assumption", kFALSE,
      nMassBins, massRange[0], massRange[1], AliReducedVarManager::kMassV0+2);
       man->AddHistogram(classStr.Data(), "Mass_V0Gamma", "Invariant mass, #gamma conversion assumption", kFALSE,
      nMassBins, massRange[0], massRange[1], AliReducedVarManager::kMassV0+3);
       */ man->AddHistogram(classStr.Data(), "Pt", "", kFALSE, 1000, 0.0, 10.0, AliReducedVarManager::kPt);
       man->AddHistogram(classStr.Data(), "Px", "", kFALSE, 1000, 0.0, 10.0, AliReducedVarManager::kPx);
       man->AddHistogram(classStr.Data(), "Py", "", kFALSE, 1000, 0.0, 10.0, AliReducedVarManager::kPy);
       man->AddHistogram(classStr.Data(), "Pz", "", kFALSE, 1000, 0.0, 10.0, AliReducedVarManager::kPz);
       man->AddHistogram(classStr.Data(), "P", "", kFALSE, 1000, 0.0, 10.0, AliReducedVarManager::kP);
       man->AddHistogram(classStr.Data(), "Rapidity", "Rapidity", kFALSE, 240, -1.2, 1.2, AliReducedVarManager::kRap);
       man->AddHistogram(classStr.Data(), "Eta", "Pseudo-rapidity #eta", kFALSE, 240, -2.0, 2.0, AliReducedVarManager::kEta);
       man->AddHistogram(classStr.Data(), "Phi", "Azimuthal distribution", kFALSE, 315, 0.0, 6.3, AliReducedVarManager::kPhi);
       man->AddHistogram(classStr.Data(), "Theta", "#theta distribution", kFALSE, 1000, 0.0, 3.2, AliReducedVarManager::kTheta);
       man->AddHistogram(classStr.Data(), "LxyOrR", "L_{xy}/Decay Radius", kFALSE, 1000, -10.0, 10.0, AliReducedVarManager::kPairLxy);
       man->AddHistogram(classStr.Data(), "OpeningAngle", "Opening angle", kFALSE, 1000, 0.0, 3.2, AliReducedVarManager::kPairOpeningAngle);
       man->AddHistogram(classStr.Data(), "PointingAngle", "Pointing angle", kFALSE, 1000, 0.0, 3.2, AliReducedVarManager::kPairPointingAngle);
       
       continue;
    }   // end if for Pair classes of histograms 
    
    if (classStr.Contains("PairME")) {
      man->AddHistClass(classStr.Data());
      cout << "  " << classStr.Data() << endl;
      man->AddHistogram(classStr.Data(), "PairType", "Pair type", kFALSE, 4, -0.5, 3.5, AliReducedVarManager::kPairType);
      man->AddHistogram(classStr.Data(), "Mass", "Invariant mass", kFALSE,500, 0.0, 5.0, AliReducedVarManager::kMass);
      man->AddHistogram(classStr.Data(), "Pt", "", kFALSE, 100, 0.0, 10.0, AliReducedVarManager::kPt);
      man->AddHistogram(classStr.Data(), "Rapidity", "Rapidity", kFALSE, 300, -1.5, 1.5, AliReducedVarManager::kRap);
      man->AddHistogram(classStr.Data(), "Eta", "Eta", kFALSE, 240, -1.2, 1.2, AliReducedVarManager::kEta);
      man->AddHistogram(classStr.Data(), "Phi", "Azimuthal distribution", kFALSE, 630, 0.0, 6.3, AliReducedVarManager::kPhi);
      continue;
    } 

    if (classStr.Contains("MCTRUTH")) {
      man->AddHistClass(classStr.Data());
      cout << "  " << classStr.Data() << endl;
      man->AddHistogram(classStr.Data(), "Mass", "Invariant mass", kFALSE, nMassBins, massRange[0], massRange[1], AliReducedVarManager::kMassMC);
      man->AddHistogram(classStr.Data(), "PseusoproperDecayTime", "Pseudoproper decay time", kFALSE,200, -1.0, 1.0, AliReducedVarManager::kPseudoProperDecayTimeMC);
      man->AddHistogram(classStr.Data(), "Pt", "", kFALSE, 500, 0.0, 30.0, AliReducedVarManager::kPtMC);
      man->AddHistogram(classStr.Data(), "Px", "", kFALSE, 500, 0.0, 10.0, AliReducedVarManager::kPxMC);
      man->AddHistogram(classStr.Data(), "Py", "", kFALSE, 500, 0.0, 10.0, AliReducedVarManager::kPyMC);
      man->AddHistogram(classStr.Data(), "Pz", "", kFALSE, 500, 0.0, 10.0, AliReducedVarManager::kPzMC);
      man->AddHistogram(classStr.Data(), "P", "", kFALSE, 500, 0.0, 10.0, AliReducedVarManager::kPMC);
      man->AddHistogram(classStr.Data(), "Rapidity", "Rapidity", kFALSE, 240, -1.2, 1.2, AliReducedVarManager::kRapMC);
      man->AddHistogram(classStr.Data(), "Eta", "Pseudo-rapidity #eta", kFALSE, 240, -2.0, 2.0, AliReducedVarManager::kEtaMC);
      man->AddHistogram(classStr.Data(), "Phi", "Azimuthal distribution", kFALSE, 315, 0.0, 6.3, AliReducedVarManager::kPhiMC);
      man->AddHistogram(classStr.Data(), "Theta", "#theta distribution", kFALSE, 1000, 0.0, 3.2, AliReducedVarManager::kThetaMC);
      man->AddHistogram(classStr.Data(), "PDG", "PDG value", kFALSE, 300, -600, 600, AliReducedVarManager::kPdgMC);
      man->AddHistogram(classStr.Data(), "PDGmother", "PDG value jpsi mother", kFALSE, 300, -6000, 6000, AliReducedVarManager::kPdgMC+1);
      continue;
    }   // end if for Pair classes of histograms 
		
		cout << " !" << classStr.Data() << " is unknown and will be ignored" << endl;
  }  // end loop over histogram classes
  cout << endl;
}
