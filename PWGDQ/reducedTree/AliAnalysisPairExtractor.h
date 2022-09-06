//
// Definition of the DataExtractor for further TMVA analysis on Pairs
//
// Creation date: 2022/05/05
// Author: Jonathan Sigrist, j.sigrist@wwu.de

#ifndef ALIANALYSISPAIREXTRACTOR_H
#define ALIANALYSISPAIREXTRACTOR_H

#include <map>
#include <string>

#include <TTree.h>
#include <TFile.h>
#include <TMap.h>
#include <AliReducedPairInfo.h>
#include <AliReducedEventInfo.h>
#include <AliReducedTrackInfo.h>

class AliAnalysisPairExtractor {

public:
    AliAnalysisPairExtractor(AliReducedPairInfo::CandidateType type=AliReducedPairInfo::CandidateType::kJpsiToEE); // Set up trees
    AliAnalysisPairExtractor(int mother, int leg1, int leg2); // Set up trees
  virtual ~AliAnalysisPairExtractor() {}; // TODO clear trees map
  
public: 
    ULong_t extractDirectory(const TString path, const TString fileName, const TString treeName="DstTree", const TString outName="tree", const TString outDescription="", const Bool_t isMC=kFALSE, const ULong_t N=kMaxULong); // extracts all files in subdirectories up to N pairs
    ULong_t extractFile(const TString path, const TString treeName="DstTree", const TString outName="tree", const TString outDescription="", const Bool_t isMC=kFALSE, const ULong_t N=kMaxULong); // opens, extracts and closes the file automatically
    ULong_t extractTree(TTree* intree, const TString outName="tree", const TString outDescription="", const Bool_t isMC=kFALSE, const ULong_t N=kMaxULong); // extracts mc from tree into signalTree and backgroundTree
    
    void setPDG(int mother, int leg1, int leg2) {pdgMother=mother; pdgLeg1=leg1; pdgLeg2=leg2;}
    void SetUp(TString outpath); // Set outfile
    void Write(); // Write data header to outfile and close it
    
    Bool_t checkTreeIntegrity(TTree* inTree); // TODO Check tree for valid extraction format
    TTree* getOutputTree(const TString treeName, const TString treeDescription); // Returns the correct tree with treeName. Sets a new tree up, if necessary.

private:

    Int_t pdgMother, pdgLeg1, pdgLeg2;

    TFile* outfile;
    TMap trees;
    
    void createBranches(TTree* tree);
    void fillVars(AliReducedEventInfo* event, AliReducedPairInfo* pair, AliReducedTrackInfo* leg1, AliReducedTrackInfo* leg2);
    
    // Event Information
    Int_t runNo; // Run no.
    Float_t vtxX, vtxY, vtxZ; // vertex approximation in cm
    Int_t vtxN; // number of vertex contributors
    Float_t diaX, diaY, diaZ; // diamond size
    Int_t totalTracks; // total reconstructed tracks

    Float_t centV0, centSPD, centTPC, centZEMvsZDC, centV0A, centV0C, centZNA,
        centV0Mnew, centV0MnewPlus05, centV0MnewMinus05, centV0MnewPlus10, centV0MnewMinus10,
        centV0MPlus05, centV0MMinus05, centV0MPlus10, centV0MMinus10; // centrality from different sources // TODO reduce to one centrality or just use one
    
    // Pair Information
    Float_t pairP, pairPt, pairPhi, pairTheta, pairEta; // basic kinematic properties
    Float_t pairMass, pairEnergy, pairRapidity; // properties using JPsi estimation
    Float_t pairDecayRadius; // radius of the secondary vertex for V0s // measured xy-distance to primary vertex // https://arxiv.org/pdf/1205.5880.pdf
    Float_t pairPsProper; // pseudo-proper decay length (pair candidates) // https://arxiv.org/pdf/1205.5880.pdf
    Float_t pairPointingAngle; // angle between the pair momentum vector and the secondary vertex position vector
    Float_t pairChi2; // chi2 for the legs matching
    
    // Track Information Leg1
    Float_t leg1p, leg1pt, leg1phi, leg1theta, leg1eta;
    Int_t leg1charge;
    
    Float_t leg1massTracking, leg1trackLength;
    Float_t leg1dcaXY, leg1dcaZ, leg1helixX, leg1helixY, leg1helixR;
    
    Float_t leg1itsnSigElec, leg1itsnSigPion, leg1itsnSigKaon, leg1itsnSigProt;
    Float_t leg1itsSignal, leg1itsChi2;
    UInt_t leg1itsNcls, leg1itsNclsShared;
    
    Float_t leg1tpcnSigElec, leg1tpcnSigPion, leg1tpcnSigKaon, leg1tpcnSigProt;
    Float_t leg1tpcSignal, leg1tpcSignalTuned, leg1tpcChi2;
    Float_t leg1tpcdEdxQmaxIROC, leg1tpcdEdxQmaxMedOROC, leg1tpcdEdxQmaxLongOROC, leg1tpcdEdxQmaxAllOROC;
    Float_t leg1tpcdEdxQtotIROC, leg1tpcdEdxQtotMedOROC, leg1tpcdEdxQtotLongOROC, leg1tpcdEdxQtotAllOROC;
    Float_t leg1tpcActiveLength, leg1tpcGeomLength;
    UInt_t leg1tpcSignalN, leg1tpcNcls, leg1tpcNclsShared, leg1tpcRows;
    
    Float_t leg1tofnSigElec, leg1tofnSigPion, leg1tofnSigKaon, leg1tofnSigProt;
    Float_t leg1tofBeta, leg1tofTime, leg1tofMisProbab, leg1tofChi2;
    
    UInt_t leg1trdNtracklets, leg1trdNtrackletsPID;
    UInt_t leg1trdSig1DElec, leg1trdSig1DPion;
    UInt_t leg1trdSig2DElec, leg1trdSig2DPion;
    
    // Track Information Leg2
    Float_t leg2p, leg2pt, leg2phi, leg2theta, leg2eta;
    Int_t leg2charge;
    
    Float_t leg2massTracking, leg2trackLength;
    Float_t leg2dcaXY, leg2dcaZ, leg2helixX, leg2helixY, leg2helixR;
    
    Float_t leg2itsnSigElec, leg2itsnSigPion, leg2itsnSigKaon, leg2itsnSigProt;
    Float_t leg2itsSignal, leg2itsChi2;
    UInt_t leg2itsNcls, leg2itsNclsShared;
    
    Float_t leg2tpcnSigElec, leg2tpcnSigPion, leg2tpcnSigKaon, leg2tpcnSigProt;
    Float_t leg2tpcSignal, leg2tpcSignalTuned, leg2tpcChi2;
    Float_t leg2tpcdEdxQmaxIROC, leg2tpcdEdxQmaxMedOROC, leg2tpcdEdxQmaxLongOROC, leg2tpcdEdxQmaxAllOROC;
    Float_t leg2tpcdEdxQtotIROC, leg2tpcdEdxQtotMedOROC, leg2tpcdEdxQtotLongOROC, leg2tpcdEdxQtotAllOROC;
    Float_t leg2tpcActiveLength, leg2tpcGeomLength;
    UInt_t leg2tpcSignalN, leg2tpcNcls, leg2tpcNclsShared, leg2tpcRows;
    
    Float_t leg2tofnSigElec, leg2tofnSigPion, leg2tofnSigKaon, leg2tofnSigProt;
    Float_t leg2tofBeta, leg2tofTime, leg2tofMisProbab, leg2tofChi2;
    
    UInt_t leg2trdNtracklets, leg2trdNtrackletsPID;
    UInt_t leg2trdSig1DElec, leg2trdSig1DPion;
    UInt_t leg2trdSig2DElec, leg2trdSig2DPion;
    
    ClassDef(AliAnalysisPairExtractor, 42); //Analysis Task for extracting information from reduced tree candidates
};

#endif
