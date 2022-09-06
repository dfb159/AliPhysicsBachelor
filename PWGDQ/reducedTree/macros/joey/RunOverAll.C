
/*
 * TODO
 * check gsi data : ~2 Full reconstructed tracks
 * check old monte carlo : 1~2 full reconstructed tracks, no MCTruth
 * check injected mc : ~3 full reconstructed tracks + dozen of full MCTruth
 * gsi data tag14 is set on a hand full of events; 15/1000 with #tracks >= 2 :)
 * 
 * run pairing on example files
 * injected mc pairing: 4100/4500 pairs
 * get total amount of reconstructable pairs, maybe some info not detected?
 * what about multiple per event?
 * 
 * prepare maf for all files
 * => all data/mc_signal/mc_bck in one file
 * train MVA on file
 * 
 * pairing MC
 * True Positives
 * False positives
 * False Negatives (low means high efficiency)
 */

#include <iostream>
#include <array>
#include <string>


#include <TTree.h>
#include <TFile.h>
#include <AliReducedPairInfo.h>
#include <AliReducedEventInfo.h>
#include <AliReducedTrackInfo.h>


using std::cout;
using std::endl;
using std::flush;

void printTrackInfo(AliReducedTrackInfo* track) {
    cout << 
    " - ID = " << track->TrackId() <<
    ", mcFlags = " << track->GetMCFlags() <<
    

    // Track Information
    ", p = " << track->P() <<
    ", pt = " << track->Pt() <<
    ", phi = " << track->Phi() <<
    ", theta = " << track->Theta() <<
    ", eta = " << track->Eta() <<
    ", charge = " << track->Charge() <<
    
    //", massTracking = " << track->MassForTracking() <<
    //", trackLength = " << track->TrackLength() <<
    ", dcaXY = " << track->DCAxy() <<
    ", dcaZ = " << track->DCAz() <<
    //", helixX = " << track->HelixX() <<
    //", helixY = " << track->HelixY() <<
    //", helixR = " << track->HelixR() <<
    endl << "   " <<

    //", itsnSigElec = " << track->ITSnSig(0) <<
    //", itsnSigPion = " << track->ITSnSig(1) <<
    //", itsnSigKaon = " << track->ITSnSig(2) <<
    //", itsnSigProt = " << track->ITSnSig(3) <<
    //", itsSignal = " << track->ITSsignal() <<
    ", itsChi2 = " << track->ITSchi2() <<
    ", itsNcls = " << track->ITSncls() <<
    //", itsNclsShared = " << track->ITSnSharedCls() <<

    //", tpcnSigElec = " << track->TPCnSig(0) <<
    //", tpcnSigPion = " << track->TPCnSig(1) <<
    //", tpcnSigKaon = " << track->TPCnSig(2) <<
    //", tpcnSigProt = " << track->TPCnSig(3) <<
    ", tpcnSig = {" << track->TPCnSig(0) << "," << track->TPCnSig(1) << "," << track->TPCnSig(2) << "," << track->TPCnSig(3) << "}" <<
    ", tpcSignal = " << track->TPCsignal() <<
    //", tpcSignalTuned = " << track->TPCsignalTunedOnData() <<
    ", tpcChi2 = " << track->TPCchi2() <<
    ", tpcSignalN = " << (int) track->TPCsignalN() <<
    ", tpcNcls = " << (int) track->TPCncls() <<
    //", tpcNclsShared = " << (int) track->TPCnclsShared() <<
    ", tpcRows = " << (int) track->TPCCrossedRows() <<

    endl << "   " <<
    /*
    ", tpcdEdxQmaxIROC = " << track->TPCdEdxInfoQmax(0) <<
    ", tpcdEdxQmaxMedOROC = " << track->TPCdEdxInfoQmax(1) <<
    ", tpcdEdxQmaxLongOROC = " << track->TPCdEdxInfoQmax(2) <<
    ", tpcdEdxQmaxAllOROC = " << track->TPCdEdxInfoQmax(3) <<
    ", tpcdEdxQtotIROC = " << track->TPCdEdxInfoQtot(0) <<
    ", tpcdEdxQtotMedOROC = " << track->TPCdEdxInfoQtot(1) <<
    ", tpcdEdxQtotLongOROC = " << track->TPCdEdxInfoQtot(2) <<
    ", tpcdEdxQtotAllOROC = " << track->TPCdEdxInfoQtot(3) <<
    ", tpcActiveLength = " << track->TPCActiveLength() <<
    ", tpcGeomLength = " << track->TPCGeomLength() <<
    //*/
    /*
    ", tofnSigElec = " << track->TOFnSig(0) <<
    ", tofnSigPion = " << track->TOFnSig(1) <<
    ", tofnSigKaon = " << track->TOFnSig(2) <<
    ", tofnSigProt = " << track->TOFnSig(3) <<
    ", tofBeta = " << track->TOFbeta() <<
    ", tofTime = " << track->TOFtime() <<
    ", tofMisProbab = " << track->TOFmismatchProbab() <<
    ", tofChi2 = " << track->TOFchi2() <<
    //*/
    /*
    ", trdNtracklets = " << track->TRDntracklets(0) <<
    ", trdNtrackletsPID = " << track->TRDntracklets(1) <<
    ", trdSig1DElec = " << track->TRDpidLQ1D(0) <<
    ", trdSig1DPion = " << track->TRDpidLQ1D(1) <<
    ", trdSig2DElec = " << track->TRDpidLQ2D(0) <<
    ", trdSig2DPion = " << track->TRDpidLQ2D(1) <<
    //*/
    ", pdg = {" << track->MCPdg(0) << "," << track->MCPdg(1) << "," << track->MCPdg(2) << "," << track->MCPdg(3) << "}" <<
    ", label = {" << track->MCLabel(0) << "," << track->MCLabel(1) << "," << track->MCLabel(2) << "," << track->MCLabel(3) << "}" <<
    endl;
}

void RunOverFile(TString path, TString treeName, int maxRead) {
    cout << " ===== Analysing: " << path << "/" << treeName << " ===== " << endl;
    TFile* fin = TFile::Open(path.Data(), "READ");
    TTree* intree = fin->Get<TTree>(treeName.Data());

    AliReducedEventInfo* event = new AliReducedEventInfo();
    intree->SetBranchAddress("Event",&event);
    int n = intree->GetEntries();
    for(int i = 0; i < n; i++) { // for every event
        if (maxRead > 0 && i > maxRead) break;
        intree->GetEntry(i);

        TClonesArray* fTracks;
        AliReducedTrackInfo* track;
        fTracks = event->GetTracks();
        TIter nextTrack(fTracks);
        //cout << "Event " << i << ", Tracks " << fTracks->GetEntries() << endl;

        int base_mc = 0, base_re = 0, full_mc = 0, full_re = 0;
        for (int j = 0; j < fTracks->GetEntries(); j++) {
            track = (AliReducedTrackInfo*) nextTrack(); // track1 always FullTracks
            switch ((track->IsA() == AliReducedTrackInfo::Class()) + 2* track->IsMCTruth()) {
                case 0: base_re++; break;
                case 1: full_re++; break;
                case 2: base_mc++; break;
                case 3: full_mc++; break;
            }
            //if (i < -1)// || !track->IsMCTruth())
            //if (event->EventTag(14) && fTracks->GetEntries() > 1)
                printTrackInfo(track);
        }
        //*
        if (full_re > 0) {
        //if (fTracks->GetEntries() > 0) {
        //if (event->EventTag(14) && fTracks->GetEntries() > 1) {
            printf("Event: %5i, Tracks1: %4i, BaseRE: %4i, FullRE: %4i, BaseMC: %4i, FullMC: %4i", i, fTracks->GetEntries(), base_re, full_re, base_mc, full_mc);
            cout << ", Tag14: " << event->EventTag(14) << endl;
        }//*/

        /*
        TClonesArray* fTracks2;
        AliReducedBaseTrack* track2;
        fTracks2 = event->GetTracks2();
        TIter nextTrack2(fTracks2);
        //cout << "Event " << i << ", Tracks " << fTracks->GetEntries() << endl;

        base_mc = 0, base_re = 0, full_mc = 0, full_re = 0;
        for (int j = 0; j < fTracks2->GetEntries(); j++) {
            track2 = (AliReducedBaseTrack*) nextTrack2(); // track2 always BaseTracks
            switch ((track2->IsA() == AliReducedTrackInfo::Class()) + 2* track2->IsMCTruth()) {
                case 0: base_re++; break;
                case 1: full_re++; break;
                case 2: base_mc++; break;
                case 3: full_mc++; break;
            }
            //if (i == -1 || !track2->IsMCTruth()) printTrackInfo(track);
        }
        //*
        if (true || full_re > 0) {
        //if (fTracks2->GetEntries() > 0) {
            printf("Event: %5i, Tracks2: %4i, BaseRE: %4i, FullRE: %4i, BaseMC: %4i, FullMC: %4i", i, fTracks2->GetEntries(), base_re, full_re, base_mc, full_mc);
            cout << endl;
        }//*/
        //printf("Event: %5i, Tracks1: %4i, Tracks2: %4i", i, fTracks->GetEntries(), fTracks2->GetEntries()); cout << endl;        
    }
    intree->ResetBranchAddresses();

    fin->Close();
    cout << " ==========================================\n" << endl;
}

void RunOverAll() {
    int N = 10;
    /*
    RunOverFile("/home/joey/Documents/raw/17k/274690/dstTree_data.root", "DstTree", N);
    RunOverFile("/home/joey/Documents/raw/17k/274690/dstTree_mc.root", "DstTree", N);
    RunOverFile("/home/joey/Documents/raw/17k/274690/dstTree_mc_old.root", "DstTree", N);
    //*/
    //RunOverFile("/home/joey/Documents/raw/16k/256941/dstTree_data_gsi.root", "DstTree", N);
    RunOverFile("/home/joey/Documents/raw/16k/256941/dstTree_mc_inject.root", "DstTree", N);
    //RunOverFile("/home/joey/Documents/raw/16k/256941/dstTree_mc_old.root", "DstTree", N);
}
