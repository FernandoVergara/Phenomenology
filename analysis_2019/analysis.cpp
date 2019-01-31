#include <iostream>
#include "ROOTFunctions.h"
#include "DelphesFunctions.h"


//Global variables
float mt=173.21;
float mw=80.385;
double m_electron=0.511e-3;
double m_muon=105e-3;
double pi = 4.0*atan(1.0);

using namespace std;

//FUNcTION DEFINITIONS
//PHI RESCALATION
double d_phi(double phi1, double phi2) {
  double dphi = fabs(phi1-phi2);
  if(dphi>pi) dphi = 2*pi-dphi;
  return dphi;
}

// MAIN PROGRAM STARTED
int main(int argc, char *argv[]) {
  // Load shared library
  gSystem->Load("lib/libExRootAnalysis.so");
  gSystem->Load("libPhysics");
  // Root ntuplas
  TChain chain("Delphes");
  chain.Add(argv[1]);
  // Create object of class ExRootTreeReader
  ExRootTreeReader *treeReader_Delphes = new ExRootTreeReader(&chain);
  // Number of entries access
  long int numberOfEntries = treeReader_Delphes->GetEntries();

  // Get pointers to branches used in this analysis
  TClonesArray *branchJet = treeReader_Delphes->UseBranch("Jet");
  TClonesArray *branchElectron = treeReader_Delphes->UseBranch("Electron");
  TClonesArray *branchMuon = treeReader_Delphes->UseBranch("Muon");
  TClonesArray *branchMissingET = treeReader_Delphes->UseBranch("MissingET");
  TClonesArray *branchScalarHT = treeReader_Delphes->UseBranch("ScalarHT");  

  // Book histograms //////////////////////////////////////////////////////////////////////////
  TH1 *hist_RM = new TH1F("RM","RM",20,0.0,1.0);
  TH1 *hist_leading_jet = new TH1F("ljpt","ljpt",20,0.0,1000);
  TH1 *histNumJets = new TH1F("num_jets","JETS",15, 0, 15);
  TH1 *hist_count = new TH1F("cuts","cuts",20,0,20);

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////

  //GLOBAL COUNTER
  int cut_counter[15]={15*0};

  // LOOP OVER ALL EVENTS
  for ( Int_t entry = 0; entry < numberOfEntries; ++entry ){
    cut_counter[0]++;   
    // Load selected branches with data from specified event
    treeReader_Delphes->ReadEntry(entry);

    //Number of jets, muons, electrons and MET in the event
    int NumJets = branchJet->GetEntries();
    histNumJets->Fill(NumJets);
    int NumElectrons = branchElectron->GetEntries();
    int NumMuons = branchMuon->GetEntries();
    int NumMissingET = branchMissingET->GetEntries();

    //POINTERS
    MissingET *metpointer;
    Muon *Muonpointer;
    Electron *Electronpointer;
    Jet *isr_jet;
    Jet *jet_i;
    //LOCAL VARIABLES
    double delta_phi_isr_MET, delta_phi_jet_MET;
    double rm;
    int nbtag=0;
    int nj123=0, npt60=0;

    // LEPTON VETO: PT(l)<10 || |eta(l)|>2.5 
    if(NumElectrons>0 || NumMuons>0) continue; cut_counter[1]++;
    //CUT IN JET NUMBER 
    if(NumJets<4) continue;  		//>>>>>>>>>>>>>>>> N(j)>=4 <<<<<<<<<<
      metpointer = (MissingET*)branchMissingET->At(0);
      isr_jet = (Jet*)branchJet->At(0);
    //CUT IN LEADING JET PT 
      if(isr_jet->PT<100) continue; cut_counter[2]++;   //>>>> PT(ISR)<700 <<<<<<<
    //ISR JET B-TAG VETO
      if(isr_jet->BTag==1) continue; cut_counter[3]++;
    //CUT IN DELTA PHI ISR-JET & MET
      delta_phi_isr_MET = d_phi(isr_jet->Phi,metpointer->Phi);
      if(delta_phi_isr_MET<(pi-0.15)) continue; cut_counter[4]++;//>>> <0.15 <<<<<<
      //loop for jet1,2,3 PT analysis
      for(int i=1; i<4; i++){
        jet_i = (Jet*)branchJet->At(i);
        if(jet_i->PT<60) nj123++;
      }
    //CUT IN JET1,2,3 PT<60
      if(nj123>0) continue; cut_counter[5]++;
      //loop for delta phi>0.2 analysis
      for(int i=0; i<NumJets; i++){
        jet_i = (Jet*)branchJet->At(i);
	if(jet_i->PT<60) break;
        delta_phi_jet_MET = d_phi(jet_i->Phi,metpointer->Phi);
        if(delta_phi_jet_MET<0.2) npt60++;//>>>>>>>>>>>>>>>> <0.55 <<<<<<<<<
      }
    //CUT FOR DELTA PHI>0.2 IF PT>60
      if(npt60>0) continue; cut_counter[6]++;
      //loop to count BTags   
      for(int i=0; i<NumJets; i++){
        jet_i = (Jet*)branchJet->At(i);
        nbtag += jet_i->BTag;           //counts BTags in the event
      }
    //CUT IN BTAGS
      if(nbtag<1) continue; cut_counter[7]++;   //>>>>>>>> N(BTags)<2 <<<<<<<<<<<

      //////////////////////////////////////////////////////////////////////////

      //HISTOGRAMS WITH BASELINE CUTS 
      rm = metpointer->MET/isr_jet->PT;
	hist_RM->Fill(rm); if(rm>1) hist_RM->Fill(0.99);
        hist_leading_jet->Fill(isr_jet->PT);

    //END CUTS
  }//END LOOP OVER ALL EVENTS

for(int i=0;i<13;i++ ){
  for(int j=1;j<cut_counter[i]+1;j++){
    hist_count->Fill(i+0.5);
  }
}

cout<<"######################### CUTS SUMMARY ########  An-Lian Tao Code  ##"<<endl;
cout<<"                         TOTAL NUMBER OF EVENTS  [0] = "<<cut_counter[0]<<endl;
cout<<" ------------------------------------------------------------"<<endl;
cout<<"     BASELINE CUTS"<<endl;
cout<<" ------------------------------------------------------------"<<endl;
cout<<" Lepton veto                        [1] = "<<cut_counter[1]<<endl;
cout<<" leading-jet PT             (>700)  [2] = "<<cut_counter[2]<<endl;
cout<<" ISR-jet BTag veto            (=1)  [3] = "<<cut_counter[3]<<endl;
cout<<" leading-jet Phi & MET (> pi-0.15)  [4] = "<<cut_counter[4]<<endl;
cout<<" jet1,2,3 PT                 (>60)  [5] = "<<cut_counter[5]<<endl;
cout<<" delta phi jets & MET       (>0.2)  [6] = "<<cut_counter[6]<<endl;
cout<<" B-tags                       (>1)  [7] = "<<cut_counter[7]<<endl;
cout<<" ------------------------------------------------------------"<<endl;
cout<<"                              EVENTS LEFT AFTER CUTS = "<<cut_counter[7]<<endl;
cout<<"####################################################################"<<endl;

  // ROOT Program output
  TFile* hfile = new TFile(argv[2], "RECREATE");
    hist_RM->Write();
    histNumJets->Write();
    hist_leading_jet->Write();
    hist_count->Write();

hfile->Close();
}  //end of the program

