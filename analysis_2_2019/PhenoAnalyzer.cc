/*
@file PhenoAnalyzer.cc
@author Andres Florez
@author Carlos Miguel Patino
@date April 2, 2017

Code used to perform phenomenological analysis of Heavy Neutrinos in the tau channel
*/

#include "PhenoAnalyzer.h"

int main(int argc, char *argv[]) {

  //TApplication app("App",&argc, argv);
 // gSystem->Load("libDelphes.so");
  TChain chain("Delphes");
  chain.Add(argv[1]);
  TFile * HistoOutputFile = new TFile(argv[2], "RECREATE");
  int nDir = 16;
  TDirectory *theDirectory[nDir];
  theDirectory[0]  = HistoOutputFile->mkdir("No_cuts");
  // Cuts
  theDirectory[1]  = HistoOutputFile->mkdir("singleLepton_tau_noBjets");
  theDirectory[2]  = HistoOutputFile->mkdir("singleLepton_tau_jets_kinematics");
  theDirectory[3]  = HistoOutputFile->mkdir("singleLepton_tau_MET");
  theDirectory[4]  = HistoOutputFile->mkdir("singleLepton_tau_VBF_deltaEta");
  theDirectory[5]  = HistoOutputFile->mkdir("singleLepton_tau_VBF_diJetMass");
  // Single tau
  theDirectory[6]  = HistoOutputFile->mkdir("singleLepton_tau_kinematics");
  // Single muon
  theDirectory[7]  = HistoOutputFile->mkdir("singleLepton_muon_VBF_diJetMass");
  /////////////////////////////////// Nestor ///////////////////////////////////////
  // Single electron
  theDirectory[8]  = HistoOutputFile->mkdir("singleLepton_electron_VBF_diJetMass");
  //////////////////////////////////////////////////////////////////////////////////
  // muTau pair 
  theDirectory[9]  = HistoOutputFile->mkdir("DiLepton_muTau_VBF_diJetMass");
  //(MM: Added folders for muMu pairs)
  theDirectory[10]  = HistoOutputFile->mkdir("DiLepton_muMu_VBF_diJetMass");
  //(MM: Added folders for tauTau pairs)
  theDirectory[11]  = HistoOutputFile->mkdir("DiLepton_tauTau_VBF_diJetMass");   
  /////////////////////////////////// Nestor ///////////////////////////////////////
  // elecTau pair
  theDirectory[12]  = HistoOutputFile->mkdir("DiLepton_elecTau_VBF_diJetMass");   
  // elecElec pair
  theDirectory[13]  = HistoOutputFile->mkdir("DiLepton_elecElec_VBF_diJetMass");   
  // muMuMu trio
  theDirectory[14]  = HistoOutputFile->mkdir("TriLepton_muMuMu_VBF_diJetMass");   
  // elElEl trio
  theDirectory[15]  = HistoOutputFile->mkdir("TriLepton_elElEl_VBF_diJetMass");   
  //////////////////////////////////////////////////////////////////////////////////
  PhenoAnalysis BSM_analysis(chain, HistoOutputFile, theDirectory, nDir);
  
}

using namespace std;
PhenoAnalysis::PhenoAnalysis(TChain& chain, TFile* theFile, TDirectory *cdDir[], int nDir)
{
  ifstream inFile;
  inFile.open ("config.in", ios::in);
  
  if (!inFile)
    {
      cerr << "ERROR: Can't open input file: " << endl;
      exit (1);
    }
  
  string inputType = "";
  
  //This set of lines are used to open and read the "config.in" file.
  ///////////////////////////////////////////////////////////////////////
  TEnv *params = new TEnv ("config_file");
  params->ReadFile ("config.in", kEnvChange);
  
  double b_jet_pt_min       = params->GetValue ("b_jet_pt_min",      30.0);
  double DR_jet_lep_max     = params->GetValue ("DR_jet_lep_max",     0.3);
  double jet_min_pt         = params->GetValue ("jet_min_pt",          30);
  double jet_max_eta        = params->GetValue ("jet_max_eta",          5);
  double VBF_jetPt_min      = params->GetValue ("VBF_jetPt_min",      30.0);
  double tau_pt_cut         = params->GetValue ("tau_pt_cut",          20.);
  double tau_pt_cut_max     = params->GetValue ("tau_pt_cut_max",      20.);
  double gen_tau_pt_min     = params->GetValue ("gen_tau_pt_cut",       8.);
  double gen_tau2_pt_min    = params->GetValue ("gen_tau2_pt_cut",      8.);
  double tau_eta_cut        = params->GetValue ("tau_eta_cut",         2.3);
  double deltaEta_diJet_cut = params->GetValue ("deltaEta_diJet_cut",  3.8);
  double diJetmass_cut      = params->GetValue ("diJetmass_cut",     500.0);
  double tauMass_cut        = params->GetValue ("tauMass_cut",       100.0);
  double MET_cut            = params->GetValue ("MET_cut",           250.0);
  double transmass_cut      = params->GetValue ("transmass_cut",      50.0);
  double muon_pt_cut        = params->GetValue ("muon_pt_cut",        100.);
  double muon_eta_cut       = params->GetValue ("muon_eta_cut",         8.);
  double gen_muon2_pt_min   = params->GetValue ("gen_tau2_pt_cut",      8.);  
  double elec_pt_cut        = params->GetValue ("elec_pt_cut",        100.);
  double elec_eta_cut       = params->GetValue ("elec_eta_cut",         8.);
  int    isSignal           = params->GetValue ("isSignal",              0);
  double muTau_mass_input   = params->GetValue ("muTau_mass_input",    10.);
  double muMu_mass_input    = params->GetValue ("muMu_mass_input",     10.);
  double tauTau_mass_input  = params->GetValue ("tauTau_mass_input",   10.);
  double elecTau_mass_input = params->GetValue ("elecTau_mass_input",  10.);
  double elecElec_mass_input= params->GetValue ("elecElec_mass_input", 10.);
  
  crateHistoMasps(nDir);

  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  Long64_t numberOfEntries = treeReader->GetEntries();
  
  TClonesArray *branchJet = treeReader->UseBranch("Jet");
  TClonesArray *branchElectron = treeReader->UseBranch("Electron");
  TClonesArray *branchMuon = treeReader->UseBranch("Muon");
  TClonesArray *branchMissingET = treeReader->UseBranch("MissingET");
  TClonesArray *branchTrack = treeReader->UseBranch("Track");
  TClonesArray *branchGenParticle = treeReader->UseBranch("Particle");
  
  MissingET *METpointer;
  std::map<unsigned int, TLorentzVector> muTau_TLV;//(MM: is this an int or a float??)
  std::map<unsigned int, TLorentzVector> pairs_muTau_TLV;
  
  //(MM: Added for tau-tau pairs)
  std::map<unsigned int, TLorentzVector> tauTau_TLV;
  std::map<unsigned int, TLorentzVector> pairs_tauTau_TLV;
  
  //(MM: Added for mu-mu pairs)
  std::map<unsigned int, TLorentzVector> muMu_TLV;
  std::map<unsigned int, TLorentzVector> pairs_muMu_TLV;
 
  /////////////////////////////////// Nestor ///////////////////////////////////////
  // Added for elec-tau pairs
  std::map<unsigned int, TLorentzVector> elecTau_TLV;
  std::map<unsigned int, TLorentzVector> pairs_elecTau_TLV;

  // Added for elec-elec pairs
  std::map<unsigned int, TLorentzVector> elecElec_TLV;
  std::map<unsigned int, TLorentzVector> pairs_elecElec_TLV;

  // Added for mu-mu-mu trio
  std::map<unsigned int, TLorentzVector> muMuMu_TLV;
  std::map<unsigned int, TLorentzVector> trio_muMuMu_TLV;

  // Added for elec-elec-elec trio
  std::map<unsigned int, TLorentzVector> elElEl_TLV;
  std::map<unsigned int, TLorentzVector> trio_elElEl_TLV;
  //////////////////////////////////////////////////////////////////////////////////
 
  for(Int_t entry = 0; entry < numberOfEntries; ++entry){
    
    treeReader->ReadEntry(entry);
    int pass_cuts[nDir];
    TLorentzVector jetLeadingVec(0., 0., 0., 0.); //tc stands for tau channel
    TLorentzVector jetSleadingVec(0., 0., 0., 0.);
    
    TLorentzVector Jet_1(0., 0., 0., 0.);
    TLorentzVector Jet_2(0., 0., 0., 0.);
    
    TLorentzVector jet_i(0., 0., 0., 0.);
    TLorentzVector elec_i(0., 0., 0., 0.);
    TLorentzVector muon_i(0., 0., 0., 0.);
    TLorentzVector tmp_tlv (0., 0., 0., 0.); //Vector used to swap vectors
    
    //Tau Vectors (MM: changed from 3 to 4)
    TLorentzVector Tau1HadTLV (0., 0., 0., 0.);
    TLorentzVector Tau2HadTLV (0., 0., 0., 0.);
    TLorentzVector Tau3HadTLV (0., 0., 0., 0.);
    TLorentzVector Tau4HadTLV (0., 0., 0., 0.);
    
    //Muon Vectors (MM: changed from 3 to 4)
    TLorentzVector Muon1HadTLV (0., 0., 0., 0.);
    TLorentzVector Muon2HadTLV (0., 0., 0., 0.);
    TLorentzVector Muon3HadTLV (0., 0., 0., 0.);
    TLorentzVector Muon4HadTLV (0., 0., 0., 0.);
    
  /////////////////////////////////// Nestor ///////////////////////////////////////
    // Electron Vectors
    TLorentzVector Elec1HadTLV (0., 0., 0., 0.);
    TLorentzVector Elec2HadTLV (0., 0., 0., 0.);
    TLorentzVector Elec3HadTLV (0., 0., 0., 0.);
    TLorentzVector Elec4HadTLV (0., 0., 0., 0.);
  //////////////////////////////////////////////////////////////////////////////////
 
    TLorentzVector TauHadTLV_gen (0., 0., 0., 0.);
    TLorentzVector Tau1HadTLV_gen (0., 0., 0., 0.);
    TLorentzVector Tau2HadTLV_gen (0., 0., 0., 0.);
    TLorentzVector Tau3HadTLV_gen (0., 0., 0., 0.);
    
    
    double gen_tau3_charge = 0.;
    
    double tau1_charge = 0.;
    double tau2_charge = 0.;
    double tau3_charge = 0.;
  /////////////////////////////////// Nestor ///////////////////////////////////////
    double tau4_charge = 0.;
  //////////////////////////////////////////////////////////////////////////////////
    
    double muon1_charge = 0.;
    double muon2_charge = 0.;
    double muon3_charge = 0.;
  /////////////////////////////////// Nestor ///////////////////////////////////////
    double muon4_charge = 0.;

    double elec1_charge = 0.;
    double elec2_charge = 0.;
    double elec3_charge = 0.;
    double elec4_charge = 0.;
  //////////////////////////////////////////////////////////////////////////////////
    
    Track *track_tau1;
    vector<TLorentzVector> jetsList;
    
    bool fill_tau1 = false;
    bool fill_tau2 = false;
    bool fill_tau3 = false;
    bool fill_tau4 = false;//(MM:added)
    
    METpointer = (MissingET*) branchMissingET->At(0);
    double MET = METpointer->MET;
    double Met_phi = METpointer->Phi;
    double tau_transmass = 0.;
    double muon_transmass = 0.;
    double elec_transmass = 0.;
    int njets_counter = 0;
    int ntau_counter = 0;
    int nmuon_counter = 0;
  /////////////////////////////////// Nestor ///////////////////////////////////////
    int nelec_counter = 0;
  //////////////////////////////////////////////////////////////////////////////////
    double DiJetMass_final = 500.;
    int nBJets = 0;
    
    //////////////////Tau Channel///////////
    //Search for taus and bjets
    for (int j = 0; j < branchJet->GetEntriesFast(); j++){

      Jet *jet = (Jet*) branchJet->At(j);
      
      if((jet->BTag == 1) && (jet->PT > b_jet_pt_min)){nBJets++;}
      
      //Tau search
      if((jet->TauTag == 1) && (jet->PT > 20.0)){
	
        ntau_counter++;
        double tau_energy = calculateE(jet->Eta, jet->PT, jet->Mass);
	
        if((fill_tau1 == false) && (fill_tau2 == false) && (fill_tau3 == false) && (fill_tau4==false)){
	  
          Tau1HadTLV.SetPtEtaPhiE(jet->PT, jet->Eta, jet->Phi, tau_energy);
          tau1_charge = jet->Charge;
          fill_tau1 = true;
          continue;
        }
	
        if((fill_tau1 == true) && (fill_tau2 == false) && (fill_tau3 == false) && (fill_tau4==false)){
	    Tau2HadTLV.SetPtEtaPhiE(jet->PT, jet->Eta, jet->Phi, tau_energy);
	    tau2_charge = jet->Charge;
	    fill_tau2 = true;
	    continue;
	  }
	
        if((fill_tau1 == true) && (fill_tau2 == true) && (fill_tau3 == false) && (fill_tau4==false)){
	    Tau3HadTLV.SetPtEtaPhiE(jet->PT, jet->Eta, jet->Phi, tau_energy);
	    tau3_charge = jet->Charge;
	    fill_tau3 = true;
	    continue;
	  }
        //(MM: added extra condition for fourth tau)
        if((fill_tau1 == true) && (fill_tau2 == true) && (fill_tau3 == true) && (fill_tau4==false)){
	    Tau4HadTLV.SetPtEtaPhiE(jet->PT, jet->Eta, jet->Phi, tau_energy);
	    tau4_charge = jet->Charge;
	    fill_tau4 = true;
	    continue;
	  }
      }
    }
    
    bool tau1_muon_overlap = false;
    bool tau2_muon_overlap = false;
    bool tau3_muon_overlap = false;
    bool tau4_muon_overlap = false;//(MM:added)
    
    bool tau1_elec_overlap = false;
    bool tau2_elec_overlap = false;
    bool tau3_elec_overlap = false;
    bool tau4_elec_overlap = false;//(MM:added)
    
    
    bool fill_muon1 = false;
    bool fill_muon2 = false;
    bool fill_muon3 = false;
    bool fill_muon4 = false;//(MM:added)
    
    //Check if taus overlap with muons
    
    for (int muo = 0; muo < branchMuon->GetEntriesFast(); muo++){
      
      if(tau1_muon_overlap && tau2_muon_overlap && tau3_muon_overlap && tau4_muon_overlap){break;}
      
      Muon *muon = (Muon*) branchMuon->At(muo);
      double muon_energy = calculateE(muon->Eta, muon->PT, 0.1056583715);
      muon_i.SetPtEtaPhiE(muon->PT, muon->Eta, muon->Phi, muon_energy);
      
      double DR_tau1_muon = Tau1HadTLV.DeltaR(muon_i);
      double DR_tau2_muon = Tau2HadTLV.DeltaR(muon_i);
      double DR_tau3_muon = Tau3HadTLV.DeltaR(muon_i);
      double DR_tau4_muon = Tau4HadTLV.DeltaR(muon_i);//(MM: added 4 tau-muon overlap)
      
      //(MM: Added DR 4 tau-muon overlap)
      if (DR_tau4_muon < DR_jet_lep_max){
        tau4_muon_overlap = true;
        Tau4HadTLV.SetPtEtaPhiE(0.,0.,0.,0.);
        ntau_counter--;
      }
      
      if (DR_tau3_muon < DR_jet_lep_max){
        tau3_muon_overlap = true;
        Tau3HadTLV.SetPtEtaPhiE(0.,0.,0.,0.);
        ntau_counter--;
      }
      
      if(DR_tau2_muon < DR_jet_lep_max){
        tau2_muon_overlap = true;
        Tau2HadTLV.SetPtEtaPhiE(0.,0.,0.,0.);
        ntau_counter--;
      }
      
      if(DR_tau1_muon < DR_jet_lep_max){
        tau1_muon_overlap = true;
        Tau1HadTLV.SetPtEtaPhiE(0.,0.,0.,0.);
        ntau_counter--;
      }
      
      if((muon->PT > 8.0) && (abs(muon->Eta) < 2.5)){
	if((fill_muon1 == false) && (fill_muon2 == false) && (fill_muon3 == false) && (fill_muon4 == false)){
	  Muon1HadTLV = muon_i;
	  fill_muon1 = true;
	  muon1_charge = muon->Charge;
	  nmuon_counter++;
	  continue;
	}
       if((fill_muon1 == true) && (fill_muon2 == false) && (fill_muon3 == false) && (fill_muon4 == false)){
         if (muon_i.DeltaR(Muon1HadTLV) > 0.3){
           Muon2HadTLV = muon_i;
           fill_muon2 = true;
           muon2_charge = muon->Charge;
           nmuon_counter++;
           continue;
         }
       }
       if((fill_muon1 == true) && (fill_muon2 == true) && (fill_muon3 == false) && (fill_muon4 == false)){
         if ((muon_i.DeltaR(Muon2HadTLV) > 0.3) && (muon_i.DeltaR(Muon1HadTLV) > 0.3)){
           Muon3HadTLV = muon_i;
           fill_muon3 = true;
           muon3_charge = muon->Charge;
           nmuon_counter++;
           continue;
         }
       }
       
       //(MM: Added conditional for 4 muons)
       if((fill_muon1 == true) && (fill_muon2 == true) && (fill_muon3 == true) && (fill_muon4 == false)){
         if ((muon_i.DeltaR(Muon2HadTLV) > 0.3) && (muon_i.DeltaR(Muon1HadTLV) > 0.3)){
           Muon4HadTLV = muon_i;
           fill_muon4 = true;
           muon3_charge = muon->Charge;
           nmuon_counter++;
           continue;
         }
       }
      }
      
    } // for (int muo = 0; muo < branchMuon->GetEntriesFast(); muo++)
    

  /////////////////////////////////// Nestor ///////////////////////////////////////
    bool fill_elec1 = false;
    bool fill_elec2 = false;
    bool fill_elec3 = false;
    bool fill_elec4 = false;
  //////////////////////////////////////////////////////////////////////////////////

    //Check if taus overlap with electrons
    
    for (int el = 0; el < branchElectron->GetEntriesFast(); el++){
      
      if(tau1_elec_overlap && tau2_elec_overlap && tau3_elec_overlap && tau4_elec_overlap){break;}
      
      Electron *elec = (Electron*) branchElectron->At(el);
      double elec_energy = calculateE(elec->Eta, elec->PT, 0.000510998902);
      elec_i.SetPtEtaPhiE(elec->PT, elec->Eta, elec->Phi, elec_energy);
      
      double DR_tau1_elec = Tau1HadTLV.DeltaR(elec_i);
      double DR_tau2_elec = Tau2HadTLV.DeltaR(elec_i);
      double DR_tau3_elec = Tau3HadTLV.DeltaR(elec_i);
      double DR_tau4_elec = Tau4HadTLV.DeltaR(elec_i);//(MM: Added for fourth tau electron overlap)
      
      //(MM: conditional for tau 4 electron overlap)
      if (DR_tau4_elec < DR_jet_lep_max){
        tau4_elec_overlap = true;
        Tau4HadTLV.SetPtEtaPhiE(0.,0.,0.,0.);
        ntau_counter--;
      }
      
      if (DR_tau3_elec < DR_jet_lep_max){
        tau3_elec_overlap = true;
        Tau3HadTLV.SetPtEtaPhiE(0.,0.,0.,0.);
        ntau_counter--;
      }
      
      if(DR_tau2_elec < DR_jet_lep_max){
        tau2_elec_overlap = true;
        Tau2HadTLV.SetPtEtaPhiE(0.,0.,0.,0.);
        ntau_counter--;
      }
      
      if(DR_tau1_elec < DR_jet_lep_max){
        tau1_elec_overlap = true;
        Tau1HadTLV.SetPtEtaPhiE(0.,0.,0.,0.);
        ntau_counter--;
      }
      
  /////////////////////////////////// Nestor ///////////////////////////////////////
      if((elec->PT > 8.0) && (abs(elec->Eta) < 2.5)){
	if((fill_elec1 == false) && (fill_elec2 == false) && (fill_elec3 == false) && (fill_elec4 == false)){
	  Elec1HadTLV = elec_i;
	  fill_elec1 = true;
	  elec1_charge = elec->Charge;
	  nelec_counter++;
	  continue;
	}
       if((fill_elec1 == true) && (fill_elec2 == false) && (fill_elec3 == false) && (fill_elec4 == false)){
         if (elec_i.DeltaR(Elec1HadTLV) > 0.3){
           Elec2HadTLV = elec_i;
           fill_elec2 = true;
           elec2_charge = elec->Charge;
           nelec_counter++;
           continue;
         }
       }
       if((fill_elec1 == true) && (fill_elec2 == true) && (fill_elec3 == false) && (fill_elec4 == false)){
         if ((elec_i.DeltaR(Elec2HadTLV) > 0.3) && (elec_i.DeltaR(Elec1HadTLV) > 0.3)){
           Elec3HadTLV = elec_i;
           fill_elec3 = true;
           elec3_charge = elec->Charge;
           nelec_counter++;
           continue;
         }
       }
       
       if((fill_elec1 == true) && (fill_elec2 == true) && (fill_elec3 == true) && (fill_elec4 == false)){
         if ((elec_i.DeltaR(Elec3HadTLV) > 0.3) && (elec_i.DeltaR(Elec2HadTLV) > 0.3) && (elec_i.DeltaR(Elec1HadTLV) > 0.3)){
           Elec4HadTLV = elec_i;
           fill_elec4 = true;
           elec4_charge = elec->Charge;
           nelec_counter++;
           continue;
         }
       }
      }
  //////////////////////////////////////////////////////////////////////////////////
      
    }
    
    
    //Check if jets overlap with taus
    int n_jets = 0;
    
    for (int l = 0; l < branchJet->GetEntriesFast(); l++){
      
      Jet *jet = (Jet*) branchJet->At(l);
      
      if((jet->PT > jet_min_pt) && (jet->TauTag == 0) && (jet->BTag == 0)){
	
        double jet_i_energy = calculateE(jet->Eta, jet->PT, jet->Mass);
        jet_i.SetPtEtaPhiE(jet->PT, jet->Eta, jet->Phi, jet_i_energy);
	
        double DR_tau1_jet = Tau1HadTLV.DeltaR(jet_i);
        double DR_tau2_jet = Tau2HadTLV.DeltaR(jet_i);
        double DR_tau3_jet = Tau3HadTLV.DeltaR(jet_i);
        double DR_tau4_jet = Tau4HadTLV.DeltaR(jet_i);//(MM: Added tau4 jet overlap)
	
        double DR_muon1_jet = Muon1HadTLV.DeltaR(jet_i);
        double DR_muon2_jet = Muon2HadTLV.DeltaR(jet_i);
        double DR_muon3_jet = Muon3HadTLV.DeltaR(jet_i);
        double DR_muon4_jet = Muon4HadTLV.DeltaR(jet_i);//(MM: Added muon4 jet overlap)
	
  /////////////////////////////////// Nestor ///////////////////////////////////////
        double DR_elec1_jet = Elec1HadTLV.DeltaR(jet_i);
        double DR_elec2_jet = Elec2HadTLV.DeltaR(jet_i);
        double DR_elec3_jet = Elec3HadTLV.DeltaR(jet_i);
        double DR_elec4_jet = Elec4HadTLV.DeltaR(jet_i);
  //////////////////////////////////////////////////////////////////////////////////
	
        if((DR_tau1_jet < DR_jet_lep_max)){
          Tau1HadTLV.SetPtEtaPhiE(0.,0.,0.,0.);
          ntau_counter--;
        }
	
        if(DR_tau2_jet < DR_jet_lep_max){
          Tau2HadTLV.SetPtEtaPhiE(0.,0.,0.,0.);
          ntau_counter--;
        }
	
        if(DR_tau3_jet < DR_jet_lep_max){
          Tau3HadTLV.SetPtEtaPhiE(0.,0.,0.,0.);
          ntau_counter--;
        }
        //(MM: tau4 jet condition overlap)
        if(DR_tau4_jet < DR_jet_lep_max){
          Tau4HadTLV.SetPtEtaPhiE(0.,0.,0.,0.);
          ntau_counter--;
        }
	
        // remove overlaps of jets with muons
        if((DR_muon1_jet < DR_jet_lep_max)){
          Muon1HadTLV.SetPtEtaPhiE(0.,0.,0.,0.);
          nmuon_counter--;
        }
	
        if(DR_muon2_jet < DR_jet_lep_max){
          Muon2HadTLV.SetPtEtaPhiE(0.,0.,0.,0.);
          nmuon_counter--;
        }
	
        if(DR_muon3_jet < DR_jet_lep_max){
          Muon3HadTLV.SetPtEtaPhiE(0.,0.,0.,0.);
          nmuon_counter--;
        }
        //(MM: muon4 jet condition overlap)
        if(DR_muon4_jet < DR_jet_lep_max){
          Muon4HadTLV.SetPtEtaPhiE(0.,0.,0.,0.);
          nmuon_counter--;
        }
	
  /////////////////////////////////// Nestor ///////////////////////////////////////
        // remove overlaps of jets with electrons
        if((DR_elec1_jet < DR_jet_lep_max)){
          Elec1HadTLV.SetPtEtaPhiE(0.,0.,0.,0.);
          nelec_counter--;
        }
	
        if(DR_elec2_jet < DR_jet_lep_max){
          Elec2HadTLV.SetPtEtaPhiE(0.,0.,0.,0.);
          nelec_counter--;
        }
	
        if(DR_elec3_jet < DR_jet_lep_max){
          Elec3HadTLV.SetPtEtaPhiE(0.,0.,0.,0.);
          nelec_counter--;
        }
        
        if(DR_elec4_jet < DR_jet_lep_max){
          Elec4HadTLV.SetPtEtaPhiE(0.,0.,0.,0.);
          nelec_counter--;
        }
  //////////////////////////////////////////////////////////////////////////////////

        if((n_jets <= 6) && (jet->TauTag == 0) && (jet->BTag == 0) && (abs(jet->Eta) < 5.0)){
          jetsList.push_back(jet_i);
          n_jets++;
        }
      }
    }
    
    
    //Order taus by pt
    for(int pt_order = 0; pt_order < 4; pt_order++){ 
      //(MM: comparison with tau4)
      if(Tau3HadTLV.Pt() < Tau4HadTLV.Pt()){
        tmp_tlv = Tau3HadTLV;
        Tau3HadTLV = Tau4HadTLV;
        tau3_charge = tau4_charge;
        Tau4HadTLV = tmp_tlv;
      }     
      
      if(Tau2HadTLV.Pt() < Tau3HadTLV.Pt()){
        tmp_tlv = Tau2HadTLV;
        Tau2HadTLV = Tau3HadTLV;
        tau2_charge = tau3_charge;
        Tau3HadTLV = tmp_tlv;
      }
      if(Tau1HadTLV.Pt() < Tau2HadTLV.Pt()){
        tmp_tlv = Tau1HadTLV;
        Tau1HadTLV = Tau2HadTLV;
        tau1_charge = tau2_charge;
        Tau2HadTLV = tmp_tlv;
      }
    }
    
    //Order muons by pt
    for(int pt_order = 0; pt_order < 4; pt_order++){
      //(MM: comparison with muon4)
      if(Muon3HadTLV.Pt() < Muon4HadTLV.Pt()){
        tmp_tlv = Muon4HadTLV;
        Muon3HadTLV = Muon4HadTLV;
        muon3_charge = muon4_charge;
        Muon4HadTLV = tmp_tlv;
      }
      if(Muon2HadTLV.Pt() < Muon3HadTLV.Pt()){
        tmp_tlv = Muon2HadTLV;
        Muon2HadTLV = Muon3HadTLV;
        muon2_charge = muon3_charge;
        Muon3HadTLV = tmp_tlv;
      }
      if(Muon1HadTLV.Pt() < Muon2HadTLV.Pt()){
        tmp_tlv = Muon1HadTLV;
        Muon1HadTLV = Muon2HadTLV;
        muon1_charge = muon2_charge;
        Muon2HadTLV = tmp_tlv;
      }
    }
    
  /////////////////////////////////// Nestor ///////////////////////////////////////
    //Order electrons by pt
    for(int pt_order = 0; pt_order < 4; pt_order++){
      if(Elec3HadTLV.Pt() < Elec4HadTLV.Pt()){
        tmp_tlv = Elec4HadTLV;
        Elec3HadTLV = Elec4HadTLV;
        elec3_charge = elec4_charge;
        Elec4HadTLV = tmp_tlv;
      }
      if(Elec2HadTLV.Pt() < Elec3HadTLV.Pt()){
        tmp_tlv = Elec2HadTLV;
        Elec2HadTLV = Elec3HadTLV;
        elec2_charge = elec3_charge;
        Elec3HadTLV = tmp_tlv;
      }
      if(Elec1HadTLV.Pt() < Elec2HadTLV.Pt()){
        tmp_tlv = Elec1HadTLV;
        Elec1HadTLV = Elec2HadTLV;
        elec1_charge = elec2_charge;
        Elec2HadTLV = tmp_tlv;
      }
    }
  //////////////////////////////////////////////////////////////////////////////////
    
    // GEN INFO
    bool fill_tau1_gen = false;
    bool fill_tau2_gen = false;
    bool fill_tau3_gen = false;
    
    for(int tau_i = 0; tau_i < branchGenParticle->GetEntriesFast(); tau_i++){
      // Look for the particle PDGID
      GenParticle *tau_gen = (GenParticle*) branchGenParticle->At(tau_i);
      if (abs(tau_gen->PID) == 15){
	double tau_energy = calculateE(tau_gen->Eta, tau_gen->PT, tau_gen->Mass);
	TauHadTLV_gen.SetPtEtaPhiE(tau_gen->PT, tau_gen->Eta, tau_gen->Phi, tau_energy);
	if(fill_tau1_gen==false){
	  Tau1HadTLV_gen = TauHadTLV_gen; 
	  fill_tau1_gen=true; 
	  continue;
	}
	if((fill_tau1_gen==true) && (fill_tau2_gen==false)){ 
	  double delta_R_gen1 = Tau1HadTLV_gen.DeltaR(TauHadTLV_gen);
	  if (delta_R_gen1 > 0.3) {
	    Tau2HadTLV_gen = TauHadTLV_gen;
	    fill_tau2_gen=true; 
	  }
	  continue;
	}
	if((fill_tau2_gen==true) && (fill_tau3_gen==false)){
	  double delta_R_gen2 = Tau2HadTLV_gen.DeltaR(TauHadTLV_gen);
	  if (delta_R_gen2 > 0.3) {
	    Tau3HadTLV_gen = TauHadTLV_gen;
	    fill_tau3_gen=true;
	  }
	  continue;
	}
      }
    }
    
    double tau1_track_DR_min = 999.;
    double tau1_track_DR;
    //Search tau track
    for(int i = 0; i < branchTrack->GetEntriesFast(); i++){
      Track *track = (Track*) branchTrack->At(i);
      tau1_track_DR = calculate_deltaR(Tau1HadTLV, track);
      if(tau1_track_DR < tau1_track_DR_min){
        tau1_track_DR_min = tau1_track_DR;
        track_tau1 = track;
      }
    }
    
    int dijet_index1 = 0;
    int dijet_index2 = 0;
    
    //Search DiJetMass
    for(int k = 0; k < jetsList.size(); k++){
      
      if (jetsList.size() < 2){break;}
      Jet_1 = jetsList[k];
      
      if ((Jet_1.Pt() < VBF_jetPt_min) || (abs(Jet_1.Eta()) > 5.0)){continue;}
      
      for (int sj = k + 1; sj < jetsList.size(); sj++){
	
        if (sj != k){
	  
          Jet_2 = jetsList[sj];
          if ((Jet_2.Pt() < VBF_jetPt_min) || (abs(Jet_2.Eta()) > 5.0)){continue;}
          double DiJetMass = (Jet_1+Jet_2).M();
          if (DiJetMass > DiJetMass_final){
            DiJetMass_final = DiJetMass;
            jetLeadingVec = Jet_1;
            jetSleadingVec = Jet_2;
            dijet_index1 = k;
            dijet_index2 = sj;
          }
        }
      }
    }
    
    if(jetLeadingVec.Pt() < jetSleadingVec.Pt()){
      tmp_tlv = jetLeadingVec;
      jetLeadingVec = jetSleadingVec;
      jetSleadingVec = tmp_tlv;
    }
    
    if(jetsList.size() > 2){
      
      //Remove the diJet pair from the jets list
      tmp_tlv = jetsList[jetsList.size() - 1];
      jetsList[jetsList.size() - 1] = jetsList[dijet_index1];
      jetsList[dijet_index1] = tmp_tlv;
      
      tmp_tlv = jetsList[jetsList.size() - 2];
      jetsList[jetsList.size() -2] = jetsList[dijet_index2];
      jetsList[dijet_index2] = tmp_tlv;
      
      jetsList.pop_back();
      jetsList.pop_back();
      
    }
    
    //Check for jets pt condition
    int jet_pt_condition = 0;
    for (int i = 0; i < jetsList.size(); i++) {
      if((jetsList[i].Pt() > jet_min_pt) && (abs(jetsList[i].Eta()) < 5.0)){
        jet_pt_condition++;
      }
    }
    
    double delta_eta_diJet = abs(jetLeadingVec.Eta()-jetSleadingVec.Eta());
    double tauMass = (Tau1HadTLV + Tau2HadTLV).M();
    tau_transmass = TMath::Sqrt(TMath::Abs(2*Tau1HadTLV.Pt()*MET*(1-TMath::Cos(normalizedDphi(Tau1HadTLV.Phi() - Met_phi)))));
    muon_transmass = TMath::Sqrt(TMath::Abs(2*Muon1HadTLV.Pt()*MET*(1-TMath::Cos(normalizedDphi(Muon1HadTLV.Phi() - Met_phi)))));
    elec_transmass = TMath::Sqrt(TMath::Abs(2*Elec1HadTLV.Pt()*MET*(1-TMath::Cos(normalizedDphi(Elec1HadTLV.Phi() - Met_phi)))));
    double ht = 0.;
    double st = 0.;
    
      
    ht += jetLeadingVec.Pt() + jetSleadingVec.Pt();
      
    for (int i = 0; i < jetsList.size(); i++) {
	
      ht += jetsList[i].Pt();
	
    }
    st += ht + Tau1HadTLV.Pt() + Tau2HadTLV.Pt() + Tau3HadTLV.Pt() + jetLeadingVec.Pt() + jetSleadingVec.Pt() + MET;
    
    muTau_TLV[0]=Tau1HadTLV;
    muTau_TLV[1]=Tau2HadTLV;
    muTau_TLV[2]=Tau3HadTLV;
    muTau_TLV[3]=Tau4HadTLV;//(MM: Added the other tau)
    muTau_TLV[4]=Muon1HadTLV;    
    muTau_TLV[5]=Muon2HadTLV;
    muTau_TLV[6]=Muon3HadTLV;
    muTau_TLV[7]=Muon4HadTLV;//(MM: Added the other muon)
    
    double muTau_mass_i = 0.;
    
    for (int i = 0; i<4; i++){
      if((muTau_TLV[i].Pt() > tau_pt_cut) && (abs(muTau_TLV[i].Eta()) < tau_eta_cut)){
        for(int j = 4; j<8; j++){
  	  if((muTau_TLV[j].Pt() > muon_pt_cut) && (abs(muTau_TLV[j].Eta() < muon_eta_cut))){
	    double muTau_mass = (muTau_TLV[i]+muTau_TLV[j]).M();
	    if(muTau_mass > muTau_mass_i){
              muTau_mass_i = muTau_mass;
 	      pairs_muTau_TLV[0] = muTau_TLV[i];
	      pairs_muTau_TLV[1] = muTau_TLV[j];
	    }
	  }
        }
      }
    }
    //(MM: mu-mu pairs)
    muMu_TLV[0]=Muon1HadTLV;
    muMu_TLV[1]=Muon2HadTLV;
    muMu_TLV[2]=Muon3HadTLV;
    muMu_TLV[3]=Muon4HadTLV;
    
    double muMu_mass_i = 0.;
    
    for (int i = 0; i<4; i++){
      if((muMu_TLV[i].Pt() > muon_pt_cut) && (abs(muMu_TLV[i].Eta()) < muon_eta_cut)){
        for(int j = i+1; j<4; j++){
  	  if((muMu_TLV[j].Pt() > muon_pt_cut) && (abs(muMu_TLV[j].Eta()) < muon_eta_cut) && (i!=j)){
	    double muMu_mass = (muMu_TLV[i]+muMu_TLV[j]).M();
	    if(muMu_mass > muMu_mass_i){
              muMu_mass_i = muMu_mass;
 	      pairs_muMu_TLV[0] = muMu_TLV[i];
	      pairs_muMu_TLV[1] = muMu_TLV[j];
	    }
	  }
        }
      }
    }
    //(MM: Tau-Tau pairs)
    tauTau_TLV[0]=Tau1HadTLV;
    tauTau_TLV[1]=Tau2HadTLV;
    tauTau_TLV[2]=Tau3HadTLV;
    tauTau_TLV[3]=Tau4HadTLV;
    
    double tauTau_mass_i = 0.;
    
    for (int i = 0; i<4; i++){
      if((tauTau_TLV[i].Pt() > tau_pt_cut) && (abs(tauTau_TLV[i].Eta()) < tau_eta_cut)){
        for(int j = i+1; j<4; j++){
  	  if((tauTau_TLV[j].Pt() > tau_pt_cut) && (abs(tauTau_TLV[j].Eta()) < tau_eta_cut) && (i!=j)){
	    double tauTau_mass = (tauTau_TLV[i]+tauTau_TLV[j]).M();
	    if(tauTau_mass > tauTau_mass_i){
              tauTau_mass_i = tauTau_mass;
 	      pairs_tauTau_TLV[0] = tauTau_TLV[i];
	      pairs_tauTau_TLV[1] = tauTau_TLV[j];
	    }
	  }
        }
      }
    }
    
  /////////////////////////////////// Nestor ///////////////////////////////////////
    // elec - tau pairs
    elecTau_TLV[0]=Elec1HadTLV;
    elecTau_TLV[1]=Elec2HadTLV;
    elecTau_TLV[2]=Elec3HadTLV;
    elecTau_TLV[3]=Elec4HadTLV;
    elecTau_TLV[4]=Tau1HadTLV;
    elecTau_TLV[5]=Tau2HadTLV;
    elecTau_TLV[6]=Tau3HadTLV;
    elecTau_TLV[7]=Tau4HadTLV;
    
    double elecTau_mass_i = 0.;
    
    for (int i = 0; i<4; i++){
      if((elecTau_TLV[i].Pt() > elec_pt_cut) && (abs(elecTau_TLV[i].Eta()) < elec_eta_cut)){
        for(int j = 4; j<8; j++){
  	  if((elecTau_TLV[j].Pt() > tau_pt_cut) && (abs(elecTau_TLV[j].Eta()) < tau_eta_cut) && (i!=j)){
	    double elecTau_mass = (elecTau_TLV[i]+elecTau_TLV[j]).M();
	    if(elecTau_mass > elecTau_mass_i){
              elecTau_mass_i = elecTau_mass;
 	      pairs_elecTau_TLV[0] = elecTau_TLV[i];
	      pairs_elecTau_TLV[1] = elecTau_TLV[j];
	    }
	  }
        }
      }
    }
    
    // elec - elec pairs
    elecElec_TLV[0]=Elec1HadTLV;
    elecElec_TLV[1]=Elec2HadTLV;
    elecElec_TLV[2]=Elec3HadTLV;
    elecElec_TLV[3]=Elec4HadTLV;
    
    double elecElec_mass_i = 0.;
    
    for (int i = 0; i<4; i++){
      if((elecElec_TLV[i].Pt() > elec_pt_cut) && (abs(elecElec_TLV[i].Eta()) < elec_eta_cut)){
        for(int j = i+1; j<4; j++){
  	  if((elecElec_TLV[j].Pt() > elec_pt_cut) && (abs(elecElec_TLV[j].Eta()) < elec_eta_cut) && (i!=j)){
	    double elecElec_mass = (elecElec_TLV[i]+elecElec_TLV[j]).M();
	    if(elecElec_mass > elecElec_mass_i){
              elecElec_mass_i = elecElec_mass;
 	      pairs_elecElec_TLV[0] = elecElec_TLV[i];
	      pairs_elecElec_TLV[1] = elecElec_TLV[j];
	    }
	  }
        }
      }
    }

    // mu-mu-mu trio
    muMuMu_TLV[0]=Muon1HadTLV;
    muMuMu_TLV[1]=Muon2HadTLV;
    muMuMu_TLV[2]=Muon3HadTLV;
    muMuMu_TLV[3]=Muon4HadTLV;
    
    double muMuMu_mass_i = 0.;
    
    for (int i = 0; i<4; i++){
      if((muMuMu_TLV[i].Pt() > muon_pt_cut) && (abs(muMuMu_TLV[i].Eta()) < muon_eta_cut)){
        for(int j = i+1; j<4; j++){
  	  if((muMuMu_TLV[j].Pt() > muon_pt_cut) && (abs(muMuMu_TLV[j].Eta()) < muon_eta_cut) && (i!=j)){
            for(int k = j+1; k<4; k++){
	      double muMuMu_mass = (muMuMu_TLV[i]+muMuMu_TLV[j]+muMuMu_TLV[k]).M();
	      if(muMuMu_mass > muMuMu_mass_i){
                muMuMu_mass_i = muMuMu_mass;
 	        trio_muMuMu_TLV[0] = muMuMu_TLV[i];
	        trio_muMuMu_TLV[1] = muMuMu_TLV[j];
	        trio_muMuMu_TLV[2] = muMuMu_TLV[k];
	      }
            }
	  }
        }
      }
    }
    
    // elec-elec-elec trio
    elElEl_TLV[0]=Elec1HadTLV;
    elElEl_TLV[1]=Elec2HadTLV;
    elElEl_TLV[2]=Elec3HadTLV;
    elElEl_TLV[3]=Elec4HadTLV;
    
    double elElEl_mass_i = 0.;
    
    for (int i = 0; i<4; i++){
      if((elElEl_TLV[i].Pt() > elec_pt_cut) && (abs(elElEl_TLV[i].Eta()) < elec_eta_cut)){
        for(int j = i+1; j<4; j++){
  	  if((elElEl_TLV[j].Pt() > elec_pt_cut) && (abs(elElEl_TLV[j].Eta()) < elec_eta_cut) && (i!=j)){
            for(int k = j+1; k<4; k++){
	      double elElEl_mass = (elElEl_TLV[i]+elElEl_TLV[j]+elElEl_TLV[k]).M();
	      if(elElEl_mass > elElEl_mass_i){
                elElEl_mass_i = elElEl_mass;
 	        trio_elElEl_TLV[0] = elElEl_TLV[i];
	        trio_elElEl_TLV[1] = elElEl_TLV[j];
	        trio_elElEl_TLV[2] = elElEl_TLV[k];
	      }
            }
	  }
        }
      }
    }
    
  //////////////////////////////////////////////////////////////////////////////////
    
  //////////////////////////////////////////////////////////////////////////////////
    
    //////// Apply cuts /////////
    // Events with no cuts
    pass_cuts[0] = 1;
    
    // tau: Number of bjets
    if((nBJets == 0) ){
      pass_cuts[1] = 1;
    }
    // tau: dijet kinematics
    if((pass_cuts[1] == 1) && (jetLeadingVec.Pt() > jet_min_pt) && (jetSleadingVec.Pt() > jet_min_pt ) && 
       (abs(jetLeadingVec.Eta()) <jet_max_eta ) && (abs(jetSleadingVec.Eta()) <jet_max_eta )){
      pass_cuts[2] = 1;
    }
    //tau: MET cut
    if ((pass_cuts[2] == 1) && (MET > MET_cut) ){
      pass_cuts[3] = 1;
    }
    // tau: deltaEta cut
    if ((pass_cuts[3] == 1) && ( delta_eta_diJet > deltaEta_diJet_cut)){
      pass_cuts[4] = 1;
    }
    // tau: diJet mass cut 
    if((pass_cuts[4] == 1) && ( DiJetMass_final > diJetmass_cut)){
      pass_cuts[5] = 1;
    }
    // tau: Events with 1 tau with kinematic cuts
    if ((pass_cuts[5] == 1) && (ntau_counter==1) && (Tau1HadTLV.Pt() > tau_pt_cut) && (Tau1HadTLV.Pt() < tau_pt_cut_max) && 
        (abs(Tau1HadTLV.Eta()) < tau_eta_cut) ){
      pass_cuts[6] = 1;
    }
    //muon: Events with 1 muon  with kinematic cuts
    if ((pass_cuts[5] == 1) && (nmuon_counter==1) && (Muon1HadTLV.Pt() > muon_pt_cut)  && (abs(Muon1HadTLV.Eta()) < muon_eta_cut) ){
      pass_cuts[7] = 1;
    }
    //electron: Events with 1 electron with kinematic cuts
    if ((pass_cuts[5] == 1) && (nelec_counter==1) && (Elec1HadTLV.Pt() > elec_pt_cut)  && (abs(Elec1HadTLV.Eta()) < elec_eta_cut) ){
      pass_cuts[8] = 1;
    }
    //muTau kinematics
    if ((pass_cuts[5] == 1) && (nmuon_counter==1) && (ntau_counter==1) && (pairs_muTau_TLV[0].Pt() >  tau_pt_cut) && (pairs_muTau_TLV[1].Pt()> muon_pt_cut) && 
    ( abs(pairs_muTau_TLV[0].Eta()) < tau_eta_cut) && ( abs(pairs_muTau_TLV[1].Eta()) < muon_eta_cut) && (muTau_mass_i > muTau_mass_input)){
      pass_cuts[9] = 1;
    }
    //muMu kinematics
    if ((pass_cuts[5] == 1) && (nmuon_counter==2) && (pairs_muMu_TLV[0].Pt() >  muon_pt_cut) && (pairs_muMu_TLV[1].Pt()> muon_pt_cut) && 
    ( abs(pairs_muMu_TLV[0].Eta()) < muon_eta_cut) && ( abs(pairs_muMu_TLV[1].Eta()) < muon_eta_cut) && (muMu_mass_i > muMu_mass_input)){
      pass_cuts[10] = 1;
    }
    //tauTau kinematics
    if ((pass_cuts[5] == 1)  && (ntau_counter==2) && (pairs_tauTau_TLV[0].Pt() >  tau_pt_cut) && (pairs_tauTau_TLV[1].Pt()> tau_pt_cut) && 
    ( abs(pairs_tauTau_TLV[0].Eta()) < tau_eta_cut) && ( abs(pairs_tauTau_TLV[1].Eta()) < tau_eta_cut) && (tauTau_mass_i > tauTau_mass_input)){
      pass_cuts[11] = 1;
    }
    //elecTau kinematics
    if ((pass_cuts[5] == 1) && (nelec_counter==1) && (ntau_counter==1) && (pairs_elecTau_TLV[0].Pt() >  elec_pt_cut) && (pairs_elecTau_TLV[1].Pt()> tau_pt_cut) && 
    ( abs(pairs_elecTau_TLV[0].Eta()) < elec_eta_cut) && ( abs(pairs_elecTau_TLV[1].Eta()) < tau_eta_cut) && (elecTau_mass_i > elecTau_mass_input)){
      pass_cuts[12] = 1;
    }
    //elecElec kinematics
    if ((pass_cuts[5] == 1)  && (nelec_counter==2) && (pairs_elecElec_TLV[0].Pt() >  elec_pt_cut) && (pairs_elecElec_TLV[1].Pt()> elec_pt_cut) && 
    ( abs(pairs_elecElec_TLV[0].Eta()) < elec_eta_cut) && ( abs(pairs_elecElec_TLV[1].Eta()) < elec_eta_cut) && (elecElec_mass_i > elecElec_mass_input)){
      pass_cuts[13] = 1;
    }
    //muMuMu kinematics
    if ((pass_cuts[5] == 1)  && (nmuon_counter==3) && (trio_muMuMu_TLV[0].Pt() >  muon_pt_cut) && (trio_muMuMu_TLV[1].Pt()> muon_pt_cut) && 
    (trio_muMuMu_TLV[2].Pt() > muon_pt_cut) && ( abs(trio_muMuMu_TLV[0].Eta()) < muon_eta_cut) && ( abs(trio_muMuMu_TLV[1].Eta()) < muon_eta_cut)    && ( abs(trio_muMuMu_TLV[2].Eta()) < muon_eta_cut)){
      pass_cuts[14] = 1;
    }
    //elElEl kinematics
    if ((pass_cuts[5] == 1)  && (nelec_counter==3) && (trio_elElEl_TLV[0].Pt() >  elec_pt_cut) && (trio_elElEl_TLV[1].Pt()> elec_pt_cut) && 
    (trio_elElEl_TLV[2].Pt() > elec_pt_cut) && ( abs(trio_elElEl_TLV[0].Eta()) < elec_eta_cut) && ( abs(trio_elElEl_TLV[1].Eta()) < elec_eta_cut)    && ( abs(trio_elElEl_TLV[2].Eta()) < elec_eta_cut)){
      pass_cuts[15] = 1;
    }

    //Fill histograms (MM: No sé como llenar los histogramas de masa...)
    for (int i = 0; i < nDir; i++){
      
      _hmap_Nevents[i]->Fill(0.0);
      
      if (pass_cuts[i] == 1){
	
        _hmap_n_tau[i]->Fill(ntau_counter);
        _hmap_n_jets[i]->Fill(n_jets);
        _hmap_Nevents[i]->Fill(1.0);
	
        if(jetLeadingVec.Pt() > 1.0){
	  _hmap_lead_jet_pT[i]->Fill(jetLeadingVec.Pt());
	  _hmap_lead_jet_eta[i]->Fill(jetLeadingVec.Eta());
	  _hmap_lead_jet_phi[i]->Fill(jetLeadingVec.Phi());
	  _hmap_slead_jet_pT[i]->Fill(jetSleadingVec.Pt());
	  _hmap_slead_jet_eta[i]->Fill(jetSleadingVec.Eta());
	  _hmap_slead_jet_phi[i]->Fill(jetSleadingVec.Phi());
        }
	if(Tau1HadTLV.Pt() > 1.0){
	  _hmap_tau1_pT[i]->Fill(Tau1HadTLV.Pt());
	  _hmap_tau1_eta[i]->Fill(Tau1HadTLV.Eta());
	  _hmap_tau1_phi[i]->Fill(Tau1HadTLV.Phi());
        }
        if(Tau2HadTLV.Pt() > 1.0){
	  _hmap_tau2_pT[i]->Fill(Tau2HadTLV.Pt());
	  _hmap_tau2_eta[i]->Fill(Tau2HadTLV.Eta());
	  _hmap_tau2_phi[i]->Fill(Tau2HadTLV.Phi());
          _hmap_Delta_pT[i]->Fill(abs(Tau1HadTLV.Pt() -Tau2HadTLV.Pt()));
        }
	if(Muon1HadTLV.Pt() > 1.0){
	  _hmap_muon1_pT[i]->Fill(Muon1HadTLV.Pt());
	  _hmap_muon1_eta[i]->Fill(Muon1HadTLV.Eta());
	  _hmap_muon1_phi[i]->Fill(Muon1HadTLV.Phi());
        }
        if(Muon2HadTLV.Pt() > 1.0){
	  _hmap_muon2_pT[i]->Fill(Muon2HadTLV.Pt());
	  _hmap_muon2_eta[i]->Fill(Muon2HadTLV.Eta());
	  _hmap_muon2_phi[i]->Fill(Muon2HadTLV.Phi());
        }
	if(Elec1HadTLV.Pt() > 1.0){
	  _hmap_elec1_pT[i]->Fill(Elec1HadTLV.Pt());
	  _hmap_elec1_eta[i]->Fill(Elec1HadTLV.Eta());
	  _hmap_elec1_phi[i]->Fill(Elec1HadTLV.Phi());
        }
        if(Elec2HadTLV.Pt() > 1.0){
	  _hmap_elec2_pT[i]->Fill(Elec2HadTLV.Pt());
	  _hmap_elec2_eta[i]->Fill(Elec2HadTLV.Eta());
	  _hmap_elec2_phi[i]->Fill(Elec2HadTLV.Phi());
        }
        if (isSignal == 1){
          if(Tau1HadTLV_gen.Pt() > 1.0){
            _hmap_Gentau1_pT[i]->Fill(Tau1HadTLV_gen.Pt());
            _hmap_Gentau1_eta[i]->Fill(Tau1HadTLV_gen.Eta());
            _hmap_Gentau1_phi[i]->Fill(Tau1HadTLV_gen.Phi());
          }
          if(Tau2HadTLV_gen.Pt() > 1.0){
            _hmap_Gentau2_pT[i]->Fill(Tau2HadTLV_gen.Pt());
            _hmap_Gentau2_eta[i]->Fill(Tau2HadTLV_gen.Eta());
            _hmap_Gentau2_phi[i]->Fill(Tau2HadTLV_gen.Phi());
            _hmap_GenDelta_pT[i]->Fill(abs(Tau1HadTLV_gen.Pt() -Tau2HadTLV_gen.Pt()));
          }
        } else{
       
          if(Muon1HadTLV.Pt() > 1.0){
            _hmap_Gentau1_pT[i]->Fill(Muon1HadTLV.Pt());
            _hmap_Gentau1_eta[i]->Fill(Muon1HadTLV.Eta());
            _hmap_Gentau1_phi[i]->Fill(Muon1HadTLV.Phi());
          }
          if(Muon2HadTLV.Pt() > 1.0){
            _hmap_Gentau2_pT[i]->Fill(Muon2HadTLV.Pt());
            _hmap_Gentau2_eta[i]->Fill(Muon2HadTLV.Eta());
            _hmap_Gentau2_phi[i]->Fill(Muon2HadTLV.Phi());
            _hmap_GenDelta_pT[i]->Fill(abs(Muon1HadTLV.Pt() -Muon2HadTLV.Pt()));
          }

        }
        if(tauMass > 0){_hmap_tauMass[i]->Fill(tauMass);}
        if(MET > 0){_hmap_MET[i]->Fill(MET);}
        if(tau_transmass > 0){_hmap_tau_transmass[i]->Fill(tau_transmass);}
        if(muon_transmass > 0){_hmap_muon_transmass[i]->Fill(muon_transmass);}
        if(elec_transmass > 0){_hmap_elec_transmass[i]->Fill(elec_transmass);}
	if(ht > 0.0){_hmap_ht[i]->Fill(ht);}
        if(st > 0.0){_hmap_st[i]->Fill(st);}
	if(DiJetMass_final > 100.0){_hmap_dijet_mass[i]->Fill(DiJetMass_final);}
	if(delta_eta_diJet > 0){_hmap_dijet_deltaEta[i]->Fill(delta_eta_diJet);}
      }
    }
  }// end entry loop for tau channel

  theFile->cd();
  for (int d = 0; d < nDir; d++)
    {
      cdDir[d]->cd();
      _hmap_Nevents[d]->Write();
      _hmap_n_jets[d]->Write();
      _hmap_n_tau[d]->Write();
      _hmap_lead_jet_pT[d]->Write();
      _hmap_lead_jet_eta[d]->Write();
      _hmap_lead_jet_phi[d]->Write();
      _hmap_slead_jet_pT[d]->Write();
      _hmap_slead_jet_eta[d]->Write();
      _hmap_slead_jet_phi[d]->Write();
      _hmap_tau1_pT[d]->Write();
      _hmap_tau1_eta[d]->Write();
      _hmap_tau1_phi[d]->Write();
      _hmap_tau2_pT[d]->Write();
      _hmap_tau2_eta[d]->Write();
      _hmap_tau2_phi[d]->Write();
      _hmap_muon1_pT[d]->Write();
      _hmap_muon1_eta[d]->Write();
      _hmap_muon1_phi[d]->Write();
      _hmap_muon2_pT[d]->Write();
      _hmap_muon2_eta[d]->Write();
      _hmap_muon2_phi[d]->Write();
      _hmap_elec1_pT[d]->Write();
      _hmap_elec1_eta[d]->Write();
      _hmap_elec1_phi[d]->Write();
      _hmap_elec2_pT[d]->Write();
      _hmap_elec2_eta[d]->Write();
      _hmap_elec2_phi[d]->Write();
      _hmap_Delta_pT[d]->Write();
      //Gen/////
      _hmap_Gentau1_pT[d]->Write();
      _hmap_Gentau1_eta[d]->Write();
      _hmap_Gentau1_phi[d]->Write();
      _hmap_Gentau2_pT[d]->Write();
      _hmap_Gentau2_eta[d]->Write();
      _hmap_Gentau2_phi[d]->Write();
      _hmap_st[d]->Write();
      _hmap_GenDelta_pT[d]->Write();
       ////////
      _hmap_tauMass[d]->Write();
      _hmap_MET[d]->Write();
      _hmap_tau_transmass[d]->Write();
      _hmap_muon_transmass[d]->Write();
      _hmap_elec_transmass[d]->Write();
      _hmap_ht[d]->Write();
      _hmap_dijet_mass[d]->Write();
      _hmap_dijet_deltaEta[d]->Write();

    }
  theFile->Close();
}

PhenoAnalysis::~PhenoAnalysis()
{
  // do anything here that needs to be done at desctruction time
}

double PhenoAnalysis::calculateE(double eta, double pt, double mass){

  double theta = TMath::ATan(TMath::Exp(-eta));
  double sin_theta = TMath::Sin(2*theta);
  double p= pt/sin_theta;
  double e = sqrt(pow(p, 2) + pow(mass, 2));

  return e;
}

double PhenoAnalysis::calculate_deltaR(TLorentzVector vector,  Track* track){

  double eta1 = vector.Eta();
  double phi1 = vector.Phi();
  double eta2 = track->Eta;
  double phi2 = track->Phi;
  double deltaR = sqrt(pow(eta1-eta2,2) + pow(phi1-phi2,2));
  return deltaR;
}

double PhenoAnalysis::normalizedDphi(double phi){
  const double PI  = 3.141592653589793238463;
  double twoPI = 2.0*PI;
  if ( phi < -PI ){phi += twoPI;}
  if ( phi > PI ){phi = twoPI-phi;}
  else phi = TMath::Abs(phi);
  return phi;
}
void PhenoAnalysis::crateHistoMasps (int directories)
{
  for (int i = 0; i < directories; i++)
    {
      _hmap_Nevents[i]       = new TH1F("Nevents", "Nevents", 3,0,3);
      _hmap_lead_jet_pT[i]   = new TH1F("jetLeadPt",    "p_{T} Leading Jet", 200, 0., 2000.);
      _hmap_lead_jet_eta[i]  = new TH1F("jetLeadEta",   "#eta Leading Jet", 50, -5.0, 5.0);
      _hmap_lead_jet_phi[i]  = new TH1F("jetLeadPhi",   "#phi Leading Jet", 70, -3.6, 3.6);
      _hmap_slead_jet_pT[i]   = new TH1F("jetSleadPt",    "p_{T} Sub-leading Jet", 200, 0., 2000.);
      _hmap_slead_jet_eta[i]  = new TH1F("jetSleadEta",   "#eta Sub-leading Jet", 50, -5.0, 5.0);
      _hmap_slead_jet_phi[i]  = new TH1F("jetSleadPhi",   "#phi Sub-leading Jet", 70, -3.6, 3.6);
      _hmap_n_jets[i]        = new TH1F("nJets",         "N(jet)", 7, 0, 7);
      _hmap_n_tau[i]         = new TH1F("nTaus",          "N(#tau)", 4, 0, 4);
      _hmap_tau1_pT[i]       = new TH1F("tau1Pt",        "p_{T}(#tau_{1})", 200, 0., 2000.);
      _hmap_tau1_eta[i]      = new TH1F("tau1Eta",       "#eta(#tau_{1})", 50, -3.5, 3.5);
      _hmap_tau1_phi[i]      = new TH1F("tau1Phi",       "#phi(#tau_{1})", 70, -3.6, 3.6);
      _hmap_tau2_pT[i]       = new TH1F("tau2Pt",        "p_{T}(#tau_{2})", 200, 0., 2000.);
      _hmap_tau2_eta[i]      = new TH1F("tau2Eta",       "#eta(#tau_{2})", 50, -3.5, 3.5);
      _hmap_tau2_phi[i]      = new TH1F("tau2Phi",       "#phi(#tau_{2})", 70, -3.6, 3.6);
      _hmap_muon1_pT[i]       = new TH1F("muon1Pt",        "p_{T}(#muon_{1})", 200, 0., 2000.);
      _hmap_muon1_eta[i]      = new TH1F("muon1Eta",       "#eta(#muon_{1})", 50, -3.5, 3.5);
      _hmap_muon1_phi[i]      = new TH1F("muon1Phi",       "#phi(#muon_{1})", 70, -3.6, 3.6);
      _hmap_muon2_pT[i]       = new TH1F("muon2Pt",        "p_{T}(#muon_{2})", 200, 0., 2000.);
      _hmap_muon2_eta[i]      = new TH1F("muon2Eta",       "#eta(#muon_{2})", 50, -3.5, 3.5);
      _hmap_muon2_phi[i]      = new TH1F("muon2Phi",       "#phi(#muon_{2})", 70, -3.6, 3.6);
      _hmap_elec1_pT[i]       = new TH1F("elec1Pt",        "p_{T}(#elec_{1})", 200, 0., 2000.);
      _hmap_elec1_eta[i]      = new TH1F("elec1Eta",       "#eta(#elec_{1})", 50, -3.5, 3.5);
      _hmap_elec1_phi[i]      = new TH1F("elec1Phi",       "#phi(#elec_{1})", 70, -3.6, 3.6);
      _hmap_elec2_pT[i]       = new TH1F("elec2Pt",        "p_{T}(#elec_{2})", 200, 0., 2000.);
      _hmap_elec2_eta[i]      = new TH1F("elec2Eta",       "#eta(#elec_{2})", 50, -3.5, 3.5);
      _hmap_elec2_phi[i]      = new TH1F("elec2Phi",       "#phi(#elec_{2})", 70, -3.6, 3.6);
      _hmap_tauMass[i]       = new TH1F("tauMass", "m(#t)", 100, 0, 1000);
      _hmap_MET[i]           = new TH1F("MET", "MET", 100, 0, 1000);
      _hmap_tau_transmass[i] = new TH1F("tauTransMass", "Tau Transverse Mass", 100, 0, 1000);
      _hmap_muon_transmass[i] = new TH1F("muonTransMass", "Muon Transverse Mass", 100, 0, 1000);
      _hmap_elec_transmass[i] = new TH1F("elecTransMass", "Electron Transverse Mass", 100, 0, 1000);
      _hmap_ht[i]            = new TH1F("HT", "H_{T}", 100, 0, 5000);
      _hmap_st[i]            = new TH1F("ST", "S_{T}", 100, 0, 5000);
      _hmap_dijet_mass[i]     = new TH1F("diJetMass", "diJet_Mass", 100, 0, 5000);
      _hmap_dijet_deltaEta[i] = new TH1F("diJetDeltaEta", "diJet_deltaEta", 160, 0, 8);
      // Gen Taus
      _hmap_Gentau1_pT[i]       = new TH1F("Gentau1Pt",        "p_{T}(#tau_{1})", 200, 0., 2000.);
      _hmap_GenDelta_pT[i]      = new TH1F("GenDeltaPt",        "#Delta p_{T}", 1000, 0., 1000.);
      _hmap_Delta_pT[i]      = new TH1F("DeltaPt",        "#Delta p_{T}", 1000, 0., 1000.);
      _hmap_Gentau1_eta[i]      = new TH1F("Gentau1Eta",       "#eta(#tau_{1})", 50, -3.5, 3.5);
      _hmap_Gentau1_phi[i]      = new TH1F("Gentau1Phi",       "#phi(#tau_{1})", 70, -3.6, 3.6);
      _hmap_Gentau2_pT[i]       = new TH1F("Gentau2Pt",        "p_{T}(#tau_{2})", 200, 0., 2000.);
      _hmap_Gentau2_eta[i]      = new TH1F("Gentau2Eta",       "#eta(#tau_{2})", 50, -3.5, 3.5);
      _hmap_Gentau2_phi[i]      = new TH1F("Gentau2Phi",       "#phi(#tau_{2})", 70, -3.6, 3.6);
    }
}
