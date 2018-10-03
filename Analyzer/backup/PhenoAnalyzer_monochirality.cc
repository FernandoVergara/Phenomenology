////////////////////////////////////////////////////////////////
//                                                            //
// Author: Andrés Flórez, Universidad de los Andes, Colombia  //
//                                                            //
////////////////////////////////////////////////////////////////


#include <iostream>
#include "ROOTFunctions.h"
#include "PhenoAnalyzer.h"
#include "DelphesFunctions.h"

int main(int argc, char *argv[]) {
  
  TChain chain("Delphes");
  chain.Add(argv[1]);
  TFile * HistoOutputFile = new TFile(argv[2], "RECREATE");
  int nDir = 5;
  TDirectory *theDirectory[nDir+1];
  theDirectory[0]  = HistoOutputFile->mkdir("No_cuts");
  theDirectory[1]  = HistoOutputFile->mkdir("After_bjetpt");
  theDirectory[2]  = HistoOutputFile->mkdir("After_bjeteta");
  theDirectory[3]  = HistoOutputFile->mkdir("After_lepton");
  theDirectory[4]  = HistoOutputFile->mkdir("After_passed_w");
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
  
  double DR_jet_elec_max  = params->GetValue ("DR_jet_elec_max", 0.); 
  double DR_elec_muon_max = params->GetValue ("DR_elec_muon_max", 0.);
  double jet_pt_min       = params->GetValue ("jet_pt_min", 0.);
  double jet_eta_max      = params->GetValue ("jet_eta_max", 0.);
  double b_pt_min         = params->GetValue ("b_pt_min", 0.);
  double b_eta_max        = params->GetValue ("b_eta_max", 0.);
  double elec_pt_min      = params->GetValue ("elec_pt_min", 0.);
  double muon_pt_min      = params->GetValue ("muon_pt_min", 0.);
  double met_min          = params->GetValue ("met_min", 0.);

 
  crateHistoMasps(nDir);
  
  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  Long64_t numberOfEntries = treeReader->GetEntries();
  
  TClonesArray *branchJet = treeReader->UseBranch("Jet");
  TClonesArray *branchElectron = treeReader->UseBranch("Electron");
  TClonesArray *branchMuon = treeReader->UseBranch("Muon");
  TClonesArray *branchMissingET = treeReader->UseBranch("MissingET"); 
  MissingET *METpointer; 
  
  
  for(Int_t entry = 0; entry < numberOfEntries; ++entry)
    {
      treeReader->ReadEntry(entry);
      TLorentzVector muon_i(0., 0., 0., 0.);
      TLorentzVector jet_i(0., 0., 0., 0.);
      TLorentzVector jet_1(0., 0., 0., 0.);
      TLorentzVector jet_2(0., 0., 0., 0.);
      TLorentzVector elec_i(0., 0., 0., 0.);
      TLorentzVector b_vec(0.,0.,0.,0.);
      TLorentzVector Tau_i(0.,0.,0.,0.);
      
      int pass_cuts[nDir] = {0};
      
      vector<int> tau_index; 
      
      bool is_b_jet = false;
      bool is_jet_elec_overlap = false;
      bool is_jet_muon_overlap = false;
      bool filled_first_jet = false;
      bool passed_mass_w = false;
      
      
      int bjet_counter=0;
      int elec_counter=0;
      int muon_counter=0;
      int tau_counter=0;
      int jet_counter = 0;
      
      double bpt=0.; 
      double jet_met_dphi=0.;
      float  e_pt_min = 10.0;
      float  m_pt_min = 10.0;
      
      METpointer = (MissingET*) branchMissingET->At(0);
      double MET = METpointer->MET;
      double MET_phi = METpointer->Phi;
      
      for (int j = 0; j < branchJet->GetEntriesFast(); j++)
	{ 
	  bool is_jet_elec_overlap = false;
          // Find any jet passing pT and eTa
	  Jet *jet = (Jet*) branchJet->At(j);
	  double jet_energy = calculateE(jet->Eta, jet->PT, jet->Mass);
	  if((jet->PT > 1.0) && (abs(jet->Eta) < 2.5)){jet_i.SetPtEtaPhiE(jet->PT, jet->Eta, jet->Phi, jet_energy);}
	  
          // Find if there are electrons overlapping with the jets
	  for (int el = 0; el < branchElectron->GetEntriesFast(); el++){
	    Electron *elec = (Electron*) branchElectron->At(el);
	    double elec_energy = calculateE(elec->Eta, elec->PT, 0.000510998902);
	    elec_i.SetPtEtaPhiE(elec->PT, elec->Eta, elec->Phi, elec_energy);
	    double DR_jet_elec = jet_i.DeltaR(elec_i);
	    if ((jet_i.Pt() > 1.0) && (DR_jet_elec < DR_jet_elec_max) && (elec->PT > e_pt_min) && (abs(elec->Eta) < 2.5)){
	      is_jet_elec_overlap = true;
              elec_i.SetPtEtaPhiE(0.0, 0.0, 0.0, 0.0);
	      break;
	    }
       

	  }
          // Find if there are muons overlapping with the jets
	  is_jet_muon_overlap = false; 
	  for (int muo = 0; muo < branchMuon->GetEntriesFast(); muo++){
	    Muon *muon = (Muon*) branchMuon->At(muo);
	    double muon_energy = calculateE(muon->Eta, muon->PT, 0.105658369);
	    muon_i.SetPtEtaPhiE(muon->PT, muon->Eta, muon->Phi, muon_energy);
	    double DR_jet_muon = jet_i.DeltaR(muon_i);
	    if ( (jet_i.Pt() > 1.0) && ( DR_jet_muon < DR_jet_elec_max) && (muon->PT > m_pt_min) && (abs(muon->Eta) < 2.5)){ 
	      is_jet_muon_overlap = true;
              muon_i.SetPtEtaPhiE(0.0, 0.0, 0.0, 0.0);
	     cout << is_jet_muon_overlap << endl;
               // break;
	    }
	  }
	  // If there no jets overlapping with electrons or muons, then fill the info for the b-jets....
	  if((!is_jet_elec_overlap) && (!is_jet_muon_overlap)){
	    if ((jet->PT > 1.0) && (jet->BTag == 1)){
	      b_vec.SetPtEtaPhiE(jet->PT, jet->Eta, jet->Phi, jet_energy);           
	      bjet_counter++;
	    }
	    // select electron
	    for (int el = 0; el < branchElectron->GetEntriesFast(); el++){
	      Electron *elec = (Electron*) branchElectron->At(el);
	      double elec_energy = calculateE(elec->Eta, elec->PT, 0.000510998902);
	      if ((elec->PT > e_pt_min) && (abs(elec->Eta) < 2.5)){
		e_pt_min = elec->PT;
		elec_i.SetPtEtaPhiE(elec->PT, elec->Eta, elec->Phi, elec_energy);
		elec_counter++;
                }
              
	      // select muon
	      for (int muo = 0; muo < branchMuon->GetEntriesFast(); muo++){
                Muon *muon = (Muon*) branchMuon->At(muo);
                double muon_energy = calculateE(muon->Eta, muon->PT, 0.105658369);
                double DR_elec_muon = elec_i.DeltaR(muon_i);
                if ((DR_elec_muon > DR_elec_muon_max) && (muon->PT > m_pt_min) && (abs(muon->Eta) < 2.5)){
                  m_pt_min = muon->PT;
                  muon_i.SetPtEtaPhiE(muon->PT, muon->Eta, muon->Phi, muon_energy);
		  muon_counter++;
                 }
              }
	     } 
	    
            if ((jet->TauTag == 1) && (jet->BTag == 0)){
	      if((jet->PT > 1.0) && (abs(jet->Eta) < 2.5))
		tau_index.push_back(j);
	        tau_counter++;
	      for(int i = 0; i < tau_index.size(); i++)
		{
		  int tau_i = tau_index.at(i);
		  Jet *jet = (Jet*) branchJet->At(tau_i);
		  double jet_energy = calculateE(jet->Eta, jet->PT, jet->Mass);
		  Tau_i.SetPtEtaPhiE(jet->PT, jet->Eta, jet->Phi, jet_energy);
		}		
	    }
	    
	    if ((jet->TauTag == 0) && ((jet->BTag == 0))){
              double jet_energy = calculateE(jet->Eta, jet->PT, jet->Mass);
	      if((filled_first_jet == false) && (jet->PT > jet_pt_min)){
		jet_1.SetPtEtaPhiE(jet->PT, jet->Eta, jet->Phi, jet_energy);
		jet_pt_min = jet->PT;
		jet_counter++;
                filled_first_jet = true;
	      }
              if ((filled_first_jet == true) && (jet->PT > jet_pt_min)){
                jet_2.SetPtEtaPhiE(jet->PT, jet->Eta, jet->Phi, jet_energy);  
                float mass_w = (jet_1 + jet_2).M();
                if ((mass_w < 100.0) && (mass_w > 60.0)){passed_mass_w = true;} 
                jet_counter++;
              }
              
	    }
	    
	  } //is not electron or muon overlap
	

  //cout << elec_counter<<endl;
	  
	}    //finaliza el ciclo sobre las entradas. 
      
       //cout << elec_pt_min << " " << muon_pt_min <<endl; 
      //cuts delphes level	
      pass_cuts[0]=1;
      if(b_vec.Pt() > b_pt_min){pass_cuts[1] = 1;}
      if((TMath::Abs(b_vec.Eta()) < jet_eta_max) && (pass_cuts[1] == 1)){pass_cuts[2] = 1;}
      if(((elec_i.Pt() > elec_pt_min) || (muon_i.Pt() > muon_pt_min)) && (pass_cuts[2] == 1)){pass_cuts[3] = 1;}
      
      if ((pass_cuts[3] == 1) && (passed_mass_w)){pass_cuts[4] = 1;}
      
      //Highest lepton Pt delphes-level
      if(elec_i.Pt() > muon_i.Pt()){
	bpt = b_vec.Pt()/(b_vec.Pt()+elec_i.Pt());
      }else{if(muon_i.Pt() > elec_i.Pt())
	  bpt = b_vec.Pt()/(b_vec.Pt()+muon_i.Pt());
      }
      
      jet_met_dphi=normalizedDphi(jet_i.Phi()-MET_phi);
      
      // save some important histograms
      for (int i = 0; i < nDir; i++){
        if (pass_cuts[i] == 1){
	  _hmap_N_events[i]->Fill(1.);  
	  _hmap_N_jets[i]->Fill(jet_counter);
	  _hmap_N_tau[i]->Fill(tau_counter);
	  _hmap_N_bjet[i]->Fill(bjet_counter);
	  // _hmap_N_elec[i]->Fill(elec_counter);
	  // _hmap_N_muon[i]->Fill(muon_counter);
          
	  if(jet_i.Pt()>0.){ 
	    _hmap_jet_pT[i]->Fill(jet_i.Pt());
	    _hmap_jet_eta[i]->Fill(jet_i.Eta());
	    _hmap_jet_phi[i]->Fill(jet_i.Phi());
	  }
	  
	  if(Tau_i.Pt()>0.){
            _hmap_tau1_pT[i]->Fill(Tau_i.Pt());
            _hmap_tau1_eta[i]->Fill(Tau_i.Eta());
            _hmap_tau1_phi[i]->Fill(Tau_i.Phi());
	  }
	  
	  if(b_vec.Pt()>0.){
	    _hmap_bjet_ratio_pT[i]->Fill(bpt);	
	    _hmap_bjet_pT[i]->Fill(b_vec.Pt());
	    _hmap_bjet_eta[i]->Fill(b_vec.Eta());
	    _hmap_bjet_phi[i]->Fill(b_vec.Phi());
	    _hmap_bjet_energy[i]->Fill(b_vec.E()*0.002);
	  }
	  
	  if (elec_i.Pt() > 0.0) { 
	    _hmap_elec_pT[i]->Fill(elec_i.Pt());
	    _hmap_elec_eta[i]->Fill(elec_i.Eta());
	    _hmap_elec_phi[i]->Fill(elec_i.Phi());
	    _hmap_N_elec[i]->Fill(elec_counter);
	  }
	  if (muon_i.Pt() > 0.0) {
	    _hmap_N_muon[i]->Fill(muon_counter);
	    _hmap_muon_pT[i]->Fill(muon_i.Pt());
	    _hmap_muon_eta[i]->Fill(muon_i.Eta());
	    _hmap_muon_phi[i]->Fill(muon_i.Phi());
	  }
	  _hmap_met[i]->Fill(MET);
	  _hmap_met_jetmetdphi[i]->Fill(jet_met_dphi,MET);
	} 
      }
      
    }
  
  
  theFile->cd();
  for (int d = 0; d < nDir; d++)
    {
      cdDir[d]->cd();
      _hmap_N_events[d]->Write();
      _hmap_N_jets[d]->Write();
      _hmap_jet_pT[d]->Write();
      _hmap_jet_eta[d]->Write();
      _hmap_jet_phi[d]->Write();
      _hmap_N_tau[d]->Write();
      _hmap_tau1_pT[d]->Write();
      _hmap_tau1_eta[d]->Write();
      _hmap_tau1_phi[d]->Write();
      _hmap_N_bjet[d]->Write();
      _hmap_bjet_ratio_pT[d]->Write();
      _hmap_bjet_pT[d]->Write();
      _hmap_bjet_eta[d]->Write();
      _hmap_bjet_phi[d]->Write();
      _hmap_bjet_energy[d]->Write();
      _hmap_N_elec[d]->Write();
      _hmap_elec_pT[d]->Write();
      _hmap_elec_eta[d]->Write();
      _hmap_elec_phi[d]->Write();
      _hmap_N_muon[d]->Write();
      _hmap_muon_pT[d]->Write();
      _hmap_muon_eta[d]->Write();
      _hmap_muon_phi[d]->Write();
      _hmap_met[d]->Write();
      _hmap_met_jetmetdphi[d]->Write();    
    }
  
  theFile->Close();
  
}

PhenoAnalysis::~PhenoAnalysis()
{
  // do anything here that needs to be done at desctruction time
}

double PhenoAnalysis::calculateE(double eta, double pt, double mass){
  
  double theta = TMath::ATan(TMath::Exp(-eta)); 
  double cos_theta = TMath::Cos(2*theta);
  double p= pt/cos_theta;
  double e = sqrt(pow(p, 2) + pow(mass, 2));
  
  return e;
  
}


double PhenoAnalysis::normalizedDphi(double phi){
  double twoPI = 2.0*TMath::Pi();
  if ( phi < -TMath::Pi() ){phi += twoPI;}
  if ( phi > TMath::Pi() ){phi = twoPI-phi;}
  else phi = abs(phi); 
  return phi;
}



void PhenoAnalysis::crateHistoMasps (int directories)
{
  for (int i = 0; i < directories; i++)
    {
      _hmap_N_events[i]		= new TH1F("N_events", "N_events", 2,-0.5,1.5); 
      
      _hmap_N_jets[i]		= new TH1F("N_jets",  "N_jets", 11, -0.5, 10.5);    
      _hmap_jet_pT[i]		= new TH1F("jet_pT",  "P_{T} jet ", 100, 0., 1000.);
      _hmap_jet_eta[i]		= new TH1F("jet_eta", "#eta jet ", 50, -5.0, 5.0);
      _hmap_jet_phi[i]		= new TH1F("jet_phi", "#phi jet ", 70, -3.6, 3.6);
      
      _hmap_N_tau[i]		= new TH1F("N_tau",  "N_tau", 5, -0.5, 4.5);
      _hmap_tau1_pT[i]		= new TH1F("tau1_pT",  "p_{T}(#tau_{1})", 100, 0., 1000.);
      _hmap_tau1_eta[i]		= new TH1F("tau1_eta", "#eta(#tau_{1})", 50, -5.0, 5.0);
      _hmap_tau1_phi[i]		= new TH1F("tau1_phi", "#phi(#tau_{1})", 72, -3.6, 3.6);
      
      _hmap_N_bjet[i]		= new TH1F("N_bjet",    "N_bjet", 5, -0.5, 4.5);
      _hmap_bjet_ratio_pT[i]	= new TH1F("bjet_ratio_pT", "p_{T} bjet", 50, 0, 1);
      _hmap_bjet_pT[i]		= new TH1F("bjet_pT",  "P_{T} bjet ", 100, 0., 1000.);
      _hmap_bjet_eta[i]		= new TH1F("bjet_eta", "#eta bjet ", 50, -5.0, 5.0);
      _hmap_bjet_phi[i]		= new TH1F("bjet_phi", "#phi bjet ", 70, -3.6, 3.6);
      _hmap_bjet_energy[i]	= new TH1F("bjet_energy", "Energy b quark", 50, 0, 2);
      
      _hmap_N_elec[i]           = new TH1F("N_elec",    "N_elec", 11, -0.5, 10.5);
      _hmap_elec_pT[i]          = new TH1F("elec_pT", "p_{T} e ", 100, 0., 1000.);
      _hmap_elec_eta[i]         = new TH1F("elec_eta","#eta e ", 50, -5.0, 5.0);
      _hmap_elec_phi[i]         = new TH1F("elec_phi", "#phi e", 70, -3.6, 3.6);      
      
      _hmap_N_muon[i]           = new TH1F("N_muon",    "N_muon", 11, -0.5, 10.5);
      _hmap_muon_pT[i]          = new TH1F("muon_pT", "p_{T} #mu ", 100, 0., 1000.);
      _hmap_muon_eta[i]         = new TH1F("muon_eta","#eta #mu ", 50, -5.0, 5.0);
      _hmap_muon_phi[i]         = new TH1F("muon_phi", "#phi #mu", 70, -3.6, 3.6); 
      
      _hmap_met[i]		= new TH1F("met", "MET", 200, 0. , 2000.);	
      _hmap_met_jetmetdphi[i]	= new TH2F("met_jetmetdphi", "#Delta #phi(jet,MET), MET", 32, 0., TMath::Pi(), 50, 0., 500.);
    }
}

